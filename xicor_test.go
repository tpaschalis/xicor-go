package xicor

import (
	"fmt"
	"math"
	"math/rand"
	"reflect"
	"testing"
	"time"
)

var anscombesQuartet = map[string][]float64{
	"x_1": {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5},
	"x_2": {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5},
	"x_3": {0, 0, 0, 0, 1, 1, 2, 3},
	"x_4": {8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8},
	"y_1": {8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68},
	"y_2": {9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74},
	"y_3": {0, 1, 2, 3, 4, 5, 6, 7},
	"y_4": {6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89},
}

// Test correctness of results -- results compared with initial code by Sourav Chatterjee and Python port by `czbiohub`
func TestXi(t *testing.T) {
	rand.Seed(21)

	// Asymptotic p-value calculation
	xi1, pval1, _ := New(
		anscombesQuartet["x_1"],
		anscombesQuartet["y_1"],
		WithAsymptoticPvalue(),
	).Pvalue()
	assertEpsilon(t, xi1, 0.275)
	assertEpsilon(t, pval1, 0.07841556)

	xi2, pval2, _ := New(
		anscombesQuartet["x_2"],
		anscombesQuartet["y_2"],
		WithAsymptoticPvalue(),
	).Pvalue()
	assertEpsilon(t, xi2, 0.6)
	assertEpsilon(t, pval2, 0.0010040217037570187)

	xi3, pval3, _ := New(
		anscombesQuartet["x_3"],
		anscombesQuartet["y_3"],
		WithAsymptoticPvalue(),
	).Pvalue()
	assertEpsilon(t, xi3, 0.38095238095238093)
	assertEpsilon(t, pval3, 0.04989192742513937)

	xi4, pval4, _ := New(
		anscombesQuartet["x_4"],
		anscombesQuartet["y_4"],
		WithAsymptoticPvalue(),
	).Pvalue()
	assertEpsilon(t, xi4, 0.19999999999999996)
	assertEpsilon(t, pval4, 0.1515801165640982)

	// Permutation p-value calculation
	xi1, pval1, _ = New(
		anscombesQuartet["x_1"],
		anscombesQuartet["y_1"],
		WithPermutationPvalue(1000),
	).Pvalue()
	assertEpsilon(t, xi1, 0.2749999999999999)

	wantPval := 0.059
	if abs(pval1-wantPval) > 0.02 {
		t.Errorf("failed float assertion: got:%v, want:%v", pval1, wantPval)
	}

	// Without ties
	xi1, pval1, _ = New(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithoutTies(),
	).Pvalue()
	assertEpsilon(t, xi1, 0.7272727)
	assertEpsilon(t, pval1, 0.000138257)

	// With asymptotic pvals
	xi1, pval1, _ = New(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithAsymptoticPvalue(),
	).Pvalue()
	assertEpsilon(t, xi1, 0.7272727)
	assertEpsilon(t, pval1, 0.0001879616)

	// With permutations
	xi1, pval1, _ = New(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithPermutationPvalue(1000),
	).Pvalue()
	assertEpsilon(t, xi1, 0.7272727)
	assertEpsilon(t, pval1, 0)

	wantPval = 0.051425
	xi1, pval1, _ = New(
		anscombesQuartet["x_1"],
		anscombesQuartet["y_1"],
		WithPermutationPvalue(200_000),
	).Pvalue()
	assertEpsilon(t, xi1, 0.275)
	if abs(pval1-wantPval) > 0.01 {
		t.Errorf("failed float assertion: got:%v, want:%v", pval1, wantPval)
	}

	// 1000-length array test
	xi1, pval1, _ = New(xx, yy, WithAsymptoticPvalue()).Pvalue()
	assertEpsilon(t, xi1, 0.001335001) // R value 0.001335001
	assertEpsilon(t, pval1, 0.4733904) // R value 0.4733904
}

func TestXiErrors(t *testing.T) {
	x := []float64{1, 2, 3}
	y := []float64{5, 6}

	xi := New(x, y)
	_, err := xi.Correlation()
	if err.Error() != "xicor: mismatched size of input vectors" {
		t.Errorf("didn't receive the correct error when providing input of different lengths: %v", err)
	}

	_, _, err = xi.Pvalue()
	if err.Error() != "xicor: mismatched size of input vectors" {
		t.Errorf("didn't receive the correct error when providing input of different lengths: %v", err)
	}

	xi = New(x, x)
	xi.Method = "invalid method"
	_, _, err = xi.Pvalue()
	if err.Error() != "xicor: invalid p-value calculation method; use either 'asymptotic' or 'permutation'" {
		t.Errorf("didn't receive the correct error when providing an invalid p-value calculation method: %v", err)
	}

	xi.Method = MethodAsymptotic
	xi.WantPvalue = false

	_, _, err = xi.Pvalue()
	if err.Error() != "xicor: trying to calculate the p-value on an object where `Xi.WantPvalues=false`" {
		t.Errorf("didn't receive the correct error when asking for p-value with WantPvalue=false: %v", err)
	}
}

// Test helpers

func TestRemoveNaNs(t *testing.T) {
	a := []float64{0, 1, math.NaN(), 3, 4, math.NaN(), 6}
	b := []float64{8, math.NaN(), 6, 5, 4, math.NaN(), 2}

	gotA, gotB := removeNaNs(a, b)
	if len(gotA) != len(gotB) {
		t.Error("removeNaNs returned slices of incorrect lengths")
	}

	wantA := []float64{0, 3, 4, 6}
	if !reflect.DeepEqual(gotA, wantA) {
		t.Errorf("removeNaNs did not correctly trim slice A ; got: %v, want: %v", gotA, wantA)
	}

	wantB := []float64{8, 5, 4, 2}
	if !reflect.DeepEqual(gotB, wantB) {
		t.Errorf("removeNaNs did not correctly trim slice B; got: %v, want: %v", gotB, wantB)
		fmt.Printf("%#v\n%T\n", b, b)
		fmt.Printf("%#v\n%T\n", wantB, wantB)
	}
}

func TestRank_RND_Max(t *testing.T) {
	rand.Seed(time.Now().Unix())
	input := []float64{0., 2., 3., 2.}

	want1 := []float64{1., 2., 4., 3.}
	want2 := []float64{1., 3., 4., 2.}

	got := rankRND(input)
	if !reflect.DeepEqual(got, want1) &&
		!reflect.DeepEqual(got, want2) {
		t.Errorf("wrong result for rankRND with input:%v, got:%v", input, got)
	}

	got = rankMax(input)
	want := []float64{1., 3., 4., 3.}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("wrong result for rankMax with input:%v, got:%v, want:%v", input, got, want)
	}
}

func TestArgsort(t *testing.T) {
	input := []float64{3., 1., 2.}

	got := argsort(input)
	want := []int{1, 2, 0}

	if !reflect.DeepEqual(got, want) {
		t.Errorf("wrong result for argsort with input:%v, got:%v, want:%v", input, got, want)
	}
}

func TestArgsortWithTies(t *testing.T) {
	rand.Seed(time.Now().Unix())
	input := []float64{5, 10, 2, 99, 5, 2, 8, 17, 5}
	got := argsort(input)

	// [2 5 0 8 4 6 1 7 3]
	// [5 2 0 8 4 6 1 7 3]
	// [2 5 4 0 8 6 1 7 3]
	// [2 5 0 8 4 6 1 7 3]
	// [2 5 4 8 0 6 1 7 3]
	// [5 2 8 0 4 6 1 7 3]
	// ties (i.e. indexes {2, 5} and {0, 4, 8}) should be shuffled

	tie1 := got[0:2]
	tie2 := got[2:5]
	rest := got[5:]
	if !(contains(tie1, 2) && contains(tie1, 5)) {
		t.Error("first tie block should contain 2 and 5")
	}

	if !(contains(tie2, 0) && contains(tie2, 4) && contains(tie2, 8)) {
		t.Error("second tie should contain 0, 4 and 8")
	}

	wantRest := []int{6, 1, 7, 3}
	if !reflect.DeepEqual(rest, wantRest) {
		t.Error("rest of the argsort result should be [1 7 3]")
	}
}

func TestAbs(t *testing.T) {
	val := 0.
	if abs(val) != abs(-val) {
		t.Errorf("wrong result for absolute value with input %v", val)
	}

	val = 101.
	if abs(val) != abs(-val) {
		t.Errorf("wrong result for absolute value with input %v", val)
	}

	val = -42.
	if abs(val) != abs(-val) {
		t.Errorf("wrong result for absolute value with input %v", val)
	}
}

func TestMinMax(t *testing.T) {
	input := []int{5, 5, 5}
	gotMin, gotMax := min(input), max(input)
	wantMin, wantMax := 5, 5

	if gotMin != wantMin {
		t.Errorf("incorrect minimum value for input: %v, got: %v, want: %v", input, gotMin, wantMin)
		t.Errorf("incorrect maximum value for input: %v, got: %v, want: %v", input, gotMax, wantMax)
	}

	input = []int{-4, 0, 8, 16, -99, 4, 7}
	gotMin, gotMax = min(input), max(input)
	wantMin, wantMax = -99, 16

	if gotMin != wantMin {
		t.Errorf("incorrect minimum value for input: %v, got: %v, want: %v", input, gotMin, wantMin)
		t.Errorf("incorrect maximum value for input: %v, got: %v, want: %v", input, gotMax, wantMax)
	}
}

func TestMean(t *testing.T) {
	input := []float64{}
	got := mean(input)
	if !math.IsNaN(got) {
		t.Errorf("expected NaN result for the mean of an empty list, got %v", got)
	}

	input = []float64{1, 1, 2, 3, 5, 8, 13, 21}
	got = mean(input)
	want := 6.75
	if got != want {
		t.Errorf("wrong result for mean with input: %v, got:%v, want:%v", input, got, want)
	}
}

func TestCumsum(t *testing.T) {
	input := []float64{}
	got := cumsum(input)
	want := []float64{}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("wrong result for cumsum with input: %v, got:%v, want:%v", input, got, want)
	}

	input = []float64{1, 2, 3, 4, 5, 6}
	got = cumsum(input)
	want = []float64{1, 3, 6, 10, 15, 21}

	if !reflect.DeepEqual(got, want) {
		t.Errorf("wrong result for cumsum with input: %v, got:%v, want:%v", input, got, want)
	}
}

func TestSum(t *testing.T) {
	input := []float64{}
	got := sum(input)
	if got != 0 {
		t.Errorf("sum produced an incorrect result for input: %v", input)
	}

	input = []float64{1, 2, 3, 4, 5, 6}
	got = sum(input)
	if got != 21 {
		t.Errorf("sum produced an incorrect result for input: %v", input)
	}

	input = []float64{1, 2, 3, math.NaN()}
	got = sum(input)
	if !math.IsNaN(got) {
		t.Error("expected NaN result when summing with a NaN", input)
	}
}

// Benchmarks

var corr, pvalue float64

func BenchmarkCorrelation(b *testing.B) {
	b.ReportAllocs()
	x := make([]float64, 100_000)
	y := make([]float64, 100_000)

	for i := range x {
		x[i] = rand.NormFloat64()
		y[i] = rand.NormFloat64()
	}

	var c float64
	xi := New(x, y)
	for i := 0; i < b.N; i++ {
		c, _ = xi.Correlation()
	}
	corr = c
}

func BenchmarkPvalues(b *testing.B) {
	b.ReportAllocs()
	x := make([]float64, 100_000)
	y := make([]float64, 100_000)

	for i := range x {
		x[i] = rand.NormFloat64()
		y[i] = rand.NormFloat64()
	}

	var c, p float64
	xi := New(x, y)
	for i := 0; i < b.N; i++ {
		c, p, _ = xi.Pvalue()
	}
	corr = c
	pvalue = p
}

// Test helpers

func contains(s []int, e int) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func assertEpsilon(t *testing.T, got, want float64) {
	epsilon := 0.00001
	if abs(got-want) > epsilon {
		t.Errorf("failed float assertion: got:%v, want:%v", got, want)
	}
}
