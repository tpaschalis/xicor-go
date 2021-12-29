package xicor

import (
	"math/rand"
	"reflect"
	"testing"
	"time"
)

var anscombes_quartet = map[string][]float64{
	"x_1": {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5},
	"x_2": {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5},
	"x_3": {0, 0, 0, 0, 1, 1, 2, 3},
	"x_4": {8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8},
	"y_1": {8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68},
	"y_2": {9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74},
	"y_3": {0, 1, 2, 3, 4, 5, 6, 7},
	"y_4": {6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.5, 5.56, 7.91, 6.89},
}

func TestRank(t *testing.T) {
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

func TestCumsum(t *testing.T) {
	input := []float64{1, 2, 3, 4, 5, 6}
	got := cumsum(input)
	want := []float64{1, 3, 6, 10, 15, 21}

	if !reflect.DeepEqual(got, want) {
		t.Errorf("wrong result for cumsum with input: %v, got:%v, want:%v", input, got, want)
	}
}

func TestXi(t *testing.T) {
	rand.Seed(21)
	//rand.Seed(time.Now().Unix())

	// Asymptotic p-value calculation
	xi_1, pval_1, _ := NewFloat64(
		anscombes_quartet["x_1"],
		anscombes_quartet["y_1"],
		WithAsymptoticPvals(),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.275)
	assertEpsilon(t, pval_1, 0.07841556)
	// TODO this value isn't correct
	//assertEpsilon(t, pval_1, 0.059)

	xi_2, pval_2, _ := NewFloat64(
		anscombes_quartet["x_2"],
		anscombes_quartet["y_2"],
		WithAsymptoticPvals(),
	).Pvalues()
	assertEpsilon(t, xi_2, 0.6)
	assertEpsilon(t, pval_2, 0.0010040217037570187)

	xi_3, pval_3, _ := NewFloat64(
		anscombes_quartet["x_3"],
		anscombes_quartet["y_3"],
		WithAsymptoticPvals(),
	).Pvalues()
	assertEpsilon(t, xi_3, 0.38095238095238093)
	assertEpsilon(t, pval_3, 0.04989192742513937)

	xi_4, pval_4, _ := NewFloat64(
		anscombes_quartet["x_4"],
		anscombes_quartet["y_4"],
		WithAsymptoticPvals(),
	).Pvalues()
	assertEpsilon(t, xi_4, 0.19999999999999996)
	assertEpsilon(t, pval_4, 0.1515801165640982)

	// Permutation p-value calculation
	xi_1, pval_1, _ = NewFloat64(
		anscombes_quartet["x_1"],
		anscombes_quartet["y_1"],
		WithPermutationPvals(1000),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.2749999999999999)

	wantPval := 0.059
	if abs(pval_1-wantPval) > 0.02 {
		t.Errorf("failed float assertion: got:%v, want:%v", pval_1, wantPval)
	}
}

func TestCompareWithR(t *testing.T) {

	// Without ties
	xi_1, pval_1, _ := NewFloat64(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithoutTies(),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.7272727)
	assertEpsilon(t, pval_1, 0.000138257)

	// With asymptotic pvals
	xi_1, pval_1, _ = NewFloat64(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithAsymptoticPvals(),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.7272727)
	assertEpsilon(t, pval_1, 0.0001879616)

	// With permutations
	xi_1, pval_1, _ = NewFloat64(
		[]float64{1., 2., 3., 4., 5., 6., 7., 8., 9., 10.},
		[]float64{5., 6., 7., 8., 9., 10., 11., 12., 13., 14.},
		WithPermutationPvals(1000),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.7272727)
	assertEpsilon(t, pval_1, 0)

	wantPval := 0.051425
	xi_1, pval_1, _ = NewFloat64(
		anscombes_quartet["x_1"],
		anscombes_quartet["y_1"],
		WithPermutationPvals(200_000),
	).Pvalues()
	assertEpsilon(t, xi_1, 0.275)
	if abs(pval_1-wantPval) > 0.01 {
		t.Errorf("failed float assertion: got:%v, want:%v", pval_1, wantPval)
	}
}

func assertEpsilon(t *testing.T, got, want float64) {
	epsilon := 0.00001
	if abs(got-want) > epsilon {
		t.Errorf("failed float assertion: got:%v, want:%v", got, want)
	}
}
