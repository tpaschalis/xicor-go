package xicor

import (
	"errors"
	"math"
	"math/rand"
	"sort"
)

// Xi is the main struct exposed by the `xicor` package.
// It is used to receive the input data as well as configure some options on how to calculate the correlation coefficient and p-values.
// In general, it is safe to just supply x and y, and leave other parameters to their default values.
type Xi struct {
	X, Y       []float64
	WantPvalue bool
	Nperms     int
	Method     string
	DataTies   bool

	// variables reused for p-values calculation
	n    float64
	f    []float64
	cval float64
}

// MethodAsymptotic employs the 'asympotic theory' to calculate the p-value.
var MethodAsymptotic = "asymptotic"

// MethodPermutation employes `NPerms` permutations to estimate the p-value. As per Dr. Sourav, "usually there is no need for the permutation test, the asymptotic theory is good enough".
var MethodPermutation = "permutation"

// New creates a `Xi` object which can be used to calculate the correlation coefficient along with the p-value. It receives the input datasets, as well as a number of functional options to configure the runtime behavior.
func New(x, y []float64, options ...func(*Xi)) *Xi {
	res := &Xi{
		X:          x,
		Y:          y,
		WantPvalue: true,
		Nperms:     1000,
		Method:     "asymptotic",
		DataTies:   true,
	}

	for _, o := range options {
		o(res)
	}

	return res
}

// WithAsymptoticPvalue makes sure that the asymptotic theory will be used to calculate the p-value.
func WithAsymptoticPvalue() func(*Xi) {
	return func(d *Xi) {
		d.WantPvalue = true
		d.Method = MethodAsymptotic
	}
}

// WithPermutationPvalue makes sure that the p-value will be estimated using `Nperms` permutations.
func WithPermutationPvalue(nperms int) func(*Xi) {
	return func(d *Xi) {
		d.WantPvalue = true
		d.Method = MethodPermutation
		d.Nperms = nperms
	}
}

// WithoutTies informs the algorithm that there are no ties in the data, and uses some simpler theory to calculate the p-value. There is no harm in leaving DataTies to `true` even if there are no ties.
func WithoutTies() func(*Xi) {
	return func(d *Xi) {
		d.WantPvalue = true
		d.DataTies = false
	}
}

// Correlation calculates and returns the correlation coefficient for the input data vectors `X` and `Y` along with an error.
func (d *Xi) Correlation() (float64, error) {
	if len(d.X) != len(d.Y) {
		return 0, errors.New("xicor: mismatched size of input vectors")
	}

	// x, y are the data vectors
	// Find and Remove N/A pair values
	removeNaNs(d.X, d.Y)

	// Factor variables should be converted to integers here
	// https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/factor
	//	if (factor == T) {
	//	if (!is.numeric(x)) x = as.numeric(factor(x))
	//	if (!is.numeric(y)) y = as.numeric(factor(y))
	//}

	// Sample Size
	d.n = float64(len(d.X))

	// PI is the rank vector for x, with ties broken at random
	pi := rankRND(d.X)

	// f[i] is number of j s.t. y[j] <= y[i], divided by n.
	d.f = rankMax(d.Y)
	for i := range d.f {
		d.f[i] = d.f[i] / d.n
	}

	// g[i] is number of j s.t. y[j] >= y[i], divided by n.
	ym := make([]float64, len(d.Y))
	for i := range d.Y {
		ym[i] = -d.Y[i]
	}
	g := rankMax(ym)
	for i := range g {
		g[i] = g[i] / d.n
	}

	// order of the x's, ties broken at random.
	ord := argsort(pi)

	// Rearrange f according to ord.
	ford := make([]float64, len(d.f))
	for i, val := range ord {
		_ = val
		ford[i] = d.f[ord[i]]
	}

	// xi is calculated in the next lines
	diffs := make([]float64, 0)
	for i := 0; i < int(d.n)-1; i++ {
		diffs = append(diffs, abs(ford[i]-ford[i+1]))
	}
	A1 := mean(diffs) * (d.n - 1) / (2 * d.n)

	muls := make([]float64, 0)
	for _, val := range g {
		muls = append(muls, val*(1-val))
	}
	d.cval = mean(muls)

	xi := 1 - A1/d.cval

	return xi, nil
}

// Pvalue calculates and returns the correlation coefficient and p-value for the input data vectors `X` and `Y` along with an error.
func (d *Xi) Pvalue() (float64, float64, error) {
	xi, err := d.Correlation()
	if err != nil {
		return 0, 0, err
	}
	var pval float64

	if d.Method != "asymptotic" && d.Method != "permutation" {
		return 0, 0, errors.New("xicor: invalid p-value calculation method; use either 'asymptotic' or 'permutation'")
	}
	if !d.WantPvalue {
		return 0, 0, errors.New("xicor: trying to calculate the p-value on an object where `Xi.WantPvalues=false`")
	}

	// If there are no data ties, we can use some simpler theory to calculate the theoretical P-value
	if d.DataTies == false {
		pval = 1 - pnorm(math.Sqrt(d.n)*xi/math.Sqrt(2./5.))
		return xi, pval, nil
	}

	// If there are ties in the input data, the algorithm employs the more elaborated theory for calculating the P-value
	// There is no harm in setting DataTies to true and using the fancy P-value calculation even if there are no ties
	if d.Method == "asymptotic" {
		q := make([]float64, len(d.f))
		copy(q, d.f)
		sort.Float64s(q)

		ind := make([]float64, int(d.n))
		ind2 := make([]float64, int(d.n))

		for i := 0; i < int(d.n); i++ {
			ind[i] = float64(i) + 1.
			ind2[i] = 2*d.n - 2*ind[i] + 1
		}

		asl := make([]float64, int(d.n))
		csl := make([]float64, int(d.n))

		for i := 0; i < int(d.n); i++ {
			asl[i] = ind2[i] * q[i] * q[i]
			csl[i] = ind2[i] * q[i]
		}
		a := mean(asl) / d.n
		c := mean(csl) / d.n
		cq := cumsum(q)

		m := make([]float64, int(d.n))
		msq := make([]float64, int(d.n))
		for i := 0; i < int(d.n); i++ {
			m[i] = (cq[i] + (d.n-ind[i])*q[i]) / d.n
			msq[i] = m[i] * m[i]
		}
		b := mean(msq)

		v := (a - 2*b + c*c) / (d.cval * d.cval)

		pval = 1 - pnorm(math.Sqrt(d.n)*xi/math.Sqrt(v))
		return xi, pval, nil
	}

	// If permutation test is to be used for calculating P-value:
	if d.Method == "permutation" {
		r := make([]float64, d.Nperms)
		for i := 0; i < d.Nperms; i++ {
			x1 := make([]float64, int(d.n))
			for i := 0; i < int(d.n); i++ {
				x1[i] = rand.Float64()
			}
			xinew, _, _ := New(x1, d.Y, WithAsymptoticPvalue()).Pvalue()
			r[i] = xinew
		}
		ps := make([]float64, d.Nperms)
		for i := 0; i < d.Nperms; i++ {
			if r[i] > xi {
				ps[i] = 1.0
			} else {
				ps[i] = 0.0
			}
		}
		pval = mean(ps)
	}

	return xi, pval, nil
}

func removeIdx(a []int, i int) []int {
	return append(a[:i], a[i+1:]...)
}

func removeNaNs(x, y []float64) ([]float64, []float64) {
	nans := make(map[int]struct{})
	for i, xv := range x {
		if math.IsNaN(xv) {
			nans[i] = struct{}{}
		}
	}
	for j, yv := range y {
		if math.IsNaN(yv) {
			nans[j] = struct{}{}
		}
	}

	var newX, newY []float64
	for i := 0; i < len(x); i++ {
		if _, ok := nans[i]; ok {
			continue
		}
		newX = append(newX, x[i])
		newY = append(newY, y[i])
	}

	return newX, newY
}

func rankRND(a []float64) []float64 {
	ns := make([]float64, len(a))
	copy(ns, a)
	sort.Float64s(ns)

	idx := make(map[float64][]int)

	for i, val := range ns {
		idx[val] = append(idx[val], i+1)
	}

	res := make([]float64, len(a))
	for i := range res {
		pool := idx[a[i]]
		selectedIdx := rand.Intn(len(pool))
		res[i] = float64(pool[selectedIdx])
		idx[a[i]] = removeIdx(idx[a[i]], selectedIdx)
	}

	return res
}

func rankMax(a []float64) []float64 {
	ns := make([]float64, len(a))
	copy(ns, a)
	sort.Float64s(ns)

	idx := make(map[float64][]int)
	for i, val := range ns {
		idx[val] = append(idx[val], i+1)
	}

	res := make([]float64, len(a))
	for i := range res {
		pool := idx[a[i]]
		res[i] = float64(pool[len(pool)-1])
	}

	return res
}

type fsl struct {
	sort.Float64Slice
	idx []int
}

func (s *fsl) Swap(i, j int) {
	s.Float64Slice.Swap(i, j)
	s.idx[i], s.idx[j] = s.idx[j], s.idx[i]
}

// Initial approach from https://stackoverflow.com/a/31141540
func argsort(a []float64) []int {

	indexes := make([]int, len(a))
	for i := range indexes {
		indexes[i] = i
	}

	s := &fsl{Float64Slice: a, idx: indexes}
	sort.Sort(s)

	// Lets identify ties, which we should randomize
	var currentStreak []int
	for i := 0; i < len(a)-1; i++ {
		currentStreak = append(currentStreak, i)
		if a[i] == a[i+1] {
			// we're on a streak
			continue
		} else {
			// a streak was broken
			if len(currentStreak) > 1 {
				// we need to randomize the indices
				first := currentStreak[0]
				last := currentStreak[len(currentStreak)-1]
				shuffle(indexes[first : last+1])
			}
			// reset the streak
			currentStreak = []int{}
		}
	}

	return indexes
}

func shuffle(a []int) {
	rand.Shuffle(len(a), func(i, j int) { a[i], a[j] = a[j], a[i] })
	return
}

func abs(a float64) float64 {
	if a >= 0 {
		return a
	}
	return -a
}

func min(a []int) int {
	res := a[0]
	for _, val := range a {
		if val < res {
			res = val
		}
	}
	return res
}

func max(a []int) int {
	res := a[0]
	for _, val := range a {
		if val > res {
			res = val
		}
	}
	return res
}

func mean(a []float64) float64 {
	sum := 0.
	for _, val := range a {
		sum += val
	}
	return sum / float64(len(a))
}

func cumsum(a []float64) []float64 {
	res := make([]float64, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = sum(a[:i+1])
	}

	return res
}

func sum(a []float64) float64 {
	var sum float64
	for _, val := range a {
		sum += val
	}

	return sum
}

// From https://en.wikipedia.org/wiki/Normal_distribution
func pnorm(a float64) float64 {
	var val, sum, res float64
	sum = a
	val = a
	for i := 1; i <= 100; i++ {
		val = (val * a * a / (2*float64(i) + 1))
		sum = sum + val
	}
	res = 0.5 + (sum/math.Sqrt(2*math.Pi))*math.Exp(-(a*a)/2)
	return res
}
