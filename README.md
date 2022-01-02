# xicor-go

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## What is xicor?
`xicor-go` is a Go implementation of the "xi" correlation metric described in the paper [_"A new coefficient of correlation"_](https://arxiv.org/pdf/1909.10140.pdf) by Sourav Chatterjee (arxiv.org/abs/1909.10140).

The current implementation is based off the [R code](https://statweb.stanford.edu/~souravc/xi.R) mentioned in the paper, and influenced by a [Python port](https://github.com/czbiohub/xicor) by user [czbiohub](https://github.com/czbiohub/xicor). There is also a [mirror](https://github.com/cran/XICOR) of the CRAN R package hosted on GitHub.

The package provides the `Xi` struct which receives the input data along with a number of functional options and is used to calculate the correlation coefficient and p-value.

<!--
There is also a dev branch which utilizes [Go generics](https://go.dev/doc/tutorial/generics) introduced in Go 1.18, which will be merged into the canonical release once Go 1.18 is released in February.
-->

## Installation
For a project utilizing Go modules, all you need to do is
```
go get github.com/tpaschalis/xicor-go
```

## Usage

```go
import "github.com/tpaschalis/xicor-go"

func main() {
	x := []float64{0, 1, 2, 3, 4}
	y := []float64{5, 6, 7, 8, 9}

	// Create a new Xi object and get the correlation
	xi, err := xicor.New(x, y).Correlation()

	// Create a new Xi object and obtain the correlation along with its p-value
	// Use functional options to define how to perform the calculation
	xi, pvalue, err := xicor.New(
		x,
		y,
		xicor.WithPermutationPvalue(1000),
		xicor.WithoutTies(),
	).Pvalue()

	// You can also use the Xi object directly
	data := &xicor.Xi{
		X:          x,
		Y:          y,
		WantPvalue: true,
		Method:     "asymptotic",
		DataTies:   true,
	}

	xi, pvalue, err = data.Pvalue()
	fmt.Println(xi, pvalue, err)
}
```

## Current status
I'm working towards a more stable and performant v0.0.1 release; the focus is on:
- Validating correctness of results by comparing against original R code (current tests haven't produced any inconsistency yet)
- Add support for categorical variables (a la the `factor` function used in R)
- Use dev branch to rewrite some of the code using Go 1.18 generics
- Run through a profiler to find and eliminate bottlenecks

## Benchmarks
The current pre-v0.0.1 version of the code runs through 100k randomized samples in about 0.12 seconds.

While there's a lot of low-hanging fruit to improve performance, I'll get to it once v0.0.1 is near.

## Contributing
If you have an idea, or would like to discuss and contribute an improvement, you can reach out in the repo [Issues](https://github.com/tpaschalis/xicor-go/issues) and open a [Pull Request](https://github.com/tpaschalis/xicor-go/pulls)

## License
This code in this repo is licensed under the [MIT license](LICENSE).
