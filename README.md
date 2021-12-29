# xicor-go

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## What is xicor?
`xicor-go` is a Go implementation of the "xi" correlation metric described in the paper [_"A new coefficient of correlation"_](https://arxiv.org/pdf/1909.10140.pdf) by Sourav Chatterjee (arxiv.org/abs/1909.10140).

The current implementation is based off the [R code](https://statweb.stanford.edu/~souravc/xi.R) mentioned in the paper, and influenced by a [Python port](https://github.com/czbiohub/xicor) by user [czbiohub](https://github.com/czbiohub/xicor). There is also a [mirror](https://github.com/cran/XICOR) of the CRAN R package hosted on GitHub.

The current package contains an implementation which uses the `Xi` struct and utilizes the [Go generics](https://go.dev/doc/tutorial/generics) introduced in Go 1.18, and a `XiFloat64` struct which works only for 64-bit slices. This distinction was made to gauge the performance difference between the two, as well as introduce float-specific performance improvements in the future!

## Installation
For a project utilizing Go modules, all you need to do is
```
go get github.com/tpaschalis/xicor-go
```

## Usage

```go
import "github.com/tpaschalis/xicor-go"

// Create a new Xi object and get the correlation
xi, err := xicor.NewFloat64(X, Y).Correlation()

// Create a new Xi object and obtain the correlation along with its p-value
// Use functional options to define
xi, pvalue, err := xicor.NewFloat64(
	X,
	Y,
	WithPermutationPvals(1000),
	WithoutTies(),
).Pvalues()

// You can also use the Xi object directly
data := &XiFloat64{
	X:         x,
	Y:         y,
	WantPvals: true,
	Nperms:    1000,
	Method:    "asymptotic",
	DataTies:  true,
}

xi, pvalue, err = data.Pvalues()
```

## Benchmarks
Coming soon!

## Contributing
If you have an idea, or would like to discuss and contribute an improvement, you can reach out in the repo [Issues](https://github.com/tpaschalis/xicor-go/issues) and open a [Pull Request](https://github.com/tpaschalis/xicor-go/pulls)

## License
This code present in this repo is licensed under the [MIT license](LICENSE).
