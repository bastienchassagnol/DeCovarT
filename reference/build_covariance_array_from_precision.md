# Build cell-type-specific covariances from precision slices

For each cell type \\j\\, inverts \\\boldsymbol{\Omega}\_j\\ to
\\\boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}\_j^{-1}\\ (via
[`qr.solve()`](https://rdrr.io/r/base/qr.html)) and stacks the result as
a \\G\times G\times J\\ array. If a slice is not numerically SPD, or if
the inverse fails Cholesky, a further uniform spectral shift is applied
to that precision (off-diagonal support unchanged). Precisions need not
be shared across cell types.

## Usage

``` r
build_covariance_array_from_precision(precision_array)
```

## Arguments

- precision_array:

  Symmetric positive-definite array
  \\(\boldsymbol{\Omega}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\.

## Value

Numeric array of dimension \\G\times G\times J\\.

## Examples

``` r
A <- generate_random_network_skeleton(6L, "hub")
W <- assign_iid_signed_weights(A)
Omega <- build_normalised_precision(W, precision_shift = 0.2)
arr <- array(Omega, dim = c(6, 6, 1))
dim(build_covariance_array_from_precision(arr))
#> [1] 6 6 1
```
