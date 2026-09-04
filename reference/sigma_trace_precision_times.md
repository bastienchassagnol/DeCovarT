# Trace of the precision times a matrix, \\\operatorname{tr}(\boldsymbol{\Theta}(\boldsymbol{p})\mathbf{S})\\

Computes \\\operatorname{tr}(\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
\mathbf{S})\\ through solves only, never materialising the dense
precision. This is the trace term that appears in the analytic gradient
and Hessian.

## Usage

``` r
sigma_trace_precision_times(backend, p, S)
```

## Arguments

- backend:

  A
  [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  object.

- p:

  Numeric proportions vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- S:

  Numeric \\G\times G\\ matrix (typically a cell-type covariance slice
  \\\boldsymbol{\Sigma}\_j\\).

## Value

Scalar trace.

## See also

[`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md)

## Examples

``` r
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
cov <- new_decovart_covariance(Sigma, "dense")
sigma_trace_precision_times(cov, c(0.6, 0.4), Sigma[, , 1])
#> [1] 3.846154
```
