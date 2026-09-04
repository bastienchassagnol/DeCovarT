# Log-determinant of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Log-determinant of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

## Usage

``` r
sigma_logdet(backend, p)
```

## Arguments

- backend:

  A
  [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  object.

- p:

  Numeric proportions vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

## Value

Scalar \\\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\.

## See also

[`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md),
[`sigma_quadform()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_quadform.md)

## Examples

``` r
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
sigma_logdet(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4))
#> [1] -1.307853
```
