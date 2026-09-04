# Precision solve \\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{B}\\

Precision solve \\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{B}\\

## Usage

``` r
sigma_solve(backend, p, b)
```

## Arguments

- backend:

  A
  [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  object.

- p:

  Numeric proportions vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- b:

  Numeric vector or matrix of right-hand sides.

## Value

\\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{b}\\, shaped like `b`.

## See also

[`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md)

## Examples

``` r
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
sigma_solve(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4), c(1, 2))
#> [1] 1.923077 3.846154
```
