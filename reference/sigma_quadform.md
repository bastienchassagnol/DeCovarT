# Mahalanobis quadratic form \\\mathbf{r}^{\mathsf{T}}\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{r}\\

Mahalanobis quadratic form
\\\mathbf{r}^{\mathsf{T}}\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{r}\\

## Usage

``` r
sigma_quadform(backend, p, r)
```

## Arguments

- backend:

  A
  [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  object.

- p:

  Numeric proportions vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- r:

  Numeric residual vector.

## Value

Scalar quadratic form.

## See also

[`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md)

## Examples

``` r
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
sigma_quadform(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4), c(1, 2))
#> [1] 9.615385
```
