# Observed information of one single-cell fit

At a fitted composition, the covariance of the ILR coordinate
\\\hat{\boldsymbol{\rho}}\\ is the inverse of the negative profile
Hessian
[`hessian_profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md).
The delta method carries it to \\\hat{\boldsymbol{p}}\\ with
\\\mathbf{J}\_{\psi}\\. The two technical scales use the other Schur
complement of the joint observed Hessian in
\\(\boldsymbol{\rho},c,\phi)\\: estimating the composition widens their
variances, and \\I\_{c\phi}\\ is not assumed to make that complement
diagonal.

Both matrices are observed information, not expected Fisher information.
They treat the reference moments as known.

## Usage

``` r
vcov_single_cell(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale = NULL,
  variance_scale = NULL,
  scales_known = FALSE,
  warn = TRUE
)
```

## Arguments

- p:

  Fitted composition. A boundary coordinate makes the ILR chart
  undefined.

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- mean_scale, variance_scale:

  Scales at `p`. `NULL` recomputes the closed forms.

- scales_known:

  If `TRUE`, the scales are treated as known: the composition covariance
  inverts the fixed-scale observed Hessian
  [`hessian_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_single_cell_ilr.md)
  at `mean_scale` and `variance_scale`, and `vcov_nuisance` is `NA`.
  This is the covariance of a `"fixed_scales"` fit. For a `"two_stage"`
  fit it ignores the uncertainty of the start scales, so its intervals
  are too narrow.

- warn:

  If `TRUE`, warn when either information matrix is not positive
  definite.

## Value

A list with `rho`, `vcov_rho`, `vcov` (delta method for
\\\boldsymbol{p}\\), `mean_scale`, `variance_scale` and `vcov_nuisance`
(order \\c\\, \\\phi\\). Singular blocks are filled with `NA`.

## See also

[`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)

## Examples

``` r
mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
vcov_single_cell(c(0.7, 0.3), y, mu, sigma)
#> Warning: The profile information is not positive definite.
#> Warning: The nuisance information is not positive definite.
#> $rho
#> [1] 0.5991301
#> 
#> $vcov_rho
#>      [,1]
#> [1,]   NA
#> 
#> $vcov
#>      [,1] [,2]
#> [1,]   NA   NA
#> [2,]   NA   NA
#> 
#> $mean_scale
#> [1] 20.00234
#> 
#> $variance_scale
#> [1] 0.7497659
#> 
#> $vcov_nuisance
#>                mean_scale variance_scale
#> mean_scale             NA             NA
#> variance_scale         NA             NA
#> 
```
