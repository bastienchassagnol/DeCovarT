# Profile-likelihood intervals for one single-cell bulk sample

Inverts the likelihood-ratio test for one coordinate of the profile
log-likelihood. The two scales are concentrated out at every trial
value. An endpoint that lands on a face is not a Wilks interval; use
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)
for a boundary null.

## Usage

``` r
confint_profile_single_cell(
  y,
  mean_signature_matrix,
  Sigma,
  level = 0.95,
  celltypes = colnames(mean_signature_matrix),
  variance = c("ml", "reml"),
  epsilon = 1e-08,
  itmax = 100L
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- level:

  Confidence level.

- celltypes:

  Cell types to profile. Defaults to every column.

- variance:

  `"ml"` (default) uses the \\\chi^2_1\\ cut-off. `"reml"` uses
  \\G\log\\1+F\_{1,\nu;1-\alpha}/\nu\\\\ with \\\nu=G-J\\, the exact
  cut-off of the profile likelihood ratio when a variance scale is
  profiled out of a linear model.

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of Newton iterations.

## Value

A matrix with columns `estimate`, `lower` and `upper`.

## See also

[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md),
[`fit_firth_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_single_cell.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md),
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)

## Examples

``` r
mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
colnames(mu) <- c("ct1", "ct2")
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
confint_profile_single_cell(y, mu, sigma)
#>     estimate     lower     upper
#> ct1 0.705328 0.7025836 0.7080680
#> ct2 0.294672 0.2919320 0.2974164
```
