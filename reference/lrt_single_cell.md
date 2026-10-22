# Likelihood-ratio test for one single-cell composition

Compares the profiled maximum with the maximum on `fixed` (passed as
`null_value`). An interior null uses Wilks' chi-square. A null on a face
uses the chi-bar-square mixture of
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md).
The two technical scales are profiled out under both hypotheses, so the
degrees of freedom count composition constraints only.

## Usage

``` r
lrt_single_cell(
  y,
  mean_signature_matrix,
  Sigma,
  null_value,
  n_boundary = NULL,
  boundary_tol = 1e-06,
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

- null_value:

  Named hypothesised ratios.

- n_boundary, boundary_tol:

  As in
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md).

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of Newton iterations.

## Value

A one-row data frame, as for
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md).

## See also

[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md),
[`fit_firth_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_single_cell.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
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
y <- 20 * drop(mu %*% c(1, 0)) + c(0.2, -0.1, 0.05, 0.05)
lrt_single_cell(y, mu, sigma, null_value = c(ct2 = 0))
#>     statistic n_boundary df_interior p_value loglik_full loglik_null
#> ct2         0          1           0     0.5    8.127704    8.127704
#>                     calibration
#> ct2 chi-bar-square (Self-Liang)
```
