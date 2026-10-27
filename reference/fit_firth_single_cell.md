# Firth fit of one single-cell bulk sample

Maximises the profile log-likelihood plus half the log determinant of
the observed profile information. The penalty uses the negative profile
Hessian, not an expected Fisher matrix. Scales stay concentrated out, so
the penalty does not add a third free scale.

## Usage

``` r
fit_firth_single_cell(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 1e-08,
  itmax = 80L
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of iterations.

## Value

A list with `coefficients`, `rho`, `loglik`, `penalty`,
`penalised_loglik`, `mean_scale`, `variance_scale` and `convergence`.

## See also

[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md)

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
fit_firth_single_cell(y, mu, sigma)$coefficients
#>       ct1       ct2 
#> 0.7053266 0.2946734 
```
