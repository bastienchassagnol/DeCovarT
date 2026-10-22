# Confidence intervals by mapping a Wald region through the ILR chart

Forms a Wald interval for the unconstrained ILR coordinate
\\\boldsymbol{\rho}\\ from the expected Fisher information, then maps it
to the simplex by the inverse chart \\\boldsymbol{\psi}\\. For two cell
types the map is monotone, so the endpoints of the interval for \\\rho\\
are transformed. For more than two cell types the joint normal draw
\\\boldsymbol{\rho}^{(b)}\\ is transformed and the interval is the
marginal quantile of \\\boldsymbol{\psi}(\boldsymbol{\rho}^{(b)})\\.
Every transformed draw is a composition, so the interval lies in
\\(0,1)\\.

This is the default of
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md).
The symmetric delta-method interval on \\\boldsymbol{p}\\ remains
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).
Columns of `y` are replicates that share one composition: their number
multiplies the Fisher information.

## Usage

``` r
confint_ilr_backtransform(
  y,
  mean_signature_matrix,
  Sigma,
  p = NULL,
  level = 0.95,
  n_draws = 2000L,
  seed = NULL,
  epsilon = 10^-8,
  itmax = 500L
)
```

## Arguments

- y:

  Numeric vector of length \\G\\, or \\G\times N\\ matrix of replicate
  bulk profiles that share one composition; their log-likelihoods are
  summed (the pooled, population-level fit).

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- p:

  Composition at which to evaluate the information. `NULL` uses the
  pooled maximum likelihood estimate.

- level:

  Confidence level.

- n_draws:

  Monte Carlo size used when there are three or more cell types. Ignored
  for two cell types.

- seed:

  Optional seed for the Monte Carlo map. The previous random-number
  state is restored. `NULL` leaves the state unchanged and draws from
  the current stream.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list with `estimate`, `interval` on the simplex, `rho`, `rho_interval`
(component-wise Wald bounds for \\\boldsymbol{\rho}\\), `vcov_rho`,
`vcov`, `mapping` (`"endpoint"` or `"monte_carlo"`) and `n_obs`. On a
simplex face the interval is `NA`.

## See also

[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md),
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
mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
colnames(mu) <- paste0("ct", 1:2)
Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
y <- drop(mu %*% c(0.6, 0.4))
confint_ilr_backtransform(y, mu, Sigma)$interval
#>          2.5%     97.5%
#> ct1 0.5514710 0.6453814
#> ct2 0.3546186 0.4485290
```
