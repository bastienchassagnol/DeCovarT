# Confidence intervals for one single-cell bulk sample

The composition and the two technical scales are separated by the Schur
complement of the observed Hessian. `method` chooses the interval for
the composition. The default is the inverse ILR chart applied to a Wald
interval for \\\boldsymbol{\rho}\\, using the profile information
\\\\-\nabla^{2}\_{\boldsymbol{\rho}}\ell^{\ast}\\^{-1}\\.

Intervals for \\(c,\phi)\\ in `nuisance$interval` are Wald intervals on
the log scale, mapped back by the exponential, from the Schur complement
that removes \\\boldsymbol{\rho}\\. The rows `platform_scale_a`
(\\a=\phi/c\\) and `cell_count_n` (\\n=c^{2}/\phi\\, a number of cells
only when \\a\\ is known) follow by the delta method on \\\log
a=\log\phi-\log c\\ and \\\log n=2\log c-\log\phi\\, so all four
intervals are positive. `nuisance$interval_wald` is the symmetric normal
interval of \\(c,\phi)\\ on the original scale. Bootstrap methods
replace both with percentile intervals.

- `ilr`:

  Default. Inverse ILR chart.

- `wald`:

  Delta method on the composition. Endpoints can leave \\\[0,1\]\\.

- `firth`:

  Wald interval at the penalised profile estimate. `interval_ilr` maps
  that estimate through the inverse chart.

- `profile`:

  Interior likelihood-ratio interval. A face null is
  [`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md).

- `bootstrap`:

  Parametric percentile interval.

- `reference`:

  Donor bootstrap, or a within-type cell resample that regenerates the
  bulk. See
  [`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md).

These intervals treat the reference covariance as known. They do not add
a sandwich or higher-order correction for that plug-in.

`variance = "reml"` corrects for the \\\nu=G-J\\ residual degrees of
freedom of the variance scale: covariances are inflated by \\G/\nu\\,
Wald and ILR quantiles come from \\t\_{\nu}\\ (a multivariate \\t\\ for
the ILR map when \\J\>2\\), \\\phi\\ gets the chi-square interval
\\\[G\hat\phi/\chi^2\_{\nu;1-\alpha/2},
G\hat\phi/\chi^2\_{\nu;\alpha/2}\]\\, and the profile interval uses the
cut-off of
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md).
The maximum-likelihood \\\hat\phi\\ is biased by about \\\nu/G\\, so
with few genes (\\G=10\\, \\J=2\\) the `"ml"` intervals cover about 0.84
to 0.93 at the 0.95 level; the bootstrap methods are unaffected by
`variance`.

## Usage

``` r
confint_single_cell(
  y,
  mean_signature_matrix = NULL,
  Sigma = NULL,
  method = c("ilr", "wald", "firth", "profile", "bootstrap", "reference"),
  level = 0.95,
  variance = c("ml", "reml"),
  n_draws = 2000L,
  n_boot = 200L,
  seed = NULL,
  reference_profiles = NULL,
  donor_ids = NULL,
  resample = NULL,
  regenerate_bulk = NULL,
  ...
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- method:

  Interval for the composition. See Description.

- level:

  Confidence level.

- variance:

  `"ml"` (default) or `"reml"`; see Description.

- n_draws:

  Monte Carlo size for the ILR map when \\J\>2\\.

- n_boot:

  Bootstrap size.

- seed:

  Optional seed, restored on exit, for the ILR Monte Carlo map and for
  the bootstrap methods.

- reference_profiles, donor_ids, resample, regenerate_bulk:

  Passed to
  [`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md)
  when `method = "reference"`.

- ...:

  `epsilon` and `itmax` for the selected method.

## Value

A list with `method`, `estimate` (the fitted composition), `interval`
and `nuisance`. Bootstrap results also include the replicate matrices.

## See also

[`vcov_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_single_cell.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md),
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md),
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
fit <- confint_single_cell(y, mu, sigma)
fit$interval
#>          2.5%     97.5%
#> ct1 0.7032075 0.7074397
#> ct2 0.2925603 0.2967925
fit$nuisance$interval
#>                          2.5%        97.5%
#> mean_scale       1.998420e+01 2.008519e+01
#> variance_scale   2.652748e-02 4.241083e-01
#> platform_scale_a 1.324079e-03 2.116881e-02
#> cell_count_n     9.464154e+02 1.513110e+04
```
