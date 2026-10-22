# Choose a confidence interval for a pooled DeCovarT composition

One entry point for the interior intervals of the population convolution
(covariance weights \\p_j^2\\, no technical scales). The default,
`method = "ilr"`, is
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md).

- `ilr`:

  Wald interval for \\\boldsymbol{\rho}\\, mapped by the inverse ILR
  chart.

- `wald`:

  Delta method on \\\boldsymbol{p}\\. Endpoints can leave \\\[0,1\]\\.

- `firth`:

  Wald interval at the Jeffreys-penalised estimate
  ([`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md)).
  `interval_ilr` is the same estimate mapped through the inverse chart.

- `profile`:

  Likelihood-ratio interval
  ([`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)).
  This is the interior calibration. A null on a face is
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
  not an interval.

- `bootstrap`:

  Parametric percentile interval
  ([`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)).

- `reference`:

  Resample the reference, then refit
  ([`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)).
  Two or more donors in every cell type use the donor bootstrap and keep
  the bulk fixed. Otherwise columns are resampled within each cell type
  and the bulk is regenerated. Cell-type labels are not shuffled; that
  check is
  [`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md).

Columns of `y` share one composition.
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
instead estimates a separate composition per column, and
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
is the per-column delta method.

## Usage

``` r
confint_decovart(
  y,
  mean_signature_matrix = NULL,
  Sigma = NULL,
  method = c("ilr", "wald", "firth", "profile", "bootstrap", "reference"),
  level = 0.95,
  p = NULL,
  n_draws = 2000L,
  n_boot = 500L,
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

  Numeric vector of length \\G\\, or \\G\times N\\ matrix of replicate
  bulk profiles that share one composition; their log-likelihoods are
  summed (the pooled, population-level fit).

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- method:

  Interval to compute. See Description.

- level:

  Confidence level.

- p:

  Composition at which to evaluate the information. `NULL` uses the
  pooled maximum likelihood estimate.

- n_draws:

  Monte Carlo size used when there are three or more cell types. Ignored
  for two cell types.

- n_boot:

  Bootstrap size for `method = "bootstrap"` or `"reference"`.

- seed:

  Optional seed for the Monte Carlo map. The previous random-number
  state is restored. `NULL` leaves the state unchanged and draws from
  the current stream.

- reference_profiles, donor_ids:

  Passed to
  [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  when `method = "reference"`.

- resample:

  `"donors"` or `"cells"`. `NULL` chooses from `donor_ids` as in
  Description.

- regenerate_bulk:

  Whether the reference bootstrap redraws the bulk. `NULL` follows the
  automatic rule in Description.

- ...:

  Passed to the selected method (`epsilon`, `itmax`, `celltypes`,
  `null_value`, `ridge`, `n_replicates`).

## Value

A list with `method`, `estimate` and `interval`. Other elements depend
on `method`.

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md),
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
confint_decovart(y, mu, Sigma)$interval
#>          2.5%     97.5%
#> ct1 0.5514710 0.6453814
#> ct2 0.3546186 0.4485290
confint_decovart(y, mu, Sigma, method = "wald")$interval
#>          2.5%     97.5%
#> ct1 0.5522501 0.6464257
#> ct2 0.3535743 0.4477499
```
