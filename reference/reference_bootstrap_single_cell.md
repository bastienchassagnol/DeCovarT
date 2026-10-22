# Reference bootstrap for the single-cell mixture

Resamples the reference that produced
\\(\tilde{\boldsymbol{\mu}},\boldsymbol{\Sigma}\_j)\\ and refits
[`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md).
`resample = NULL` uses donors when every cell type has at least two
donor labels, and otherwise resamples columns within each cell type. The
within-type resample regenerates the bulk from the fitted single-cell
law, at the original \\(\hat{\boldsymbol{p}},\hat c,\hat\phi)\\ and the
original moments, and refits it with the resampled moments, so that both
the bulk noise and the reference error enter the interval.

Shuffling cell-type names is not a reference bootstrap. It does not
redraw cells or donors, and the likelihood is equivariant under a
relabeling, so the draw would not be a source of uncertainty. That check
remains
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
for the population model.

## Usage

``` r
reference_bootstrap_single_cell(
  y,
  reference_profiles,
  donor_ids = NULL,
  resample = NULL,
  n_boot = 100L,
  level = 0.95,
  regenerate_bulk = NULL,
  ridge = 1e-06,
  epsilon = 1e-06,
  itmax = 40L
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- reference_profiles:

  Named list of gene-by-cell (or gene-by-replicate) matrices, one per
  cell type.

- donor_ids:

  Optional list of donor labels, aligned with `reference_profiles`.

- resample:

  `"donors"`, `"cells"`, or `NULL` for the automatic rule above.

- n_boot:

  Number of replicates.

- level:

  Confidence level.

- regenerate_bulk:

  `NULL` regenerates the bulk only for the within-type resample. `TRUE`
  or `FALSE` overrides that rule.

- ridge:

  Ridge passed to the sample covariance, relative to its mean diagonal.

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of Newton iterations.

## Value

A list with `resample`, `regenerate_bulk`, `estimates`, `interval`,
`n_failed`, `scales`, `scale_interval` and `fit`.

## See also

[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)

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
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md),
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)

## Examples

``` r
set.seed(1)
mu <- cbind(ct1 = c(2, 2, 0.2, 0.2), ct2 = c(0.2, 0.2, 3, 3))
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
refs <- lapply(seq_len(2), function(j) {
  t(MASS::mvrnorm(n = 6, mu = mu[, j], Sigma = sigma[, , j]))
})
names(refs) <- colnames(mu)
y <- 20 * drop(mu %*% c(0.7, 0.3)) + rnorm(4, sd = 0.2)
reference_bootstrap_single_cell(y, refs, n_boot = 4)$resample
#> [1] "cells"
```
