# Reference-based bootstrap for signature and composition uncertainty

The frequentist DeCovarT likelihood treats \\\boldsymbol{\mu}\\ and each
\\\boldsymbol{\Sigma}\_j\\ as known. In a labelled, reference-based
model those moments come from purified profiles whose *cell-type names
stay attached*. This bootstrap resamples the experimental units that
generated the reference (donors, by default), re-estimates the plug-in
moments and refits (Efron 1979) .

## Usage

``` r
reference_bootstrap_decovart(
  y,
  reference_profiles,
  donor_ids = NULL,
  method = c("donors", "cells"),
  n_boot = 199L,
  level = 0.95,
  regenerate_bulk = FALSE,
  p = NULL,
  n_replicates = NULL,
  ridge = 1e-06,
  epsilon = 10^-8,
  itmax = 500L
)
```

## Arguments

- y:

  Numeric vector of length \\G\\, or \\G\times N\\ matrix of replicate
  bulk profiles that share one composition; their log-likelihoods are
  summed (the pooled, population-level fit).

- reference_profiles:

  Named list of \\G\times n_j\\ matrices, one purified sample matrix per
  cell type (at least two columns each). The list order fixes the column
  order of the signature.

- donor_ids:

  Donor (or other clustering) labels: a list named like
  `reference_profiles`, with one label per purified column. Required for
  `method = "donors"`.

- method:

  One of `"donors"` (default) or `"cells"`.

- n_boot:

  Number of bootstrap replicates.

- level:

  Confidence level for percentile intervals.

- regenerate_bulk:

  If `TRUE`, simulate a new bulk from the fitted model and refit it with
  the bootstrapped moments, instead of keeping the observed
  \\\boldsymbol{Y}\\ fixed.

- p:

  Composition to simulate from. Defaults to the MLE computed from `y`
  (or the restricted MLE under `null_value`).

- n_replicates:

  Number of bulk columns simulated per replicate; defaults to the number
  of observed columns.

- ridge:

  Non-negative ridge multiple for the sample covariance of each cell
  type. The default `1e-6` is relative to the mean diagonal.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list with `method`, `estimates` (\\J\times\\ `n_boot`), `interval`
(percentile bounds), `n_failed`, `mean_signature_matrix` / `Sigma` from
the original reference, and `p_hat` from the original plug-in fit.

## Details

Permuting gene order, or permuting cell-type labels of an already
averaged signature, is **not** a bootstrap: the first destroys the
gene-wise pairing of bulk and reference, and the second is algebraic
equivariance rather than a source of uncertainty (see
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
and
[`vignette("theory-population-inference")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.md)).
The two `method` options instead resample, with replacement, units that
actually vary.

- `donors`:

  Cluster bootstrap: resample donor identifiers *within each cell type*,
  take every purified column of the sampled donors (a donor drawn twice
  is included twice), rebuild
  \\(\boldsymbol{\mu},\boldsymbol{\Sigma}\_j)\\, and refit. This is the
  default: it targets biological / reference-population uncertainty
  rather than technical cell-level noise alone.

- `cells`:

  Resample purified columns independently within each cell type (with
  replacement) and rebuild the plug-in moments. Use this when donor
  labels are unavailable.

The observed bulk is kept fixed unless `regenerate_bulk = TRUE`, in
which case each replicate also draws \\\boldsymbol{Y}^{(b)}\\ from the
fitted model, at `p` (the observed MLE by default) and the moments of
the full reference, and refits it with the bootstrapped moments
\\S^{(b)}\\. Drawing \\\boldsymbol{Y}^{(b)}\\ from \\S^{(b)}\\ instead
would make each refit consistent with its own moments, so the reference
error would not reach \\\hat{\boldsymbol{p}}^{(b)}\\: with 60 cells per
type, that variant covered a known composition 76\\ level, against 92\\

When a cell type has fewer purified columns than genes, the sample
covariance is singular. A ridge `ridge * mean(diag(S))` is then added on
the diagonal so that \\\boldsymbol{\Sigma}(\boldsymbol{p})\\ remains
factorisable. Set `ridge = 0` to disable that guard (and accept possible
Cholesky failures).

## References

Efron B (1979). “Bootstrap Methods: Another Look at the Jackknife.” *The
Annals of Statistics*, **7**(1), 1–26.
[doi:10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552) .

## See also

[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)

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
[`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md),
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)

## Examples

``` r
set.seed(1)
mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
colnames(mu) <- paste0("ct", 1:2)
Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
refs <- lapply(seq_len(2), function(j) {
  t(MASS::mvrnorm(n = 8, mu = mu[, j], Sigma = Sigma[, , j]))
})
names(refs) <- colnames(mu)
donor_ids <- lapply(refs, function(x) {
  rep(paste0("d", 1:2), length.out = ncol(x))
})
y <- drop(mu %*% c(0.6, 0.4))
reference_bootstrap_decovart(
  y,
  refs,
  donor_ids = donor_ids,
  n_boot = 15
)$interval
#>          2.5%     97.5%
#> ct1 0.6016025 0.6234908
#> ct2 0.3765092 0.3983975
```
