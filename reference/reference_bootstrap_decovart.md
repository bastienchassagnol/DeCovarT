# Reference-based bootstrap for signature and composition uncertainty

The frequentist DeCovarT likelihood treats \\\boldsymbol{\mu}\\ and each
\\\boldsymbol{\Sigma}\_j\\ as known. In a labelled, reference-based
model those moments come from purified profiles whose *cell-type names
stay attached*. This bootstrap resamples the experimental units that
generated the reference (donors, by default), or redraws compositions
from a Dirichlet law, then re-estimates the plug-in moments and refits
(Efron 1979) .

## Usage

``` r
reference_bootstrap_decovart(
  y = NULL,
  reference_profiles,
  bulk_expression = NULL,
  cell_type_labels = NULL,
  donor_ids = NULL,
  method = c("donors", "cells", "dirichlet"),
  n_boot = 199L,
  level = 0.95,
  dirichlet_alpha = 1,
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

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- reference_profiles:

  Named list of \\G\times n_j\\ matrices, one purified sample matrix per
  cell type. Names must match the intended column order of the
  signature. Alternatively, a single \\G\times n\_{\mathrm{ref}}\\
  matrix together with `cell_type_labels`.

- bulk_expression:

  Numeric matrix \\\boldsymbol{Y}\in \mathcal{M}\_{G\times N}\\ of
  replicate bulk profiles sharing one composition, used instead of `y`
  for the pooled (population) test.

- cell_type_labels:

  Optional character or factor vector of length
  `ncol(reference_profiles)` when `reference_profiles` is a matrix.

- donor_ids:

  Donor (or other clustering) labels aligned with the purified columns.
  A named list of vectors matching `reference_profiles`, or a vector of
  length `ncol(reference_profiles)` when the reference is a single
  matrix. Required for `method = "donors"`.

- method:

  One of `"donors"` (default), `"cells"`, or `"dirichlet"`.

- n_boot:

  Number of bootstrap replicates.

- level:

  Confidence level for percentile intervals.

- dirichlet_alpha:

  Positive Dirichlet concentration: a scalar recycled to \\J\\, or a
  length-\\J\\ vector. Ignored unless `method = "dirichlet"`. The
  default `1` is uniform on the simplex.

- regenerate_bulk:

  If `TRUE` (`donors` / `cells` only), simulate a new bulk from the
  bootstrapped moments instead of keeping the observed
  \\\boldsymbol{Y}\\ fixed.

- p:

  Composition to simulate from. Defaults to the MLE computed from `y` /
  `bulk_expression`.

- n_replicates:

  Number of bulk columns simulated per replicate; defaults to the number
  of observed columns.

- ridge:

  Non-negative ridge multiple for the sample covariance of each cell
  type. The default `1e-6` is relative to the mean diagonal.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `method`, `estimates` (\\J\times\\ `n_boot`), `interval`
(percentile bounds), `mean_signature_matrix` / `Sigma` from the original
reference, and `p_hat` from the original plug-in fit when a bulk is
supplied. For `method = "dirichlet"`, `p_simulated` stores the drawn
compositions.

## Details

Permuting gene order, or permuting cell-type labels of an already
averaged signature, is **not** a bootstrap: the first destroys the
gene-wise pairing of bulk and reference, and the second is algebraic
equivariance rather than a source of uncertainty (see
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
and
[`vignette("DeCovarT-MLE-properties")`](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-MLE-properties.md)).
The three `method` options instead resample units that actually vary.

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

- `dirichlet`:

  Hold the original plug-in moments fixed, draw
  \\\boldsymbol{p}^{(b)}\sim\mathrm{Dirichlet}(\boldsymbol{\alpha})\\,
  simulate bulk profiles from the Gaussian convolution, and refit. This
  is a composition-sweep diagnostic, not a confidence interval for one
  observed bulk.

For `donors` and `cells` the observed bulk is kept fixed unless
`regenerate_bulk = TRUE`, in which case each replicate also draws
\\\boldsymbol{Y}^{(b)}\\ from the bootstrapped moments at `p` (the
observed MLE by default). That is the pipeline reference units \\\to
S^{(b)}\to Y^{(b)}\to\hat{\boldsymbol{p}}^{(b)}\\.

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
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

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
#> ct1 0.6016023 0.6234908
#> ct2 0.3765092 0.3983977
```
