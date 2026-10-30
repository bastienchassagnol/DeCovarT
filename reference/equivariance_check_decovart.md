# Permutation-equivariance check for a labelled reference

Software check, not a bootstrap and not a test of the bulk-to-signature
match. In a *reference-based* model the cell-type names are anchored by
the signature columns, so shuffling those names (or the gene axis) is
not a valid uncertainty procedure. What *is* required is that the
estimator be equivariant: if
\\\boldsymbol{\mu}^{\star}=\boldsymbol{\mu}Q\\ for a permutation matrix
\\Q\\, the fitted composition should satisfy
\\\hat{\boldsymbol{p}}^{\star}\approx Q^{\top}\hat{\boldsymbol{p}}\\ and
the reconstructed bulk should be unchanged.

## Usage

``` r
equivariance_check_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  perm = rev(seq_len(ncol(mean_signature_matrix))),
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

- perm:

  Integer permutation of `1:J`. Defaults to reversing the column order
  so the check is deterministic.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list with `perm`, `p_hat`, `p_star` (fit on the relabelled reference),
`p_expected` (\\\hat{\boldsymbol{p}}\\ reordered by `perm`),
`max_abs_diff`, and `loglik_diff`.

## Details

Columns of \\\boldsymbol{\mu}\\ and matching slices of
\\\boldsymbol{\Sigma}\_j\\ are reordered by `perm` while the *names*
stay in the original order, matching the convention of
[`vignette("01-2-population-inference")`](https://bastienchassagnol.github.io/DeCovarT/articles/01-2-population-inference.md).
The labelled MLE on the relabelled reference is then compared with
\\\hat{\boldsymbol{p}}\_{\mathrm{perm}}\\.

Use
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
for percentile intervals that resample the experimental units of the
reference (donors or cells).

## See also

[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)

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
y <- drop(mu %*% c(0.7, 0.3))
equivariance_check_decovart(y, mu, Sigma)$max_abs_diff
#> [1] 5.551115e-17
```
