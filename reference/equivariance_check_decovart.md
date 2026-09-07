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
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  perm = NULL,
  epsilon = 10^-8,
  itmax = 500L
)
```

## Arguments

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- bulk_expression:

  Numeric matrix \\\boldsymbol{Y}\in \mathcal{M}\_{G\times N}\\ of
  replicate bulk profiles sharing one composition, used instead of `y`
  for the pooled (population) test.

- perm:

  Integer permutation of `1:J`. Defaults to reversing the column order
  so the check is deterministic.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `perm`, `p_hat`, `p_star` (fit on the relabelled reference),
`p_expected` (\\\hat{\boldsymbol{p}}\\ reordered by `perm`),
`max_abs_diff`, and `loglik_diff`.

## Details

Columns of \\\boldsymbol{\mu}\\ and matching slices of
\\\boldsymbol{\Sigma}\_j\\ are reordered by `perm` while the *names*
stay in the original order, matching the convention of
[`vignette("theory-DeCovarT-MLE-properties")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.md).
The labelled MLE on the relabelled reference is then compared with
\\\hat{\boldsymbol{p}}\_{\mathrm{perm}}\\.

Use
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
for percentile intervals that resample the experimental units of the
reference (donors or cells) or redraw compositions from a Dirichlet law.

## See also

[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

## Examples

``` r
mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
colnames(mu) <- paste0("ct", 1:2)
Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
y <- drop(mu %*% c(0.7, 0.3))
equivariance_check_decovart(y, mu, Sigma)$max_abs_diff
#> [1] 0
```
