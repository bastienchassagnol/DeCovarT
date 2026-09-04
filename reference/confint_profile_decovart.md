# Profile-likelihood confidence intervals for cellular ratios

Inverts the likelihood-ratio test of
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
coordinate by coordinate: the interval collects every \\c\\ with
\\2\\\ell(\hat{\boldsymbol{p}})-\ell^{\mathrm{prof}}\_j(c)\\
\le\chi^{2}\_{1,1-\alpha}\\. Because likelihood ratios are invariant
under reparametrisation, the interval respects \\0\le p_j\le 1\\ without
the delta-method linearisation that makes Wald intervals unreliable near
a simplex face.

## Usage

``` r
confint_profile_decovart(
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  level = 0.95,
  celltypes = NULL,
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

- level:

  Confidence level \\1-\alpha\\.

- celltypes:

  Cell types to profile; defaults to all.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A matrix with one row per cell type and columns `estimate`, `lower`,
`upper`.

## See also

[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
for the Wald analogue,
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
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
y <- drop(mu %*% c(0.6, 0.4))
confint_profile_decovart(y, mu, Sigma)
#>     estimate     lower     upper
#> ct1 0.599338 0.5530429 0.6473202
#> ct2 0.400662 0.3526798 0.4469571
```
