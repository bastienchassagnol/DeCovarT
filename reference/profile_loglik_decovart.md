# Profile log-likelihood of one cellular ratio

Returns \\\ell^{\mathrm{prof}}\_j(c)=\max\\\ell(\boldsymbol{p}):
\boldsymbol{p}\in\Delta^{J-1},\\p_j=c\\\\, i.e. the restricted maximum
of
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)
for a single coordinate. Profiling concentrates out the remaining
\\J-2\\ free proportions instead of relying on a quadratic (Wald)
approximation, so it is reliable even when \\\hat{p}\_j\\ sits close to
a simplex face.

## Usage

``` r
profile_loglik_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  celltype,
  value,
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

- celltype:

  Cell-type name or integer index \\j\\.

- value:

  Numeric vector of candidate ratios \\c\in\[0,1\]\\.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

Numeric vector of profile log-likelihoods, named by `value`.

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

## Examples

``` r
mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
colnames(mu) <- paste0("ct", 1:3)
Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
y <- drop(mu %*% c(0.5, 0.3, 0.2))
profile_loglik_decovart(y, mu, Sigma, "ct3", c(0, 0.2, 0.4))
#>       0.0       0.2       0.4 
#> 0.4623643 0.9941606 0.8435155 
```
