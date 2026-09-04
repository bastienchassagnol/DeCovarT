# Refit DeCovarT from several random starts to probe multimodality

The realised DeCovarT log-likelihood is not globally concave, so a
successful convergence code, a small relative-distance-to-minimum and a
negative-definite local Hessian together establish a **local** maximum
only. Refitting from independent Dirichlet\\(1,\ldots,1)\\ starts and
comparing the attained log-likelihoods is the cheapest direct probe of a
second, equally good or better mode.

## Usage

``` r
multistart_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  n_starts = 5L,
  solver = deconvolute_ratios_Marquardt_Levenberg,
  loglik_tol = 1e-04,
  dirichlet_alpha = 1,
  ...
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

- n_starts:

  Number of Dirichlet restarts (in addition to the equi-balanced start).

- solver:

  Solver accepting `initial_p` and `return_model`; defaults to
  [`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md).

- loglik_tol:

  Two starts count as reaching the same mode when their log-likelihoods
  differ by less than this amount.

- dirichlet_alpha:

  Concentration for those random starts (default `1`). See
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md).

- ...:

  Passed to `solver` (e.g. `epsilon`, `itmax`).

## Value

A list with `coefficients` (best fit), `loglik`, `loglik_range`,
`starts` (matrix of per-start estimates), `logliks`, and `multimodal`.

## See also

[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

## Examples

``` r
mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
colnames(mu) <- paste0("ct", 1:2)
Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
y <- drop(mu %*% c(0.6, 0.4))
set.seed(1)
multistart_decovart(y, mu, Sigma, n_starts = 3L)$multimodal
#> [1] FALSE
```
