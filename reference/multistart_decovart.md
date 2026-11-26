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
  solver = deconvolute_ratios_Newton_Raphson,
  loglik_tol = 1e-04,
  dirichlet_alpha = 1,
  initial_p = NULL,
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

  Number of Dirichlet restarts (in addition to the first start).

- solver:

  Solver accepting `y`, `mean_signature_matrix`, `Sigma` and
  `initial_p`, and returning `coefficients` (population solvers, called
  with `return_model = TRUE`) or `p` (as
  [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md))
  together with `loglik`. Defaults to
  [`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md).
  With `solver = deconvolute_single_cell` the single-cell Newton fit is
  restarted; pass `optimiser` through `...` to change it.

- loglik_tol:

  Two starts count as reaching the same mode when their log-likelihoods
  differ by less than this amount.

- dirichlet_alpha:

  Concentration for those random starts (default `1`). See
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md).

- initial_p:

  First start: `NULL` (default) for the barycentre, a length-\\J\\
  numeric vector (for example the estimate of another method), or
  `"qp"`; see
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md).

- ...:

  Passed to `solver` (e.g. `epsilon`, `itmax`).

## Value

A list with `coefficients` (best fit), `loglik`, `loglik_range`,
`starts` (matrix of per-start estimates), `logliks`, `multimodal`,
`best_start` (index of the best start, 1 for the first) and `fit` (the
full solver output of the best start).

## See also

[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)

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
set.seed(1)
multistart_decovart(y, mu, Sigma, n_starts = 3L)$multimodal
#> [1] FALSE
# First start supplied by the user, e.g. another method's estimate
multistart_decovart(y, mu, Sigma, n_starts = 2L, initial_p = c(0.5, 0.5))$loglik
#> [1] 0.9812718
```
