# Parametric bootstrap for one single-cell bulk sample

Draws \\\boldsymbol{Y}^{(b)}\sim\mathcal{N}(\hat
c\tilde{\boldsymbol{\mu}}
\hat{\boldsymbol{p}},\hat\phi\boldsymbol{S}\_W(\hat{\boldsymbol{p}}))\\
and refits the profile model. Percentile intervals are returned for the
composition and for the two scales. The reference moments stay fixed.

With `mean_scale` and `variance_scale`, the scales are known: every
replicate is drawn at those scales and refitted with them held fixed.
`n_cells_sdlog > 0` then propagates an uncertain cell count: each
replicate multiplies both scales (\\c=an\\ and \\\phi=a^2n\\ are linear
in \\n\\) by a log-normal factor of mean one and that log-scale standard
deviation before drawing, and still refits at the supplied scales.

## Usage

``` r
bootstrap_single_cell(
  y,
  mean_signature_matrix,
  Sigma,
  n_boot = 200L,
  level = 0.95,
  epsilon = 1e-06,
  itmax = 40L,
  optimiser = c("newton", "lbfgsb", "marquardt"),
  mean_scale = NULL,
  variance_scale = NULL,
  n_cells_sdlog = 0
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- n_boot:

  Number of replicates.

- level:

  Confidence level.

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of iterations.

- optimiser:

  Passed to
  [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)
  for the fit and every refit, so the bootstrap matches the estimator
  being assessed.

- mean_scale, variance_scale:

  Optional fixed scales. Supply both or neither.

- n_cells_sdlog:

  Log-scale standard deviation of the cell count, used with known scales
  only.

## Value

A list with `estimates`, `interval`, `n_failed` (refits that errored,
left as `NA` columns), `scales`, `scale_interval` and `fit`. `scales`
and `scale_interval` have rows `mean_scale` (\\c\\), `variance_scale`
(\\\phi\\), `platform_scale_a` (\\a=\phi/c\\) and `cell_count_n`
(\\n=c^{2}/\phi\\).

## See also

[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
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
[`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md),
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)

## Examples

``` r
mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
colnames(mu) <- c("ct1", "ct2")
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
set.seed(1)
bootstrap_single_cell(y, mu, sigma, n_boot = 5)$interval
#>          2.5%     97.5%
#> ct1 0.7045195 0.7066303
#> ct2 0.2933697 0.2954805
```
