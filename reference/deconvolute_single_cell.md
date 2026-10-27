# Single-cell deconvolution of one bulk sample

Fits one bulk sample from the linear single-cell mixture
\\\boldsymbol{y}\sim\mathcal{N}(c\tilde{\boldsymbol{\mu}}\boldsymbol{p},
\phi\boldsymbol{S}\_W(\boldsymbol{p}))\\.

The start is
[`single_cell_nuisance_scales()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_nuisance_scales.md):
nonnegative least squares with a free scale (DeconRNASeq's quadratic
programme without the simplex equality), then the closed-form scales at
that composition.

With `profile = TRUE` (default) the two scales are profiled out: for
every trial composition they are recomputed in closed form, so the
weight matrix \\\boldsymbol{S}\_W(\boldsymbol{p})^{-1}\\ is rebuilt at
every iterate, and Newton-Raphson runs on the profile log-likelihood of
the \\J-1\\ ILR coordinates
([`profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md)).
Its limit is the joint maximum-likelihood estimate of
\\(\boldsymbol{p},c,\phi)\\, reached in one loop.

With `profile = FALSE` the scales stay at their start values (a
two-stage plug-in fit), which is not the joint maximum. Supplying both
`mean_scale` and `variance_scale` fixes them at those values.

With `optimiser = "newton"` (default) every fit uses the analytic score
and observed Hessian with the PORT trust region of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html): the step
maximises the quadratic model inside a ball whose radius grows or
shrinks with the ratio of actual to predicted increase; there is no line
search and no fixed step length. With `optimiser = "lbfgsb"` the same
objective is maximised by L-BFGS-B on the closure
\\\boldsymbol{p}=\boldsymbol{w}/\mathbf{1}^{\top} \boldsymbol{w}\\,
\\\boldsymbol{w}\ge0\\, as in
[`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md);
it can return exact zeros, where `rho` and the covariances are `NA`. The
returned `vcov` is the delta-method image of the inverse negative
profile Hessian; intervals are
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md).

## Usage

``` r
deconvolute_single_cell(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 1e-08,
  itmax = 100L,
  profile = TRUE,
  initial_p = NULL,
  mean_scale = NULL,
  variance_scale = NULL,
  estimand = c("cellular", "transcript"),
  cell_size = NULL,
  optimiser = c("newton", "lbfgsb")
)
```

## Arguments

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- epsilon:

  Relative function and step tolerances passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- itmax:

  Maximum number of iterations.

- profile:

  If `TRUE`, profile the two scales out; if `FALSE`, keep them at their
  start values.

- initial_p:

  Optional open-simplex start. `NULL` uses the free-scale least-squares
  start.

- mean_scale, variance_scale:

  Optional fixed scales. Supply both or neither.

- estimand:

  `cellular` (default) interprets the fitted composition as a cell
  fraction, which is valid when the mean keeps the cell sizes \\S_j\\.
  `transcript` interprets it as a transcript fraction; pass `cell_size`
  to convert it back to a cell fraction.

- cell_size:

  Optional \\S_j\\. Defaults to the column sums of the mean, which equal
  the cell sizes for a `cellular_raw` reference.

- optimiser:

  `"newton"` (trust-region Newton in ILR coordinates) or `"lbfgsb"`
  (L-BFGS-B on the closed simplex).

## Value

A list with `composition` (the fitted simplex coordinate), `p`, `q`,
`rho`, `estimand`, `fit_type` (`"profile"`, `"two_stage"` or
`"fixed_scales"`), `mean_scale`, `variance_scale`, `platform_scale`,
`n_cells`, `cell_size`, `ols_scale`, `loglik`, `convergence` (0 from
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) means success),
`iterations`, `p_start`, `message`, `vcov_rho`, `vcov` and
`vcov_nuisance` (see
[`vcov_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_single_cell.md)).
A singular information matrix is returned as `NA` without failing the
fit.

## See also

[`single_cell_nuisance_scales()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_nuisance_scales.md),
[`profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md),
[`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md),
[`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md)

## Examples

``` r
genes <- paste0("g", 1:4)
cts <- paste0("ct", 1:2)
mu <- matrix(
  c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
  nrow = 4,
  dimnames = list(genes, cts)
)
sigma <- array(0, c(4, 4, 2), dimnames = list(genes, genes, cts))
sigma[,, 1] <- diag(0.05, 4)
sigma[,, 2] <- diag(0.05, 4)
set.seed(1)
sim <- simulate_bulk_from_cells(
  mu,
  sigma,
  p = c(0.8, 0.2),
  n_cells = 2000,
  platform_scale = 1
)
fit <- deconvolute_single_cell(sim$Y[, 1], mu, sigma)
round(fit$p, 3)
#>   ct1   ct2 
#> 0.801 0.199 
```
