# Parametric bootstrap for DeCovarT proportions and boundary tests

The conditional law of the bulk profile is fully specified, so the
asymptotic calibrations of
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
can be replaced by simulation. Given plug-in moments and a composition
\\\boldsymbol{p}\\, the function draws \\\boldsymbol{Y}^{(b)}\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\, refits DeCovarT, and returns the
empirical distribution of \\\hat{\boldsymbol{p}}^{(b)}\\ together with
percentile intervals.

## Usage

``` r
bootstrap_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  p = NULL,
  null_value = NULL,
  n_boot = 500L,
  level = 0.95,
  n_replicates = NULL,
  epsilon = 10^-8,
  itmax = 500L,
  optimiser = c("newton", "lbfgsb")
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

- p:

  Composition to simulate from. Defaults to the MLE computed from `y`
  (or the restricted MLE under `null_value`).

- null_value:

  Optional numeric vector of hypothesised ratios named by cell types;
  switches to the restricted bootstrap (see Details).

- n_boot:

  Number of bootstrap replicates.

- level:

  Confidence level for percentile intervals.

- n_replicates:

  Number of bulk columns simulated per replicate; defaults to the number
  of observed columns.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- optimiser:

  Solver for the unrestricted fits: `"newton"` (trust-region Newton in
  ILR coordinates) or `"lbfgsb"` (L-BFGS-B on the closed simplex, as in
  [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)),
  so that the bootstrap distribution matches the estimator being
  assessed. Restricted fits under `null_value` always use the Newton
  driver.

## Value

A list with `estimates` (\\J\times\\ `n_boot` matrix), `interval`
(percentile bounds), `n_failed` (refits that errored, left as `NA`
columns), `p_simulated`, and, when `null_value` is supplied, `statistic`
(observed \\D\\), `null_statistics`, and `p_value`.

## Details

Supplying `null_value` switches to a restricted bootstrap: replicates
are generated under the restricted MLE and the likelihood-ratio
statistic is recomputed on each, giving a Monte Carlo p-value that does
not rely on the \\50{:}50\\ asymptotic weights. For a composite null the
remaining nuisance proportions are profiled at their restricted
estimate, so the procedure is a plug-in (approximate) rather than an
exact bootstrap.

Reference uncertainty is **not** propagated here: \\\boldsymbol{\mu}\\
and \\\boldsymbol{\Sigma}\_j\\ are treated as known, exactly as in the
frequentist likelihood. To resample the purified observations that
produced those plug-in moments, use
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md).

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)

Other decovart_inference:
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
bootstrap_decovart(y, mu, Sigma, n_boot = 20)$interval
#>          2.5%     97.5%
#> ct1 0.5685610 0.6251911
#> ct2 0.3748089 0.4314390
```
