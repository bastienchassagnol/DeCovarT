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
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  p = NULL,
  null_value = NULL,
  n_boot = 500L,
  level = 0.95,
  n_replicates = NULL,
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

- p:

  Composition to simulate from. Defaults to the MLE computed from `y` /
  `bulk_expression`.

- null_value:

  Named numeric vector of hypothesised ratios, e.g. `c(celltype_2 = 0)`.

- n_boot:

  Number of bootstrap replicates.

- level:

  Confidence level for percentile intervals.

- n_replicates:

  Number of bulk columns simulated per replicate; defaults to the number
  of observed columns.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `estimates` (\\J\times\\ `n_boot` matrix), `interval`
(percentile bounds), `p_simulated`, and, when `null_value` is supplied,
`statistic` (observed \\D\\), `null_statistics`, and `p_value`.

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
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
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
set.seed(1)
bootstrap_decovart(y, mu, Sigma, n_boot = 20)$interval
#>          2.5%     97.5%
#> ct1 0.5647987 0.6266518
#> ct2 0.3733482 0.4352013
```
