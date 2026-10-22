# Firth-penalised DeCovarT fit

Maximises \\\ell(\boldsymbol{\rho})+\tfrac12\log\det
I\_{\boldsymbol{\rho}}\\ on the ILR chart. The penalty is Jeffreys'
prior for the expected Fisher information of one bulk column; extra
replicate columns add only a constant, so they do not change the penalty
gradient. The returned Wald covariance is that information multiplied by
the number of columns and evaluated at the penalised estimate.

The penalty is differentiated numerically. The log-likelihood score
stays analytic. This is not an exact removal of the \\O(N^{-1})\\ bias:
the convolution is a curved exponential family.

## Usage

``` r
fit_firth_decovart(
  y,
  mean_signature_matrix,
  Sigma,
  p_start = NULL,
  epsilon = 1e-08,
  itmax = 200L
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

- p_start:

  Start on the simplex. `NULL` uses the unpenalised MLE, pulled into the
  open simplex when that MLE lies on a face.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

## Value

A list with `coefficients`, `rho`, `loglik` (unpenalised), `penalty`,
`penalised_loglik`, `vcov`, `vcov_rho` and `convergence`.

## References

Firth D (1993). “Bias Reduction of Maximum Likelihood Estimates.”
*Biometrika*, **80**(1), 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27) .

## See also

[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md)

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
fit_firth_decovart(y, mu, Sigma)$coefficients
#>       ct1       ct2 
#> 0.5986526 0.4013474 
```
