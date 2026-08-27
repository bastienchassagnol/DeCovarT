# Boundary and stationarity diagnostics for one DeCovarT fit

Separates three claims that a bare optimiser return code conflates: the
solver stopped, the iterate is genuinely stationary with the right local
curvature, and the estimate sits close to a simplex face. None of them
implies global uniqueness, which the DeCovarT log-likelihood does not
have for a single observed sample (see
[`vignette("DeCovarT-MLE-properties")`](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-MLE-properties.md)).

## Usage

``` r
boundary_diagnostics(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  boundary_tol = 1e-08,
  score_tol = 1e-04
)
```

## Arguments

- p:

  Estimated proportions of length \\J\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- boundary_tol:

  Threshold on \\\min_j\hat{p}\_j\\ below which the estimate is flagged
  as near-boundary.

- score_tol:

  Threshold on the ALR score norm below which the iterate counts as
  stationary.

## Value

A one-row data frame of diagnostics.

## Details

Reported fields are the ALR score norm
\\\lVert\nabla\_{\boldsymbol{\rho}}\ell\rVert\\, the largest eigenvalue
\\\lambda\_{\max}(\mathbf{H}\_{\boldsymbol{\rho}})\\ (negative at a
local maximum), `boundary_distance` \\=\min_j\hat{p}\_j\\, and the flags
`near_boundary` and `local_maximum`.

`boundary_tol` is a **statistical** warning threshold for Wald / ALR
linearisation, deliberately much larger than the machine-precision guard
that decides whether a logarithm is representable. A fit with
\\\min_j\hat{p}\_j\ll 1\\ is not evidence of optimiser failure: the
solver may be correctly approaching a genuine boundary optimum. It *is*
evidence that interior Wald intervals should be replaced by the profile
or boundary-calibrated tests of
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
and
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md).

## See also

[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
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
boundary_diagnostics(c(0.6, 0.4), y, mu, Sigma)
#>   boundary_distance near_boundary score_norm max_eigenvalue local_maximum
#> 1               0.4         FALSE  0.2769231      -100.2504         FALSE
```
