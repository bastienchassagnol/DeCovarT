# Likelihood-ratio test for cellular ratios

Tests \\H_0:p_j=c_j,\\j\in A\\ against the unrestricted simplex with the
profile likelihood-ratio statistic \$\$ D = 2\bigl\\
\ell(\hat{\boldsymbol{p}}) -\ell(\hat{\boldsymbol{p}}^{(0)}) \bigr\\,
\$\$ where \\\hat{\boldsymbol{p}}^{(0)}\\ is the restricted MLE of
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md).
Unlike the Wald interval of
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
\\D\\ is invariant under the additive log-ratio reparametrisation,
because a smooth bijection leaves likelihood ratios unchanged.

## Usage

``` r
lrt_decovart(
  y = NULL,
  mean_signature_matrix,
  Sigma,
  null_value,
  bulk_expression = NULL,
  n_boundary = NULL,
  boundary_tol = 1e-06,
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

- null_value:

  Named numeric vector of hypothesised ratios, e.g. `c(celltype_2 = 0)`.

- bulk_expression:

  Numeric matrix \\\boldsymbol{Y}\in \mathcal{M}\_{G\times N}\\ of
  replicate bulk profiles sharing one composition, used instead of `y`
  for the pooled (population) test.

- n_boundary:

  Optional integer override for the number of active boundary
  constraints; `NULL` infers it from `null_value`.

- boundary_tol:

  Distance from 0 or 1 below which a null value counts as being on the
  boundary.

- epsilon, itmax:

  Relative convergence tolerance and iteration budget passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A one-row data frame with columns `statistic`, `n_boundary`,
`df_interior`, `p_value`, `loglik_full`, `loglik_null`, and
`calibration`.

## Details

**Interior null.** If every tested ratio lies strictly inside \\(0,1)\\,
\\D\\ is asymptotically \\\chi^{2}\_{\|A\|}\\ (Wilks 1938) .

**Boundary null.** A null such as \\p_j=0\\ lets the parameter move in
one direction only, so the local parameter space is a half-line rather
than a vector space and Wilks' theorem fails. The limit is the
chi-bar-square mixture of
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md):
for a single active constraint,
\\D\rightsquigarrow\tfrac{1}{2}\chi^{2}\_{0}+\tfrac{1}{2}\chi^{2}\_{1}\\
Chernoff (1954) , generalised to several active constraints by Self and
Liang (1987) . Boundary calibration is selected automatically when any
tested value is within `boundary_tol` of 0 or 1; override with
`n_boundary`.

**Replication.** Both calibrations are *population* statements: they
need replicate bulk samples that genuinely share one composition
\\\boldsymbol{p}\\. A single \\G\\-vector supplies no replication (genes
are dependent through \\\boldsymbol{\Sigma}(\boldsymbol{p})\\, and
adding bulk samples with distinct compositions adds parameters as fast
as observations). Use `bulk_expression` with several replicate columns,
or calibrate by
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md).

## References

Wilks SS (1938). “The Large-Sample Distribution of the Likelihood Ratio
for Testing Composite Hypotheses.” *The Annals of Mathematical
Statistics*, **9**(1), 60–62. ISSN 0003-4851.
[doi:10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360) .

Self SG, Liang K (1987). “Asymptotic Properties of Maximum Likelihood
Estimators and Likelihood Ratio Tests under Nonstandard Conditions.”
*Journal of the American Statistical Association*, **82**(398), 605–610.
ISSN 0162-1459.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)
.

## See also

[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

## Examples

``` r
mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
colnames(mu) <- paste0("ct", 1:3)
Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
y <- drop(mu %*% c(0.5, 0.5, 0))
lrt_decovart(y, mu, Sigma, null_value = c(ct3 = 0))
#>     statistic n_boundary df_interior   p_value loglik_full loglik_null
#> ct3 0.2785736          1           0 0.2988188    0.832434   0.6931472
#>                     calibration
#> ct3 chi-bar-square (Self-Liang)
```
