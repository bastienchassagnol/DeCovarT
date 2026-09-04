# Fit the DeCovarT Gaussian-convolution model

Maximum-likelihood estimator of cellular proportions
\\\boldsymbol{p}\_{\cdot i}\\ under the multivariate Gaussian
convolution \$\$ \boldsymbol{y}\_{\cdot
i}\mid(\boldsymbol{\mu},\boldsymbol{\Sigma}\_j, \boldsymbol{p}\_{\cdot
i}) \sim \mathcal{N}\_{G}\bigl( \boldsymbol{\mu}\boldsymbol{p}\_{\cdot
i},\\ \boldsymbol{\Sigma}(\boldsymbol{p}\_{\cdot i}) \bigr), \qquad
\boldsymbol{\Sigma}(\boldsymbol{p})
=\sum\_{j=1}^{J}p_j^{2}\boldsymbol{\Sigma}\_j. \$\$ This is a **variance
/ likelihood specification**, not ordinary least squares: cell-type
covariances enter the residual law and cannot be written as extra
columns of \\\boldsymbol{\mu}\\. There is therefore no `formula` / `lm`
interface, and no [`predict()`](https://rdrr.io/r/stats/predict.html)
method for forecasting bulk expression.

The returned object is an S3 class `decovart_fit` with accessors
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). Residuals are
\\\boldsymbol{Y}-\boldsymbol{\mu}\hat{\boldsymbol{P}}\\; they are
**not** OLS residuals. Goodness of fit is the convolution
log-likelihood, not residual sum of squares.

## Usage

``` r
fit_decovart(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  method = c("Marquardt-Levenberg", "L-BFGS-B", "Newton-Raphson"),
  epsilon = 10^-4,
  itmax = 200,
  standardise = FALSE,
  scaled = FALSE,
  n_starts = 0L,
  boundary_tol = 1e-08
)

# S3 method for class 'decovart_fit'
print(x, ...)

# S3 method for class 'decovart_fit'
summary(object, ...)

# S3 method for class 'summary.decovart_fit'
print(x, ...)

# S3 method for class 'decovart_fit'
coef(object, ...)

# S3 method for class 'decovart_fit'
fitted(object, ...)

# S3 method for class 'decovart_fit'
residuals(object, ...)

# S3 method for class 'decovart_fit'
vcov(object, ...)

# S3 method for class 'decovart_fit'
nobs(object, ...)

# S3 method for class 'decovart_fit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'decovart_fit'
plot(x, ...)
```

## Arguments

- signature_matrix:

  Mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (rownames = genes, colnames = cell types).

- bulk_expression:

  Bulk matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\ (rownames =
  genes, colnames = samples).

- true_ratios:

  Optional ground-truth proportions (\\J\\ vector or \\J\times N\\
  matrix). Stored on the fit; not used by the MLE.

- Sigma:

  Array \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times
  J}\\ of cell-type covariances.

- method:

  Optimiser; one of `"Marquardt-Levenberg"`, `"L-BFGS-B"`,
  `"Newton-Raphson"` (case-insensitive). These three maps already land
  on the simplex (ILR or \\p/\sum p\\); they do **not** call
  [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md).

- epsilon, itmax:

  Absolute convergence tolerance and iteration budget, in the same roles
  as `reltol` / `maxit` in
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- standardise:

  If `TRUE`, apply a **gene-wise** affine z-score computed once from
  \\\boldsymbol{\mu}\\ to bulk, means and covariances (see Details).
  Cell-type-wise or sample-wise transforms are not supported.

- scaled:

  Deprecated. `TRUE` (log2 mixing) always errors.

- n_starts:

  Number of additional random Dirichlet starts per sample. With
  `n_starts > 0` the best log-likelihood is kept and
  [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  records the spread of attained optima, the only direct probe of
  multimodality (the realised log-likelihood is not globally concave).

- boundary_tol:

  Threshold on \\\min_j\hat{p}\_j\\ below which a sample is flagged
  `near_boundary` by
  [`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md).
  Wald intervals are unreliable there; prefer
  [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)
  or
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md).

- x, object:

  A `decovart_fit`.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)
  (plot method) or ignored.

- parm:

  Kept for compatibility with
  [`stats::confint()`](https://rdrr.io/r/stats/confint.html); all
  simplex coordinates are always returned (subsetting by `parm` is not
  implemented). Removing this formal would break the S3 method contract
  with the generic.

- level:

  Confidence level \\1-\alpha\\ (default `0.95`). Wald intervals use
  asymptotic normality of the MLE with standard errors from
  [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
  /
  [`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
  (see Details of `fit_decovart()` and of
  [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)):
  \\\hat{p}\_j \pm z\_{1-\alpha/2}\\\mathrm{SE}\_j\\.

## Value

An object of class `decovart_fit`.

## Details

**Standardisation.** CIBERSORT requires non-negative expression, no
missing values, and a non-log linear scale (Newman et al. 2015) . A
logarithm is concave, so Jensen's inequality shifts first and second
moments and
\\\log(\boldsymbol{\mu}\boldsymbol{p})\neq(\log\boldsymbol{\mu})\boldsymbol{p}\\.
A gene-wise affine map \\\boldsymbol{x}\mapsto
D^{-1}(\boldsymbol{x}-\boldsymbol{m})\\ with the **same**
\\D,\boldsymbol{m}\\ on \\\boldsymbol{Y}\\, \\\boldsymbol{\mu}\\ and
\\\boldsymbol{\Sigma}\_j^\star=D^{-1}\boldsymbol{\Sigma}\_j D^{-1}\\
leaves the MLE of \\\boldsymbol{p}\\ unchanged (equivariance).

**Wald covariance.** Let
\\\boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\.
The expected Fisher information of the unconstrained mean–covariance map
(multivariate normal; see e.g. the Wikipedia entry *Fisher information*,
multivariate normal) is \$\$ I(\boldsymbol{p})\_{jk} =
\boldsymbol{\mu}\_{\cdot j}^{\top} \boldsymbol{\Theta}(\boldsymbol{p})
\boldsymbol{\mu}\_{\cdot k} + 2 p_j p_k\\ \mathrm{tr}\bigl(
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_j
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_k \bigr). \$\$
Cramer–Rao gives \\\mathrm{Var}(\hat{\boldsymbol{z}})\succeq
I\_{\boldsymbol{z}}^{-1}\\ in ILR coordinates, with
\\I\_{\boldsymbol{z}} =\mathbf{J}\_{\boldsymbol{\psi}}^{\top}
I(\boldsymbol{p}) \mathbf{J}\_{\boldsymbol{\psi}}\\ and
\\\mathbf{J}\_{\boldsymbol{\psi}}
=\mathbf{S}(\boldsymbol{p})\mathbf{V}\\
([`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)).
The delta method maps the bound back to the simplex: \$\$
\mathrm{Var}(\hat{\boldsymbol{p}}) = \mathbf{J}\_{\boldsymbol{\psi}}
I\_{\boldsymbol{z}}^{-1} \mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}.
\$\$

## References

Newman A, Liu C, Green M, Gentles A, Feng W, Xu Y, Hoang C, Diehn M,
Alizadeh A (2015). “Robust Enumeration of Cell Subsets from Tissue
Expression Profiles.” *Nature methods*, **12**.
[doi:10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337) .

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md),
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md),
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md),
`coef.decovart_fit()`, `vcov.decovart_fit()`, `confint.decovart_fit()`

## Examples

``` r
toy <- readRDS(system.file(
  "extdata", "toy_deconvolution.rds",
  package = "DeCovarT"
))
fit <- fit_decovart(
  signature_matrix = toy$signature_matrix,
  bulk_expression = toy$bulk_expression,
  Sigma = toy$Sigma,
  itmax = 40
)
coef(fit)
#>      sample_1  sample_2
#> ct1 0.4461301 0.7062033
#> ct2 0.5538699 0.2937967
nobs(fit)
#> [1] 2
#> attr(,"n_genes")
#> [1] 2
#> attr(,"n_celltypes")
#> [1] 2
#> attr(,"n_samples")
#> [1] 2
```
