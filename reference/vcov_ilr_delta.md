# Cramer–Rao / ILR delta-method covariance of \\\hat{\boldsymbol{p}}\\

Maps the expected Fisher information of unconstrained proportions
through the isometric log-ratio (ILR) chart and back to the simplex via
the delta method.

Let \\\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{z})\\
with Jacobian
\\\mathbf{J}\_{\boldsymbol{\psi}}=\mathbf{S}(\boldsymbol{p})
\mathbf{V}\\
([`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)).
Fisher information transforms as \$\$ I\_{\boldsymbol{z}} =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}} I(\boldsymbol{p})
\mathbf{J}\_{\boldsymbol{\psi}} = \mathbf{V}^{\mathsf{T}}
\mathbf{S}(\boldsymbol{p}) I(\boldsymbol{p}) \mathbf{S}(\boldsymbol{p})
\mathbf{V}. \$\$ The first-order delta method then yields \$\$
\mathrm{Var}(\hat{\boldsymbol{p}}) \approx
\mathbf{J}\_{\boldsymbol{\psi}} I\_{\boldsymbol{z}}^{-1}
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}. \$\$ This simplex
covariance is invariant to orthogonal rotations of \\\mathbf{V}\\. The
construction is undefined on the simplex boundary (the log-ratio chart
blows up); the function then returns `NA` with a warning.

Wald intervals in
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
and
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
all use this expected-Fisher ILR map. The observed Hessian
[`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md)
is for Newton steps only; it is not the Wald covariance.

## Usage

``` r
vcov_ilr_delta(p, mean_signature_matrix, Sigma, warn = TRUE)
```

## Arguments

- p:

  Numeric proportions on the open simplex.

- mean_signature_matrix:

  Mean signature \\\boldsymbol{\mu}\\ (\\G\times J\\).

- Sigma:

  Cell-type covariances \\G\times G\times J\\.

- warn:

  If `FALSE`, skip the boundary / singularity warnings (Monte Carlo
  loops). Default `TRUE`.

## Value

Symmetric \\J\times J\\ asymptotic covariance of
\\\hat{\boldsymbol{p}}\\, or a matrix of `NA` if the bound is undefined
/ singular.

## See also

[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md),
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)
