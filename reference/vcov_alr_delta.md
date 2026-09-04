# Cramer–Rao / ALR delta-method covariance of \\\hat{\boldsymbol{p}}\\

Maps the expected Fisher information of unconstrained proportions
through the additive log-ratio (ALR) chart and back to the simplex via
the delta method.

Let \\\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})\\ with
Jacobian \\\mathbf{J}\_{\boldsymbol{\psi}}
=\partial\boldsymbol{\psi}/\partial\boldsymbol{\rho}^{\top}\\
([`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md)).
Fisher information transforms as the covariant quadratic form \$\$
I\_{\boldsymbol{\rho}} = \mathbf{J}\_{\boldsymbol{\psi}}^{\top}
I(\boldsymbol{p}) \mathbf{J}\_{\boldsymbol{\psi}}. \$\$ Under a regular
large-sample regime the MLE in ALR coordinates is asymptotically normal,
\\\hat{\boldsymbol{\rho}} \overset{a}{\sim}
\mathcal{N}(\boldsymbol{\rho}\_{0}, I\_{\boldsymbol{\rho}}^{-1})\\. The
first-order delta method then yields the same simplex covariance as
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
when both charts are transformed correctly: \$\$
\mathrm{Var}(\hat{\boldsymbol{p}}) \approx
\mathbf{J}\_{\boldsymbol{\psi}} I\_{\boldsymbol{\rho}}^{-1}
\mathbf{J}\_{\boldsymbol{\psi}}^{\top}. \$\$ This helper is kept for
reference-invariance checks;
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
uses
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md).
The construction is undefined on the simplex boundary (the ALR chart
blows up); the function then returns `NA` with a warning. See also
<https://en.wikipedia.org/wiki/Delta_method> and
<https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution>.

## Usage

``` r
vcov_alr_delta(p, mean_signature_matrix, Sigma)
```

## Arguments

- p:

  Numeric proportions on the open simplex.

- mean_signature_matrix:

  Mean signature \\\boldsymbol{\mu}\\ (\\G\times J\\).

- Sigma:

  Cell-type covariances \\G\times G\times J\\.

## Value

Symmetric \\J\times J\\ asymptotic covariance of
\\\hat{\boldsymbol{p}}\\, or a matrix of `NA` if the bound is undefined
/ singular.

## See also

[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md),
[`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md)
