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
\mathcal{N}(\boldsymbol{\rho}\_{0}, I\_{\boldsymbol{\rho}}^{-1})\\
(Cramer–Rao / asymptotic normality of MLEs; law of large numbers for the
score). The first-order delta method then yields the simplex covariance
used by
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
and the Wald standard errors in
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md):
\$\$ \mathrm{Var}(\hat{\boldsymbol{p}}) \approx
\mathbf{J}\_{\boldsymbol{\psi}} I\_{\boldsymbol{\rho}}^{-1}
\mathbf{J}\_{\boldsymbol{\psi}}^{\top}. \$\$ Diagonal square roots of
this matrix are the asymptotic standard errors; Wald intervals at level
\\1-\alpha\\ are \\\hat{p}\_j \pm z\_{1-\alpha/2}\\\mathrm{SE}\_j\\ with
\\z\_{q}=\Phi^{-1}(q)\\. The construction is undefined on the simplex
boundary (ALR chart blows up); the function then returns `NA` with a
warning. See also <https://en.wikipedia.org/wiki/Delta_method> and
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

[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md),
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
[`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md)
