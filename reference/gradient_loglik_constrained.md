# Constrained gradient via the chain rule

Returns \$\$ \nabla\_{\boldsymbol{\rho}}\ell =
\bigl(\nabla\_{\boldsymbol{p}}\ell\bigr)^{\mathsf{T}}
\mathbf{J}\_{\boldsymbol{\psi}}(\boldsymbol{\rho}), \$\$ i.e.
first-order chain rule for \\\ell\circ\boldsymbol{\psi}\\.

## Usage

``` r
gradient_loglik_constrained(rho, y, mean_signature_matrix, Sigma)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Numeric vector in \\\mathbb{R}^{J-1}\\.
