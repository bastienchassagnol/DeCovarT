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

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
gradient_loglik_constrained(additive_log_ratio(p), y, mu, Sigma)
#>            [,1]
#> [1,] -0.1846154
```
