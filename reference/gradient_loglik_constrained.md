# Constrained gradient via the chain rule

Returns \$\$ \nabla\_{\boldsymbol{z}}\ell =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\nabla\_{\boldsymbol{p}}\ell = \mathbf{V}^{\mathsf{T}}
\mathbf{S}(\boldsymbol{p}) \nabla\_{\boldsymbol{p}}\ell, \$\$ i.e. the
first-order chain rule for \\\ell\circ\boldsymbol{\psi}\\ in ILR
coordinates.

## Usage

``` r
gradient_loglik_constrained(z, y, mean_signature_matrix, Sigma, V = NULL)
```

## Arguments

- z:

  Numeric vector \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- V:

  Optional ILR basis; see
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md).

## Value

Numeric vector in \\\mathbb{R}^{J-1}\\.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
gradient_loglik_constrained(isometric_log_ratio(p), y, mu, Sigma)
#> [1] -0.2610856
```
