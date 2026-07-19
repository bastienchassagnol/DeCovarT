# Constrained log-likelihood \\\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{\psi}(\boldsymbol{\rho}))\\

Composes
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
with
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md),
so that optimisation may be performed over
\\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

## Usage

``` r
loglik_multivariate_constrained(rho, y, mean_signature_matrix, Sigma)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Scalar log-likelihood on the constrained manifold.
