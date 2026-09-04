# Constrained log-likelihood \\\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{\psi}(\boldsymbol{z}))\\

Composes
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
with
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md),
so that optimisation may be performed over unconstrained ILR coordinates
\\\boldsymbol{z}\in\mathbb{R}^{J-1}\\.

## Usage

``` r
loglik_multivariate_constrained(z, y, mean_signature_matrix, Sigma, V = NULL)
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

Scalar log-likelihood on the constrained manifold.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
loglik_multivariate_constrained(isometric_log_ratio(p), y, mu, Sigma)
#> [1] 0.6539265
```
