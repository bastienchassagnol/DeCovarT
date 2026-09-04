# Constrained Hessian of \\\ell\circ\boldsymbol{\psi}\\

Second-order chain rule in ILR coordinates \$\$
\mathbf{H}\_{\boldsymbol{z}} =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\mathbf{H}\_{\boldsymbol{p}} \mathbf{J}\_{\boldsymbol{\psi}}
+\sum\_{i=1}^{J} \frac{\partial\ell}{\partial p_i}\\
\mathbf{H}\_{\psi_i}. \$\$ The second summand cannot be dropped away
from stationarity. At an interior KKT point
\\\nabla_p\ell=\lambda\mathbf{1}\\ it vanishes because
\\\sum_i\mathbf{H}\_{\psi_i}=\mathbf{0}\\.

## Usage

``` r
hessian_loglik_constrained(z, y, mean_signature_matrix, Sigma, V = NULL)
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

Symmetric matrix in \\\mathcal{M}\_{(J-1)\times(J-1)}\\.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
hessian_loglik_constrained(isometric_log_ratio(p), y, mu, Sigma)
#>          [,1]
#> [1,] -2.51645
```
