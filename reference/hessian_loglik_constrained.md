# Constrained Hessian of \\\ell\circ\boldsymbol{\psi}\\

Second-order chain rule \$\$ \mathbf{H}\_{\boldsymbol{\rho}} =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\mathbf{H}\_{\boldsymbol{p}} \mathbf{J}\_{\boldsymbol{\psi}}
+\sum\_{i=1}^{J} \frac{\partial\ell}{\partial p_i}\\
\frac{\partial^{2}p_i}{\partial\boldsymbol{\rho}\partial\boldsymbol{\rho}^{\mathsf{T}}}.
\$\$

## Usage

``` r
hessian_loglik_constrained(rho, y, mean_signature_matrix, Sigma)
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

Symmetric matrix in \\\mathcal{M}\_{(J-1)\times(J-1)}\\.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
hessian_loglik_constrained(additive_log_ratio(p), y, mu, Sigma)
#>           [,1]
#> [1,] -1.258225
```
