# Hessian \\\mathbf{H}\\ of the unconstrained log-likelihood

Analytic Hessian \\\mathbf{H}\in\mathcal{M}\_{J\times J}\\ with entries
\\\mathbf{H}\_{i,j}=\partial^{2}\ell/(\partial p_i\partial p_j)\\,
matching the matrix formulae of the article (quadratic forms in
\\\boldsymbol{\Theta}\\, \\\boldsymbol{\Sigma}\_i\\,
\\\boldsymbol{\mu}\_{\cdot i}\\ and residual
\\\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}\\).

## Usage

``` r
hessian_loglik_unconstrained(p, y, mean_signature_matrix, Sigma)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

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

Symmetric numeric matrix \\\mathbf{H}\\.

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
hessian_loglik_unconstrained(p, y, mu, Sigma)
#>           [,1]      [,2]
#> [1,] -1697.041 -1685.207
#> [2,] -1685.207 -1702.959
```
