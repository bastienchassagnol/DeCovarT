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

## Details

The log-determinant contributes
\\-\delta\_{ij}\mathrm{Tr}(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_j)
+2p_ip_j\mathrm{Tr}(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_i
\boldsymbol{\Theta}\boldsymbol{\Sigma}\_j)\\, i.e. half the coefficients
of the pre-2.3.0 objective, which used
\\-\log\det\boldsymbol{\Sigma}(\boldsymbol{p})\\. Residual terms are
unchanged. Taking expectations under
\\\boldsymbol{r}\sim\mathcal{N}\_G(\boldsymbol{0},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\ cancels the determinant and
residual traces, leaving exactly
\\\mathbb{E}\[-\mathbf{H}\]=I(\boldsymbol{p})\\ of
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md).

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
p <- c(0.6, 0.4)
y <- drop(mu %*% p)
hessian_loglik_unconstrained(p, y, mu, Sigma)
#>           [,1]      [,2]
#> [1,] -1698.521 -1688.757
#> [2,] -1688.757 -1701.479
```
