# Log-likelihood of one single-cell bulk sample

Gaussian log-likelihood of
\\\boldsymbol{y}\sim\mathcal{N}(c\tilde{\boldsymbol{\mu}}\boldsymbol{p},
\phi\boldsymbol{S}\_W(\boldsymbol{p}))\\ with
\\\boldsymbol{S}\_W(\boldsymbol{p})=\sum_j
p_j\boldsymbol{\Sigma}\_{W,j}\\. The scales \\c\\ and \\\phi\\ are held
fixed. Up to \\-\tfrac{G}{2}\log(2\pi)\\, \$\$ \ell =
-\tfrac{G}{2}\log\phi
-\tfrac{1}{2}\log\det\boldsymbol{S}\_W(\boldsymbol{p}) -\tfrac{1}{2\phi}
\boldsymbol{r}^{\top}\boldsymbol{S}\_W(\boldsymbol{p})^{-1}\boldsymbol{r},
\$\$
\\\boldsymbol{r}=\boldsymbol{y}-c\tilde{\boldsymbol{\mu}}\boldsymbol{p}\\.
The population objective
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
is a different model: its covariance is \\\sum_j
p_j^2\boldsymbol{\Sigma}\_j\\.

## Usage

``` r
loglik_single_cell(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
)
```

## Arguments

- p:

  Composition of length \\J\\. It does not have to sum to 1 (numerical
  derivatives step off the simplex).

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- mean_scale:

  Positive plug-in \\c=an\\.

- variance_scale:

  Positive plug-in \\\phi=a^2 n\\.

## Value

Scalar log-likelihood.

## Examples

``` r
mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
p <- c(0.7, 0.3)
y <- 20 * drop(mu %*% p)
loglik_single_cell(p, y, mu, sigma, 20, 4)
#> [1] 3.218876
```
