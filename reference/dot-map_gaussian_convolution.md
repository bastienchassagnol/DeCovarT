# Maximum a posteriori purified profiles under a Gaussian convolution

For fixed plug-in
\\\boldsymbol{\zeta}=(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_j\\)\\ and
a bulk vector \\\boldsymbol{y}\\, returns cell-type-specific MAP
estimates of the **latent** purified profiles \\\boldsymbol{x}\_{\cdot
j}\\ in the additive Gaussian model of the article. This is the Bayesian
counterpart to the frequentist plug-in that replaces
\\\boldsymbol{x}\_{\cdot j}\\ by \\\boldsymbol{\mu}\_{\cdot j}\\ when
estimating proportions alone.

## Usage

``` r
.map_gaussian_convolution(y, mean_signature_matrix, Sigma)
```

## Arguments

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

List of length \\J\\ with MAP vectors in \\\mathbb{R}^{G}\\.

## Examples

``` r
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
y <- drop(mu %*% c(0.6, 0.4))
.map_gaussian_convolution(y, mu, Sigma)
#> [[1]]
#>     ct1
#> g1  9.4
#> g2 11.6
#> 
#> [[2]]
#>     ct2
#> g1 11.4
#> g2  9.6
#> 
```
