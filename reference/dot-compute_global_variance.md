# Bulk mixture covariance \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Assembles the conditional covariance of the Gaussian convolution \$\$
\boldsymbol{\Sigma}(\boldsymbol{p})
=\sum\_{j=1}^{J}p_j^{2}\\\boldsymbol{\Sigma}\_j, \$\$ stored as slices
of the array `Sigma`.

## Usage

``` r
.compute_global_variance(p, Sigma)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Symmetric matrix
\\\boldsymbol{\Sigma}(\boldsymbol{p})\in\mathcal{M}\_{G\times G}\\.

## Examples

``` r
p <- c(0.6, 0.4)
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
.compute_global_variance(p, Sigma)
#>      [,1] [,2]
#> [1,] 0.52 0.00
#> [2,] 0.00 0.52
```
