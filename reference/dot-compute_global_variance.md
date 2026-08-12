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
