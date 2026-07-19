# Maximum a posteriori purified profiles under a Gaussian convolution

For fixed plug-in
\\\boldsymbol{\zeta}=(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_j\\)\\ and
a bulk vector \\\boldsymbol{y}\\, returns cell-type-specific MAP
estimates of the latent purified profiles in the additive Gaussian model
of the article.

## Usage

``` r
.map_gaussian_convolution(y, mean_signature_matrix, Sigma)
```

## Arguments

- y:

  Bulk vector \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Mean matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\.

- Sigma:

  Array \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times
  J}\\.

## Value

List of length \\J\\ with MAP vectors in \\\mathbb{R}^{G}\\.
