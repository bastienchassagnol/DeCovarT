# Simulate bulk mixtures from a Gaussian convolution

Draws purified profiles
\\\boldsymbol{x}\_j\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\_{\cdot j},
\boldsymbol{\Sigma}\_j)\\ independently for each cell type and forms
bulk observations by the linear mixture \$\$ \boldsymbol{y}\_{\cdot i}
=\boldsymbol{\mu}^{(i)}\boldsymbol{p}
=\sum\_{j=1}^{J}p_j\\\boldsymbol{x}\_j^{(i)}, \$\$ matching the
article's conditional model
\\\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\ with
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2}\boldsymbol{\Sigma}\_j\\.

## Usage

``` r
simulate_bulk_mixture(
  signature_matrix,
  Sigma,
  p = rep(1/ncol(signature_matrix), ncol(signature_matrix)),
  n = 500
)
```

## Arguments

- signature_matrix:

  Mean matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\.

- Sigma:

  Array of covariances
  \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\.

- p:

  Proportion vector \\\boldsymbol{p}\in\Delta^{J-1}\\ (default:
  uniform).

- n:

  Number of bulk samples \\n\\.

## Value

A list with:

- `mean_signature_matrix`: array \\(x\_{gji})\in\mathcal{M}\_{G\times
  J\times N}\\ of simulated purified profiles;

- `Y`: matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\ whose
  columns are bulk vectors \\\boldsymbol{y}\_{\cdot i}\\.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`benchmark_bivariate_gaussian_convolutions()`](https://bastienchassagnol.github.io/DeCovarT/reference/benchmark_bivariate_gaussian_convolutions.md)
