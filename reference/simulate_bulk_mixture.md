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
\\\boldsymbol{\Sigma}(\boldsymbol{p})= \sum_j
p_j^{2}\boldsymbol{\Sigma}\_j\\.

Equivalently, stacking the purified draws into the three-way array
\\\mathcal{X}=(x\_{gji})\in\mathcal{M}\_{G\times J\times N}\\, the bulk
matrix is the mode-2 tensor–vector contraction \$\$ \boldsymbol{Y}
=\mathcal{X}\times\_{2}\boldsymbol{p}, \qquad
y\_{gi}=\sum\_{j=1}^{J}x\_{gji}\\p\_{j} \quad(g=1,\ldots,G;\\
i=1,\ldots,N), \$\$ which for each sample recovers the matrix–vector
product \\\boldsymbol{y}\_{\cdot i}=\boldsymbol{X}\_{\cdot\cdot
i}\\\boldsymbol{p}\\ with \\\boldsymbol{X}\_{\cdot\cdot
i}\in\mathcal{M}\_{G\times J}\\.

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

- `mean_signature_matrix`: array
  \\\mathcal{X}=(x\_{gji})\in\mathcal{M}\_{G\times J\times N}\\ of
  simulated purified profiles;

- `Y`: matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\ whose
  columns are bulk vectors \\\boldsymbol{y}\_{\cdot i}\\, obtained as
  \\\mathcal{X}\times\_{2}\boldsymbol{p}\\.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`benchmark_bivariate_gaussian_convolutions()`](https://bastienchassagnol.github.io/DeCovarT/reference/benchmark_bivariate_gaussian_convolutions.md)
