# Simulate bulk mixtures from a multivariate Gaussian convolution

For each bootstrap sample \\i=1,\ldots,N\\, draws **latent** purified
profiles \\\boldsymbol{x}\_{\cdot j}^{(i)}\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\_{\cdot j}, \boldsymbol{\Sigma}\_j)\\
independently for each cell type \\j=1,\ldots,J\\, then forms the bulk
by the linear convolution \$\$ \boldsymbol{y}\_{\cdot i}
=\boldsymbol{X}^{(i)}\boldsymbol{p}
=\sum\_{j=1}^{J}p_j\\\boldsymbol{x}\_{\cdot j}^{(i)}, \$\$ matching the
article's conditional model
\\\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p}))\\ with
\\\boldsymbol{\Sigma}(\boldsymbol{p})= \sum_j
p_j^{2}\boldsymbol{\Sigma}\_j\\.

The mean signature \\\boldsymbol{\mu}\\ is shared across samples; only
the latent draws \\\boldsymbol{x}\_{\cdot j}^{(i)}\\ vary with \\i\\.
Stacking those draws into the three-way array
\\\mathcal{X}=(x\_{gji})\in\mathcal{M}\_{G\times J\times N}\\, the bulk
matrix is the mode-2 tensor–vector contraction \$\$ \boldsymbol{Y}
=\mathcal{X}\times\_{2}\boldsymbol{p}, \qquad
y\_{gi}=\sum\_{j=1}^{J}x\_{gji}\\p\_{j} \quad(g=1,\ldots,G;\\
i=1,\ldots,N), \$\$ which for each sample recovers
\\\boldsymbol{y}\_{\cdot i}=\boldsymbol{X}^{(i)}\\\boldsymbol{p}\\ with
\\\boldsymbol{X}^{(i)}\in\mathcal{M}\_{G\times J}\\.

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

  Mean signature \\\boldsymbol{\mu}=(\mu\_{gj})\in\mathcal{M}\_{G\times
  J}\\ (shared across samples; not the latent profiles).

- Sigma:

  Array of covariances
  \\(\boldsymbol{\Sigma}\_j)\_{j}\in\mathcal{M}\_{G\times G\times J}\\.

- p:

  Proportion vector \\\boldsymbol{p}\in\Delta^{J-1}\\ (default:
  uniform).

- n:

  Number of bulk / bootstrap samples \\N\\.

## Value

A list with:

- `latent_profiles`: array
  \\\mathcal{X}=(x\_{gji})\in\mathcal{M}\_{G\times J\times N}\\ of
  **unobserved** cell-type-specific draws (one \\G\times J\\ slice per
  sample \\i\\);

- `Y`: matrix \\\boldsymbol{Y}\in\mathcal{M}\_{G\times N}\\ whose
  columns are bulk vectors \\\boldsymbol{y}\_{\cdot i}\\, obtained as
  \\\mathcal{X}\times\_{2}\boldsymbol{p}\\.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
sim <- simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 5)
dim(sim$Y)
#> [1] 2 5
```
