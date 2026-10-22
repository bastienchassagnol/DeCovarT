# Simulate a bulk sample as a sum of single cells

Draws cell counts \\m\_{\cdot i}\sim Multinomial(n_i,p\_{\cdot i})\\,
then the exact Gaussian sum of those cells, \\T\_{ji}\sim
N(m\_{ji}\mu\_{\cdot j},m\_{ji}\Sigma\_{W,j})\\, and scales the total by
the platform factor: \\\boldsymbol{y}\_{\cdot
i}=a_i\sum_j\boldsymbol{T}\_{ji}\\. A sum of \\m\\ independent Gaussian
cells is Gaussian, so there is no per-cell loop.

This is the single-cell mixture (covariance linear in the cell counts).
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
remains the population-profile convolution, whose covariance is weighted
by \\p_j^2\\.

Pass the cellular moments (`generative_mean`, `generative_sigma` from
[`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md)),
not a library-normalised reference. The returned estimand is the
realised cell fraction \\\boldsymbol{m}\_{\cdot i}/n_i\\.

## Usage

``` r
simulate_bulk_from_cells(
  mean_signature_matrix,
  Sigma,
  p,
  n_cells,
  platform_scale = 1,
  n_samples = NULL,
  seed = NULL
)
```

## Arguments

- mean_signature_matrix:

  Per-cell means \\\tilde{\boldsymbol{\mu}}\\ (\\G\times J\\).

- Sigma:

  Per-cell covariances \\\boldsymbol{\Sigma}\_{W,j}\\ (\\G\times G\times
  J\\).

- p:

  Cell-fraction vector of length \\J\\, or a \\J\times N\\ matrix. Each
  column must sum to 1.

- n_cells:

  Cells per bulk sample, \\n_i\\. Recycled to \\N\\. The multinomial
  size is the rounded value; the fitted \\\hat n\\ need not be an
  integer.

- platform_scale:

  Capture scale \\a_i\\. Use `1` when the bulk and the single-cell
  platform share their scale, and `2` when the bulk captures twice as
  much transcript per cell. Recycled to \\N\\.

- n_samples:

  Number of bulk samples when `p` is a vector. Ignored when `p` already
  has one column per sample, unless it matches that number.

- seed:

  Optional seed, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the draws.

## Value

A list with `Y` (\\G\times N\\), nominal `p`, realised `p_realised`
(\\\boldsymbol{m}/n\\), `counts` (\\\boldsymbol{m}\\), transcript
fractions `q`, `cell_size` \\S_j=\mathbf{1}^{\top}
\tilde{\boldsymbol{\mu}}\_{\cdot j}\\, `platform_scale` \\a\\, `n_cells`
\\n\\, `mean_scale` \\c=an\\ and `variance_scale` \\\phi=a^{2}n\\.

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md),
[`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)

## Examples

``` r
genes <- paste0("g", 1:4)
cts <- paste0("ct", 1:2)
mu <- matrix(
  c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
  nrow = 4,
  dimnames = list(genes, cts)
)
sigma <- array(0, c(4, 4, 2), dimnames = list(genes, genes, cts))
sigma[,, 1] <- diag(0.05, 4)
sigma[,, 2] <- diag(0.05, 4)
set.seed(1)
sim <- simulate_bulk_from_cells(
  mu,
  sigma,
  p = c(0.8, 0.2),
  n_cells = 2000,
  platform_scale = 2,
  n_samples = 2
)
colSums(sim$counts)
#> sample_1 sample_2 
#>     2000     2000 
```
