# Cell-size conventions for a single-cell reference

Turns nonnegative gene-by-cell-type moments into the reference used by
the single-cell mixture. Column \\j\\ of the supplied mean is normalised
to a composition \\\boldsymbol{\theta}\_{\cdot j}\\. Two conventions are
available; the default is the raw-covariance, cell-fraction model of the
article.

- `cellular_raw` (default). Mean \\\tilde{\boldsymbol{\mu}}\_{\cdot
  j}=S_j\boldsymbol{\theta}\_{\cdot j}\\ with the covariance left as
  supplied. Keeping \\S_j\\ in the mean makes the fitted composition a
  cell fraction, not a transcript fraction.

- `transcript_raw`. Library-normalised mean
  \\\boldsymbol{\theta}\_{\cdot j}\\ and per-transcript covariance
  \\\boldsymbol{\Sigma}\_j/S_j\\. The fitted composition is the
  transcript fraction \\\boldsymbol{q}\\. The physical generator is
  still the cellular model: `generative_mean` and `generative_sigma` are
  the moments
  [`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md)
  should draw from.

## Usage

``` r
rescale_single_cell_moments(
  moments,
  cell_size,
  reference = c("cellular_raw", "transcript_raw")
)
```

## Arguments

- moments:

  List with `mean_profiles` (\\G\times J\\) and `covariance_matrices`
  (\\G\times G\times J\\), as returned by
  [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
  with `nonnegative = TRUE`.

- cell_size:

  Positive cell sizes \\S_j\\. A scalar is repeated once per cell type.
  For two types the scenarios use \\(1, 1)\\ and \\(1, 5)\\; \\(1, 2,
  5)\\ is the same idea for three types.

- reference:

  `cellular_raw` or `transcript_raw`.

## Value

A list with the matrices passed to the estimator
(`mean_signature_matrix`, `Sigma`), the cellular generator
(`generative_mean`, `generative_sigma`), `cell_size`, `composition`
(\\\boldsymbol{\theta}\\), `reference` and `estimand` (`"cellular"` or
`"transcript"`).

## See also

[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md),
[`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md),
[`single_cell_scenario_design()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_scenario_design.md)

## Examples

``` r
mu <- cbind(c(2, 2, 0, 0), c(0, 0, 1, 1))
sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
moments <- list(mean_profiles = mu, covariance_matrices = sigma)
ref <- rescale_single_cell_moments(moments, c(1, 5))
colSums(ref$generative_mean)
#> [1] 1 5
```
