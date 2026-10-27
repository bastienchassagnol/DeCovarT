# Simulate an independent reference panel

Draws, for every cell type \\j\\, `n_individuals` \\\times\\ `n_cells`
profiles from \\\mathcal{N}\_G(\boldsymbol{\mu}\_{\cdot
j},\boldsymbol{\Sigma}\_j)\\ and labels each with its individual. With
`n_cells = 1` every column is the purified profile of one individual, as
in the population convolution; with `n_cells > 1` every individual
contributes `n_cells` cells, as in a single-cell reference. The
generative models of this package have no between-individual layer, so
individuals are exchangeable labels.

Benchmark methods that estimate the reference moments (MEAD averages
per-individual pseudo-bulks, RNA-Sieve pools the cells) receive this
panel; methods that take plug-in moments receive the true
\\(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_j\\)\\. Call it with its own
`seed`, distinct from the seed of the bulk samples, so that no bulk is
built from a profile of the panel.

## Usage

``` r
simulate_reference_replicates(
  mean_signature_matrix,
  Sigma,
  n_individuals,
  n_cells = 1L,
  seed = NULL
)
```

## Arguments

- mean_signature_matrix:

  \\G\times J\\ mean matrix with gene row names and cell-type column
  names.

- Sigma:

  \\G\times G\times J\\ covariance array.

- n_individuals:

  Number of individuals \\R\ge2\\.

- n_cells:

  Profiles per individual and cell type.

- seed:

  Optional seed; the global RNG state is restored afterwards.

## Value

A list with `reference`, a list named by cell type of
\\G\times(R\\n\_{\mathrm{cells}})\\ matrices, and `individual`, a list
named by cell type of integer vectors giving the individual of each
column.

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md)

## Examples

``` r
mu <- matrix(
  c(5, 1, 1, 1, 5, 1),
  nrow = 3,
  dimnames = list(paste0("g", 1:3), c("ct1", "ct2"))
)
sigma <- array(
  c(diag(3), diag(3)),
  dim = c(3, 3, 2),
  dimnames = list(rownames(mu), rownames(mu), colnames(mu))
)
panel <- simulate_reference_replicates(mu, sigma, n_individuals = 4, seed = 1)
dim(panel$reference$ct1)
#> [1] 3 4
panel$individual$ct1
#> [1] 1 2 3 4
```
