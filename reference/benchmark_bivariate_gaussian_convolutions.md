# Benchmark bivariate Gaussian convolutions

Wrapper reproducing the bivariate (\\G=2\\, \\J=2\\) toy study of the
article: for each scenario it builds \\\boldsymbol{\mu}\\ and
\\(\boldsymbol{\Sigma}\_j)\_{j}\\, simulates \\\boldsymbol{Y}\\ via
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
and deconvolves with the supplied algorithms. Performance is summarised
against entropy of \\\boldsymbol{p}\\ and overlap of the Gaussian
mixture.

## Usage

``` r
benchmark_bivariate_gaussian_convolutions(
  proportions = list(balanced = c(0.5, 0.5), `small unbalanced` = c(0.6, 0.4),
    `highly unbalanced` = c(0.05, 0.95)),
  signature_matrices = list(`small OVL` = matrix(c(20, 40, 40, 20), nrow = 2)),
  corr_sequence = seq(-0.8, 0.8, 0.2),
  diagonal_terms = list(homoscedastic = c(1, 1), heteroscedastic = c(1, 2)),
  deconvolution_functions = list(lsfit = list(FUN = deconvolute_ratios_lsfit,
    additional_parameters = NULL)),
  n = 200,
  scaled = FALSE,
  cores = ifelse(.Platform$OS.type == "unix", getOption("mc.cores",
    parallel::detectCores()), 1)
)
```

## Arguments

- proportions:

  List of simplex vectors \\\boldsymbol{p}\\.

- signature_matrices:

  List of mean matrices \\\boldsymbol{\mu}\in\mathcal{M}\_{2\times
  2}^{+}\\.

- corr_sequence, diagonal_terms:

  Correlation sequence and diagonal variance templates used to assemble
  \\\boldsymbol{\Sigma}\_j=
  \mathrm{D}\_{j}^{1/2}\mathbf{R}\_j\mathrm{D}\_{j}^{1/2}\\.

- deconvolution_functions:

  Named list of deconvolution callables (each with `FUN` and optional
  `additional_parameters`).

- n:

  Number of bulk replicates \\N\\ per scenario.

- scaled:

  Logical; whether to log-scale inputs before estimation.

- cores:

  Parallel workers for
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

## Value

A list with `config` (design + entropy/overlap) and `simulations`
(estimation tibble).

## Details

Designed for two cell types and two genes. Larger \\(G,J)\\ with only
bivariate observations is prone to non-identifiability.

Scenarios are enumerated with
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
and tagged with a unique `ID` via
[`dplyr::row_number()`](https://dplyr.tidyverse.org/reference/row_number.html),
so no side-effect mutation is needed while looping over the design.

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)

## Examples

``` r
set.seed(1)
out <- benchmark_bivariate_gaussian_convolutions(
  proportions = list("balanced" = c(0.5, 0.5)),
  signature_matrices = list("small" = matrix(c(20, 22, 22, 20), 2)),
  corr_sequence = 0,
  diagonal_terms = list("homoscedastic" = c(1, 1)),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  n = 2,
  cores = 1
)
#> Scenario B1_Ho: balanced, corr=(0, 0), centroids=small, variance=homoscedastic.
nrow(out$simulations)
#> [1] 2
```
