# Horizontal raincloud of Monte Carlo proportion estimates

Half-eye densities, dots, and empirical 50% / 95% intervals (Allen et
al. 2019) for the Monte Carlo sampling distribution of \\\hat p_j\\.
Default aesthetics: numeric axis is the estimation error \\\hat
p_j-p_j^{\star}\\ (so a vertical line at 0 is bias); the y-axis is the
cell type; colour and dodge group the deconvolution algorithm. Optional
`facet_rows` / `facet_cols` split scenarios (for example number of genes
versus pairwise cosine).

## Usage

``` r
plot_mc_raincloud(
  benchmark,
  quantity = c("error", "estimate"),
  facet_rows = NULL,
  facet_cols = NULL,
  .width = c(0.5, 0.95)
)
```

## Arguments

- benchmark:

  List from
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md),
  or a table from
  [`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md).

- quantity:

  `"error"` (default) or `"estimate"`.

- facet_rows, facet_cols:

  Optional column names for
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  rows and columns.

- .width:

  Passed to
  [`ggdist::stat_halfeye()`](https://mjskay.github.io/ggdist/reference/stat_halfeye.html);
  default `c(0.5, 0.95)`.

## Value

A `ggplot` object.

## Details

The inner interval is the central 50% of Monte Carlo replicates; the
outer interval is the central 95%. These are **not** confidence
intervals for \\p_j\\.

## References

Allen M, Poggiali D, Whitaker K, Marshall TR, van Langen J, Kievit RA
(2019). “Raincloud Plots: A Multi-Platform Tool for Robust Data
Visualization.” *Wellcome Open Research*, **4**, 63. ISSN 2398-502X.
[doi:10.12688/wellcomeopenres.15191.2](https://doi.org/10.12688/wellcomeopenres.15191.2)
. <https://wellcomeopenresearch.org/articles/4-63>.

## See also

[`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md),
[`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md)

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
out <- run_simulation_benchmark(
  tibble::tibble(
    n_genes = 2L,
    cosine = 0,
    true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
  ),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  n = 4,
  cores = 1
)
if (requireNamespace("ggdist", quietly = TRUE)) {
  plot_mc_raincloud(out, facet_cols = "cosine")
}
```
