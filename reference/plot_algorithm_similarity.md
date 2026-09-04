# Tile heatmap of algorithm-similarity correlations

[`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
display of
[`algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/algorithm_similarity.md),
with rows and columns ordered by average-linkage clustering of \\1-r\\.
Optional dendrogram via `ggdendro` (Suggests). This is the default for a
small correlation matrix;
[`plot_correlation_Heatmap()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_correlation_Heatmap.md)
is reserved for linked multi-omics grids.

## Usage

``` r
plot_algorithm_similarity(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL,
  dendrogram = FALSE
)
```

## Arguments

- benchmark:

  List from
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md).

- facet_rows, facet_cols:

  Optional scenario columns. When supplied, correlations are computed
  inside each facet cell.

- dendrogram:

  If `TRUE`, attach a `ggdendro` ggplot as attribute `"dendrogram"`
  (ignored when scenario facets are used).

## Value

A `ggplot` object.

## See also

[`algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/algorithm_similarity.md),
[`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md)

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
    true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
  ),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls),
    "rlm" = list(FUN = deconvolute_ratios_rlm)
  ),
  n = 4,
  cores = 1
)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_algorithm_similarity(out)
}
```
