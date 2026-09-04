# Plot deconvolution metric heatmaps

For each algorithm, visualises a selected score over the design grid of
simulated \\\boldsymbol{p}\\ (and related scenario factors). Requires
optional Suggests packages `ComplexHeatmap`, `circlize`, and `viridis`
(install `ComplexHeatmap` via Bioconductor). Use this helper only for
**linked, annotated, hierarchical** grids. For algorithm-similarity
correlation matrices prefer
[`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md)
([`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)).

## Usage

``` r
plot_correlation_Heatmap(
  distribution_metrics,
  score_variable = "model_mse",
  n_break = 20,
  uni_scale = TRUE,
  file = NULL
)
```

## Arguments

- distribution_metrics:

  Tibble of metric scores from a benchmark.

- score_variable:

  Column name of the metric to display (`"model_mse"`, `"model_rmse"`,
  `"model_coef_determination"`, `"model_coef_determination_adjusted"`,
  `"model_mae"`, `"model_cor"`, `"model_ccc"`). Matching is
  case-insensitive.

- n_break:

  Number of colour breaks.

- uni_scale:

  If `FALSE`, each panel uses its own colour scale.

- file:

  Optional PDF path. When supplied, heatmaps are drawn with
  [`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html); a missing
  `.pdf` suffix is added (G4.0).

## Value

A named list of `ComplexHeatmap` heatmap objects (one per algorithm).

## See also

[`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md),
[`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md)

## Examples

``` r
metrics <- tibble::tibble(
  correlation_celltype1 = c(0, 0, 0.5, 0.5),
  correlation_celltype2 = c(0, 0.5, 0, 0.5),
  algorithm = "nnls",
  model_mse = c(0.01, 0.02, 0.015, 0.03)
)
if (
  requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE) &&
    requireNamespace("viridis", quietly = TRUE)
) {
  ht <- plot_correlation_Heatmap(metrics, score_variable = "model_mse")
  names(ht)
}
#> [1] "nnls"
```
