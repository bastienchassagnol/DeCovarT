# Plot deconvolution metric heatmaps

For each algorithm, visualises a selected score over the design grid of
simulated \\\boldsymbol{p}\\ (and related scenario factors).

## Usage

``` r
plot_correlation_Heatmap(
  distribution_metrics,
  score_variable = "model_mse",
  n_break = 20,
  uni_scale = TRUE
)
```

## Arguments

- distribution_metrics:

  Tibble of metric scores from a benchmark.

- score_variable:

  Column name of the metric to display.

- n_break:

  Number of colour breaks.

- uni_scale:

  If `FALSE`, each panel uses its own colour scale.

## Value

A named list of `ComplexHeatmap` heatmap objects (one per algorithm).

## Examples

``` r
metrics <- tibble::tibble(
  correlation_celltype1 = c(0, 0, 0.5, 0.5),
  correlation_celltype2 = c(0, 0.5, 0, 0.5),
  algorithm = "nnls",
  model_mse = c(0.01, 0.02, 0.015, 0.03)
)
ht <- plot_correlation_Heatmap(metrics, score_variable = "model_mse")
names(ht)
#> [1] "nnls"
```
