# Tile heatmap of a scenario metric on the inner 2-by-2 design

Fig02: one page per meta-scenario (composition \\\times\\ variance
\\\times\\ CLD) with tiles on \\(\rho_1,\rho_2)\\ and one panel per
solver. Fig03: tiles on the CT1 / CT2 graph families.

## Usage

``` r
plot_bivariate_metric_tiles(metrics, title)
```

## Arguments

- metrics:

  Long table with `algorithm`, design columns, and `value`.

- title:

  Page title.

## Value

A `ggplot`.
