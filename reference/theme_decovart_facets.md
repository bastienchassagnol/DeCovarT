# Faceted ggplot2 theme (black strips, panel border)

`theme_minimal()` plus a panel border and white-on-black facet strips,
after `theme_features()` in the
[atlas-feature-selection-benchmark](https://github.com/theislab/atlas-feature-selection-benchmark/blob/b89fc0f66747062e6e1b4b35bd392b27ad035295/analysis/R/plotting.R)
plotting helpers.

## Usage

``` r
theme_decovart_facets(base_size = 11, ...)
```

## Arguments

- base_size:

  Base font size for
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- ...:

  Passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

## Value

A `ggplot2` theme object.

## See also

[`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md),
[`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md),
[`plot_bivariate_metric_tiles()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bivariate_metric_tiles.md)
