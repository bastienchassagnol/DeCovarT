# Kolmogorov-Smirnov p-values versus N(0, 1) by algorithm

One lollipop per scenario (panel) and solver, pooling whitened ILR
coordinates. Values below 0.05 reject normality at the conventional 5%
level.

## Usage

``` r
plot_ks_normality_box(df, title = NULL)
```

## Arguments

- df:

  Long whitened table.

- title:

  Page title.

## Value

A `ggplot` or cowplot grob with an inset legend.
