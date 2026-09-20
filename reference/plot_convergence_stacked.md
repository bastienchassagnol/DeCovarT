# Stacked bar counts of theoretical success and both failure modes

Stacked bar counts of theoretical success and both failure modes

## Usage

``` r
plot_convergence_stacked(df, title = NULL, icon_dir = NULL)
```

## Arguments

- df:

  Optimisation rows with `numerical_converged`, `theoretical_converged`,
  `algorithm`, and `panel`.

- title:

  Page title.

- icon_dir:

  Directory containing the three outcome PNGs (`succes_converence.png`,
  `failed_nuermical_convergence.png`, `failed_convergence.png`). When
  `ggtext` is installed, those icons are rendered in the fill legend.

## Value

A `ggplot` or cowplot grob.
