# Multi-page PDF of quantile-dot rainclouds at four corners

Same layout as
[`save_bivariate_raincloud_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_raincloud_book.md),
but
[`ggdist::stat_dotsinterval()`](https://mjskay.github.io/ggdist/reference/stat_dotsinterval.html)
replaces the bounded half-eye KDE.

## Usage

``` r
save_bivariate_dotsinterval_book(artefacts, file, data_rds = NULL)
```

## Arguments

- artefacts:

  List from
  [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md)
  with `assemble = TRUE`, or the same named pieces.

- file:

  Output PDF path.

- data_rds:

  Optional directory for ggplot `data` RDS files.
