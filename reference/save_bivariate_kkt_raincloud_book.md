# Multi-page KKT raincloud book

Multi-page KKT raincloud book

## Usage

``` r
save_bivariate_kkt_raincloud_book(artefacts, file, data_rds = NULL)
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
