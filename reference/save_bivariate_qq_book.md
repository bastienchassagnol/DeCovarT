# Multi-page Q-Q book of whitened ILR coordinates versus N(0, 1)

Multi-page Q-Q book of whitened ILR coordinates versus N(0, 1)

## Usage

``` r
save_bivariate_qq_book(artefacts, file, data_rds = NULL)
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
