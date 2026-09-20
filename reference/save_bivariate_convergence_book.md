# Multi-page stacked-bar book of solver outcomes

Multi-page stacked-bar book of solver outcomes

## Usage

``` r
save_bivariate_convergence_book(
  artefacts,
  file,
  data_rds = NULL,
  icon_dir = NULL
)
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

- icon_dir:

  Directory of outcome PNG icons.
