# Multi-page Q-Q book of joint Mahalanobis D^2 versus chi-square

Multi-page Q-Q book of joint Mahalanobis D^2 versus chi-square

## Usage

``` r
save_bivariate_chi2_book(artefacts, file, data_rds = NULL)
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
