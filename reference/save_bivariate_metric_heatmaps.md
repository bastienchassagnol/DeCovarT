# Write RMSE / MAE / Aitchison tile PDFs for the bivariate toy

Write RMSE / MAE / Aitchison tile PDFs for the bivariate toy

## Usage

``` r
save_bivariate_metric_heatmaps(artefacts, dir, data_rds = NULL)
```

## Arguments

- artefacts:

  List from
  [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md)
  with `assemble = TRUE`, or the same named pieces.

- dir:

  Output directory.

- data_rds:

  Optional directory for ggplot `data` RDS files.

## Value

Named paths.
