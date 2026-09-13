# Funky heatmap of scenario-level geometry metrics

Rows nest Shannon composition, MixSim overlap, then the four CT1/CT2
graph pairs. Circles are min–max scaled within each column. No
clustering.

## Usage

``` r
save_scenario_metrics_funkyheatmap(artefacts, file, data_rds = NULL)
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
