# 12-page PDF of RMSE/Aitchison solver dots on the correlation grid

One page per meta-scenario. Each page facets solvers on the same 9-by-9
\\(\rho_1,\rho_2)\\ factorial as the metric heatmaps (enrichplot-style
dots: colour = RMSE, size = Aitchison).

## Usage

``` r
save_bivariate_solver_dots_book(artefacts, file, data_rds = NULL)
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
