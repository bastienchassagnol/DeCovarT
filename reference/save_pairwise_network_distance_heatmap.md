# 3-by-3 Hellinger heatmaps of cell-type pairs versus MixSim overlap

Facets are pairwise comparisons (rows) by overlap (columns). Each panel
is the 2-by-2 of CT1 / CT2 graph families (CT3 is Scale-free). Hellinger
does not depend on \\p\\, so compositions are collapsed.

## Usage

``` r
save_pairwise_network_distance_heatmap(artefacts, file, data_rds = NULL)
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
