# Three-page network book (increasing MixSim overlap)

Each page: CT3 (Scale-free) centred at the top; 2-by-2 of CT1 / CT2
topologies below. Edge widths use \\\|\Sigma\_{ij}\|\\ on the skeleton
with one shared scale across pages. Edge colour is the sign of
\\\Sigma\_{ij}\\ (inhibitory precision completion can yield negative
covariances).

## Usage

``` r
save_network_topology_book(artefacts, file, data_rds = NULL, png_file = NULL)
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

- png_file:

  Optional PNG copy (vignette / README).
