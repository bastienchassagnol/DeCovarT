# Write split simulation artefacts (config, descriptors, theta, metrics)

[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
still returns a single list for tests. Scripts persist four RDS files
keyed by `ID` so the metrics object does not duplicate geometry, MixSim,
or `\theta`.

## Usage

``` r
write_simulation_artefacts(
  benchmark,
  dir,
  stem,
  config = NULL,
  rewrite_bivariate_id = FALSE
)
```

## Arguments

- benchmark:

  List from
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md).

- dir:

  Output directory.

- stem:

  File-name stem (`bivariate`, `hybrid`, …).

- config:

  Optional scenario grid (defaults to `benchmark$config`).

- rewrite_bivariate_id:

  If `TRUE`, rebuild fig02-style IDs.

## Value

Invisibly, a named list of written paths.

## See also

[`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md)
