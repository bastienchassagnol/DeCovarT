# Read split simulation artefacts and optionally reassemble a benchmark list

Read split simulation artefacts and optionally reassemble a benchmark
list

## Usage

``` r
read_simulation_artefacts(dir, stem, assemble = FALSE)
```

## Arguments

- dir:

  Output directory.

- stem:

  File-name stem (`bivariate`, `hybrid`, …).

- assemble:

  If `TRUE`, join config onto metric tables for plotting helpers that
  still expect design columns.

## Value

Named list of tibbles, or a benchmark-like list when `assemble = TRUE`.
