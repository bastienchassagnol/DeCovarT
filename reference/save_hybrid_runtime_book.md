# Solver wall-clock time (one page per composition)

Facets MixSim overlap; topology is on the x-axis; colour / fill /
grouping are the five fig03 solvers. y is elapsed seconds on a log10
scale, with a left-side rug of the raw Monte Carlo draws.

## Usage

``` r
save_hybrid_runtime_book(artefacts, file, data_rds = NULL)
```

## Arguments

- artefacts:

  Split or assembled fig03 artefacts from
  [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md).

- file:

  Output PDF path.

- data_rds:

  Optional directory for ggplot `data` RDS files.

## Value

`file`, invisibly.

## See also

[`save_hybrid_memory_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_hybrid_memory_book.md),
[`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md)

## Examples

``` r
if (requireNamespace("ggdist", quietly = TRUE)) {
  opt <- tibble::tibble(
    ID = rep("V1", 20),
    sample_id = paste0("s", 1:20),
    algorithm = rep(c("lsei", "LBFGS"), each = 10),
    elapsed_sec = runif(20, 0.01, 0.2),
    graph_ct1 = "scale_free",
    graph_ct2 = "scale_free",
    overlap_label = "low",
    proportions = "balanced"
  )
  artefacts <- list(optimisation = opt)
  tf <- withr::local_tempfile(fileext = ".pdf")
  save_hybrid_runtime_book(artefacts, tf)
}
```
