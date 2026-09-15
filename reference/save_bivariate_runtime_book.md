# Solver wall-clock time for the bivariate toy (12 pages)

One page per CLD / variance / composition. The x-axis is the four
correlation corners; colour, fill, and grouping are solvers. y is
elapsed seconds on a log10 scale, with a left-side rug.

## Usage

``` r
save_bivariate_runtime_book(artefacts, file, data_rds = NULL)
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

[`save_hybrid_runtime_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_hybrid_runtime_book.md),
[`save_bivariate_memory_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_memory_book.md)

## Examples

``` r
if (requireNamespace("ggdist", quietly = TRUE)) {
  opt <- tibble::tibble(
    ID = rep(c("A", "B", "C", "D"), each = 8),
    sample_id = paste0("s", 1:32),
    algorithm = rep(c("lsei", "LBFGS"), 16),
    elapsed_sec = runif(32, 0.01, 0.2)
  )
  cfg <- tibble::tibble(
    ID = c("A", "B", "C", "D"),
    centroids = "small_CLD",
    variance = "homoscedastic",
    proportions = "balanced",
    correlation_celltype1 = c(0, -0.8, 0.8, -0.8),
    correlation_celltype2 = c(0, -0.8, 0.8, 0.8)
  )
  artefacts <- list(config = cfg, optimisation = opt)
  tf <- withr::local_tempfile(fileext = ".pdf")
  save_bivariate_runtime_book(artefacts, tf)
}
```
