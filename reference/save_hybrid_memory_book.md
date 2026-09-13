# Solver peak memory (one page per composition)

Same layout as
[`save_hybrid_runtime_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_hybrid_runtime_book.md).
Memory is plotted in mebibytes (`memory_bytes / 2^20`) on a log10
y-axis.

## Usage

``` r
save_hybrid_memory_book(artefacts, file, data_rds = NULL)
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

[`save_hybrid_runtime_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_hybrid_runtime_book.md)

## Examples

``` r
if (requireNamespace("ggdist", quietly = TRUE)) {
  opt <- tibble::tibble(
    ID = rep("V1", 20),
    sample_id = paste0("s", 1:20),
    algorithm = rep(c("lsei", "LBFGS"), each = 10),
    memory_bytes = runif(20, 3.8e8, 4.2e8),
    graph_ct1 = "scale_free",
    graph_ct2 = "scale_free",
    overlap_label = "low",
    proportions = "balanced"
  )
  artefacts <- list(optimisation = opt)
  tf <- withr::local_tempfile(fileext = ".pdf")
  save_hybrid_memory_book(artefacts, tf)
}
```
