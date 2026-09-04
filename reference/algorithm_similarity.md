# Algorithm-similarity correlation from a Monte Carlo benchmark

Pearson correlations \\r\_{ab}=\mathrm{cor}(\hat p_a,\hat p_b)\\ across
Monte Carlo replicates (and cell types). This is **behavioural
similarity**: two solvers can correlate near 1 while remaining
systematically biased. Hierarchical order uses \\d\_{ab} = 1 -
r\_{ab}\\.

## Usage

``` r
algorithm_similarity(benchmark, facet_rows = NULL, facet_cols = NULL)
```

## Arguments

- benchmark:

  List from
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md).

- facet_rows, facet_cols:

  Optional scenario columns. When supplied, correlations are computed
  inside each facet cell.

## Value

A tibble with `algorithm_x`, `algorithm_y`, `correlation`, and any facet
columns.

## See also

[`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md),
[`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md)

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
out <- run_simulation_benchmark(
  tibble::tibble(
    true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
  ),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls),
    "rlm" = list(FUN = deconvolute_ratios_rlm)
  ),
  n = 4,
  cores = 1
)
algorithm_similarity(out)
#> # A tibble: 4 × 3
#>   algorithm_x algorithm_y correlation
#>   <chr>       <chr>             <dbl>
#> 1 nnls        nnls                  1
#> 2 rlm         nnls                  1
#> 3 nnls        rlm                   1
#> 4 rlm         rlm                   1
```
