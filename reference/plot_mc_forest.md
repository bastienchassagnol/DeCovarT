# Forest plot of ADEMP Monte Carlo summaries

Dot-and-whisker display of bias, RMSE, MAE, coverage, mean interval
width, and optimiser failure rate by algorithm and cell type (Allen et
al. 2019) . Coverage whiskers are the Wilson interval already stored on
`monte_carlo`
([`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md);
(Wilson 1927) ): they are intervals for the coverage *rate*, not for
\\p_j\\. Bias is referenced at 0; coverage at 0.95. Pairwise algorithm
contrasts (MAE differences versus a reference solver on the same Monte
Carlo replicates) can be read from the raincloud of paired errors; they
do not need a second bootstrap.

## Usage

``` r
plot_mc_forest(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL,
  metrics = c("bias", "rmse", "mae", "coverage", "mean_interval_width", "failure_rate")
)
```

## Arguments

- benchmark:

  List from
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md),
  or a table from
  [`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md).

- facet_rows, facet_cols:

  Optional column names for
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
  rows and columns.

- metrics:

  Character vector of summaries to display.

## Value

A `ggplot` object.

## References

Allen M, Poggiali D, Whitaker K, Marshall TR, van Langen J, Kievit RA
(2019). “Raincloud Plots: A Multi-Platform Tool for Robust Data
Visualization.” *Wellcome Open Research*, **4**, 63. ISSN 2398-502X.
[doi:10.12688/wellcomeopenres.15191.2](https://doi.org/10.12688/wellcomeopenres.15191.2)
. <https://wellcomeopenresearch.org/articles/4-63>.  
  
Wilson EB (1927). “Probable Inference, the Law of Succession, and
Statistical Inference.” *Journal of the American Statistical
Association*, **22**(158), 209–212. ISSN 0162-1459.
[doi:10.2307/2276774](https://doi.org/10.2307/2276774) .

## See also

[`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md),
[`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md)

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
    n_genes = 2L,
    cosine = 0,
    true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
  ),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  n = 4,
  cores = 1
)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_mc_forest(out, facet_cols = "cosine")
}
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```
