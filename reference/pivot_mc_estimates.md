# Pivot Monte Carlo proportion estimates to a long table

Turns the `optimisation` block of
[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
into one row per replicate, algorithm, and cell type, with the matching
true proportion from `theta_true`. The raincloud intervals drawn from
this table are **empirical Monte Carlo quantiles** of
\\\hat{\boldsymbol{p}}\\, not confidence intervals for
\\\boldsymbol{p}\\ (Allen et al. 2019) .

## Usage

``` r
pivot_mc_estimates(benchmark)
```

## Arguments

- benchmark:

  List returned by
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md).

## Value

A tibble with `algorithm`, `cell_type`, `estimate`, `p_true`, `error`,
`sample_id`, and any scenario metadata columns.

## References

Allen M, Poggiali D, Whitaker K, Marshall TR, van Langen J, Kievit RA
(2019). “Raincloud Plots: A Multi-Platform Tool for Robust Data
Visualization.” *Wellcome Open Research*, **4**, 63. ISSN 2398-502X.
[doi:10.12688/wellcomeopenres.15191.2](https://doi.org/10.12688/wellcomeopenres.15191.2)
. <https://wellcomeopenresearch.org/articles/4-63>.

## See also

[`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md),
[`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md)

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
    ID = "B1",
    true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
  ),
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  n = 2,
  cores = 1
)
nrow(pivot_mc_estimates(out))
#> [1] 4
```
