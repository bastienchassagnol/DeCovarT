# Simulate bulk mixtures and benchmark deconvolution algorithms

Wrapper around
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
and
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
(called inside
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)).
Each row of `scenario_config` describes one generative model
(\\\boldsymbol{\mu}\\, \\(\boldsymbol{\Sigma}\_j)\_j\\,
\\\boldsymbol{p}\\) stored in a list column `true_theta`. Scenario
builders (factorial grids, overlap summaries, etc.) should live in
analysis scripts; see `scripts/configure_bivariate_toy_scenarios.R` and
the manuscript scenario vignettes.

## Usage

``` r
run_simulation_benchmark(
  scenario_config,
  deconvolution_functions,
  n = 200,
  standardise = FALSE,
  scaled = FALSE,
  cores = NULL,
  parallel_scenarios = FALSE,
  parallel_cores = NULL
)
```

## Arguments

- scenario_config:

  Tibble or list of scenario rows. Each row must contain a `true_theta`
  list column (or list element) with at least `mu` and `sigma`. Optional
  per-row `n` overrides the default `n`.

- deconvolution_functions:

  Named list passed to
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md);
  each element has `FUN` and optional `additional_parameters` for
  [`do.call()`](https://rdrr.io/r/base/do.call.html).

- n:

  Default number of bulk replicates \\N\\ when `scenario_config` has no
  `n` column.

- standardise, scaled:

  Passed to
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

- cores:

  Workers for the per-sample loop inside
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).
  When `parallel_scenarios = TRUE`, defaults to `1` to avoid nested
  parallelism. Otherwise defaults to half of detected cores (at most
  \\\lfloor C/2\rfloor\\).

- parallel_scenarios:

  If `TRUE`, parallelise across scenario rows with `furrr` (optional
  Suggests `furrr` and `future`). Defaults to `FALSE`.

- parallel_cores:

  Maximum workers for scenario-level parallelism; defaults to half of
  detected cores.

## Value

A list with:

- `simulations`: tibble of per-sample estimates and metrics;

- `config`: tibble of scenario metadata (one row per scenario).

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)

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
scenario_config <- tibble::tibble(
  ID = "B1_Ho",
  true_theta = list(list(
    p = c(0.5, 0.5),
    mu = mu,
    sigma = Sigma
  ))
)
out <- run_simulation_benchmark(
  scenario_config = scenario_config,
  deconvolution_functions = list(
    "nnls" = list(FUN = deconvolute_ratios_nnls)
  ),
  n = 2,
  cores = 1
)
nrow(out$simulations)
#> [1] 2
```
