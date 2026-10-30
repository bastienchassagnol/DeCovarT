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
\\\boldsymbol{p}\\) stored in a list column `true_theta`. Scenario rows
are always evaluated **sequentially** to avoid nested parallelism;
sample-level workers live only in
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).
Scenario builders (factorial grids, overlap summaries, etc.) should live
in analysis scripts; see `scripts/fig02_bivariate_toy.R` and the
paper-scenario vignettes.

## Usage

``` r
run_simulation_benchmark(
  scenario_config,
  deconvolution_functions,
  n = 200,
  standardise = FALSE,
  scaled = FALSE,
  cores = 1L,
  coverage_interval = "wilson",
  verbose = FALSE,
  progress_every = 10L,
  likelihood = population_likelihood(),
  return_data = FALSE,
  reference = NULL,
  seed = NULL,
  scenarios = NULL
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
  Defaults to `1L`.

- coverage_interval:

  Coverage interval for the Monte Carlo coverage *rate*; see
  [`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md).

- verbose:

  If `TRUE`, print each scenario row and (when the grid has at most 10
  scenarios) every `progress_every` inferred samples. Large factorial
  grids log one line per scenario only, so logs stay readable.

- progress_every:

  Sample-progress interval passed to
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  when `verbose` is `TRUE`. Defaults to `10L`.

- likelihood:

  A likelihood object or its constructor, which both draws the bulks and
  scores the fits.
  [`population_likelihood()`](https://bastienchassagnol.github.io/DeCovarT/reference/population_likelihood.md)
  (default) draws every bulk from the convolution
  [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md);
  [`single_cell_likelihood()`](https://bastienchassagnol.github.io/DeCovarT/reference/population_likelihood.md)
  draws it as a sum of cells with
  [`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md),
  using `n_cells` and `platform_scale` from `true_theta` (and
  `generative_mean` / `generative_sigma` when the fitted reference
  differs from the generator), and scores every estimate against its
  realised composition.

- return_data:

  If `TRUE`, keep each scenario's bulk matrix, targets and true scales
  in `data`, so methods run outside this function can be scored on the
  same samples.

- reference:

  Optional list with `n_individuals`, `n_cells` (default `1`) and
  `seed`: with `return_data = TRUE`, an independent reference panel is
  drawn per scenario by
  [`simulate_reference_replicates()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_reference_replicates.md),
  with seed `reference$seed` plus the scenario index, for methods that
  estimate the reference moments themselves.

- seed:

  Optional integer. Scenario `i` then runs under
  `withr::with_seed(seed + i)`, so its draws depend only on its index: a
  scenario run alone (`scenarios = i`) reproduces the uninterrupted run.
  `NULL` (default) keeps the current random stream across scenarios.

- scenarios:

  Optional indices of the rows of `scenario_config` to run; the others
  are skipped. With `seed`, `scenarios = i` reproduces scenario `i` of
  the full run.

## Value

A list with:

- `regression`: `global` (per-sample composition scores) and `cell_type`
  (across-sample Pearson / F1 / spillover);

- `monte_carlo`: ADEMP summaries per cell type;

- `optimisation`: per-sample elapsed time, memory, KKT residual, and
  \\\hat{\boldsymbol{p}}\\;

- `config`: tibble of scenario metadata (one row per scenario);

- `theta_true`: list of convolution parameters (\\\boldsymbol{p}\\,
  \\\boldsymbol{\mu}\\, \\\boldsymbol{\Sigma}\_j\\) actually used to
  draw the bulk;

- `descriptors`: kept scenario statistics (composition, mean, SPD,
  network, tangent Fisher, MixSim BarOmega, pairwise Hellinger);

- `supplementary`: Jeffreys overlap, recorded separately;

- `intervals`: coverage and width of the interval methods requested in
  `deconvolution_functions`;

- `data`: per-scenario simulated data when `return_data = TRUE`, else
  `NULL`;

- `call`: the matched call
  ([`match.call()`](https://rdrr.io/r/base/match.call.html)). Analysis
  scripts persist those pieces with `write_simulation_artefacts()`
  (`scripts/benchmark_metrics/benchmark_artefacts.R`) rather than saving
  the whole list.

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md),
[`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md),
[`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md),
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
nrow(out$optimisation)
#> [1] 2
```
