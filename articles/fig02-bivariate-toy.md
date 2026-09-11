# 2.1 Bivariate toy model (G = 2 genes, J = 2 cell types)

This article is the bivariate Gaussian-convolution toy of the methods
paper (G=2, J=2). It isolates gene–gene correlation from mean
separation. The factorial grid and ADEMP pipeline live in a **single**
script, `scripts/fig02_bivariate_toy.R`. Scenario builders
(`build_bivariate_scenario_config()`,
`bivariate_toy_deconvolution_functions()`) sit at the top of that file.
Sourcing the script from a vignette or from a temp copy in tests defines
those functions without launching the 972-scenario pipeline.

[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
wraps
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
and
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md).
Scenario rows are sequential; sample-level workers live only in
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).
Reporting conventions (ADEMP, Nature Methods, raincloud / forest plots)
are in [how to build synthetic
scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ademp).

``` default
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  P["3 compositions<br/>Ba / Mo / Hi"] --> S["972 scenarios"]
  V["2 variance structures<br/>Ho homoscedastic / He heteroscedastic"] --> S
  C["2 CLD mean separations<br/>Sm small / Lg large"] --> S
  R["9 x 9 gene-gene correlations<br/>rho in -0.8 to 0.8, step 0.2"] --> S
  S --> A["7 deconvolution solvers"]
```

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  P["3 compositions<br/>Ba / Mo / Hi"] --> S["972 scenarios"]
  V["2 variance structures<br/>Ho homoscedastic / He heteroscedastic"] --> S
  C["2 CLD mean separations<br/>Sm small / Lg large"] --> S
  R["9 x 9 gene-gene correlations<br/>rho in -0.8 to 0.8, step 0.2"] --> S
  S --> A["7 deconvolution solvers"]
```

Figure 1: Bivariate toy factorial design: three compositions, two
variance structures, two CLD mean separations, a 9-by-9 gene-gene
correlation grid, and seven solvers (3 x 2 x 2 x 9 x 9 = 972 scenarios).

``` default
---
config:
  layout: elk
  theme: sandstone
---
flowchart TB
  ROOT["Bivariate toy"]
  ROOT --> Sm["Small CLD"]
  ROOT --> Lg["Large CLD"]
  Sm --> SmHo["Homoscedastic"]
  Sm --> SmHe["Heteroscedastic"]
  Lg --> LgHo["Homoscedastic"]
  Lg --> LgHe["Heteroscedastic"]
  SmHo --> SmHoBa["Balanced"]
  SmHo --> SmHoMo["Moderately unbalanced"]
  SmHo --> SmHoHi["Highly unbalanced"]
  SmHe --> SmHeBa["Balanced"]
  SmHe --> SmHeMo["Moderately unbalanced"]
  SmHe --> SmHeHi["Highly unbalanced"]
  LgHo --> LgHoBa["Balanced"]
  LgHo --> LgHoMo["Moderately unbalanced"]
  LgHo --> LgHoHi["Highly unbalanced"]
  LgHe --> LgHeBa["Balanced"]
  LgHe --> LgHeMo["Moderately unbalanced"]
  LgHe --> LgHeHi["Highly unbalanced"]
```

``` mermaid
---
config:
  layout: elk
  theme: sandstone
---
flowchart TB
  ROOT["Bivariate toy"]
  ROOT --> Sm["Small CLD"]
  ROOT --> Lg["Large CLD"]
  Sm --> SmHo["Homoscedastic"]
  Sm --> SmHe["Heteroscedastic"]
  Lg --> LgHo["Homoscedastic"]
  Lg --> LgHe["Heteroscedastic"]
  SmHo --> SmHoBa["Balanced"]
  SmHo --> SmHoMo["Moderately unbalanced"]
  SmHo --> SmHoHi["Highly unbalanced"]
  SmHe --> SmHeBa["Balanced"]
  SmHe --> SmHeMo["Moderately unbalanced"]
  SmHe --> SmHeHi["Highly unbalanced"]
  LgHo --> LgHoBa["Balanced"]
  LgHo --> LgHoMo["Moderately unbalanced"]
  LgHo --> LgHoHi["Highly unbalanced"]
  LgHe --> LgHeBa["Balanced"]
  LgHe --> LgHeMo["Moderately unbalanced"]
  LgHe --> LgHeHi["Highly unbalanced"]
```

Figure 2: Hierarchical reading of the same factorial: small versus large
CLD, then homoscedastic versus heteroscedastic variances, then balanced
/ moderately unbalanced / highly unbalanced compositions (2 x 2 x 3 = 12
meta-scenarios). Each leaf is crossed with the 9-by-9 gene-gene
correlation grid.

## Generative model

With p_1+p_2=1, only one free ILR coordinate is estimated (see
[derivatives under simplex
transforms](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md)).
Even when mean profiles are similar, gene–gene correlation can degrade
mean-only solvers, and is partly recovered once
\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}\_j
enters the likelihood.

![](figures/fig_toy_model_parameters.svg)

\(a\) Hyperparameter grid used in the bivariate manuscript study: two
centroid geometries, balanced versus unbalanced \boldsymbol{p}, a
gene–gene correlation sweep, and the first-generation plus DeCovarT
solvers.

![](figures/fig_toy_2d_density.svg)

\(b\) AI-generated 2D density of the two-gene mixture while only the
pairwise correlation changes, at fixed means, marginal variances, and
composition. The ggplot panel below recomputes the same contrast and
annotates MixSim overlap.

Figure 3: Factorial design (left) and correlation-only 2D densities
(right) for the bivariate toy. {#fig-toy-design}

### Factorial design (972 scenarios)

| Factor | Levels |
|----|----|
| Composition \boldsymbol{p} | balanced (1/2,\\1/2); moderately unbalanced (17/20,\\3/20); highly unbalanced (99/100,\\1/100) |
| Mean separation (CLD) | small: \boldsymbol{\mu}\_{\cdot 1}=(20,22), \boldsymbol{\mu}\_{\cdot 2}=(22,20); large: (20,40), (40,20) |
| Gene–gene corr. CT 1 (\rho_1) | \\-0.8,\\-0.6,\\\ldots,\\0.8\\ — 9 levels |
| Gene–gene corr. CT 2 (\rho_2) | \\-0.8,\\-0.6,\\\ldots,\\0.8\\ — 9 levels |
| Variance structure | homoscedastic \sigma^2 = (1,1); heteroscedastic \sigma^2 = (1,2) |
| **Total scenarios** | 3 \times 2 \times 9 \times 9 \times 2 = \mathbf{972} |

N = 500 Monte Carlo replicates per scenario.

**CLD** is centroid (mean-profile) Euclidean separation
\lVert\boldsymbol{\mu}\_{\cdot 1}-\boldsymbol{\mu}\_{\cdot 2}\rVert_2,
not an Aitchison distance; see [the descriptor
notes](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#nte-desc-cld).
Small CLD uses \boldsymbol{\mu}\_{\cdot 1}=(20,22),
\boldsymbol{\mu}\_{\cdot 2}=(22,20) (d=2\sqrt{2}\approx 2.83, cosine
\approx 0.995). Large CLD uses (20,40) and (40,20) (d=20\sqrt{2}\approx
28.3, cosine 0.8).

#### Scenario ID

Each row is labelled `B{index}_{variance}_{composition}_{CLD}`. The
index is the factorial row number (1–972). Example: `B1_Ho_Ba_Sm` is the
first row, homoscedastic, balanced proportions, small CLD.

| Token | Meaning |
|----|----|
| `B` | Bivariate toy (G=2, J=2) |
| `1`–`972` | Factorial row index |
| `Ho` / `He` | Homoscedastic / heteroscedastic marginal variances |
| `Ba` / `Mo` / `Hi` | Balanced / moderately unbalanced / highly unbalanced \boldsymbol{p} |
| `Sm` / `Lg` | Small / large CLD (Euclidean centroid gap) |

The slim design grid is `output/fig02/bivariate_config.rds` (no
`true_theta`). Geometry, MixSim `BarOmega`, Hellinger, and Jeffreys live
in `bivariate_descriptors.rds`. Convolution parameters are
`bivariate_theta.rds`. ADEMP metrics (`regression`, `monte_carlo`,
`optimisation`, `call`) are `bivariate_benchmark.rds`. All four files
share the `ID` key.
`read_simulation_artefacts(dir, "bivariate", assemble = TRUE)` rejoins
them for plotting.

The Monte Carlo \hat{\boldsymbol{p}} draws themselves are the cell-type
columns of `optimisation` (one row per replicate \times algorithm \times
scenario).
[`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md)
stacks those columns into long form. Interior Wald intervals for
convolution-likelihood solvers use
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
(expected Fisher information through the ILR chart), matching
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).
Coverage, `mean_model_se`, and `se_sd_ratio` (mean model SE / empirical
SD) are then finite. Mean-only solvers (NNLS, LSEI) leave those Wald
columns missing: they do not maximise the convolution likelihood, so an
ILR interval would not apply.

``` r

library(DeCovarT)
source("scripts/fig02_bivariate_toy.R")

bivariate_config <- build_bivariate_scenario_config()
nrow(bivariate_config)
```

### Log-likelihood along the simplex

For J=2 the free coordinate is p_1\in(0,1). The listing below evaluates
\ell\_{\boldsymbol{y}\mid\boldsymbol{\zeta}}(p_1,1-p_1) on one bulk
draw. It is not executed during `R CMD build`: Quarto’s Windows CLI
starts a new R process that cannot see the temporary library
([quarto-r#217](https://github.com/quarto-dev/quarto-r/issues/217)). The
same likelihood path is checked in `tests/testthat/`.

``` r

set.seed(20260828)
rho <- 0.6
R <- matrix(c(1, rho, rho, 1), nrow = 2)
Sigma <- array(c(R, R), dim = c(2, 2, 2))
dimnames(Sigma) <- list(
  rownames(mu_small),
  rownames(mu_small),
  colnames(mu_small)
)
y_one <- DeCovarT::simulate_bulk_mixture(
  mu_small,
  Sigma,
  p = p_bal,
  n = 1
)$Y[, 1]
p1_grid <- seq(0.02, 0.98, by = 0.01)
ll_grid <- vapply(
  p1_grid,
  function(p1) {
    DeCovarT::loglik_multivariate(c(p1, 1 - p1), y_one, mu_small, Sigma)
  },
  numeric(1)
)
ll_df <- data.frame(p1 = p1_grid, loglik = ll_grid)
ggplot2::ggplot(ll_df, ggplot2::aes(x = p1, y = loglik)) +
  ggplot2::geom_line(colour = "#1B4F72", linewidth = 0.8) +
  ggplot2::geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    colour = "grey30"
  ) +
  ggplot2::labs(
    x = "p1 (p2 = 1 - p1)",
    y = "Log-likelihood"
  ) +
  ggplot2::theme_bw(base_size = 11)
```

## Inference

Seven solvers are evaluated on each scenario:

| Solver | Description |
|----|----|
| `NNLS` | Non-negative least squares (mean-only baseline) |
| `DeconRNASeq` | NNLS via `lsei` ([Gong and Szustakowski 2013](#ref-gongDeconRNASeqStatisticalFramework2013)) |
| `L-BFGS-B` | Quasi-Newton on box-constrained simplex |
| `gradient` | Projected gradient descent in ILR |
| `Newton-Raphson` | Second-order in ILR |
| `Marquardt-Levenberg` | Damped Newton via `marqLevAlg` |
| `SA` | Simulated annealing |

Performance is summarised using the ADEMP framework ([Morris et al.
2019](#ref-morrisUsingSimulationStudies2019)): **bias**, **RMSE**,
**coverage** of ILR Wald intervals (Wilson intervals on the coverage
*rate*), **SE/SD**, and optimiser **failure rate**. Mean-only baselines
have no convolution Wald SE.

``` r

N_REPLICATES <- as.integer(Sys.getenv("N_REPLICATES", unset = "500"))
deconvolution_functions <- bivariate_toy_deconvolution_functions(
  itmax = 200L,
  epsilon = 1e-4
)
bivariate <- run_simulation_benchmark(
  scenario_config = bivariate_config,
  deconvolution_functions = deconvolution_functions,
  n = N_REPLICATES,
  cores = 1L
)
```

`config` on the in-memory benchmark still carries Shannon entropy of
\boldsymbol{p} and MixSim overlap for tests. Saved artefacts split that
information: `optimisation` stores per-sample \hat{\boldsymbol{p}},
elapsed time, and memory; `regression` and `monte_carlo` are the
composition and ADEMP blocks from
[`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md).

## Visualisations

Figures are written from the split RDS (no ADEMP refit). Density books
go to `output/fig02/density_visualisations/`; RMSE / MAE / Aitchison
tiles, raincloud, forest, similarity, and solver-dot books go to
`output/fig02/performance_visualisations/`. Tile heatmaps
([`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html))
show RMSE, MAE, and Aitchison distance on the (\rho_1,\rho_2) plane: one
PDF per metric, one page per meta-scenario (CLD () variance ()
composition), seven solver panels. Scenario factors are relevelled
(small then large CLD; homoscedastic then heteroscedastic; balanced,
moderately unbalanced, highly unbalanced) when artefacts are read, and
solvers follow `nnls`, `lsei`, `SA`, `gradient`, `LBFGS`,
`Newton-Raphson`, `Marquardt-Levenberg`. Density and performance books
share the four correlation corners
(\rho_1,\rho_2)\in\\(0,0),(-0.8,-0.8),(0.8,0.8),(-0.8,0.8)\\ (2\times
2\times 3=12 pages). Similarity heatmaps compute Pearson (r) of paired
() **within each corner** (not pooled across the 9-by-9 grid), with
average-linkage clustering of (1-r) and a dendrogram on the right of the
tiles. Wald forests drop NNLS, LSEI, and SA; solid whiskers are (1.96)
times the expected-Fisher Wald SE at the true composition, dashed
whiskers use the empirical Monte Carlo SD, and bold two-line labels (one
per solver, centred between cell types) report RMSE and coverage, which
are the same for (p_1) and (p_2) on the simplex. Rainclouds plot (p)
(not the error) with vertical lines at the two true proportions.
Solver-dot pages use the same 9-by-9 ((\_1,\_2)) grid as the RMSE
heatmaps: colour is mean RMSE and size is mean Aitchison distance, with
horizontal y-axis labels. ggplot `data` tables are saved under
`output/fig02/ggplot_rds/`.

#### Expected findings

With **large mean separation** and \rho_j = 0 both groups of algorithms
should perform well. As \|\rho\| increases and the means move closer,
NNLS and DeconRNASeq degrade faster than DeCovarT solvers because they
ignore covariance information. The MSE heatmap should display a
characteristic saddle: lowest error on the diagonal \rho_1 = \rho_2 and
highest near zero mean separation.

### See also

- Variance-driven hybrid (G=20, J=3):
  [§2.2](https://bastienchassagnol.github.io/DeCovarT/articles/fig03-covariance-driven.md)
- Moment generator and ADEMP reporting: [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.md)
- Regular-case MLE checks: [Appendix
  S1](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.md)

### References

Gong, Ting, and Joseph D. Szustakowski. 2013. ‘DeconRNASeq: A
Statistical Framework for Deconvolution of Heterogeneous Tissue Samples
Based on mRNA-Seq Data’. *Bioinformatics (Oxford, England)* 29.
<https://doi.org/10.1093/bioinformatics/btt090>.

Morris, Tim P., Ian R. White, and Michael J. Crowther. 2019. ‘Using
Simulation Studies to Evaluate Statistical Methods’. *Statistics in
Medicine* 38 (11): 2074–102. <https://doi.org/10.1002/sim.8086>.
