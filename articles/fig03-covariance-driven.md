# 2.2 Covariance-driven scenario (G = 20, J = 3)

``` r

library(DeCovarT)
```

The bivariate toy confirms the theoretical benefit of covariance
modelling in a transparent setting but cannot demonstrate realistic
gene-regulatory network (GRN) topology. This scenario uses G = 20 genes
and J = 3 cell types in which types 1 and 2 are **mean-collinear**
(\cos(\boldsymbol{\mu}\_{\cdot 1},\boldsymbol{\mu}\_{\cdot 2})=0.9): the
network topology encoded in \boldsymbol{\Sigma}\_j is the remaining
discriminative signal. Construction, ADEMP benchmark, and the static
network figure live in a **single** script,
`scripts/fig03_variance_covariance_driven.R` (seed `20260807`).

![](fig03-covariance-driven_files/mediabag/fig-covariance-design.svg)

Figure 1: Covariance-driven factorial: one fixed Gram (cosines 0.9 /
0.1), a 2-by-2 graph assignment on the mean-collinear types (scale-free
versus cluster SBM; type 3 held at scale-free), three target MixSim
average overlaps, three Shannon compositions (4 x 3 x 3 = 36 scenarios),
analysed by eight methods.

``` default
---
config:
  layout: elk
  theme: sandstone
---
flowchart TB
  ROOT["Covariance-driven G=20, J=3"]
  ROOT --> SfSf["CT1 SF, CT2 SF"]
  ROOT --> SfSb["CT1 SF, CT2 SBM"]
  ROOT --> SbSf["CT1 SBM, CT2 SF"]
  ROOT --> SbSb["CT1 SBM, CT2 SBM"]
  SfSf --> SfSfLo["Low overlap"]
  SfSf --> SfSfMo["Moderate overlap"]
  SfSf --> SfSfHi["High overlap"]
  SfSb --> SfSbLo["Low overlap"]
  SfSb --> SfSbMo["Moderate overlap"]
  SfSb --> SfSbHi["High overlap"]
  SbSf --> SbSfLo["Low overlap"]
  SbSf --> SbSfMo["Moderate overlap"]
  SbSf --> SbSfHi["High overlap"]
  SbSb --> SbSbLo["Low overlap"]
  SbSb --> SbSbMo["Moderate overlap"]
  SbSb --> SbSbHi["High overlap"]
```

``` mermaid
---
config:
  layout: elk
  theme: sandstone
---
flowchart TB
  ROOT["Covariance-driven G=20, J=3"]
  ROOT --> SfSf["CT1 SF, CT2 SF"]
  ROOT --> SfSb["CT1 SF, CT2 SBM"]
  ROOT --> SbSf["CT1 SBM, CT2 SF"]
  ROOT --> SbSb["CT1 SBM, CT2 SBM"]
  SfSf --> SfSfLo["Low overlap"]
  SfSf --> SfSfMo["Moderate overlap"]
  SfSf --> SfSfHi["High overlap"]
  SfSb --> SfSbLo["Low overlap"]
  SfSb --> SfSbMo["Moderate overlap"]
  SfSb --> SfSbHi["High overlap"]
  SbSf --> SbSfLo["Low overlap"]
  SbSf --> SbSfMo["Moderate overlap"]
  SbSf --> SbSfHi["High overlap"]
  SbSb --> SbSbLo["Low overlap"]
  SbSb --> SbSbMo["Moderate overlap"]
  SbSb --> SbSbHi["High overlap"]
```

Figure 2: Hierarchical reading of the covariance factorial: 2-by-2 graph
assignment on mean-collinear cell types 1 and 2 (type 3 held at
scale-free), then low / moderate / high target MixSim BarOmega (4 x 3 =
12 covariance draws). Each leaf is crossed with the three Shannon
compositions H\* = 1 / 0.5 / 0.1.

## Generative model

One fixed mean signature \boldsymbol{\mu}\in\mathbb{R}^{G\times J} is
built with
[`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
and an explicit Gram `target_gram` ([mean signature
construction](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-lr-block-mu)).
There is no marker block and no null (`equal_all`) block: uninformative
genes would be removed before deconvolution.

R = \begin{pmatrix} 1 & 0.9 & 0.1 \\ 0.9 & 1 & 0.1 \\ 0.1 & 0.1 & 1
\end{pmatrix}, \qquad \boldsymbol{\mu} = s\mathbf{Q}R^{1/2}, \quad s=10.
\tag{1}

Because \mathbf{Q}^{\mathsf{T}}\mathbf{Q}=\mathbf{I}\_J, the pairwise
cosines of the columns of \boldsymbol{\mu} equal R exactly (up to
rounding). Euclidean gaps follow \lVert\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\_{\cdot k}\rVert =s\sqrt{2(1-R\_{jk})}. The call
uses `nonnegative = TRUE` (disjoint gene-block frame) so that
\boldsymbol{\mu}\ge 0, as required by
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

| pair                     | cosine |                  Euclidean |
|:-------------------------|-------:|---------------------------:|
| celltype_1 vs celltype_2 |    0.9 |  10\sqrt{0.2}\approx 4.472 |
| celltype_1 vs celltype_3 |    0.1 | 10\sqrt{1.8}\approx 13.416 |
| celltype_2 vs celltype_3 |    0.1 | 10\sqrt{1.8}\approx 13.416 |

Table 1: Exact pairwise geometry of the 20\times 3 mean signature
(`target_gram`, `mean_scale` =10).

Each of the mean-collinear types (CT1, CT2) independently draws an
undirected skeleton from one of two generators, then completes a signed
precision with
[`build_normalised_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_normalised_precision.md)
([GGM
networks](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ggm-networks);
topology families in [the topology
table](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#tbl-topologies)):

- **Scale-free** (`graph_model = "scale_free"`): Barabási–Albert
  preferential attachment (`power = 1`, `edges_per_node = 1`).
- **Cluster / stochastic block model**
  (`graph_model = "stochastic_block_model"`): three blocks with
  probabilities (0.5,0.25,0.25), p\_{\mathrm{within}}=0.25,
  p\_{\mathrm{between}}=0.01.

Erdős–Rényi is not used. Type 3, already separated by the Gram
(R\_{13}=R\_{23}=0.1), is held at scale-free, so the topology factorial
is 2^{2}=4 assignments on CT1 and CT2. Those generator names and
parameter lists are stored on every row of the `config` element of
`output/fig03/fig03_benchmark.rds` as `graph_ct1` / `graph_ct2` /
`graph_ct3` and the list-column `graph_generators`; the true moments of
each scenario are in its `theta_true` element.

The spectral cushion u used at completion is now a **fixed** engineering
constant (`PRECISION_SHIFT_BASE = 0.2`) so that each graph has one SPD
precision. The experimental knob is a global covariance scale s:
\boldsymbol{\Sigma}\_j\mapsto s\boldsymbol{\Sigma}\_j keeps every
structural zero of \boldsymbol{\Omega}\_j because
\boldsymbol{\Omega}\_j(s)=\boldsymbol{\Omega}\_j/s.
[`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md)
searches s so that MixSim `BarOmega` matches a target in
\\0.02,0.08,0.20\\ (low / moderate / high), using Sobol Monte Carlo for
G=20 ([overlap
note](https://bastienchassagnol.github.io/DeCovarT/articles/theory-distance-covariance.html#sec-fig03-scale)).
This is the graph-constrained analogue of MixSim / FSDA average-overlap
simulation ([CRAN
MixSim](https://cran.r-project.org/web/packages/MixSim/refman/MixSim.html#MixSim);
([Melnykov et al. 2012](#ref-melnykovMixSimPackageSimulating2012);
[Riani et al. 2015](#ref-rianiSimulatingMixturesMultivariate2015))): we
do **not** call `MixSim()` as a generator, because it cannot take a
precision support as input.

Larger s also raises the covariance-information fraction
f\_{\mathrm{cov}} (mean Fisher scales as 1/s; the covariance block is
scale-invariant). See the [tangent-Fisher
callout](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#nte-desc-fisher).

| Factor                 | Levels                                   |
|------------------------|------------------------------------------|
| Graph model (CT1, CT2) | scale-free; cluster SBM                  |
| Graph model (CT3)      | scale-free (fixed)                       |
| Target `BarOmega`      | low (0.02); moderate (0.08); high (0.20) |
| Topology assignments   | 2^{2}=\mathbf{4}                         |
| Covariance draws       | 4\times 3=\mathbf{12}                    |

Table 2: Covariance factorial (graphs \times target average overlap).

[`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md)
records the realised `BarOmega`, f\_{\mathrm{cov}},
\kappa\\\boldsymbol{\Sigma}(\boldsymbol{p})\\, and pairwise AIRM
distance of the \boldsymbol{\Sigma}\_j on every row. The overlap target
is matched at the **balanced** composition; unbalanced \boldsymbol{p}
change MixSim’s MAP weights, so realised `BarOmega` is higher there.
AIRM is invariant to the global scale s (it scores topology / shape).
f\_{\mathrm{cov}} rises with s.

![](figures/fig_network_topologies.png)

Figure 3: Two generator examples on G=20 nodes (scale-free, cluster
SBM). Mean-collinear cell types 1 and 2 are each assigned one of these
families; type 3 is held at scale-free.

Two generators were considered but not adopted:

- **Toeplitz covariances**: encode only local (banded) dependencies;
  cannot model cell-type-specific sparse precision.
- **`MixSim()` as a mixture generator**: draws unconstrained means and
  covariances for a target average overlap; no graph-constrained
  precision input. DeCovarT keeps the graphs and matches `BarOmega` by
  scaling \boldsymbol{\Sigma}\_j instead.

## Inference

Compositions are one-dominant vectors from
[`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md)
targeting normalised Shannon entropy H^{\star}\in\\1,0.5,0.1\\
(balanced; moderately unbalanced; highly unbalanced). All DeCovarT
solvers start at the barycentre (`initial_p = "barycentre"`). Other
initialisation strategies are reserved for separate scenarios.

| Factor                     | Levels                                    |
|----------------------------|-------------------------------------------|
| Composition \boldsymbol{p} | H^{\star}=1; H^{\star}=0.5; H^{\star}=0.1 |
| **Total scenarios**        | 4\times 3\times 3=\mathbf{36}             |

Table 3: Composition factor.

The run size is hard-coded at the top of the script:

| Constant | Meaning | Value |
|----|----|----|
| `N_REPLICATES` | Monte Carlo bulks per scenario, N | 200 |
| `N_REFERENCE` | individuals of the reference panel, R (MEAD and RNA-Sieve) | 200 |
| `CORES` | workers for the N bulks of a scenario; the three figure scripts run at the same time | a third of half the cores (14 of 88) |
| `POSTPROCESS_ONLY` | `TRUE` redraws the figures from the saved benchmark | `FALSE` |

Table 4: Run-size constants of
`scripts/fig03_variance_covariance_driven.R`.

Parallel workers load the **installed** DeCovarT, so run
`devtools::install()` before a parallel run.

### Methods

Every method is scored on the same simulated bulks. Feature selection is
switched off everywhere: all G=20 genes enter every fit.

| Method | Function or code | Input | Intervals |
|:---|:---|:---|:---|
| `lsei` (DeconRNASeq) | [`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true means | none |
| NNLS | [`deconvolute_ratios_nnls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true means | none |
| RLM (robust linear model) | [`deconvolute_ratios_rlm()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true means | none |
| CIBERSORT (\nu-SVR) | [`deconvolute_ratios_cibersort()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true means | none |
| Newton–Raphson | [`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true moments | Wald, ILR, profile, oracle |
| L-BFGS-B | [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md) | true moments | Wald, ILR |
| Firth | [`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md) | true moments | Wald, ILR |
| MEAD | `MEAD::MEAD()` (GitHub, commit `d66fa10`) | reference panel, cohort of N bulks | native sandwich interval |
| RNA-Sieve | `rnasieve` 0.1.4 (Python) | reference panel, cohort of N bulks | released, Fisher and Godambe (see below) |
| ReDeconv | `redeconv` 1.2.1 (Python) | true means and standard deviations | none |

Table 5: Methods of the covariance-driven benchmark.

- **`lsei` and NNLS.** NNLS followed by renormalisation is least squares
  on the simplex with a **free** scale (\boldsymbol\beta=c\boldsymbol p,
  c=\mathbf 1^{\top}\boldsymbol\beta); `lsei` fixes c=1. They agree only
  when \hat c=1. In population scenarios the true scale is 1, so `lsei`
  uses information NNLS discards; both are run, to show the difference.
  On 500 bulks of a design like this one with \boldsymbol
  p=(0.84,0.08,0.08), the largest coordinate difference between the two
  has median 0.03, 95th percentile 0.11 and maximum 0.19.
- **L-BFGS-B** runs on the closure chart \boldsymbol p=\boldsymbol
  w/\mathbf 1^{\top}\boldsymbol w with the box \boldsymbol w\ge 0
  ([closure
  chart](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-generative-optimisation.html#sec-closure-box)):
  it maximises the convolution likelihood itself, reaches exact zeros,
  and its box KKT points are the simplex KKT points. A boundary estimate
  (one, two or up to J-1 zero proportions: a face, an edge or a vertex)
  has no Wald or ILR interval, so their coverage for L-BFGS-B is over
  interior estimates; `n_boundary` gives the count. Profile intervals
  refit from the bulk and are reported under Newton–Raphson only.
- **Firth** maximises the Jeffreys-penalised likelihood ([Firth
  1993](#ref-firthBiasReductionMaximum1993)); its estimate is the same
  for both optimisers, so it appears once.
- **Second-generation methods** run from their authors’ code, with gene
  selection and preprocessing filters off and the other defaults kept
  (`scripts/benchmark_metrics/second_generation_methods.R`; the Python
  ones through a `uv` project in `scripts/benchmark_metrics/python/`,
  exchanging CSV files with R). MEAD and RNA-Sieve read a reference
  panel of R purified profiles per type drawn **independently** of the
  bulks
  ([`simulate_reference_replicates()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_reference_replicates.md),
  its own seed), and are fitted on the cohort of the N bulks of a
  scenario, whereas DeCovarT fits one bulk at a time.
  - MEAD rescales every reference profile to a unit total, so its
    proportions are transcript fractions; its intervals are scored
    against both the cell and the transcript fractions.
  - The released interval code of RNA-Sieve returns \[0,1\] whenever its
    averaged Hessian has a condition number above 10^{6} (`stats.py`,
    lines 224 and 234 of version 0.1.4): it then computes \boldsymbol
    H\boldsymbol B\boldsymbol H^{\top} instead of \boldsymbol
    H^{-1}\boldsymbol B\boldsymbol H^{-\top} on a gene subsample. The
    runner also returns the two intervals of the paper, recomputed on
    all genes with the package’s own per-gene scores and averaged
    Hessian: Fisher (\boldsymbol B^{-1}/G, the empirical Fisher
    information) and Godambe (\boldsymbol H^{-1}\boldsymbol B\boldsymbol
    H^{-\top}/G), with \alpha_K=1-\sum\_{k\<K}\alpha_k. In the latest
    dry run the faulty branch never fired, so the released and Godambe
    intervals coincide; RNA-Sieve’s solver converged in all 36 scenarios
    of that run, but failed in one to three of them in earlier draws
    (recorded as a numerical failure, a grey tile in the heatmaps).
  - DECALS has no public implementation and is not run; its point
    estimate, constrained least squares solved as a quadratic programme,
    is the `lsei` baseline.

### Interval constructions

| Interval | Construction |
|:---|:---|
| `wald` | delta method on \boldsymbol p, symmetric |
| `ilr` | Wald interval in ILR coordinates mapped back to the simplex |
| `profile` | likelihood-ratio inversion; chart-invariant, inside \[0,1\] |
| `oracle` | \hat p_j\pm the Wald half-width at \boldsymbol p^{\star} (Fisher at the truth) |

Table 6: Interval menu ([population
inference](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.md)).

Coverage is reported with a Wilson Monte Carlo interval ([Morris et al.
2019](#ref-morrisUsingSimulationStudies2019)); the interval table also
gives the mean and median width and the share of intervals leaving
\[0,1\].

``` r

benchmark <- readRDS("output/fig03/fig03_benchmark.rds")
names(benchmark)
#> meta config theta_true descriptors supplementary regression
#> monte_carlo optimisation intervals data
dplyr::filter(benchmark$intervals, is_composition) |>
  dplyr::select(ID, algorithm, ci_method, parameter, coverage, mean_width)
```

## Visualisations

Every figure is redrawn from two files in `output/fig03/`:
`fig03_benchmark.rds` (run settings, design, true parameters and every
metric table) and `fig03_plot_data.rds`, a list keyed by figure file
stem that holds the data handed to each plot function.

``` r

source("scripts/visualisation/load_visualisation.R")
plot_data <- readRDS("output/fig03/fig03_plot_data.rds")
names(plot_data)
render_figure_from_data("ci_width_forest", plot_data, "output/fig03")
```

| Output | Description |
|----|----|
| `density_visualisations/mean_signature.pdf` | Unscaled \boldsymbol{\mu} heatmap |
| `density_visualisations/network_topologies.pdf` | Graph skeletons, one page per overlap target; edges coloured by \Sigma\_{ij}, nodes by SBM block with a thin black outline, scale-free hubs ringed |
| `density_visualisations/scenario_metrics.pdf` | Scenario descriptors as shapes, one geometry and palette per column group, each column min–max scaled; the conditioning of the mean signature, shared by the scenarios, is in the caption |
| `density_visualisations/scenario_metrics_table.html`, `.pdf` | The same descriptors as numbers (reactable, tinytable); the PDF marks the easiest value of each column in green and the hardest in red, and shades three columns (entropy, mean overlap, condition number) |
| `density_visualisations/latent_projections.pdf` | PCA of the pooled purified profiles, with a scree inset (share of the variance of the first ten components), and supervised `MclustDR` plane |
| `density_visualisations/pairwise_network_aware_distances.pdf` | Pairwise discrepancy of the purified laws |
| `density_visualisations/loglik_surface_p.pdf` | Log-likelihood of one bulk on the balance plane 1\mid2, \\1,2\\\mid3 |
| `density_visualisations/loglik_expected_surface_p.pdf` | Expected log-likelihood; red point = Monte Carlo mean of \hat{\boldsymbol p} |
| `density_visualisations/loglik_rgl.pdf` | 3D relative likelihood (`rgl`) |
| `performance_visualisations/heatmap_rmse.pdf`, `heatmap_aitchison.pdf` | RMSE and Aitchison tiles; grey = no finite estimate (solver failure) |
| `performance_visualisations/ci_width_forest.pdf` | Mean interval, coverage \[Wilson interval\] and width, one page per scenario and cell type; RNA-Sieve released, Fisher and Godambe rows |
| `performance_visualisations/p_hat_monte_carlo_raincloud.pdf` | Monte Carlo law of \hat{\boldsymbol{p}} |
| `performance_visualisations/p_hat_monte_carlo_quantile_dots.pdf` | The same law as quantile dots |
| `performance_visualisations/kkt_residual_raincloud.pdf` | Projected-score KKT residual |
| `performance_visualisations/convergence_theoretical_vs_numerical_failures.pdf` | Theoretical success, numerical failure, theoretical failure |
| `performance_visualisations/wald_normality_qq.pdf` | Whitened ILR estimate against N(0,1) |
| `performance_visualisations/wald_mahalanobis_chi2_qq.pdf` | D^{2} against \chi^{2}\_{J-1}, \log\_{10} axes, coloured by convergence |
| `performance_visualisations/estimate_agreement_hellinger_ward.pdf` | Mean Hellinger distance between paired estimates, Ward clustering; dashed line at the mean over method pairs |
| `performance_visualisations/rmse_colour_aitchison_size_dots.pdf` | Colour = RMSE, size = Aitchison distance |
| `performance_visualisations/runtime.pdf`, `memory.pdf` | \log\_{10} seconds and \log\_{10} MiB |

Table 7: Plot books of the covariance-driven scenario.

Redraw every book without refitting with `POSTPROCESS_ONLY <- TRUE` at
the top of `scripts/fig03_variance_covariance_driven.R`.

### Two-dimensional projections of a Gaussian convolution

With G=20 the purified Gaussians and the bulk convolution are no longer
drawable in the gene plane. `latent_projections.pdf` maps each selected
scenario (lowest and highest overlap, balanced and highly unbalanced
composition, four graph assignments) to two planes.

> **Warning 1**
>
> A J-type bulk draw is
> Y\mid\boldsymbol{p}\sim\mathcal{N}\bigl(\boldsymbol{\mu}(\boldsymbol{p}),\boldsymbol{\Sigma}(\boldsymbol{p})\bigr)
> with \boldsymbol{\mu}(\boldsymbol{p})=\sum_j
> p_j\boldsymbol{\mu}\_{\cdot j} and
> \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
> p_j^2\boldsymbol{\Sigma}\_j. Every observation is a weighted sum of
> **all** types
> ([`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)).
> A finite mixture instead draws a latent label and then one component,
> f(\boldsymbol{y})=\sum_k\pi_k\\\phi(\boldsymbol{y};\boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k),
> with mixing weights \boldsymbol{p} (not \boldsymbol{p}^2). `MclustDR`
> targets that mixture. Its \pi_k are not cell-type proportions, and its
> component densities are not the convolution law used for
> deconvolution.

#### Principal components of the purified profiles

Principal components of the pooled purified draws, coloured by type,
show what the **means** alone separate: types 1 and 2, whose means have
cosine 0.9, overlap on the first two components, while type 3 stands
apart.

#### Supervised `MclustDR` plane

The second panel fits `mclust::MclustDA(..., modelType = "EDDA")` to the
labelled purified draws, then
[`mclust::MclustDR()`](https://mclust-org.github.io/mclust/reference/MclustDR.html)
([Fraley et al. 2026](#ref-R-mclust); [Scrucca
2010](#ref-scruccaDimensionReductionModelbased2010a),
[2015](#ref-scruccaGraphicalToolsModelbased2015)) ([`MclustDR`
reference](https://mclust-org.github.io/mclust/reference/MclustDR.html)).
Its kernel combines the separation of the class means with the
differences of the class covariances, which is the signal DeCovarT
exploits. Its generalised eigenvalues measure class-relevant
variability, not ordinary PCA variance. The panels with independent
factor analyses (Thomson scores) and mixtures of common factor analysers
were dropped: they described each type, or an unlabelled mixture, but
not the contrast between types.

### Solver runtime and memory

`runtime.pdf` and `memory.pdf` use the same composition order as the
other performance books (balanced, then moderately unbalanced, then
highly unbalanced): **one page per** (). Each page facets the three
MixSim overlap levels, puts CT1/CT2 topology on the x-axis (`SF/SF`,
`SF/SBM`, `SBM/SF`, `SBM/SBM`), and maps solver to colour, fill, and
grouping. The geometry is a vertical ggdist half-eye (raincloud) of the
per-sample `optimisation$elapsed_sec` or `memory_bytes` (converted to
MiB), with a left-side `geom_rug()` of the raw Monte Carlo draws and a
log10 y-axis.

### Likelihood geometry and solver behaviour

> **Note**
>
> The numbers in this section come from the archived run of
> `output/archives/fig03/` (N=50, Wald intervals only), which still
> included Marquardt–Levenberg and the earlier L-BFGS-B that
> renormalised a box-constrained fit. They will be replaced by those of
> the production run; the mechanisms they illustrate do not depend on
> the solver list.

Cell types 1 and 2 share cosine 0.9, so the mean map
\boldsymbol{\mu}\boldsymbol{p} is close to singular in the (p_1,p_2)
plane: type 3 is the only well-separated mean. MixSim \bar\omega then
sets how much of the remaining signal lives in
\\\boldsymbol{\Sigma}\_j\\. The pack in `output/fig03/mle_explanation/`
(descriptors, ILR slices of \ell at
\boldsymbol{y}=\boldsymbol{\mu}\boldsymbol{p}^{\star}, expected-Fisher
SEs, Newton/Marquardt forests on **balanced** compositions, and a
quantile-sampled raincloud across all H^{\star}) is the numerical
reading of that geometry. As in the bivariate toy, three layers must be
kept apart: the convolution \ell, the surrogate each solver optimises,
and the chart / KKT / stopping rule.

#### Surrogates versus the convolution

`lsei` (`limSolve`) is again a convex QP for
\\\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}\\\_2^2 on the simplex
([Soetaert et al. 2026](#ref-R-limSolve); [Gong and Szustakowski
2013](#ref-gongDeconRNASeqStatisticalFramework2013); [Dessole et al.
2023](#ref-dessoleLawsonHansonAlgorithmDeviation2023); [Boyd et al.
2004](#ref-boydConvexOptimization2004)). With two nearly collinear
signature columns the QP is under-determined in the (p_1,p_2) direction;
the equality \sum p_j=1 plus bounds then put mass on a face. `CIBERSORT`
replaces squares by a \nu-SVR hinge ([Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)). The dual is a convex
quadratic programme solved by sequential working-set methods ([Fan et
al. 2005](#ref-fanWorkingSetSelection2005)). The wrapper does **not**
constrain coefficients to the simplex during the fit: negatives are
clipped to zero and the remainder is renormalised. Clipping is a
Euclidean projection onto the positive orthant, not a constrained MLE,
so a support vector that overshoots becomes a vertex of \Delta^{2}.
Neither method sees p_j^2\boldsymbol{\Sigma}\_j.

DeCovarT solvers maximise the Gaussian convolution. `Newton-Raphson` and
`Marquardt-Levenberg` do so in ILR, with a full Hessian versus Marquardt
damping plus RDM ([Marquardt
1963](#ref-marquardtAlgorithmLeastSquaresEstimation1963); [Commenges et
al. 2006](#ref-commengesNewtonLikeAlgorithmLikelihood2006); [Philipps et
al. 2021](#ref-philippsRobustEfficientOptimization2021),
[2023](#ref-R-marqLevAlg)). `L-BFGS-B` uses a limited-memory inverse
Hessian on the closure chart \boldsymbol p=\boldsymbol w/\mathbf
1^{\top}\boldsymbol w, \boldsymbol w\ge0 ([Byrd et al.
1995](#ref-byrdLimitedMemoryAlgorithm1995); [Zhu et al.
1997](#ref-zhuAlgorithm778LBFGSB1997)). The same caveats as in [the
bivariate
interpretation](https://bastienchassagnol.github.io/DeCovarT/articles/fig02-bivariate-toy.html#sec-bivariate-findings)
apply: Boyd’s local-quadratic Newton rate needs a locally strongly
concave \ell; RDM stops when further ILR steps are smaller than
statistical error; the archived box-constrained BFGS, which renormalised
a free-scale fit, could sit on a face that is not a stationary point of
\ell (the closure chart now used cannot). All of this is still N=1 per
bulk ([finite-sample
MLE](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-finite-sample)).

#### Balanced compositions (forest + Fisher)

The Newton/Marquardt forest in the pack is restricted to H^{\star}=1
(p^{\star}=(1/3,1/3,1/3)). On that slice the two second-order solvers
are interchangeable. At low overlap (`V1`, both graphs scale-free) RMSE
for p_1 is 0.113, bias 0.012, empirical SD 0.113, coverage 0.99. At high
overlap (`V7`) RMSE rises only to 0.124 and coverage stays 0.99. Median
RMSE for p_1 across the four graph assignments is 0.122 (low
\bar\omega), 0.134 (moderate), 0.147 (high); type 3, whose mean is not
collinear, stays easier at low overlap (median RMSE 0.064) and catches
up as overlap grows (0.130). Expected-Fisher ILR standard errors in
`theoretical_wald_se.csv` tell the same story: \mathrm{SE}\_{\rho_1}
goes from 0.63 (`V1`) to 0.97 (`V7`), and the condition number of the
ILR information drops from 6.8 to 1.7 — not because the problem becomes
better posed, but because *both* eigenvalues shrink and the ellipsoid
becomes a large, round disk. Graph assignment (scale-free versus cluster
SBM) moves those SEs by a smaller amount than \bar\omega does. Topology
is visible; **average overlap is the dominant knob**, matching the
design claim, for *balanced* p^{\star} and for *second-order convolution
MLEs*.

#### Boundary pile-up (raincloud, all H^{\star})

The sampled raincloud (`raincloud_sampled.csv`) includes moderately and
highly unbalanced Shannon targets. Among non-missing point estimates,
76\\ of `CIBERSORT` values lie at 0 or 1. At the scenario-wise median
(`sample_role = q50`), type 3 is exactly zero in 35 of 36 scenarios
(mean \hat p_3=0.0008). That is the clip-and-renormalise step, not a
convolution mode on a face. `lsei` sits on a bound in 42\\ of sampled
cells overall and in 58\\ of highly unbalanced cells; `LBFGS` in 26\\
and 36\\; Marquardt in 10\\ and 17\\; Newton in 4\\ and 12\\. On
balanced raincloud rows Newton and Marquardt never hit \\0,1\\; LSEI
still does (21\\). Mean-only and SVR methods convert a collinear mean
map into a **hard** boundary estimate; ILR convolution methods convert
it into a **soft** interior mode, until p^{\star} itself is a rare-type
vertex.

Median absolute error of the q\_{50} p_1 fit across scenarios is 0.048
(Newton), 0.100 (`lsei`), 0.114 (`LBFGS`), 0.142 (Marquardt), 0.159
(`CIBERSORT`). The Marquardt figure is not a contradiction of the
balanced forest: on the highly unbalanced scenario `V3`
(p^{\star}\approx(0.980,0.010,0.010)) the q\_{50} estimates are 0.92
(Newton), 0.93 (`lsei`), 1.00 (`CIBERSORT`), 0.82 (`LBFGS`), and 0.78
(Marquardt). LSEI is closer to the vertex because the QP *wants* a
vertex; Marquardt is pulled toward the barycentre because RDM refuses a
long ILR walk, the same mechanism as `B811` in the bivariate pack.
Covariance information therefore helps only if the solver actually
follows the convolution ridge to the rare-type coordinate. A mean-only
vertex can beat a damped interior MLE on raw p_1 error when the truth is
already a vertex; it will still report \hat p_3=0 and cannot support a
chi-bar-square test that the type is absent ([boundary
inference](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-boundary)).

#### Ridges in the ILR slices

`ilr_loglik_profiles.csv` slices \ell(\rho_1,\rho_2) at
\boldsymbol{y}=\boldsymbol{\mu}\boldsymbol{p}^{\star}. High \bar\omega
produces ridges along the type-1/type-2 contrast (the collinear mean
direction): many (\rho_1,\rho_2) give almost the same convolution
Gaussian. Newton follows the ridge until the Hessian in the remaining
well-identified direction (type 3 versus \\1,2\\) stops it. Marquardt
stops earlier when RDM sees no inferential gain. `L-BFGS-B`, stepping in
p with only box bounds, can leave the simplex during the line search and
land near a face after p/\sum p, which is the raincloud’s 26\\ bound
rate. None of these paths is the Godambe sandwich or a Firth penalty
([Firth 1993](#ref-firthBiasReductionMaximum1993)); the Firth fit is now
benchmarked as its own method
([Firth](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-firth)).
With N=1, the honest frequentist statement for a rare type is still the
restricted LRT / parametric bootstrap of the MLE article, not a Wald
ellipse from
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md).

### See also

- Bivariate toy:
  [§2.1](https://bastienchassagnol.github.io/DeCovarT/articles/fig02-bivariate-toy.md)
- Moment generator: [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.md)
- Overlap and SPD distances: [Distances, overlap, and covariance
  information](https://bastienchassagnol.github.io/DeCovarT/articles/theory-distance-covariance.md)
- Feature-selection on this design: [Appendix
  S6](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S6-feature-selection.md)
- Identifiability of collinear means: [Appendix
  S1](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.md)

### References

Boyd, Stephen, Stephen P. Boyd, and Lieven Vandenberghe. 2004. *Convex
Optimization*. 1st edn. Cambridge University Press; Cambridge University
Press. <https://doi.org/10.1017/cbo9780511804441>.

Byrd, Richard H., Peihuang Lu, Jorge Nocedal, and Ciyou Zhu. 1995. ‘A
Limited Memory Algorithm for Bound Constrained Optimization’. *SIAM
Journal on Scientific Computing* 16 (5): 1190–208.
<https://doi.org/10.1137/0916069>.

Commenges, Daniel, Helene Jacqmin-Gadda, Cecile Proust, and Jeremie
Guedj. 2006. *A Newton-Like Algorithm for Likelihood Maximization: The
Robust-Variance Scoring Algorithm*. arXiv.
<https://doi.org/10.48550/arxiv.math/0610402>.

Dessole, Monica, Marco Dell’Orto, and Fabio Marcuzzi. 2023. ‘The
Lawson-Hanson Algorithm with Deviation Maximization: Finite Convergence
and Sparse Recovery’. *Numerical Linear Algebra with Applications* 30
(5): e2490. <https://doi.org/10.1002/nla.2490>.

Fan, Rong-En, Pai-Hsuen Chen, and Chih-Jen Lin. 2005. ‘Working Set
Selection Using Second Order Information for Training Support Vector
Machines’. *Journal of Machine Learning Research* 6 (63): 1889–918.

Firth, David. 1993. ‘Bias Reduction of Maximum Likelihood Estimates’.
*Biometrika* 80 (1): 27–38. <https://doi.org/10.1093/biomet/80.1.27>.

Fraley, Chris, Adrian E. Raftery, and Luca Scrucca. 2026. *Mclust:
Gaussian Mixture Modelling for Model-Based Clustering, Classification,
and Density Estimation*. <https://mclust-org.github.io/mclust/>.

Gong, Ting, and Joseph D. Szustakowski. 2013. ‘DeconRNASeq: A
Statistical Framework for Deconvolution of Heterogeneous Tissue Samples
Based on mRNA-Seq Data’. *Bioinformatics (Oxford, England)* 29.
<https://doi.org/10.1093/bioinformatics/btt090>.

Marquardt, Donald W. 1963. ‘An Algorithm for Least-Squares Estimation of
Nonlinear Parameters’. *Journal of the Society for Industrial and
Applied Mathematics* 11. <https://doi.org/10.1137/0111030>.

Melnykov, Volodymyr, Wei-Chen Chen, and Ranjan Maitra. 2012. ‘MixSim: An
R Package for Simulating Data to Study Performance of Clustering
Algorithms’. *Journal of Statistical Software* 51.
<https://doi.org/10.18637/jss.v051.i12>.

Morris, Tim P., Ian R. White, and Michael J. Crowther. 2019. ‘Using
Simulation Studies to Evaluate Statistical Methods’. *Statistics in
Medicine* 38 (11): 2074–102. <https://doi.org/10.1002/sim.8086>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Philipps, Viviane, Boris P. Hejblum, Mélanie Prague, Daniel Commenges,
and Cécile Proust-Lima. 2021. ‘Robust and Efficient Optimization Using a
Marquardt-Levenberg Algorithm with R Package marqLevAlg’. *The R
Journal* 13. <https://doi.org/10.32614/rj-2021-089>.

Philipps, Viviane, Cecile Proust-Lima, Melanie Prague, Boris Hejblum,
Daniel Commenges, and Amadou Diakite. 2023. *marqLevAlg: A Parallelized
General-Purpose Optimization Based on Marquardt-Levenberg Algorithm*.

Riani, Marco, Andrea Cerioli, Domenico Perrotta, and Francesca Torti.
2015. ‘Simulating Mixtures of Multivariate Data with Fixed Cluster
Overlap in FSDA Library’. *Advances in Data Analysis and Classification*
9 (4): 461–81. <https://doi.org/10.1007/s11634-015-0223-9>.

Scrucca, Luca. 2010. ‘Dimension Reduction for Model-Based Clustering’.
*Statistics and Computing* 20 (4): 471–84.
<https://doi.org/10.1007/s11222-009-9138-7>.

Scrucca, Luca. 2015. *Graphical Tools for Model-Based Mixture
Discriminant Analysis*. Https://arxiv.org/abs/1508.01695v1.
<https://doi.org/10.1007/s11634-013-0147-1>.

Soetaert, Karline, Karel Van den Meersche, and Dick van Oevelen. 2026.
*limSolve: Solving Linear Inverse Models*.

Zhu, Ciyou, Richard H. Byrd, Peihuang Lu, and Jorge Nocedal. 1997.
‘Algorithm 778: L-BFGS-B: Fortran Subroutines for Large-Scale
Bound-Constrained Optimization’. *ACM Transactions on Mathematical
Software* 23 (4): 550–60. <https://doi.org/10.1145/279232.279236>.
