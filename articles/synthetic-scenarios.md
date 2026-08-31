# Simulating synthetic pseudo-bulk mixtures for benchmarking

## Overview

This vignette documents the synthetic first- and second-order moment
generator
[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
and how its outputs feed bulk mixture simulation. Mean profiles follow
an AutoGeneS-inspired construction ([Aliee and Theis
2021](#ref-alieeAutoGeneSAutomaticGene2021)) in which a single dial \rho
(`target_cosine`) sets pairwise collinearity between cell-type columns
of \boldsymbol{\mu}, while a fixed scale s (`mean_scale`) sets column
norms. Prefer varying \rho across scenarios when precision weights
already control second-order dependence.

Each cell type carries its own graph-constrained signed precision
\boldsymbol{\Omega}\_j, mapped to
\boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}\_j^{-1}. An undirected
Gaussian Markov simulation separates four layers
([Equation 12](#eq-ggm-pipeline)): graph G_j, signed weights W(G_j), SPD
precision \Omega_j, and latent Gaussians (an optional observation layer
can add Poisson–log-normal or zero-inflated counts). The **graph** and
the **precision** must not be conflated: a binary adjacency is almost
never itself positive definite. DeCovarT draws G_j via
[`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md),
completes \Omega_j\succ 0 with
[`build_normalised_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_normalised_precision.md),
and stacks the inverted slices with
[`build_covariance_array_from_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_covariance_array_from_precision.md).
[Section 4](#sec-ggm-networks) develops the topology, weight, and SPD
design in detail.

## Pipeline API

The exported entry point orchestrates internal helpers and returns
moments that can be passed to
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md).
For end-to-end benchmarking, see [Manuscript synthetic simulation
scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-manuscript-scenarios.html#sec-simulation-benchmark);
scenario grids live in `scripts/configure_bivariate_toy_scenarios.R`.
Direct solver calls use
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart LR
  A["Graph G"] --> B["Weights W(G)"]
  B --> C["Precision Ω ≻ 0"]
  C --> D["Σ = Ω^{-1}"]
  E["Means μ_j\n(AutoGeneS-style ρ)"] --> F["simulate_hierarchical_grn_moments()"]
  D --> F
  F --> G["simulate_bulk_mixture()"]
  G --> H["deconvolute_ratios()"]
```

Figure 1: Pipeline from synthetic moment generation to deconvolution.

``` r

library(DeCovarT)

set.seed(42)
moments <- simulate_hierarchical_grn_moments(
  n_genes = 40L,
  n_celltypes = 3L,
  mean_scale = 10,
  target_cosine = 0.1,
  precision_shift = 0.1,
  precision_scale = 0.3,
  prop_inhibitory = 0.5,
  graph_model = "scale_free"
)

moments$objectives
#> $mean_abs_cosine
#> [1] 0.3315284
#> 
#> $sum_euclidean_distance
#> [1] 34.68786

bulk <- simulate_bulk_mixture(
  signature_matrix = moments$mean_profiles,
  Sigma = moments$covariance_matrices,
  p = c(0.5, 0.3, 0.2),
  n = 20
)
dim(bulk$Y)
#> [1] 40 20
dim(bulk$latent_profiles) # G x J x N latent draws; input mu is shared
#> [1] 40  3 20
```

## Generate mean expression profiles

### Mean signature: `generate_mean_signature_matrix()`

The construction uses one symbol per quantity:

- G: number of genes
- J: number of cell types, indexed by j=1,\ldots,J
- \boldsymbol{\mu}\_{\cdot j}\in\mathbb{R}^{G}: column j of the mean
  signature
- s: column Euclidean scale (`mean_scale`; default 10 in the
  nine-scenario grid)
- \rho\in\[0,1\]: target pairwise cosine (`target_cosine`)
- \boldsymbol{u}=G^{-1/2}\mathbf{1}\_{G}: shared unit direction
- \boldsymbol{v}\_{j}: \ell_2-normalised indicator of the gene block
  assigned to type j (a disjoint partition of \\1,\ldots,G\\)

The private directions are mutually orthogonal:
\boldsymbol{v}\_{j}^{\mathsf{T}}\boldsymbol{v}\_{k}=0 for j\neq k. Each
column is the normalised blend

\tilde{\boldsymbol{\mu}}\_{\cdot j} = \sqrt{\rho}\\\boldsymbol{u}
+\sqrt{1-\rho}\\\boldsymbol{v}\_{j}, \qquad \boldsymbol{\mu}\_{\cdot j}
= s\\ \frac{\tilde{\boldsymbol{\mu}}\_{\cdot j}}{
\\\tilde{\boldsymbol{\mu}}\_{\cdot j}\\\_2 }, \qquad j=1,\ldots,J.
\tag{1}

Stacking the un-normalised columns yields a **rank-one shared term plus
a block-sparse private term**:

\tilde{\boldsymbol{M}} = \bigl\[ \tilde{\boldsymbol{\mu}}\_{\cdot 1}
\mid \cdots \mid \tilde{\boldsymbol{\mu}}\_{\cdot J} \bigr\] =
\sqrt{\rho}\\ \boldsymbol{u}\mathbf{1}\_{J}^{\mathsf{T}} +
\sqrt{1-\rho}\\ \mathbf{V}, \tag{2}

where
\mathbf{V}=\bigl\[\boldsymbol{v}\_{1}\mid\cdots\mid\boldsymbol{v}\_{J}\bigr\]
has **disjoint supports** across columns (each \boldsymbol{v}\_j is
non-zero only on block j). The shared part is rank one; the private part
is block-diagonal in the gene-partition sense. After column
normalisation and scaling,

\boldsymbol{\mu} = s\\\tilde{\boldsymbol{M}}\\ \mathbf{D}^{-1}, \qquad
\mathbf{D}=\mathrm{diag}\bigl( \\\tilde{\boldsymbol{\mu}}\_{\cdot
1}\\\_2, \ldots, \\\tilde{\boldsymbol{\mu}}\_{\cdot J}\\\_2 \bigr),
\tag{3}

so \boldsymbol{\mu} is a low-rank-plus-block construction followed by
per-column calibration. This matches synthetic signature designs that
combine a common baseline with type-specific marker programmes ([Aliee
and Theis 2021](#ref-alieeAutoGeneSAutomaticGene2021); [Schelker et al.
2017](#ref-schelkerEstimationImmuneCell2017); [Ba et al.
2026](#ref-baWhenLessNot2026)).

Expanding the un-normalised inner product for j\neq k gives the exact
identity

\tilde{\boldsymbol{\mu}}\_{\cdot j}^{\mathsf{T}}
\tilde{\boldsymbol{\mu}}\_{\cdot k} = \rho + \sqrt{\rho(1-\rho)}\\
\bigl( \boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{j} +
\boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{k} \bigr), \tag{4}

because \boldsymbol{u}^{\mathsf{T}}\boldsymbol{u}=1 and
\boldsymbol{v}\_{j}^{\mathsf{T}}\boldsymbol{v}\_{k}=0. Block indicators
are **not** orthogonal to the shared direction: if block j contains L_j
genes then

\boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{j} = \sqrt{\frac{L_j}{G}},
\tag{5}

so [Equation 4](#eq-mean-inner) becomes

\tilde{\boldsymbol{\mu}}\_{\cdot j}^{\mathsf{T}}
\tilde{\boldsymbol{\mu}}\_{\cdot k} = \rho + \sqrt{\rho(1-\rho)}\\
c\_{jk}, \qquad c\_{jk} = \sqrt{\tfrac{L_j}{G}} + \sqrt{\tfrac{L_k}{G}}.
\tag{6}

Column-wise renormalisation in [Equation 1](#eq-mean-cosine) maps this
to a realised cosine that is **not** exactly \rho for finite G and J:

\cos\bigl( \boldsymbol{\mu}\_{\cdot j}, \boldsymbol{\mu}\_{\cdot k}
\bigr) = \frac{ \rho + \sqrt{\rho(1-\rho)}\\c\_{jk} }{ 1 +
\sqrt{\rho(1-\rho)}\\c\_{jk} }, \qquad j\neq k. \tag{7}

#### Asymptotic cosine control

> **Note 1: The realised cosine is systematically larger than \rho**
>
> For \rho\in(0,1) and finite J, the shared direction \boldsymbol{u} is
> **not** orthogonal to the private blocks, so
> \boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}\_{j}=\sqrt{L_j/G}\>0. That
> strictly positive cross term inflates every pairwise cosine relative
> to the dial \rho. Increasing J (equal blocks, G/J bounded away from
> zero) drives the bias to zero like J^{-1/2}, but the sign of the bias
> does not flip. AutoGeneS uses a related shared-plus-private geometry
> when it penalises inter-type correlation ([Aliee and Theis
> 2021](#ref-alieeAutoGeneSAutomaticGene2021)); the closed-form
> O(J^{-1/2}) remainder below is the finite-(G,J) identity for
> [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md).

> **Note**
>
> **Theorem 1 (Asymptotic cosine control)** For equal blocks (L_j\approx
> G/J) we have c\_{jk}\approx 2/\sqrt{J}. Writing the bias relative to
> the target,
>
> \cos\bigl( \boldsymbol{\mu}\_{\cdot j}, \boldsymbol{\mu}\_{\cdot k}
> \bigr) - \rho = \frac{ (1-\rho)\sqrt{\rho(1-\rho)}\\c\_{jk} }{ 1 +
> \sqrt{\rho(1-\rho)}\\c\_{jk} } = O\\\bigl(J^{-1/2}\bigr) = o(1)
> \quad\text{as }J\to\infty \tag{8}
>
> (with G/J bounded away from zero). Equivalently,
>
> \lim\_{J\to\infty} \cos\bigl( \boldsymbol{\mu}\_{\cdot j},
> \boldsymbol{\mu}\_{\cdot k} \bigr) = \rho, \tag{9}
>
> so the cross terms in [Equation 4](#eq-mean-inner) become negligible
> relative to \rho and the realised pairwise cosine tracks the dial
> asymptotically. The scale s sets column norms without changing angles:
> \\\boldsymbol{\mu}\_{\cdot j}-\boldsymbol{\mu}\_{\cdot k}\\\_2\propto
> s at fixed \rho. We keep s=10 across scenarios and only vary \rho.

[Figure 3](#fig-cosine-convergence) illustrates the finite-sample
approach to \rho=0.5 on a 3\times 3 grid of (J,G).
[Figure 4](#fig-cosine-bias-asymptotic) isolates the **positive**
remainder as a function of J: realised cosine always sits above the
target (except at the trivial endpoints \rho\in\\0,1\\ up to
discretisation), and the gap shrinks as more cell types split the gene
panel into smaller private blocks.

[Figure 2](#fig-mean-signature-schematic) gives a toy schematic for J=3
cell types: a uniform shared component, three orthogonal private blocks,
and the resulting mean directions.

![](figures/fig_mean_signature_shared_private.png)

Figure 2: Shared-plus-private mean signature for J=3 cell types. **A.**
Gene-by-cell-type heatmap: each type is high on its private gene block
and shares a grey baseline elsewhere. **B.** The shared unit direction
\boldsymbol{u} blends with mutually orthogonal private indicators
\boldsymbol{v}\_j at weights \sqrt{\rho} and \sqrt{1-\rho}; columns are
then \ell_2-normalised and scaled to \boldsymbol{\mu}\_{\cdot j}.
Regenerated for this vignette (two-panel layout; the previous toy PCA
panel is omitted).

#### Visualising the cosine construction

![](synthetic-scenarios_files/figure-html/fig-cosine-convergence-1.png)

Figure 3: Realised mean absolute pairwise cosine vs target \rho=0.5 on a
3\times 3 grid: J\in\\3,10,20\\ nested with G\in\\20,100,400\\ (equal
gene blocks). Thick black line: target; coloured bands: harder (red,
small J) to easier (green, large J); arrows: residual to the target.

![](synthetic-scenarios_files/figure-html/fig-cosine-bias-asymptotic-1.png)

Figure 4: Positive finite-J bias of the realised pairwise cosine at
target \rho=0.5 and G=400 genes (equal blocks). Points: empirical mean
absolute cosine from
[`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md).
Solid curve: closed form [Equation 7](#eq-realised-cosine). Dashed line:
target \rho. The remainder is strictly positive and shrinks like
J^{-1/2} ([Theorem 1](#thm-cosine-asymptotics)).

[Figure 5](#fig-cosine-geometry) visualises the same three operating
points (\rho\in\\0,\\0.3,\\1\\) for G=2 genes and J=2 cell types. Cell
type 1 is fixed along the positive gene-1 axis; only the **angular**
position of cell type 2 is varied so that \rho=0 is exactly orthogonal.
Independent Gaussian noise is added *only in this vignette chunk* around
each centroid (N=20 i.i.d. draws per population) to emulate within-type
dispersion under an independence assumption, without altering the
mean-signature generator itself. Coloured arrows mark each type’s mean
direction; the shaded wedge is the angle between the two means.

This illustration is inspired by panel B in the AutoGeneS paper ([Aliee
and Theis 2021](#ref-alieeAutoGeneSAutomaticGene2021)), but focuses here
on the cosine-angle control induced by [Equation 1](#eq-mean-cosine)
(without modelling centroid-distance optimisation). With only two genes
the cross terms in [Equation 4](#eq-mean-inner) are large, so the
realised cosine from the generator is a warped but strictly increasing
map of \rho; the panels below enforce the target angle geometrically
after generation.

![](synthetic-scenarios_files/figure-html/fig-cosine-geometry-1.png)

Figure 5: Cosine-control toy for J=2 cell types and G=2 genes. Cell type
1 is fixed on the positive gene-1 axis; only the angular position of
cell type 2 varies (so \rho=0 is orthogonal). Coloured arrows: mean
directions; shaded wedge: angle between means. Labels report the
realised cosine. Points are N=20 i.i.d. Gaussian draws around each
centroid (stars). Styling inspired by panel B of `AutoGeneS` (Aliee and
Theis ([2021](#ref-alieeAutoGeneSAutomaticGene2021))).

> **Note 2: Related signature constructions**
>
> Blending a common baseline with type-specific markers is standard in
> deconvolution benchmarking. **`AutoGeneS`** ([Aliee and Theis
> 2021](#ref-alieeAutoGeneSAutomaticGene2021)) selects genes by jointly
> minimising inter-type correlation and maximising centroid distance;
> its simulations use a similar shared-plus-private logic to control
> signature geometry. Schelker and colleagues ([Schelker et al.
> 2017](#ref-schelkerEstimationImmuneCell2017)) build tumour-derived
> reference gene expression profiles from single-cell RNA-seq clusters
> and marker genes, then simulate bulk mixtures by combining those
> columns—empirically a shared baseline plus distinct marker programmes
> per immune cell type.
>
> The **`DICEPro`** framework ([Ba et al. 2026](#ref-baWhenLessNot2026))
> targets incomplete reference matrices in supervised deconvolution. Its
> simulation engine generates synthetic signature matrices from
> multivariate Gaussian (or Poisson) models with a prescribed
> correlation structure; a **`bloc=TRUE`** option enforces **block
> sparsity** by zeroing expression outside type-specific gene
> blocks—closely analogous to the private indicators \boldsymbol{v}\_j
> in [Equation 2](#eq-lr-block). DICEPro then optimises reference
> signatures when cell types are missing from the matrix, complementing
> the controlled \rho-dial we use here for method comparison under known
> ground truth.

### Objectives: `compute_mean_profile_objectives()`

For columns of \boldsymbol{\mu},

\overline{\lvert\cos\rvert} = \frac{2}{J(J-1)} \sum\_{1\le j\<k\le J}
\Biggl\| \frac{ \boldsymbol{\mu}\_{\cdot j}^{\mathsf{T}}
\boldsymbol{\mu}\_{\cdot k} }{ \\\boldsymbol{\mu}\_{\cdot j}\\\_2\\
\\\boldsymbol{\mu}\_{\cdot k}\\\_2 } \Biggr\|,

D\_{\mathrm{Euc}} = \sum\_{1\le j\<k\le J} \\\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\_{\cdot k}\\\_2.

AutoGeneS treats the first as a quantity to minimise and the second as a
quantity to maximise ([Aliee and Theis
2021](#ref-alieeAutoGeneSAutomaticGene2021)). Cosine is preferred to
Pearson correlation so that collinear but unequally scaled profiles
remain penalised.

> **Note 3: Perspectives: a single geometric score for
> \boldsymbol{\mu}**
>
> The present API exposes two AutoGeneS-style objectives. A natural
> refinement is to replace them by one scalar that jointly rewards large
> column norms (and hence Euclidean separation) and penalises alignment.
>
> Two candidates are immediate. The Gram determinant
>
> V(\boldsymbol{\mu}) = \sqrt{ \det\bigl(
> \boldsymbol{\mu}^{\mathsf{T}}\boldsymbol{\mu} \bigr) } \tag{10}
>
> is the volume of the parallelepiped spanned by the columns of
> \boldsymbol{\mu}. It vanishes under exact collinearity and grows when
> columns lengthen or become more orthogonal—precisely the desired MOO
> trade-off in one number.
>
> Equivalently, the reciprocal condition number
>
> \kappa_2(\boldsymbol{\mu})^{-1} =
> \frac{\sigma\_{\min}(\boldsymbol{\mu})}{\sigma\_{\max}(\boldsymbol{\mu})},
> \tag{11}
>
> available in R via `kappa(mean_profiles, exact = TRUE)` (or
> `1 / kappa(...)`), measures how ill-posed linear recovery of
> \boldsymbol{p} from
> \boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p} is. Maximising
> \kappa_2(\boldsymbol{\mu})^{-1} (or minimising
> \kappa_2(\boldsymbol{\mu})) collapses cosine and distance into a
> deconvolution-centric criterion without a Pareto front.
>
> A practical next step is to expose an optional
> `mean_objective = c("autogenes", "volume", "condition")` in
> [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md),
> keep the constructive (\rho,s) generator for controllable scenarios,
> and report the chosen scalar alongside the existing pairwise
> diagnostics.

## Undirected Gaussian Markov network generation

G \longrightarrow W(G) \longrightarrow \Omega \succ 0 \longrightarrow
X_i \sim \mathcal{N}\_p(\mu_i,\Omega^{-1}) \tag{12}

Here G is an **undirected** graph, W(G) a symmetric weighted matrix with
the same off-diagonal support, \Omega a strictly positive-definite
precision, and X_i an expression profile. This section covers
**undirected** Gaussian Markov / precision-graph simulation only.
[`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md)
draws G;
[`build_normalised_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_normalised_precision.md)
completes it to \Omega\succ 0 by an affine spectral shift;
[`build_covariance_array_from_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_covariance_array_from_precision.md)
maps each slice \boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}\_j^{-1} (cell
types need not share one network).

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart LR
  G["Undirected graph G"] --> W["Signed weights W(G)"]
  W --> Omega["SPD precision Ω"]
  Omega --> Latent["Latent Gaussian X"]
  Mu["Mean design μ"] --> Latent
  Latent --> Obs["Observation layer"]
```

Figure 6: Undirected simulation pipeline from topology to observations
(optional observation layer dashed).

### Graph structure generators

\rho\_{jk\mid -\\j,k\\} =
-\frac{\Omega\_{jk}}{\sqrt{\Omega\_{jj}\Omega\_{kk}}} \tag{13}

[Table 1](#tbl-topologies) compares the main **undirected** topology
families used in GGM benchmarks, with R entry points.
[Figure 7](#fig-topologies) sketches the six families most often used in
random-like network structure generation.

| Structure | Definition and controls | R entry points | Advantages | Limitations |
|----|----|----|----|----|
| **Erdős–Rényi** | Each edge present independently with probability q; q=d\_{\mathrm{avg}}/(p-1) targets average degree | `huge::huge.generator(graph="random")` ([Jiang et al. 2026](#ref-R-huge)); `BDgraph::graph.sim(graph="random")` ([Mohammadi and Wit 2025](#ref-R-BDgraph)); DeCovarT `graph_model = "erdos_renyi"` | Exact sparsity control in expectation; few nuisance parameters | Narrow degree distribution; little modular organisation |
| **Band / AR** | Edge j–k if 1\le\lvert j-k\rvert\le b | `huge.generator(graph="band")` ([Jiang et al. 2026](#ref-R-huge)); `BDgraph::bdgraph.sim(graph="AR(1)")` / `"AR(2)"` ([Mohammadi and Wit 2025](#ref-R-BDgraph)) | Local dependence and controlled max degree | Artificial gene ordering; no hubs or modules |
| **Hub / star** | g groups; one centre linked to remaining members | `huge.generator(graph="hub")` ([Jiang et al. 2026](#ref-R-huge)); `graph.sim(graph="hub")` / `"star"` ([Mohammadi and Wit 2025](#ref-R-BDgraph)); DeCovarT `graph_model = "hub"` | Stress-tests high-degree regulators | Pure stars omit cross-talk |
| **Scale-free** | Preferential attachment (e.g. Barabási–Albert with m edges per new node) | `huge.generator(graph="scale-free")` ([Jiang et al. 2026](#ref-R-huge)); DeCovarT `graph_model = "scale_free"` | Heterogeneous degrees | Weak modularity under plain BA growth |
| **Cluster / SBM** | Within-block edge probability exceeds between-block | `huge.generator(graph="cluster")` ([Jiang et al. 2026](#ref-R-huge)); DeCovarT `graph_model = "stochastic_block_model"` | Pathway / module structure | Equal blocks can be unrealistically regular |
| **Small-world / lattice / circle** | High clustering with short paths, or fixed neighbourhoods | DeCovarT `graph_model = "small_world"` (Watts–Strogatz); `BDgraph::graph.sim(graph="smallworld")` ([Mohammadi and Wit 2025](#ref-R-BDgraph)) | Local modules and spatial organisation | Narrow degrees; needs a defensible ordering |

Table 1: Topology families for undirected GGM simulation (R-focused).

![](figures/fig_network_topologies_six.png)

Figure 7: Six undirected topology families used in GGM simulation:
Erdős–Rényi, Barabási–Albert (preferential attachment / scale-free),
stochastic block (cluster), band / AR, star / hub, and small-world
(Watts–Strogatz).

Not every topology is equally useful as a **gene-regulatory** prior. The
three simulation studies in [Table 3](#tbl-studies), together with
classical network biology, suggest the following reading.

| Family | Biological reading | Why it matters for DeCovarT |
|----|----|----|
| **Star / hub** | Master-regulator TF with many targets; one-to-many signalling | Stress-tests recovery of a few high-degree hubs that dominate partial correlations ([Wu and Luo 2022](#ref-wuEstimatingHeterogeneousGene2022); [Federico et al. 2023](#ref-federicoStructureLearningGene2023)) |
| **Scale-free (BA)** | Preferential attachment produces a heavy-tailed degree distribution and several hubs ([Barabási and Albert 1999](#ref-barabasiEmergenceScalingRandom1999)) | Standard high-dimensional GGM stress test used by SILGGM via `huge` ([Zhang et al. 2018](#ref-zhangSILGGMExtensivePackage2018)); useful for testing degree heterogeneity, but not a default empirical model of gene regulation |
| **Cluster / SBM** | Co-regulated modules / pathways with dense within-block links | Matches pathway organisation; BLGGM’s global block partition is SBM-like before finer within-block motifs ([Wu and Luo 2022](#ref-wuEstimatingHeterogeneousGene2022)) |
| **Small-world** | High local clustering with short paths (ring + shortcuts) | Local co-expression neighbourhoods plus long-range regulatory shortcuts; useful when modularity is soft rather than hard blocks ([Federico et al. 2023](#ref-federicoStructureLearningGene2023)) |
| **Band / AR** | Ordered local dependence only | Useful numerical control, but gene indices rarely encode a natural order—limited biological fidelity ([Zhang et al. 2018](#ref-zhangSILGGMExtensivePackage2018)) |
| **Erdős–Rényi** | Homogeneous random wiring | Null / baseline sparsity; little pathway or hub structure ([Zhang et al. 2018](#ref-zhangSILGGMExtensivePackage2018)) |

Table 2: Biological interpretability of undirected topology families.

> **Note 4: The scale-free hypothesis: a useful stress test, not a
> biological default**
>
> The phrase *scale-free network* usually means that the degree
> distribution follows a power law, \Pr(K=k)\propto k^{-\alpha}, at
> least above a lower cut-off. This statistical statement is distinct
> from the Barabási–Albert (BA) growth mechanism. Preferential
> attachment can generate a power law, but observing a heavy-tailed
> degree distribution does not establish preferential attachment; other
> mechanisms can produce similar tails, and variants of preferential
> attachment need not produce a pure power law ([Lima-Mendez and van
> Helden 2009](#ref-lima-mendezPowerfulLawPower2009); [Broido and
> Clauset 2019](#ref-broidoScalefreeNetworksAre2019)).
>
> Whether biological networks are scale-free has therefore remained
> contested. Early claims often relied on approximately straight log–log
> degree plots. Such plots are not goodness-of-fit tests: binning, the
> plotting scale, a small number of hubs, and incomplete or biased
> sampling can all make non-power-law distributions appear linear. The
> underlying graph representation matters as well. For example, pool
> metabolites can become artificial hubs in metabolic graphs, while
> protein-interaction networks depend strongly on the assay and sampling
> scheme ([Lima-Mendez and van Helden
> 2009](#ref-lima-mendezPowerfulLawPower2009)). A power law should
> instead be fitted to the discrete degree data, tested for
> plausibility, and compared on the same support with alternatives such
> as log-normal, exponential, stretched-exponential, and
> power-law-with-cut-off distributions.
>
> The evidence is not uniformly negative. Using adjacency-spectrum
> distributions rather than degree alone, Takahashi and colleagues
> classified eight protein–protein interaction networks as scale-free
> ([Takahashi et al.
> 2012](#ref-takahashiDiscriminatingDifferentClasses2012)). Importantly,
> this was a relative model-selection result: BA-like networks fitted
> better than the two candidate alternatives, Erdős–Rényi and
> Watts–Strogatz networks. It did not show that a power law was an
> adequate absolute model, that preferential attachment generated the
> observed networks, or that unobserved interactions would preserve the
> classification. The authors explicitly noted both the restricted
> candidate set and possible sampling artefacts.
>
> A broader test of 928 real-world network data sets reached a more
> cautious conclusion ([Broido and Clauset
> 2019](#ref-broidoScalefreeNetworksAre2019)). Across all domains, only
> 4% met the strongest scale-free criterion, whereas a log-normal
> distribution fitted as well as or better than a power law for 88% of
> degree distributions. Among the biological networks, 63% showed
> neither direct nor indirect evidence of scale-free structure, although
> 6% met the strongest criterion, principally metabolic networks. These
> proportions depend on the network corpus and the operational
> definition of scale-freeness, but they reject a universal scale-free
> law rather than excluding scale-free structure from every biological
> system.

[Table 3](#tbl-studies) summarises how recent simulation studies combine
topology, precision construction, and mean / observation models.

| Source | Graph structures | Precision construction | Mean / observation | Purpose |
|----|----|----|----|----|
| PLNnetwork ([Chiquet et al. 2018](#ref-chiquetVariationalInferenceSparse2018)) | Erdős–Rényi, preferential attachment, affiliation | \Omega = vG + \operatorname{diag}(\lvert\lambda\_{\min}(vG)\rvert + u) | Latent log-abundance XB; compositional / multinomial counts | Count-network estimators under compositionality |
| SILGGM ([Zhang et al. 2018](#ref-zhangSILGGMExtensivePackage2018)) | Band, hub, Erdős–Rényi, scale-free via `huge` | Package covariance / precision simulation | Centred Gaussian; focus on precision entries | Large-p inferential and computational stress tests |
| BLGGM ([Wu and Luo 2022](#ref-wuEstimatingHeterogeneousGene2022)) | Block mixtures of dense, circle, star, signed modules | Module matrices with structured signed off-diagonals | Cell-type-specific \mu_k, \Omega_k; expression-dependent dropout | Joint clustering, heterogeneous networks, zero inflation |

Table 3: Literature survey of undirected graph → precision → observation
designs.

Several conclusions follow.

1.  **Band** is not intrinsically a random-graph model: its support is
    usually deterministic once the bandwidth is fixed. Hub constructions
    are often partly deterministic as well. Erdős–Rényi, preferential
    attachment, and stochastic block models are genuinely stochastic.
2.  A positive off-diagonal entry in \Omega encodes a **negative**
    partial correlation ([Equation 13](#eq-partial-cor)). Uniformly
    positive precision weights therefore induce uniformly negative
    edgewise partial correlations; random or signed biological weights
    are needed when activation-like and inhibition-like associations
    matter.
3.  **BLGGM is a two-scale hybrid** ([Wu and Luo
    2022](#ref-wuEstimatingHeterogeneousGene2022)), shown in
    [Figure 8](#fig-blggm-two-pronged):
    - **Global partition.** Genes are assigned to blocks (a
      stochastic-block–like cut); between-block edges are sparse.
    - **Local motifs.** Inside each block the topology can be a dense
      clique, a circle, a star / hub, or a signed dense module, so
      distinct blocks can stand for distinct regulatory regimes
      (co-expression programmes, cyclic motifs, master-regulator hubs).
      Signs may differ across modules or cell types.

    The construction remains fully undirected, but it is richer than a
    flat Erdős–Rényi draw. DeCovarT’s exported generators currently
    expose the three families most useful as standalone benchmarks:
    `scale_free`, `stochastic_block_model`, and `small_world`.

![](figures/fig_blggm_two_pronged.png)

Figure 8: Two-pronged BLGGM-style topology: a global stochastic-block
partition of genes, then a local motif (dense module, star, or circle)
inside each block ([Wu and Luo
2022](#ref-wuEstimatingHeterogeneousGene2022)).

### Weight design and random signs

Topology generators return a binary undirected support E. **Signs and
magnitudes** are a separate design layer W(G) in
[Equation 12](#eq-ggm-pipeline): they are assigned *after* (or jointly
with) the skeleton, then completed to an SPD precision. Remember that

\operatorname{sign}(\rho\_{jk\mid -\\j,k\\}) =
-\operatorname{sign}(\Omega\_{jk}) \tag{14}

see [Equation 13](#eq-partial-cor), so a positive precision weight is an
inhibitory partial correlation. [Table 4](#tbl-weights) contrasts three
standard undirected strategies.

| Approach | What is randomised | How \Omega is obtained | When to use |
|----|----|----|----|
| **i.i.d. edge signs** | For each \\j,k\\\in E, draw s\_{jk}\in\\-1,+1\\ (e.g. fair coin or \mathbb{P}(+)=\pi) and magnitude m\_{jk}\>0; set W\_{jk}=W\_{kj}=s\_{jk}m\_{jk} | Support-preserving SPD map of W (spectral shift, diagonal dominance, …) | Simple signed ER / hub / SBM benchmarks; explicit control of the inhibitory fraction |
| **Partial-correlation scaling** | Fill a symmetric matrix R with R\_{jk}=0 off E and R\_{jk}\in(-1,1) (signed) on E; ensure \rho(R)\<1 | \Omega = D^{1/2}(I-R)D^{1/2} ([Equation 17](#eq-partial-scale)) | Direct control of partial-correlation signs and strengths |
| **G-Wishart** (`BDgraph::rgwish()`) | Continuous random \Omega on the fixed graph zeros | \Omega\sim W_G(b,D) already SPD ([Mohammadi and Wit 2025](#ref-R-BDgraph), [2019](#ref-BDgraph2019)) | Heterogeneous random weights without hand-tuned magnitudes; Bayesian-style replicates |

Table 4: Three undirected strategies for signed / random edge weights
after a fixed skeleton.

#### i.i.d. signs on a fixed skeleton

Draw G (ER, band, hub, …). Independently for each undirected edge,

W\_{jk}=W\_{kj}=s\_{jk}\\m\_{jk}, \qquad s\_{jk}\in\\-1,+1\\, \quad
m\_{jk}\>0, \tag{15}

then apply a support-preserving SPD completion (next section). Uniform
positive loadings W=vG are the special case s\_{jk}\equiv +1 used by
PLN-style and many `huge` simulations ([Chiquet et al.
2018](#ref-chiquetVariationalInferenceSparse2018))—all partial
correlations then share the same sign. DeCovarT uses i.i.d. signs with
inhibitory fraction `prop_inhibitory` and magnitude v
(`precision_scale`) baked into \boldsymbol{W} before the spectral SPD
completion.

#### Partial-correlation scaling and why I-R

In a Gaussian Markov model the **matrix of partial correlations** (with
ones on the diagonal) relates to the precision by a diagonal congruence.
Write the off-diagonal partial correlations as a symmetric matrix R with

R\_{jj}=0, \qquad R\_{jk} = \begin{cases} \rho\_{jk\mid -\\j,k\\} &
\\j,k\\\in E,\\ 0 & \\j,k\\\notin E, \end{cases} \tag{16}

and choose magnitudes so that the spectral radius satisfies \rho(R)\<1.
Then

\Omega = D^{1/2}(I-R)D^{1/2} \tag{17}

for a positive diagonal D (often D=I after standardising margins).

#### G-Wishart draws on a fixed graph

Given adjacency G, `BDgraph::rgwish()` samples \Omega\sim W_G(b,D)
already positive definite and with \Omega\_{jk}=0 whenever (j,k)\notin E
([Mohammadi and Wit 2025](#ref-R-BDgraph), [2019](#ref-BDgraph2019)).
Signs and magnitudes arise from the continuous distribution rather than
from an explicit coin-flip layer. This is the most “automatic”
signed-weight generator once the skeleton is fixed; the trade-off is
less direct control of the inhibitory fraction than
[Equation 15](#eq-iid-signs) or [Equation 17](#eq-partial-scale).

### Positive-definiteness strategies

[Table 5](#tbl-pd) contrasts constructions that turn a weighted support
W (or a partial-correlation design) into \Omega\succ 0. DeCovarT’s
[`build_normalised_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_normalised_precision.md)
implements the uniform spectral shift below, with diagonal cushion u
(`precision_shift`):

\boldsymbol{\Omega} = \boldsymbol{W} + \bigl(
\lvert\lambda\_{\min}(\boldsymbol{W})\rvert + u \bigr) \mathbf{I}.

| Construction | Formula | Support preserved? | Role |
|----|----|:--:|----|
| Spectral / PLN loading | \Omega=W+(\lvert\lambda\_{\min}(W)\rvert+u)I | Yes | DeCovarT default; same affine shift as `huge` / PLN ([Chiquet et al. 2018](#ref-chiquetVariationalInferenceSparse2018)) |
| Partial-correlation scaling | \Omega=D^{1/2}(I-R)D^{1/2}, \rho(R)\<1 | Yes | Signs and strengths set on the partial-correlation scale ([Equation 17](#eq-partial-scale)) |
| Graph Laplacian / M-matrix | \Omega=L_G+\varepsilon I | Yes | Off-diagonals non-positive ⇒ all edgewise partial correlations positive |
| G-Wishart | \Omega\sim W_G(b,D) | Yes | Heterogeneous random weights; `BDgraph::rgwish()` ([Mohammadi and Wit 2025](#ref-R-BDgraph)) |

Table 5: Support-preserving maps from a weighted graph to \Omega\succ 0.
Eigenvalue clipping (reconstruct Q\Lambda\_{\mathrm{clip}}Q^{\top})
fills structural zeros, so it is omitted from the generative truth;
nearest-PD projections are repair tools, not support-preserving
simulators.

#### Why the uniform spectral shift is attractive

If W=Q\Lambda Q^{\top}, then ([Equation 18](#eq-spectral-shift))

W+cI = Q(\Lambda+cI)Q^{\top}. \tag{18}

Every eigenvalue increases by c, while every off-diagonal of W is
unchanged, so the edge set

E=\\(j,k):j\<k,\\ \Omega\_{jk}\neq 0\\ \tag{19}

is identical before and after loading—and **signs of** W are preserved.

The floor \varepsilon also controls the condition number

\kappa(\Omega)=\frac{\lambda\_{\max}(\Omega)}{\lambda\_{\min}(\Omega)}.
\tag{20}

Large loading yields a well-conditioned but weakly dependent network;
loading only slightly above -\lambda\_{\min} strengthens dependence but
can make estimation fragile. Benchmarks should report
partial-correlation distributions and \kappa(\Omega), not only the
adjacency.

### Cell-type covariances: `build_covariance_array_from_precision()`

By default each cell type has its own precision (and thus covariance),

\boldsymbol{\Sigma}\_j = \boldsymbol{\Omega}\_j^{-1}, \qquad
j=1,\ldots,J,

stacked as an array in \mathcal{M}\_{G\times G\times J}. Shared networks
remain available by supplying the same adjacency (or the same graph
model seed path) for every type.

> **Note 5: Why the bulk covariance is \sum_j
> p_j^{2}\boldsymbol{\Sigma}\_j**
>
> If each latent profile is \boldsymbol{x}\_{\cdot
> j}\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\_{\cdot
> j},\boldsymbol{\Sigma}\_j) independently, the bulk
> \boldsymbol{y}=\sum_j p_j\boldsymbol{x}\_{\cdot j} is an **affine**
> image of a stacked Gaussian vector. Affine images of multivariate
> Gaussians remain Gaussian, with mean \boldsymbol{\mu}\boldsymbol{p}
> and covariance \sum_j p_j^{2}\boldsymbol{\Sigma}\_j (the p_j^{2}
> factor is the square of the linear coefficient, not a mixture-weight
> identity). That is the conditional law used in the article (Gaussian
> convolution with covariance \sum_j p_j^{2}\boldsymbol{\Sigma}\_j; see
> the log-likelihood derivation in `article/main.tex`) and by
> [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md).
> Downstream bulk simulation therefore uses
>
> \boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p}) \sim
> \mathcal{N}\_{G}\bigl(\boldsymbol{\mu}\boldsymbol{p}, \textstyle\sum_j
> p_j^{2}\boldsymbol{\Sigma}\_j\bigr).

### Imbalanced cellular proportions

Topology and mean design fix the *reference* moments \boldsymbol{\zeta}.
Mixture difficulty also depends on how uneven the true composition
\boldsymbol{p} is. DeCovarT summarises that imbalance with the
normalised Shannon entropy
[`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md),

H^{\star}(\boldsymbol{p}) = \frac{-\sum\_{j:p\_{j}\>0}p\_{j}\log
p\_{j}}{\log J} \in\[0,1\], \tag{21}

so H^{\star}=0 for a Dirac mass on one type and H^{\star}=1 for the
uniform vector on J types. Zero masses are dropped only inside the sum;
the normaliser still uses the full panel size J.

A practical grid therefore crosses network / mean axes with a few
composition designs, for example:

| Design | Example \boldsymbol{p} (up to permutation) | Typical H^{\star} |
|:---|:---|:---|
| Pure / near-pure | (1,0,\ldots) or (0.95,0.05,0,\ldots) | near 0 |
| One dominant type | (0.7,0.15,0.15) | low–moderate |
| Two-way balance | (0.5,0.5) or (0.4,0.4,0.2) | moderate–high |
| Uniform | (1/J,\ldots,1/J) | 1 |

Table 6: Composition designs scored by normalised Shannon entropy.

The Bioconductor package `SimBu` ([Dietrich 2024](#ref-R-SimBu)) offers
a complementary pseudo-bulk toolkit with six named fraction scenarios
(`even`, `random`, `mirror_db`, `weighted`, `pure`, `custom`) when
aggregating single-cell profiles; see the [SimBu “Simulate pseudo bulk
datasets”](https://bioconductor.org/packages/release/bioc/vignettes/SimBu/inst/doc/SimBu.html)
section. [Figure 9](#fig-simbu-entropy) shows five of those designs for
J=6 types as stacked compositions, labelled by
H^{\star}(\boldsymbol{p}). High-entropy (`even`) mixtures spread mass
across types; low-entropy (`pure`) mixtures collapse onto one type.
DeCovarT’s Gaussian convolution uses the same idea at the *moment*
layer: supply `p` (or a J\times N matrix of sample-wise ratios) to
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
and report H^{\star} alongside overlap or condition-number diagnostics.
The hybrid J=3 manuscript scenario uses the imbalanced design
\boldsymbol{p}=(0.4,0.4,0.2) ([Manuscript synthetic simulation
scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-manuscript-scenarios.html#sec-scenario-grid)).

![](synthetic-scenarios_files/figure-html/fig-simbu-entropy-1.png)

Figure 9: Illustrative J=6 compositions spanning `SimBu`-style fraction
scenarios, ordered by decreasing normalised Shannon entropy
H^{\star}(\boldsymbol{p}). Each bar is one design; `geom_label` reports
H^{\star}. `mirror_db` is omitted (it copies an empirical atlas rather
than a fixed simplex vector).

## Conclusions

A reproducible synthetic study should keep three design layers separate,
then score both the reference geometry and the mixture composition:

1.  **Mean signatures.** Build \boldsymbol{\mu} with
    [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md):
    hold `mean_scale` s fixed and dial only `target_cosine` \rho. Report
    the realised pairwise cosine ([Theorem 1](#thm-cosine-asymptotics)),
    Euclidean separation, and optionally \kappa_2(\boldsymbol{\mu}).
2.  **Precision graphs.** Draw an undirected skeleton (`scale_free`,
    `stochastic_block_model`, or `small_world`), assign i.i.d. signed
    weights (`prop_inhibitory`), and complete \Omega\succ 0 by the
    uniform spectral shift. Report \kappa(\Omega) and the
    partial-correlation sign mix.
3.  **Cellular compositions.** Choose \boldsymbol{p} on a Shannon grid
    from uniform (H^{\star}=1) to near-pure (H^{\star}\approx 0),
    matching the `SimBu` fraction vocabulary when comparing to
    pseudo-bulk tools ([Figure 9](#fig-simbu-entropy)).
4.  **Bulk draws.** Simulate \boldsymbol{Y} from the Gaussian
    convolution
    \boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p})\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p},\sum_j
    p_j^{2}\boldsymbol{\Sigma}\_j) and pass the same moments to every
    solver.
5.  **Hybrid stress test.** The manuscript scenario (two mean-collinear
    types told apart only by topology, plus an orthogonal third type) is
    documented in [Manuscript synthetic simulation
    scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-manuscript-scenarios.html#sec-scenario-grid).

See [Note 6](#nte-benchmark-spec) for a compact factorial checklist.

> **Important 6: Recommended benchmark specification**
>
> Use one pipeline for every topology
> ([Equation 22](#eq-benchmark-pipe), [Figure 6](#fig-ggm-pipeline)):
>
> \text{topology} \rightarrow \text{signed weights} \rightarrow
> \text{SPD precision} \rightarrow \text{mean design} \rightarrow
> \text{latent Gaussian} \rightarrow \text{observation model} \tag{22}
>
> | Axis | Suggested levels |
> |----|----|
> | Dimension (genes) | G \in \\100,500,1000\\ |
> | Topology | ER, band, hub, scale-free, SBM, small-world |
> | Expected degree | \approx 2,4,8 (keep graphs sparse) |
> | Edge weights | Constant; i.i.d. signed; partial-correlation scaling; G-Wishart |
> | Conditioning | Report \kappa(\Omega) and partial-correlation summaries |
> | Composition | Pure / weighted / balanced / uniform ([Section 4.5](#sec-imbalanced)); record H^{\star} |
>
> Table 7: Compact factorial axes for a reproducible undirected GGM
> simulation study.
>
> Prefer the uniform spectral shift ([Equation 18](#eq-spectral-shift))
> or partial-correlation scaling ([Equation 17](#eq-partial-scale)) when
> support and signs must be exact; add the mean layer **independently**
> of graph generation. The end-to-end **hybrid multi-topology reference
> scenario** (two mean-collinear cell types distinguished only by
> network topology, a third orthogonal type, NSGA-II panel curation, and
> the frequentist solver comparison on imbalanced mixtures) is
> documented in [Manuscript synthetic simulation
> scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-manuscript-scenarios.html#sec-scenario-grid).
> Edge-case numerical checks (gene-wise z-score equivariance, collinear
> signatures, small bulk perturbations, random ALR starts) live in the
> [appendix simulation
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT_appendix_simualtion_frameworks.md).

## Perspectives on high-dimensional or more biologically realistic designs

Beyond the families in [Table 1](#tbl-topologies):

- **Degree-prescribed / configuration** graphs separate degree
  heterogeneity from preferential-attachment *mechanisms*.
- **Spatial / geometric** graphs suit spatial transcriptomics or
  chromatin neighbourhoods.
- **Differential networks:**
  - Start from a shared base G_0
  - Control the differential shift \Omega_k=\Omega_0+\Delta_k
  - Ensure the SPD step
  - Alternatively, vary the *support, weight, and sign* separately
    ([Federico et al. 2023](#ref-federicoStructureLearningGene2023); [Wu
    and Luo 2022](#ref-wuEstimatingHeterogeneousGene2022)).
- **Other distributions,** using the GGM layer as a parameter:
  - **Latent-variable GGMs** — sparse \Omega_S plus low-rank confounding
    Bf_i.
  - **Nonparanormal / Gaussian-copula** margins on a latent Gaussian
    graph.
  - **Poisson–log-normal / compositional** observation layers on latent
    log-abundances ([Chiquet et al.
    2018](#ref-chiquetVariationalInferenceSparse2018)).
  - **Zero-inflated mixture GGMs** with cell-type-specific
    (\mu_k,\Omega_k) and expression-dependent dropout ([Wu and Luo
    2022](#ref-wuEstimatingHeterogeneousGene2022)).

## References

Aliee, Hananeh, and Fabian J. Theis. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems* 12. <https://doi.org/10.1016/j.cels.2021.05.006>.

Ba, Kalidou, Rodolphe Thiébaut, Xavier Hinaut, and Boris Hejblum. 2026.
*When Less Is Not More: DICEPro Mitigates the Impact of Incomplete
Reference Matrices on Cellular Frequency Deconvolution*. bioRxiv.
<https://doi.org/10.64898/2026.06.17.732876>.

Barabási, Albert-László, and Réka Albert. 1999. ‘Emergence of Scaling in
Random Networks’. *Science* 286.
<https://doi.org/10.1126/science.286.5439.509>.

Broido, Anna D., and Aaron Clauset. 2019. ‘Scale-Free Networks Are
Rare’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-019-08746-5>.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2018.
*Variational Inference for Sparse Network Reconstruction from Count
Data*. arXiv. <https://doi.org/10.48550/arxiv.1806.03120>.

Dietrich, Alexander. 2024. *SimBu: Bias-Aware Simulation of Bulk RNA-Seq
Data with Variable Cell-Type Composition*.
<https://doi.org/10.18129/B9.bioc.SimBu>.

Federico, Anthony, Joseph Kern, Xaralabos Varelas, and Stefano Monti.
2023. ‘Structure Learning for Gene Regulatory Networks’. *PLOS
Computational Biology* 19.
<https://doi.org/10.1371/journal.pcbi.1011118>.

Jiang, Haoming, Xinyu Fei, Han Liu, et al. 2026. *Huge: High-Dimensional
Undirected Graph Estimation*. <https://github.com/Gatech-Flash/huge>.

Lima-Mendez, Gipsi, and Jacques van Helden. 2009. ‘The Powerful Law of
the Power Law and Other Myths in Network Biology1’. *Molecular
BioSystems(MBS)* 5. <https://doi.org/10.1039/b908681a>.

Mohammadi, Reza, and Ernst Wit. 2019. ‘BDgraph: An R Package for
Bayesian Structure Learning in Graphical Models’. *Journal of
Statistical Software* 89 (3): 1–30.
<https://doi.org/10.18637/jss.v089.i03>.

Mohammadi, Reza, and Ernst Wit. 2025. *BDgraph: Bayesian Structure
Learning in Graphical Models Using Birth-Death MCMC*.
<https://www.uva.nl/profile/a.mohammadi>.

Schelker, Max, Sonia Feau, Jinyan Du, et al. 2017. ‘Estimation of Immune
Cell Content in Tumour Tissue Using Single-Cell RNA-seq Data’. *Nature
Communications* 8 (1): 2032.
<https://doi.org/10.1038/s41467-017-02289-3>.

Takahashi, Daniel Yasumasa, João Ricardo Sato, Carlos Eduardo Ferreira,
and André Fujita. 2012. ‘Discriminating Different Classes of Biological
Networks by Analyzing the Graphs Spectra Distribution’. *PLOS ONE* 7.
<https://doi.org/10.1371/journal.pone.0049949>.

Wu, Qiuyu, and Xiangyu Luo. 2022. ‘Estimating Heterogeneous Gene
Regulatory Networks from Zero-Inflated Single-Cell Expression Data’.
*The Annals of Applied Statistics* 16.
<https://doi.org/10.1214/21-aoas1582>.

Zhang, Rong, Zhao Ren, and Wei Chen. 2018. ‘SILGGM: An Extensive R
Package for Efficient Statistical Inference in Large-Scale Gene
Networks’. *PLOS Computational Biology* 14.
<https://doi.org/10.1371/journal.pcbi.1006369>.
