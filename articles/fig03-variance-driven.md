# §2.2 Variance-driven hybrid scenario (G = 50, J = 3)

``` r

library(DeCovarT)
```

The bivariate toy confirms the theoretical benefit of covariance
modelling in a transparent setting but cannot demonstrate realistic
gene-regulatory network (GRN) topology. This scenario uses G = 50 genes
and J = 3 cell types in which types 1 and 2 are **near mean-collinear**
(\rho\_{12} \approx 0.97): the network topology encoded in
\boldsymbol{\Sigma}\_j is the discriminative signal. Construction, ADEMP
benchmark, and the static network figure live in a **single** script,
`scripts/fig03_variance_driven.R` (seed `20260807`).

> **Script:** `Rscript scripts/fig03_variance_driven.R` (full, n = 200)
> or `N_REPLICATES=2 Rscript scripts/fig03_variance_driven.R` (smoke).
> Moments: `data/synthetic_networks/true_grn_moments.rds`. Figures:
> `output/fig03/`. ADEMP reporting: [how to build synthetic
> scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.html#sec-ademp).

## Generative model

Types 1 and 2 share a high local cosine on 30 genes and differ only
through \boldsymbol{\Omega}; type 3 has a compact marker block; 10 genes
are a null (`equal_all`) block.

\underbrace{\mathcal{G}\_{12}}\_{30\text{ genes}} \\ \cup\\
\underbrace{\mathcal{G}\_{3}}\_{10\text{ genes}} \\ \cup\\
\underbrace{\mathcal{G}\_{\mathrm{eq}}}\_{10\text{ genes}}, \qquad
\boldsymbol{\Omega}\_j = \mathrm{blockdiag}\bigl(
\boldsymbol{\Omega}\_j^{(\mathcal{G}\_{12})},\\
\boldsymbol{\Omega}\_j^{(\mathcal{G}\_{3})},\\
\boldsymbol{\Omega}\_j^{(\mathcal{G}\_{\mathrm{eq}})} \bigr).

\boldsymbol{\mu} uses two calls to
[`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
([mean signature
construction](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.html#sec-lr-block-mu)):
\rho\_{12}=0.95 on \mathcal{G}\_{12}, then \rho\_{3}=0.05 on a marker
pool for type 3. Local topologies follow the BLGGM two-scale pattern
([Table 1](#tbl-hybrid-topology-design); [GGM
networks](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.html#sec-ggm-networks)).

| Gene block | \# genes | Cell type 1 | Cell type 2 | Cell type 3 |
|----|----|----|----|----|
| `shared_12_vs_3` | 30 | SBM (hub-like modules) | star (one driver) | Erdős–Rényi |
| `marker_3` | 10 | Erdős–Rényi | Erdős–Rényi | scale-free |
| `equal_all` | 10 | Erdős–Rényi | Erdős–Rényi | Erdős–Rényi |

Table 1: Local topology per (gene block, cell type); block-diagonal
precision support.

| pair                     | cosine | Euclidean |
|:-------------------------|-------:|----------:|
| celltype_1 vs celltype_2 |  0.973 |     2.765 |
| celltype_1 vs celltype_3 |  0.562 |    11.219 |
| celltype_2 vs celltype_3 |  0.562 |    11.219 |

Table 2: Realised pairwise geometry of the 50\times 3 mean signature
(seed 20260807).

| cell type | topology | \lambda\_{\min} | \lambda\_{\max} | \kappa(\Omega) | prop. inhib. |
|:---|:---|---:|---:|---:|---:|
| celltype_1 | SBM | 0.100 | 1.685 | 16.846 | 0.500 |
| celltype_2 | star | 0.100 | 3.331 | 33.311 | 0.489 |
| celltype_3 | scale-free | 0.100 | 2.096 | 20.961 | 0.500 |

Table 3: Precision spectra of each cell type’s \boldsymbol{\Omega}\_j
(seed 20260807).

Types 1 and 2 remain near-collinear on the full G=50 vector (cosine
\approx 0.97). Type 3 is only *moderately* separated (cosine \approx
0.56) because the shared \mathcal{G}\_{12} baseline and the flat
`equal_all` block inflate every pairwise inner product—the finite-J
cosine bias of [asymptotic cosine
control](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.html#sec-asymptotic-cosine).
[Table 3](#tbl-topology-diagnostics) shows that types 1 and 2 still
differ in \boldsymbol{\Omega}.

![](figures/fig_network_topologies.png)

Figure 1: Cell-type-specific block topologies: SBM/hub-like modules
(type 1) and a star (type 2) on the 30 `shared_12_vs_3` genes;
scale-free wiring (type 3) on `marker_3`. The `equal_all` block is
Erdős–Rényi in every type. Edge colour: precision sign (red inhibitory,
teal activatory).

### Related literature (`muscat` / `scDD`)

Pseudobulk tools such as `muscat` ([Crowell et al.
2020](#ref-crowellMuscatDetectsSubpopulationspecific2020)) aggregate
single-cell counts across samples. `scDD` ([Korthauer et al.
2016](#ref-korthauerStatisticalApproachIdentifying2016)) partitions
genes into EE, DE, DP, DM, and DB patterns. The hybrid scenario is a
network-level analogue of three of these ([Fig. 2](#fig-dd-taxonomy)):
`equal_all` is EE; `marker_3` is DE for type 3; `shared_12_vs_3` has
**no mean shift between types 1 and 2**, only a precision-topology
contrast (DM-like for deconvolution, not scDD’s original single-gene
modality test).

![](figures/fig_dd_taxonomy_ee_de_dm.jpg)

Figure 2: Schematic EE / DE / DM-like taxonomy for the three gene
blocks.

Two generators were considered but not adopted:

- **Toeplitz covariances**: encode only local (banded) dependencies;
  cannot model cell-type-specific sparse precision.
- **`MixSim` overlap** ([Melnykov et al.
  2012](#ref-melnykovMixSimPackageSimulating2012)): designed for
  isotropic or full-rank Gaussian clusters; no cell-type-specific
  precision input.

## Inference

| Factor | Levels |
|----|----|
| Composition \boldsymbol{p} | balanced (1/3,1/3,1/3); moderately unbalanced (0.5,0.3,0.2); highly unbalanced (0.7,0.2,0.1) |
| Algorithms | NNLS, DeconRNASeq, Marquardt–Levenberg |
| **Total scenarios** | 3 \times 3 = \mathbf{9} (per fixed GRN realisation) |

N = 200 Monte Carlo replicates per scenario (smoke test: n = 2). All
solvers start from the open simplex via
[`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md)
(barycentre by default; Dirichlet or QP via `initial_p`).

| Algorithm | Function | Constraint |
|:---|:---|:---|
| NNLS | `deconvolute_ratios_nnls` | [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md) |
| QP (`DeconRNASeq`-style) | `deconvolute_ratios_deconrnaseq` | simplex equality / inequality |
| Marquardt–Levenberg | `deconvolute_ratios_Marquardt_Levenberg` | ILR \to open simplex |

Table 4: Solvers used in the variance-driven comparison.

The shipped catalogue (CIBERSORT-style \nu-SVR,
[`MASS::rlm`](https://rdrr.io/pkg/MASS/man/rlm.html), GLS,
Newton–Raphson, L-BFGS-B, BFGS, simulated annealing) is listed in
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).
A temporary solver smoke test on this GRN lives in
`tests/manual/deconvolute_simulated_scenario.R` (Rbuildignored).

``` r

N_REPLICATES <- as.integer(Sys.getenv("N_REPLICATES", unset = "200"))
ALGORITHMS <- c("NNLS", "DeconRNASeq", "Marquardt-Levenberg")
# Full pipeline: Rscript scripts/fig03_variance_driven.R
```

## Visualisations

| Output | Description |
|----|----|
| `output/fig03/fig03_raincloud.pdf` | Raincloud of RMSE by algorithm and composition |
| `output/fig03/fig03_forest.pdf` | ADEMP forest plot |
| `output/fig03/fig03_metric_dots.pdf` | Metric dot plot across compositions |

#### Expected findings

Because CT 1 and CT 2 share near-identical mean profiles (\rho\_{12}
\approx 0.97), NNLS and DeconRNASeq cannot separate them — their RMSE
for that pair should sit near the variance floor. DeCovarT
(Marquardt–Levenberg) uses the block-diagonal precision structure; its
advantage grows as the composition becomes unbalanced and the rare type
(CT 3 at 10%) requires second-order information.

### See also

- Bivariate toy:
  [§2.1](https://bastienchassagnol.github.io/DeCovarT/articles/fig02-bivariate-toy.md)
- Moment generator: [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.md)
- Feature-selection on this design: [Appendix
  S6](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S6-feature-selection.md)
- Identifiability of collinear means: [Appendix
  S1](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.md)

### References

Crowell, Helena L., Charlotte Soneson, Pierre-Luc Germain, et al. 2020.
‘Muscat Detects Subpopulation-Specific State Transitions from
Multi-Sample Multi-Condition Single-Cell Transcriptomics Data’. *Nature
Communications* 11: 6077. <https://doi.org/10.1038/s41467-020-19894-4>.

Korthauer, Keegan D., Li-Fang Chu, Michael A. Newton, et al. 2016. ‘A
Statistical Approach for Identifying Differential Distributions in
Single-Cell RNA-seq Experiments’. *Genome Biology* 17: 222.
<https://doi.org/10.1186/s13059-016-1077-y>.

Melnykov, Volodymyr, Wei-Chen Chen, and Ranjan Maitra. 2012. ‘MixSim: An
R Package for Simulating Data to Study Performance of Clustering
Algorithms’. *Journal of Statistical Software* 51.
<https://doi.org/10.18637/jss.v051.i12>.
