# §2.2 Variance-driven hybrid scenario (G = 50, J = 3)

``` r

library(DeCovarT)
```

The bivariate toy confirms the theoretical benefit of covariance
modelling in a transparent setting but cannot demonstrate realistic
gene-regulatory network (GRN) topology. This scenario uses G = 50 genes
and J = 3 cell types in which types 1 and 2 are **mean-collinear**
(\cos(\boldsymbol{\mu}\_{\cdot 1},\boldsymbol{\mu}\_{\cdot 2})=0.9): the
network topology encoded in \boldsymbol{\Sigma}\_j is the remaining
discriminative signal. Construction, ADEMP benchmark, and the static
network figure live in a **single** script,
`scripts/fig03_variance_driven.R` (seed `20260807`).

> **Script:** `Rscript scripts/fig03_variance_driven.R` (full, n = 50)
> or `N_REPLICATES=2 Rscript scripts/fig03_variance_driven.R` (smoke).
> Figures: `output/fig03/`. ADEMP reporting: [how to build synthetic
> scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ademp).
> Scenario geometry: [scenario
> descriptors](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-scenario-descriptors).

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

Table 1: Exact pairwise geometry of the 50\times 3 mean signature
(`target_gram`, `mean_scale` =10).

Each cell type independently draws an undirected skeleton from one of
three generators (`stochastic_block_model`, `erdos_renyi`, `hub` /
star), then completes a signed precision with
[`build_normalised_precision()`](https://bastienchassagnol.github.io/DeCovarT/reference/build_normalised_precision.md)
([GGM
networks](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ggm-networks)).
The spectral cushion u (`precision_shift`) sets
\lambda\_{\min}(\boldsymbol{\Omega}\_j)=u after the shift, so smaller u
worsens \kappa(\boldsymbol{\Omega}\_j)=\lambda\_{\max}/u and, after
inversion and convolution,
\kappa\\\boldsymbol{\Sigma}(\boldsymbol{p})\\.

| Factor                      | Levels                                 |
|-----------------------------|----------------------------------------|
| Graph model (per cell type) | SBM, Erdős–Rényi, hub/star             |
| Precision shift u           | well (0.5); moderate (0.1); ill (0.02) |
| Topology assignments        | 3^3=\mathbf{27}                        |
| Covariance draws            | 27\times 3=\mathbf{81}                 |

Table 2: Covariance factorial (graphs \times precision cushion).

[`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md)
records the realised \kappa\\\boldsymbol{\Sigma}(\boldsymbol{p})\\ and
the covariance-information fraction f\_{\mathrm{cov}} on every row (see
the [descriptors
section](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-scenario-descriptors),
especially the tangent-Fisher callout). Those columns are **measured**,
not additional design factors: u is the knob; f\_{\mathrm{cov}} reports
how much of the Fisher information is second-order.

![](figures/fig_network_topologies.png)

Figure 1: Three generator examples on G=50 nodes (SBM, Erdős–Rényi,
hub/star). Each cell type in the factorial is assigned one of these
families independently.

Two generators were considered but not adopted:

- **Toeplitz covariances**: encode only local (banded) dependencies;
  cannot model cell-type-specific sparse precision.
- **`MixSim` overlap** ([Melnykov et al.
  2012](#ref-melnykovMixSimPackageSimulating2012)): designed for
  isotropic or full-rank Gaussian clusters; no cell-type-specific
  precision input.

## Inference

Compositions are one-dominant vectors from
[`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md)
targeting Pielou evenness H^{\star}\in\\1,0.5,0.1\\ (balanced;
moderately unbalanced; highly unbalanced). All DeCovarT solvers start at
the barycentre (`initial_p = "barycentre"`). Other initialisation
strategies are reserved for separate scenarios.

| Factor | Levels |
|----|----|
| Composition \boldsymbol{p} | H^{\star}=1; H^{\star}=0.5; H^{\star}=0.1 |
| Algorithms | DeconRNASeq (`lsei`), CIBERSORT, L-BFGS-B, Newton–Raphson, Marquardt–Levenberg |
| **Total scenarios** | 27\times 3\times 3=\mathbf{243} |

N = 50 Monte Carlo replicates per scenario (smoke test:
`N_REPLICATES=2`).

| Algorithm | Function | Constraint |
|:---|:---|:---|
| QP (`DeconRNASeq`-style) | `deconvolute_ratios_deconrnaseq` | simplex equality / inequality |
| CIBERSORT-style \nu-SVR | `deconvolute_ratios_cibersort` | non-negative then [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md) |
| L-BFGS-B | `deconvolute_ratios_L_BFGS_B` | ILR \to open simplex |
| Newton–Raphson | `deconvolute_ratios_Newton_Raphson` | ILR / `nlminb` |
| Marquardt–Levenberg | `deconvolute_ratios_Marquardt_Levenberg` | ILR \to open simplex |

Table 3: Solvers used in the variance-driven comparison.

Mean-only baselines (LSEI, CIBERSORT) ignore \boldsymbol{\Sigma}\_j. The
three DeCovarT maps use the same convolution likelihood. The shipped
catalogue (robust linear model, GLS, BFGS, simulated annealing) is
listed in
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

``` r

N_REPLICATES <- as.integer(Sys.getenv("N_REPLICATES", unset = "50"))
ALGORITHMS <- c(
  "lsei",
  "cibersort",
  "LBFGS",
  "Newton-Raphson",
  "Marquardt-Levenberg"
)
# Full pipeline: Rscript scripts/fig03_variance_driven.R
```

## Visualisations

| Output | Description |
|----|----|
| `output/fig03/fig03_raincloud.pdf` | Raincloud of Monte Carlo errors, faceted by composition and \kappa label |
| `output/fig03/fig03_forest.pdf` | ADEMP forest plot |
| `output/fig03/fig03_metric_dots.pdf` | Metric dot plot (RMSE, MAE, coverage) |

#### Expected findings

Because CT 1 and CT 2 share cosine 0.9, LSEI and CIBERSORT cannot
separate them from the mean signature alone. DeCovarT (L-BFGS-B,
Newton–Raphson, Marquardt–Levenberg) uses \boldsymbol{\Sigma}\_j; its
advantage should grow where f\_{\mathrm{cov}} is large, where u is small
(ill-conditioned but still informative second-order structure), and
where the composition is highly unbalanced (H^{\star}=0.1).

### See also

- Bivariate toy:
  [§2.1](https://bastienchassagnol.github.io/DeCovarT/articles/fig02-bivariate-toy.md)
- Moment generator: [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.md)
- Feature-selection on this design: [Appendix
  S6](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S6-feature-selection.md)
- Identifiability of collinear means: [Appendix
  S1](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.md)

### References

Melnykov, Volodymyr, Wei-Chen Chen, and Ranjan Maitra. 2012. ‘MixSim: An
R Package for Simulating Data to Study Performance of Clustering
Algorithms’. *Journal of Statistical Software* 51.
<https://doi.org/10.18637/jss.v051.i12>.
