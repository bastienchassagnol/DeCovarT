# DeCovarT perspectives

> **Scope**
>
> Outlook on extending DeCovarT beyond bulk transcriptomic
> deconvolution. Sections below cover Scheffé-type mixture extensions,
> lineage-aware hierarchical deconvolution, and compositional
> validation; future topics (spatial transcriptomics, multi-modal /
> multi-omics) will be added here. Compositional reparametrisation is
> implemented in
> [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
> and documented numerically in the [softmax / ALR
> derivatives](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md)
> vignette.

Bulk transcriptomic deconvolution is one instance of a **mixture inverse
problem**: a composite signal (tissue expression, a spectrum, a
chromatogram) is observed, and we seek the proportions \boldsymbol{p} of
underlying components. Because \boldsymbol{p} is **compositional**
(parts of a whole), estimates must lie on the simplex
\Delta^{J-1}=\\\boldsymbol{p}\in\[0,1\]^J:\sum_j p_j=1\\. ([Chassagnol
et al. 2023](#ref-chassagnolDeCovarTMultidimensionalProbalistic2023))
model each reference profile as multivariate Gaussian and recover
\boldsymbol{p} by constrained maximum likelihood; classical
deconvolution tools such as `CIBERSORT` assume gene-wise independence
instead ([Newman et al. 2015](#ref-newmanRobustEnumerationCell2015)).

## Extending the generative model

### Mixture design and first-order structure

In designed mixture experiments (chemistry, food science, formulation),
component proportions are explicit regressors with no free intercept
because \sum_j p_j=1. A Scheffé cubic model for three components reads

y = \beta_1 p_1 + \beta_2 p_2 + \beta_3 p_3 + \beta\_{12} p_1 p_2 +
\beta\_{13} p_1 p_3 + \beta\_{23} p_2 p_3 + \beta\_{123} p_1 p_2 p_3,
\tag{1}

with every p_j\in\[0,1\] and \sum_j p_j=1 ([Brown et al.
2015](#ref-brownGeneralBlendingModels2015); [Aitchison
1982](#ref-aitchisonStatisticalAnalysisCompositional1982)). DeCovarT’s
current release implements only the **first-order** part of
[Eq. 1](#eq-scheffe-cubic) in a probabilistic setting: for cell type j,
reference mean \boldsymbol{\mu}\_j and covariance \boldsymbol{\Sigma}\_j
are given, and the bulk mixture satisfies

\boldsymbol{Y}\mid\boldsymbol{p} \sim \mathcal{N}\\\left(
\sum\_{j=1}^{J} p_j\\\boldsymbol{\mu}\_j,\\ \sum\_{j=1}^{J}
p_j^{2}\\\boldsymbol{\Sigma}\_j \right). \tag{2}

[Eq. 2](#eq-gaussian-convolution) is the Gaussian-convolution likelihood
optimised in
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md);
it extends mean-only linear deconvolution by propagating gene–gene
covariance, not by adding Scheffé interaction monomials.

### Interaction and synergy terms

Synergy or antagonism between components corresponds to **second-order**
terms in [Eq. 1](#eq-scheffe-cubic): a positive \beta\_{12} means the
joint presence of components 1 and 2 raises the response above the
additive prediction \beta_1 p_1+\beta_2 p_2. A natural generative
extension would enrich the mean and, optionally, the covariance:

\mathbb{E}\[\boldsymbol{Y}\] = \sum_j p_j\\\boldsymbol{\mu}\_j +
\sum\_{j\<k} f(p_j,p_k)\\\boldsymbol{\mu}\_{jk} + \cdots, \tag{3}

\boldsymbol{\Sigma}(\boldsymbol{p}) = \sum_j
p_j^{2}\\\boldsymbol{\Sigma}\_j + \sum\_{j\<k}
g(p_j,p_k)\\\boldsymbol{\Sigma}\_{jk} + \cdots, \tag{4}

with f,g often taken as p_j p_k. Allowing \operatorname{Cov}(\boldsymbol
{X}\_j,\boldsymbol{X}\_k)\neq\boldsymbol{0} relaxes DeCovarT’s
independence assumption between components and captures non-additive
mixing. Such a model remains compatible with the same simplex constraint
on \boldsymbol{p} and could still be fitted by constrained MLE, at the
cost of many additional interaction parameters and careful guard against
overfitting.

### Lineage-aware hierarchical deconvolution

Flat deconvolution treats every cell type as a sibling on one simplex.
**Hierarchical** methods instead estimate proportions at several
granularities linked by a cell-type tree: broad compartments first, then
closely related subtypes whose mean signatures are collinear at the leaf
level.

We reuse the DeCovarT set-up from the
[article](https://bastienchassagnol.github.io/DeCovarT/article/main.pdf):
gene g=1,\ldots,G; cell type or tree node j; sample i=1,\ldots,N
(independence across i). Slices follow the usual dot notation
(\boldsymbol{\mu}\_{\cdot j}, \boldsymbol{y}\_{\cdot i}).

- \boldsymbol{y}=(y\_{gi})\in\mathbb{R}\_{+}^{G\times N}: bulk
  expression; \boldsymbol{y}\_{\cdot i}\in\mathbb{R}\_{+}^{G} is sample
  i.
- \boldsymbol{\mu}=(\mu\_{gj})\in\mathbb{R}^{G\times J}: mean signature;
  \boldsymbol{\mu}\_{\cdot j}\in\mathbb{R}^{G} is the cross-sample mean
  profile of type j (as in `MuSiC` ([Wang et al.
  2019](#ref-wangBulkTissueCell2019))).
- \boldsymbol{p}=(p\_{ji})\in\\\]0,1\[^{J\times N}: unknown relative
  proportions; \boldsymbol{p}\_{\cdot i}\in\\\]0,1\[^{J} is the
  composition of sample i.

\boldsymbol{y}\_{\cdot i}=\boldsymbol{\mu}\\\boldsymbol{p}\_{\cdot i},
\tag{5}

\sum\_{j=1}^{J} p\_{ji}=1, \qquad p\_{ji}\ge 0 \qquad (i=1,\ldots,N).
\tag{6}

p\_{j,i}=\sum\_{c\in\mathcal{C}(j)} p\_{c,i} \qquad (i=1,\ldots,N).
\tag{7}

At each tree level the bulk mixture satisfies
[Eq. 5](#eq-hierarchical-linear) with each composition on the unit
simplex ([Eq. 6](#eq-simplex-constraint)). For a parent node j with
child set \mathcal{C}(j), coherent hierarchical estimates require
[Eq. 7](#eq-parent-child).

#### Hard versus soft lineage constraints

Despite diverse implementations, hierarchical bulk mRNA deconvolution
methods fall into **two constraint regimes**:

- **Hard constraints** — child lineage ratios are forced to sum exactly
  to the parent estimate ([Eq. 7](#eq-parent-child)). `MuSiC` and `SCDC`
  (tree-guided mode) first fix a parent total \hat{p}\_{j,i}, then fit
  children q\_{c,i} on the within-parent simplex so that
  \hat{p}\_{c,i}=q\_{c,i}\\\hat{p}\_{j,i} ([Wang et al.
  2019](#ref-wangBulkTissueCell2019); [Dong et al.
  2021](#ref-dongSCDCBulkGene2021)). `HIDE` applies the same equality by
  renormalising child estimates after each top-down split ([Völkl et al.
  2025](#ref-voelklHIDEHierarchicalCelltype2025)); `Rectangle`
  propagates coarse cluster fractions as hard bounds on fine-grained
  estimates ([Eder et al. 2026](#ref-ederRectangleRobustScalable2026)).
- **Soft constraints** — parent–child coherence is encouraged but not
  enforced exactly. `HiDecon` adds a squared mismatch penalty to the
  deconvolution loss ([Huang et al.
  2024](#ref-huangAccurateEstimationRare2024)); `Kassandra` trains
  separate predictors at predefined hierarchy levels without an explicit
  parent–child equality ([Zaitsev et al.
  2022](#ref-zaitsevPreciseReconstructionTME2022)).

Hard splitting after fixing the parent ([Eq. 8](#eq-hard-child-split)):

\hat{p}\_{c,i}=q\_{c,i}\\\hat{p}\_{j,i}, \quad q\_{c,i}\ge 0,\\
\sum\_{c\in\mathcal{C}(j)} q\_{c,i}=1. \tag{8}

`HIDE` renormalisation ([Eq. 9](#eq-hide-renorm)):

\hat{p}\_{c,i}\leftarrow \hat{p}\_{c,i}\\
\frac{\hat{p}\_{j,i}}{\sum\_{c'\in\mathcal{C}(j)}\hat{p}\_{c',i}}.
\tag{9}

`HiDecon` soft penalty ([Eq. 10](#eq-soft-penalty)):

\mathcal{L}\_{\mathrm{tot}} = \mathcal{L}\_{\mathrm{deconv}} +
\lambda\sum\_{j}\sum\_{i=1}^{N} \Bigl(
p\_{j,i}-\sum\_{c\in\mathcal{C}(j)} p\_{c,i} \Bigr)^{2}. \tag{10}

Within each regime, methods further differ by whether they reuse a
**single reference signature** at every granularity or **reselect /
reweight genes** locally at each internal node
([Figure 1](#fig-hierarchical-taxonomy)).

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  A["Lineage-aware bulk mRNA deconvolution"] --> B{"Parent–child linkage"}
  B -->|"Hard equality (Eq. parent–child)"| C["MuSiC, SCDC, HIDE, Rectangle"]
  B -->|"Soft / uncoupled levels"| D["HiDecon, Kassandra"]
  C --> E{"Genes at each node"}
  D --> F{"Genes at each node"}
  E -->|"Fixed reference matrix"| G["SCDC ensemble"]
  E -->|"Node-local reweighting"| H["MuSiC, HIDE, Rectangle"]
  F -->|"Shared matrix + penalty"| I["HiDecon"]
  F -->|"Level-specific signatures"| J["Kassandra"]
```

Figure 1: Taxonomy of lineage-aware bulk mRNA deconvolution:
parent–child constraint (hard vs soft) and reference-gene strategy
(fixed matrix vs node-local weights).

`xCell` 2.0 ([Angel et al. 2025](#ref-angelXCell20Robust2025)) uses Cell
Ontology ancestry to build signatures and control spillover between
related labels, but its outputs are enrichment scores rather than a
compositional tree; it sits outside the hard/soft fraction taxonomy
above.

| method | constraint | lineage constraint | gene strategy |
|----|----|----|----|
| MuSiC | Hard | Fix \\\hat{p}\_{j,i}\\; \\\hat{p}\_{c,i}=q\_{c,i}\hat{p}\_{j,i}\\, \\\sum_c q\_{c,i}=1\\ | Node-local consistent genes |
| SCDC | Hard (tree-guided) | MuSiC step; \\\hat{p}\_{c,i}=q\_{c,i}\hat{p}\_{j,i}\\ | Full reference; optional node split |
| HIDE | Hard | Fit then renormalise: \\\hat{p}\_{c,i}\leftarrow\hat{p}\_{c,i}\hat{p}\_{j,i}/\sum\_{c'}\hat{p}\_{c',i}\\ | Parent-specific learned weights |
| Rectangle | Hard (multiscale) | Coarse \\\boldsymbol{p}^{(\mathrm{coarse})}\\ bounds fine estimates | Clustered + direct signatures |
| HiDecon | Soft | \\\mathcal{L}\_{\mathrm{tot}}=\mathcal{L}\_{\mathrm{deconv}}+\lambda\sum\_{j,i}(p\_{j,i}-\sum_c p\_{c,i})^2\\ | Shared \$\boldsymbol{\mu}\$; joint penalised fit |
| Kassandra | Soft / uncoupled | Separate models; no parent–child equality | Level-specific trained signatures |
| xCell 2.0 | None (scores) | Ontology spillover control only | Ontology-aware signature generation |

Table 1: Lineage-aware bulk mRNA deconvolution: constraint regime,
enforcing equation, and gene strategy (2019–2026 census).

[Table 1](#tbl-hierarchical-methods) and
[Figure 1](#fig-hierarchical-taxonomy) summarise the landscape aligned
with [Eq. 5](#eq-hierarchical-linear)–[Eq. 7](#eq-parent-child).
**Hard** methods (`MuSiC`, `HIDE`, `Rectangle`) are easier to audit node
by node but inherit upstream bias; **soft** coupling (`HiDecon`) can
stabilise rare subtypes when the tree is credible. **Node-local gene
reweighting** (`MuSiC`, `HIDE`, `Rectangle`) targets collinear siblings;
**fixed or level-specific signatures** (`SCDC`, `Kassandra`) trade
flexibility for simpler reference construction.

A natural DeCovarT extension would combine
[Eq. 2](#eq-gaussian-convolution) with either
[Eq. 8](#eq-hard-child-split) or [Eq. 10](#eq-soft-penalty), node-local
marker panels, and ALR optimisation ([Sec. 2.1](#sec-alr-reparam)) when
a curated lineage tree is supplied—benchmarked at broad, intermediate,
and terminal resolutions against `HiDecon`, `HIDE`, tree-guided `MuSiC`,
and `Rectangle`.

## Compositional geometry and validation

### ALR reparametrisation and links to the maths vignette

Any estimator of \boldsymbol{p} must respect compositional geometry
([Aitchison 1982](#ref-aitchisonStatisticalAnalysisCompositional1982);
[Pawlowsky-Glahn and Buccianti
2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)). DeCovarT maps
unconstrained coordinates \boldsymbol{\rho}\in\mathbb{R}^{J-1} to the
open simplex via the **additive logistic** (inverse ALR) map implemented
in
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md):

p_j=\frac{e^{\rho_j}}{\sum\_{k\<J}e^{\rho_k}+1}\\(j\<J), \qquad
p_J=\frac{1}{\sum\_{k\<J}e^{\rho_k}+1}. \tag{11}

The inverse
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_log_ratio.md)
recovers \rho_j=\log(p_j/p_J). The forward and inverse maps, their
Jacobians, and Hessians used in the constrained likelihood are derived
in the companion vignette [Simplex maps and softmax / ALR
derivatives](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.html#eq-additive-logistic-intro)
(`@eq-additive-logistic-intro` and `@eq-additive-log-ratio-intro`
there). Optimisation therefore runs in \boldsymbol{\rho}-space while
reported estimates satisfy \sum_j p_j=1.

Extensions on the same footing include **ILR** or **CLR** coordinates
([Pawlowsky-Glahn and Buccianti
2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)) and Bayesian
formulations: a Dirichlet or logistic-normal prior on \boldsymbol{p}
pairs naturally with the Gaussian bulk likelihood (the experimental MAP
`CTS` path in `R/03_04_DeCovarT_estimate_CTS_MAP_Bayesian.R` is one
starting point). Regardless of coordinate system, the workflow is the
same: optimise in an unconstrained parameterisation, then map back to
the simplex for interpretation.

### Validation metrics

Extended models should be checked on **designed mixtures** with known
\boldsymbol{p}: semi-synthetic bulk profiles (as in [synthetic
scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/synthetic-scenarios.md)),
chemical standard blends, or mixture-design benchmarks with and without
injected interaction signals. Because \boldsymbol{p} is compositional,
error is often summarised with **Aitchison distance** after a log-ratio
transform ([Aitchison
1982](#ref-aitchisonStatisticalAnalysisCompositional1982)), or with MAE
/ RMSE and Pearson correlation on raw proportions (common in
deconvolution benchmarks). Comparing against mean-only baselines
(`CIBERSORT`, NNLS, second-generation single-cell-informed tools)
quantifies the gain from modelling covariance in
[Eq. 2](#eq-gaussian-convolution) and, eventually, interaction terms in
[Eq. 3](#eq-interaction-mean)–[Eq. 4](#eq-interaction-cov).

### R tooling for mixture experimental designs

The archived CRAN package `mixexp` formerly provided constrained mixture
regions, extreme-vertex designs, and ternary contour plots; archived
releases remain on CRAN but there is no maintained one-package
replacement. For validation workflows that require **designed simplex
points** rather than random \boldsymbol{p}, two practical options are:

- **`AlgDesign`** ([Wheeler 2025](#ref-R-AlgDesign)) — lattice and
  simplex-centroid mixture designs via `genMixture()` / `gen.mixture()`,
  plus D-, A-, and I-optimal designs; the closest general substitute for
  classical mixture DOE.
- **`skpr`** ([Morgan-Wall and Khoury 2025](#ref-R-skpr)) — optimal
  mixture designs from a candidate set under the \sum_j p_j=1
  constraint; preferable when optimality criteria matter more than fixed
  simplex-centroid lattices.

A modular workflow for Scheffé-type benchmarks would combine `AlgDesign`
(or `skpr`) to propose \boldsymbol{p}, DeCovarT (or baselines) to invert
the mixture, and standard compositional metrics above. Constrained
extreme-vertex designs under linear bounds remain less well served in
current R packages than the original `mixexp` workflow; custom simplex
sampling may still be required for those regions.

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Angel, Almog, Loai Naom, Shir Nabet-Levy, and Dvir Aran. 2025. ‘xCell
2.0: Robust Algorithm for Cell Type Proportion Estimation Predicts
Response to Immune Checkpoint Blockade’. *Genome Biology* 26 (1): 335.
<https://doi.org/10.1186/s13059-025-03784-3>.

Brown, L., A. N. Donev, and A. C. Bissett. 2015. ‘General Blending
Models for Data from Mixture Experiments’. *Technometrics* 57 (4):
449–56. <https://doi.org/10.1080/00401706.2014.947003>.

Chassagnol, Bastien, Grégory Nuel, and Etienne Becht. 2023. *DeCovarT, a
Multidimensional Probabilistic Model for the Deconvolution of
Heterogeneous Transcriptomic Samples*. arXiv.
<https://doi.org/10.48550/arxiv.2309.09557>.

Dong, Meichen, Aatish Thennavan, Eugene Urrutia, et al. 2021. ‘SCDC:
Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA
Sequencing References’. *Briefings in Bioinformatics* 22.
<https://doi.org/10.1093/bib/bbz166>.

Eder, Bernhard, Irene Rigato, Alexander Dietrich, et al. 2026.
*Rectangle: Robust and Scalable Multiscale Deconvolution Informed by
Single-Cell RNA Sequencing Data*. bioRxiv.
<https://doi.org/10.64898/2026.07.07.736950>.

Huang, Penghui, Manqi Cai, Xinghua Lu, Chris McKennan, and Jiebiao Wang.
2024. ‘Accurate Estimation of Rare Cell-Type Fractions from Tissue Omics
Data via Hierarchical Deconvolution’. *The Annals of Applied Statistics*
18. <https://doi.org/10.1214/23-aoas1829>.

Morgan-Wall, Tyler, and George Khoury. 2025. *Skpr: Design of
Experiments Suite: Generate and Evaluate Optimal Designs*.
<https://CRAN.R-project.org/package=skpr>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Völkl, Dennis, Malte Mensching-Buhr, Thomas Sterr, et al. 2025. ‘HIDE:
Hierarchical Cell-Type Deconvolution’. *Bioinformatics* 41
(Supplement_1): i207–16.
<https://doi.org/10.1093/bioinformatics/btaf179>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.

Wheeler, Bob. 2025. *AlgDesign: Algorithmic Experimental Design*.
<https://github.com/jvbraun/AlgDesign>.

Zaitsev, Aleksandr, Maksim Chelushkin, Daniiar Dyikanov, et al. 2022.
‘Precise Reconstruction of the TME Using Bulk RNA-seq and a Machine
Learning Algorithm Trained on Artificial Transcriptomes’. *Cancer Cell*
40. <https://doi.org/10.1016/j.ccell.2022.07.006>.
