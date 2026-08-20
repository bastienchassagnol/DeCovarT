# DeCovarT perspectives

> **Scope**
>
> Outlook on extending DeCovarT beyond the closed-reference Gaussian
> convolution. Sections cover Scheffé-type mixture structure,
> alternative observation laws and Bayesian CTS inference, sample-level
> covariates, incomplete references, isoforms, RNA–cell uncoupling,
> weighted / generalised least squares, lineage and archetypes,
> time-resolved composition, ensembles, spatial transcriptomics, and
> multi-omics. Compositional reparametrisation is implemented in
> [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
> and documented numerically in the [softmax / ALR
> derivatives](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md)
> vignette. The same catalogue is tracked in [GitHub issue
> 5](https://github.com/bastienchassagnol/DeCovarT/issues/5).

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
  deconvolution loss ([Huang, Cai, Lu, et al.
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
marker panels, and ALR optimisation ([Sec. 5.1](#sec-alr-reparam)) when
a curated lineage tree is supplied—benchmarked at broad, intermediate,
and terminal resolutions against `HiDecon`, `HIDE`, tree-guided `MuSiC`,
and `Rectangle`. Spatial `HIDF` applies the same tree idea to spots with
an `RCTD` back-end ([Zou et al.
2026](#ref-zouHidfIntegratingTreeStructured2026)); `scDETECT` uses a
lineage prior for differential expression rather than for \boldsymbol{p}
([Xu et al. 2025](#ref-xuScdetectNovelStatisticalModel2025)). Archetypes
and cell states (type versus shared gene-expression programmes) are
taken up in [Sec. 4.1](#sec-archetypes).

## Statistical perspectives

The present estimator is a **closed-reference, continuous, frequentist**
convolution: \boldsymbol{\mu}\_{\cdot j} and \boldsymbol{\Sigma}\_{j}
are plug-in moments, and only \boldsymbol{p}\_{\cdot i} is unknown. The
subsections below change the observation law, the reference, or the
loss, while keeping the simplex constraint
([Eq. 6](#eq-simplex-constraint)).

### Alternative distributions

Split first by the **support of \boldsymbol{y}** (integer counts versus
continuous intensities), then by **frequentist plug-in versus Bayesian
joint inference** of \boldsymbol{p}\_{\cdot i} and sample-level
cell-type-specific (CTS) profiles \boldsymbol{x}\_{\cdot j,i}.

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  A["Bulk observation law"] --> B["Discrete counts"]
  A --> C["Continuous intensities"]
  B --> D["Multinomial / Dirichlet: ISOpureR, BayesPrism"]
  B --> E["Poisson / PLN / ZIPLN: DeconV, Chiquet"]
  C --> F["Univariate Gaussian: DSection, DeMix, BayICE"]
  C --> G["Multivariate Gaussian convolution: DeCovarT"]
  C --> H["Log-normal convolution: BLADE"]
  D --> I["Bayesian joint p and CTS"]
  G --> J["Frequentist plug-in mu, Sigma"]
  H --> I
```

Figure 2: Taxonomy of observation models for bulk deconvolution,
relative to DeCovarT’s multivariate Gaussian convolution.

#### Discrete counts

A multinomial (or Dirichlet–multinomial) model treats
\boldsymbol{y}\_{\cdot i} as integer reads that sum to a library size.
`ISOpureR` implements a hierarchical multinomial–Dirichlet purification
of mixed tumour profiles ([Anghel et al.
2015](#ref-anghelISOpureRImplementationComputational2015)). `BayesPrism`
places a multinomial likelihood on bulk counts with a patient-derived
scRNA-seq prior and jointly infers \boldsymbol{p}\_{\cdot i} and
sample-level CTS expression ([Chu et al.
2022](#ref-chuCellTypeGeneExpression2022)). `DeconV` instead uses
Poisson aggregation from single-cell to bulk (sums of independent
Poissons remain Poisson) and returns interval estimates of
\boldsymbol{p}\_{\cdot i} ([Gynter et al.
2023](#ref-gynterDeconvProbabilisticCellType2023)).

Over-dispersed, correlated counts are the domain of the multivariate
**Poisson–log-normal (PLN)** family: a latent Gaussian vector induces
dependence, and a Poisson observation layer yields integer reads
([Chiquet et al. 2021](#ref-chiquetPoissonLognormalModelVersatile2021),
[2018](#ref-chiquetVariationalInferenceSparse2018)). **ZIPLN** adds a
Bernoulli zero-inflation layer for dropout ([Batardière et al.
2025](#ref-batardiereZeroInflationMultivariatePoisson2025)). A
DeCovarT-style extension would replace [Eq. 2](#eq-gaussian-convolution)
by a PLN (or ZIPLN) convolution on the latent log-abundance scale,
keeping the ALR map on \boldsymbol{p}.

#### Continuous intensities: frequentist versus Bayesian CTS

DSection is the univariate ancestor of DeCovarT: a Gaussian convolution
that is heteroscedastic across genes, with Normal / Gamma / Dirichlet
priors and a Gibbs–Metropolis sampler ([Erkkilä et al.
2010](#ref-erkkilaProbabilisticAnalysisGene2010)). DeMix and DeMixT
remain frequentist convolutions that treat tumour or stromal profiles as
latent ([Ahn et al. 2013](#ref-ahnDeMixDeconvolutionMixed2013); [Wang et
al. 2018](#ref-wangTranscriptomeDeconvolutionHeterogeneous2018)).
DeCovarT keeps that convolution but makes it **multivariate** with
sparse \boldsymbol{\Sigma}\_{j}.

`BLADE` replaces the Gaussian with a **gene-wise log-normal**
convolution, jointly purifying CTS profiles by variational inference
([Andrade Barbosa et al.
2021](#ref-andrade-barbosaBayesianLogNormalDeconvolution2021)). `bMIND`
estimates sample-level CTS with an scRNA-seq-derived prior ([Wang et al.
2020](#ref-wangBayesianEstimationCelltypespecific2020)). The
experimental MAP path in `R/03_04_DeCovarT_estimate_CTS_MAP_Bayesian.R`
is the corresponding DeCovarT starting point: recover
\boldsymbol{x}\_{\cdot j,i} given \boldsymbol{y}\_{\cdot i} and
\boldsymbol{p}\_{\cdot i} under [Eq. 2](#eq-gaussian-convolution) rather
than under independent genes.

Zhang et al. benchmark joint \boldsymbol{p} + CTS engines and report
**BayesPrism with DWLS gene weights** as the strongest combination on
pseudobulk and real bulk data ([Zhang et al.
2026](#ref-zhangIntegratedInferenceCellularCompositions2026)). That is a
weighted-likelihood analogue of [Sec. 2.6](#sec-robust-gls).

### Sample-level covariates

A condition, tissue, batch, sex or age vector \boldsymbol{z}\_{i} is
indexed by sample i, whereas \boldsymbol{\mu} has genes in rows and cell
types in columns. Covariates therefore cannot be appended as extra
columns of the signature. Two generative extensions are distinct ([Fan
et al. 2022](#ref-fanMusic2CellTypeDeconvolution2022)).

**State model.** Condition alters the reference distribution,

\boldsymbol{x}\_{\cdot j,i}\mid\boldsymbol{z}\_{i} \sim\mathcal{N}\_G
\bigl(\boldsymbol{\mu}\_{j}(\boldsymbol{z}\_{i}),\boldsymbol{\Sigma}\_{j}(\boldsymbol{z}\_{i})\bigr),
\tag{11}

for example
\boldsymbol{\mu}\_{j}(\boldsymbol{z}\_{i})=\boldsymbol{\mu}\_{j}+B\_{j}\boldsymbol{z}\_{i}.
The bulk law is then [Eq. 2](#eq-gaussian-convolution) with those
condition-dependent moments. `MuSiC2` addresses the same mismatch by
iteratively dropping cell-type-specific DE genes between the bulk
condition and the scRNA-seq reference; it remains closed-reference and
cannot invent novel types ([Fan et al.
2022](#ref-fanMusic2CellTypeDeconvolution2022)).

**Composition model.** Condition alters \boldsymbol{p}\_{\cdot i}
through the existing ALR coordinates \boldsymbol{\rho}\_{i},

\boldsymbol{\rho}\_{i}=B\boldsymbol{z}\_{i}+\boldsymbol{u}\_{i}, \qquad
\boldsymbol{p}\_{\cdot i}=\psi(\boldsymbol{\rho}\_{i}), \tag{12}

with \psi the additive logistic map ([Eq. 16](#eq-alr-forward)). This is
the setting in which a future
[`predict()`](https://rdrr.io/r/stats/predict.html) method on new
\boldsymbol{z} would be meaningful.

A known sequencing exposure s\_{i} is a **multiplicative** scale,
\boldsymbol{y}\_{\cdot
i}\sim\mathcal{N}\_G(s\_{i}\boldsymbol{\mu}\boldsymbol{p}\_{\cdot
i},s\_{i}^{2}\Sigma(\boldsymbol{p}\_{\cdot i})), not an additive `lm`
offset. Cell-type transcriptome size r\_{j} is a different parameter
([Sec. 2.5](#sec-uncoupling)). A gene-wise affine standardisation
applied identically to \boldsymbol{y}, \boldsymbol{\mu} and every
\boldsymbol{\Sigma}\_{j} leaves \hat{\boldsymbol{p}} unchanged in exact
arithmetic; column-wise [`scale()`](https://rdrr.io/r/base/scale.html)
on \boldsymbol{\mu} does not.

### Incomplete references

Closed-reference DeCovarT assumes every abundant type appears as a
column of \boldsymbol{\mu}. An unrepresented population contributes a
structured gene vector \boldsymbol{u}\_{i}, not a scalar intercept:

\boldsymbol{y}\_{\cdot i} =\boldsymbol{\mu}\\\boldsymbol{p}\_{\cdot
i}^{\mathrm{known}}
+p\_{u,i}\boldsymbol{u}\_{i}+\boldsymbol{\varepsilon}\_{i}, \qquad
\sum\_{j}p\_{ji}+p\_{u,i}=1. \tag{13}

`DICEPro` shows that supervised engines remain stable while most types
are present and then collapse as \boldsymbol{\mu} becomes incomplete,
especially when the missing type is abundant; it wraps existing methods
by adjusting reference signatures ([Ba et al.
2026](#ref-baWhenLessNot2026)). `BayICE` is a univariate-Gaussian
Bayesian semi-reference model that estimates known types **and one
unknown type**, with spike-and-slab selection of genes and columns of
\boldsymbol{\mu} ([Tai et al.
2021](#ref-taiBayiceBayesianHierarchicalModel2021)). Montierth et al.
convolve negative binomials for known components plus an unknown tumour
component, assuming proportions shift means but not dispersions
([Montierth et al.
2025](#ref-montierthDeconvolutionSparsecountRNA2025)). `CDState`
reconstructs malignant **states** rather than a single tumour column
([Kraft et al. 2026](#ref-kraftCdstateResolvesMalignantCell2026)).
`EPIC` already includes an uncharacterised compartment in its signature
([Racle et al. 2017](#ref-racleSimultaneousEnumerationCancer2017)).
Semi-CAM uses partial marker information ([Dong et al.
2020](#ref-dongSemiCAMSemisupervisedDeconvolution2020)); `BLEND`
automates reference-panel selection ([Huang, Cai, McKennan, et al.
2024](#ref-huangBLENDProbabilisticCellular2024)).

### Isoform-level observations

Gene-level \boldsymbol{y}\in\mathbb{R}\_{+}^{G} discards differential
isoform usage. Expanding the observation to transcripts t=1,\ldots,T
(short-read quantification or long-read) yields a taller signature
\boldsymbol{\mu}^{\mathrm{iso}}\in\mathbb{R}^{T\times J}. `IsoDeconvMM`
estimates \boldsymbol{p} from isoform-level expression, even from a
single gene, by exploiting differential isoform usage rather than
gene-level DE; full-length long-read or spatial RNA-seq are the natural
sources of cell-type- and isoform-specific references when droplet
scRNA-seq misses isoforms ([Heiling et al.
2023](#ref-heilingEstimatingCellTypeComposition2023)).

### RNA fraction versus cell fraction

DeCovarT, like most linear engines, estimates the **RNA mass fraction**
of type j, not the cytometric cell fraction, unless transcriptome sizes
r\_{j} (mean transcripts per cell) are homogeneous. With \hat p\_{j} the
RNA-scale estimate,

\hat p\_{j}^{\ast} =\frac{\hat p\_{j}/r\_{j}}{\sum\_{k=1}^{J}\hat
p\_{k}/r\_{k}}. \tag{14}

**Post-correction** treats r\_{j} as measured (or proxied). `EPIC` and
`quanTIseq` rescale \hat{\boldsymbol{p}} using kit-based mRNA content or
housekeeping-gene / proteasome-subunit surrogates ([Racle et al.
2017](#ref-racleSimultaneousEnumerationCancer2017); [Finotello et al.
2019](#ref-finotello_etal19)). `MuSiC` uses average cell-type library
size as a size proxy and warns that TPM can discard the information
needed to recover cell fractions ([Wang et al.
2019](#ref-wangBulkTissueCell2019)). `ReDeconv` separates global
library-size normalisation from cell-type transcriptome size ([Lu et al.
2025](#ref-luTranscriptomeSizeMatters2025)).

**Joint estimation** treats r\_{j} as unknown. `MMAD` estimates
extraction efficiencies by non-linear conjugate gradients, so the
mixture is no longer linear ([Liebner et al.
2014](#ref-liebnerMMADMicroarrayMicrodissection2014)). A DeCovarT
analogue would introduce r\_{j} inside the mean
\sum\_{j}p\_{ji}r\_{j}\boldsymbol{\mu}\_{\cdot j} (with a
product-simplex constraint on (p\_{j}r\_{j})), rather than assuming
r\_{j} known.

### Robust regression, GLS, and gene weights

Let W=\mathrm{diag}(w\_{1},\ldots,w\_{G}) be gene-specific precision
weights (DWLS dampening, voom mean–variance weights, or inverse
leverage). Weighted least squares minimises (\boldsymbol{y}\_{\cdot
i}-\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i})^{\top}
W(\boldsymbol{y}\_{\cdot i}-\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i})
([Tsoucas et al. 2019](#ref-tsoucasAccurateEstimationCelltype2019); [Law
et al. 2014](#ref-lawVoomPrecisionWeights2014)). CIBERSORT’s \nu-SVR
replaces squared error by a robust hinge ([Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)).

**Generalised least squares** uses a single residual covariance \Sigma,
not a cell-type tensor. The equi-balanced approximation p\_{j}=1/J of
the DeCovarT covariance is
\Sigma=J^{-2}\sum\_{j}\boldsymbol{\Sigma}\_{j}, and

\hat{\boldsymbol{p}}^{\mathrm{GLS}}
=(\boldsymbol{\mu}^{\top}\Theta\boldsymbol{\mu})^{-1}\boldsymbol{\mu}^{\top}\Theta\boldsymbol{y}\_{\cdot
i}, \qquad \Theta=\Sigma^{-1}, \tag{15}

after which the simplex is imposed by ALR or projection. WLS is the
diagonal special case \Theta=W. DeCovarT is strictly richer:
\Sigma(\boldsymbol{p})=\sum\_{j}p\_{ji}^{2}\boldsymbol{\Sigma}\_{j}
depends on the unknown coefficients.

A direct extension, following the BayesPrism–DWLS construction ([Zhang
et al. 2026](#ref-zhangIntegratedInferenceCellularCompositions2026)), is
to ingest W into the convolution: transform
\boldsymbol{y}^{\star}=W^{1/2}\boldsymbol{y}\_{\cdot i},
\boldsymbol{\mu}^{\star}=W^{1/2}\boldsymbol{\mu},
\boldsymbol{\Sigma}\_{j}^{\star}=W^{1/2}\boldsymbol{\Sigma}\_{j}W^{1/2}
and run the existing solver — weighted GLS with cell-type-specific
networks retained.

`cellGeometry` estimates uncertainty from gene-wise heteroscedasticity
and the cross-sample variance of each gene, i.e. a **global**
gene/sample covariance rather than \boldsymbol{\Sigma}\_{j} ([Lau et al.
2026](#ref-lauCellGeometryUltrafastSinglecell2026)). `Unico` deconvolves
a 2-D bulk matrix into a 3-D sample \times feature \times cell-type
tensor and explicitly models cell-type-level covariances, including on
methylation ([Chen et al. 2025](#ref-chenUnicoUnifiedModelCell2025)).
Interaction monomials p\_{j}p\_{k} in the mean
([Eq. 3](#eq-interaction-mean)) are the Scheffé / general-linear-model
route; `DecOT` instead changes the discrepancy to an optimal-transport
loss ([Liu et al. 2022](#ref-liuDecotBulkDeconvolutionOptimal2022)).

### Time-resolved composition

When J\>G the linear map is under-determined. Elastic-net methods such
as `DCQ` (`glmnet`) select the **support of cell types that change**
between biological states, rather than recovering a full simplex
snapshot ([Altboum et al.
2014](#ref-altboumDigitalCellQuantification2014); [Friedman et al.
2025](#ref-R-glmnet)). Static benchmarks that score these methods as
poor deconvolution miss that target ([Jin and Liu
2021](#ref-jinBenchmarkRNAseqDeconvolution2021)). `DCQ` tracked 213
immune subtypes across ten influenza time points and reported changes in
about 70 types; `ImmQuant` packages a similar pipeline ([Frishberg et
al. 2016](#ref-frishbergImmQuantUserfriendlyTool2016)).

A DeCovarT trajectory would put a temporal prior on ALR coordinates,
e.g. \boldsymbol{\rho}\_{i}(t)=\boldsymbol{\rho}\_{i}(t-1)+\boldsymbol{\eta}\_{t},
as in `ChronoStrain`’s Bayesian strain trajectories ([Kim et al.
2025](#ref-kimLongitudinalProfilingLowAbundance2025)). `scTREND`
supplies an annotation-free single-cell hazard clock for organoid /
gastruloid time courses ([Yuki et al.
2026](#ref-yukiSctrendAnnotationFreeSingle2026)). `SpaDecoder` adds
space–time kernels on spots ([Lobo et al.
2026](#ref-loboSpatiotemporalCellTypeDeconvolution2026)).

### Ensembles

Two operations are easily confused. **Several references, one
algorithm:** `MuSiC` weights subjects and `SCDC` weights scRNA-seq
studies ([Wang et al. 2019](#ref-wangBulkTissueCell2019); [Dong et al.
2021](#ref-dongSCDCBulkGene2021)). **Several algorithms, one consensus
\hat{\boldsymbol{p}}:** `EnsDeconv` fuses 11 bulk engines plus variation
in references, markers and normalisation by cell-type-specific robust
regression, on 4,937 samples with measured fractions ([Cai et al.
2022](#ref-caiRobustAccurateEstimation2022)). `EnDecon` runs 14 spatial
(and bulk) methods and forms a weighted-median consensus, down-weighting
discordant engines ([Tu et al.
2023](#ref-tuEndeconCellTypeDeconvolution2023)). The DREAM community
assessment found that a mean-rank ensemble beat every individual method
([White et al.
2024](#ref-whiteCommunityAssessmentMethodsDeconvolve2024)).

## Spatial transcriptomics

Each spot (or pixel) s is a local mixture \boldsymbol{y}(s) with its own
\boldsymbol{p}(s)\in\Delta^{J-1}. The Gaussian convolution still applies
per location; space enters as dependence among neighbouring
\boldsymbol{p}(s). Sequencing-based surveys classify probabilistic
assumptions and provide practical benchmarks ([Saqib and Kim
2025](#ref-saqibPixelsCellTypes2025); [Gaspard-Boulinc et al.
2025](#ref-gaspard-boulincCelltypeDeconvolutionMethods2025); [Li et al.
2023](#ref-liComprehensiveBenchmarkingPractical2023)). Imaging-based
reconstruction goes the other way: `HistoMap` generates single-cell-like
profiles from bulk with a variational autoencoder and maps them onto
histological coordinates with an H-ViT ([He et al.
2026](#ref-heHistomapReconstructingSpatiallyResolved2026)); `HEDeST`
pairs H&E with spot-level \hat{\boldsymbol{p}}(s) ([Gortana et al.
2026](#ref-gortanaHedestIntegrativeApproachEnhance2026)). `SpaDecoder`
aligns slices, infers neighbourhoods, and uses 3-D Gaussian kernels,
with an explicit simplex constraint and a penalty on the number of types
per spot ([Lobo et al.
2026](#ref-loboSpatiotemporalCellTypeDeconvolution2026)). Under
**cell-type mismatch**, missing types are mostly absorbed by the most
similar column of \boldsymbol{\mu} ([Mahamune et al.
2025](#ref-mahamuneSystematicEvaluationRobustnessDeconvolution2025)) —
the spatial counterpart of [Sec. 2.3](#sec-semi-reference).

## Multi-modal and multi-omics

Bulk deconvolution is one arm of a multi-scale tumour-immune map that
also includes scRNA-seq and spatial assays ([Sun et al.
2026](#ref-sunMultiScaleTranscriptomicsRedefining2026)). `DECODE` trains
a shared deconvolution architecture across omics ([Zhao et al.
2026](#ref-zhaoDecodeDeepLearningBased2026)). Yao et al. reconstruct
cell-type-specific regulatory processes from paired ATAC-seq and bulk
RNA-seq, benchmarking against CIBERSORTx ([Yao et al.
2026](#ref-yaoHighResolutionReconstructionCell2026)).
Proteomics-constrained deconvolution uses protein abundance as a
regulariser on transcriptomic \boldsymbol{p} / CTS ([Işık et al.
2026](#ref-isikProteomicsConstrainedDeconvolutionReveals2026)). `Unico`
already treats expression and DNA methylation under one tensor model
([Chen et al. 2025](#ref-chenUnicoUnifiedModelCell2025)). The HADACA3
community benchmark stresses the limits of naive multimodal
concatenation ([Barbot and Richard
2026](#ref-barbotPromisesLimitsMultimodal2026)).

### Archetypes, states, and potency

A finite J-column signature is a piecewise approximation of within-type
continua. **Archetypal analysis** places cells on a simplicial polytope
whose vertices are extreme programmes, avoiding NMF/ICA rotational
ambiguity ([Hart et al. 2015](#ref-hartInferringBiologicalTasks2015);
[Crowley et al. 2026](#ref-crowleyParetoOptimalityRevealsAtlas2026)).
`ACTION` separates transcriptional identity (cell type) from activity
states (shared gene-expression programmes) and chooses k automatically
([Mohammadi et al.
2020](#ref-mohammadiAMultiresolutionFrameworkCharacterize2020));
`scAAnet` extends that construction with a VAE and a zero-inflated
negative-binomial reconstruction ([Wang and Zhao
2022](#ref-wangNonLinearArchetypalAnalysis2022)). `tissueResolver`
builds a virtual tissue without predefined labels, aiming at fine states
within types ([Simeth et al.
2024](#ref-simethVirtualTissueExpressionAnalysis2024)). Potency and
differentiation scores, as in hypothalamic–pituitary organoids, replace
discrete labels by a continuous axis ([Asano et al.
2024](#ref-asanoADeepLearningApproach2024)). Macrophage M1/M2
polarisation is the standard biological example of a
microenvironment-dependent continuum rather than two extra columns of
\boldsymbol{\mu}.

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
p_J=\frac{1}{\sum\_{k\<J}e^{\rho_k}+1}. \tag{16}

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

Ahn, Jaeil, Ying Yuan, Giovanni Parmigiani, et al. 2013. ‘DeMix:
Deconvolution for Mixed Cancer Transcriptomes Using Raw Measured Data’.
*Bioinformatics (Oxford, England)* 29.
<https://doi.org/10.1093/bioinformatics/btt301>.

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Altboum, Zeev, Yael Steuerman, Eyal David, et al. 2014. ‘Digital Cell
Quantification Identifies Global Immune Cell Dynamics During Influenza
Infection’. *Molecular Systems Biology* 10.
<https://doi.org/10.1002/msb.134947>.

Andrade Barbosa, Bárbara, Saskia D. van Asten, Ji Won Oh, et al. 2021.
‘Bayesian Log-Normal Deconvolution for Enhanced in Silico
Microdissection of Bulk Gene Expression Data’. *Nature Communications*
12 (1). <https://doi.org/10.1038/s41467-021-26328-2>.

Angel, Almog, Loai Naom, Shir Nabet-Levy, and Dvir Aran. 2025. ‘xCell
2.0: Robust Algorithm for Cell Type Proportion Estimation Predicts
Response to Immune Checkpoint Blockade’. *Genome Biology* 26 (1): 335.
<https://doi.org/10.1186/s13059-025-03784-3>.

Anghel, Catalina V, Gerald Quon, Syed Haider, et al. 2015. ‘ISOpureR: An
R Implementation of a Computational Purification Algorithm of Mixed
Tumour Profiles’. *BMC Bioinformatics* 16.
<https://doi.org/10.1186/s12859-015-0597-x>.

Asano, Tomoyoshi, Hidetaka Suga, Hirohiko Niioka, et al. 2024. ‘A Deep
Learning Approach to Predict Differentiation Outcomes in
Hypothalamic-Pituitary Organoids’. *Communications Biology* 7 (1).
<https://doi.org/10.1038/s42003-024-07109-1>.

Ba, Kalidou, Rodolphe Thiébaut, Xavier Hinaut, and Boris Hejblum. 2026.
*When Less Is Not More: DICEPro Mitigates the Impact of Incomplete
Reference Matrices on Cellular Frequency Deconvolution*. bioRxiv.
<https://doi.org/10.64898/2026.06.17.732876>.

Barbot, Hugo, and Magali Richard. 2026. ‘On the Promises and Limits of
Multimodal Integration for Deconvolution: The HADACA3 Benchmark’.
*NeurIPS*.

Batardière, Bastien, Julien Chiquet, François Gindraud, and Mahendra
Mariadassou. 2025. ‘Zero-Inflation in the Multivariate Poisson Lognormal
Family’. *Statistics and Computing* 35 (6).
<https://doi.org/10.1007/s11222-025-10729-0>.

Brown, L., A. N. Donev, and A. C. Bissett. 2015. ‘General Blending
Models for Data from Mixture Experiments’. *Technometrics* 57 (4):
449–56. <https://doi.org/10.1080/00401706.2014.947003>.

Cai, Manqi, Molin Yue, Tianmeng Chen, et al. 2022. ‘Robust and Accurate
Estimation of Cellular Fraction from Tissue Omics Data via Ensemble
Deconvolution’. *Bioinformatics* 38.
<https://doi.org/10.1093/bioinformatics/btac279>.

Chassagnol, Bastien, Grégory Nuel, and Etienne Becht. 2023. *DeCovarT, a
Multidimensional Probabilistic Model for the Deconvolution of
Heterogeneous Transcriptomic Samples*. arXiv.
<https://doi.org/10.48550/arxiv.2309.09557>.

Chen, Zeyuan Johnson, Elior Rahmani, and Eran Halperin. 2025. ‘Unico: A
Unified Model for Cell-Type Resolution Genomics from Heterogeneous Omics
Data’. *Genome Biology* 26 (1).
<https://doi.org/10.1186/s13059-025-03776-3>.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2018.
*Variational Inference for Sparse Network Reconstruction from Count
Data*. arXiv. <https://doi.org/10.48550/arxiv.1806.03120>.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2021. ‘The
Poisson-Lognormal Model as a Versatile Framework for the Joint Analysis
of Species Abundances’. *Frontiers in Ecology and Evolution* 9.
<https://doi.org/10.3389/fevo.2021.588292>.

Chu, Tinyi, Zhong Wang, Dana Pe’er, and Charles G. Danko. 2022. ‘Cell
Type and Gene Expression Deconvolution with BayesPrism Enables Bayesian
Integrative Analysis Across Bulk and Single-Cell RNA Sequencing in
Oncology’. *Nature Cancer* 3 (4): 505–17.
<https://doi.org/10.1038/s43018-022-00356-3>.

Crowley, George, Uri Alon, and Stephen R. Quake. 2026. ‘Pareto
Optimality Reveals an Atlas of Cellular Archetypes’. *Proceedings of the
National Academy of Sciences* 123 (11).
<https://doi.org/10.1073/pnas.2530194123>.

Dong, Li, Avinash Kollipara, Toni Darville, Fei Zou, and Xiaojing Zheng.
2020. ‘Semi-CAM: A Semi-Supervised Deconvolution Method for Bulk
Transcriptomic Data with Partial Marker Gene Information’. *Scientific
Reports* 10. <https://doi.org/10.1038/s41598-020-62330-2>.

Dong, Meichen, Aatish Thennavan, Eugene Urrutia, et al. 2021. ‘SCDC:
Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA
Sequencing References’. *Briefings in Bioinformatics* 22.
<https://doi.org/10.1093/bib/bbz166>.

Eder, Bernhard, Irene Rigato, Alexander Dietrich, et al. 2026.
*Rectangle: Robust and Scalable Multiscale Deconvolution Informed by
Single-Cell RNA Sequencing Data*. bioRxiv.
<https://doi.org/10.64898/2026.07.07.736950>.

Erkkilä, Timo, Saara Lehmusvaara, Pekka Ruusuvuori, Tapio Visakorpi,
Ilya Shmulevich, and Harri Lähdesmäki. 2010. ‘Probabilistic Analysis of
Gene Expression Measurements from Heterogeneous Tissues’.
*Bioinformatics* 26. <https://doi.org/10.1093/bioinformatics/btq406>.

Fan, Jiaxin, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, and Rui
Xiao. 2022. ‘MuSiC2: Cell-Type Deconvolution for Multi-Condition Bulk
RNA-seq Data’. *Briefings in Bioinformatics* 23 (6).
<https://doi.org/10.1093/bib/bbac430>.

Finotello, Francesca, Clemens Mayer, Christina Plattner, et al. 2019.
‘Molecular and Pharmacological Modulators of the Tumor Immune Contexture
Revealed by Deconvolution of RNA-seq Data’. *Genome Medicine* 11.
<https://doi.org/10.1186/s13073-019-0638-6>.

Friedman, Jerome, Trevor Hastie, Rob Tibshirani, et al. 2025. *Glmnet:
Lasso and Elastic-Net Regularized Generalized Linear Models*.
<https://glmnet.stanford.edu>.

Frishberg, Amit, Avital Brodt, Yael Steuerman, and Irit Gat-Viks. 2016.
‘ImmQuant: A User-Friendly Tool for Inferring Immune Cell-Type
Composition from Gene-Expression Data’. *Bioinformatics* 32.
<https://doi.org/10.1093/bioinformatics/btw535>.

Gaspard-Boulinc, Lucie C., Luca Gortana, Thomas Walter, Emmanuel
Barillot, and Florence M. G. Cavalli. 2025. ‘Cell-Type Deconvolution
Methods for Spatial Transcriptomics’. *Nature Reviews Genetics* 26.
<https://doi.org/10.1038/s41576-025-00845-y>.

Gortana, Luca, Loïc Chadoutaud, Raphaël Bourgade, Emmanuel Barillot, and
Thomas Walter. 2026. *HEDeST: An Integrative Approach to Enhance Spatial
Transcriptomic Deconvolution with Histology*. bioRxiv.
<https://doi.org/10.64898/2026.01.06.697922>.

Gynter, Artur, Dimitri Meistermann, Harri Lähdesmäki, and Helena
Kilpinen. 2023. *DeconV: Probabilistic Cell Type Deconvolution from Bulk
RNA-sequencing Data*. bioRxiv.
<https://doi.org/10.1101/2023.12.07.570524>.

Hart, Yuval, Hila Sheftel, Jean Hausser, et al. 2015. ‘Inferring
Biological Tasks Using Pareto Analysis of High-Dimensional Data’.
*Nature Methods* 12 (3): 233–35. <https://doi.org/10.1038/nmeth.3254>.

He, Jia, Yong Cao, Yan Liu, et al. 2026. ‘HistoMap: Reconstructing
Spatially Resolved Single-Cell Profiles from Bulk RNA-Seq to Decipher
the Immune-Excluded Microenvironment in Colon Cancer’. *International
Journal of Molecular Sciences* 27 (12): 5259.
<https://doi.org/10.3390/ijms27125259>.

Heiling, Hillary M., Douglas R. Wilson, Naim U. Rashid, Wei Sun, and
Joseph G. Ibrahim. 2023. ‘Estimating Cell Type Composition Using Isoform
Expression One Gene at a Time’. *Biometrics* 79 (2): 854–65.
<https://doi.org/10.1111/biom.13614>.

Huang, Penghui, Manqi Cai, Xinghua Lu, Chris McKennan, and Jiebiao Wang.
2024. ‘Accurate Estimation of Rare Cell-Type Fractions from Tissue Omics
Data via Hierarchical Deconvolution’. *The Annals of Applied Statistics*
18. <https://doi.org/10.1214/23-aoas1829>.

Huang, Penghui, Manqi Cai, Chris McKennan, and Jiebiao Wang. 2024.
*BLEND: Probabilistic Cellular Deconvolution with Automated Reference
Selection*. bioRxiv. <https://doi.org/10.1101/2024.08.02.606458>.

Işık, Esra Büşra, Michael J. Haley, Ali Hussein Al-Anbaki, et al. 2026.
*Proteomics-Constrained Deconvolution Reveals Spatial Cell-Type Programs
in Tumours*. bioRxiv. <https://doi.org/10.64898/2026.06.01.729268>.

Jin, Haijing, and Zhandong Liu. 2021. ‘A Benchmark for RNA-seq
Deconvolution Analysis Under Dynamic Testing Environments’. *Genome
Biology* 22. <https://doi.org/10.1186/s13059-021-02290-6>.

Kim, Younhun, Colin J. Worby, Sawal Acharya, et al. 2025. ‘Longitudinal
Profiling of Low-Abundance Strains in Microbiomes with ChronoStrain’.
*Nature Microbiology* 10 (5): 1184–97.
<https://doi.org/10.1038/s41564-025-01983-z>.

Kraft, Agnieszka, Josephine Yates, Florian Barkmann, and Valentina
Boeva. 2026. ‘CDState Resolves Malignant Cell Heterogeneity from Bulk
Tumor RNA-Sequencing Data’. *Cancer Research*, ahead of print.
<https://doi.org/10.1158/0008-5472.can-25-4102>.

Lau, Rachel, Cankut Çubuk, Athina Spiliopoulou, et al. 2026.
*cellGeometry: Ultra-Fast Single-Cell Deconvolution of Bulk RNA-Seq
Using a Geometric Solution*. bioRxiv.
<https://doi.org/10.64898/2026.01.24.701240>.

Law, Charity W., Yunshun Chen, Wei Shi, and Gordon K. Smyth. 2014.
‘Voom: Precision Weights Unlock Linear Model Analysis Tools for RNA-seq
Read Counts’. *Genome Biology* 15.
<https://doi.org/10.1186/gb-2014-15-2-r29>.

Li, Haoyang, Juexiao Zhou, Zhongxiao Li, et al. 2023. ‘A Comprehensive
Benchmarking with Practical Guidelines for Cellular Deconvolution of
Spatial Transcriptomics’. *Nature Communications* 14.
<https://doi.org/10.1038/s41467-023-37168-7>.

Liebner, David A., Kun Huang, and Jeffrey D. Parvin. 2014. ‘MMAD:
Microarray Microdissection with Analysis of Differences Is a
Computational Tool for Deconvoluting Cell Type-Specific Contributions
from Tissue Samples’. *Bioinformatics (Oxford, England)* 30.
<https://doi.org/10.1093/bioinformatics/btt566>.

Liu, Gan, Xiuqin Liu, and Liang Ma. 2022. ‘DecOT: Bulk Deconvolution
with Optimal Transport Loss Using a Single-Cell Reference’. *Frontiers
in Genetics* 13. <https://doi.org/10.3389/fgene.2022.825896>.

Lobo, Macrina Maria, Ziqi Zhang, and Xiuwei Zhang. 2026. *Spatiotemporal
Cell Type Deconvolution Leveraging Tissue Structure*. bioRxiv.
<https://doi.org/10.64898/2026.02.10.705204>.

Lu, Songjian, Jiyuan Yang, Lei Yan, et al. 2025. ‘Transcriptome Size
Matters for Single-Cell RNA-seq Normalization and Bulk Deconvolution’.
*Nature Communications* 16 (1).
<https://doi.org/10.1038/s41467-025-56623-1>.

Mahamune, Utkarsh M., Aldo Jongejan, Antoine H. C. van Kampen, Lisa G.
M. van Baarsen, and Perry D. Moerland. 2025. *Systematic Evaluation of
Robustness of Deconvolution Methods for Spatial Transcriptomics Data in
Case of Cell Type Mismatch*. bioRxiv.
<https://doi.org/10.1101/2025.08.12.669903>.

Mohammadi, Shahin, Jose Davila-Velderrain, and Manolis Kellis. 2020. ‘A
Multiresolution Framework to Characterize Single-Cell State Landscapes’.
*Nature Communications* 11 (1).
<https://doi.org/10.1038/s41467-020-18416-6>.

Montierth, Matthew D., Hao Yan, Liyang Xie, et al. 2025. *Deconvolution
of Sparse-count RNA Sequencing Data for Tumor Cells Using Embedded
Negative Binomial Distributions*. bioRxiv.
<https://doi.org/10.1101/2025.11.21.689822>.

Morgan-Wall, Tyler, and George Khoury. 2025. *Skpr: Design of
Experiments Suite: Generate and Evaluate Optimal Designs*.
<https://CRAN.R-project.org/package=skpr>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Racle, Julien, Kaat de Jonge, Petra Baumgaertner, Daniel E Speiser, and
David Gfeller. 2017. ‘Simultaneous Enumeration of Cancer and Immune Cell
Types from Bulk Tumor Gene Expression Data’. *eLife* 6.
<https://doi.org/10.7554/elife.26476>.

Saqib, Jahanzeb, and Junil Kim. 2025. ‘From Pixels to Cell Types: A
Comprehensive Review of Computational Methods for Spatial
Transcriptomics Deconvolution’. *Genomics & Informatics* 23.
<https://doi.org/10.1186/s44342-025-00055-2>.

Simeth, Jakob, Paul Hüttl, Marian Schön, et al. 2024. ‘Virtual Tissue
Expression Analysis’. *Bioinformatics* 40 (12).
<https://doi.org/10.1093/bioinformatics/btae709>.

Sun, Jing, Yingxue Xiao, Lingling Xie, et al. 2026. ‘Multi-Scale
Transcriptomics Redefining the Tumor Immune Microenvironment’. *BioTech*
15 (1): 7. <https://doi.org/10.3390/biotech15010007>.

Tai, An-Shun, George C. Tseng, and Wen-Ping Hsieh. 2021. ‘BayICE: A
Bayesian Hierarchical Model for Semireference-Based Deconvolution of
Bulk Transcriptomic Data’. *The Annals of Applied Statistics* 15 (1).
<https://doi.org/10.1214/20-aoas1376>.

Tsoucas, Daphne, Rui Dong, Haide Chen, Qian Zhu, Guoji Guo, and
Guo-Cheng Yuan. 2019. ‘Accurate Estimation of Cell-Type Composition from
Gene Expression Data’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-019-10802-z>.

Tu, Jia-Juan, Hui-Sheng Li, Hong Yan, and Xiao-Fei Zhang. 2023.
‘EnDecon: Cell Type Deconvolution of Spatially Resolved Transcriptomics
Data via Ensemble Learning’. *Bioinformatics* 39 (1).
<https://doi.org/10.1093/bioinformatics/btac825>.

Völkl, Dennis, Malte Mensching-Buhr, Thomas Sterr, et al. 2025. ‘HIDE:
Hierarchical Cell-Type Deconvolution’. *Bioinformatics* 41
(Supplement_1): i207–16.
<https://doi.org/10.1093/bioinformatics/btaf179>.

Wang, Jiebiao, Kathryn Roeder, and Bernie Devlin. 2020. *Bayesian
Estimation of Cell-Type-Specific Gene Expression Per Bulk Sample with
Prior Derived from Single-Cell Data*. openRxiv; openRxiv.
<https://doi.org/10.1101/2020.08.05.238949>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.

Wang, Yuge, and Hongyu Zhao. 2022. ‘Non-Linear Archetypal Analysis of
Single-Cell RNA-seq Data by Deep Autoencoders’. *PLOS Computational
Biology* 18 (4): e1010025.
<https://doi.org/10.1371/journal.pcbi.1010025>.

Wang, Zeya, Shaolong Cao, Jeffrey S. Morris, et al. 2018. ‘Transcriptome
Deconvolution of Heterogeneous Tumor Samples with Immune Infiltration’.
*iScience* 9. <https://doi.org/10.1016/j.isci.2018.10.028>.

Wheeler, Bob. 2025. *AlgDesign: Algorithmic Experimental Design*.
<https://github.com/jvbraun/AlgDesign>.

White, Brian S., Aurélien de Reyniès, Aaron M. Newman, et al. 2024.
‘Community Assessment of Methods to Deconvolve Cellular Composition from
Bulk Gene Expression’. *Nature Communications* 15 (1).
<https://doi.org/10.1038/s41467-024-50618-0>.

Xu, Yuhan, Weiwei Zhang, and Hao Wu. 2025. ‘scDETECT: A Novel
Statistical Model Accounting for Cell Type Correlation in Single-Cell
RNA-seq Differential Expression Analysis’. *Briefings in Bioinformatics*
26 (5). <https://doi.org/10.1093/bib/bbaf556>.

Yao, Li, Sagar R. Shah, Abdullah Ozer, et al. 2026. ‘High-Resolution
Reconstruction of Cell-Type-Specific Transcriptional Regulatory
Processes from Bulk Sequencing Samples’. *Nature Biotechnology*, ahead
of print. <https://doi.org/10.1038/s41587-026-03218-w>.

Yuki, Shintaro, Chikara Mizukoshi, Ko Abe, and Teppei Shimamura. 2026.
*scTREND: An Annotation-Free Single-Cell Time-Resolved and
Condition-Dependent Hazard Model*. bioRxiv.
<https://doi.org/10.64898/2026.01.26.701686>.

Zaitsev, Aleksandr, Maksim Chelushkin, Daniiar Dyikanov, et al. 2022.
‘Precise Reconstruction of the TME Using Bulk RNA-seq and a Machine
Learning Algorithm Trained on Artificial Transcriptomes’. *Cancer Cell*
40. <https://doi.org/10.1016/j.ccell.2022.07.006>.

Zhang, Ze, Xu Wang, Fan Hong, Pei Yu, Shengbao Suo, and Ye-Guang Chen.
2026. ‘Integrated Inference of Cellular Compositions and Gene Expression
Programs by Deconvolution’. *Cell Regeneration* 15 (1).
<https://doi.org/10.1186/s13619-026-00299-5>.

Zhao, Tianyi, Renjie Liu, Yuzhi Sun, et al. 2026. ‘DECODE: Deep
Learning-Based Common Deconvolution Framework for Various Omics Data’.
*Nature Methods* 23 (3): 596–608.
<https://doi.org/10.1038/s41592-026-03007-y>.

Zou, Zhiyi, Yuting Bai, Bo Wang, Wanwan Shi, Xiao Liang, and Jiawei Luo.
2026. ‘HIDF: Integrating Tree‐structured scRNA‐seq Heterogeneity for
Hierarchical Deconvolution of Spatial Transcriptomics’. *Advanced
Science* 13 (12). <https://doi.org/10.1002/advs.202514073>.
