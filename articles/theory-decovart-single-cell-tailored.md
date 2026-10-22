# DeCovarT tailored to single-cell references: reference moments, scales and four alignment strategies

> **Scope**
>
> Released DeCovarT fits the convolution \boldsymbol{y}\_{\cdot
> i}\mid\boldsymbol{p}\_{\cdot i}
> \sim\mathcal{N}\_G(\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i}, \sum_j
> p\_{ji}^{2}\boldsymbol{\Sigma}\_j) of the [derivatives
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-unconstrained).
> That law is exact when each reference column is one random
> **population profile**, as for the FACS-sorted references it was
> designed for. A single-cell reference supplies a different object:
> many cells per type, counted by a UMI protocol, with a within-type
> cell-to-cell covariance, and a bulk sample that is a physical sum of
> cells measured on another platform.
>
> This vignette rebuilds DeCovarT for that reference. The measurement
> and expression models of UMI counts, and the published methods that
> use second-order moments, are reviewed in the [distribution-aware
> deconvolution
> review](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.md);
> this vignette uses their conclusions.
>
> 1.  [Sec. 1](#sec-reference) builds the mean and the covariance layers
>     of a reference from cells and replicates.
> 2.  [Sec. 2](#sec-one-formula) shows that mixture, convolution and
>     correlated cells are three cases of one covariance identity.
> 3.  [Sec. 3](#sec-general) gives one score, Hessian and Fisher
>     information for any Gaussian whose mean and covariance depend on
>     parameters.
> 4.  [Sec. 4](#sec-strategies) derives four alignment strategies: cell
>     counts, averaged mixture, compositional bulk, and plug-in scales
>     estimated before the deconvolution.
> 5.  [Sec. 5](#sec-uncoupling) converts RNA fractions into cell
>     fractions after the fit.
>
> The final model for real data, written as a generalised linear model,
> is in the [MLE and inference
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-final-choice).
> The simplex charts of the derivatives vignette apply unchanged to the
> proportion block of every formulation. Nothing here is implemented in
> [`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md);
> the equations are the specification for a later release.

> **Note 1: Notation**
>
> Indices: gene g=1,\ldots,G; cell type j=1,\ldots,J (a second type is
> \ell); bulk sample i=1,\ldots,N; independent biological replicate
> (donor) of the single-cell reference r=1,\ldots,R; cell c, with
> C\_{jr} cells of type j in replicate r and C_j in total.
>
> | Symbol | Meaning |
> |----|----|
> | x\_{gc}, \boldsymbol{x}\_{\cdot c} | Observed UMI count of gene g in cell c |
> | \boldsymbol{m}\_{\cdot c} | True molecule content of cell c |
> | s_c=s^{\mathrm{bio}}\_cs^{\mathrm{tech}}\_c | Cell size factor: biological (cell volume, cycle) times technical (capture, depth); \mathbb{E}\[s_c\]=1 |
> | \nu=\operatorname{CV}^{2}(s_c) | Extrinsic noise, 1+\nu=(1+\nu^{\mathrm{bio}})(1+\nu^{\mathrm{tech}}) |
> | b_g | Capture efficiency of gene g in one protocol |
> | \tilde b_g=b_g^{\mathrm{bulk}}/b_g^{\mathrm{ref}}, \tilde{\mathbf{B}}=\operatorname{diag}(\tilde{\boldsymbol{b}}) | Gene-wise platform ratio between the bulk and the reference protocols |
> | \boldsymbol{\lambda}\_{\cdot c}, \boldsymbol{\Sigma}^{\lambda}\_j | Size-free expression of cell c and its covariance inside type j |
> | \boldsymbol{y}\_{\cdot i}\in\mathbb{R}^{G} | Raw bulk profile of sample i |
> | m\_{ji}\in\mathbb{N}, n_i=\sum_j m\_{ji} | Cells of type j in bulk sample i, total cells (`MuSiC` writes m_j^{k}, `RNA-Sieve` n) |
> | \tilde{\boldsymbol{y}}\_{\cdot i}=\boldsymbol{y}\_{\cdot i}/n_i | Bulk profile per cell |
> | \boldsymbol{p}\_{\cdot i}=\boldsymbol{m}\_{\cdot i}/n_i\in\Delta^{J-1} | Cell-number fractions |
> | \boldsymbol{q}\_{\cdot i}, q\_{ji}=p\_{ji}S_j/\sum\_\ell p\_{\ell i}S\_\ell | RNA (transcript) fractions |
> | \tilde{\boldsymbol{\mu}}\_{\cdot j}=C_j^{-1}\sum\_{c\in j}\boldsymbol{x}\_{\cdot c}=S_j\boldsymbol{\theta}\_{\cdot j} | Per-cell mean of type j on raw UMI counts, without library normalisation: cell size times relative abundance, \mathbf{1}^{\top}\boldsymbol{\theta}\_{\cdot j}=1; \tilde{\boldsymbol{\mu}} is the G\times J matrix of these columns |
> | \boldsymbol{\mu}\_{\cdot j} | Mean of one sorted population profile, as in released DeCovarT |
> | S_j=\mathbf{1}^{\top}\tilde{\boldsymbol{\mu}}\_{\cdot j} | Mean transcriptome size (RNA content) of one cell of type j |
> | \boldsymbol{\Sigma}\_{W,j} | Within-replicate cell-to-cell covariance of type j |
> | \boldsymbol{\Sigma}\_{B,j} | Between-replicate covariance of the type-j population mean |
> | \boldsymbol{M}\_{jr}, \boldsymbol{B}\_{jr} | Replicate population mean and its deviation \boldsymbol{M}\_{jr}-\tilde{\boldsymbol{\mu}}\_{\cdot j} |
> | a_i\>0 | Capture ratio between the bulk library and the reference cells (depth, capture), shared by all genes |
> | c_i=a_in_i, \phi_i=a_i^{2}n_i | Mean scale and covariance scale of bulk sample i ([Sec. 4.2.2](#sec-av-scale)) |
> | \boldsymbol{D}\_i=\boldsymbol{\Sigma}\_{\mathrm{obs},i} | Diagonal bulk measurement noise, \operatorname{diag}(\tau\_{1i}^{2},\ldots,\tau\_{Gi}^{2}) |
> | \boldsymbol{T}\_{ji}=\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic} | Total expression contributed by type j to sample i |
> | \boldsymbol{S}\_W(\boldsymbol{v})=\sum_jv_j\boldsymbol{\Sigma}\_{W,j}, \boldsymbol{S}\_B(\boldsymbol{v})=\sum_jv_j^{2}\boldsymbol{\Sigma}\_{B,j} | Linear and quadratic covariance mixtures |
> | \varrho_j | Cell-to-cell correlation inside type j; \boldsymbol{\rho} stays the ILR coordinate |
> | \boldsymbol{\vartheta}, \boldsymbol{\eta} | Generic parameter vector; nuisance scales of [Sec. 4.4](#sec-strategy-plugin) |
> | \boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k,\boldsymbol{\mu}\_{kl},\boldsymbol{\Sigma}\_{kl} | First and second partial derivatives of the mean and covariance in \vartheta_k,\vartheta_l |
>
> A tilde marks a quantity rescaled to a reference unit: per cell for
> expression (\tilde{\boldsymbol{y}}\_{\cdot i},
> \tilde{\boldsymbol{\mu}}\_{\cdot j}), or relative to the reference
> protocol for capture (\tilde{\boldsymbol{b}}). The capture efficiency
> b_g of one protocol multiplies every cell measured with it, so it
> cancels from any contrast inside that protocol; between a single-cell
> reference and a bulk library only the ratio \tilde b_g survives. That
> ratio is the gene-wise platform effect: `RCTD` estimates it on a
> pseudo-bulk of all spots as e^{\hat\gamma_g}, with
> \hat\gamma_g=\log\bar Y_g-\log\sum_j\mu\_{gj}\hat W_j for the pooled
> count \bar Y_g and average type weights \hat W_j, then renormalises
> the reference means to \bar\mu\_{gj}=\mu\_{gj}e^{\hat\gamma_g} ([Cable
> et al. 2022](#ref-cableRobustDecompositionCell2022a)). The ratio
> \lambda_g of `MEAD` and the gene-wise slope of `Bisque` are the same
> nuisance.
>
> When one bulk sample is fixed, the index i is dropped inside
> derivations and restored in the model statements. The single-cell
> papers use j for the gene (Sarkar and Stephens; Fang and Pachter;
> `DECALS`; `Bisque`) or for the subject (`MuSiC`); [the notation map of
> the
> review](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#tbl-notation-map)
> translates every reviewed method.
>
> All plate diagrams follow one symbol key: white squares are fixed
> quantities, blue circles latent random variables, grey circles
> observed variables, dotted circles deterministic functions of their
> parents. Text beside a node is an equation, not a variable. When a
> paper builds increasingly complex models, its panels run from left to
> right.

![](theory-decovart-single-cell-tailored_files/mediabag/fig-four-strategies.svg)

Figure 1: Four strategies to align a single-cell reference with a bulk
sample. Strategies 1–3 differ by the scale on which the bulk is observed
and estimate the scales jointly with the composition; Strategy 4
estimates the technical and biological scales first and plugs them into
the regression or likelihood. Italic lines name the reviewed methods
that already sit on each branch.

## From cells to reference moments

The signature that enters [the released
convolution](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-gaussian-convolution)
is a **population profile**: one mean vector and one covariance per cell
type. Single-cell counts are not that object. This section states how to
build the mean and how the covariance splits, then how to place the
resulting profiles on the same scale as bulk sample i.

![Schematic of three covariance layers: linear within-cell variability,
quadratic between-sample variability, and diagonal bulk measurement
noise.](images/covariance_three_variability_factors.png)

Figure 2: Within-sample cell noise, between-replicate biological
covariance, and residual bulk noise.

![MuSiC model assumptions A1 same population, A2 conserved cell-size
ratios, A3 conserved library-size ratios, and A3-prime for absolute
cell-type proportions.](images/music_model_assumptions.png)

Figure 3: MuSiC assumptions for library-size normalisation: same
population, conserved cell-size ratios, conserved library-size ratios,
and the stronger A3-prime hypothesis for absolute cell-type fractions.

[Figure 2](#fig-three-variability) is the target decomposition. Panel A
is a purified population profile (quadratic in p\_{ji}). Panel B is one
single-cell reference (linear in p\_{ji}). Panel C is the hierarchical
model once R independent replicates exist.

### Mean: average cells, then average replicates

For a fixed type j, the per-cell mean does not depend on how many cells
were captured. With raw counts on an additive scale,

\bar{\boldsymbol{X}}\_{jr} =
\frac{1}{C\_{jr}}\sum\_{c=1}^{C\_{jr}}\boldsymbol{X}\_{jrc}, \qquad
\hat{\tilde{\boldsymbol{\mu}}}\_{\cdot j} =
\frac{1}{R}\sum\_{r=1}^{R}\bar{\boldsymbol{X}}\_{jr}. \tag{1}

Divide by the **cell count** C\_{jr}, not by the library size.
Library-size weights w_c=L_c/\sum_d L_d are functions of
\boldsymbol{X}\_c itself. Their probability limit is
\tilde{\boldsymbol{\mu}}\_{\cdot
j}+\boldsymbol{\Sigma}\_{W,j}\mathbf{1}/S_j, so they target an
RNA-weighted cell and destroy the exact Gaussian scale that leaves the
precision graph unchanged. Under i.i.d. cells the arithmetic mean also
has the largest effective sample size among deterministic convex
weights.

Pool cells across replicates only after the within-replicate mean. A
replicate with 10\\000 cells must not outweigh a replicate with 500
cells. `MuSiC` is built on that multi-subject split: cross-subject means
and variances, with genes weighted by cross-subject and cross-cell
consistency ([Wang et al. 2019](#ref-wangBulkTissueCell2019)).

When every cell is scaled to a common total, the product
S_j\boldsymbol{\theta}\_{\cdot j} collapses to the arithmetic mean of
the raw counts. `MuSiC` writes the signature as cell size times relative
abundance and then recovers that same per-cell mean. Keep S_j and
\boldsymbol{\theta}\_{\cdot j} as separate factors
([Sec. 1.3](#sec-alignment)). Do not library-weight the raw vectors and
call the result S_j.

> **Warning 2: Summing cells is the wrong signature, and the right
> pseudobulk**
>
> Summing is the correct total for a physical pool of exactly C\_{jr}
> cells, \boldsymbol{T}\_{jr}=\sum_c\boldsymbol{X}\_{jrc}, with mean
> C\_{jr}\tilde{\boldsymbol{\mu}}\_{\cdot j} and covariance
> C\_{jr}\boldsymbol{\Sigma}\_{W,j}. Using \boldsymbol{T}\_{jr} as a
> DeCovarT column inserts the arbitrary capture count into
> \boldsymbol{p}. A type with ten times more reference cells looks ten
> times more abundant.
>
> The same sum is the right input for a different question. Differential
> expression of one cell type across donors needs one count vector per
> donor, then a bulk count model such as `DESeq2` ([Love et al.
> 2014](#ref-loveModeratedEstimationFold2014a)). There the estimand is a
> relative fold change of averaged counts, and the donor is the
> independent unit. Cells from one donor are subsamples. Treating them
> as replicates inflates false discoveries ([Squair et al.
> 2021](#ref-squairConfrontingFalseDiscoveries2021)). That analysis does
> not estimate \boldsymbol{p}\_{\cdot i}.

### Three reference regimes

Write the cell as a sum of a replicate effect and a cell effect,

\boldsymbol{X}\_{jrc} = \tilde{\boldsymbol{\mu}}\_{\cdot
j}+\boldsymbol{B}\_{jr}+\boldsymbol{\varepsilon}\_{jrc}, \tag{2}

with
\boldsymbol{B}\_{jr}\sim\mathcal{N}(\boldsymbol{0},\boldsymbol{\Sigma}\_{B,j})
and
\boldsymbol{\varepsilon}\_{jrc}\sim\mathcal{N}(\boldsymbol{0},\boldsymbol{\Sigma}\_{W,j}),
independent across j given \boldsymbol{p}\_{\cdot i}. The bulk profile
per cell of a tissue of n_i cells is

\tilde{\boldsymbol{y}}\_{\cdot i} = \sum_j
p\_{ji}(\tilde{\boldsymbol{\mu}}\_{\cdot j}+\boldsymbol{B}\_{jr(i)}) +
\frac{1}{n_i}\sum_j\sum\_{c=1}^{m\_{ji}}\boldsymbol{\varepsilon}\_{jic}.
\tag{3}

#### FACS-sorted population profiles

Historically the reference was one bulk profile per purified population:
surface-marker sorting, then a pooled microarray or RNA-seq library. The
random object is the whole population vector
\boldsymbol{Z}\_{ij}\sim\mathcal{N}(\boldsymbol{\mu}\_{\cdot
j},\boldsymbol{\Sigma}\_{B,j}), and the bulk is the convolution

\boldsymbol{y}\_{\cdot i} = \sum_j p\_{ji}\boldsymbol{Z}\_{ij}, \qquad
\operatorname{Cov}(\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot
i}) = \sum_j p\_{ji}^{2}\\\boldsymbol{\Sigma}\_{B,j}. \tag{4}

The square appears because one random profile is multiplied by p\_{ji},
not because p\_{ji} cells were added. This is the covariance in [the
released
convolution](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-gaussian-convolution).
If several types are sorted from the same donors, cross-type blocks
\boldsymbol{\Sigma}\_{B,j\ell} add terms p\_{ji}p\_{\ell i}. Before
droplet scRNA-seq, such population profiles came from FACS- or
MACS-purified bulk libraries, already with coexpression among the genes
of a purified type ([Szaniszlo et al.
2004](#ref-szaniszloGettingRightCells); [Novershtern et al.
2011](#ref-novershternDenselyInterconnectedTranscriptional2011)).

#### One single-cell replicate (R=1)

R=1 identifies \boldsymbol{M}\_{j1} and \boldsymbol{\Sigma}\_{W,j}. It
does not identify \boldsymbol{\Sigma}\_{B,j}: extra cells shrink
\boldsymbol{\Sigma}\_{W,j}/C\_{j1} and leave the single draw
\boldsymbol{B}\_{j1} untouched. The coherent bulk covariance is the
**mixture / physical sum**, linear in p\_{ji}:

\operatorname{Cov}(\tilde{\boldsymbol{y}}\_{\cdot
i}\mid\boldsymbol{p}\_{\cdot i}) = \frac{1}{n_i}\sum_j
p\_{ji}\\\boldsymbol{\Sigma}\_{W,j} +
\boldsymbol{\Sigma}\_{\mathrm{obs},i}. \tag{5}

Here n_i is the number of cells behind the bulk library, in practice an
effective number ([Sec. 4.2.2](#sec-av-scale)). The number of reference
cells C\_{j1} is not n_i. A finite-cell reference mean adds a separate
quadratic term \sum_j p\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/C\_{j1}.
That term is Monte Carlo error of \hat{\tilde{\boldsymbol{\mu}}}\_{\cdot
j}, not biological between-donor covariance. Do not drop
\boldsymbol{\Sigma}\_{W,j} into the p\_{ji}^{2} slot of [the released
convolution](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-gaussian-convolution).

For a rare type, p\_{ji}=0.01 makes p\_{ji}^{2} a hundred times smaller
than p\_{ji}. A within-cell covariance placed in the quadratic slot is
almost silenced.

#### Enough biological replicates

The smallest R that identifies a variance component is R\ge 2. The
replicate means satisfy

\bar{\boldsymbol{X}}\_{jr} \sim \mathcal{N}\\\left(
\boldsymbol{\mu}\_{\cdot j},\\
\boldsymbol{\Sigma}\_{B,j}+\frac{\boldsymbol{\Sigma}\_{W,j}}{C\_{jr}}
\right). \tag{6}

A method-of-moments contrast is then available. Let
\boldsymbol{S}\_{B,j} be the sample covariance of
\\\bar{\boldsymbol{X}}\_{jr}\\\_{r=1}^{R} and \bar{C}\_j a typical cell
count. Then

\mathbb{E}\[\boldsymbol{S}\_{B,j}\] = \boldsymbol{\Sigma}\_{B,j} +
\overline{\boldsymbol{\Sigma}\_{W,j}/C\_{jr}}, \tag{7}

so
\hat{\boldsymbol{\Sigma}}\_{B,j}=\boldsymbol{S}\_{B,j}-\widehat{\boldsymbol{\Sigma}}\_{W,j}/\bar{C}\_j,
projected onto the positive semidefinite cone. An unstructured G\times G
covariance has G(G+1)/2 free entries. Transcriptome-scale G makes that
matrix unestimable for the R of a normal single-cell study. What R\ge 2
buys is the **scalar separation** of the two layers, or a sparse
precision for \boldsymbol{\Sigma}\_{B,j}, not a dense G\times G matrix.
Reference cohorts in current single-cell deconvolution are small; `MEAD`
states this explicitly ([Xie et al.
2026](#ref-xieStatisticalInferenceCell2026)). Cells from one donor share
genetic and environmental background, so they are subsamples rather than
independent biological units. Zimmerman, Espeland and Langefeld document
that dependence across cell types, and they treat donor as the random
effect when contrasting conditions ([Zimmerman et al.
2021](#ref-zimmermanPracticalSolutionPseudoreplication2021)). That
within-versus-between split is the same hierarchy as
\boldsymbol{\Sigma}\_{W,j} versus \boldsymbol{\Sigma}\_{B,j} when a
single-cell reference is aligned with bulk.

With that split, the bulk law per cell is

\tilde{\boldsymbol{y}}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i} \sim
\mathcal{N}\\\Big( \sum_j p\_{ji}\tilde{\boldsymbol{\mu}}\_{\cdot j},\\
\sum_j p\_{ji}^{2}\\\boldsymbol{\Sigma}\_{B,j} + \frac{1}{n_i}\sum_j
p\_{ji}\\\boldsymbol{\Sigma}\_{W,j} +
\boldsymbol{\Sigma}\_{\mathrm{obs},i} \Big). \tag{8}

[Sec. 4.2](#sec-averaged) puts it back on the scale of the bulk library.

\boldsymbol{\Sigma}\_{\mathrm{obs},i} is diagonal. Gene–gene correlation
is already inside \boldsymbol{\Sigma}\_{B,j} and
\boldsymbol{\Sigma}\_{W,j}. Heteroscedastic diagonal entries are
allowed. Off-diagonal residual correlation would double-count the
network.

![](theory-decovart-single-cell-tailored_files/mediabag/fig-three-regimes.svg)

Figure 4: Which covariance terms the reference design identifies. FACS
profiles supply only the quadratic layer. One single-cell replicate
supplies only the linear layer. R\ge 2 replicates separate both, plus a
diagonal bulk residual.

### Aligning reference profiles with bulk sample i

[Eq. 1](#eq-mean-two-stage) is a per-cell mean. Bulk sample i is a
library. Three scales sit between them. Only two belong inside the
estimator.

**Sample scales a_i and n_i.** On the raw molecule scale,
\boldsymbol{y}\_{\cdot
i}=a_in_i\tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot i} plus noise:
n_i is the number of cells and a_i the capture ratio between the bulk
library and the reference cells (depth and capture), one positive number
per bulk sample shared by all genes. The mean identifies only
c_i=a_in_i; [Sec. 4.2.2](#sec-av-scale) shows that the covariance
identifies the second scalar. Neither is a composition parameter.

**Cell size S_j.** RNA content differs by cell type,
\tilde{\boldsymbol{\mu}}\_{\cdot j}=S_j\boldsymbol{\theta}\_{\cdot j},
so

\boldsymbol{y}\_{\cdot i} = a_in_i\sum_j p\_{ji}
S_j\boldsymbol{\theta}\_{\cdot j} + \boldsymbol{\varepsilon}\_i. \tag{9}

`MuSiC` states the same product: cell fraction, average cell size,
relative abundance ([Wang et al. 2019](#ref-wangBulkTissueCell2019)). If
S_j is measured, insert it. Do not re-estimate it. Free S_j and
\boldsymbol{p}\_{\cdot i} are not separately identifiable from
\boldsymbol{y}\_{\cdot i} alone, because only the products p\_{ji}S_j
enter [Eq. 9](#eq-music-mean). The RNA fraction

q\_{ji} = \frac{p\_{ji}S_j}{\sum\_\ell p\_{\ell i}S\_\ell} \tag{10}

is what a normalised bulk identifies. \boldsymbol{q}\_{\cdot
i}=\boldsymbol{p}\_{\cdot i} when S_j is constant across j. That
constant-efficiency hypothesis is assumption A3' of `MuSiC`: the ratio
of average library size equals the ratio of average cell size,
equivalently a common capture efficiency \gamma_j=L_j/S_j
([Figure 3](#fig-music-assumptions)). Under A3', the per-cell mean of
raw counts \tilde{\boldsymbol{\mu}}\_{\cdot j} carries the true size
ratios and the fitted weights are cytometry fractions
\boldsymbol{p}\_{\cdot i}. A library-normalised signature
\boldsymbol{\theta}\_{\cdot j} returns RNA fractions
\boldsymbol{q}\_{\cdot i} instead; and if A3' fails, the UMI totals S_j
mix cell size with type-specific capture. Both cases need the
RNA-to-cell correction of [Sec. 5](#sec-uncoupling), which also compares
`EPIC`, `quanTIseq`, `MMAD` and `ReDeconv`.

**Gene-wise platform ratio \tilde{\boldsymbol{b}}.** A platform shift
can be gene-specific, y\_{gi}\approx\tilde b_g a_in_i\sum_j
p\_{ji}\tilde\mu\_{gj}. `Bisque` learns such a map from paired bulk and
single-cell measurements ([Jew et al.
2020](#ref-jewAccurateEstimationCell2020)), `RCTD` from a pseudo-bulk of
all spots ([Note 1](#nte-notation)). Fitting \tilde{\boldsymbol{b}} on
the bulk sample that is being deconvolved lets the scale absorb any
residual, and \boldsymbol{p}\_{\cdot i} is then unidentified. Estimate
\tilde{\boldsymbol{b}} on an external paired cohort or a pool of many
samples, or do not estimate it.

#### Normalisation that preserves the linear mixture

The moment identities above are for **raw counts on an additive scale**.
They do not survive a logarithm: \log\boldsymbol{y}\_{\cdot i} is not
\sum_j p\_{ji}\log\boldsymbol{\mu}\_{\cdot j}.

The comparison DeCovarT needs is the **same gene on two technologies**,
not two genes inside one sample. Gene length is a property of the assay.
It cancels for a within-gene contrast and should not be divided out of
the reference unless the bulk was length-normalised in the same way.
`ReDeconv` treats length as a bulk-only TPM/RPKM correction and leaves
UMI single-cell counts uncorrected for length ([Lu et al.
2025](#ref-luTranscriptomeSizeMatters2025)).

| Procedure | What it estimates | Use for \tilde{\boldsymbol{\mu}}\_{\cdot j} |
|----|----|----|
| Arithmetic mean of raw counts | Per-cell mean | Default signature |
| CPM / CP10K | Library-size relative profile | Erases S_j |
| TPM / RPKM / FPKM | Length-adjusted composition | Breaks cross-platform per-gene alignment |
| TMM (`edgeR`) or median-of-ratios (`DESeq2`) | Effective library size for differential expression | Donor-level DE, not the signature ([Robinson et al. 2010](#ref-robinsonEdgeRBioconductorPackage2010); [Love et al. 2014](#ref-loveModeratedEstimationFold2014a)) |
| `SCTransform` | Regularised negative-binomial residuals | Clustering and integration, not an additive mixture mean ([Hafemeister and Satija 2019](#ref-hafemeisterNormalizationVarianceStabilization2019)) |
| CLTS (`ReDeconv`) | Counts rescaled by a linearised transcriptome size | Keeps S_j for deconvolution ([Lu et al. 2025](#ref-luTranscriptomeSizeMatters2025)) |

CPM and CP10K divide by the cell total. That is a relative profile
\boldsymbol{\theta}\_{\cdot j} with S_j discarded. TMM and the `DESeq2`
size factor answer a different question: an effective depth for
between-sample differential expression, under a hypothesis that most
genes are stable ([Robinson et al.
2010](#ref-robinsonEdgeRBioconductorPackage2010); [Love et al.
2014](#ref-loveModeratedEstimationFold2014a)). `SCTransform` models
depth inside a negative binomial; the residuals are not an additive
convolution mean ([Hafemeister and Satija
2019](#ref-hafemeisterNormalizationVarianceStabilization2019)).

## Convolution, mixture and correlated cells

Let
\boldsymbol{X}\sim\mathcal{N}\_G(\boldsymbol{\mu},\boldsymbol{\Sigma}).
Two statements look contradictory,

\operatorname{Cov}(a\boldsymbol{X})=a^{2}\boldsymbol{\Sigma}
\qquad\text{versus}\qquad
\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)=a\boldsymbol{\Sigma},
\tag{11}

and both are correct. The difference lies entirely in the dependence
between the a summands, not in two rules for Gaussian vectors.

### Characteristic-function proof

The characteristic function of \boldsymbol{X} is
\phi\_{\boldsymbol{X}}(\boldsymbol{t})
=\exp(i\boldsymbol{t}^{\top}\boldsymbol{\mu}
-\tfrac12\boldsymbol{t}^{\top}\boldsymbol{\Sigma}\boldsymbol{t}).

**Scaling** substitutes a\boldsymbol{t} for \boldsymbol{t}:

\phi\_{a\boldsymbol{X}}(\boldsymbol{t})
=\phi\_{\boldsymbol{X}}(a\boldsymbol{t})
=\exp\bigl(i\boldsymbol{t}^{\top}(a\boldsymbol{\mu})
-\tfrac12\boldsymbol{t}^{\top}(a^{2}\boldsymbol{\Sigma})\boldsymbol{t}\bigr),
\qquad
a\boldsymbol{X}\sim\mathcal{N}\_G(a\boldsymbol{\mu},a^{2}\boldsymbol{\Sigma}).
\tag{12}

**Independent convolution** multiplies a identical factors:

\phi\_{\boldsymbol{X}\_1+\cdots+\boldsymbol{X}\_a}(\boldsymbol{t})
=\phi\_{\boldsymbol{X}}(\boldsymbol{t})^{a}
=\exp\bigl(i\boldsymbol{t}^{\top}(a\boldsymbol{\mu})
-\tfrac12\boldsymbol{t}^{\top}(a\boldsymbol{\Sigma})\boldsymbol{t}\bigr),
\qquad
\sum\_{r}\boldsymbol{X}\_r\sim\mathcal{N}\_G(a\boldsymbol{\mu},a\boldsymbol{\Sigma}).
\tag{13}

The first route squares the quadratic exponent; the second adds a copies
of it. The density convolution (f\*f)(\boldsymbol{s})=\int
f(\boldsymbol{x})f(\boldsymbol{s}-\boldsymbol{x})\\d\boldsymbol{x} is
the law of \boldsymbol{X}\_1+\boldsymbol{X}\_2 only when the two draws
are independent. Saying “a convolution of Gaussians” therefore already
asserts independent summands.

### The master identity and its two endpoints

For arbitrary random vectors with
\operatorname{Cov}(\boldsymbol{X}\_r)=\boldsymbol{\Sigma} and
cross-covariances
\boldsymbol{C}\_{rs}=\operatorname{Cov}(\boldsymbol{X}\_r,\boldsymbol{X}\_s),

\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)
=\sum\_{r,s}\operatorname{Cov}(\boldsymbol{X}\_r,\boldsymbol{X}\_s)
=a\boldsymbol{\Sigma}+\sum\_{r\neq s}\boldsymbol{C}\_{rs}. \tag{14}

Under **equicorrelation**,
\boldsymbol{C}\_{rs}=\varrho\boldsymbol{\Sigma} for r\neq s, the a(a-1)
off-diagonal blocks collapse to

\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)
=a\bigl\[1+(a-1)\varrho\bigr\]\boldsymbol{\Sigma}, \qquad
-\frac{1}{a-1}\le\varrho\le 1 . \tag{15}

Setting \varrho=1 (the same realised vector reused a times) gives
a^{2}\boldsymbol{\Sigma}; setting \varrho=0 (fresh independent draws)
gives a\boldsymbol{\Sigma}. [Table 1](#tbl-dependence-regimes) lists the
cases. A sum of merely *marginally* Gaussian vectors need not be
Gaussian; [Eq. 15](#eq-equicorrelation) is a Gaussian law only when the
stacked vector
(\boldsymbol{X}\_1^{\top},\ldots,\boldsymbol{X}\_a^{\top})^{\top} is
**jointly** Gaussian, which is the assumption made whenever an
intermediate \varrho is used below.

| Construction | Dependence | Mean | Covariance | Law |
|----|----|----|----|----|
| Same \boldsymbol{X} repeated | \varrho=1 | a\boldsymbol{\mu} | a^{2}\boldsymbol{\Sigma} | \mathcal{N}\_G |
| a i.i.d. copies | \varrho=0 | a\boldsymbol{\mu} | a\boldsymbol{\Sigma} | \mathcal{N}\_G |
| Correlated copies | general \boldsymbol{C}\_{rs} | a\boldsymbol{\mu} | a\boldsymbol{\Sigma}+\sum\_{r\neq s}\boldsymbol{C}\_{rs} | Gaussian if jointly Gaussian |
| Equicorrelated copies | \boldsymbol{C}\_{rs}=\varrho\boldsymbol{\Sigma} | a\boldsymbol{\mu} | a\[1+(a-1)\varrho\]\boldsymbol{\Sigma} | Gaussian if jointly Gaussian |

Table 1: Dependence regimes for a sum of a Gaussian vectors with common
marginal covariance. The first two rows are the endpoints \varrho=1 and
\varrho=0 of the last.

The biological reading is the one that matters for a single-cell
reference. If \boldsymbol{\Sigma}\_j is the covariance **between cells**
of type j, summing m_j independent cells yields
m_j\boldsymbol{\Sigma}\_j. If \boldsymbol{\Sigma}\_j is the covariance
of **one latent type-level profile** that is multiplied by a weight, the
contribution is p_j^{2}\boldsymbol{\Sigma}\_j. Released DeCovarT is the
second case. The random-effects decomposition of
[Sec. 1.2](#sec-three-regimes),

\boldsymbol{X}\_{jc}=\boldsymbol{M}\_j+\boldsymbol{\varepsilon}\_{jc},
\qquad
\operatorname{Cov}\Bigl(\sum\_{c=1}^{m_j}\boldsymbol{X}\_{jc}\Bigr)
=m_j\boldsymbol{\Sigma}\_{W,j}+m_j^{2}\boldsymbol{\Sigma}\_{B,j},
\tag{16}

contains both endpoints at once: independent cell fluctuations
accumulate linearly, a shared type-level fluctuation is amplified
quadratically. Equicorrelated cells are the special case
\boldsymbol{\Sigma}\_{B,j}=\varrho_j\boldsymbol{\Sigma}\_j and
\boldsymbol{\Sigma}\_{W,j}=(1-\varrho_j)\boldsymbol{\Sigma}\_j
([Sec. 2.4](#sec-persp-equicorrelation)).

### Mixture of cells versus convolution of populations

The identity [Eq. 14](#eq-master-identity) fixes the covariance once the
dependence between summands is fixed. Deconvolution models differ by a
second choice: whether the cell type of each summand is random or given.
[Figure 5](#fig-dag-dependence) lays out the four cases from no
dependence to full dependence.

![](theory-decovart-single-cell-tailored_files/mediabag/fig-dag-dependence.svg)

Figure 5: Four generative stories for a bulk profile \boldsymbol{y},
ordered by the dependence between cells. (a) A mixture: each of n cells
draws its type label from \mathrm{Cat}(\boldsymbol{p}). (b) Independent
cells with fixed counts m_j. (c) Equicorrelated cells sharing a
type-level effect \boldsymbol{B}\_j. (d) One population profile per
type, scaled by p_j. The note under each \boldsymbol{y} is its
covariance; \boldsymbol{\delta}\_j=\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\boldsymbol{p}.

**Mixture of cells.** If the type label z_c of each of n cells is drawn
from \mathrm{Cat}(\boldsymbol{p}), one cell follows the Gaussian mixture
\sum_jp_j\mathcal{N}\_G(\boldsymbol{\mu}\_{\cdot
j},\boldsymbol{\Sigma}\_j), with

\mathbb{E}\[\boldsymbol{X}\_c\]=\boldsymbol{\mu}\boldsymbol{p}, \qquad
\operatorname{Cov}(\boldsymbol{X}\_c) =\sum_jp_j\boldsymbol{\Sigma}\_j
+\sum_jp_j(\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\boldsymbol{p})(\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\boldsymbol{p})^{\top} \tag{17}

by the law of total covariance. The sum of n independent such cells has
n times that covariance and is Gaussian only in the central-limit sense.
This is the `RNA-Sieve` likelihood, written gene by gene
([`RNA-Sieve`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-rnasieve)).
The second term is the spread of the type means; it is present because
the labels are random, and it disappears once the counts \boldsymbol{m}
are conditioned on.

**Convolution of independent cells.** With the counts given, the sum is
an exact Gaussian with covariance \sum_jm_j\boldsymbol{\Sigma}\_j:
linear weights. This is Strategy 1 ([Sec. 4.1](#sec-cell-counts)) and,
after division by n, the linear layer of Strategy 2.

**Convolution of population profiles.** If each type contributes one
random profile scaled by its weight, the covariance is
\sum_jp_j^{2}\boldsymbol{\Sigma}\_j: quadratic weights. This is released
DeCovarT, the `DECALS` covariance, and the cross-subject layer of
`MuSiC`.

The equicorrelated model sits between the last two.

### Correlated cells inside a type

[Sec. 4.1](#sec-cell-counts) assumed independent cells within a type.
Cells of one type in one tissue share micro-environment, cell cycle
phase and paracrine signals, so a positive within-type correlation is
plausible. The jointly Gaussian equicorrelated model of
[Eq. 15](#eq-equicorrelation) with a type-specific \varrho_j gives

\operatorname{Cov}(\boldsymbol{T}\_{ji})
=m\_{ji}\bigl\[1+(m\_{ji}-1)\varrho_j\bigr\]\boldsymbol{\Sigma}\_{W,j}
=m\_{ji}(1-\varrho_j)\boldsymbol{\Sigma}\_{W,j}+m\_{ji}^{2}\varrho_j\boldsymbol{\Sigma}\_{W,j}.
\tag{18}

The second form shows that equicorrelation is the random-effects model
[Eq. 16](#eq-random-effects) with proportional layers,
\boldsymbol{\Sigma}\_{B,j}=\varrho_j\boldsymbol{\Sigma}\_{W,j} and a
within covariance deflated to (1-\varrho_j)\boldsymbol{\Sigma}\_{W,j}.
The covariance derivative in \varrho_j is
\boldsymbol{\Sigma}\_{\varrho_j}=m\_{ji}(m\_{ji}-1)\boldsymbol{\Sigma}\_{W,j},
so [Theorem 1](#thm-score-general) applies. Identifiability is the
obstacle: with one bulk sample, \varrho_j and m\_{ji} enter only through
the scalar m\_{ji}\[1+(m\_{ji}-1)\varrho_j\] in front of
\boldsymbol{\Sigma}\_{W,j} and cannot be separated. Several samples with
different compositions, or an external estimate of \varrho_j from
spatial neighbourhoods of the same type, are needed. The lower bound
\varrho_j\ge-1/(m\_{ji}-1) is the positive-semidefiniteness constraint
and shrinks to zero for abundant types, so negative within-type
correlation is a rare-type phenomenon only. The claim that a sum of
equicorrelated cells is Gaussian rests on **joint** Gaussianity of the
stacked cells, which is an assumption and not a consequence of Gaussian
marginals.

## One score and Hessian for every formulation

All models below are Gaussian with a parametrised mean and covariance.
Write \boldsymbol{y}\mid\boldsymbol{\vartheta}\sim
\mathcal{N}\_G\bigl(\boldsymbol{\mu}(\boldsymbol{\vartheta}),
\boldsymbol{\Sigma}(\boldsymbol{\vartheta})\bigr),
\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}(\boldsymbol{\vartheta}),
\boldsymbol{\Omega}=\boldsymbol{\Sigma}^{-1}, and up to a constant

\ell(\boldsymbol{\vartheta})
=\tfrac12\log\det\boldsymbol{\Omega}(\boldsymbol{\vartheta})
-\tfrac12\boldsymbol{r}^{\top}\boldsymbol{\Omega}(\boldsymbol{\vartheta})\boldsymbol{r}.
\tag{19}

The identities of the derivatives vignette ([matrix
calculus](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-matrix-calculus))
give the following.

**Theorem 1 (Score, Hessian and expected information of a parametrised
Gaussian)** With
\boldsymbol{\mu}\_k=\partial\boldsymbol{\mu}/\partial\vartheta_k,
\boldsymbol{\Sigma}\_k=\partial\boldsymbol{\Sigma}/\partial\vartheta_k,
\boldsymbol{\mu}\_{kl}=\partial^{2}\boldsymbol{\mu}/\partial\vartheta_k\partial\vartheta_l
and
\boldsymbol{\Sigma}\_{kl}=\partial^{2}\boldsymbol{\Sigma}/\partial\vartheta_k\partial\vartheta_l,

\frac{\partial\ell}{\partial\vartheta_k}
=\underbrace{-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k)}\_{\text{determinant}}
+\underbrace{\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{r}}\_{\text{mean
residual}}
+\underbrace{\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}}\_{\text{covariance
quadratic}}, \tag{20}

\begin{aligned}
\frac{\partial^{2}\ell}{\partial\vartheta_k\partial\vartheta_l}
&=\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l)
-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl})
-\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_l
+\boldsymbol{\mu}\_{kl}^{\top}\boldsymbol{\Omega}\boldsymbol{r} \\
&\quad
-\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{\mu}\_l^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl}\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r},
\end{aligned} \tag{21}

I(\boldsymbol{\vartheta})\_{kl}
=\mathbb{E}\Bigl\[-\frac{\partial^{2}\ell}{\partial\vartheta_k\partial\vartheta_l}\Bigr\]
=\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_l
+\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l).
\tag{22}

*Proof*. The determinant term follows from
\partial\log\det\boldsymbol{\Omega}/\partial\vartheta_k
=-\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k) and
\partial\boldsymbol{\Omega}/\partial\vartheta_k
=-\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}.
Differentiating
-\tfrac12\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{r} uses
\partial\boldsymbol{r}/\partial\vartheta_k=-\boldsymbol{\mu}\_k and the
symmetry of \boldsymbol{\Omega}. Differentiating
[Eq. 20](#eq-score-general) once more and using
\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}
=\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}
(transpose of a scalar) gives [Eq. 21](#eq-hess-general). For
[Eq. 22](#eq-fisher-general) substitute
\mathbb{E}\[\boldsymbol{r}\]=\boldsymbol{0} and
\mathbb{E}\[\boldsymbol{r}\boldsymbol{r}^{\top}\]=\boldsymbol{\Sigma}:
the \boldsymbol{\Sigma}\_{kl} traces cancel and
\mathbb{E}\[\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}\]
=\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l).

> **Tip 3: Check against released DeCovarT**
>
> Take \boldsymbol{\vartheta}=\boldsymbol{p},
> \boldsymbol{\mu}(\boldsymbol{p})=\boldsymbol{\mu}\boldsymbol{p} and
> \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
> p_j^{2}\boldsymbol{\Sigma}\_j. Then
> \boldsymbol{\mu}\_{p_j}=\boldsymbol{\mu}\_{\cdot j},
> \boldsymbol{\mu}\_{p_jp\_\ell}=\boldsymbol{0}, the covariance
> derivative is \boldsymbol{\Sigma}\_{p_j}=2p_j\boldsymbol{\Sigma}\_j
> with \boldsymbol{\Sigma}\_j the type covariance,
> \boldsymbol{\Sigma}\_{p_jp_j}=2\boldsymbol{\Sigma}\_j and
> \boldsymbol{\Sigma}\_{p_jp\_\ell}=\boldsymbol{0} for j\neq\ell.
> Substituting into [Eq. 20](#eq-score-general) to
> [Eq. 22](#eq-fisher-general) recovers the [unconstrained
> gradient](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-score-unconstrained),
> [Hessian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-hess-unconstrained)
> and [expected Fisher
> information](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-fisher-p)
> of the derivatives vignette term by term, including the factor
> 2p_jp\_\ell in front of the covariance trace.

Every model below is therefore specified by its tuple
(\boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k,\boldsymbol{\mu}\_{kl},\boldsymbol{\Sigma}\_{kl}).
When a parameter block lives on the simplex, the ILR chain rule
([constrained score and
Hessian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-chain-rule))
is applied to that block only; the Jacobian of the full map is block
diagonal with identity blocks for the scalar parameters (c_i, \phi_i,
\tau\_{gi}^{2}). Positivity of a scalar s\>0 is handled by the log chart
s=e^{t}, which multiplies the corresponding row of the score by s and
adds the score itself to the diagonal of the Hessian.

### Compound parametrisations and nuisance scales

[Theorem 1](#thm-score-general) is the general form of the Fisher
information of a multivariate normal whose mean and covariance depend on
a parameter,
I\_{kl}=\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Sigma}^{-1}\boldsymbol{\mu}\_l
+\tfrac12\operatorname{tr}(\boldsymbol{\Sigma}^{-1}\boldsymbol{\Sigma}\_k\boldsymbol{\Sigma}^{-1}\boldsymbol{\Sigma}\_l)
([Fisher information of the multivariate
normal](https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution)).
Every model below is **compound**: \boldsymbol{\vartheta} is itself a
function \boldsymbol{\vartheta}=\boldsymbol{h}(\boldsymbol{t}) of the
free coordinates \boldsymbol{t} that are optimised (ILR coordinates, log
scales, counts and sizes). With the Jacobian
\mathbf{J}=\partial\boldsymbol{\vartheta}/\partial\boldsymbol{t}^{\top},

\nabla\_{\boldsymbol{t}}\ell=\mathbf{J}^{\top}\nabla\_{\boldsymbol{\vartheta}}\ell,
\qquad \frac{\partial^{2}\ell}{\partial t_a\partial t_b}
=\bigl(\mathbf{J}^{\top}\mathbf{H}\_{\boldsymbol{\vartheta}}\mathbf{J}\bigr)\_{ab}
+\sum_k\frac{\partial\ell}{\partial\vartheta_k}\frac{\partial^{2}\vartheta_k}{\partial
t_a\partial t_b}, \qquad
I(\boldsymbol{t})=\mathbf{J}^{\top}I(\boldsymbol{\vartheta})\mathbf{J}.
\tag{23}

The second Hessian term is the score contraction of the map Hessian; it
vanishes in expectation, which is why the information transforms without
it.

Most strategies below carry nuisance scales
\boldsymbol{\eta}=(c_i,\phi_i,\boldsymbol{\tau}^{2},\ldots) next to the
composition. Partition the information into a \boldsymbol{p} block and
an \boldsymbol{\eta} block, with cross block
I\_{p\eta,jl}=\boldsymbol{\mu}\_{p_j}^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_{\eta_l}
+\tfrac12\operatorname{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{p_j}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{\eta_l})
from [Eq. 22](#eq-fisher-general). When \boldsymbol{\eta} is estimated
jointly, the information left for \boldsymbol{p} is the Schur complement

I\_{p\mid\eta}=I\_{pp}-I\_{p\eta}I\_{\eta\eta}^{-1}I\_{\eta p}, \qquad
\operatorname{Var}(\hat{\boldsymbol{p}})\approx
I\_{p\mid\eta}^{-1}\succeq I\_{pp}^{-1}. \tag{24}

The loss I\_{pp}-I\_{p\mid\eta} is zero when the cross block vanishes,
for example when a scale enters the covariance only and its trace term
is orthogonal to those of \boldsymbol{p}.
[Sec. 4.4](#sec-strategy-plugin) uses [Eq. 24](#eq-schur-information) to
compare joint and plug-in estimation of the scales.

## Four strategies to align a single-cell reference with bulk

[Figure 1](#fig-four-strategies) is the map of the strategies.
Strategies 1–3 differ by the scale on which the bulk is observed, and
estimate every scale jointly with the composition. Strategy 4 keeps one
of those likelihoods but estimates the technical and biological scales
first, from other data, and plugs them in. The reviewed methods sit on
these branches as follows: `RNA-Sieve` and `DWLS` estimate counts or a
cell number (Strategy 1); `DECALS`, `MEAD` and DeCovarT estimate a
composition with a nuisance scale (Strategy 2); `MuSiC`, `Bisque` and
`DECALS` work on library-normalised or FPKM bulk (Strategy 3); and
`MuSiC`’s S_j, `RCTD`’s platform, `Bisque`’s map and `MEAD`’s reference
scale are plug-ins (Strategy 4).

### Strategy 1: cell-count convolution

The first formulation keeps the bulk on an absolute molecule scale and
asks how many cells of each type it contains. It is the cell-level
reading of `MuSiC`’s equation (2), which in our indices reads
y\_{gi}=\sum_j m\_{ji}S\_{j}\theta\_{gj} with subject-specific sizes and
abundances ([Wang et al. 2019](#ref-wangBulkTissueCell2019)). With
\tilde{\boldsymbol{\mu}}\_{\cdot j}=S_j\boldsymbol{\theta}\_{\cdot j}
the mean below is that identity before `MuSiC` divides by the library
total. The estimand is \boldsymbol{m}\_{\cdot i}\in\mathbb{N}^{J},
proportions follow by division, and no simplex constraint is imposed.

#### Layer 1: independent cells, linear covariance

Assume one single-cell replicate (R=1) and, conditional on the counts,
independent cells

\boldsymbol{X}\_{jic}\overset{\text{i.i.d.}}{\sim}
\mathcal{N}\_G(\tilde{\boldsymbol{\mu}}\_{\cdot
j},\boldsymbol{\Sigma}\_{W,j}), \qquad c=1,\ldots,m\_{ji}, \tag{25}

independent across types. The bulk is the physical sum

\boldsymbol{y}\_{\cdot i} =\sum\_{j=1}^{J}\boldsymbol{T}\_{ji}
=\sum\_{j=1}^{J}\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}. \tag{26}

Applying [Eq. 14](#eq-master-identity) with
\boldsymbol{C}\_{rs}=\boldsymbol{0} inside each type and independence
across types,

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i}
\sim\mathcal{N}\_G\Bigl( \tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot
i},\\ \boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i}) \Bigr), \qquad
\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i})=\sum_j
m\_{ji}\boldsymbol{\Sigma}\_{W,j}. \tag{27}

Each term has one meaning.
\tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot i} is the expected
molecule count of a tissue made of m\_{ji} average cells of each type.
m\_{ji}\boldsymbol{\Sigma}\_{W,j} is the accumulated cell-to-cell
scatter of type j: independent deviations add, so the weight is linear.
There is no \boldsymbol{\Sigma}\_{B,j} because one replicate cannot
identify it. The law is exact under Gaussian cells. Under non-Gaussian
cells it is the central-limit approximation of `RNA-Sieve`, applied here
to the full vector rather than gene by gene ([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)).

Treat \boldsymbol{m}\_{\cdot i} as a point of (0,\infty)^{J} for the
moment. The tuple of [Theorem 1](#thm-score-general) is

\boldsymbol{\mu}\_{m_j}=\tilde{\boldsymbol{\mu}}\_{\cdot j}, \qquad
\boldsymbol{\Sigma}\_{m_j}=\boldsymbol{\Sigma}\_{W,j}, \qquad
\boldsymbol{\mu}\_{m_jm\_\ell}=\boldsymbol{0}, \qquad
\boldsymbol{\Sigma}\_{m_jm\_\ell}=\boldsymbol{0}. \tag{28}

**Theorem 2 (Score, Hessian and information for cell counts)** With
\boldsymbol{\Omega}=\boldsymbol{S}\_W(\boldsymbol{m})^{-1} and
\boldsymbol{r}=\boldsymbol{y}-\tilde{\boldsymbol{\mu}}\boldsymbol{m},

\frac{\partial\ell}{\partial m_j}
=-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j})
+\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{\Omega}\boldsymbol{r}
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j}\boldsymbol{\Omega}\boldsymbol{r},
\tag{29}

\frac{\partial^{2}\ell}{\partial m_j\partial m\_\ell}
=\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,\ell})
-\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{\Omega}\tilde{\boldsymbol{\mu}}\_{\cdot \ell}
-\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,\ell}\boldsymbol{\Omega}\boldsymbol{r}
-\tilde{\boldsymbol{\mu}}\_{\cdot
\ell}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j}\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,\ell}\boldsymbol{\Omega}\boldsymbol{r},
\tag{30}

I(\boldsymbol{m})\_{j\ell} =\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{\Omega}\tilde{\boldsymbol{\mu}}\_{\cdot \ell}
+\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,j}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,\ell}).
\tag{31}

Compared with released DeCovarT the Hessian loses two terms, because
\boldsymbol{\Sigma}(\boldsymbol{m}) is **linear** in the parameter:
there is no \mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl})
and no second residual quadratic. The covariance derivative is a
constant matrix, so one Cholesky factorisation of
\boldsymbol{S}\_W(\boldsymbol{m}) per iteration serves objective, score
and Hessian, as in
[`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-numerical-speedups).

Proportions and their uncertainty follow by the delta method. With \hat
n=\mathbf{1}^{\top}\hat{\boldsymbol{m}} and
\hat{\boldsymbol{p}}=\hat{\boldsymbol{m}}/\hat n,

\frac{\partial\boldsymbol{p}}{\partial\boldsymbol{m}^{\top}}
=\frac{1}{n}\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr),
\qquad \operatorname{Var}(\hat{\boldsymbol{p}}) \approx\frac{1}{n^{2}}
\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr)
I(\boldsymbol{m})^{-1}
\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr)^{\top}.
\tag{32}

That covariance is singular along \mathbf{1}, as it must be for a
composition, without any log-ratio chart. A type with \hat m\_{ji}=0 is
an ordinary point of the parameter space, not a boundary of a logarithm;
the structural-zero machinery of the [derivatives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-structural-zeros)
is not needed here.

> **Important 4: Where the GRN precision goes**
>
> The reproducibility book estimates one sparse precision per type from
> single cells and exports \\\hat{\tilde{\boldsymbol{\mu}}}\_{\cdot
> j},\hat{\boldsymbol{\Sigma}}\_j,\hat{\boldsymbol{\Omega}}\_j\\ on the
> linear count scale ([GGM
> chapter](https://github.com/bastienchassagnol/DeCovarT_reproducibility/blob/main/docs/02-grn-graphical-lasso.qmd)).
> The rows of that Gaussian table are **cells** from one experiment, so
> the estimate is (\tilde{\boldsymbol{\mu}}\_{\cdot
> j},\boldsymbol{\Sigma}\_{W,j},\boldsymbol{\Omega}\_{W,j}), not a
> population-profile covariance. Its slot is the linear term
> m\_{ji}\boldsymbol{\Sigma}\_{W,j} of [Eq. 27](#eq-cc-law), and the
> per-cell average of C_j such cells has covariance
> \boldsymbol{\Sigma}\_{W,j}/C_j and precision
> C_j\boldsymbol{\Omega}\_{W,j}: the zero pattern survives averaging.
> The precision of the **sum** \boldsymbol{S}\_W(\boldsymbol{m}) is
> dense in general even when every \boldsymbol{\Omega}\_{W,j} is sparse.
> The structure-aware backends of the derivatives vignette act on
> \boldsymbol{\Sigma}, not on \boldsymbol{\Omega}, so a dense Cholesky
> is the default here unless the \boldsymbol{\Sigma}\_{W,j} share a
> block or low-rank structure.

#### Layer 2: a sample-wide capture scale

Bulk and reference are measured on different platforms. Let a_i\>0 be
the sample-wide ratio of molecule capture between the bulk library and
the reference cells, shared by all genes. The whole realised signal is
amplified,

\boldsymbol{y}\_{\cdot i}
=a_i\sum_j\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}, \qquad
\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl(
a_i\tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot i},\\
a_i^{2}\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i}) \Bigr). \tag{33}

The two scalings of [Sec. 2](#sec-one-formula) appear side by side:
m\_{ji} enters linearly because m\_{ji} independent cells are convolved;
a_i enters squared because one realised fluctuation is amplified a_i
times. Reparametrise to effective counts \boldsymbol{w}\_{\cdot
i}=a_i\boldsymbol{m}\_{\cdot i}:

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{w}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl( \tilde{\boldsymbol{\mu}}\boldsymbol{w}\_{\cdot
i},\\ a_i\boldsymbol{S}\_W(\boldsymbol{w}\_{\cdot i}) \Bigr). \tag{34}

The mean identifies \boldsymbol{w}\_{\cdot i} alone; the covariance
magnitude then identifies a_i, and \boldsymbol{m}\_{\cdot
i}=\boldsymbol{w}\_{\cdot i}/a_i. That is the route by which `RNA-Sieve`
infers its total cell number n: the mean grows like an and the variance
like a^{2}n, so their ratio separates the capture ratio from the count
([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)). Its
supplement states the caveat that applies here as well: n is physically
meaningful within one protocol, or when relative amplification factors
are known, and loses that meaning when the cross-protocol scale is
unknown. In [Eq. 34](#eq-cc-w-law) the same caveat reads: a_i is
identified only if \boldsymbol{\Sigma}\_{W,j} is on the molecule scale
of the reference cells, which is where the GRN estimates of
[Note 4](#nte-sparse-precision-slot) live.

In the (\boldsymbol{w},a) chart the tuple is

\boldsymbol{\mu}\_{w_j}=\tilde{\boldsymbol{\mu}}\_{\cdot j},\quad
\boldsymbol{\mu}\_{a}=\boldsymbol{0},\quad
\boldsymbol{\Sigma}\_{w_j}=a\boldsymbol{\Sigma}\_{W,j},\quad
\boldsymbol{\Sigma}\_{a}=\boldsymbol{S}\_W(\boldsymbol{w}),\quad
\boldsymbol{\Sigma}\_{w_ja}=\boldsymbol{\Sigma}\_{W,j},\quad
\boldsymbol{\Sigma}\_{w_jw\_\ell}=\boldsymbol{\Sigma}\_{aa}=\boldsymbol{0},
\tag{35}

and all second mean derivatives vanish. The scale score is

\frac{\partial\ell}{\partial a}
=-\tfrac12\\\mathrm{tr}\bigl(\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})\bigr)
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})\boldsymbol{\Omega}\boldsymbol{r}
=\frac{1}{2a}\bigl(\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{r}-G\bigr),
\tag{36}

because
\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})=a^{-1}\mathbf{I}\_G
when \boldsymbol{D}\_i=\boldsymbol{0}. Setting it to zero gives the
closed form \hat
a=\boldsymbol{r}^{\top}\boldsymbol{S}\_W(\boldsymbol{w})^{-1}\boldsymbol{r}/G
at fixed \boldsymbol{w}: the scale is the mean squared Mahalanobis
residual, and it can be profiled out. The Fisher block for a is
I\_{aa}=G/(2a^{2}), the usual variance-parameter information, so the
relative standard error of \hat a is \sqrt{2/G}. With G in the hundreds
that is a few percent, which is the “fairly weak dependence” `RNA-Sieve`
reports for n.

#### Layer 3: diagonal bulk measurement noise

Add technical bulk noise that is independent across genes,

\boldsymbol{y}\_{\cdot i}
=a_i\sum_j\boldsymbol{T}\_{ji}+\boldsymbol{\epsilon}\_i, \qquad
\boldsymbol{\epsilon}\_i\sim\mathcal{N}\_G(\boldsymbol{0},\boldsymbol{D}\_i),
\qquad \boldsymbol{\Sigma}(\boldsymbol{w},a,\boldsymbol{\tau}^{2})
=a\boldsymbol{S}\_W(\boldsymbol{w})+\boldsymbol{D}\_i . \tag{37}

\boldsymbol{D}\_i is diagonal by construction: gene to gene correlation
is already carried by the \boldsymbol{\Sigma}\_{W,j}, and an
off-diagonal residual would count the network twice
([Eq. 8](#eq-three-layer)). The new derivatives are
\boldsymbol{\Sigma}\_{\tau_g^{2}}=\boldsymbol{e}\_g\boldsymbol{e}\_g^{\top}
and all second derivatives involving \boldsymbol{\tau}^{2} vanish, so

\frac{\partial\ell}{\partial\tau_g^{2}}
=-\tfrac12\\\Omega\_{gg}+\tfrac12\\(\boldsymbol{\Omega}\boldsymbol{r})\_g^{2}.
\tag{38}

The closed form of [Eq. 36](#eq-cc-score-a) no longer holds, because
\boldsymbol{\Omega}\boldsymbol{S}\_W is not a multiple of the identity.
Identifiability now rests on structure rather than on magnitude: a_i
inflates every entry of a_i\boldsymbol{S}\_W(\boldsymbol{w}),
\boldsymbol{D}\_i inflates the diagonal only. The **off-diagonal**
entries of the bulk covariance, which are the gene network, identify
a_i; the diagonal then identifies \boldsymbol{D}\_i. A gene-wise
deconvolution has no off-diagonal entries and cannot separate the two;
`MuSiC` and `RNA-Sieve` accordingly keep a single gene-wise variance.
With one bulk column and G free \tau\_{gi}^{2} the model is
over-parametrised; constrain \boldsymbol{D}\_i=\tau_i^{2}\mathbf{I}\_G,
or
\tau\_{gi}^{2}=\tau_i^{2}\\(\tilde{\boldsymbol{\mu}}\boldsymbol{w})\_g^{\gamma}
with a mean–variance exponent, or estimate \boldsymbol{D} from technical
bulk replicates.

#### Layer 4: a nested hierarchical model when several single-cell replicates are available

With R\ge 2 biological replicates in the reference, write
\boldsymbol{X}\_{jrc}=\tilde{\boldsymbol{\mu}}\_{\cdot
j}+\boldsymbol{B}\_{jr}+\boldsymbol{\varepsilon}\_{jrc} as in
[Eq. 2](#eq-cell-hierarchy), with
\boldsymbol{B}\_{jr}\sim\mathcal{N}(\boldsymbol{0},\boldsymbol{\Sigma}\_{B,j})
shared by all cells of type j in replicate r. A new bulk donor i draws
its own \boldsymbol{B}\_{ji}, so

\boldsymbol{T}\_{ji} =m\_{ji}(\tilde{\boldsymbol{\mu}}\_{\cdot
j}+\boldsymbol{B}\_{ji})
+\sum\_{c=1}^{m\_{ji}}\boldsymbol{\varepsilon}\_{jic}, \qquad
\operatorname{Cov}(\boldsymbol{T}\_{ji})
=m\_{ji}\boldsymbol{\Sigma}\_{W,j}+m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j},
\tag{39}

and the full five-layer law is

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl(
a_i\tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot i},\\
a_i^{2}\bigl\[\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
i})+\boldsymbol{S}\_B(\boldsymbol{m}\_{\cdot i})\bigr\]
+\boldsymbol{D}\_i \Bigr). \tag{40}

The quadratic term is DeCovarT’s original covariance written in cell
counts. In the \boldsymbol{m} chart the new derivatives are

\boldsymbol{\Sigma}\_{m_j}
=a^{2}\bigl(\boldsymbol{\Sigma}\_{W,j}+2m_j\boldsymbol{\Sigma}\_{B,j}\bigr),
\qquad \boldsymbol{\Sigma}\_{m_jm_j}=2a^{2}\boldsymbol{\Sigma}\_{B,j},
\qquad \boldsymbol{\Sigma}\_{m_jm\_\ell}=\boldsymbol{0}\\ (j\neq\ell),
\tag{41}

so the two Hessian terms that vanished in [Sec. 4.1.1](#sec-cc-linear)
return, carried by \boldsymbol{\Sigma}\_{B,j} alone.
\boldsymbol{\Sigma}\_{B,j} is estimated by the method-of-moments
contrast of [Sec. 1.2.3](#sec-enough-replicates). A dense G\times G
between-replicate covariance needs many more replicates than a
single-cell study provides; with R around ten, independent donors
justify a **diagonal** \boldsymbol{\Sigma}\_{B,j} whose G marginal
variances are estimable, and the gene network stays inside
\boldsymbol{\Sigma}\_{W,j}.

#### Strategy 1 and the reviewed methods

Three reviewed methods sit on this branch; their models are in the
review
([`MuSiC`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-music),
[`RNA-Sieve`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-rnasieve),
[`DWLS`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-dwls)).
Only the contrast with [Eq. 40](#eq-cc-five-layer) is kept here.

- **`MuSiC`** has the mean of [Eq. 27](#eq-cc-law) before division by
  the library total, so it never estimates counts. Its variance carries
  quadratic weights p\_{ji}^{2} on a cross-subject variance: the
  diagonal of the \boldsymbol{S}\_B layer, plus a residual in the role
  of \boldsymbol{D}\_i. The cell scatter \boldsymbol{\Sigma}\_{W,j} only
  selects genes ([Wang et al. 2019](#ref-wangBulkTissueCell2019)).
- **`RNA-Sieve`** has a variance linear in n, as here, but for a
  randomly labelled cell: it adds the spread of the type means
  ([Eq. 17](#eq-mixture-cell)), which [Eq. 27](#eq-cc-law) drops by
  conditioning on the counts. Genes are independent, and the noisy
  reference means add the errors-in-variables term listed in
  [Sec. 6](#sec-summary).
- **`DWLS`** weights genes by
  1/(\tilde{\boldsymbol{\mu}}\hat{\boldsymbol{p}})\_g^{2} plus a
  dampening constant: [Eq. 37](#eq-cc-noise-law) with
  \boldsymbol{\Sigma}\_{W,j}=\boldsymbol{0} and
  \boldsymbol{D}\_i\propto\operatorname{diag}(\tilde{\boldsymbol{\mu}}\boldsymbol{p})^{2},
  a multiplicative noise without a count scale, whose weights are cell
  fractions only if every cell carries the same RNA ([Tsoucas et al.
  2019](#ref-tsoucasAccurateEstimationCelltype2019)).

### Strategy 2: averaged mixture on the simplex

The second formulation divides both sides by the number of cells and
targets the composition \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1}. The
reference is already a per-cell average ([mean
signature](#sec-mean-signature)). The bulk is not: a library is a sum,
and the number of cells behind it is unknown and usually far from the
number of reference cells. Each layer below adds one way of closing that
gap.

#### Layer 1: the per-cell average and the missing n_i

Divide [Eq. 26](#eq-cc-sum) by n_i:

\tilde{\boldsymbol{y}}\_{\cdot i}
=\frac{1}{n_i}\sum_j\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}, \qquad
\tilde{\boldsymbol{y}}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},n_i
\sim\mathcal{N}\_G\Bigl( \tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot
i},\\ \frac{1}{n_i}\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i}) \Bigr).
\tag{42}

The mean is the DeCovarT mean. The covariance is the linear mixture of
[Sec. 1.2.2](#sec-one-replicate): averaging n_i independent cells
shrinks the cell-level covariance by 1/n_i. The quadratic weights also
appear algebraically,
\operatorname{Cov}(p\_{ji}\bar{\boldsymbol{X}}\_{ji})
=p\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/m\_{ji}=p\_{ji}\boldsymbol{\Sigma}\_{W,j}/n_i,
but the covariance of the type-specific average depends on its own cell
count, which turns the square into a linear weight.

n_i is unobserved. Three ways to supply it:

- **Variance route (`RNA-Sieve`).** Estimate n_i jointly with
  \boldsymbol{p}\_{\cdot i} from the covariance magnitude, as in
  1.  This needs \boldsymbol{\Sigma}\_{W,j} on the molecule scale of the
      bulk.
- **Pairing route (`Bisque`).** A gene-wise map from observed bulk to a
  pseudo-bulk of paired individuals, whose compositions are counted in
  their single-cell data, is learned once ([Jew et al.
  2020](#ref-jewAccurateEstimationCell2020)). Its multivariate Gaussian
  version takes the observed bulk of a paired sample i' as an affine
  image of [Eq. 27](#eq-cc-law), \boldsymbol{y}\_{\cdot
  i'}\sim\mathcal{N}\_G\bigl(\boldsymbol{o}+\tilde{\mathbf{B}}\tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot
  i'}, \tilde{\mathbf{B}}\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
  i'})\tilde{\mathbf{B}}+\boldsymbol{D}\bigr), with an intercept
  \boldsymbol{o}; (\boldsymbol{o},\tilde{\mathbf{B}},\boldsymbol{D}) are
  estimated on the paired cohort, and n_i for an unpaired sample follows
  from the fitted mean scale ([Sec. 4.4.4](#sec-persp-affine)).
- **External counts.** Nuclei counts, DNA content or tissue mass give
  n_i directly, as a known offset.

With \boldsymbol{\vartheta}=(\boldsymbol{p},\log n) the tuple is

\boldsymbol{\mu}\_{p_j}=\tilde{\boldsymbol{\mu}}\_{\cdot j},\quad
\boldsymbol{\Sigma}\_{p_j}=n^{-1}\boldsymbol{\Sigma}\_{W,j},\quad
\boldsymbol{\Sigma}\_{\log
n}=-n^{-1}\boldsymbol{S}\_W(\boldsymbol{p}),\quad
\boldsymbol{\Sigma}\_{p_j\log n}=-n^{-1}\boldsymbol{\Sigma}\_{W,j},\quad
\boldsymbol{\Sigma}\_{\log n\log
n}=n^{-1}\boldsymbol{S}\_W(\boldsymbol{p}), \tag{43}

with \boldsymbol{\Sigma}\_{p_jp\_\ell}=\boldsymbol{0} and all second
mean derivatives zero. The simplex is imposed on the \boldsymbol{p}
block by the ILR chart,
\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{\rho}), with
Jacobian \mathbf{J}\_\psi=\mathbf{S}(\boldsymbol{p})\mathbf{V} and the
score contraction of the map Hessian exactly as in the [chain
rule](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-chain-rule);
the \log n row is untouched. In particular the stationarity condition
\mathbf{V}^{\top}\mathbf{S}(\boldsymbol{p})\nabla\_{\boldsymbol{p}}\ell=\boldsymbol{0}
is the same equation with \nabla\_{\boldsymbol{p}}\ell taken from
[Eq. 20](#eq-score-general) and the tuple [Eq. 43](#eq-av-tuple).

#### Layer 2: a global scale, per sample or shared

On the raw scale the bulk is \boldsymbol{y}\_{\cdot
i}=a_in_i\tilde{\boldsymbol{y}}\_{\cdot i}, with the capture ratio a_i
of 1:

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},a_i,n_i
\sim\mathcal{N}\_G\Bigl(
a_in_i\tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot i},\\
a_i^{2}n_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i}) \Bigr)
=\mathcal{N}\_G\Bigl( c_i\tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot
i},\\ \phi_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i}) \Bigr), \qquad
c_i=a_in_i,\quad\phi_i=a_i^{2}n_i . \tag{44}

It is [Eq. 33](#eq-cc-scale-law) with \boldsymbol{m}\_{\cdot
i}=n_i\boldsymbol{p}\_{\cdot i}. The two scales enter only through the
mean scale c_i and the covariance scale \phi_i, and
(a_i,n_i)\mapsto(c_i,\phi_i) is a bijection of (0,\infty)^{2}, with
a_i=\phi_i/c_i and n_i=c_i^{2}/\phi_i. Here the simplex constraint buys
identifiability: because \mathbf{1}^{\top}\boldsymbol{p}\_{\cdot i}=1,
the mean pins c_i from first moments alone, and \phi_i is read from the
covariance magnitude. In Strategy 1 the mean fixes only the product
a_i\boldsymbol{m}\_{\cdot i}. In the (\boldsymbol{p},c,\phi) chart the
tuple is

\boldsymbol{\mu}\_{p_j}=c\tilde{\boldsymbol{\mu}}\_{\cdot j},\quad
\boldsymbol{\mu}\_{c}=\tilde{\boldsymbol{\mu}}\boldsymbol{p},\quad
\boldsymbol{\mu}\_{p_jc}=\tilde{\boldsymbol{\mu}}\_{\cdot j},\quad
\boldsymbol{\Sigma}\_{p_j}=\phi\boldsymbol{\Sigma}\_{W,j},\quad
\boldsymbol{\Sigma}\_{\phi}=\boldsymbol{S}\_W(\boldsymbol{p}),\quad
\boldsymbol{\Sigma}\_{p_j\phi}=\boldsymbol{\Sigma}\_{W,j}, \tag{45}

and every other derivative is zero: the mean is bilinear in
(c,\boldsymbol{p}) and the covariance bilinear in (\phi,\boldsymbol{p}).
At fixed \boldsymbol{p}, \hat c is the generalised least-squares
coefficient of \boldsymbol{y} on \tilde{\boldsymbol{\mu}}\boldsymbol{p}
and \hat\phi the mean squared Mahalanobis residual, as in
[Eq. 36](#eq-cc-score-a), so both scales can be profiled out. When N
bulk samples share one protocol, a **shared** capture ratio a ties the
two scales, \phi_i=ac_i, and the pooled log-likelihood \sum_i\ell_i
couples all samples through the a row of a block arrow-shaped Hessian.
This pooling is what `MuSiC` achieves with its sample-level normalising
constant and `RCTD` with a known depth offset
([`RCTD`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-rctd)).

**A Poisson reading of the scales (our remark).** The mean
a_in_i\sum_jp\_{ji}S_j\boldsymbol{\theta}\_{\cdot j} has the
multiplicative structure of a Poisson measurement model: `RCTD` writes
y\_{gi}\sim\mathrm{Poisson}(N_i\sum_jw\_{ji}\bar\mu\_{gj}) with
w\_{ji}=p\_{ji}e^{\alpha_i} ([Cable et al.
2022](#ref-cableRobustDecompositionCell2022a)), and Sarkar and Stephens
x\_{gc}\sim\mathrm{Poisson}(s_cb_g\lambda\_{gc}) ([Sarkar and Stephens
2021](#ref-sarkarSeparatingMeasurementExpression2021)). Here the rate is
carried by cell type rather than by gene: an average type-j cell
contributes S_j\boldsymbol{\theta}\_{\cdot j}, multiplied by its number
of cells and by the capture ratio. The variance does not follow the
Poisson mean–variance law; it comes from the cell-to-cell covariance
\boldsymbol{\Sigma}\_{W,j}. The scale (n in `RNA-Sieve`, N_ie^{\alpha_i}
in `RCTD`) is the total number of cells only when bulk and reference
share the platform, or for paired samples with similar biological and
technical size factors ([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)); otherwise
n_i is an effective number of reference-equivalent cells.

#### Layer 3: a residual bulk term

Adding
\boldsymbol{\epsilon}\_i\sim\mathcal{N}\_G(\boldsymbol{0},\boldsymbol{D}\_i)
gives

\boldsymbol{\Sigma}(\boldsymbol{p},\phi,\boldsymbol{\tau}^{2})
=\phi\boldsymbol{S}\_W(\boldsymbol{p})+\boldsymbol{D}\_i, \tag{46}

with
\boldsymbol{\Sigma}\_{\tau_g^{2}}=\boldsymbol{e}\_g\boldsymbol{e}\_g^{\top}
as in [Eq. 38](#eq-cc-score-tau). \phi_i and \boldsymbol{D}\_i both
inflate the variance and are separated by structure, not by magnitude:
the network entries of \boldsymbol{S}\_W(\boldsymbol{p}) identify
\phi_i, the diagonal excess identifies \boldsymbol{D}\_i. If the
reference network is weak, a free \boldsymbol{D}\_i absorbs \phi_i, and
n_i becomes unidentifiable; constrain \boldsymbol{D}\_i as in
[Sec. 4.1.3](#sec-cc-noise), fix n_i externally, or drop
\boldsymbol{D}\_i, as the [final
model](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-final-choice)
does.

#### Layer 4: within and between covariance from R replicates

With R\ge 2 replicates, multiply [Eq. 39](#eq-cc-replicate-type) by a_i
and set \boldsymbol{m}\_{\cdot i}=n_i\boldsymbol{p}\_{\cdot i}: the
between term
a_i^{2}m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}=c_i^{2}p\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}
keeps its quadratic weight, the within term becomes
\phi_ip\_{ji}\boldsymbol{\Sigma}\_{W,j}. The four-layer law is
[Eq. 8](#eq-three-layer) on the scale of the bulk library,

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},c_i,\phi_i
\sim\mathcal{N}\_G\Bigl(
c_i\tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot i},\\
c_i^{2}\boldsymbol{S}\_B(\boldsymbol{p}\_{\cdot
i})+\phi_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i})
+\boldsymbol{D}\_i \Bigr), \tag{47}

and the covariance derivatives in \boldsymbol{p} are

\boldsymbol{\Sigma}\_{p_j}
=2c^{2}p_j\boldsymbol{\Sigma}\_{B,j}+\phi\boldsymbol{\Sigma}\_{W,j},
\qquad \boldsymbol{\Sigma}\_{p_jp_j}=2c^{2}\boldsymbol{\Sigma}\_{B,j},
\qquad \boldsymbol{\Sigma}\_{p_jp\_\ell}=\boldsymbol{0}\\ (j\neq\ell).
\tag{48}

Setting n_i\to\infty at fixed c_i=1 (so \phi_i\to 0),
\boldsymbol{D}\_i=\boldsymbol{0}, returns the released DeCovarT
derivatives with \boldsymbol{\mu}=\tilde{\boldsymbol{\mu}} on the bulk
scale and \boldsymbol{\Sigma}\_j=\boldsymbol{\Sigma}\_{B,j}. Released
DeCovarT is therefore the large-n_i, noise-free limit of this
formulation, which is why its \boldsymbol{\Sigma}\_j must be a
between-replicate covariance to carry a sampling interpretation. The
diagonal-\boldsymbol{\Sigma}\_{B,j} remark of
[Sec. 4.1.4](#sec-cc-replicates) applies: a diagonal quadratic layer
plus a sparse-precision linear layer has a dense precision, so the
backends of the derivatives vignette see only the structure of
\boldsymbol{\Sigma}.

#### Strategy 2, DECALS and MEAD

`DECALS` and `MEAD` sit on this branch; their models are in the review
([`DECALS`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-decals),
[`MEAD`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-mead)).
In our indices the `DECALS` bulk covariance is
\operatorname{Cov}(\boldsymbol{y}\_{\cdot
i})=\sum_jp\_{ji}^{2}\boldsymbol{\Sigma}\_j (their Section 2.3), the
quadratic layer \boldsymbol{S}\_B(\boldsymbol{p}) and nothing else ([Cai
et al. 2024](#ref-caiStatisticalInferenceCelltype2024)). Three contrasts
with [Eq. 47](#eq-av-four-layer) remain.

- **Source of \boldsymbol{\Sigma}\_j.** `DECALS` regresses residual
  outer products on \hat p\_{ji}^{2} across bulk samples, which needs
  N\ge J samples of varying composition; DeCovarT plugs in reference
  estimates. With R=1, that regression is a route to
  \boldsymbol{\Sigma}\_{B,j}.
- **Estimator.** `DECALS` keeps constrained least squares for
  \hat{\boldsymbol{p}} and uses \boldsymbol{\Sigma}\_j only in a
  sandwich variance; DeCovarT maximises the Gaussian likelihood, whose
  determinant and covariance-quadratic terms no least-squares criterion
  has ([GLS
  competitor](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-extensions-roadmap.html#sec-robust-gls)).
- **Finite-sample correction.** The `DECALS` correction targets its
  covariance estimator, not \hat{\boldsymbol{p}} ([Firth
  penalisation](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-firth)).

`MEAD` estimates \boldsymbol{\beta}\_i=\gamma_i\boldsymbol{p}\_{\cdot
i}, the (c_i,\boldsymbol{p}\_{\cdot i}) pair of
[Eq. 44](#eq-av-scale-law), and normalises at the end ([Xie et al.
2026](#ref-xieStatisticalInferenceCell2026)). Its random platform ratio
\lambda_g, our \tilde b_g, is a mean-proportional diagonal inflation
shared by all samples, and its reference covariance is the
errors-in-variables term of [Sec. 6](#sec-summary); neither
\boldsymbol{S}\_W nor \boldsymbol{S}\_B enters its estimating equation.

> **Important 5: Strategies 1 and 2 are one model with one redundant
> coordinate**
>
> Start from the five-layer law [Eq. 40](#eq-cc-five-layer) and set
> \boldsymbol{m}\_{\cdot i}=n_i\boldsymbol{p}\_{\cdot i}:
>
> \begin{aligned} a_i\tilde{\boldsymbol{\mu}}\boldsymbol{m}\_{\cdot i}
> &=c_i\\\tilde{\boldsymbol{\mu}}\boldsymbol{p}\_{\cdot i}, \\
> a_i^{2}\bigl\[\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
> i})+\boldsymbol{S}\_B(\boldsymbol{m}\_{\cdot i})\bigr\]
> &=\phi_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot
> i})+c_i^{2}\boldsymbol{S}\_B(\boldsymbol{p}\_{\cdot i}), \end{aligned}
> \tag{49}
>
> with c_i=a_in_i and \phi_i=a_i^{2}n_i, which is
> [Eq. 47](#eq-av-four-layer) exactly. The map (\boldsymbol{m}\_{\cdot
> i},a_i)\mapsto(\boldsymbol{p}\_{\cdot i},c_i,\phi_i) is a bijection
> from (0,\infty)^{J}\times(0,\infty) onto
> \mathring{\Delta}^{J-1}\times(0,\infty)^{2}, so the two strategies are
> reparametrisations of the same Gaussian family. The cellular ratio
> p\_{ji} is the count m\_{ji} averaged by n_i, and the quadratic
> DeCovarT weight is the between-population layer written per cell:
> m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}/n_i^{2}=p\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}.
>
> What differs is how the redundant degree of freedom is spent. Strategy
> 1 has J+1 free parameters (\boldsymbol{m},a) and the mean identifies
> only a\boldsymbol{m}; the count scale comes from the covariance.
> Strategy 2 has (J-1)+2 free parameters (\boldsymbol{p},c,\phi); the
> simplex removes one coordinate and lets the mean identify c outright,
> leaving \phi to the covariance. The integer constraint of
> **?@sec-cc-integer** exists only in Strategy 1, because
> n_i\boldsymbol{p}\_{\cdot i} is not constrained to the lattice once
> n_i is a continuous nuisance. Choose by the question: counts and a
> scale calibrated to the reference platform (Strategy 1), or a
> composition with n_i as a nuisance (Strategy 2).

### Strategy 3: library-size normalisation gives a composition

The third formulation divides the bulk by its own total rather than by a
cell count, as CPM, TPM, RPKM and FPKM do. Let
\boldsymbol{T}\_i=\sum_j\sum_c\boldsymbol{X}\_{jic} be the unnormalised
vector and

\boldsymbol{Z}\_i=\frac{\boldsymbol{T}\_i}{\mathbf{1}^{\top}\boldsymbol{T}\_i},
\qquad \mathbf{1}^{\top}\boldsymbol{Z}\_i=1 . \tag{50}

\boldsymbol{Z}\_i is a random ratio. Even when \boldsymbol{T}\_i is
exactly Gaussian, \boldsymbol{Z}\_i is not; its support is the
composition hyperplane and a Gaussian description lives on the
(G-1)-dimensional tangent space. With
\boldsymbol{t}\_i=\mathbb{E}\[\boldsymbol{T}\_i\],
\boldsymbol{V}\_i=\operatorname{Cov}(\boldsymbol{T}\_i),
L_i=\mathbf{1}^{\top}\boldsymbol{t}\_i and
\boldsymbol{z}\_i=\boldsymbol{t}\_i/L_i, the delta method gives

\boldsymbol{J}\_i=\frac{1}{L_i}\bigl(\mathbf{I}\_G-\boldsymbol{z}\_i\mathbf{1}^{\top}\bigr),
\qquad
\boldsymbol{Z}\_i\approx\mathcal{N}\bigl(\boldsymbol{z}\_i,\boldsymbol{J}\_i\boldsymbol{V}\_i\boldsymbol{J}\_i^{\top}\bigr),
\qquad \mathbf{1}^{\top}\boldsymbol{J}\_i=\boldsymbol{0}. \tag{51}

The covariance is singular in \mathbb{R}^{G}. Closure also **creates**
correlation. For a diagonal
\boldsymbol{V}\_i=\operatorname{diag}(\boldsymbol{v}),

\bigl(\boldsymbol{J}\_i\boldsymbol{V}\_i\boldsymbol{J}\_i^{\top}\bigr)\_{gg'}
=\frac{1}{L_i^{2}}\Bigl(
v_g\delta\_{gg'}-z_gv\_{g'}-z\_{g'}v_g+z_gz\_{g'}\sum_hv_h \Bigr),
\tag{52}

so independent genes become negatively correlated after normalisation,
and the gene network estimated on raw counts is no longer the network of
the normalised bulk. The target also changes: with
\tilde{\boldsymbol{\mu}}\_{\cdot j}=S_j\boldsymbol{\theta}\_{\cdot j}
the normalised mean is \sum_jq\_{ji}\boldsymbol{\theta}\_{\cdot j} with
RNA fractions q\_{ji}=m\_{ji}S_j/\sum\_\ell m\_{\ell i}S\_\ell, which
equal cell fractions only when S_j is constant
([Sec. 5](#sec-uncoupling)). This vignette does not develop Strategy 3
further: it removes the count scale that Strategies 1 and 2 estimate,
breaks the additive mixture that the convolution needs ([normalisation
that preserves the linear mixture](#sec-norm-linear)), and is the regime
in which `MuSiC` and `DECALS` operate by choice.

### Strategy 4: plug-in scales before deconvolution

Strategies 1–3 estimate the scales from the bulk sample that is being
deconvolved. Strategy 4 estimates them **before** the deconvolution,
from auxiliary data, and holds them fixed in the regression or
likelihood. It is the strategy that most published methods use in
practice, usually without saying so.

![](theory-decovart-single-cell-tailored_files/mediabag/fig-dag-plugin.svg)

Figure 6: Joint versus plug-in estimation of nuisance scales
\boldsymbol{\eta} (size factors, cell sizes, platform, sample scale).
(a) Strategies 1–3 maximise the bulk likelihood in
\boldsymbol{p}\_{\cdot i} and \boldsymbol{\eta}\_i together. (b)
Strategy 4 estimates \hat{\boldsymbol{\eta}} from auxiliary data, then
fits \boldsymbol{p}\_{\cdot i} with \hat{\boldsymbol{\eta}} held fixed.

#### What is plugged in

| Scale | Plug-in estimator | Source | Methods |
|----|----|----|----|
| Technical size factor s^{\mathrm{tech}}\_c, extrinsic noise \nu^{\mathrm{tech}} | [the Poisson-gene size factor](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-poisson-size-factor) on Poisson genes; spike-ins | reference cells | Fang and Pachter; [the denoised within-type covariance](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-denoised-within) |
| Cell size S_j | mean UMI total per type; CLTS; RNA-content assay | reference, or external assay | `MuSiC`, `ReDeconv`, `EPIC`, `quanTIseq` |
| Reference scale a^{\mathrm{ref}}\_r | grand mean of the reference matrix | reference alone | `MEAD` (\hat\gamma_j) |
| Gene-wise platform ratio \tilde{\boldsymbol{b}} | regression on paired samples; moment matching; pseudo-bulk MLE | paired cohort, or all bulks | `Bisque`, `RCTD` stage 1 |
| Sample scale c_i | observed library size as an offset | bulk total | `RCTD` (N_i) |
| Type moments \tilde{\boldsymbol{\mu}}\_{\cdot j}, \boldsymbol{\Sigma}\_{W,j} | reference averages and covariances | reference | every method except `RNA-Sieve` and `MEAD` |

Table 2: Quantities that the reviewed methods estimate before the
deconvolution and then hold fixed.

With the plug-ins, the bulk law of Strategy 2 becomes

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot
i};\hat{\boldsymbol{\eta}} \sim\mathcal{N}\_G\Bigl( \hat
c_i\\\hat{\tilde{\mathbf{B}}}\sum_jp\_{ji}\hat
S_j\hat{\boldsymbol{\theta}}\_{\cdot j},\\
\hat{\tilde{\mathbf{B}}}\bigl\[\hat
c_i^{2}\boldsymbol{S}\_B(\boldsymbol{p}\_{\cdot
i})+\hat\phi_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot
i})\bigr\]\hat{\tilde{\mathbf{B}}} +\boldsymbol{D}\_i\Bigr), \tag{53}

with
\hat{\tilde{\mathbf{B}}}=\operatorname{diag}(\hat{\tilde{\boldsymbol{b}}}),
plugged-in sizes \hat S_j and relative abundances
\hat{\boldsymbol{\theta}}\_{\cdot j}, and \boldsymbol{S}\_W built from
[the denoised within-type
covariance](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-denoised-within).
The tuple of [Theorem 1](#thm-score-general) is that of
[Sec. 4.2.4](#sec-av-replicates) with \tilde{\boldsymbol{\mu}}\_{\cdot
j} replaced by \hat c_i\hat{\tilde{\mathbf{B}}}\hat
S_j\hat{\boldsymbol{\theta}}\_{\cdot j} and every covariance block
sandwiched by \hat{\tilde{\mathbf{B}}}; only \boldsymbol{p}\_{\cdot i}
(and, if kept free, \boldsymbol{\tau}^{2}) remains in the score.

#### Why plug in: identifiability

Some scales are not identified from the bulk sample at all, so a plug-in
is not a convenience but the only route.

- **Cell size and composition.** Only p\_{ji}S_j enters the mean
  ([Eq. 9](#eq-music-mean)). Fitting with a wrong \hat S_j makes the
  estimator converge to \tilde p\_{ji}\propto p\_{ji}S_j/\hat S_j: the
  cell fractions are rescaled by the ratio of true to plugged-in sizes.
  That is `MEAD`’s cell-size rescaling and `MuSiC`’s assumption A3'.
- **Gene-wise platform and composition.** A free \tilde{\boldsymbol{b}}
  fitted on the sample being deconvolved absorbs any residual, and with
  perfect markers it absorbs a type-specific factor (`MEAD`, Theorem 1).
  `RCTD` and `Bisque` estimate it on all bulks or on a paired cohort;
  `MEAD` declines to estimate it and treats it as noise.
- **Technical versus biological size factors.** One cell’s counts
  identify only the product s^{\mathrm{bio}}\_cs^{\mathrm{tech}}\_c;
  spike-ins or ambient RNA separate them ([Fang and Pachter
  2025](#ref-fangExtrinsicBiologicalStochasticity2025)).

#### What plugging in costs: uncertainty

Holding \hat{\boldsymbol{\eta}} fixed treats it as exact. Let
\hat{\boldsymbol{p}} solve the bulk score
\boldsymbol{U}(\boldsymbol{p};\hat{\boldsymbol{\eta}})=\boldsymbol{0},
with \hat{\boldsymbol{\eta}} estimated from independent auxiliary data
with covariance \boldsymbol{V}\_{\eta}. A first-order expansion in both
arguments (the standard two-step M-estimation argument) gives

\operatorname{Var}(\hat{\boldsymbol{p}}) \approx
I\_{pp}^{-1}+I\_{pp}^{-1}I\_{p\eta}\\\boldsymbol{V}\_{\eta}\\I\_{\eta
p}I\_{pp}^{-1}, \tag{54}

with the blocks of [Eq. 24](#eq-schur-information). The naive plug-in
variance I\_{pp}^{-1} omits the second term and under-covers whenever
the cross block I\_{p\eta} is non-zero; `MuSiC` and `DECALS` treat the
reference means as exact in that sense, whereas `RNA-Sieve` and `MEAD`
propagate reference noise. Compare with joint estimation. By the
Woodbury identity,

I\_{p\mid\eta}^{-1}
=I\_{pp}^{-1}+I\_{pp}^{-1}I\_{p\eta}\bigl(I\_{\eta\eta}-I\_{\eta
p}I\_{pp}^{-1}I\_{p\eta}\bigr)^{-1}I\_{\eta p}I\_{pp}^{-1}, \tag{55}

so joint estimation is [Eq. 54](#eq-plugin-variance) with
\boldsymbol{V}\_{\eta} equal to the variance of \hat{\boldsymbol{\eta}}
estimated from the bulk itself. The comparison is therefore explicit: a
plug-in is more precise than joint estimation exactly when the auxiliary
data know the scale better than the bulk does, and it is biased, not
merely noisier, when the auxiliary data measure a different scale
(another protocol, another tissue, another condition).

#### Affine platform maps

1 used one scalar a_i. Capture efficiency in single-cell protocols is
gene dependent and technology dependent, so the alignment between
reference and bulk is in general an affine map. Let
\boldsymbol{o}\in\mathbb{R}^{G} be a background intercept and
\tilde{\mathbf{B}}\_j=\operatorname{diag}(\tilde{\boldsymbol{b}}\_j) a
gene-wise, type-wise platform ratio that also carries the sample scale.
Affine invariance of the Gaussian family gives exactly

\boldsymbol{y}\_{\cdot i}
=\boldsymbol{o}+\sum_j\tilde{\mathbf{B}}\_j\boldsymbol{T}\_{ji}+\boldsymbol{\epsilon}\_i
\sim\mathcal{N}\_G\Bigl(
\boldsymbol{o}+\sum_jm\_{ji}\tilde{\mathbf{B}}\_j\tilde{\boldsymbol{\mu}}\_{\cdot
j},\\
\sum_jm\_{ji}\tilde{\mathbf{B}}\_j\boldsymbol{\Sigma}\_{W,j}\tilde{\mathbf{B}}\_j+\boldsymbol{D}\_i
\Bigr), \tag{56}

and the tuple of [Theorem 1](#thm-score-general) follows by replacing
\tilde{\boldsymbol{\mu}}\_{\cdot j} with
\tilde{\mathbf{B}}\_j\tilde{\boldsymbol{\mu}}\_{\cdot j} and
\boldsymbol{\Sigma}\_{W,j} with
\tilde{\mathbf{B}}\_j\boldsymbol{\Sigma}\_{W,j}\tilde{\mathbf{B}}\_j.
Three special cases are already in the literature.
\tilde{\mathbf{B}}\_j=a_i\mathbf{I}\_G, \boldsymbol{o}=\boldsymbol{0}
is 1. A gene- and subject-specific intercept \boldsymbol{o} with
\tilde{\mathbf{B}}\_j=\mathbf{I}\_G is what `MuSiC` adds to its equation
(8) to absorb protocol bias ([Wang et al.
2019](#ref-wangBulkTissueCell2019)).
\tilde{\mathbf{B}}\_j=\tilde{\mathbf{B}} shared across types, with a
per-gene intercept, is the `Bisque` transformation learned on paired
individuals ([Jew et al. 2020](#ref-jewAccurateEstimationCell2020)) and
the platform ratio of `MEAD` and `RCTD` ([Sec. 1.3](#sec-alignment)).

Two warnings carry over. A free \tilde{\boldsymbol{b}} fitted on the
sample being deconvolved is not identified from \boldsymbol{p}\_{\cdot
i}; it must come from a paired cohort, a pseudo-bulk of all samples, or
an external measurement model of the capture probabilities b_g of both
protocols ([Sarkar and Stephens
2021](#ref-sarkarSeparatingMeasurementExpression2021)). A type-wise
\tilde{\boldsymbol{b}}\_j is weaker still, and a candidate only when
sorted bulk profiles of type j exist on the bulk platform. Any such
plug-in is protocol specific, so \tilde{\mathbf{B}}\_j and S_j are
stored with the reference, never re-estimated per bulk sample; the
measurement-model estimates of S_j are in [Sec. 5](#sec-uncoupling).

## RNA fraction versus cell fraction

A linear engine fitted on a library-normalised signature
\boldsymbol{\theta}\_{\cdot j}, or on a normalised bulk, estimates the
**RNA mass fraction** \boldsymbol{q}\_{\cdot i}, not the cytometric cell
fraction \boldsymbol{p}\_{\cdot i}, unless transcriptome sizes S_j (mean
transcripts per cell; [Note 1](#nte-notation)) are homogeneous. Fitted
on the raw per-cell means \tilde{\boldsymbol{\mu}}\_{\cdot j}, it
returns \boldsymbol{p}\_{\cdot i} directly, provided their UMI totals
S_j reflect RNA content. With \hat q\_{j} an RNA-scale estimate,

\hat p\_{j} =\frac{\hat q\_{j}/S\_{j}}{\sum\_{\ell=1}^{J}\hat
q\_{\ell}/S\_{\ell}}. \tag{57}

**Post-correction** treats S_j as measured (or proxied). `EPIC` and
`quanTIseq` rescale \hat{\boldsymbol{q}} using kit-based mRNA content or
housekeeping-gene / proteasome-subunit surrogates ([Racle et al.
2017](#ref-racleSimultaneousEnumerationCancer2017); [Finotello et al.
2019](#ref-finotello_etal19)). `MuSiC` uses average cell-type library
size as a size proxy and warns that TPM can discard the information
needed to recover cell fractions ([Wang et al.
2019](#ref-wangBulkTissueCell2019)). `ReDeconv` separates global
library-size normalisation from cell-type transcriptome size ([Lu et al.
2025](#ref-luTranscriptomeSizeMatters2025)).

**Joint estimation** treats S_j as unknown. `MMAD` estimates extraction
efficiencies by non-linear conjugate gradients, so the mixture is no
longer linear ([Liebner et al.
2014](#ref-liebnerMMADMicroarrayMicrodissection2014)). A DeCovarT
analogue would introduce S_j inside the mean
\sum\_{j}p\_{ji}S\_{j}\boldsymbol{\theta}\_{\cdot j} (with a
product-simplex constraint on (p\_{j}S\_{j})), rather than assuming S_j
known. If S_j is observed, plug it in. Do not estimate a free gene-wise
platform ratio \tilde{\boldsymbol{b}} in the same fit
([Sec. 1.3](#sec-alignment)).

**Bayesian transcript totals.** With no RNA-content assay and no
cytometry ground truth, S_j has to come from a measurement model on the
single-cell counts. `bayNorm` takes the observed count of gene g in cell
c as a binomial draw from an unobserved original count, with
cell-specific capture efficiency b_c, and a negative binomial prior fit
by empirical Bayes ([Tang et al.
2020](#ref-tangBayNormBayesianGene2020)). The posterior is a shifted
negative binomial for that original count. Summing the posterior means
over genes is that cell’s inferred transcript total. S_j is the mean of
those totals inside type j. A single mean capture efficiency scales
every cell, so it cancels in the ratios of [Eq. 57](#eq-rna-to-cell);
the absolute molecule number still needs an external anchor (the paper
uses smFISH for that check). `Sanity` is the expression-model
counterpart of the split in [measurement and expression
models](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-measurement):
a Poisson measurement model and a prior on the expression state, with
the posterior reported as a denoised activity and an error bar, not as a
molecule count ([Breda et al.
2021](#ref-bredaBayesianInferenceGene2021)). Use the `bayNorm` total as
S_j. Use `Sanity` when \boldsymbol{\theta}\_{\cdot j} itself must be
corrected for sampling noise before it enters [the `MuSiC`
mean](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-single-cell-tailored.html#eq-music-mean).
The observed library size is a capture total, not this S_j.

## Summary and critical comparison

|  | Strategy 1: cell counts | Strategy 2: averaged mixture | Strategy 3: compositional | Strategy 4: plug-in scales | Released DeCovarT |
|----|----|----|----|----|----|
| Target | \boldsymbol{m}\_{\cdot i}\in\mathbb{N}^{J}, then \boldsymbol{p}=\boldsymbol{m}/n | \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1}, n_i nuisance | RNA fractions \boldsymbol{q}\_{\cdot i} | \boldsymbol{p}\_{\cdot i} given \hat{\boldsymbol{\eta}} | \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1} |
| Covariance weights | linear m\_{ji} in \boldsymbol{\Sigma}\_{W,j}; quadratic m\_{ji}^{2} in \boldsymbol{\Sigma}\_{B,j} | linear \phi_ip\_{ji}; quadratic c_i^{2}p\_{ji}^{2} | delta-method image of either | as Strategy 2, rescaled by \hat{\tilde{\boldsymbol{b}}}, \hat c_i | quadratic p\_{ji}^{2} only |
| Constraint | none; lattice | simplex via ILR | closure, singular | simplex via ILR | simplex via ILR |
| Scale identification | mean fixes a\boldsymbol{m}; covariance fixes a | mean fixes c; covariance fixes \phi | scale removed | from auxiliary data | c_i=1 |
| Reviewed methods | `RNA-Sieve`, `DWLS`, `MuSiC` eq. (2) | `DECALS`, `MEAD`, DeCovarT | `MuSiC`, `Bisque`, `DECALS` | `MuSiC` (S_j), `RCTD`, `Bisque`, `MEAD` (\hat\gamma) | – |
| Advantages | counts are the physical quantity; zeros are ordinary points; linear covariance, one Cholesky per iteration | DeCovarT is its noise-free limit; simplex identifies a; direct comparison with `DECALS`, `MEAD` | matches CPM/TPM pipelines; no scale to estimate | breaks non-identifiability of S_j and \tilde{\boldsymbol{b}}; uses external information | implemented; exact Gaussian; sparse \boldsymbol{\Omega}\_j as declared input |
| Drawbacks | a needs \boldsymbol{\Sigma}\_{W,j} on the reference molecule scale; integer search for rare types | \phi and \boldsymbol{D} separated only through off-diagonal structure | induced correlations; \boldsymbol{q}\neq\boldsymbol{p} unless S_j constant | biased if the auxiliary scale differs; variance needs [Eq. 54](#eq-plugin-variance) | \boldsymbol{\Sigma}\_j must be a between-replicate covariance; silences \boldsymbol{\Sigma}\_{W,j} for rare types |

Table 3: The four formulations and the released model.

[Table 3](#tbl-strategies) points to Strategy 2, with Strategy 1 as its
count-scale chart. For one UMI donor, one bulk sample and no covariates,
the model kept is its linear layer, [Eq. 44](#eq-av-scale-law), on raw
per-cell means, with both scales c_i,\phi_i estimated or profiled and no
\boldsymbol{S}\_B or \boldsymbol{D}\_i layer: \boldsymbol{S}\_B needs
R\ge 2 donors, and a bulk residual cannot be separated from the scales
without further assumptions. Its quasi-likelihood form and inference are
the [final
choice](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-final-choice)
of the inference vignette. With R\ge 2 donors the four layers of
[Sec. 4.2.4](#sec-av-replicates) apply, with S_j and
\tilde{\boldsymbol{b}} plugged in as in Strategy 4 and the within-type
covariance cleaned as in [the denoised within-type
covariance](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#eq-denoised-within).
Four further alternatives stay inside the multivariate Gaussian regime
and deserve a place in a simulation study before any is implemented.

1.  **Errors in variables for the signature.** The reference mean is
    itself an average of C_j cells, so
    \hat{\tilde{\boldsymbol{\mu}}}\_{\cdot
    j}\sim\mathcal{N}(\tilde{\boldsymbol{\mu}}\_{\cdot
    j},\boldsymbol{\Sigma}\_{W,j}/C_j). Marginalising it adds
    c_i^{2}\sum_jp\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/C_j to
    [Eq. 47](#eq-av-four-layer): a quadratic within-type term that is
    Monte Carlo error of the signature, not biology, and that vanishes
    as the reference grows. `RNA-Sieve` and `MEAD` carry it
    ([Erdmann-Pham et al.
    2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021); [Xie
    et al. 2026](#ref-xieStatisticalInferenceCell2026)). It is the one
    quadratic term that a one-replicate reference does identify.
2.  **Residual regression for \boldsymbol{\Sigma}\_{B,j}.** When R=1 but
    N\ge J bulk samples share the reference, the `DECALS` moment
    regression of residual outer products on \hat p\_{ji}^{2} estimates
    a between-population covariance without single-cell replicates. A
    diagonal restriction keeps it estimable.
3.  **Multi-sample pooling.** Shared a across samples
    ([Sec. 4.2.2](#sec-av-scale)), shared \boldsymbol{D} across
    technical replicates, or shared n across samples of the same tissue
    mass turn weakly identified per-sample scalars into pooled ones.
4.  **Compositional Gaussian on the bulk side.** If a pipeline insists
    on normalised bulk, model \operatorname{clr}(\boldsymbol{Z}\_i) or
    \operatorname{ilr}(\boldsymbol{Z}\_i) as Gaussian (a logistic-normal
    bulk) with mean the log-ratio image of the mixture ([Aitchison
    1982](#ref-aitchisonStatisticalAnalysisCompositional1982);
    [Pawlowsky-Glahn and Buccianti
    2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)). The mean
    map is no longer linear in \boldsymbol{p}, so the analytic score of
    [Theorem 1](#thm-score-general) still applies but
    \boldsymbol{\mu}\_{kl}\neq\boldsymbol{0}.

Outside the Gaussian regime, a Poisson–log-normal read-level mixture
(`RCTD`,
[`RCTD`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.html#sec-rctd))
is the natural competitor.

Strategy 4 is orthogonal to that choice: any of the laws above can be
fitted with plug-in scales. The rule that follows from
[Sec. 4.4](#sec-strategy-plugin) is to plug in what the bulk cannot
identify (S_j, \tilde{\boldsymbol{b}}, s^{\mathrm{tech}}\_c), to
estimate jointly what it identifies well (c_i and \phi_i on the simplex
chart), and in both cases to report [Eq. 54](#eq-plugin-variance) rather
than the naive I\_{pp}^{-1}.

## References

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Breda, Jérémie, Mihaela Zavolan, and Erik van Nimwegen. 2021. ‘Bayesian
Inference of Gene Expression States from Single-Cell RNA-Seq Data’.
*Nature Biotechnology* 39 (8): 1008–16.
<https://doi.org/10.1038/s41587-021-00875-x>.

Cable, Dylan M., Evan Murray, Luli S. Zou, et al. 2022. ‘Robust
Decomposition of Cell Type Mixtures in Spatial Transcriptomics’. *Nature
Biotechnology* 40 (4): 517–26.
<https://doi.org/10.1038/s41587-021-00830-w>.

Cai, Biao, Emma Jingfei Zhang, Hongyu Li, Chang Su, and Hongyu Zhao.
2024. ‘Statistical Inference of Cell-Type Proportions Estimated from
Bulk Expression Data’. *Journal of the American Statistical Association*
119 (548): 2521–32. <https://doi.org/10.1080/01621459.2024.2382435>.

Erdmann-Pham, Dan D., Jonathan Fischer, Justin Hong, and Yun S. Song.
2021. ‘Likelihood-Based Deconvolution of Bulk Gene Expression Data Using
Single-Cell References’. *Genome Research* 31 (10): 1794–806.
<https://doi.org/10.1101/gr.272344.120>.

Fang, Meichen, and Lior Pachter. 2025. *Extrinsic Biological
Stochasticity and Technical Noise Normalization of Single-Cell RNA
Sequencing Data*. bioRxiv preprint; bioRxiv.
<https://doi.org/10.1101/2025.05.11.653373>.

Finotello, Francesca, Clemens Mayer, Christina Plattner, et al. 2019.
‘Molecular and Pharmacological Modulators of the Tumor Immune Contexture
Revealed by Deconvolution of RNA-seq Data’. *Genome Medicine* 11.
<https://doi.org/10.1186/s13073-019-0638-6>.

Hafemeister, Christoph, and Rahul Satija. 2019. ‘Normalization and
Variance Stabilization of Single-Cell RNA-Seq Data Using Regularized
Negative Binomial Regression’. *Genome Biology* 20 (1): 296.
<https://doi.org/10.1186/s13059-019-1874-1>.

Jew, Brandon, Marcus Alvarez, Elior Rahmani, et al. 2020. ‘Accurate
Estimation of Cell Composition in Bulk Expression Through Robust
Integration of Single-Cell Information’. *Nature Communications* 11.
<https://doi.org/10.1038/s41467-020-15816-6>.

Liebner, David A., Kun Huang, and Jeffrey D. Parvin. 2014. ‘MMAD:
Microarray Microdissection with Analysis of Differences Is a
Computational Tool for Deconvoluting Cell Type-Specific Contributions
from Tissue Samples’. *Bioinformatics (Oxford, England)* 30.
<https://doi.org/10.1093/bioinformatics/btt566>.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. ‘Moderated
Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2’.
*Genome Biology* 15. <https://doi.org/10.1186/s13059-014-0550-8>.

Lu, Songjian, Jiyuan Yang, Lei Yan, et al. 2025. ‘Transcriptome Size
Matters for Single-Cell RNA-seq Normalization and Bulk Deconvolution’.
*Nature Communications* 16 (1).
<https://doi.org/10.1038/s41467-025-56623-1>.

Novershtern, Noa, Aravind Subramanian, Lee N. Lawton, et al. 2011.
‘Densely Interconnected Transcriptional Circuits Control Cell States in
Human Hematopoiesis’. *Cell* 144 (2): 296–309.
<https://doi.org/10.1016/j.cell.2011.01.004>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Racle, Julien, Kaat de Jonge, Petra Baumgaertner, Daniel E Speiser, and
David Gfeller. 2017. ‘Simultaneous Enumeration of Cancer and Immune Cell
Types from Bulk Tumor Gene Expression Data’. *eLife* 6.
<https://doi.org/10.7554/elife.26476>.

Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. 2010. ‘edgeR:
A Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data’. *Bioinformatics* 26.
<https://doi.org/10.1093/bioinformatics/btp616>.

Sarkar, Abhishek, and Matthew Stephens. 2021. ‘Separating Measurement
and Expression Models Clarifies Confusion in Single-Cell RNA Sequencing
Analysis’. *Nature Genetics* 53 (6): 770–77.
<https://doi.org/10.1038/s41588-021-00873-4>.

Squair, Jordan W., Matthieu Gautier, Claudia Kathe, et al. 2021.
‘Confronting False Discoveries in Single-Cell Differential Expression’.
*Nature Communications* 12 (1).
<https://doi.org/10.1038/s41467-021-25960-2>.

Szaniszlo, Peter, Nan Wang, Mala Sinha, et al. 2004. ‘Getting the Right
Cells to the Array: Gene Expression Microarray Analysis of Cell Mixtures
and Sorted Cells’. *Cytometry. Part A* 59 (2): 191–202.
<https://doi.org/10.1002/cyto.a.20055>.

Tang, Wenhao, François Bertaux, Philipp Thomas, et al. 2020. ‘bayNorm:
Bayesian Gene Expression Recovery, Imputation and Normalization for
Single-Cell RNA-Sequencing Data’. *Bioinformatics* 36 (4): 1174–81.
<https://doi.org/10.1093/bioinformatics/btz726>.

Tsoucas, Daphne, Rui Dong, Haide Chen, Qian Zhu, Guoji Guo, and
Guo-Cheng Yuan. 2019. ‘Accurate Estimation of Cell-Type Composition from
Gene Expression Data’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-019-10802-z>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.

Xie, Dongyue, Lin Gui, and Jingshu Wang. 2026. ‘Statistical Inference
for Cell Type Deconvolution’. *Journal of the Royal Statistical Society
Series B: Statistical Methodology* 88 (4): 1449–68.
<https://doi.org/10.1093/jrsssb/qkag054>.

Zimmerman, Kip D., Mark A. Espeland, and Carl D. Langefeld. 2021. ‘A
Practical Solution to Pseudoreplication Bias in Single-Cell Studies’.
*Nature Communications* 12 (1): 738.
<https://doi.org/10.1038/s41467-021-21038-1>.
