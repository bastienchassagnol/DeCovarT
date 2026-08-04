# Feature selection for reference-based deconvolution

## 1 Scope

This vignette surveys **feature selection for reference-based
(signature-based) deconvolution**—methods that recover cellular
proportions \boldsymbol{p} from a bulk mixture using a purified
expression design matrix ([Gaspard-Boulinc et al.
2025](#ref-gaspard-boulincCelltypeDeconvolutionMethods2025); [Avila
Cobos et al.
2018](#ref-avilacobosComputationalDeconvolutionTranscriptomics2018);
[Sturm et al.
2019](#ref-sturmComprehensiveEvaluationTranscriptomebased2019)). Feature
selection itself is a classical pre-modelling filter ([Guyon and
Elisseeff 2003](#ref-guyonIntroductionVariableFeature2003)). We
deliberately exclude **marker-based** abundance scoring (enrichment /
GSEA-style modules), which returns proxies of presence rather than a
consensual gene panel suited to linear or probabilistic mixture models
([Aran et al. 2017](#ref-aranXCellDigitallyPortraying2017); [Zhang et
al. 2017](#ref-zhang_etal17)).

The focus is **second-generation** signatures built from single-cell
pseudobulks (CIBERSORTx, MuSiC, DWLS, Kassandra and related pipelines)
([Newman et al. 2019](#ref-newmanDeterminingCellType2019); [Wang et al.
2019](#ref-wangBulkTissueCell2019); [Tsoucas et al.
2019](#ref-tsoucasAccurateEstimationCelltype2019); [Zaitsev et al.
2022](#ref-zaitsevPreciseReconstructionTME2022)), where the selected
genes must discriminate **all** cell types in the signature jointly, not
only one-versus-rest markers. Correlation-aware compact panels such as
DUBStepR illustrate the same goal from a single-cell feature-selection
angle ([Ranjan et al.
2021](#ref-ranjanDUBStepRScalableCorrelationbased2021)).

Interactive signature refinement utilities (zero filtering, specificity
binning, Gini / entropy compaction) are implemented in
[DeconvExplorer](https://github.com/omnideconv/deconvExplorer)
(omnideconv ecosystem; cite as software once added to `packages.bib`).

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  A["Reference μ, Σ from sc / bulk"] --> B["Univariate"]
  A --> C["Multivariate"]
  B --> B1["One vs others"]
  B --> B2["All vs all"]
  C --> C1["One vs others\n(GGM / DiffNet)"]
  C --> C2["All vs all\n(κ, D/E-opt, overlap)"]
  B2 --> D["Consensual panel G"]
  C2 --> D
  D --> E["Perspectives:\nmechanistic / potency / archetypes"]
```

Figure 1: Layout of this vignette: univariate and multivariate
selection, each split into one-versus-others then all-versus-all,
followed by perspectives for continuous cell states.

## 2 Notation

We align notation with the DeCovarT generative model. Let G genes and J
cell types index the purified means
\boldsymbol{\mu}=(\boldsymbol{\mu}\_{.j})\_{j=1}^{J}\in\mathbb{R}^{G\times
J} and (optional) covariances \\\boldsymbol{\Sigma}\_{j}\\\_{j=1}^{J}. A
bulk observation obeys

\boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p}) \sim
\mathcal{N}\_{G}\\\left( \boldsymbol{\mu}\boldsymbol{p},\\
\sum\_{j=1}^{J}p\_{j}^{2}\boldsymbol{\Sigma}\_{j} \right), \quad
\boldsymbol{p}\in\Delta^{J-1}, \tag{1}

with \boldsymbol{\zeta}=(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_{j}\\)
and simplex \Delta^{J-1}=\\\boldsymbol{p}:p\_{j}\ge
0,\sum\_{j}p\_{j}=1\\. For a gene subset
\mathcal{G}\subseteq\\1,\ldots,G\\, write
\boldsymbol{\mu}\_{\mathcal{G}} for the corresponding rows. The
condition number of a design matrix \boldsymbol{X} is

\kappa(\boldsymbol{X}) =
\frac{\lambda\_{\max}(\boldsymbol{X})}{\lambda\_{\min}(\boldsymbol{X})},
\tag{2}

with large \kappa indicating multicollinearity and numerical instability
([Gong and Szustakowski
2013](#ref-gongDeconRNASeqStatisticalFramework2013); [Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)).

## 3 Preprocessing shared by all strategies

Before scoring genes:

1.  Aggregate **donor-level** cell-type pseudobulks (not individual
    cells as independent replicates). Pseudobulking reduces sparsity and
    avoids treating cells as independent biological replicates in
    differential testing ([Squair et al.
    2021](#ref-squairConfrontingFalseDiscoveries2021)).
2.  Estimate gene–cell-type means \mu\_{gj} and variance components on a
    **linear** scale (counts, CPM, or TPM). Log-space breaks the mixture
    linearity in [Equation 1](#eq-mixture).
3.  Filter genes with poor bulk detectability, extreme zero inflation
    after pseudobulking, strong donor instability, or systematic sc–bulk
    discordance.
4.  When the biological target is **cell number** rather than RNA mass,
    correct transcriptome-size differences across cell types ([Lu et al.
    2025](#ref-luTranscriptomeSizeMatters2025)).
5.  Apply the **same** normalisation to bulk and signature ([Jin and Liu
    2021](#ref-jinBenchmarkRNAseqDeconvolution2021); [Racle et al.
    2017](#ref-racleSimultaneousEnumerationCancer2017)).

## 4 Univariate approaches

Univariate selectors score each gene marginally (or with a simple
contrast) and ignore partial correlations among genes.

### 4.1 One versus others

Classical signature pipelines score gene g for cell type j against the
remaining types, then concatenate per-type shortlists ([Abbas et al.
2009](#ref-abbasDeconvolutionBloodMicroarray2009); [Newman et al.
2015](#ref-newmanRobustEnumerationCell2015); [Finotello et al.
2019](#ref-finotello_etal19)).

#### 4.1.1 Differential expression and fold-change ranks

Limma–voom, edgeR and DESeq2 remain the workhorses for one-versus-rest
contrasts ([Ritchie et al.
2015](#ref-ritchieLimmaPowersDifferential2015a); [Robinson et al.
2010](#ref-robinsonEdgeRBioconductorPackage2010); [Love et al.
2014](#ref-loveModeratedEstimationFold2014a)). CIBERSORT retains genes
with adjusted p\le 0.3 ordered by decreasing fold-change before a global
condition-number prune ([Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)). Abbas et al. keep genes
that exceed the second- and third-ranked cell types, then choose the
panel with smallest \kappa ([Abbas et al.
2009](#ref-abbasDeconvolutionBloodMicroarray2009)). quanTIseq /
xCell-style filters further drop non-haematopoietic contaminants
([Finotello et al. 2019](#ref-finotello_etal19); [Aran et al.
2017](#ref-aranXCellDigitallyPortraying2017)).

#### 4.1.2 ANOVA F-test

Wang et al. rank genes by the nested-model F statistic comparing a
cell-type factor to an intercept-only model ([Wang et al.
2010](#ref-wangSilicoEstimatesTissue2010)). For gene g with expression
x\_{gi} in sample i and cell-type labels z\_{i}\in\\1,\ldots,J\\,

F\_{g} = \frac{ \bigl(\mathrm{SSR}\_{0}-\mathrm{SSR}\_{1}\bigr)/(J-1) }{
\mathrm{SSR}\_{1}/(n-J) }, \tag{3}

where \mathrm{SSR}\_{0}=\sum\_{i}(x\_{gi}-\bar{x}\_{g})^{2} and
\mathrm{SSR}\_{1}=\sum\_{j}\sum\_{i:z\_{i}=j}(x\_{gi}-\bar{x}\_{gj})^{2}.
Large F\_{g} indicates that cell-type means explain a substantial share
of gene variance.

#### 4.1.3 Gini and entropy specificity

BioQC and DeconvExplorer compact signatures with Gini or entropy scores
over the row \boldsymbol{\mu}\_{g\cdot} ([Zhang et al.
2017](#ref-zhang_etal17)). With relative shares
q\_{gj}=\mu\_{gj}/\sum\_{\ell}\mu\_{g\ell},

\mathrm{Gini}(g) = \frac{J}{J-1} \left( 1-2\sum\_{j=1}^{J}
\frac{J-j+\tfrac{1}{2}}{J}\\q\_{g(j)} \right), \quad H(g) =
-\sum\_{j=1}^{J}q\_{gj}\log q\_{gj}, \tag{4}

where q\_{g(j)} are the ordered shares. High Gini (or low entropy) flags
genes concentrated in few cell types—useful for one-versus-rest
compaction, but prone to redundancy when panels are concatenated.

> **Limitation of one-versus-rest concatenation**
>
> Union of per-type markers often yields a multicollinear
> \boldsymbol{\mu}\_{\mathcal{G}} that is still too large for stable
> NNLS / SVR / GLS ([Newman et al.
> 2015](#ref-newmanRobustEnumerationCell2015); [Avila Cobos et al.
> 2018](#ref-avilacobosComputationalDeconvolutionTranscriptomics2018)).
> A global all-versus-all refinement is therefore required.

### 4.2 All versus all

#### 4.2.1 Condition-number pruning

After candidate generation, regression-based methods select \mathcal{G}
to minimise \kappa(\boldsymbol{\mu}\_{\mathcal{G}})
([Equation 2](#eq-condition-number)), optionally after projecting out
the sum-to-one direction ([Abbas et al.
2009](#ref-abbasDeconvolutionBloodMicroarray2009); [Gong and
Szustakowski 2013](#ref-gongDeconRNASeqStatisticalFramework2013);
[Newman et al. 2015](#ref-newmanRobustEnumerationCell2015); [Boldina et
al. 2022](#ref-boldinaA2SignAgnosticAlgorithms2022)). This is the
historical gold-standard metric for LM22-style panels and remains the
default in many second-generation pipelines.

#### 4.2.2 Genetic algorithms (AutoGeneS)

AutoGeneS searches gene subsets with a multi-objective genetic algorithm
that **minimises inter-population correlation** while **maximising
centroid distance** ([Aliee and Theis
2021](#ref-alieeAutoGeneSAutomaticGene2021)). Related multi-objective
genetic searches appear in network module discovery (e.g. MOGAMUN)
([Novoa-del-Toro et al.
2021](#ref-novoa-del-toroMultiobjectiveGeneticAlgorithm2021)). AutoGeneS
is explicitly all-versus-all: the loss is defined on the full J-column
signature rather than on independent one-versus-rest lists.

DeCovarT replaces that dual objective with a single **global overlap**
criterion on the Gaussian mixture
\\(p\_{j},\boldsymbol{\mu}\_{.j},\boldsymbol{\Sigma}\_{j})\\ (see
[Equation 9](#eq-overlap) below and the synthetic-scenarios vignette).

## 5 Multivariate approaches

Multivariate selectors use joint dependence—covariance, precision, or
information matrices—so that a gene is scored by how it improves
identifiability of \boldsymbol{p} under [Equation 1](#eq-mixture).

### 5.1 One versus others

Univariate DGE ignores gene–gene interactions and therefore requires
aggressive multiple-testing correction without modelling network
rewiring. Differential-network and GGM methods address that gap.

#### 5.1.1 Differential networks (INDEED and related)

INDEED combines mean (fold-change) evidence with bootstrap tests of
changes in partial correlation, ranking genes that both shift in mean
**and** rewire neighbourhood structure ([Zuo et al.
2016](#ref-zuoINDEEDIntegratedDifferential2016)). Related ideas appear
in weighted graphical lasso pipelines that incorporate PPI priors and
differential-network scores (DNS). Relative to pure DGE, these methods
favour genes whose local Markov blanket differs between a focal cell
type and the rest.

#### 5.1.2 Sparse GGMs and cooperative penalties

Cell-type precision matrices
\boldsymbol{\Omega}\_{j}=\boldsymbol{\Sigma}\_{j}^{-1} estimated by
graphical lasso are shrunk; non-zero partial correlations are typically
underestimated. A practical mitigation is to take the gLasso **support**
as topological constraints and re-estimate by constrained MLE when the
graph is chordal (otherwise the mapping is non-trivial) ([Dahl et al.
2005](#ref-dahlMaximumLikelihoodEstimation2005)). Cooperative /
sign-coherent group lasso penalties can encode pathway modules and
PPI-informed edge weights ([Chiquet et al.
2012](#ref-chiquetSparsitySigncoherentGroups2012)). PPI priors also
shrink the combinatorial space of undirected graphs explored per cell
type.

> **From one-versus-others networks to a panel**
>
> Network-based one-versus-others lists still concatenate into a
> redundant design. They are best treated as **candidate generators**
> for the all-versus-all stage below.

### 5.2 All versus all

#### 5.2.1 Variance-weighted condition number

With gene-specific uncertainties \tau\_{g}^{2} (biological + sampling +
optional sc–bulk discordance), standardise rows and evaluate
conditioning in the (J-1)-dimensional contrast space induced by
\sum\_{j}p\_{j}=1:

\widetilde{\mu}\_{gj} = \frac{\mu\_{gj}-\bar{\mu}\_{g}}{\tau\_{g}},
\qquad \kappa\\\left(
\boldsymbol{W}\_{\mathcal{G}}^{1/2}\boldsymbol{\mu}\_{\mathcal{G}}\boldsymbol{C}
\right), \tag{5}

where \boldsymbol{W}\_{\mathcal{G}}=\mathrm{diag}(1/\tau\_{g}^{2}) and
\boldsymbol{C} removes the common-expression direction. Unweighted
\kappa(\boldsymbol{\mu}\_{\mathcal{G}}) alone can favour quiet but
biased genes ([Gong and Szustakowski
2013](#ref-gongDeconRNASeqStatisticalFramework2013)).

#### 5.2.2 D- / E- / A-optimal design

Under approximately Gaussian errors
\boldsymbol{y}\_{\mathcal{G}}\mid\boldsymbol{p}\sim\mathcal{N}(\boldsymbol{\mu}\_{\mathcal{G}}\boldsymbol{p},\boldsymbol{\Sigma}\_{\mathcal{G}}),
the information for \boldsymbol{p} (after the simplex constraint) is

\mathcal{I}\_{\mathcal{G}} = \boldsymbol{\mu}\_{\mathcal{G}}^{\top}
\boldsymbol{\Sigma}\_{\mathcal{G}}^{-1} \boldsymbol{\mu}\_{\mathcal{G}}.
\tag{6}

Common criteria are

\max\_{\mathcal{G}}\log\det(\mathcal{I}\_{\mathcal{G}}+\lambda\boldsymbol{I})
\quad\text{(D)}, \qquad
\max\_{\mathcal{G}}\lambda\_{\min}(\mathcal{I}\_{\mathcal{G}})
\quad\text{(E)}, \qquad
\min\_{\mathcal{G}}\operatorname{tr}(\mathcal{I}\_{\mathcal{G}}^{-1})
\quad\text{(A)}. \tag{7}

E-optimality is closely related to protecting the hardest cell-type
contrast; unweighted condition-number minimisation is a rough E-proxy
without an explicit noise model. In practice use diagonal or shrinkage
\boldsymbol{\Sigma}\_{\mathcal{G}}.

#### 5.2.3 Pairwise Gaussian separation with redundancy control

A cheap multivariate screen scores each gene by standardised pairwise
separation

d\_{g,jl}^{2} = \frac{(\mu\_{gj}-\mu\_{gl})^{2}}
{\sigma\_{gj}^{2}+\sigma\_{gl}^{2}+\varepsilon} \tag{8}

and solves a set-cover-like objective with diminishing returns and a
redundancy penalty. Recursive / hierarchical selection (broad lineage
then subtypes) mirrors MuSiC’s motivation for correlated cell types
([Wang et al. 2019](#ref-wangBulkTissueCell2019)). Follow with
[Equation 6](#eq-fisher) or [Equation 5](#eq-weighted-kappa) for the
final panel.

#### 5.2.4 Global overlap (DeCovarT proposal)

AutoGeneS optimises correlation and centroid distance separately ([Aliee
and Theis 2021](#ref-alieeAutoGeneSAutomaticGene2021)). We instead
minimise the **average pairwise overlap** of the purified Gaussians,
which folds mean separation **and** covariance structure into one
interpretable probability: the chance of mis-assigning a draw under a
MAP rule. With MixSim-style pairwise misclassification masses
\Omega\_{jl} ([Melnykov et al. 2024](#ref-R-MixSim)),

\overline{\Omega} = \frac{2}{J(J-1)} \sum\_{1\le j\<l\le J}
\bigl(p\_{j}\Omega\_{jl}+p\_{l}\Omega\_{lj}\bigr). \tag{9}

Package helper
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
implements [Equation 9](#eq-overlap). Lower \overline{\Omega}
corresponds to better-separated purified profiles and, empirically, to
more accurate recovery of \boldsymbol{p}.

#### 5.2.5 Wrapper optimisation on pseudo-bulks and matched bulks

When matched sc/bulk aliquots exist (e.g. Kassandra-style benchmarks)
([Zaitsev et al. 2022](#ref-zaitsevPreciseReconstructionTME2022)),
evaluate candidate panels by running the **intended** deconvolution
algorithm on held-out mixtures and scoring RMSE / worst-type error /
rare-cell recall under nested **donor-level** cross-validation. Use the
one-standard-error rule to prefer the smallest panel statistically
indistinguishable from the best. Wrappers overfit the mixture simulator
if real matched bulks are omitted ([Sturm et al.
2019](#ref-sturmComprehensiveEvaluationTranscriptomebased2019); [Jin and
Liu 2021](#ref-jinBenchmarkRNAseqDeconvolution2021)).

| Stage | Examples | Role |
|:---|:---|:---|
| Univariate · 1 vs others | DE / F-test / Gini–entropy | Candidate generation |
| Univariate · all vs all | κ pruning; AutoGeneS | Panel compaction |
| Multivariate · 1 vs others | INDEED; sparse GGM | Network-aware candidates |
| Multivariate · all vs all | Weighted κ; D/E-opt; overlap | Identifiability of p |
| Wrapper | Nested CV on matched bulks | Empirical size selection |

Table 1: Compact map from selection strategy to typical role in a
reference-based pipeline.

## 6 Recommended pipeline

1.  **Reference moments.** Donor-level pseudobulks; estimate
    \boldsymbol{\mu} and variance components on a linear scale.
2.  **Reliability filter.** Drop unstable / undetectable / discordant
    genes; optionally inflate \tau\_{g}^{2} by observed sc–bulk
    discrepancy.
3.  **Candidates.** Union of stable one-versus-others DE, high
    d\_{g,jl}^{2}, and (optionally) differential-network hits → roughly
    10^{2}–10^{3} genes.
4.  **Information-optimal subset.** Greedy D- or E-optimal selection
    1.  or weighted \kappa ([Equation 5](#eq-weighted-kappa)), with
        overlap ([Equation 9](#eq-overlap)) as a secondary monitor.
5.  **Wrapper prune.** Evaluate panel sizes along the path under nested
    donor CV and matched bulks; keep the smallest competitive panel.

The distinctive contribution for DeCovarT-style second-generation models
is a **bulk-calibrated, variance-weighted information criterion plus
overlap**, exploiting paired sc/bulk designs when available.

## 7 Perspectives

### 7.1 Mechanistic priors

PPI / TF networks and pathway modules can regularise both GGM estimation
and cooperative penalties ([Chiquet et al.
2012](#ref-chiquetSparsitySigncoherentGroups2012); [Zuo et al.
2016](#ref-zuoINDEEDIntegratedDifferential2016)), shrinking the search
toward biologically plausible supports rather than purely statistical
markers.

### 7.2 Continuous states, potency and archetypes

Many tissues are better described by **continua** (differentiation
trajectories, activation gradients) than by discrete labels. In that
regime, a finite J-column signature is a piecewise approximation of an
archetype geometry ([Hart et al.
2015](#ref-hartInferringBiologicalTasks2015)) or a potency axis ([Gulati
et al. 2020](#ref-gulatiSinglecellTranscriptionalDiversity2020)): genes
should separate extreme vertices **and** resolve interior mixtures along
trajectories. Open directions include replacing hard one-versus-others
contrasts with scores along inferred latent axes, and evaluating overlap
/ Fisher information on continuous rather than categorical components.
Until such methods mature, hierarchical discrete panels (lineage →
subtype) remain a pragmatic compromise ([Wang et al.
2019](#ref-wangBulkTissueCell2019)).

## 8 Software notes

- DeCovarT helpers:
  [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md),
  [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
  [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md).
- Signature QC / compaction patterns: DeconvExplorer
  (`removePercentZeros`, `removeUnspecificGenes`, `selectGenesByScore`
  with Gini or entropy).
- Gold-standard reference-based estimators for downstream validation:
  CIBERSORT / CIBERSORTx, EPIC, quanTIseq, ABIS, MuSiC, DWLS ([Newman et
  al. 2015](#ref-newmanRobustEnumerationCell2015); [Newman et al.
  2019](#ref-newmanDeterminingCellType2019); [Racle et al.
  2017](#ref-racleSimultaneousEnumerationCancer2017); [Finotello et al.
  2019](#ref-finotello_etal19); [Monaco et al.
  2019](#ref-monacoRNASeqSignaturesNormalized2019); [Wang et al.
  2019](#ref-wangBulkTissueCell2019); [Tsoucas et al.
  2019](#ref-tsoucasAccurateEstimationCelltype2019)).

Abbas, Alexander R., Kristen Wolslegel, Dhaya Seshasayee, Zora Modrusan,
and Hilary F. Clark. 2009. ‘Deconvolution of Blood Microarray Data
Identifies Cellular Activation Patterns in Systemic Lupus
Erythematosus’. *PloS One* 4.
<https://doi.org/10.1371/journal.pone.0006098>.

Aliee, Hananeh, and Fabian J. Theis. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems* 12. <https://doi.org/10.1016/j.cels.2021.05.006>.

Aran, Dvir, Zicheng Hu, and Atul Butte. 2017. ‘xCell: Digitally
Portraying the Tissue Cellular Heterogeneity Landscape’. *Genome
Biology* 18. <https://doi.org/10.1186/s13059-017-1349-1>.

Avila Cobos, Francisco, Jo Vandesompele, Pieter Mestdagh, and Katleen De
Preter. 2018. ‘Computational Deconvolution of Transcriptomics Data from
Mixed Cell Populations’. *Bioinformatics (Oxford, England)* 34.
<https://doi.org/10.1093/bioinformatics/bty019>.

Boldina, Galina, Paul Fogel, Corinne Rocher, Charles Bettembourg, George
Luta, and Franck Augé. 2022. ‘A2Sign: Agnostic Algorithms for
Signatures—a Universal Method for Identifying Molecular Signatures from
Transcriptomic Datasets Prior to Cell-Type Deconvolution’.
*Bioinformatics* 38 (4): 1015–21.
<https://doi.org/10.1093/bioinformatics/btab773>.

Chiquet, Julien, Yves Grandvalet, and Camille Charbonnier. 2012.
‘Sparsity with Sign-Coherent Groups of Variables via the
Cooperative-Lasso’. *The Annals of Applied Statistics* 6.
<https://doi.org/10.1214/11-aoas520>.

Dahl, Joachim, Vwani Roychowdhury, and Lieven Vandenberghe. 2005.
*Maximum Likelihood Estimation of Gaussian Graphical Models: Numerical
Implementation and Topology Selection*. University of California, Los
Angeles.

Finotello, Francesca, Clemens Mayer, Christina Plattner, et al. 2019.
‘Molecular and Pharmacological Modulators of the Tumor Immune Contexture
Revealed by Deconvolution of RNA-seq Data’. *Genome Medicine* 11.
<https://doi.org/10.1186/s13073-019-0638-6>.

Gaspard-Boulinc, Lucie C., Luca Gortana, Thomas Walter, Emmanuel
Barillot, and Florence M. G. Cavalli. 2025. ‘Cell-Type Deconvolution
Methods for Spatial Transcriptomics’. *Nature Reviews Genetics* 26.
<https://doi.org/10.1038/s41576-025-00845-y>.

Gong, Ting, and Joseph D. Szustakowski. 2013. ‘DeconRNASeq: A
Statistical Framework for Deconvolution of Heterogeneous Tissue Samples
Based on mRNA-Seq Data’. *Bioinformatics (Oxford, England)* 29.
<https://doi.org/10.1093/bioinformatics/btt090>.

Gulati, Gunsagar S., Shaheen S. Sikandar, Daniel J. Wesche, et al. 2020.
‘Single-Cell Transcriptional Diversity Is a Hallmark of Developmental
Potential’. *Science* 367 (6476): 405–11.
<https://doi.org/10.1126/science.aax0249>.

Guyon, Isabelle, and André Elisseeff. 2003. ‘An Introduction to Variable
and Feature Selection’. *Journal of Machine Learning Research* 3:
1157–82.

Hart, Yuval, Hila Sheftel, Jean Hausser, et al. 2015. ‘Inferring
Biological Tasks Using Pareto Analysis of High-Dimensional Data’.
*Nature Methods* 12 (3): 233–35. <https://doi.org/10.1038/nmeth.3254>.

Jin, Haijing, and Zhandong Liu. 2021. ‘A Benchmark for RNA-seq
Deconvolution Analysis Under Dynamic Testing Environments’. *Genome
Biology* 22. <https://doi.org/10.1186/s13059-021-02290-6>.

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. ‘Moderated
Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2’.
*Genome Biology* 15. <https://doi.org/10.1186/s13059-014-0550-8>.

Lu, Songjian, Jiyuan Yang, Lei Yan, et al. 2025. ‘Transcriptome Size
Matters for Single-Cell RNA-seq Normalization and Bulk Deconvolution’.
*Nature Communications* 16 (1).
<https://doi.org/10.1038/s41467-025-56623-1>.

Melnykov, Volodymyr, Wei-Chen Chen, and Ranjan Maitra. 2024. *MixSim:
Simulating Data to Study Performance of Clustering Algorithms*.

Monaco, Gianni, Bernett Lee, Weili Xu, et al. 2019. ‘RNA-Seq Signatures
Normalized by mRNA Abundance Allow Absolute Deconvolution of Human
Immune Cell Types’. *Cell Reports* 26.
<https://doi.org/10.1016/j.celrep.2019.01.041>.

Newman, Aaron M., Chloé B. Steen, Chih Long Liu, et al. 2019.
‘Determining Cell Type Abundance and Expression from Bulk Tissues with
Digital Cytometry’. *Nature Biotechnology* 37.
<https://doi.org/10.1038/s41587-019-0114-2>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Novoa-del-Toro, Elva María, Efrén Mezura-Montes, Matthieu Vignes, et al.
2021. ‘A Multi-Objective Genetic Algorithm to Find Active Modules in
Multiplex Biological Networks’. *PLOS Computational Biology* 17 (8):
e1009263. <https://doi.org/10.1371/journal.pcbi.1009263>.

Racle, Julien, Kaat de Jonge, Petra Baumgaertner, Daniel E Speiser, and
David Gfeller. 2017. ‘Simultaneous Enumeration of Cancer and Immune Cell
Types from Bulk Tumor Gene Expression Data’. *eLife* 6.
<https://doi.org/10.7554/elife.26476>.

Ranjan, Bobby, Wenjie Sun, Jinyu Park, et al. 2021. ‘DUBStepR Is a
Scalable Correlation-Based Feature Selection Method for Accurately
Clustering Single-Cell Data’. *Nature Communications* 12 (1).
<https://doi.org/10.1038/s41467-021-26085-2>.

Ritchie, Matthew E., Belinda Phipson, Di Wu, et al. 2015. ‘Limma Powers
Differential Expression Analyses for RNA-sequencing and Microarray
Studies’. *Nucleic Acids Research* 43.
<https://doi.org/10.1093/nar/gkv007>.

Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. 2010. ‘edgeR:
A Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data’. *Bioinformatics* 26.
<https://doi.org/10.1093/bioinformatics/btp616>.

Squair, Jordan W., Matthieu Gautier, Claudia Kathe, et al. 2021.
‘Confronting False Discoveries in Single-Cell Differential Expression’.
*Nature Communications* 12 (1).
<https://doi.org/10.1038/s41467-021-25960-2>.

Sturm, Gregor, Francesca Finotello, Florent Petitprez, et al. 2019.
‘Comprehensive Evaluation of Transcriptome-Based Cell-Type
Quantification Methods for Immuno-Oncology’. *Bioinformatics (Oxford,
England)* 35. <https://doi.org/10.1093/bioinformatics/btz363>.

Tsoucas, Daphne, Rui Dong, Haide Chen, Qian Zhu, Guoji Guo, and
Guo-Cheng Yuan. 2019. ‘Accurate Estimation of Cell-Type Composition from
Gene Expression Data’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-019-10802-z>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.

Wang, Yipeng, Xiao-Qin Xia, Zhenyu Jia, et al. 2010. ‘In Silico
Estimates of Tissue Components in Surgical Samples Based on Expression
Profiling Data’. *Cancer Research* 70.
<https://doi.org/10.1158/0008-5472.can-10-0021>.

Zaitsev, Aleksandr, Maksim Chelushkin, Daniiar Dyikanov, et al. 2022.
‘Precise Reconstruction of the TME Using Bulk RNA-seq and a Machine
Learning Algorithm Trained on Artificial Transcriptomes’. *Cancer Cell*
40. <https://doi.org/10.1016/j.ccell.2022.07.006>.

Zhang, Jitao David, Klas Hatje, Gregor Sturm, et al. 2017. ‘Detect
Tissue Heterogeneity in Gene Expression Data with BioQC’. *BMC Genomics*
18. <https://doi.org/10.1186/s12864-017-3661-2>.

Zuo, Yiming, Yi Cui, Cristina Di Poto, et al. 2016. ‘INDEED: Integrated
Differential Expression and Differential Network Analysis of Omic Data
for Biomarker Discovery’. *Methods (San Diego, Calif.)* 111.
<https://doi.org/10.1016/j.ymeth.2016.08.015>.
