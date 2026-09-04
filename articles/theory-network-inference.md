# In silico inference of gene regulatory networks

> **Note 1: Notation (aligned with the DeCovarT manuscript)**
>
> | Symbol | Role |
> |:---|:---|
> | g=1,\ldots,G | genes |
> | j=1,\ldots,J | cell types |
> | i=1,\ldots,N | purified samples / bootstrap replicates |
> | x\_{gji} | latent expression of gene g in cell type j, replicate i |
> | \mathcal{X}=(x\_{gji}) | tensor \mathcal{M}\_{G\times J\times N} of **single-cell** (purified) profiles |
> | \boldsymbol{\mu}=(\mu\_{gj}) | mean signature \mathcal{M}\_{G\times J}; column \boldsymbol{\mu}\_{\cdot j} is the **aggregated** reference for type j |
> | \boldsymbol{\Sigma}\_j | cell-type covariance \mathcal{M}\_{G\times G}; precision \boldsymbol{\Theta}\_j=\boldsymbol{\Sigma}\_j^{-1} |
> | \boldsymbol{\Omega} | precision of an undirected GGM (distinct from DeCovarT’s cell-type \boldsymbol{\Theta}\_j when discussing GRN inference) |
>
> DeCovarT treats purified draws \boldsymbol{x}\_{\cdot
> j}^{(i)}\sim\mathcal{N}\_G(\boldsymbol{\mu}\_{\cdot j},
> \boldsymbol{\Sigma}\_j) as latent and forms bulk mixtures
> \boldsymbol{y}\_{\cdot i}=\sum_j p_j\\\boldsymbol{x}\_{\cdot j}^{(i)}.
> GRN inference instead targets the **conditional dependence structure**
> among genes within each purified profile—or across a pooled design
> matrix.

## General reviews

Marginal co-expression is a poor proxy for regulation. Shared upstream
TFs, cell-type programmes, and common stimuli induce correlation without
a direct edge; time lags and nonlinear kinetics can hide true edges.
Conditional graphical models target **direct** dependence via the
precision \boldsymbol{\Omega}, while causal / interventional designs are
needed for direction and mechanism ([Gene regulatory network inference
in the era of single-cell
multi-omics](https://doi.org/10.1038/s41576-023-00618-5); Svinin and
Glaab ([2025](#ref-svininCausalNetworkAnalysis2025)); [Network inference
in systems biology](https://doi.org/10.1016/j.copbio.2019.12.002)).

[Badia-i-Mompel *et al.*
(2023)](https://doi.org/10.1038/s41576-023-00618-5) review single-cell
multi-omics GRN tools and prior databases. Svinin and Glaab
([2025](#ref-svininCausalNetworkAnalysis2025)) place association,
prior-informed, and causal network analyses in one taxonomy. Together
they motivate separating undirected skeletons
([Sec. 2](#sec-grn-undirected)) from DAG / perturbation work
([Sec. 3](#sec-grn-directed)).

## Undirected probabilistic strategies

### General inference

Gaussian graphical models encode conditional independence in
\boldsymbol{\Omega}. The MLE exists only under sample-size / chordality
conditions formalised by [Uhler *et al.*
(2012)](https://doi.org/10.1214/11-AOS957): if n is too small relative
to clique size, the covariance MLE fails. Sparse \ell_1 estimation
(graphical lasso; Friedman et al.
([2008](#ref-friedmanSparseInverseCovariance2008))) is the workhorse
when p\gg n ([Chiquet 2015](#ref-chiquetContributionsSparseMethods2015);
[Menéndez et al. 2010](#ref-menendezGeneRegulatoryNetworks2010)).

Two practical workflows appear repeatedly:

1.  **Seeds then expand** — select a small discriminating gene set
    (e.g. low condition number / high centroid separation, in the spirit
    of Aliee and Theis ([2021](#ref-alieeAutoGeneSAutomaticGene2021))),
    then grow a prior-informed network with tools such as
    [NeKo](https://doi.org/10.1371/journal.pcbi.1013300).
2.  **Structure then parameters** — estimate a sparse skeleton (glasso,
    NeighbourNet, CeSpGRN), then fit continuous weights, possibly under
    hard zeros on \Theta\_{ij} for (i,j)\notin E.

Count-aware extensions replace the Gaussian likelihood by
Poisson–log-normal (PLN) or copula models ([Chiquet et al.
2021](#ref-chiquetPoissonLognormalModelVersatile2021),
[2018](#ref-chiquetVariationalInferenceSparse2018)) and
[CeSpGRN](https://doi.org/10.1093/bioinformatics/btag324), which matters
when DeCovarT moves from Gaussian convolutions to sequencing counts.

### Regularisation and topological constraints

Three standard devices:

\max\_{\boldsymbol{\Theta} \succ 0} \log\det\boldsymbol{\Theta} -
\operatorname{tr}(S\boldsymbol{\Theta}) \quad\text{s.t.}\quad
\boldsymbol{\Theta}\_{ij}=0\\\forall\\(i,j)\notin E \tag{1}

(hard covariance selection);

\max\_{\boldsymbol{\Theta} \succ 0} \log\det\boldsymbol{\Theta} -
\operatorname{tr}(S\boldsymbol{\Theta}) - \sum\_{i\neq
j}\lambda\_{ij}\\\|\boldsymbol{\Theta}\_{ij}\| \tag{2}

(edge-specific / hybrid penalties); and structured group or fused
graphical lasso for hubs, modules, or shared biological constraints
([Chiquet 2015](#ref-chiquetContributionsSparseMethods2015)). Trees and
chordal graphs remain PD-friendly ([Uhler *et al.*
(2012)](https://doi.org/10.1214/11-AOS957)). DeCovarT’s affine spectral
shift of a binary adjacency is a constructive PD completion compatible
with hard zeros off E.

### Algorithm unrolling and deep inference

Algorithm unrolling turns iterative sparse-precision solvers into
trainable layers ([Monga *et al.*
(2021)](https://doi.org/10.1109/MSP.2020.3016905);
[GLAD](https://arxiv.org/abs/1906.00271)). For hierarchical count
ecosystems, [Chaussard *et al.*
(2024)](https://hal.science/hal-04622671v1) introduce **PLN-Tree**: a
tree-structured extension of PLN with backward Markov variational
inference and identifiability results—an attractive route when taxonomy
or gene modules should inform \boldsymbol{\Omega}.

### Prior knowledge and PPI

Priors from TF–target, chromatin, and PPI databases can enter as
\lambda\_{ij} masks or as Bayesian edge probabilities ([Badia-i-Mompel
*et al.* (2023)](https://doi.org/10.1038/s41576-023-00618-5); Svinin and
Glaab ([2025](#ref-svininCausalNetworkAnalysis2025));
[NeKo](https://doi.org/10.1371/journal.pcbi.1013300)). This is the
undirected analogue of pathway contextualisation in
[Sec. 3](#sec-grn-directed): association graphs become biologically
constrained without claiming full causality.

## Directed approaches

Gaussian Bayesian networks orient edges under Markov and faithfulness
assumptions; interventions (Perturb-seq and related CRISPR screens)
supply the environments needed for stronger causal claims ([Peters *et
al.* (2016)](https://doi.org/10.1111/rssb.12167); Svinin and Glaab
([2025](#ref-svininCausalNetworkAnalysis2025))). Pathway tools (NPA,
SignalingProfiler, LIONESS/PANDA, NLBayes) contextualise expression on
signed graphs but often stop short of a full weighted GRN MLE—treat them
as priors or downstream scorers rather than drop-in replacements for
\boldsymbol{\Omega}.

> **Note 2: From factorised Gaussian BNs to a joint multivariate
> normal**
>
> Consider a directed acyclic graph on G genes with topological order
> 1,\ldots,G. A **linear Gaussian Bayesian network** specifies, for each
> gene i,
>
> X_i \mid \mathrm{Pa}(i) \sim \mathcal{N}\\\left( \mu_i +
> \sum\_{j\in\mathrm{Pa}(i)} \beta\_{ij} X_j,\\ \sigma_i^{2} \right),
> \tag{3}
>
> where \mathrm{Pa}(i) is the parent set in the DAG and \beta\_{ij}=0
> when j\not\in\mathrm{Pa}(i). Equivalently,
>
> X_i = \mu_i + \sum\_{j=1}^{G} b\_{ij} X_j + \varepsilon_i, \qquad
> \varepsilon_i \sim \mathcal{N}(0,\sigma_i^{2}), \qquad b\_{ij}=0
> \text{ unless } j\in\mathrm{Pa}(i). \tag{4}
>
> Collect coefficients in a strictly upper-triangular matrix B=(b\_{ij})
> (after ordering genes along the DAG) and noise variances
> D=\mathrm{diag}(\sigma_1^{2},\ldots,\sigma_G^{2}). Then
> [Eq. 4](#eq-gbn-structural) is
> \boldsymbol{X}=B\boldsymbol{X}+\boldsymbol{\varepsilon} with
> \boldsymbol{\varepsilon}\sim\mathcal{N}\_G(0,D), hence
>
> (\boldsymbol{I}\_G-B)\boldsymbol{X}=\boldsymbol{\mu}+\boldsymbol{\varepsilon}.
> \tag{5}
>
> When (\boldsymbol{I}\_G-B) is invertible (always true for a DAG with
> consistent ordering),
>
> \boldsymbol{X}\sim\mathcal{N}\_G(\boldsymbol{\mu},\boldsymbol{\Sigma}),
> \qquad \boldsymbol{\Sigma} =(\boldsymbol{I}\_G-B)^{-1} D
> (\boldsymbol{I}\_G-B)^{-\mathsf{T}}. \tag{6}
>
> The joint precision is
>
> \boldsymbol{\Omega}=\boldsymbol{\Sigma}^{-1}
> =(\boldsymbol{I}\_G-B)^{\mathsf{T}} D^{-1}(\boldsymbol{I}\_G-B),
> \tag{7}
>
> which is generally **not** sparse even when the DAG skeleton is
> sparse: moralisation and conditioning introduce fill-in. That is why
> undirected GGMs ([Sec. 2](#sec-grn-undirected)) and directed BNs
> answer different questions: [Eq. 6](#eq-gbn-joint-cov) is the directed
> generative model; [Eq. 1](#eq-hard-cov-selection) targets conditional
> independences in an undirected graph.
>
> In R, `bnlearn::gbn2mvnorm()` ([Scutari 2026](#ref-R-bnlearn)) returns
> (\boldsymbol{\mu}, \boldsymbol{\Sigma}) from a fitted Gaussian BN—the
> directed counterpart of taking
> \boldsymbol{\Sigma}=\boldsymbol{\Omega}^{-1} in an undirected GGM. For
> DeCovarT simulation, the same construction can populate
> cell-type-specific \boldsymbol{\Sigma}\_j when a DAG is the ground
> truth for purified profiles \boldsymbol{x}\_{\cdot j}^{(i)}, before
> bulk convolution with \boldsymbol{p}.

## Benchmarks

Comparative reverse-engineering studies and DREAM-style challenges
remain the default yardsticks for undirected and directed GRN callers
([Penfold & Wild (2011)](https://doi.org/10.1098/rsfs.2011.0053);
[Network inference in systems
biology](https://doi.org/10.1016/j.copbio.2019.12.002); Menéndez et al.
([2010](#ref-menendezGeneRegulatoryNetworks2010))). Emulation tests for
TF–target function complement topological accuracy when gold-standard
DAGs are incomplete. Any DeCovarT-related network estimator should
report against these protocols before biological claims.

## Perspectives

### Multi-omics

Collaborative / multi-layer glasso and tensor-fusion GNNs couple
transcriptomic layers with metabolomic, methylation, or miRNA views.
They generalise the single-precision story to coupled
\boldsymbol{\Omega} blocks—relevant when references cease to be
expression-only.

### Boolean meets Bayesian

Noisy logic models (e.g. NLBayes) combine Boolean pathway structure with
Bayesian TF-activity inference. They are powerful for regulon-level
questions but are not substitutes for a dense weighted precision used
inside DeCovarT’s Gaussian convolution.

### Dynamics

Static snapshots average over regulatory delays. Time-series reviews
([Natalia *et al.*
(2023)](https://doi.org/10.1371/journal.pcbi.1011254)),
[Dictys](https://doi.org/10.1038/s41592-023-01971-3), and related
dynamic biomarker methods target \boldsymbol{\Omega}\_t along
trajectories. Real-time or densely time-stamped assays would further
reduce confounding between co-expression and causation.

### Graph neural networks

OOD-robust GNNs and invariant causal prediction ([Peters *et al.*
(2016)](https://doi.org/10.1111/rssb.12167)) address batch / cell-type
shift when message-passing replaces classical glasso. Unrolling
([Sec. 2.3](#sec-grn-unrolling)) remains preferable when
interpretability of each layer as an optimisation step matters for
statistical GRN recovery.

## Link to DeCovarT simulation

[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
currently draws a random-graph adjacency and completes it to a PD
precision. The literature above suggests three near-term upgrades, in
increasing ambition:

1.  replace the random skeleton by a prior-informed or glasso skeleton
    on real purified profiles;
2.  use edge-specific \lambda\_{ij} (PPI / TF) instead of a uniform
    `precision_scale`;
3.  move from a shared Gaussian \boldsymbol{\Omega} to PLN / PLN-Tree or
    cell-type-specific copula models when count structure dominates.

## References

Aliee, Hananeh, and Fabian J. Theis. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems* 12. <https://doi.org/10.1016/j.cels.2021.05.006>.

Chiquet, Julien. 2015. ‘Contributions to Sparse Methods for Complex Data
Analysis’. PhD thesis, Université Paris-Saclay / AgroParisTech.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2018.
*Variational Inference for Sparse Network Reconstruction from Count
Data*. arXiv. <https://doi.org/10.48550/arxiv.1806.03120>.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2021. ‘The
Poisson-Lognormal Model as a Versatile Framework for the Joint Analysis
of Species Abundances’. *Frontiers in Ecology and Evolution* 9.
<https://doi.org/10.3389/fevo.2021.588292>.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. ‘Sparse
Inverse Covariance Estimation with the Graphical Lasso’. *Biostatistics
(Oxford, England)* 9. <https://doi.org/10.1093/biostatistics/kxm045>.

Menéndez, Patricia, Yiannis A. I. Kourmpetis, Cajo J. F. ter Braak, and
Fred A. van Eeuwijk. 2010. ‘Gene Regulatory Networks from Multifactorial
Perturbations Using Graphical Lasso: Application to the DREAM4
Challenge’. *PLOS ONE* 5.
<https://doi.org/10.1371/journal.pone.0014147>.

Scutari, Marco. 2026. *Bnlearn: Bayesian Network Structure Learning,
Parameter Learning and Inference*. <https://www.bnlearn.com/>.

Svinin, Gleb, and Enrico Glaab. 2025. ‘Causal Network Analysis of Omics
Data Using Prior Knowledge Databases’. *Briefings in Bioinformatics* 26.
<https://doi.org/10.1093/bib/bbaf654>.
