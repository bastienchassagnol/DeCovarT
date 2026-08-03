# In silico inference of gene regulatory networks

## Overview

This vignette curates methodological reading for **Step 1.b** of
DeCovarT: inferring gene regulatory network (GRN) structure that can
inform cell-type precision matrices () (and hence (=^{-1})). The living
checklist lives in [GitHub issue
\#1](https://github.com/bastienchassagnol/DeCovarT/issues/1); here we
summarise the same layout with bibliographic cross-references.

Pandoc `@citekey` citations below resolve only against the Zotero/BBT
export in `inst/REFERENCES.bib`. Papers not yet in that library are
linked by DOI (Markdown) until they are added to Zotero and re-exported.

DeCovarT’s synthetic pipeline already builds a graph-constrained () from
a random adjacency ([synthetic scenarios
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/synthetic-scenarios.md)).
Empirical inference must decide whether that graph is undirected (GGM /
glasso), directed (Gaussian BN), prior-informed, dynamic, or multi-omic.

## i) General reviews — correlation, causation, direct vs indirect effects

Marginal co-expression is a poor proxy for regulation. Shared upstream
TFs, cell-type programmes, and common stimuli induce correlation without
a direct edge; time lags and nonlinear kinetics can hide true edges.
Conditional graphical models target **direct** dependence via the
precision (), while causal / interventional designs are needed for
direction and mechanism ([Gene regulatory network inference in the era
of single-cell multi-omics](https://doi.org/10.1038/s41576-023-00618-5);
Svinin and Glaab ([2025](#ref-svininCausalNetworkAnalysis2025));
[Network inference in systems
biology](https://doi.org/10.1016/j.copbio.2019.12.002)).

[Badia-i-Mompel *et al.*
(2023)](https://doi.org/10.1038/s41576-023-00618-5) review single-cell
multi-omics GRN tools and prior databases. Svinin and Glaab
([2025](#ref-svininCausalNetworkAnalysis2025)) place association,
prior-informed, and causal network analyses in one taxonomy (see figures
in issue \#1). Together they motivate separating undirected skeletons
(**ii**) from DAG / perturbation work (**iii**).

## ii) Undirected probabilistic strategies

### ii.a) General inference

Gaussian graphical models (GGMs) encode conditional independence in ().
The MLE exists only under sample-size / chordality conditions formalised
by [Uhler *et al.* (2012)](https://doi.org/10.1214/11-AOS957): if (n) is
too small relative to clique size, the covariance MLE fails. Sparse
(\_1) estimation (graphical lasso; Friedman et al.
([2008](#ref-friedmanSparseInverseCovariance2008))) is the workhorse
when (pn) ([Chiquet 2015](#ref-chiquetContributionsSparseMethods2015);
[Menéndez et al. 2010](#ref-menendezGeneRegulatoryNetworks2010)).

Two practical workflows appear repeatedly:

1.  **Seeds then expand** — select a small discriminating gene set
    (e.g. low condition number / high centroid separation, in the spirit
    of Aliee and Theis ([2021](#ref-alieeAutoGeneSAutomaticGene2021))),
    then grow a prior-informed network with tools such as
    [NeKo](https://doi.org/10.1371/journal.pcbi.1013300).
2.  **Structure then parameters** — estimate a sparse skeleton (glasso,
    NeighbourNet, CeSpGRN), then fit continuous weights, possibly under
    hard zeros on (\_{ij}) for ((i,j)E).

Count-aware extensions replace the Gaussian likelihood by
Poisson–log-normal (PLN) or copula models
([**chiquetPoissonLognormalModelVersatile2021?**](#ref-chiquetPoissonLognormalModelVersatile2021);
[Chiquet et al. 2018](#ref-chiquetVariationalInferenceSparse2018)), and
[CeSpGRN](https://doi.org/10.1093/bioinformatics/btag324), which matters
when DeCovarT moves from Gaussian convolutions to sequencing counts.

### ii.b) Regularisation and topological constraints

Three standard devices:

\max\_{\Theta \succ 0} \log\det\Theta - \operatorname{tr}(S\Theta)
\quad\text{s.t.}\quad \Theta\_{ij}=0\\\forall\\(i,j)\notin E

(hard covariance selection);

\max\_{\Theta \succ 0} \log\det\Theta - \operatorname{tr}(S\Theta) -
\sum\_{i\neq j}\lambda\_{ij}\\\|\Theta\_{ij}\|

(edge-specific / hybrid penalties); and structured group or fused
graphical lasso for hubs, modules, or shared biological constraints
([Chiquet 2015](#ref-chiquetContributionsSparseMethods2015)). Trees and
chordal graphs remain PD-friendly ([Uhler *et al.*
(2012)](https://doi.org/10.1214/11-AOS957)). DeCovarT’s affine spectral
shift of a binary adjacency is a constructive PD completion compatible
with hard zeros off (E).

### ii.c) Algorithm unrolling and deep / structured variational inference

Algorithm unrolling turns iterative sparse-precision solvers into
trainable layers ([Monga *et al.*
(2021)](https://doi.org/10.1109/MSP.2020.3016905);
[GLAD](https://arxiv.org/abs/1906.00271)). For hierarchical count
ecosystems, [Chaussard *et al.*
(2024)](https://hal.science/hal-04622671v1) introduce **PLN-Tree**: a
tree-structured extension of PLN with backward Markov variational
inference and identifiability results—an attractive route when taxonomy
or gene modules should inform ().

### ii.d) Prior knowledge and PPI as precision weights

Priors from TF–target, chromatin, and PPI databases can enter as
(\_{ij}) masks or as Bayesian edge probabilities ([Badia-i-Mompel *et
al.* (2023)](https://doi.org/10.1038/s41576-023-00618-5); Svinin and
Glaab ([2025](#ref-svininCausalNetworkAnalysis2025));
[NeKo](https://doi.org/10.1371/journal.pcbi.1013300)). This is the
undirected analogue of pathway contextualisation in **iii**: association
graphs become biologically constrained without claiming full causality.

## iii) Directed approaches — DAGs, causality, perturbations

Gaussian Bayesian networks orient edges under Markov and faithfulness
assumptions; interventions (Perturb-seq and related CRISPR screens)
supply the environments needed for stronger causal claims ([Peters *et
al.* (2016)](https://doi.org/10.1111/rssb.12167); Svinin and Glaab
([2025](#ref-svininCausalNetworkAnalysis2025))). Pathway tools (NPA,
SignalingProfiler, LIONESS/PANDA, NLBayes) contextualise expression on
signed graphs but often stop short of a full weighted GRN MLE—treat them
as priors or downstream scorers rather than drop-in replacements for ().

> **From factorised Gaussian BNs to a joint multivariate normal**
>
> For each node, (X_i(i)(i+\*\_{ij}X_j,,\_i^2)). Stacking these linear
> structural equations yields a unique joint (\_G(,)). In R,
> `bnlearn::gbn2mvnorm()` ([Scutari 2026](#ref-R-bnlearn)) returns that
> mean and covariance from a fitted Gaussian BN—the directed counterpart
> of (=^{-1}) in an undirected GGM.

## iv) Benchmarks

Comparative reverse-engineering studies and DREAM-style challenges
remain the default yardsticks for undirected and directed GRN callers
([Penfold & Wild (2011)](https://doi.org/10.1098/rsfs.2011.0053);
[Network inference in systems
biology](https://doi.org/10.1016/j.copbio.2019.12.002); Menéndez et al.
([2010](#ref-menendezGeneRegulatoryNetworks2010))). Emulation tests for
TF–target function complement topological accuracy when gold-standard
DAGs are incomplete. Any DeCovarT-related network estimator should
report against these protocols before biological claims.

## v) Perspectives

### v.i) Multi-omics

Collaborative / multi-layer glasso and tensor-fusion GNNs couple
transcriptomic layers with metabolomic, methylation, or miRNA views.
They generalise the single-precision story to coupled () blocks—relevant
when references cease to be expression-only.

### v.ii) Boolean meets Bayesian

Noisy logic models (e.g. NLBayes) combine Boolean pathway structure with
Bayesian TF-activity inference. They are powerful for regulon-level
questions but are not substitutes for a dense weighted precision used
inside DeCovarT’s Gaussian convolution.

### v.iii) Dynamics and real-time capture

Static snapshots average over regulatory delays. Time-series reviews
([Natalia *et al.*
(2023)](https://doi.org/10.1371/journal.pcbi.1011254)),
[Dictys](https://doi.org/10.1038/s41592-023-01971-3), and related
dynamic biomarker methods target (\_t) along trajectories. Real-time or
densely time-stamped assays would further reduce confounding between
co-expression and causation.

### v.iv) Graph neural networks

OOD-robust GNNs and invariant causal prediction ([Peters *et al.*
(2016)](https://doi.org/10.1111/rssb.12167)) address batch / cell-type
shift when message-passing replaces classical glasso. Unrolling
(**ii.c**) remains preferable when interpretability of each layer as an
optimisation step matters for statistical GRN recovery.

## Link to DeCovarT simulation

[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
currently draws a random-graph adjacency and completes it to a PD
precision. The literature above suggests three near-term upgrades, in
increasing ambition:

1.  replace the random skeleton by a prior-informed or glasso skeleton
    on real purified profiles;
2.  use edge-specific (\_{ij}) (PPI / TF) instead of a uniform
    `precision_scale`;
3.  move from a shared Gaussian () to PLN / PLN-Tree or
    cell-type-specific copula models when count structure dominates.

## References

Aliee, Hananeh, and Fabian J. Theis. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems* 12. <https://doi.org/10.1016/j.cels.2021.05.006>.

Chiquet, Julien. 2015. ‘Contributions to Sparse Methods for Complex Data
Analysis’. PhD thesis.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2018.
*Variational Inference for Sparse Network Reconstruction from Count
Data*. arXiv. <https://doi.org/10.48550/arxiv.1806.03120>.

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
