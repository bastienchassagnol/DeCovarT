# In silico inference of gene regulatory networks

## Overview

This vignette curates methodological reading for **Step 1.b** of
DeCovarT: inferring gene regulatory network (GRN) structure that can
inform cell-type precision matrices () (and hence (=^{-1})). The living
checklist lives in [GitHub issue
\#1](https://github.com/bastienchassagnol/DeCovarT/issues/1); here we
summarise the same layout with bibliographic cross-references.

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
direction and mechanism ([Badia-i-Mompel et al.
2023](#ref-badia_etal23); [Nuzzo et al. 2025](#ref-nuzzo_etal26);
[Saint-Antoine and Singh 2020](#ref-saint_antoine_singh20)).

Badia-i-Mompel et al. ([2023](#ref-badia_etal23)) review single-cell
multi-omics GRN tools and prior databases. Nuzzo et al.
([2025](#ref-nuzzo_etal26)) place association, prior-informed, and
causal network analyses in one taxonomy (see figures in issue \#1).
Together they motivate separating undirected skeletons (**ii**) from DAG
/ perturbation work (**iii**).

## ii) Undirected probabilistic strategies

### ii.a) General inference

Gaussian graphical models (GGMs) encode conditional independence in ().
The MLE exists only under sample-size / chordality conditions formalised
by Uhler ([2012](#ref-uhler_etal12)): if (n) is too small relative to
clique size, the covariance MLE fails. Sparse (\_1) estimation
(graphical lasso; Friedman et al. ([2008](#ref-friedman_etal08))) is the
workhorse when (pn) ([Chiquet 2015](#ref-chiquet15); [Menéndez et al.
2010](#ref-menendez_etal10)).

Two practical workflows appear repeatedly:

1.  **Seeds then expand** — select a small discriminating gene set
    (e.g. low condition number / high centroid separation, in the spirit
    of Aliee, Hananeh and Theis, Fabian J.
    ([2021](#ref-aliee_theis21))), then grow a prior-informed network
    with tools such as NeKo ([Ruscone et al.
    2025](#ref-ruscone_etal25)).
2.  **Structure then parameters** — estimate a sparse skeleton (glasso,
    NeighbourNet, CeSpGRN), then fit continuous weights, possibly under
    hard zeros on (\_{ij}) for ((i,j)E).

Count-aware extensions replace the Gaussian likelihood by
Poisson–log-normal (PLN) or copula models ([Chiquet et al.
2018](#ref-chiquet_etal19_pln); [Zhang et al.
2025](#ref-zhang_etal25_cespgrn)), which matters when DeCovarT moves
from Gaussian convolutions to sequencing counts.

### ii.b) Regularisation and topological constraints

Three standard devices:

\\ \max\_{\Theta \succ 0} \log\det\Theta - \operatorname{tr}(S\Theta)
\quad\text{s.t.}\quad \Theta\_{ij}=0\\\forall\\(i,j)\notin E \\

(hard covariance selection);

\\ \max\_{\Theta \succ 0} \log\det\Theta - \operatorname{tr}(S\Theta) -
\sum\_{i\neq j}\lambda\_{ij}\\\|\Theta\_{ij}\| \\

(edge-specific / hybrid penalties); and structured group or fused
graphical lasso for hubs, modules, or shared biological constraints
([Chiquet 2015](#ref-chiquet15)). Trees and chordal graphs remain
PD-friendly ([Uhler 2012](#ref-uhler_etal12)). DeCovarT’s affine
spectral shift of a binary adjacency is a constructive PD completion
compatible with hard zeros off (E).

### ii.c) Algorithm unrolling and deep / structured variational inference

Algorithm unrolling turns iterative sparse-precision solvers into
trainable layers ([Monga et al. 2021](#ref-monga_etal21_unrolling);
[Shrivastava et al. 2020](#ref-shrivastava_etal19_glad)). For
hierarchical count ecosystems, Chaussard et al.
([2024](#ref-chaussard_etal24)) introduce **PLN-Tree**: a
tree-structured extension of PLN with backward Markov variational
inference and identifiability results—an attractive route when taxonomy
or gene modules should inform ().

### ii.d) Prior knowledge and PPI as precision weights

Priors from TF–target, chromatin, and PPI databases can enter as
(\_{ij}) masks or as Bayesian edge probabilities ([Badia-i-Mompel et al.
2023](#ref-badia_etal23); [Nuzzo et al. 2025](#ref-nuzzo_etal26);
[Ruscone et al. 2025](#ref-ruscone_etal25)). This is the undirected
analogue of pathway contextualisation in **iii**: association graphs
become biologically constrained without claiming full causality.

## iii) Directed approaches — DAGs, causality, perturbations

Gaussian Bayesian networks orient edges under Markov and faithfulness
assumptions; interventions (Perturb-seq and related CRISPR screens)
supply the environments needed for stronger causal claims ([Peters et
al. 2016](#ref-peters_etal16); [Nuzzo et al. 2025](#ref-nuzzo_etal26)).
Pathway tools (NPA, SignalingProfiler, LIONESS/PANDA, NLBayes)
contextualise expression on signed graphs but often stop short of a full
weighted GRN MLE—treat them as priors or downstream scorers rather than
drop-in replacements for ().

> **From factorised Gaussian BNs to a joint multivariate normal**
>
> For each node, (X_i(i)(*i+*{j(i)}\_{ij}X_j,,\_i^2)). Stacking these
> linear structural equations yields a unique joint (\_G(,)). In R,
> `bnlearn::gbn2mvnorm()` returns that mean and covariance from a fitted
> Gaussian BN ([Scutari 2024](#ref-bnlearn))—the directed counterpart of
> (=^{-1}) in an undirected GGM.

## iv) Benchmarks

Comparative reverse-engineering studies and DREAM-style challenges
remain the default yardsticks for undirected and directed GRN callers
([Penfold and Wild 2011](#ref-penfold_wild11); [Saint-Antoine and Singh
2020](#ref-saint_antoine_singh20); [Menéndez et al.
2010](#ref-menendez_etal10)). Emulation tests for TF–target function
complement topological accuracy when gold-standard DAGs are incomplete.
Any DeCovarT-related network estimator should report against these
protocols before biological claims.

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
([Navarro et al. 2023](#ref-natalia_etal23_ts)), Dictys ([Wang et al.
2023](#ref-wang_etal23_dictys)), and related dynamic biomarker methods
target (\_t) along trajectories. Real-time or densely time-stamped
assays would further reduce confounding between co-expression and
causation.

### v.iv) Graph neural networks

OOD-robust GNNs and invariant causal prediction ([Peters et al.
2016](#ref-peters_etal16)) address batch / cell-type shift when
message-passing replaces classical glasso. Unrolling (**ii.c**) remains
preferable when interpretability of each layer as an optimisation step
matters for statistical GRN recovery.

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

Aliee, Hananeh and Theis, Fabian J. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems*, ahead of print.
<https://doi.org/10.1016/j.cels.2021.05.006>.

Badia-i-Mompel, Pau, Lorna Wessels, Sophia Müller-Dott, et al. 2023.
‘Gene Regulatory Network Inference in the Era of Single-Cell
Multi-Omics’. *Nature Reviews Genetics* 24: 739–54.
<https://doi.org/10.1038/s41576-023-00618-5>.

Chaussard, Alexandre, Anna Bonnet, Élisabeth Gassiat, and Sylvain Le
Corff. 2024. ‘Tree-Based Variational Inference for Poisson Log-Normal
Models’. <https://hal.science/hal-04622671v1>.

Chiquet, Julien. 2015. ‘Contributions to Sparse Methods for Complex Data
Analysis’. PhD thesis, Université d’Évry-Val-d’Essonne.
<https://hal.science/tel-01288976/>.

Chiquet, Julien, Mahendra Mariadassou, and Stéphane Robin. 2018.
‘Variational Inference for Probabilistic Poisson PCA’. *The Annals of
Applied Statistics* 12 (4): 2674–98.
<https://doi.org/10.1214/18-AOAS1177>.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. ‘Sparse
Inverse Covariance Estimation with the Graphical Lasso’. *Biostatistics*
9 (3): 432–41. <https://doi.org/10.1093/biostatistics/kxm045>.

Menéndez, Patricia, Yiannis A. I. Kourmpetis, Cajo J. F. ter Braak, and
Fred A. van Eeuwijk. 2010. ‘Gene Regulatory Networks from Multifactorial
Perturbations Using Graphical Lasso’. *PLOS ONE* 5 (12): e14147.
<https://doi.org/10.1371/journal.pone.0014147>.

Monga, Vishal, Yuelong Li, and Yonina C. Eldar. 2021. ‘Algorithm
Unrolling: Interpretable, Efficient Deep Learning for Signal and Image
Processing’. *IEEE Signal Processing Magazine* 38 (2): 18–44.
<https://doi.org/10.1109/MSP.2020.3016905>.

Navarro, Jeremie et al. 2023. ‘From Time-Series Transcriptomics to Gene
Regulatory Networks: A Review on Inference Methods’. *PLOS Computational
Biology*, ahead of print.
<https://doi.org/10.1371/journal.pcbi.1011254>.

Nuzzo, Angelo et al. 2025. ‘Causal Network Analysis of Omics Data Using
Prior Knowledge Databases’. *Briefings in Bioinformatics* 26 (6).
<https://doi.org/10.1093/bib/bbaf654>.

Penfold, Christopher A., and David L. Wild. 2011. ‘How to Infer Gene
Networks from Expression Profiles’. *Interface Focus*, ahead of print.
<https://doi.org/10.1098/rsfs.2011.0053>.

Peters, Jonas, Peter Bühlmann, and Nicolai Meinshausen. 2016. ‘Causal
Inference by Using Invariant Prediction: Identification and Confidence
Intervals’. *Journal of the Royal Statistical Society: Series B* 78 (5):
947–1012. <https://doi.org/10.1111/rssb.12167>.

Ruscone, Marco, Eirini Tsirvouli, Andrea Checcoli, et al. 2025. ‘NeKo: A
Tool for Automatic Network Construction from Prior Knowledge’. *PLOS
Computational Biology*, ahead of print.
<https://doi.org/10.1371/journal.pcbi.1013300>.

Saint-Antoine, Michael M., and Abhyudai Singh. 2020. ‘Network Inference
in Systems Biology: Recent Developments, Challenges, and Applications’.
*Current Opinion in Biotechnology* 63: 89–98.
<https://doi.org/10.1016/j.copbio.2019.12.002>.

Scutari, Marco. 2024. *bnlearn: Bayesian Network Structure Learning,
Parameter Learning and Inference*. <https://www.bnlearn.com/>.

Shrivastava, Harsh et al. 2020. ‘GLAD: Learning Sparse Graph Recovery’.
*International Conference on Learning Representations*.
<https://arxiv.org/abs/1906.00271>.

Uhler, Caroline. 2012. ‘Geometry of Maximum Likelihood Estimation in
Gaussian Graphical Models’. *The Annals of Statistics* 40 (1): 238–61.
<https://doi.org/10.1214/11-AOS957>.

Wang, Lingfei et al. 2023. ‘Dictys: Dynamic Gene Regulatory Network
Dissects Developmental Continuum with Single-Cell Multiomics’. *Nature
Methods* 20: 1368–78. <https://doi.org/10.1038/s41592-023-01971-3>.

Zhang, Shiyu et al. 2025. ‘CeSpGRN: Inferring Cell-Specific Gene
Regulatory Networks from Single-Cell Multi-Omics and Spatial Data’.
*Bioinformatics*, ahead of print.
<https://doi.org/10.1093/bioinformatics/btag324>.
