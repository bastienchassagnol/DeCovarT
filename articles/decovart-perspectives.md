# DeCovarT perspectives

> **Scope**
>
> Outlook on extending DeCovarT beyond bulk transcriptomic
> deconvolution. Future sections (spatial transcriptomics, multi-modal /
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
p_J=\frac{1}{\sum\_{k\<J}e^{\rho_k}+1}. \tag{5}

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

Brown, L., A. N. Donev, and A. C. Bissett. 2015. ‘General Blending
Models for Data from Mixture Experiments’. *Technometrics* 57 (4):
449–56. <https://doi.org/10.1080/00401706.2014.947003>.

Chassagnol, Bastien, Grégory Nuel, and Etienne Becht. 2023. *DeCovarT, a
Multidimensional Probabilistic Model for the Deconvolution of
Heterogeneous Transcriptomic Samples*. arXiv.
<https://doi.org/10.48550/arxiv.2309.09557>.

Morgan-Wall, Tyler, and George Khoury. 2025. *Skpr: Design of
Experiments Suite: Generate and Evaluate Optimal Designs*.
<https://CRAN.R-project.org/package=skpr>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Wheeler, Bob. 2025. *AlgDesign: Algorithmic Experimental Design*.
<https://github.com/jvbraun/AlgDesign>.
