# DeCovarT as a mixture model: links with designed mixture experiments

> **Scope**
>
> Bulk deconvolution is one instance of a problem that appears wherever
> a measured signal is a blend of known ingredients. This vignette
> places DeCovarT in the vocabulary of **mixture experiments**, where
> the proportions are regressors on the simplex, and lists the modelling
> devices of that field (Scheffé polynomials, interaction terms,
> designed mixture points) that could be borrowed. Extensions specific
> to transcriptomics are in the [extensions
> roadmap](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-extensions-roadmap.md).

Write g=1,\ldots,G for the measured features (genes), j=1,\ldots,J for
the components (cell types), \boldsymbol{p}\in\Delta^{J-1} for the
proportions, \boldsymbol{\mu}\_j and \boldsymbol{\Sigma}\_j for the mean
and covariance of component j.

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

## Mixture design and first-order structure

> **Mixture experiments and compositional data**
>
> DeCovarT—and bulk deconvolution more generally—sits inside the wider
> class of **mixture inverse problems** in formulation science,
> spectroscopy, and designed experiments where component fractions are
> explicit regressors on the simplex. R offers only a handful of
> packages for **designed mixture points** (as opposed to random
> \boldsymbol{p} on \Delta^{J-1}):
>
> - **`AlgDesign`** ([Wheeler 2025](#ref-R-AlgDesign)) — lattice and
>   simplex-centroid mixture designs via `genMixture()` /
>   `gen.mixture()`, plus D-, A-, and I-optimal designs; the closest
>   general substitute for classical mixture DOE.
> - **`skpr`** ([Morgan-Wall and Khoury 2025](#ref-R-skpr)) — optimal
>   mixture designs from a candidate set under \sum_j p_j=1; preferable
>   when optimality criteria matter more than fixed simplex-centroid
>   lattices.
>
> A modular Scheffé-type benchmark would combine `AlgDesign` (or `skpr`)
> to propose \boldsymbol{p}, DeCovarT (or baselines) to invert the
> mixture, and compositional metrics from the simulation vignettes.

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

The present generative framework keeps a **linear** mean map
\boldsymbol{y}=\boldsymbol{\mu}\boldsymbol{p}. When the regression of a
compositional response on Euclidean predictors is itself nonlinear,
flexible simplicial nonparametric models ([Tsagris et al.
2023](#ref-tsagrisFlexibleNonparametricRegression2023)) and folded
\alpha-distributions on the simplex ([Tsagris and Stewart
2020](#ref-tsagrisFoldedModelCompositional2020)) are a natural outlook;
they are not used in the current convolution.

## Other mixture inverse problems

The same structure, a signal equal to a proportion-weighted sum of
component profiles, recurs outside genomics. Three settings are close
enough to DeCovarT that their methods are worth importing, and their
literature should be added to the Zotero library before it is cited
here.

| Field | Signal | Components | Second-order information |
|----|----|----|----|
| Designed mixture experiments | a measured property of a blend | ingredients with fixed fractions | replicate variance of each ingredient’s response |
| Spectral unmixing (chemometrics, remote sensing) | a spectrum or a pixel | pure-component spectra (“endmembers”) | variability of each endmember across samples |
| Receptor-based source apportionment (atmospheric chemistry) | chemical species in a filter sample | source profiles | measurement uncertainty of each source profile |

Table 1: Mixture inverse problems with the same algebra as bulk
deconvolution. References for the last two rows are to be imported into
Zotero.

The common lesson is the one DeCovarT formalises: when the component
profiles are themselves uncertain, a weighted least-squares fit should
weight each feature by the variance the profiles induce at the current
proportions, \sum_jp_j^{2}\boldsymbol{\Sigma}\_j for independent
components, rather than by a fixed residual variance.

## References

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

Tsagris, Michail, Abdulaziz Alenazi, and Connie Stewart. 2023. ‘Flexible
Non-Parametric Regression Models for Compositional Response Data with
Zeros’. *Statistics and Computing* 33 (5): 106.
<https://doi.org/10.1007/s11222-023-10277-5>.

Tsagris, Michail, and Connie Stewart. 2020. ‘A Folded Model for
Compositional Data Analysis’. *Australian & New Zealand Journal of
Statistics* 62 (2): 249–77. <https://doi.org/10.1111/anzs.12289>.

Wheeler, Bob. 2025. *AlgDesign: Algorithmic Experimental Design*.
<https://github.com/jvbraun/AlgDesign>.
