# Articles

### External resources

- [External resources: manuscript, presentations, and
  preprints](https://bastienchassagnol.github.io/DeCovarT/articles/decovart-paper-slides.md):

### Getting started

- [Simulating synthetic pseudo-bulk mixtures for
  benchmarking](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.md):
- [Distances, overlap, and covariance
  information](https://bastienchassagnol.github.io/DeCovarT/articles/theory-distance-covariance.md):
- [Feature selection for reference-based
  deconvolution](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S6-feature-selection.md):
- [In silico inference of gene regulatory
  networks](https://bastienchassagnol.github.io/DeCovarT/articles/theory-network-inference.md):

### Paper: numerical simulation study

The two simulation scenarios of the DeCovarT article. **§ 2.1**
Bivariate toy (J = 2, G = 2): factorial sweep of gene–gene correlation,
mean distance, and composition. **§ 2.2** Covariance-driven hybrid (J =
3, G = 20): fixed Gram (cosines 0.9 / 0.1), 2^2 scale-free / cluster-SBM
assignments on the mean-collinear types \times three target MixSim
average overlaps (precision zeros kept), Shannon H^{\star} compositions.
**§ 2.4** Single-cell reference, two lineages (J = 4): sibling cosine
0.8, lineage stochastic block models with scale-free or hub children,
2^5 crossed factors (composition, sibling overlap, cells per bulk, cell
sizes, genes) at platform scale 5, plus one scenario at scale 1.

- [2.1 Bivariate toy model (G = 2 genes, J = 2 cell
  types)](https://bastienchassagnol.github.io/DeCovarT/articles/fig02-bivariate-toy.md):
- [2.2 Covariance-driven scenario (G = 20, J =
  3)](https://bastienchassagnol.github.io/DeCovarT/articles/fig03-covariance-driven.md):
- [2.4 Single-cell scenario: four cell types in two lineages (G = 20 or
  200, J =
  4)](https://bastienchassagnol.github.io/DeCovarT/articles/fig05-single-cell-hierarchical.md):

### Theory: population-level references

Gaussian convolution of population profiles: the generative model, its
derivatives and the two solvers, then the MLE properties and the
confidence intervals.

- [Population-level DeCovarT: generative model, derivatives and
  optimisation on the
  simplex](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-generative-optimisation.md):
- [Population-level DeCovarT inference: MLE properties, confidence
  intervals and a GLM
  formulation](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.md):

### Theory: single-cell references

Mixture of independent cells with two technical scales: reference
moments, alignment strategies and the profile fit, then the profile
information, the intervals and the REML calibration.

- [Single-cell DeCovarT: reference moments, alignment strategies and the
  profile
  fit](https://bastienchassagnol.github.io/DeCovarT/articles/theory-single-cell-generative-optimisation.md):
- [Single-cell DeCovarT inference: profile information, intervals and
  the REML
  calibration](https://bastienchassagnol.github.io/DeCovarT/articles/theory-single-cell-inference.md):

### Mathematics: review and perspectives

- [Deconvolution with reference distributions: a review of methods that
  use more than the
  mean](https://bastienchassagnol.github.io/DeCovarT/articles/theory-review-distribution-aware-deconvolution.md):
- [A Bayesian DeCovarT: latent purified profiles, the Wiener filter and
  variational
  EM](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-bayesian-wiener.md):
- [DeCovarT extensions roadmap: lineages, covariates, modalities and
  beyond](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-extensions-roadmap.md):
- [DeCovarT as a mixture model: links with designed mixture
  experiments](https://bastienchassagnol.github.io/DeCovarT/articles/theory-mixture-regression-applications.md):

### Supplementary simulation appendix

Regular-case MLE checks in the open simplex, identifiability edge cases,
and covariance / scaling / misspecification studies.

- [Appendix S1 — Regular-case MLE properties and
  identifiability](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.md):
- [Appendix S2 — Covariance Inversion
  Strategies](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S2-covariance-inversion.md):
- [Appendix S3 — Solver
  Scalability](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S3-scaling.md):
- [Appendix S4 — Covariance Modelling
  Assumptions](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S4-covariance-modeling.md):
- [Appendix S5 — Model Mis-specification (Description
  Only)](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S5-misspecification.md):
