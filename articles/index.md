# Articles

### External resources

- [External resources: manuscript, presentations, and
  preprints](https://bastienchassagnol.github.io/DeCovarT/articles/decovart-paper-slides.md):

### Theory 1 · Population-level references

Gaussian convolution of population profiles: the generative model, its
derivatives and the two solvers, then the MLE properties and the
confidence intervals.

- [1.1 Population-level DeCovarT: generative model, derivatives and
  optimisation on the
  simplex](https://bastienchassagnol.github.io/DeCovarT/articles/01-1-population-generative-optimisation.md):
- [1.2 Population-level DeCovarT inference: MLE properties, confidence
  intervals and a GLM
  formulation](https://bastienchassagnol.github.io/DeCovarT/articles/01-2-population-inference.md):

### Theory 2 · Single-cell references

Mixture of independent cells with two technical scales: reference
moments, alignment strategies and the profile fit, then the profile
information, the intervals and the REML calibration.

- [2.1 Single-cell DeCovarT: reference moments, alignment strategies and
  the profile
  fit](https://bastienchassagnol.github.io/DeCovarT/articles/02-1-single-cell-generative-optimisation.md):
- [2.2 Single-cell DeCovarT inference: profile information, intervals
  and the REML
  calibration](https://bastienchassagnol.github.io/DeCovarT/articles/02-2-single-cell-inference.md):

### Theory 3 · Estimating the reference: genes and covariances

Which genes enter the reference, and how each cell type’s covariance is
estimated before it is plugged into the likelihood.

- [3.1 Gene selection for the
  reference](https://bastienchassagnol.github.io/DeCovarT/articles/03-1-feature-selection.md):
- [3.2 Cell-type covariance for the DeCovarT
  plug-in](https://bastienchassagnol.github.io/DeCovarT/articles/03-2-covariance-network-estimation.md):
- [3.3 Precision networks: alternatives to the covariance
  plug-in](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.md):

### Theory 4 · Benchmark design

Synthetic mean and covariance scenarios (ADEMP), and the distances and
overlap scores that describe them.

- [4.1 Simulating synthetic pseudo-bulk mixtures for
  benchmarking](https://bastienchassagnol.github.io/DeCovarT/articles/04-1-synthetic-scenarios.md):
- [4.2 Distances, overlap, and covariance
  information](https://bastienchassagnol.github.io/DeCovarT/articles/04-2-distances-overlap-metrics.md):

### Theory 5 · Review and perspectives

- [5.1 Deconvolution with reference distributions: a review of methods
  that use more than the
  mean](https://bastienchassagnol.github.io/DeCovarT/articles/05-1-review-distribution-aware.md):
- [5.2 A Bayesian DeCovarT: latent purified profiles, the Wiener filter
  and variational
  EM](https://bastienchassagnol.github.io/DeCovarT/articles/05-2-bayesian-wiener.md):
- [5.3 DeCovarT extensions roadmap: lineages, covariates, modalities and
  beyond](https://bastienchassagnol.github.io/DeCovarT/articles/05-3-extensions-roadmap.md):
- [5.4 DeCovarT as a mixture model: links with designed mixture
  experiments](https://bastienchassagnol.github.io/DeCovarT/articles/05-4-mixture-experiments.md):

### Simulation study (article §2)

The three simulation scenarios of the DeCovarT article. **§2.2**
Bivariate toy (J = 2, G = 2): factorial sweep of gene–gene correlation,
mean distance, and composition. **§2.3** Covariance-driven hybrid (J =
3, G = 20): fixed Gram (cosines 0.9 / 0.1), 2^2 scale-free / cluster-SBM
assignments on the mean-collinear types \times three target MixSim
average overlaps (precision zeros kept), Shannon H^{\star} compositions.
**§2.4** Single-cell reference, two lineages (J = 4): sibling cosine
0.8, lineage stochastic block models with scale-free or hub children,
2^6 = 64 crossed scenarios (composition, sibling overlap, cells per
bulk, cell sizes, genes, platform scale 1 or 5).

- [§2.2 Bivariate toy model (G = 2 genes, J = 2 cell
  types)](https://bastienchassagnol.github.io/DeCovarT/articles/R1-bivariate-toy.md):
- [§2.3 Covariance-driven scenario (G = 20, J =
  3)](https://bastienchassagnol.github.io/DeCovarT/articles/R2-covariance-driven.md):
- [§2.4 Single-cell scenario: four cell types in two lineages (G = 20 or
  200, J =
  4)](https://bastienchassagnol.github.io/DeCovarT/articles/R3-single-cell-hierarchical.md):

### Supplementary appendix

Regular-case MLE checks in the open simplex, identifiability edge cases,
covariance / scaling / misspecification studies, and a synthetic
scenario for comparing gene-screening strategies.

- [Appendix S1 — Regular-case MLE properties and
  identifiability](https://bastienchassagnol.github.io/DeCovarT/articles/S1-identifiability.md):
- [Appendix S2 — Covariance Inversion
  Strategies](https://bastienchassagnol.github.io/DeCovarT/articles/S2-covariance-inversion.md):
- [Appendix S3 — Solver
  Scalability](https://bastienchassagnol.github.io/DeCovarT/articles/S3-scaling.md):
- [Appendix S4 — Covariance Modelling
  Assumptions](https://bastienchassagnol.github.io/DeCovarT/articles/S4-covariance-modeling.md):
- [Appendix S5 — Model Mis-specification (Description
  Only)](https://bastienchassagnol.github.io/DeCovarT/articles/S5-misspecification.md):
- [Appendix S6 — Synthetic scenario for comparing gene-screening
  strategies](https://bastienchassagnol.github.io/DeCovarT/articles/S6-feature-screening-scenario.md):
