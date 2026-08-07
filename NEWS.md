# DeCovarT 2.0.0

* **Release** aligned with the GitHub `v2.0.0` tag and the first CRAN
  submission, focused on numerical stability and optimiser cost in
  `R/03_03_DeCovarT_estimate_ratios_frequentist.R` (documented in
  `vignette("softmax-alr-derivatives")`, section *Numerical speed-ups and
  solver safeguards*). The CRAN tarball ships the use-cases and
  softmax/ALR vignettes; remaining articles stay on the pkgdown site.

* **Newton–Raphson evaluation budget.** Removed an erroneous
  `eval.max = 1` control passed to `stats::nlminb()` in
  `deconvolute_ratios_Newton_Raphson()`, which capped the *entire* run at a
  single objective evaluation and silently returned the equi-balanced start.

* **Cached Cholesky factorisation of $\boldsymbol{\Sigma}(\boldsymbol{p})$.**
  New internal `.sigma_p_factorisation()` assembles the mixture covariance
  once per trial point, shares `chol()` / `chol2inv()` across the
  log-likelihood, gradient and Hessian, and uses a Cholesky-based
  log-determinant (verified against `numDeriv`).

* **Analytic gradient for L-BFGS-B.** `deconvolute_ratios_L_BFGS_B()` now
  passes the unconstrained analytic score instead of finite differences,
  with finite-penalty guards when box-constrained line searches drive
  $\boldsymbol{\Sigma}(\boldsymbol{p})$ near-singular.

* **`marqLevAlg` Hessian sign under maximisation.** When
  `minimize = FALSE`, `marqLevAlg` only negates `fn` / `gr`; the analytic
  Hessian must be negated by hand so the Commenges et al. RDM criterion
  can fire. Reported upstream as
  [VivianePhilipps/marqLevAlgParallel#3](https://github.com/VivianePhilipps/marqLevAlgParallel/issues/3).

* **Documentation.** Expanded use-cases and feature-selection vignettes;
  CRAN-oriented `DESCRIPTION` (MIT licence, arXiv method reference
  <doi:10.48550/arXiv.2309.09557>).

# DeCovarT 1.0.0

* First **stabilised** release of DeCovarT (semver `1.0.0`), marking the
  end of the exploratory `0.x` series. The core ALR-based frequentist
  deconvolution API, analytic score equations, and simulation tooling are
  now treated as a stable public surface for downstream use and paper
  reproduction.

* **Analytic derivatives, numerically verified.** The unconstrained and
  constrained log-likelihood gradients and Hessians (and the additive
  logistic Jacobian / Hessian) are derived in closed form and unit-tested
  against `numDeriv` finite differences. Earlier faulty numerical /
  finite-difference-only paths that could mislead Marquardt–Levenberg and
  related optimisers have been replaced by these verified analytic
  expressions (see `vignette("softmax-alr-derivatives")` and
  `tests/testthat/test-03_03_DeCovarT.R`).

* **Simulation framework beyond the bivariate toy.** Synthetic
  first- and second-order moments now support arbitrary gene panels
  ($G\gg 2$) and multi-type mixtures: AutoGeneS-style mean signatures
  with a cosine dial, graph-constrained precisions (ER, hub/star,
  scale-free, SBM, small-world), and bulk convolution via
  `simulate_bulk_mixture()`. A hybrid multi-topology reference scenario
  ($G=50$, $J=3$) in `scripts/generate_random_markov_network.R` and
  `vignette("synthetic-scenarios")` stresses feature selection under
  EE / DE / differential-modality-like blocks (scDD / `muscat`-inspired
  taxonomy).

* Feature-selection metrics (`compute_average_overlap()`,
  `compute_average_jeffreys()`, `compute_glmnet_gene_scores()`) and a
  shared `check_true_theta()` validator complete the end-to-end
  simulation → pre-screen → NSGA-II refinement loop documented in
  `vignette("feature-selection")`.

# DeCovarT 0.1.0

* First official release of DeCovarT: bulk transcriptomic deconvolution
  that accounts for gene-gene covariance in purified reference
  populations.
