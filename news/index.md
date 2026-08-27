# Changelog

## DeCovarT 2.3.0

- **Log-likelihood evaluation.**
  [`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
  now uses the cached Cholesky factor with
  [`backsolve()`](https://rdrr.io/r/base/backsolve.html) for the
  Mahalanobis term, matching `mvtnorm::dmvnorm(..., log = TRUE)` up to
  ((2)), rather than a dense (^({}){-1}). The inverse remains cached for
  the analytic score and Hessian.

- **Corrected log-likelihood (breaking numerical change).**
  [`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
  implemented instead of the Gaussian , doubling the determinant
  contribution. The objective now matches
  `mvtnorm::dmvnorm(..., log = TRUE)` up to , and holds exactly, so the
  objective and the Wald / Fisher inference are finally consistent.
  [`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md)
  and
  [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md)
  halve their determinant terms accordingly; residual terms and
  [`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
  are unchanged. Reported log-likelihood values, AIC and
  likelihood-ratio statistics all change; point estimates move only
  slightly. The article and the derivatives vignette were updated in
  step.

- **Likelihood-ratio and boundary inference.** New
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
  [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
  [`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
  [`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md),
  [`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md)
  and
  [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md).
  Testing the absence of a cell type puts the null on a simplex face,
  where Wilks’ theorem fails and ALR Wald intervals are undefined; the
  null law is the chi-bar-square mixture of Chernoff (1954) and Self and
  Liang (1987), with a restricted parametric bootstrap when the binomial
  weights cannot be justified. Profile intervals are
  reparametrisation-invariant and stay inside .
  [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  resamples the experimental units of a *labelled* reference: donors
  within each cell type (default), cells within each cell type, or
  Dirichlet draws of the composition. Gene-order and cell-type-label
  shuffles of an averaged signature are not inferential procedures;
  [`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
  only verifies that permuting signature columns relabels .

- **Boundary and multimodality diagnostics.** New
  [`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md)
  and
  [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
  wired into `fit_decovart(n_starts =, boundary_tol =)` and stored in
  `fit$diagnostics`. The realised log-likelihood is not globally concave
  (a worked indefinite-Hessian counterexample is in the new vignette),
  so `converged`, `local_maximum`, `near_boundary` and `multimodal` are
  reported as separate claims. A near-boundary estimate is flagged as a
  caution about Wald linearisation, not as optimiser failure.

- **Numerical stability of
  [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md).**
  Now evaluated through the log-sum-exp shift instead of
  `exp(rho) / sum(exp(rho))`, which overflowed to `NaN` for and could
  make non-factorisable when the MLE approached a simplex face.

- **New vignette**
  [`vignette("DeCovarT-MLE-properties")`](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-MLE-properties.md):
  identifiability, existence, failure of uniqueness for one observed
  sample, population consistency through the Kullback–Leibler divergence
  and the Fisher information metric, affine (but not logarithmic)
  equivariance, the four strata of the simplex and their tangent cones,
  and why replication is what makes the asymptotics meaningful.

## DeCovarT 2.2.3

- **Optional heatmap stack.** Moved `ComplexHeatmap`, `circlize`, and
  `viridis` from Imports to Suggests.
  [`plot_correlation_Heatmap()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_correlation_Heatmap.md)
  guards with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and a
  clear install message. Removed `Additional_repositories` (Bioconductor
  is a standard CRAN dependency source; that field was causing the
  availability `? ?` NOTE).

## DeCovarT 2.2.2

- **CRAN resubmission.** Quarto Markdown tables in the two CRAN
  vignettes (no `tinytable` during vignette rebuild on Windows);
  absolute pkgdown URLs for articles excluded from the tarball; Rd
  figure widths in pixels; `cran-comments` clarifies DESCRIPTION
  Hunspell false positives.

## DeCovarT 2.2.0

- **Model-fit wrapper.** Added
  [`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
  an S3 wrapper in the style of standard R model fits (`print`,
  `summary`, `coef`, `fitted`, `residuals`, `vcov`, `nobs`, `confint`,
  `plot`). This makes DeCovarT more versatile and easier for the R
  statistical community to adopt.
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  remains the multi-algorithm benchmark tibble. There is no `formula` /
  [`predict()`](https://rdrr.io/r/stats/predict.html) method: the
  estimator is a variance model, not OLS for bulk expression.

- **[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  solvers.** Marquardt–Levenberg, L-BFGS-B and Newton–Raphson only; ALR
  maps already land on the simplex (no
  [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  on these three).

- **Gene-wise `standardise`.** Affine z-score of bulk, means and
  covariances (equivariant MLE). Log2 mixing (`scaled = TRUE`) is
  rejected (Jensen / CIBERSORT linear space).

- **Wald `vcov`.** Expected Fisher information of the multivariate
  normal convolution, mapped to the simplex by the ALR delta method.

- **Methods paper on the site.** pkgdown embeds `article/main.pdf` via
  the [embedpdf](https://github.com/jmgirard/embedpdf) Quarto extension
  (`vignette` *DeCovarT methods paper*).

## DeCovarT 2.0.0

- **Release** aligned with the GitHub `v2.0.0` tag and the first CRAN
  submission, focused on numerical stability and optimiser cost in
  `R/03_03_DeCovarT_estimate_ratios_frequentist.R` (documented in
  [`vignette("generative-model-derivatives")`](https://bastienchassagnol.github.io/DeCovarT/articles/generative-model-derivatives.md),
  section *Numerical speed-ups and solver safeguards*). The CRAN tarball
  ships the use-cases and softmax/ALR vignettes; remaining articles stay
  on the pkgdown site.

- **Newton–Raphson evaluation budget.** Removed an erroneous
  `eval.max = 1` control passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) in
  [`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md),
  which capped the *entire* run at a single objective evaluation and
  silently returned the equi-balanced start.

- **Cached Cholesky factorisation of
  \boldsymbol{\Sigma}(\boldsymbol{p}).** New internal
  [`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
  assembles the mixture covariance once per trial point, shares
  [`chol()`](https://rdrr.io/r/base/chol.html) /
  [`chol2inv()`](https://rdrr.io/r/base/chol2inv.html) across the
  log-likelihood, gradient and Hessian, and uses a Cholesky-based
  log-determinant (verified against `numDeriv`).

- **Analytic gradient for L-BFGS-B.**
  [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  now passes the unconstrained analytic score instead of finite
  differences, with finite-penalty guards when box-constrained line
  searches drive \boldsymbol{\Sigma}(\boldsymbol{p}) near-singular.

- **`marqLevAlg` Hessian sign under maximisation.** When
  `minimize = FALSE`, `marqLevAlg` only negates `fn` / `gr`; the
  analytic Hessian must be negated by hand so the Commenges et al. RDM
  criterion can fire. Reported upstream as
  [VivianePhilipps/marqLevAlgParallel#3](https://github.com/VivianePhilipps/marqLevAlgParallel/issues/3).

- **Documentation.** Expanded use-cases and feature-selection vignettes;
  CRAN-oriented `DESCRIPTION` (MIT licence, arXiv method reference
  <doi:10.48550/arXiv.2309.09557>).

## DeCovarT 1.0.0

- First **stabilised** release of DeCovarT (semver `1.0.0`), marking the
  end of the exploratory `0.x` series. The core ALR-based frequentist
  deconvolution API, analytic score equations, and simulation tooling
  are now treated as a stable public surface for downstream use and
  paper reproduction.

- **Analytic derivatives, numerically verified.** The unconstrained and
  constrained log-likelihood gradients and Hessians (and the additive
  logistic Jacobian / Hessian) are derived in closed form and
  unit-tested against `numDeriv` finite differences. Earlier faulty
  numerical / finite-difference-only paths that could mislead
  Marquardt–Levenberg and related optimisers have been replaced by these
  verified analytic expressions (see
  [`vignette("generative-model-derivatives")`](https://bastienchassagnol.github.io/DeCovarT/articles/generative-model-derivatives.md)
  and `tests/testthat/test-03_03_DeCovarT.R`).

- **Simulation framework beyond the bivariate toy.** Synthetic first-
  and second-order moments now support arbitrary gene panels (G\gg 2)
  and multi-type mixtures: AutoGeneS-style mean signatures with a cosine
  dial, graph-constrained precisions (ER, hub/star, scale-free, SBM,
  small-world), and bulk convolution via
  [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md).
  A hybrid multi-topology reference scenario (G=50, J=3) in
  `scripts/generate_random_markov_network.R` and
  [`vignette("synthetic-scenarios")`](https://bastienchassagnol.github.io/DeCovarT/articles/synthetic-scenarios.md)
  stresses feature selection under EE / DE / differential-modality-like
  blocks (scDD / `muscat`-inspired taxonomy).

- Feature-selection metrics
  ([`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
  [`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md),
  [`compute_glmnet_gene_scores()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_glmnet_gene_scores.md))
  and a shared
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md)
  validator complete the end-to-end simulation → pre-screen → NSGA-II
  refinement loop documented in
  [`vignette("feature-selection")`](https://bastienchassagnol.github.io/DeCovarT/articles/feature-selection.md).

## DeCovarT 0.1.0

- First official release of DeCovarT: bulk transcriptomic deconvolution
  that accounts for gene-gene covariance in purified reference
  populations.
