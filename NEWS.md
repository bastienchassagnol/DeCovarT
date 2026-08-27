# DeCovarT 2.3.0

* **Log-likelihood evaluation.** `loglik_multivariate()` now uses the
  cached Cholesky factor with `backsolve()` for the Mahalanobis term,
  matching `mvtnorm::dmvnorm(..., log = TRUE)` up to
  \(\tfrac{G}{2}\log(2\pi)\), rather than a dense
  \(\boldsymbol{r}^{\mathsf{T}}\boldsymbol{\Sigma}^{-1}\boldsymbol{r}\).
  The inverse remains cached for the analytic score and Hessian.

* **Corrected log-likelihood (breaking numerical change).**
  `loglik_multivariate()` implemented
  \eqn{-\log\det\boldsymbol{\Sigma}(\boldsymbol{p})} instead of the
  Gaussian \eqn{-\tfrac{1}{2}\log\det\boldsymbol{\Sigma}(\boldsymbol{p})},
  doubling the determinant contribution. The objective now matches
  `mvtnorm::dmvnorm(..., log = TRUE)` up to \eqn{\tfrac{G}{2}\log(2\pi)},
  and \eqn{\mathbb{E}[-\mathbf{H}]=I(\boldsymbol{p})} holds exactly, so
  the objective and the Wald / Fisher inference are finally consistent.
  `gradient_loglik_unconstrained()` and
  `hessian_loglik_unconstrained()` halve their determinant terms
  accordingly; residual terms and `expected_fisher_unconstrained()` are
  unchanged. Reported log-likelihood values, AIC and likelihood-ratio
  statistics all change; point estimates move only slightly. The article
  and the derivatives vignette were updated in step.

* **Likelihood-ratio and boundary inference.** New
  `lrt_decovart()`, `confint_profile_decovart()`,
  `profile_loglik_decovart()`, `restricted_mle_decovart()`,
  `chi_bar_square_pvalue()` and `bootstrap_decovart()`. Testing the
  absence of a cell type puts the null on a simplex face, where Wilks'
  theorem fails and ALR Wald intervals are undefined; the null law is the
  chi-bar-square mixture of Chernoff (1954) and Self and Liang (1987),
  with a restricted parametric bootstrap when the binomial weights cannot
  be justified. Profile intervals are reparametrisation-invariant and stay
  inside \eqn{[0,1]}. `reference_bootstrap_decovart()` resamples the
  experimental units of a *labelled* reference: donors within each cell
  type (default), cells within each cell type, or Dirichlet draws of the
  composition. Gene-order and cell-type-label shuffles of an averaged
  signature are not inferential procedures;
  `equivariance_check_decovart()` only verifies that permuting signature
  columns relabels \eqn{\hat{\boldsymbol{p}}}.

* **Boundary and multimodality diagnostics.** New
  `boundary_diagnostics()` and `multistart_decovart()`, wired into
  `fit_decovart(n_starts =, boundary_tol =)` and stored in
  `fit$diagnostics`. The realised log-likelihood is not globally concave
  (a worked indefinite-Hessian counterexample is in the new vignette), so
  `converged`, `local_maximum`, `near_boundary` and `multimodal` are
  reported as separate claims. A near-boundary estimate is flagged as a
  caution about Wald linearisation, not as optimiser failure.

* **Numerical stability of `additive_logistic()`.** Now evaluated through
  the log-sum-exp shift instead of `exp(rho) / sum(exp(rho))`, which
  overflowed to `NaN` for \eqn{\rho_i\gtrsim 710} and could make
  \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} non-factorisable when the MLE
  approached a simplex face.

* **New vignette** `vignette("DeCovarT-MLE-properties")`: identifiability,
  existence, failure of uniqueness for one observed sample, population
  consistency through the Kullback--Leibler divergence and the Fisher
  information metric, affine (but not logarithmic) equivariance, the four
  strata of the simplex and their tangent cones, and why replication is
  what makes the asymptotics meaningful.

# DeCovarT 2.2.3

* **Optional heatmap stack.** Moved `ComplexHeatmap`, `circlize`, and
  `viridis` from Imports to Suggests. `plot_correlation_Heatmap()`
  guards with `requireNamespace()` and a clear install message.
  Removed `Additional_repositories` (Bioconductor is a standard CRAN
  dependency source; that field was causing the availability `? ?` NOTE).

# DeCovarT 2.2.2

* **CRAN resubmission.** Quarto Markdown tables in the two CRAN vignettes
  (no `tinytable` during vignette rebuild on Windows); absolute pkgdown
  URLs for articles excluded from the tarball; Rd figure widths in
  pixels; `cran-comments` clarifies DESCRIPTION Hunspell false positives.

# DeCovarT 2.2.0

* **Model-fit wrapper.** Added `fit_decovart()`, an S3 wrapper in the
  style of standard R model fits (`print`, `summary`, `coef`, `fitted`,
  `residuals`, `vcov`, `nobs`, `confint`, `plot`). This makes DeCovarT
  more versatile and easier for the R statistical community to adopt.
  `deconvolute_ratios()` remains the multi-algorithm benchmark tibble.
  There is no `formula` / `predict()` method: the estimator is a
  variance model, not OLS for bulk expression.

* **`fit_decovart()` solvers.** Marquardt--Levenberg, L-BFGS-B and
  Newton--Raphson only; ALR maps already land on the simplex (no
  `repair_simplex()` on these three).

* **Gene-wise `standardise`.** Affine z-score of bulk, means and
  covariances (equivariant MLE). Log2 mixing (`scaled = TRUE`) is
  rejected (Jensen / CIBERSORT linear space).

* **Wald `vcov`.** Expected Fisher information of the multivariate
  normal convolution, mapped to the simplex by the ALR delta method.

* **Methods paper on the site.** pkgdown embeds `article/main.pdf` via
  the [embedpdf](https://github.com/jmgirard/embedpdf) Quarto extension
  (`vignette` *DeCovarT methods paper*).

# DeCovarT 2.0.0

* **Release** aligned with the GitHub `v2.0.0` tag and the first CRAN
  submission, focused on numerical stability and optimiser cost in
  `R/03_03_DeCovarT_estimate_ratios_frequentist.R` (documented in
  `vignette("generative-model-derivatives")`, section *Numerical speed-ups and
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
  expressions (see `vignette("generative-model-derivatives")` and
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
