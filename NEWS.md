# DeCovarT (development version)

* **Scenario 3 (variance-driven).** Fixed Gram
  $\cos(\mu_1,\mu_2)=0.9$, $\cos(\mu_j,\mu_3)=0.1$; $3^3$ graph
  assignments (SBM / Erdős–Rényi / hub) $\times$ three precision
  cushions; compositions from `composition_from_entropy()` at
  $H^{\star}\in\{1,0.5,0.1\}$; solvers LSEI, CIBERSORT, L-BFGS-B,
  Newton–Raphson, Marquardt–Levenberg ($n=50$). Descriptor docs use
  one callout per family in the synthetic-scenarios vignette.
  `kappa_sigma_p` is now returned beside the reciprocal ratio.
  `generate_mean_signature_matrix(nonnegative = TRUE)` uses a
  disjoint-support frame so signatures stay nonnegative. Graph
  precisions are completed by a uniform spectral shift and extra-loaded
  until both $\Omega$ and $\Sigma=\Omega^{-1}$ admit a Cholesky factor
  (support and signs of $W$ unchanged).

* **Terminal UI.** Optional Suggests `cli` formats messages, warnings
  and errors (`cli_alert_*` / `cli_abort`, with base `message()` /
  `stop()` fallback). `run_simulation_benchmark(verbose = TRUE)` logs
  each scenario row; grids of 10 or fewer scenarios also tick every 10
  inferred samples.

* **Documentation layout.** Paper scenarios live in the
  `fig02-bivariate-toy` and `fig03-variance-driven` articles (one script
  each: `scripts/fig02_bivariate_toy.R`,
  `scripts/fig03_variance_driven.R`). Regular-case MLE checks and
  identifiability sit in Appendix S1.   ADEMP / Nature Methods reporting
  is at the end of
  `vignettes/theory-synthetic-scenarios-mean-covariance.qmd`.

* **Monte Carlo figures.** `pivot_mc_estimates()`, `plot_mc_raincloud()`
  (`ggplot2` + `ggdist` horizontal rainclouds of
  $\hat p_j-p_j^{\star}$), `plot_mc_forest()` (ADEMP summaries with
  Wilson whiskers on the coverage rate), `plot_algorithm_similarity()`
  (`geom_tile` + `hclust(1-r)`), and `plot_mc_metric_dots()` (faceted
  min-max scores; no default composite). Optional Suggests: `ggdist`,
  `ggdendro`. `ComplexHeatmap` remains for linked multi-omics grids
  only.

* **Solver starts.** `starting_simplex()` (used by all ILR maps) accepts
  the barycentre (default), a Dirichlet draw (`initial_p = "dirichlet"`,
  `dirichlet_alpha = 1` uniform on the simplex; $\alpha>1$ centre-biased;
  $\alpha<1$ face-biased; several independent draws for multistarts), or
  a mean-only simplex QP (`initial_p = "qp"` via
  `deconvolute_ratios_deconrnaseq()`). `multistart_decovart()` uses the
  Dirichlet path.

* **Fixed-covariance GLS competitor.** `deconvolute_ratios_gls()` wraps
  `MASS::lm.gls()` with a known $G\times G$ residual covariance;
  `fixed_gls_covariance()` builds the global-diagonal $W$ at a declared
  $p$ (default $1/J$). Unconstrained OLS (`deconvolute_ratios_lsfit`)
  has been removed.

* **ILR / Helmert simplex chart.** Frequentist solvers and
  `vcov.decovart_fit()` now use isometric log-ratio coordinates on a
  Helmert basis (`helmert_basis()`, `isometric_logistic()`,
  `isometric_log_ratio()`, `jacobian_isometric_logistic()`,
  `hessian_isometric_logistic()`, `vcov_ilr_delta()`). No cell type is
  pinned as a reference; changing the ILR basis by an orthogonal rotation
  leaves \(\hat{\boldsymbol{p}}\), \(\ell(\hat{\boldsymbol{p}})\) and
  \(\mathrm{Var}(\hat{\boldsymbol{p}})\) unchanged. Additive log-ratio maps
  remain exported for the derivatives vignette appendix and for
  reference-invariance checks (`vcov_alr_delta()`). Restricted MLE
  assembly treats a fully constrained face without calling ILR on an
  empty coordinate, and an unrestricted fit that sits near a face is
  compared with the exact-face restricted MLE so the reported point can
  be the closed-simplex supremum.

* **Gram-matrix mean signatures.** `generate_mean_signature_matrix()`
  realises a target cosine Gram via the symmetric square root
  \(\boldsymbol{\mu}=s\,QR^{1/2}\) (`equicorrelation_gram()` or
  `target_gram`). Pairwise cosines match \(R\) exactly; the previous
  shared-plus-private blend is documented as an appendix construction
  that does not hit \(\rho\) at finite \(J\).

* **Identifiability.** Full column rank of \(\boldsymbol{\mu}\) is
  sufficient but not necessary: identical means remain identifiable when
  the cell-type covariances differ. The article and MLE vignette now
  state injectivity of \(p\mapsto(\mu p,\Sigma(p))\).

* **Benchmark metrics and sample-level parallelism.**
  `compute_benchmark_metrics()` now returns a three-block list:
  composition / regression scores (global TV, RMSE, angular distance,
  SDID, MaxAE; cell-type Pearson, presence F1, and false-positive
  mass), ADEMP Monte Carlo summaries, and optimisation / runtime
  diagnostics (simplex KKT residual, numerical versus theoretical
  convergence, per-sample elapsed time, and `ps` PSS memory).
  Coverage of the Monte Carlo coverage *rate* uses a Wilson interval by
  default (`coverage_mc_interval()`; Wald and Agresti--Coull optional).
  `run_simulation_benchmark()` also returns `theta_true` (the
  convolution parameters), `descriptors` /
  `supplementary` from `describe_simulation_scenario()`, and
  `call` (`match.call()`). MixSim `BarOmega` and averaged pairwise
  Hellinger are kept in `descriptors`; Jeffreys overlap is
  supplementary.
  `deconvolute_ratios()` iterates bulk columns with `furrr` when
  `cores > 1`, using L'Ecuyer-CMRG streams (`furrr_options(seed = TRUE)`).
  Scenario rows are sequential, so workers are never nested. There is
  no composite global score.

* **Structure-aware covariance backends.** New
  `new_decovart_covariance()` wraps
  \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}
  with a declared `structure` (`"dense"`, `"block"`, `"band"`, `"sparse"`,
  `"diag_lowrank"`) and exposes operators `sigma_logdet()`,
  `sigma_solve()`, `sigma_quadform()` and `sigma_trace_precision_times()`
  that return the log-determinant and solves without materialising a dense
  precision. Block Cholesky (disconnected modules / stochastic block
  model), banded Cholesky (band / AR), sparse Cholesky with a cached
  symbolic ordering (Erdős–Rényi, scale-free, small-world), and the
  Woodbury / matrix-determinant-lemma path for diagonal-plus-low-rank
  covariances (shared regulatory programs; `Sigma = W W^T + Psi`) each
  match the dense default to machine precision.
  `covariance_structure_from_graph_model()` maps a network topology to the
  recommended backend, and `.sigma_p_factorisation()` gains an optional
  `backend` argument. Band and sparse backends use the imported
  `Matrix` package. See the new
  "Structure-aware covariance backends" section of
  `vignette("theory-decovart-generative-model")` and
  `scripts/supp_S2_covariance_inversion.R`.

# DeCovarT 2.3.1

* **CRAN resubmission.** Tests no longer write under the package tree.
  The bivariate toy scenario is built in memory by
  `new_bivariate_toy_scenario()` in `tests/testthat/helper.R`; the
  former `fixtures/make-useful-things.R` writer and its golden RDS files
  are gone. CRAN Quarto vignettes do not execute `library(DeCovarT)`
  during rebuild (Windows CLI subprocess cannot see `R CMD build`'s
  temporary library).

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
