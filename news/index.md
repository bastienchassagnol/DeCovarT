# Changelog

## DeCovarT 2.5.0

### Breaking changes

- **Likelihood objects replace the `model` string**, in the manner of
  `glm(family = binomial)`.
  [`population_likelihood()`](https://bastienchassagnol.github.io/DeCovarT/reference/population_likelihood.md)
  and
  [`single_cell_likelihood()`](https://bastienchassagnol.github.io/DeCovarT/reference/population_likelihood.md)
  return an object of class `decovart_likelihood` holding the simulator,
  log-likelihood, score, Wald standard errors, curvature check,
  parameter rows and interval menu of each model.
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  and
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
  take `likelihood =` (the object or its constructor) in place of
  `model = "population" | "single_cell"`.
- `benchmark_likelihood()`, exported in 2.4.0, is removed; its functions
  are fields of the likelihood objects and take `(p, y, mu, Sigma)`.

### Boundary estimates

- L-BFGS-B works on the closed orthant and returns exact zeros: one zero
  proportion (a face), two (an edge), up to J - 1 (a vertex, one cell
  type making up the whole sample). Newton-Raphson in ILR coordinates
  never does. At any such estimate the Wald, ILR and oracle intervals
  are undefined and are `NA`, the scale rows of the single-cell model
  included.
- The figure 5 crash at G = 200 was a shape error, not a failure of an
  interval derivation: on a boundary estimate the empty single-cell
  interval had the J composition rows only, while interior ones also had
  the four scale rows. Every bulk now has the rows of
  `likelihood$parameters()`.
- The optimisation table gains `n_zero` (zero proportions per bulk) and
  the interval table `n_boundary` (bulks with a zero proportion): Wald
  and ILR coverage is over interior estimates only.
- The population simulator now draws from `generative_mean` and
  `generative_sigma` when a scenario sets them, as the single-cell one
  already did; population scenarios do not set them, so their bulks are
  unchanged.
- Profile intervals refit from the bulk and do not depend on the solver;
  the figure scripts now request them under Newton-Raphson only.

## DeCovarT 2.4.0

- **Benchmark interval tables tolerate a missing parameter row.** When
  one bulk’s interval omitted a parameter (seen with L-BFGS-B at G = 200
  in figure 5),
  [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  stopped with a subscript error. That bulk now counts as having no
  interval for the parameter. A sample whose optimiser returns no finite
  estimate is counted separately in `n_numerical_failures`, and its
  convergence record is `"numerical failure"`.

- **Exported for scripts, so none needs `:::` or
  [`get()`](https://rdrr.io/r/base/get.html).**
  [`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md),
  [`hellinger_gaussian()`](https://bastienchassagnol.github.io/DeCovarT/reference/hellinger_gaussian.md)
  (formerly `.hellinger_gaussian()`),
  [`kkt_residual()`](https://bastienchassagnol.github.io/DeCovarT/reference/kkt_residual.md)
  (formerly `.kkt_residual()`) and `benchmark_likelihood()` (formerly
  `.benchmark_model()`).

- **No package-level mutable state.** `cli` and `withr` move to Imports.
  The plain-text fallback of the `.ui_*` wrappers is removed and the
  package calls `cli` directly. The single-slot Cholesky cache of the
  covariance mixture is removed. Estimates are unchanged. Newton fits
  take 1.1 to 1.4 times as long at G = 20 (about 0.1 s per 200 fits); at
  G = 200 repeated timings ranged from 0.6 to 1.5 times, because the
  identity check on the covariance array costs about as much as the
  factorisation it saves. The qrng-missing note is shown once through
  `cli`’s own frequency control; its
  [`runif()`](https://rdrr.io/r/stats/Uniform.html) fallback no longer
  changes the caller’s random state. `.with_seed()` delegates to
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html).

- New test `test-linear-algebra.R` fails when BLAS / LAPACK return wrong
  Cholesky factors or inverses. This happens with OpenBLAS 0.3.20 on
  Intel Xeon w9-3575X (set `OPENBLAS_CORETYPE=SkylakeX` before R
  starts); it was the cause of the fig05 failure at G = 200.

- **Simulation benchmarks rebuilt (figures 2, 3 and the new figure 5).**

  - [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
    now runs L-BFGS-B on the closure chart p = w / sum(w), w \>= 0, and
    maximises the convolution likelihood itself: box KKT points are
    simplex KKT points and exact zeros are reachable. The former fit
    renormalised a box-constrained surrogate, whose faces were not
    stationary (KKT residual about 0.1). A stick-breaking box was tried
    first and rejected: once a stick is exhausted the later coordinates
    have no gradient, and the solver stopped at vertices up to 2.19
    log-likelihood units below the maximum (0.0008 at most for the
    closure chart).
  - `deconvolute_single_cell(optimiser = "lbfgsb")` maximises the
    profile log-likelihood on the same chart;
    [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)
    and
    [`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md)
    gain `optimiser`, so bootstrap refits use the optimiser of the
    estimate.
  - [`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md):
    the cell size now scales the covariance quadratically, S_j^2
    Sigma_j^theta, as it scales every gene of a cell; the overlap of two
    types of equal size no longer depends on that size.
  - `scale_covariances_to_overlap(pairs = )` calibrates the mean overlap
    of selected pairs of types (siblings of a lineage).
  - New
    [`simulate_reference_replicates()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_reference_replicates.md):
    an independent reference panel of individuals (and cells) per type,
    for benchmark methods that estimate the reference moments
    themselves.
  - `run_simulation_benchmark(model = "single_cell", return_data = TRUE, reference = )`
    simulates bulks from cells, scores the realised fractions and the
    two scales, and returns the bulks and a reference panel. Each
    deconvolution entry accepts a `ci` list; the new `intervals` table
    reports, per method and interval construction (Wald, ILR, profile,
    bootstrap, oracle, and the single-cell REML variants), the coverage
    with its Wilson Monte Carlo interval, the mean and median width and
    the share of intervals leaving \[0, 1\].
  - Figure scripts: `scripts/fig02_bivariate_toy.R`,
    `scripts/fig03_variance_covariance_driven.R` and the new
    `scripts/fig05_single_cell_hierarchical.R` (four cell types in two
    lineages, 32 + 1 scenarios) share one layout: run size from
    `DECOVART_N_REPLICATES`, `DECOVART_N_BOOT`, `DECOVART_N_REFERENCE`,
    `DECOVART_CORES` and `DECOVART_POSTPROCESS_ONLY`;
    `figXX_benchmark.rds` with every result and `figXX_plot_data.rds`
    keyed by figure file stem. Figures 3 and 5 add MEAD (R package),
    RNA-Sieve and ReDeconv (Python, through a `uv` project in
    `scripts/benchmark_metrics/python/`) with gene filters off; DECALS
    has no public code and is not run.
  - Single-cell intervals and bootstraps also cover the platform scale a
    = phi / c and the cell count n = c^2 / phi (rows `platform_scale_a`
    and `cell_count_n`), formed on the log scale so both stay positive;
    the benchmark scores them against the true a and n.
  - [`generate_random_network_skeleton()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_random_network_skeleton.md)
    accepts `block_sizes` for a stochastic block model with fixed,
    contiguous blocks, and returns the block of every gene as attribute
    `"block"`. Figures 3 and 5 use three blocks of equal size; in figure
    5 the two siblings of a lineage share one skeleton (block model for
    types 1-2, scale-free for types 3-4) and differ only by 5% of sign
    flips on the lineage’s edge weights.
  - Figures 2, 3 and 5 run NNLS next to `lsei`; figures 3 and 5 add the
    robust linear model
    ([`deconvolute_ratios_rlm()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)).
    RNA-Sieve is scored with its released interval and with the paper’s
    Fisher and Godambe intervals recomputed on all genes.
  - `scripts/auxiliary/` is split into `scripts/benchmark_metrics/` and
    `scripts/visualisation/`; the vignette figure scripts move to
    `vignettes/scripts/`; the bibliography tools in
    `scripts/curate_bibliography/` are a `uv` project.
  - Theory vignettes are reorganised into population and single-cell
    references, each split into generative model and optimisation, then
    inference: `theory-population-generative-optimisation`,
    `theory-population-inference`,
    `theory-single-cell-generative-optimisation` and the new
    `theory-single-cell-inference`. New vignette
    `fig05-single-cell-hierarchical`.

- **Fixes to the interval and test machinery, found by simulation.**
  `scripts/benchmark_metrics/check_confidence_interval_coverage.R`
  checks the information matrices against independent computations and
  the coverage of every interval method (300 replicates).

  - [`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)
    and
    [`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)
    used BFGS from the barycentre; with pooled columns the first step
    could jump to a vertex, where the likelihood is flat, and stop there
    (an estimate of (1, 0) for a true (0.7, 0.3)). Both now run the
    trust-region Newton driver on the face with the pooled analytic
    Hessian (`.fit_on_face()`). This affects
    [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
    [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
    [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)
    and the pooled
    [`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md)
    menu.
  - A profile interval whose upper end reaches the face returned the
    estimate instead of 1; a failed refit now gives `NA` with a warning.
  - A null `p_j = 1` fixes J - 1 coordinates on faces; the
    chi-bar-square tests now count them.
  - `vcov_single_cell(scales_known = TRUE)` gives the covariance of a
    fixed-scale or two-stage fit instead of the profile covariance.
  - The reference bootstraps regenerate the bulk from the fitted model
    and refit it with the resampled moments; drawing it from the
    resampled moments hid the reference error (coverage 0.76 instead of
    0.92).
  - Bootstraps report `n_failed`.
  - `confint_single_cell(variance = "reml")` and
    `confint_profile_single_cell(variance = "reml")` correct for the G -
    J residual degrees of freedom of the variance scale (covariance
    inflated by G / (G - J), Student-t quantiles, chi-square interval
    for phi, exact linear-model profile cut-off). With G = 10 the
    default `"ml"` intervals cover 0.84 to 0.93 at the 0.95 level.

- **Smaller, flatter internals.** `equicorrelation_gram()` is no longer
  exported: `generate_mean_signature_matrix(target_cosine = )` builds
  the same Gram matrix, and its generated means are unchanged bit for
  bit. The gene frames, the Cholesky test, the two precision-inversion
  loops, the one-line metrics and the convergence helpers of the
  benchmark are merged or inlined. Every solver returns one convergence
  record (`code`, `iterations`, `message`, `rdm`);
  [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  loses the unused `lower` and `upper` arguments. A seeded gene frame
  restores the global random state. Simulated covariances are now
  inverted through the Cholesky factor of the precision (exactly
  symmetric, three to five times faster than QR, no less accurate);
  simulated datasets change in their last bits, and through
  [`MASS::mvrnorm()`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) in
  their draws, so benchmarks must be rerun.

- **Confidence intervals for both mixture models.**
  [`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md)
  chooses an interval for the population convolution. The default maps a
  Wald interval for the ILR coordinate through the inverse chart
  ([`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md)).
  `wald` is the delta method, `firth` is a Jeffreys-penalised Wald
  interval, `profile` inverts the likelihood ratio, and `bootstrap` /
  `reference` are the existing resamples.
  [`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  is unchanged and remains the per-column delta method.
  [`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
  uses the same names for the profile likelihood. Intervals for the
  composition come from the Schur complement that removes the two
  scales; intervals for the scales come from the complementary Schur
  complement, on the log scale by default.
  [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)
  now returns `vcov`, `vcov_rho` and `vcov_nuisance`. Higher-order and
  sandwich corrections for a plug-in covariance are not included. With
  several donors in every cell type the reference bootstrap resamples
  donors and keeps the bulk fixed; otherwise it resamples columns within
  each type and regenerates the bulk. Cell-type labels are not shuffled.

- **Single-cell mixture, separate from the population convolution.**
  Cell-level moments and the bulk sum of cells are
  [`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md)
  and
  [`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md).
  [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)
  starts from a free-scale non-negative least-squares fit, then profiles
  the mean and variance scales out in closed form and runs
  Newton-Raphson on the profile log-likelihood of the ILR coordinates
  ([`profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md),
  exact score, Schur-complement Hessian). One loop reaches the joint
  maximum; the former scale/composition alternation needed up to 95
  outer updates. `profile = FALSE` keeps the two-stage plug-in fit for
  comparison. The fit requires `G > J`. Default reference: raw cell
  covariance and cell fractions. The population log-likelihood is
  unchanged (`p^2` covariance weights).

- **Shared ILR chart and Newton driver.** The simplex chart
  ([`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md),
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md),
  [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md),
  …) is in `R/02_00_ilr_chart.R`; the cached covariance factorisation
  and the `nlminb` trust-region driver are in
  `R/02_00_likelihood_core.R`. Files specific to the `p^2` convolution
  carry a `population_` tag (`R/01_01_population_toy_simulation.R`,
  `R/02_02_population_decovart_estimate_ratios.R`,
  `R/02_03_population_decovart_fit.R`). The Newton driver now sets
  `abs.tol = 0`: PORT’s absolute-function test assumes a nonnegative
  objective, and a negative log-likelihood can cross zero.

- **Benchmark code split into run, metrics and visualisation.**
  `R/04_01_run_benchmark.R` fits every algorithm to every sample;
  `R/04_02_run_simulation_benchmark.R` runs the scenario grid;
  `R/04_03_performance_metrics_global.R` scores the whole composition
  vector; `R/04_04_performance_metrics_celltype.R` scores each cellular
  ratio across replicates (Pearson, presence F1, bias, coverage). Plots
  move to `R/05_01_monte_carlo_plots.R` and
  `R/05_02_algorithm_similarity.R`, with a third call graph. Entropy,
  MixSim overlap, Jeffreys divergence and the () parser move to the
  scenario files (`01_03`, `01_04`);
  [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  to the estimation file.

- **Removed from the package.** `plot_correlation_Heatmap()` (its input
  columns no longer exist in benchmark output), `check_true_theta()` (a
  wrapper of the internal parser), `spd_affine_invariant_distance()`
  (folded into
  [`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)),
  `compute_glmnet_gene_scores()` (moved, with its cross-validated
  variant, to `scripts/experimental/compute_glmnet_gene_scores.R`), and
  the unexported file writers. Scores that do not apply in a benchmark
  table are now `NA` instead of the string `"irrelevant"`. The package
  no longer writes files (srr G4.0 is now not applicable).

- **Script helpers out of `R/`.** `write_simulation_artefacts()`,
  `read_simulation_artefacts()`, the fig02 / fig03 factor orders and
  display labels, and the expected-Fisher Wald refresh move to
  `scripts/auxiliary/benchmark_artefacts.R`. The figure-book scripts
  become `scripts/visualisation/05_03` to `05_06`, loaded together by
  `scripts/visualisation/load_visualisation.R`. The fig02, fig03 and
  regen scripts now source that loader: before, they called the book
  functions without sourcing the files that define them.

- **Bug fix: Hessian for four or more cell types.**
  [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md)
  filled its lower triangle with `H[lower.tri(H)] <- H[upper.tri(H)]`,
  which pairs the wrong entries as soon as (J): the returned matrix was
  neither symmetric nor the Hessian. Newton–Raphson (the
  [`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  default), Marquardt–Levenberg,
  [`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md)
  and the eigenvalues of
  [`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md)
  were affected for (J); the score and the Fisher (Wald) covariance were
  not. The score, Hessian and expected Fisher information are now
  vectorised around one shared information block (about six times
  faster), and a (J=4) test checks the Hessian against `numDeriv`.

- **Dense covariance only; structured backends move to
  `scripts/experimental/`.** `new_decovart_covariance()`,
  `sigma_logdet()`, `sigma_solve()`, `sigma_quadform()`,
  `sigma_trace_precision_times()` and
  `covariance_structure_from_graph_model()` leave the package, and
  [`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
  loses its `backend` argument. The appendix S2 script sources
  `scripts/experimental/covariance_backends.R`. `Matrix` is no longer
  imported.

- **Fixed-covariance GLS competitor removed.** `fixed_gls_covariance()`
  and `deconvolute_ratios_gls()` leave the package; the S4 script and
  the S1 vignette compute the mean-only GLS fit locally.

- **Inference split by geometry.** `R/03_01_interior_intervals.R` (Wald
  and profile-likelihood intervals), `R/03_02_boundary_tests.R`
  (restricted MLE on faces, chi-bar-square, likelihood-ratio tests,
  boundary diagnostics) and `R/03_03_resampling_inference.R` (parametric
  and reference bootstraps, equivariance check).
  [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  moves next to
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md).

- **Simpler inference interface.**
  [`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)
  is now the single fitting engine: `y` may be a vector or a (GN) matrix
  of replicates, and `fixed = NULL` gives the unrestricted MLE on the
  closed simplex. The `bulk_expression` argument of
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
  [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
  [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
  [`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
  and
  [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  is removed: pass replicates as the columns of `y`. Constrained ratios
  must be named by cell types (integer indices are no longer accepted).
  [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  takes a named list of purified matrices only (`cell_type_labels`
  removed). `profile_loglik_decovart()` and the internal
  `.inner_product()` are removed.

- **File layout and dropped ALR helpers.** Source files now follow the
  article: `01_` simulations, `02_` optimisation and estimation, `03_`
  intervals, `04_` benchmarks. `additive_logistic()`,
  `additive_log_ratio()` and their Jacobian and Hessian leave the
  package; solvers use the isometric log-ratio chart. Figure scripts
  that still draw an ALR surface keep those two maps under
  `scripts/auxiliary/`.

- **Confidence intervals: regularity, geometry and the three tests.**
  The MLE vignette and the article appendix now state the regularity
  conditions of the large-sample approximation (interior truth, fixed
  dimension, replication over biological samples, not genes), add the
  small-dispersion limit in the number of cells for one bulk sample,
  compare Wald, likelihood-ratio and score tests and their inversion
  into intervals, treat nuisance scales through the envelope theorem and
  parameter orthogonality, and describe higher-order corrections
  (modified likelihood root, Bartlett). The logit-scale Wald route and
  the variational sandwich are removed. The article now splits each
  model into optimisation and inference, moves the chart choice, the
  single-cell fit and the robust variants to the appendix, derives
  Fisher scoring as iteratively reweighted least squares, and groups the
  perspectives into modelling, optimisation and estimating-equation
  inference (GEE, GMM, HC corrections).

- **[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  loses `method = "dirichlet"`.** The Dirichlet composition sweep
  measured performance across random mixtures, not uncertainty for one
  observed bulk, so it is removed together with its `dirichlet_alpha`
  argument and `p_simulated` output. The function now always needs `y`
  or `bulk_expression` and resamples donors (default) or cells within
  each cell type, with replacement. Random Dirichlet starting points of
  the solvers are unchanged.

- **Fig02 performance books.** Twelve meta-scenario pages (small then
  large CLD; homoscedastic then heteroscedastic; balanced / moderately /
  highly unbalanced), each a 2-by-2 of correlation corners ((0,0)),
  ((-0.8,-0.8)), ((0.8,0.8)), ((-0.8,0.8)). Similarity heatmaps cluster
  paired (p) within each corner and draw a dendrogram on the **right**,
  flush with the tiles. Wald forests omit NNLS, LSEI, and SA. Whiskers
  use the expected-Fisher SE at the true composition
  (`confint.decovart_fit` / `vcov_ilr_delta`; constant per scenario)
  together with the empirical Monte Carlo SD; bold two-line `geom_label`
  boxes (one per solver, vertically centred between the two cell types)
  report RMSE and coverage, which coincide for (p_1) and (p_2) on the
  unit simplex. Rainclouds plot (p) with vertical lines at the true
  proportions on a 30-inch canvas and translucent slabs. A 12-page
  solver-dot book uses the same 9-by-9 correlation factorial as the RMSE
  heatmaps (enrichplot-style: colour RMSE, size Aitchison; horizontal
  y-axis labels). ggplot `data` for each book (not rgl) is written to
  `output/fig02/ggplot_rds/`. Scenario factors are relevelled with
  `forcats` (Suggests) when artefacts are read.

- **Fig02 density books.** Kernel-density rasters fill the gene-axis
  window (no inner white frame). Bulk pages drop the 95% Gaussian
  ellipses with a (^2_2) footnote. ILR profiles use () so unbalanced
  compositions are not cropped, and plot the relative likelihood L/\max
  L on a log10 y-axis with `annotation_logticks()`. Page titles follow
  the same CLD / variance / composition order as the performance books.
  Log-likelihood outputs split into `loglik_surface_p.pdf` (ratio
  lattice), `loglik_ilr_profile.pdf` (ILR () profile with MLE), and
  `loglik_rgl.pdf` plus interactive `loglik_rgl.html`.

- **Fig03 covariance-driven grid.** G=20, J=3. Graph generators are
  scale-free and cluster (stochastic block model) only: a 2^{2}
  assignment on the mean-collinear types, type 3 held at scale-free,
  times three *target MixSim average overlaps* (global covariance scale,
  precision zeros kept; 4\times 3\times 3=36 scenarios).
  `hybrid_config.rds` stores `graph_generators` (model + `graph_params`
  per cell type) alongside `true_theta`, plus realised `BarOmega`,
  f\_{\mathrm{cov}}, and pairwise AIRM distance. Convolution-likelihood
  solvers attach interior simplex standard errors from
  [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
  inside
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  (the same object as
  [`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  /
  [`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)).
  Mean-only solvers (NNLS, LSEI) leave coverage, mean model SE, and
  `se_sd_ratio` missing. The ALR helper `vcov_alr_delta()` has been
  removed.

- **Fig03 latent-space and resource books.**
  `save_hybrid_latent_projection_book()` writes a 16-page 2-D book
  (extreme MixSim overlap () extreme Shannon composition () four CT1/CT2
  topologies). The top row is independent Thomson
  [`factanal()`](https://rdrr.io/r/stats/factanal.html) scores per type;
  the bottom row is a shared MCFA plane on convolution draws, MCFA on an
  unsupervised mixture with weights (), and supervised `MclustDR`.
  Runtime and memory rainclouds (`save_hybrid_runtime_book()`,
  `save_hybrid_memory_book()`) use one page per composition, facet
  overlap, and log10 y with a rug of Monte Carlo draws. Fig02 has the
  matching pair `save_bivariate_runtime_book()` /
  `save_bivariate_memory_book()` (correlation corners on the x-axis).

- **Overlap in moderate dimension.**
  [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
  uses MixSim Davies quadrature for G\<4 and stratified Sobol Monte
  Carlo
  ([`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md),
  default 10^4 draws per component) when G\ge 4.
  [`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md)
  calibrates a global covariance scale to a target average overlap.
  Affine-invariant Riemannian distance of the \Sigma_j is
  `riemannian_sigma`. Theory:
  `vignettes/theory-distance-covariance.qmd`.

- **Faceted plots.**
  [`theme_decovart_facets()`](https://bastienchassagnol.github.io/DeCovarT/reference/theme_decovart_facets.md)
  (black strips, panel border, after the
  atlas-feature-selection-benchmark `theme_features`) is the default for
  forest, raincloud, metric tiles, and similarity heatmaps. Fig02 PDFs
  are written under `output/fig02/density_visualisations/` and
  `output/fig02/performance_visualisations/` at a wide canvas and 320
  dpi.

- **Simulation artefacts.** `write_simulation_artefacts()` /
  `read_simulation_artefacts()` persist the design grid, scenario
  descriptors (including Jeffreys), `\theta`, and ADEMP metrics as four
  RDS files keyed by `ID`. The written config uses the canonical design
  columns (`proportions`, `centroids`) and drops `true_theta` (stored in
  the theta artefact) and `scenario_idx` (encoded in `ID`). Fig02 IDs
  are `B{idx}_{Ho|He}_{Ba|Mo|Hi}_{Sm|Lg}`.

- **Scenario descriptors.**
  [`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md)
  now reports `max_cosine`, `mean_euclidean` (centroid / CLD gap),
  unweighted pairwise Hellinger (`hellinger`; the distance is symmetric)
  plus optional `hellinger_weighted`, and keeps Jeffreys as a
  supplementary column on the descriptors artefact.

- **Bivariate figures.** Purified and bulk 2-D densities fill the
  gene-axis window with a
  [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html) raster
  (viridis). Purified pages keep red-circle / green-triangle centroids
  and exact 95% Gaussian ellipses (`gaussian_confidence_ellipse()`);
  bulk pages keep centroids only. Bulk pages collect a shared legend
  with
  [`cowplot::get_legend()`](https://wilkelab.org/cowplot/reference/get_legend.html).
  Log-likelihood books evaluate ((p_1,p_2)) on a ratio lattice (not gene
  space): ggplot rasters versus a true
  [`rgl::persp3d()`](https://dmurdoch.github.io/rgl/dev/reference/persp3d.html)
  surface that marks the simulation proportions.
  `plot_mc_raincloud(include_dots = FALSE, max_rows = …)` avoids the
  `quadprog` long-vector crash on the 972-scenario grid.

- **Negative bulk.** `.prepare_deconvolution_inputs()` warns on negative
  expression instead of aborting.
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
  no longer floors Gaussian draws at 0.

- **Scenario 3 (variance-driven).** Fixed Gram \cos(\mu_1,\mu_2)=0.9,
  \cos(\mu_j,\mu_3)=0.1; 2^3 graph assignments (Erdős–Rényi /
  scale-free) \times three precision cushions; compositions from
  [`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md)
  at H^{\star}\in\\1,0.5,0.1\\; solvers LSEI, CIBERSORT, L-BFGS-B,
  Newton–Raphson, Marquardt–Levenberg (n=50, G=10). Descriptor docs use
  one callout per family in the synthetic-scenarios vignette.
  `kappa_sigma_p` is now returned beside the reciprocal ratio.
  `generate_mean_signature_matrix(nonnegative = TRUE)` uses a
  disjoint-support frame so signatures stay nonnegative. Graph
  precisions are completed by a uniform spectral shift and extra-loaded
  until both \Omega and \Sigma=\Omega^{-1} admit a Cholesky factor
  (support and signs of W unchanged). Failed \hat p entries are treated
  as absent when scoring presence F1.
  `simulate_bulk_mixture(truncate_negative = TRUE)` floors bulk
  coordinates at 0 (used by
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)).

- **Terminal UI.** Optional Suggests `cli` formats messages, warnings
  and errors (`cli_alert_*` / `cli_abort`, with base
  [`message()`](https://rdrr.io/r/base/message.html) /
  [`stop()`](https://rdrr.io/r/base/stop.html) fallback).
  `run_simulation_benchmark(verbose = TRUE)` logs each scenario row;
  grids of 10 or fewer scenarios also tick every 10 inferred samples.

- **Documentation layout.** Paper scenarios live in the
  `fig02-bivariate-toy` and `fig03-variance-driven` articles (one script
  each: `scripts/fig02_bivariate_toy.R`,
  `scripts/fig03_variance_driven.R`). Regular-case MLE checks and
  identifiability sit in Appendix S1. ADEMP / Nature Methods reporting
  is at the end of
  `vignettes/theory-synthetic-scenarios-mean-covariance.qmd`.

- **Monte Carlo figures.**
  [`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md),
  [`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md)
  (`ggplot2` + `ggdist` horizontal rainclouds of \hat p_j-p_j^{\star}),
  [`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md)
  (ADEMP summaries with Wilson whiskers on the coverage rate),
  [`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md)
  (`geom_tile` + `hclust(1-r)`), and
  [`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md)
  (faceted min-max scores; no default composite). Optional Suggests:
  `ggdist`, `ggdendro`. `ComplexHeatmap` remains for linked multi-omics
  grids only.

- **Solver starts.**
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md)
  (used by all ILR maps) accepts the barycentre (default), a Dirichlet
  draw (`initial_p = "dirichlet"`, `dirichlet_alpha = 1` uniform on the
  simplex; \alpha\>1 centre-biased; \alpha\<1 face-biased; several
  independent draws for multistarts), or a mean-only simplex QP
  (`initial_p = "qp"` via
  [`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)).
  [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  uses the Dirichlet path.

- **Fixed-covariance GLS competitor.** `deconvolute_ratios_gls()` wraps
  [`MASS::lm.gls()`](https://rdrr.io/pkg/MASS/man/lm.gls.html) with a
  known G\times G residual covariance; `fixed_gls_covariance()` builds
  the global-diagonal W at a declared p (default 1/J). Unconstrained OLS
  (`deconvolute_ratios_lsfit`) has been removed.

- **ILR / Helmert simplex chart.** Frequentist solvers and
  [`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  now use isometric log-ratio coordinates on a Helmert basis
  ([`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md),
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md),
  [`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md),
  [`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md),
  [`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md),
  [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)).
  No cell type is pinned as a reference; changing the ILR basis by an
  orthogonal rotation leaves (), (()) and (()) unchanged. Additive
  log-ratio maps remain exported for the derivatives vignette appendix.
  Restricted MLE assembly treats a fully constrained face without
  calling ILR on an empty coordinate, and an unrestricted fit that sits
  near a face is compared with the exact-face restricted MLE so the
  reported point can be the closed-simplex supremum.

- **Gram-matrix mean signatures.**
  [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
  realises a target cosine Gram via the symmetric square root
  (=s,QR^{1/2}) (`equicorrelation_gram()` or `target_gram`). Pairwise
  cosines match (R) exactly; the previous shared-plus-private blend is
  documented as an appendix construction that does not hit () at finite
  (J).

- **Identifiability.** Full column rank of () is sufficient but not
  necessary: identical means remain identifiable when the cell-type
  covariances differ. The article and MLE vignette now state injectivity
  of (p(p,(p))).

- **Benchmark metrics and sample-level parallelism.**
  [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  now returns a three-block list: composition / regression scores
  (global TV, RMSE, angular distance, SDID, MaxAE; cell-type Pearson,
  presence F1, and false-positive mass), ADEMP Monte Carlo summaries,
  and optimisation / runtime diagnostics (simplex KKT residual,
  numerical versus theoretical convergence, per-sample elapsed time, and
  `ps` PSS memory). Coverage of the Monte Carlo coverage *rate* uses a
  Wilson interval by default
  ([`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md);
  Wald and Agresti–Coull optional).
  [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
  also returns `theta_true` (the convolution parameters), `descriptors`
  / `supplementary` from
  [`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md),
  and `call` ([`match.call()`](https://rdrr.io/r/base/match.call.html)).
  MixSim `BarOmega` and averaged pairwise Hellinger are kept in
  `descriptors`; Jeffreys overlap is supplementary.
  [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  iterates bulk columns with `furrr` when `cores > 1`, using
  L’Ecuyer-CMRG streams (`furrr_options(seed = TRUE)`). Scenario rows
  are sequential, so workers are never nested. There is no composite
  global score.

- **Structure-aware covariance backends.** New
  `new_decovart_covariance()` wraps with a declared `structure`
  (`"dense"`, `"block"`, `"band"`, `"sparse"`, `"diag_lowrank"`) and
  exposes operators `sigma_logdet()`, `sigma_solve()`,
  `sigma_quadform()` and `sigma_trace_precision_times()` that return the
  log-determinant and solves without materialising a dense precision.
  Block Cholesky (disconnected modules / stochastic block model), banded
  Cholesky (band / AR), sparse Cholesky with a cached symbolic ordering
  (Erdős–Rényi, scale-free, small-world), and the Woodbury /
  matrix-determinant-lemma path for diagonal-plus-low-rank covariances
  (shared regulatory programs; `Sigma = W W^T + Psi`) each match the
  dense default to machine precision.
  `covariance_structure_from_graph_model()` maps a network topology to
  the recommended backend, and
  [`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
  gains an optional `backend` argument. Band and sparse backends use the
  imported `Matrix` package. See the new “Structure-aware covariance
  backends” section of `vignette("theory-decovart-generative-model")`
  and `scripts/supp_S2_covariance_inversion.R`.

## DeCovarT 2.3.1

- **CRAN resubmission.** Tests no longer write under the package tree.
  The bivariate toy scenario is built in memory by
  `new_bivariate_toy_scenario()` in `tests/testthat/helper.R`; the
  former `fixtures/make-useful-things.R` writer and its golden RDS files
  are gone. CRAN Quarto vignettes do not execute
  [`library(DeCovarT)`](https://github.com/bastienchassagnol/DeCovarT)
  during rebuild (Windows CLI subprocess cannot see `R CMD build`’s
  temporary library).

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
  `profile_loglik_decovart()`,
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

- **Numerical stability of `additive_logistic()`.** Now evaluated
  through the log-sum-exp shift instead of `exp(rho) / sum(exp(rho))`,
  which overflowed to `NaN` for and could make non-factorisable when the
  MLE approached a simplex face.

- **New vignette** `vignette("DeCovarT-MLE-properties")`:
  identifiability, existence, failure of uniqueness for one observed
  sample, population consistency through the Kullback–Leibler divergence
  and the Fisher information metric, affine (but not logarithmic)
  equivariance, the four strata of the simplex and their tangent cones,
  and why replication is what makes the asymptotics meaningful.

## DeCovarT 2.2.3

- **Optional heatmap stack.** Moved `ComplexHeatmap`, `circlize`, and
  `viridis` from Imports to Suggests. `plot_correlation_Heatmap()`
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
  `vignette("generative-model-derivatives")`, section *Numerical
  speed-ups and solver safeguards*). The CRAN tarball ships the
  use-cases and softmax/ALR vignettes; remaining articles stay on the
  pkgdown site.

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
  `vignette("generative-model-derivatives")` and
  `tests/testthat/test-03_03_DeCovarT.R`).

- **Simulation framework beyond the bivariate toy.** Synthetic first-
  and second-order moments now support arbitrary gene panels (G\gg 2)
  and multi-type mixtures: AutoGeneS-style mean signatures with a cosine
  dial, graph-constrained precisions (ER, hub/star, scale-free, SBM,
  small-world), and bulk convolution via
  [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md).
  A hybrid multi-topology reference scenario (G=50, J=3) in
  `scripts/generate_random_markov_network.R` and
  `vignette("synthetic-scenarios")` stresses feature selection under EE
  / DE / differential-modality-like blocks (scDD / `muscat`-inspired
  taxonomy).

- Feature-selection metrics
  ([`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
  [`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md),
  `compute_glmnet_gene_scores()`) and a shared `check_true_theta()`
  validator complete the end-to-end simulation → pre-screen → NSGA-II
  refinement loop documented in `vignette("feature-selection")`.

## DeCovarT 0.1.0

- First official release of DeCovarT: bulk transcriptomic deconvolution
  that accounts for gene-gene covariance in purified reference
  populations.
