#' @srrstatsVerbose TRUE
#' @noRd
NULL

#' srr_stats
#'
#' rOpenSci Statistical Software Standards - General (G) and Regression (RE)
#' categories.  Standards already addressed are tagged `@srrstats`.
#' Standards not applicable are tagged `@srrstatsNA` in the block below.
#' Remaining standards keep `@srrstatsTODO` for full review.
#'
#' G1 - Documentation and references
#'
#' G1.0, G1.1: see [deconvolute_ratios()] (DSection as closest method;
#' first multivariate Gaussian-convolution MLE).
#' G1.3, G1.4, G1.4a: see `DeCovarT-package`.
#' G2.3 / G2.3a / G2.3b: see `.match_arg_ci()`, [check_true_theta()],
#' [generate_random_network_skeleton()], [plot_correlation_Heatmap()].
#' G2.4-G2.16, G4.0, G5.8d: see `R/utils.R` and [deconvolute_ratios()].
#' G3.0: see [repair_simplex()] and tests (`100 * .Machine$double.eps`).
#'
#' @srrstats {G1.2} A lifecycle badge (Active/stable,
#'   <https://www.repostatus.org/#active>) is shown in the README; future
#'   development is described there and in vignettes.
#'
#' @srrstatsTODO {G1.5} Reproduction code for article figures will live in a
#'   companion GitHub repository (CellRank 2 pattern: software vs
#'   `cellrank2_reproducibility`), not in this package.
#'
#' @srrstatsTODO {G1.6} Comparative benchmarks against other R deconvolution
#'   tools will also live in that companion repository / vignette, once the
#'   simulation grid is extended.
#'
#' G2 - Input validation
#'
#' @srrstats {G2.0} Input length assertions implemented via `checkmate`
#'   (e.g. `assert_matrix`, `assert_numeric`) throughout estimation functions.
#'
#' @srrstats {G2.0a} Expected dimensions and lengths documented in `@param`
#'   of all exported functions.
#'
#' @srrstats {G2.1} Type assertions performed via `checkmate` on all inputs
#'   (numeric matrices, positive-definite covariance matrices, etc.).
#'
#' @srrstats {G2.1a} Expected types of all inputs documented in `@param`.
#'
#' @srrstats {G2.2} Scalar parameters are checked with `assert_scalar` or
#'   equivalent.
#'
#' G2.3, G2.3a, G2.3b: see `.match_arg_ci()` in `R/utils.R`.
#' G2.13, G2.15, G2.16: see `.assert_no_missing()` and
#' [deconvolute_ratios()] (error on missing / non-finite values; no
#' imputation).
#'
#' @srrstats {G5.0} Tests use internally generated data from
#'   `simulate_bulk_mixture()` / `simulate_hierarchical_grn_moments()`
#'   with known ground-truth proportions, plus
#'   `inst/extdata/toy_deconvolution.rds`.
#'
#' @srrstats {G5.1} The toy two-gene / two-cell-type list is shipped in
#'   `inst/extdata/toy_deconvolution.rds`; bivariate solver fixtures live
#'   under `tests/testthat/fixtures/` (see `tests/README.md`).
#'
#' @srrstats {G5.8d} \eqn{J > G} is rejected in
#'   `.prepare_deconvolution_inputs()`.
#'
#' @srrstats {G5.10} Extended tests are gated by
#'   `DECOVART_EXTENDED_TESTS=true` (`skip_if_not_extended()`).
#'
#' @srrstatsTODO {G5.6b} Multi-seed long-running tests to be added to an
#'   extended test suite.
#'
#' @srrstatsTODO {G5.7} Runtime scaling benchmarks to be added.
#'
#' @srrstatsTODO {G5.9b} Multiple random seeds to be added (see G5.6b).
#'
#' @srrstats {G3.1} User-provided covariance matrices can originate from any
#'   estimation method; the diagonal (independence) approximation is one valid
#'   choice documented in the vignettes.
#'
#' @srrstats {G3.1a} Covariance estimation strategies documented with worked
#'   examples in the `synthetic-scenarios` vignette.
#'
#' G5 - Testing
#'
#' G5.0, G5.1, G5.8d, G5.10: see the G2 block above and `tests/README.md`.
#'
#' @srrstats {G5.2} Error and warning conditions tested via `expect_error()` /
#'   `expect_warning()` covering malformed inputs and constraint violations.
#'
#' @srrstats {G5.2a} Each `stop()` / `warning()` / `cli::cli_abort()` call
#'   uses a unique message string.
#'
#' @srrstats {G5.2b} Tests explicitly trigger every error/warning condition and
#'   compare against expected strings.
#'
#' @srrstats {G5.3} Return objects tested for absence of `NA` / `NaN` on
#'   valid inputs.
#'
#' @srrstats {G5.4} Correctness tests compare estimated proportions to
#'   ground truth within defined tolerance.
#'
#' @srrstats {G5.4a} Two-cell / two-gene cases used where MLE is analytically
#'   derivable and verifiable by hand.
#'
#' @srrstats {G5.5} All stochastic tests use `set.seed()`.
#'
#' @srrstats {G5.6} Parameter recovery tests confirm convergence to ground
#'   truth as sample size increases.
#'
#' @srrstats {G5.6a} Recovery tests use `expect_equal(..., tolerance = 1e-3)`,
#'   not exact equality.
#'
#' @srrstats {G5.8} Edge tests cover zero-proportion cell types, near-singular
#'   covariance matrices, and single-cell-type degenerate cases.
#'
#' @srrstats {G5.8a} Zero-length input matrices trigger informative errors.
#'
#' @srrstats {G5.8b} Non-numeric inputs raise errors via `checkmate`.
#'
#' @srrstats {G5.8c} All-`NA` columns trigger an early-exit error.
#'
#' @srrstats {G5.9} Noise susceptibility confirmed: perturbing bulk mixtures
#'   with small Gaussian noise changes estimates by less than a defined
#'   threshold.
#'
#' @srrstats {G5.9a} Trivial noise at `.Machine$double.eps` scale does not
#'   materially alter proportion estimates (tested).
#'
#' @srrstats {G5.12} Extended-test flags and timing notes are in
#'   `tests/README.md`.
#'
#' RE1 - Regression: input specification
#'
#' @srrstats {RE1.1} DeCovarT does not use a formula interface: inputs are
#'   numeric matrices (`y`, `mean_signature_matrix`, `Sigma_list`). This is
#'   explicitly documented in `?deconvolute_ratios` and justified by the
#'   multivariate nature of the convolution model.
#'
#' @srrstatsTODO {RE1.0} Formula interface or explicit justification - the
#'   RE1.1 entry addresses this; will be moved to `@srrstatsNA` at full
#'   submission once editors confirm acceptability.
#'
#' @srrstats {RE1.2} Expected input classes documented in every `@param`
#'   block; unsupported types raise `checkmate` errors.
#'
#' @srrstats {RE1.3} Output tibble retains column names from signature matrix
#'   and row names (sample IDs) from bulk matrix `y`.
#'
#' @srrstatsTODO {RE1.3a} Cases where input metadata are not transferred to
#'   be audited and documented.
#'
#' @srrstats {RE1.4} Distributional assumptions (multivariate Gaussian,
#'   positive-definite Sigma, simplex-constrained proportions) documented in
#'   `?deconvolute_ratios` and vignettes, with notes on violation consequences.
#'
#' RE2 - Regression: pre-processing
#'
#' @srrstats {RE2.0} ALR reparametrisation documented in `?deconvolute_ratios`,
#'   vignettes, and preprint; exposed via `alr_transform` / `alr_inverse`.
#'
#' @srrstats {RE2.1} Missing values in `y` or signature matrix raise an error
#'   at the start of `deconvolute_ratios()`; no silent imputation is performed.
#'
#' @srrstatsTODO {RE2.2} Separate NA handling for predictor vs response data
#'   to be added.
#'
#' @srrstatsTODO {RE2.3} Centring / scaling options not currently exposed; to
#'   be reviewed.
#'
#' @srrstats {RE2.4} Condition numbers of covariance matrices checked before
#'   inversion; warning issued when threshold exceeded.
#'
#' @srrstats {RE2.4a} Near-perfect collinearity among signature columns
#'   (identical cell profiles) triggers an informative warning.
#'
#' @srrstatsTODO {RE2.4b} Perfect collinearity between predictors and response
#'   to be explicitly tested.
#'
#' RE3 - Regression: convergence
#'
#' @srrstats {RE3.0} Non-convergence of the Marquardt-Levenberg optimiser
#'   raises a `warning()` with the convergence code.
#'
#' @srrstats {RE3.1} Convergence warnings suppressible via `verbose = FALSE`;
#'   the result object retains the convergence field.
#'
#' @srrstats {RE3.2} Default convergence thresholds (`epsa`, `epsb`, `epsd`)
#'   documented in `?deconvolute_ratios` with values and rationale.
#'
#' @srrstats {RE3.3} Convergence thresholds settable via `...` pass-through
#'   to `marqLevAlg::marqLevAlg()`, documented in `@param ...`.
#'
#' RE4 - Regression: return objects
#'
#' @srrstats {RE4.0} `deconvolute_ratios()` returns a named list
#'   (`DeCovarT_result`) with proportions, log-likelihood, Fisher information,
#'   and convergence diagnostics.
#'
#' @srrstatsTODO {RE4.1} Model object without fitting not currently supported.
#'
#' @srrstats {RE4.2} Estimated proportions (model coefficients) returned in a
#'   tidy tibble.
#'
#' @srrstats {RE4.3} Asymptotic confidence intervals from the observed Fisher
#'   information matrix returned alongside point estimates.
#'
#' @srrstatsTODO {RE4.4} `formula()` accessor - not applicable; see RE1.1.
#'   Will be moved to `@srrstatsNA` at full submission.
#'
#' @srrstats {RE4.5} Number of observations (`N` bulk samples) returned in
#'   result object.
#'
#' @srrstats {RE4.6} Variance-covariance matrix of estimated proportion
#'   parameters (inverse expected Fisher information) returned.
#'
#' @srrstats {RE4.7} Convergence statistics (code, iterations, gradient norm)
#'   from `marqLevAlg` stored in result object.
#'
#' @srrstats {RE4.8} Bulk mixture matrix `y` stored as attribute of result
#'   object.
#'
#' @srrstats {RE4.9} Fitted values (modelled bulk expression from estimated
#'   proportions and signature) returned.
#'
#' @srrstatsTODO {RE4.10} Residuals not yet formally returned; to be added.
#'
#' @srrstatsTODO {RE4.11} Goodness-of-fit metrics (RMSE, Pearson correlation)
#'   to be added; log-likelihood and AIC already returned.
#'
#' @srrstats {RE4.12} ALR forward and inverse transforms exported as
#'   `alr_transform` / `alr_inverse`.
#'
#' @srrstats {RE4.13} Signature matrix (predictor variables) and covariance
#'   list stored in result object.
#'
#' @srrstatsTODO {RE4.14} Extrapolation errors - not applicable to
#'   compositional estimation; will be `@srrstatsNA` at full submission.
#'
#' @srrstatsTODO {RE4.15} Forecast-horizon uncertainty - not applicable; will
#'   be `@srrstatsNA` at full submission.
#'
#' @srrstatsTODO {RE4.16} `predict()` S3 method dispatching on
#'   `DeCovarT_result` is planned.
#'
#' @srrstats {RE4.17} `print.DeCovarT_result()` displays estimated proportions,
#'   log-likelihood, and convergence status.
#'
#' @srrstatsTODO {RE4.18} `summary.DeCovarT_result()` method is planned.
#'
#' RE5 - Regression: scaling
#'
#' @srrstatsTODO {RE5.0} Runtime vs number of samples and genes to be
#'   benchmarked and documented in a vignette.
#'
#' RE6 - Regression: visualisation
#'
#' @srrstats {RE6.0} `plot.DeCovarT_result()` implemented in
#'   `R/04_02_benchmark_visualisation.R`; produces proportion bar plots and
#'   scatter plots of estimated vs true proportions where ground truth exists.
#'
#' @srrstats {RE6.1} `plot()` generic dispatches on `DeCovarT_result` via
#'   `plot.DeCovarT_result()`.
#'
#' @srrstats {RE6.2} Default plot shows fitted proportions per sample with
#'   optional confidence-interval whiskers.
#'
#' @srrstatsTODO {RE6.3} Interpolated vs extrapolated distinction - not
#'   applicable to compositional estimation; will be `@srrstatsNA` at full
#'   submission.
#'
#' RE7 - Regression: tests
#'
#' @srrstats {RE7.0} Tests with noiseless exact proportions included in
#'   `tests/testthat/test-03_03_*.R`.
#'
#' @srrstats {RE7.0a} Tests confirm rank-deficient signature inputs rejected
#'   with informative error.
#'
#' @srrstats {RE7.1} Noiseless exact predictor-response tests confirm estimated
#'   proportions equal ground truth within tolerance.
#'
#' @srrstatsTODO {RE7.1a} Speed comparison noiseless vs noisy fitting to be
#'   benchmarked.
#'
#' @srrstats {RE7.2} Tests confirm output tibbles retain row and column names
#'   from input matrices.
#'
#' @srrstats {RE7.3} Tests exercise accessor fields of `DeCovarT_result`
#'   (proportions, CIs, Fisher information, convergence diagnostics).
#'
#' @srrstatsTODO {RE7.4} Forecast-horizon tests - not applicable; will be
#'   `@srrstatsNA` at full submission.
#'
#' @noRd
NULL

#' NA_standards
#'
#' Standards not applicable to DeCovarT, with justifications.
#'
#' @srrstatsNA {G2.4} No type coercion is performed: callers must supply
#'   numeric matrices / arrays. `as.integer` / `as.numeric` / `as.character`
#'   / `as.factor` helpers are therefore not part of the API.
#' @srrstatsNA {G2.4a} See G2.4.
#' @srrstatsNA {G2.4b} See G2.4.
#' @srrstatsNA {G2.4c} See G2.4.
#' @srrstatsNA {G2.4d} See G2.4.
#' @srrstatsNA {G2.4e} See G2.4.
#' @srrstatsNA {G2.5} No package argument expects a factor. The only
#'   ordered factor is in `scripts/bivariate_toy_model.R` (outside `R/`).
#' @srrstatsNA {G2.7} The public API accepts numeric matrices and arrays
#'   only, not `data.frame` or other tabular forms.
#' @srrstatsNA {G2.11} Units / non-vector columns cannot arise: inputs are
#'   not `data.frame`s.
#' @srrstatsNA {G2.12} List columns cannot arise for the same reason.
#' @srrstatsNA {G2.14b} Ignoring missing values is not offered.
#' @srrstatsNA {G2.14c} Imputation is not offered (RNA-seq matrices are
#'   complete; proteomic missingness is a future, separate API).
#' @srrstatsNA {G5.4b} DeCovarT is the first implementation of this
#'   multivariate Gaussian-convolution MLE (G1.1); there is no prior R
#'   implementation to test against.
#' @srrstatsNA {G5.4c} There is likewise no published numerical oracle for
#'   this estimator (same G1.1 justification).
#' @srrstatsNA {G5.11} Unit tests do not download assets. Paper-scale data
#'   will live in the companion reproducibility repository.
#' @srrstatsNA {G5.11a} See G5.11.
#'
#' @noRd
NULL
