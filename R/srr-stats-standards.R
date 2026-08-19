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
#' @srrstats {RE1.2} Expected input classes documented in every `@param`
#'   block; unsupported types raise `checkmate` / `.assert_numeric_array`
#'   errors.
#'
#' @srrstats {RE1.3} Gene rownames of \(Y\) and \(\mu\) identify genes;
#'   colnames of \(Y\) identify samples; colnames of \(\mu\) identify
#'   cell types. `fit_decovart()` keeps those dimnames on `coef()`,
#'   `fitted()` and `residuals()`. `deconvolute_ratios()` labels
#'   proportion columns from \(\mu\).
#'
#' @srrstats {RE1.4} Distributional assumptions (multivariate Gaussian,
#'   positive-definite Sigma, simplex-constrained proportions) documented in
#'   `?fit_decovart`, `?deconvolute_ratios` and vignettes.
#'
#' RE2 - Regression: pre-processing
#'
#' @srrstats {RE2.0} ALR reparametrisation documented in `?fit_decovart`,
#'   `?additive_logistic`, and `vignette("softmax-alr-derivatives")`.
#'
#' @srrstats {RE2.1} Missing values in `y` or the signature raise an error
#'   in `.prepare_deconvolution_inputs()`; no silent imputation.
#'
#' @srrstats {RE2.3} `standardise = TRUE` applies one gene-wise affine
#'   z-score (centre/scale from \(\mu\)) to bulk, means and covariances.
#'   Log2 mixing (`scaled = TRUE`) is rejected (Jensen / CIBERSORT).
#'
#' @srrstats {RE2.4} Condition of \(\Sigma(\boldsymbol{p})\) is handled in
#'   `.sigma_p_factorisation()`; collinear signature columns warn.
#'
#' @srrstats {RE2.4a} Identical or rank-deficient signature columns
#'   trigger `.warn_collinear_signature()`.
#'
#' @srrstats {RE2.4b} Perfect collinearity among predictors (duplicate
#'   columns, or a parent equal to a combination of children) is tested
#'   in `tests/testthat/test-03_05_decovart_fit.R` and the simulation
#'   appendix vignette.
#'
#' RE3 - Regression: convergence
#'
#' @srrstats {RE3.0} Non-convergence of Marquardt--Levenberg raises a
#'   `warning()` with the `istop` / RDM code.
#'
#' @srrstats {RE3.1} Convergence diagnostics are stored on `decovart_fit`
#'   (`$convergence`); warnings remain suppressible with
#'   `suppressWarnings()`.
#'
#' @srrstats {RE3.2} Default `epsilon = 1e-4`, `itmax = 200` documented
#'   in `?fit_decovart` (optim-style).
#'
#' @srrstats {RE3.3} Both knobs are explicit arguments of `fit_decovart()`
#'   and of the three native solvers.
#'
#' RE4 - Regression: return objects
#'
#' @srrstats {RE4.0} `fit_decovart()` returns class `decovart_fit` with
#'   proportions, log-likelihood, Fisher `vcov`, and convergence.
#'   `deconvolute_ratios()` remains a benchmark tibble of many algorithms.
#'
#' @srrstats {RE4.2} `coef(fit)` is \(\hat{\boldsymbol{P}}\) (\(J\times N\)).
#'
#' @srrstats {RE4.3} `confint(fit)` is the ALR delta-method Wald interval.
#'
#' @srrstats {RE4.5} `nobs(fit)` is \(N\) (one MVN observation per bulk
#'   sample), with attributes `n_genes` and `n_celltypes`.
#'
#' @srrstats {RE4.6} `vcov(fit)` is the Cramer--Rao / expected-Fisher
#'   bound mapped to the simplex.
#'
#' @srrstats {RE4.7} Optimiser stop codes and iteration counts are stored
#'   on `$convergence` (this tag is **convergence**, not prediction).
#'
#' @srrstats {RE4.8} Observed bulk `Y` is `$bulk_expression`.
#'
#' @srrstats {RE4.9} `fitted(fit)` is \(\boldsymbol{\mu}\hat{\boldsymbol{P}}\).
#'
#' @srrstats {RE4.10} `residuals(fit)` is \(\boldsymbol{Y}-\hat{\boldsymbol{Y}}\).
#'   These are convolution residuals, not OLS residuals; GoF is the MLE
#'   log-likelihood (RE4.11).
#'
#' @srrstats {RE4.11} `summary(fit)` reports log-likelihood and AIC, not
#'   least-squares \(R^{2}\).
#'
#' @srrstats {RE4.12} ALR maps: `additive_logistic` / `additive_log_ratio`.
#'
#' @srrstats {RE4.13} Signature and `Sigma` are stored on the fit.
#'
#' @srrstats {RE4.17} `print.decovart_fit()` shows proportions and logLik.
#'
#' @srrstats {RE4.18} `summary.decovart_fit()` adds Wald SEs and AIC.
#'
#' RE5 - Regression: scaling
#'
#' @srrstatsTODO {RE5.0} Runtime vs \(G\), \(J\), overlap, CPM / log2
#'   normalisation, and tolerance will be reported in a later benchmark
#'   paper / issue, not in this package release.
#'
#' RE6 - Regression: visualisation
#'
#' @srrstats {RE6.0} `plot.decovart_fit()` compares observed and fitted
#'   bulk expression after \(\hat{\boldsymbol{p}}\) has been inferred.
#'
#' @srrstats {RE6.1} The `plot()` generic dispatches on `decovart_fit`.
#'
#' @srrstats {RE6.2} Default plot is observed vs fitted bulk profiles
#'   with the identity line (not estimated vs true proportions).
#'
#' RE7 - Regression: tests
#'
#' @srrstats {RE7.0} Tests with noiseless exact proportions in
#'   `tests/testthat/test-03_03_*.R` and `test-03_05_decovart_fit.R`.
#'
#' @srrstats {RE7.0a} Rank-deficient signatures warn
#'   (`.warn_collinear_signature()`); \(J>G\) errors.
#'
#' @srrstats {RE7.1} Noiseless tests recover ground-truth proportions
#'   within numerical tolerance.
#'
#' @srrstatsTODO {RE7.1a} Speed comparison noiseless vs noisy fitting,
#'   and scaling with \(G,J\), deferred with RE5.0.
#'
#' @srrstats {RE7.2} Tests confirm dimnames of \(Y\) and \(\mu\) on
#'   `decovart_fit` accessors.
#'
#' @srrstats {RE7.3} Tests exercise `coef`, `fitted`, `residuals`,
#'   `vcov`, `nobs`, `print`, `summary`, `plot`.
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
#' @srrstatsNA {RE1.0} No `formula` interface. DeCovarT is not a linear
#'   model for predicting bulk transcriptomes: \(\boldsymbol{\Sigma}(p)\) is
#'   a variance / likelihood specification, not extra regressors. A
#'   formula / `lm` API would misrepresent the estimator.
#' @srrstatsNA {RE1.1} There is no formula to convert (see RE1.0).
#' @srrstatsNA {RE1.3a} Dimnames that exist on \(Y\) and \(\mu\) are
#'   always transferred; there is no supported path that silently drops
#'   metadata.
#' @srrstatsNA {RE2.2} Predictors and response are complete numeric
#'   matrices with a single missing-data policy (error). Separate NA
#'   handling is not relevant for quantified bulk RNA-seq.
#' @srrstatsNA {RE4.1} An unfitted `lm`-style model object is not
#'   meaningful: the convolution has no formula and is not estimated by
#'   least squares.
#' @srrstatsNA {RE4.4} `formula()` is not implemented (see RE1.0).
#' @srrstatsNA {RE4.14} Forecast / extrapolation errors do not apply:
#'   the target is a simplex of mixing weights, not a future bulk
#'   expression.
#' @srrstatsNA {RE4.15} There is no forecast horizon.
#' @srrstatsNA {RE4.16} `predict()` is not implemented: DeCovarT does
#'   not forecast bulk transcriptomic values (new groups, interpolation
#'   or extrapolation).
#' @srrstatsNA {RE6.3} Interpolated vs extrapolated predictions do not
#'   arise (no `predict()`).
#' @srrstatsNA {RE7.4} Forecast-horizon tests do not apply (see RE4.15).
#'
#' @noRd
NULL
