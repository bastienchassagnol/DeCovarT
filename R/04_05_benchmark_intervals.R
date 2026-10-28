# Coverage table of the interval menu: turns the per-bulk interval
# matrices of .likelihood_intervals() (04_00_likelihoods.R), or intervals
# computed outside the package, into coverage and width per method and
# parameter. It does not depend on the model.

#' Coverage and width of every interval method over the Monte Carlo samples
#'
#' @param intervals List over samples of the output of
#'   [.likelihood_intervals()]: per sample, a named list of two-column
#'   matrices (lower, upper), one per method.
#' @param p_hat,p_true \eqn{J\times N} estimates and targets.
#' @param nuisance_hat,nuisance_true Optional \eqn{2\times N} matrices of
#'   the estimated and true scales (`mean_scale`, `variance_scale`); the
#'   platform scale \eqn{a=\phi/c} and the cell count \eqn{n=c^{2}/\phi}
#'   are derived from them and scored as well.
#' @return A tibble with one row per algorithm, interval method and
#'   parameter. `n_intervals` counts finite bounds.
#'   `n_numerical_failures` counts samples whose estimate for that
#'   parameter is not finite (the optimiser returned no MLE). A finite
#'   estimate with no interval row is neither: it is a missing interval.
#'   `n_boundary` counts samples whose estimate has a zero proportion (face,
#'   edge or vertex), where Wald, ILR and oracle intervals are undefined.
#' @noRd
.interval_table <- function(
  intervals,
  p_hat,
  p_true,
  nuisance_hat = NULL,
  nuisance_true = NULL,
  coverage_interval = "wilson",
  algorithm = NA_character_
) {
  if (length(intervals) == 0L || all(vapply(intervals, is.null, TRUE))) {
    return(tibble::tibble())
  }
  derive <- function(scales) {
    if (is.null(scales)) NULL else .derived_scales(scales)
  }
  estimate <- rbind(p_hat, derive(nuisance_hat))
  truth <- rbind(p_true, derive(nuisance_true))
  methods <- unique(unlist(lapply(intervals, names)))
  composition <- rownames(p_hat)
  n_boundary <- sum(colSums(p_hat <= 0, na.rm = TRUE) > 0)
  purrr::map_dfr(methods, function(method) {
    # Same rows for every bulk; a bulk without a row has an NA interval.
    parameters <- intersect(
      rownames(truth),
      unique(unlist(lapply(intervals, function(x) rownames(x[[method]]))))
    )
    padded <- lapply(intervals, function(x) .pad_bounds(x[[method]], parameters))
    purrr::map_dfr(parameters, function(parameter) {
      bound <- function(side) {
        vapply(padded, function(m) m[parameter, side], numeric(1))
      }
      lower <- bound("lower")
      upper <- bound("upper")
      target <- truth[parameter, ]
      estimate_i <- estimate[parameter, ]
      numerical_failure <- !is.finite(estimate_i)
      finite <- is.finite(lower) & is.finite(upper)
      width <- upper - lower
      rate <- coverage_mc_interval(
        (lower <= target & target <= upper)[finite],
        method = coverage_interval
      )
      tibble::tibble(
        algorithm = algorithm,
        ci_method = method,
        parameter = parameter,
        is_composition = parameter %in% composition,
        n_numerical_failures = sum(numerical_failure),
        n_boundary = n_boundary,
        n_intervals = sum(finite),
        coverage = rate$coverage,
        coverage_lower = rate$lower,
        coverage_upper = rate$upper,
        mean_width = if (any(finite)) mean(width[finite]) else NA_real_,
        median_width = if (any(finite)) {
          stats::median(width[finite])
        } else {
          NA_real_
        },
        mean_lower = if (any(finite)) mean(lower[finite]) else NA_real_,
        mean_upper = if (any(finite)) mean(upper[finite]) else NA_real_,
        mean_estimate = mean(estimate[parameter, ], na.rm = TRUE),
        mean_truth = mean(target, na.rm = TRUE),
        outside_unit = if (parameter %in% composition && any(finite)) {
          mean((lower < 0 | upper > 1)[finite])
        } else {
          NA_real_
        }
      )
    })
  })
}
