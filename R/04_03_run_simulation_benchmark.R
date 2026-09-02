#' Run one generative scenario: simulate bulk mixtures and deconvolve
#'
#' @param true_theta Named list with `mu`, `sigma` or `Theta`, and `p`.
#' @param n_samples Number of bulk replicates \eqn{N}.
#' @param deconvolution_functions Named list for [deconvolute_ratios()].
#' @param standardise,scaled Passed to [deconvolute_ratios()].
#' @param cores Workers for the per-sample loop inside [deconvolute_ratios()].
#' @param scenario_meta Tibble row of scenario metadata (excluding `true_theta`
#'   and per-row `n`).
#'
#' @return List with `regression`, `monte_carlo`, `optimisation`, and
#'   `config` tibbles.
#' @keywords internal
#' @noRd
.run_one_simulation_scenario <- function(
  true_theta,
  n_samples,
  deconvolution_functions,
  standardise = FALSE,
  scaled = FALSE,
  cores = 1L,
  scenario_meta = NULL
) {
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (is.null(Sigma)) {
    stop("`true_theta` must contain `sigma` for simulation.", call. = FALSE)
  }
  p <- true_theta$p

  simulated_data <- simulate_bulk_mixture(
    signature_matrix = mu,
    Sigma = Sigma,
    p = p,
    n = n_samples
  )

  estimated_ratios <- suppressWarnings(deconvolute_ratios(
    signature_matrix = mu,
    bulk_expression = simulated_data$Y,
    true_ratios = p,
    Sigma = Sigma,
    deconvolution_functions = deconvolution_functions,
    standardise = standardise,
    scaled = scaled,
    cores = cores
  ))

  .attach_scenario_meta <- function(tbl) {
    if (
      nrow(tbl) > 0L &&
        !is.null(scenario_meta) &&
        ncol(scenario_meta) > 0L
    ) {
      meta_rep <- scenario_meta[rep(1L, nrow(tbl)), , drop = FALSE]
      tbl <- dplyr::bind_cols(meta_rep, tbl)
    }
    tbl
  }

  config <- scenario_meta
  if (is.null(config)) {
    config <- tibble::tibble()
  }
  if (nrow(config) == 0L) {
    config <- tibble::tibble(true_parameters = list(as.list(true_theta)))
  } else {
    config$true_parameters <- list(as.list(true_theta))
  }
  config$nobservations <- n_samples

  list(
    regression = list(
      global = .attach_scenario_meta(estimated_ratios$regression$global),
      cell_type = .attach_scenario_meta(
        estimated_ratios$regression$cell_type
      )
    ),
    monte_carlo = .attach_scenario_meta(estimated_ratios$monte_carlo),
    optimisation = .attach_scenario_meta(estimated_ratios$optimisation),
    config = config
  )
}

#' Simulate bulk mixtures and benchmark deconvolution algorithms
#'
#' @description
#' Wrapper around [simulate_bulk_mixture()], [deconvolute_ratios()], and
#' [compute_benchmark_metrics()] (called inside `deconvolute_ratios()`).
#' Each row of `scenario_config` describes one generative model
#' (\eqn{\boldsymbol{\mu}}, \eqn{(\boldsymbol{\Sigma}_j)_j},
#' \eqn{\boldsymbol{p}}) stored in a list column `true_theta`. Scenario
#' rows are always evaluated **sequentially** to avoid nested
#' parallelism; sample-level workers live only in
#' [deconvolute_ratios()]. Scenario builders (factorial grids, overlap
#' summaries, etc.) should live in analysis scripts; see
#' `scripts/configure_bivariate_toy_scenarios.R` and the manuscript
#' scenario vignettes.
#'
#' @param scenario_config Tibble or list of scenario rows. Each row must
#'   contain a `true_theta` list column (or list element) with at least
#'   `mu` and `sigma`. Optional per-row `n` overrides the default `n`.
#' @param deconvolution_functions Named list passed to
#'   [deconvolute_ratios()]; each element has `FUN` and optional
#'   `additional_parameters` for `do.call()`.
#' @param n Default number of bulk replicates \eqn{N} when `scenario_config`
#'   has no `n` column.
#' @param standardise,scaled Passed to [deconvolute_ratios()].
#' @param cores Workers for the per-sample loop inside
#'   [deconvolute_ratios()]. Defaults to `1L`.
#'
#' @return A list with:
#' * `regression`: `global` (per-sample composition scores) and
#'   `cell_type` (across-sample Pearson / F1 / spillover);
#' * `monte_carlo`: ADEMP summaries per cell type;
#' * `optimisation`: per-sample elapsed time, memory, KKT residual, and
#'   \eqn{\hat{\boldsymbol{p}}};
#' * `config`: tibble of scenario metadata (one row per scenario).
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' scenario_config <- tibble::tibble(
#'   ID = "B1_Ho",
#'   true_theta = list(list(
#'     p = c(0.5, 0.5),
#'     mu = mu,
#'     sigma = Sigma
#'   ))
#' )
#' out <- run_simulation_benchmark(
#'   scenario_config = scenario_config,
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   n = 2,
#'   cores = 1
#' )
#' nrow(out$optimisation)
#' @importFrom rlang .data
#' @export
#' @seealso [simulate_bulk_mixture()], [deconvolute_ratios()],
#'   [compute_benchmark_metrics()]
run_simulation_benchmark <- function(
  scenario_config,
  deconvolution_functions,
  n = 200,
  standardise = FALSE,
  scaled = FALSE,
  cores = 1L
) {
  if (!is.data.frame(scenario_config)) {
    if (!is.list(scenario_config)) {
      stop(
        "`scenario_config` must be a tibble or list of scenario rows.",
        call. = FALSE
      )
    }
    scenario_config <- dplyr::bind_rows(scenario_config)
  }
  if (!"true_theta" %in% names(scenario_config)) {
    stop("`scenario_config` must contain a `true_theta` column.", call. = FALSE)
  }
  if (nrow(scenario_config) == 0L) {
    stop("`scenario_config` must have at least one row.", call. = FALSE)
  }

  has_n_col <- "n" %in% names(scenario_config)
  meta_cols <- setdiff(
    names(scenario_config),
    c("true_theta", if (has_n_col) "n")
  )

  .run_row <- function(i) {
    row <- scenario_config[i, , drop = FALSE]
    true_theta <- row$true_theta[[1L]]
    if (!is.list(true_theta)) {
      stop(
        "`true_theta` entries must be lists with `mu` and `sigma`.",
        call. = FALSE
      )
    }
    n_samples <- if (has_n_col && !is.na(row$n[1L])) {
      as.integer(row$n[1L])
    } else {
      as.integer(n)
    }
    scenario_meta <- if (length(meta_cols) > 0L) {
      row[, meta_cols, drop = FALSE]
    } else {
      NULL
    }
    .run_one_simulation_scenario(
      true_theta = true_theta,
      n_samples = n_samples,
      deconvolution_functions = deconvolution_functions,
      standardise = standardise,
      scaled = scaled,
      cores = cores,
      scenario_meta = scenario_meta
    )
  }

  scenario_results <- lapply(seq_len(nrow(scenario_config)), .run_row)

  list(
    regression = list(
      global = dplyr::bind_rows(
        purrr::map(scenario_results, \(x) x$regression$global)
      ),
      cell_type = dplyr::bind_rows(
        purrr::map(scenario_results, \(x) x$regression$cell_type)
      )
    ),
    monte_carlo = dplyr::bind_rows(
      purrr::map(scenario_results, "monte_carlo")
    ),
    optimisation = dplyr::bind_rows(
      purrr::map(scenario_results, "optimisation")
    ),
    config = dplyr::bind_rows(purrr::map(scenario_results, "config"))
  )
}
