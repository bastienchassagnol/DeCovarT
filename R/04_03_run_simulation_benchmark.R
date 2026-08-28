#' Default parallel worker count (half of detected cores)
#'
#' @return Integer in \eqn{[1,\lfloor C/2\rfloor]}.
#' @keywords internal
#' @noRd
.default_parallel_cores <- function() {
  detected <- parallel::detectCores()
  as.integer(max(1L, floor(detected / 2L)))
}

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
#' @return List with `simulations` and `config` tibbles.
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

  simulations <- estimated_ratios
  if (
    nrow(simulations) > 0L &&
      !is.null(scenario_meta) &&
      ncol(scenario_meta) > 0L
  ) {
    meta_rep <- scenario_meta[rep(1L, nrow(simulations)), , drop = FALSE]
    simulations <- dplyr::bind_cols(meta_rep, simulations)
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

  list(simulations = simulations, config = config)
}

#' Simulate bulk mixtures and benchmark deconvolution algorithms
#'
#' @description
#' Wrapper around [simulate_bulk_mixture()], [deconvolute_ratios()], and
#' [compute_benchmark_metrics()] (called inside `deconvolute_ratios()`).
#' Each row of `scenario_config` describes one generative model
#' (\eqn{\boldsymbol{\mu}}, \eqn{(\boldsymbol{\Sigma}_j)_j},
#' \eqn{\boldsymbol{p}}) stored in a list column `true_theta`. Scenario
#' builders (factorial grids, overlap summaries, etc.) should live in
#' analysis scripts; see `scripts/configure_bivariate_toy_scenarios.R` and
#' the manuscript scenario vignettes.
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
#'   [deconvolute_ratios()]. When `parallel_scenarios = TRUE`, defaults to
#'   `1` to avoid nested parallelism. Otherwise defaults to half of detected
#'   cores (at most \eqn{\lfloor C/2\rfloor}).
#' @param parallel_scenarios If `TRUE`, parallelise across scenario rows
#'   with `furrr` (optional Suggests `furrr` and `future`). Defaults to
#'   `FALSE`.
#' @param parallel_cores Maximum workers for scenario-level parallelism;
#'   defaults to half of detected cores.
#'
#' @return A list with:
#' * `simulations`: tibble of per-sample estimates and metrics;
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
#' nrow(out$simulations)
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
  cores = NULL,
  parallel_scenarios = FALSE,
  parallel_cores = NULL
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

  parallel_cores <- if (is.null(parallel_cores)) {
    .default_parallel_cores()
  } else {
    parallel_cores
  }
  if (isTRUE(parallel_scenarios)) {
    cores <- if (is.null(cores)) 1L else cores
  } else {
    cores <- if (is.null(cores)) parallel_cores else cores
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

  if (isTRUE(parallel_scenarios)) {
    .check_suggested_package("furrr", "run_simulation_benchmark")
    .check_suggested_package("future", "run_simulation_benchmark")
    old_plan <- future::plan(
      future::multisession,
      workers = parallel_cores
    )
    on.exit(future::plan(old_plan), add = TRUE)
    scenario_results <- furrr::future_map(
      seq_len(nrow(scenario_config)),
      .run_row,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    scenario_results <- lapply(seq_len(nrow(scenario_config)), .run_row)
  }

  list(
    simulations = dplyr::bind_rows(
      purrr::map(scenario_results, "simulations")
    ),
    config = dplyr::bind_rows(purrr::map(scenario_results, "config"))
  )
}
