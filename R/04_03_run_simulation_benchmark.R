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
  scenario_meta = NULL,
  coverage_interval = "wilson",
  verbose = FALSE,
  progress_every = 10L
) {
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (is.null(Sigma)) {
    .ui_abort("{.arg true_theta} must contain {.field sigma} for simulation.")
  }
  p <- true_theta$p

  described <- describe_simulation_scenario(
    true_theta = true_theta,
    adjacency = true_theta$adjacency
  )

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
    cores = cores,
    verbose = verbose,
    progress_every = progress_every,
    coverage_interval = coverage_interval
  ))

  .attach_scenario_meta <- function(tbl) {
    if (
      nrow(tbl) > 0L &&
        !is.null(scenario_meta) &&
        ncol(scenario_meta) > 0L
    ) {
      extra <- setdiff(names(scenario_meta), names(tbl))
      if (length(extra) > 0L) {
        meta_rep <- scenario_meta[
          rep(1L, nrow(tbl)),
          extra,
          drop = FALSE
        ]
        tbl <- dplyr::bind_cols(meta_rep, tbl)
      }
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
    config = config,
    theta_true = described$theta_true,
    descriptors = .attach_scenario_meta(described$descriptors),
    supplementary = .attach_scenario_meta(described$supplementary)
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
#' `scripts/fig02_bivariate_toy.R` and the paper-scenario vignettes.
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
#' @param coverage_interval Coverage interval for the Monte Carlo
#'   coverage *rate*; see [coverage_mc_interval()].
#' @param verbose If `TRUE`, print each scenario row and (when the grid
#'   has at most 10 scenarios) every `progress_every` inferred samples.
#'   Large factorial grids log one line per scenario only, so logs stay
#'   readable.
#' @param progress_every Sample-progress interval passed to
#'   [deconvolute_ratios()] when `verbose` is `TRUE`. Defaults to `10L`.
#'
#' @return A list with:
#' * `regression`: `global` (per-sample composition scores) and
#'   `cell_type` (across-sample Pearson / F1 / spillover);
#' * `monte_carlo`: ADEMP summaries per cell type;
#' * `optimisation`: per-sample elapsed time, memory, KKT residual, and
#'   \eqn{\hat{\boldsymbol{p}}};
#' * `config`: tibble of scenario metadata (one row per scenario);
#' * `theta_true`: list of convolution parameters
#'   (\eqn{\boldsymbol{p}}, \eqn{\boldsymbol{\mu}}, \eqn{\boldsymbol{\Sigma}_j})
#'   actually used to draw the bulk;
#' * `descriptors`: kept scenario statistics (composition, mean, SPD,
#'   network, tangent Fisher, MixSim BarOmega, pairwise Hellinger);
#' * `supplementary`: Jeffreys overlap, recorded separately;
#' * `call`: the matched call ([match.call()]).
#' There is no composite global score: each metric answers a different
#' question.
#'
#' @examplesIf .Platform$OS.type != "windows"
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
#'   [compute_benchmark_metrics()], [describe_simulation_scenario()],
#'   [coverage_mc_interval()], [plot_mc_raincloud()], [plot_mc_forest()]
run_simulation_benchmark <- function(
  scenario_config,
  deconvolution_functions,
  n = 200,
  standardise = FALSE,
  scaled = FALSE,
  cores = 1L,
  coverage_interval = "wilson",
  verbose = FALSE,
  progress_every = 10L
) {
  call <- match.call()
  if (!is.data.frame(scenario_config)) {
    if (!is.list(scenario_config)) {
      .ui_abort(
        "{.arg scenario_config} must be a tibble or list of scenario rows."
      )
    }
    scenario_config <- dplyr::bind_rows(scenario_config)
  }
  if (!"true_theta" %in% names(scenario_config)) {
    .ui_abort(
      "{.arg scenario_config} must contain a {.field true_theta} column."
    )
  }
  if (nrow(scenario_config) == 0L) {
    .ui_abort("{.arg scenario_config} must have at least one row.")
  }

  n_scen <- nrow(scenario_config)
  n_algo <- length(deconvolution_functions)
  sample_verbose <- isTRUE(verbose) && n_scen <= 10L
  if (isTRUE(verbose)) {
    .ui_info(
      "Benchmark: {.val {n_scen}} scenarios · {.val {n_algo}} algorithms · n = {.val {n}}."
    )
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
      .ui_abort(
        "{.arg true_theta} entries must be lists with {.field mu} and {.field sigma}."
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
    if (isTRUE(verbose)) {
      label <- .ui_scenario_label(row)
      .ui_info(
        "Scenario {.val {i}}/{.val {n_scen}} · {label} · n = {.val {n_samples}}."
      )
    }
    .run_one_simulation_scenario(
      true_theta = true_theta,
      n_samples = n_samples,
      deconvolution_functions = deconvolution_functions,
      standardise = standardise,
      scaled = scaled,
      cores = cores,
      scenario_meta = scenario_meta,
      coverage_interval = coverage_interval,
      verbose = sample_verbose,
      progress_every = progress_every
    )
  }

  scenario_results <- lapply(seq_len(nrow(scenario_config)), .run_row)
  if (isTRUE(verbose)) {
    .ui_success("Benchmark complete ({.val {n_scen}} scenarios).")
  }

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
    config = dplyr::bind_rows(purrr::map(scenario_results, "config")),
    theta_true = purrr::map(scenario_results, "theta_true"),
    descriptors = dplyr::bind_rows(
      purrr::map(scenario_results, "descriptors")
    ),
    supplementary = dplyr::bind_rows(
      purrr::map(scenario_results, "supplementary")
    ),
    call = call
  )
}
