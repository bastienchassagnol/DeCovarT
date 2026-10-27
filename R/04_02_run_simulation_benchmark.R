#' Run one generative scenario: simulate bulk mixtures and deconvolve
#'
#' @param true_theta Named list with `mu`, `sigma` or `Theta`, and `p`.
#' @param n_samples Number of bulk replicates \eqn{N}.
#' @param deconvolution_functions Named list for [deconvolute_ratios()].
#' @param standardise,scaled Passed to [deconvolute_ratios()].
#' @param cores Workers for the per-sample loop inside [deconvolute_ratios()].
#' @param scenario_meta Tibble row of scenario metadata (excluding `true_theta`
#'   and per-row `n`).
#' @param model,reference,reference_seed,return_data See
#'   [run_simulation_benchmark()].
#'
#' @return List with `regression`, `monte_carlo`, `optimisation`,
#'   `intervals` and `config` tibbles, and `data` when `return_data` is
#'   `TRUE`.
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
  progress_every = 10L,
  model = "population",
  reference = NULL,
  reference_seed = NULL,
  return_data = FALSE
) {
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (is.null(Sigma)) {
    .ui_abort("{.arg true_theta} must contain {.field sigma} for simulation.")
  }
  p <- true_theta$p
  generative_mean <- if (is.null(true_theta$generative_mean)) {
    mu
  } else {
    true_theta$generative_mean
  }
  generative_sigma <- if (is.null(true_theta$generative_sigma)) {
    Sigma
  } else {
    true_theta$generative_sigma
  }

  described <- describe_simulation_scenario(
    true_theta = true_theta,
    adjacency = true_theta$adjacency
  )

  nuisance_true <- NULL
  if (identical(model, "single_cell")) {
    simulated_data <- simulate_bulk_from_cells(
      generative_mean,
      generative_sigma,
      p = p,
      n_cells = true_theta$n_cells,
      platform_scale = if (is.null(true_theta$platform_scale)) {
        1
      } else {
        true_theta$platform_scale
      },
      n_samples = n_samples
    )
    # The target of a single-cell bulk is its realised composition.
    true_ratios <- simulated_data$p_realised
    nuisance_true <- rbind(
      mean_scale = simulated_data$mean_scale,
      variance_scale = simulated_data$variance_scale
    )
  } else {
    simulated_data <- simulate_bulk_mixture(
      signature_matrix = mu,
      Sigma = Sigma,
      p = p,
      n = n_samples,
      truncate_negative = FALSE
    )
    true_ratios <- p
  }

  estimated_ratios <- suppressWarnings(deconvolute_ratios(
    signature_matrix = mu,
    bulk_expression = simulated_data$Y,
    true_ratios = true_ratios,
    Sigma = Sigma,
    deconvolution_functions = deconvolution_functions,
    standardise = standardise,
    scaled = scaled,
    cores = cores,
    verbose = verbose,
    progress_every = progress_every,
    coverage_interval = coverage_interval,
    model = model,
    nuisance_true = nuisance_true
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
    intervals = .attach_scenario_meta(estimated_ratios$intervals),
    config = config,
    theta_true = described$theta_true,
    descriptors = .attach_scenario_meta(described$descriptors),
    supplementary = .attach_scenario_meta(described$supplementary),
    data = if (isTRUE(return_data)) {
      list(
        bulk = simulated_data$Y,
        truth = .as_p_matrix(true_ratios, colnames(mu), n_samples),
        nuisance = nuisance_true,
        reference = if (!is.null(reference)) {
          simulate_reference_replicates(
            generative_mean,
            generative_sigma,
            n_individuals = reference$n_individuals,
            n_cells = if (is.null(reference$n_cells)) 1L else reference$n_cells,
            seed = reference_seed
          )
        }
      )
    }
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
#' @param model `"population"` (default) draws every bulk from the
#'   convolution [simulate_bulk_mixture()]; `"single_cell"` draws it as a
#'   sum of cells with [simulate_bulk_from_cells()], using `n_cells` and
#'   `platform_scale` from `true_theta` (and `generative_mean` /
#'   `generative_sigma` when the fitted reference differs from the
#'   generator), and scores every estimate against its realised
#'   composition.
#' @param return_data If `TRUE`, keep each scenario's bulk matrix,
#'   targets and true scales in `data`, so methods run outside this
#'   function can be scored on the same samples.
#' @param reference Optional list with `n_individuals`, `n_cells` (default
#'   `1`) and `seed`: with `return_data = TRUE`, an independent reference
#'   panel is drawn per scenario by [simulate_reference_replicates()], with
#'   seed `reference$seed` plus the scenario index, for methods that
#'   estimate the reference moments themselves.
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
#' * `intervals`: coverage and width of the interval methods requested in
#'   `deconvolution_functions`;
#' * `data`: per-scenario simulated data when `return_data = TRUE`, else
#'   `NULL`;
#' * `call`: the matched call ([match.call()]).
#' Analysis scripts persist those pieces with `write_simulation_artefacts()`
#' (`scripts/benchmark_metrics/benchmark_artefacts.R`) rather than saving the whole
#' list.
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
  progress_every = 10L,
  model = c("population", "single_cell"),
  return_data = FALSE,
  reference = NULL
) {
  model <- match.arg(model)
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
      "Benchmark: {.val {n_scen}} scenarios \u00b7 {.val {n_algo}} algorithms \u00b7 n = {.val {n}}."
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
        "Scenario {.val {i}}/{.val {n_scen}} \u00b7 {label} \u00b7 n = {.val {n_samples}}."
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
      progress_every = progress_every,
      model = model,
      reference = reference,
      reference_seed = if (is.null(reference$seed)) {
        NULL
      } else {
        reference$seed + i
      },
      return_data = return_data
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
    intervals = dplyr::bind_rows(purrr::map(scenario_results, "intervals")),
    data = if (isTRUE(return_data)) purrr::map(scenario_results, "data"),
    call = call
  )
}
