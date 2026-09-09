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
    n = n_samples,
    truncate_negative = FALSE
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
#' Scripts should persist those pieces with
#' [write_simulation_artefacts()] rather than saving the whole list.
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
#'   [coverage_mc_interval()], [plot_mc_raincloud()], [plot_mc_forest()],
#'   [write_simulation_artefacts()]
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
      progress_every = progress_every
    )
  }

  scenario_results <- lapply(seq_len(nrow(scenario_config)), .run_row)
  if (isTRUE(verbose)) {
    .ui_success("Benchmark complete ({.val {n_scen}} scenarios).")
  }

  list(
    regression = list(
      global = .bind_metrics_rows(
        purrr::map(scenario_results, \(x) x$regression$global)
      ),
      cell_type = .bind_metrics_rows(
        purrr::map(scenario_results, \(x) x$regression$cell_type)
      )
    ),
    monte_carlo = .bind_metrics_rows(
      purrr::map(scenario_results, "monte_carlo")
    ),
    optimisation = .bind_metrics_rows(
      purrr::map(scenario_results, "optimisation")
    ),
    config = .bind_metrics_rows(purrr::map(scenario_results, "config")),
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

#' Redundant design columns dropped from slim scenario tables
#'
#' @keywords internal
#' @noRd
.redundant_scenario_columns <- function() {
  c(
    "scenario_idx",
    "rho_ct1",
    "rho_ct2",
    "proportion_name",
    "centroid"
  )
}

#' Columns that belong on the descriptors artefact (not the design grid)
#'
#' @keywords internal
#' @noRd
.descriptor_metric_columns <- function() {
  c(
    "n_genes",
    "n_celltypes",
    "h_star",
    "n_eff",
    "n_active",
    "min_active",
    "concentration",
    "mean_abs_cosine",
    "min_cosine",
    "max_cosine",
    "mean_euclidean",
    "kappa_mu",
    "gram_volume",
    "lambda_min_sigma_p",
    "kappa_sigma_p",
    "kappa_sigma_reciprocal",
    "lambda_min_it",
    "kappa_it",
    "f_cov",
    "f_cov_max",
    "network_density",
    "network_mean_degree",
    "hoyer_abs_correlation",
    "mixsim_baromega",
    "hellinger",
    "hellinger_weighted",
    "jeffreys"
  )
}

#' Encode a bivariate-toy scenario ID
#'
#' Pattern `B{index}_{Ho|He}_{Ba|Mo|Hi}_{Sm|Lg}`.
#'
#' @keywords internal
#' @noRd
.encode_bivariate_id <- function(
  scenario_idx,
  variance,
  proportions,
  centroids
) {
  var_code <- ifelse(variance == "homoscedastic", "Ho", "He")
  p_code <- dplyr::case_when(
    proportions == "balanced" ~ "Ba",
    proportions == "moderately unbalanced" ~ "Mo",
    proportions == "highly unbalanced" ~ "Hi",
    TRUE ~ "Xx"
  )
  cld_code <- ifelse(
    grepl("small", centroids, ignore.case = TRUE),
    "Sm",
    "Lg"
  )
  paste0(
    "B",
    as.integer(scenario_idx),
    "_",
    var_code,
    "_",
    p_code,
    "_",
    cld_code
  )
}

#' Drop duplicated design columns from a scenario-tagged table
#'
#' Aliases (`scenario_idx`, `rho_ct1` / `rho_ct2`, `proportion_name`,
#' `centroid`) are removed after copying them onto the canonical names
#' `proportions` and `centroids` when those are missing.
#'
#' @param tbl A tibble that may contain redundant aliases.
#' @return The same table without alias columns.
#' @export
slim_scenario_table <- function(tbl) {
  tbl <- tibble::as_tibble(tbl)
  if ("proportion_name" %in% names(tbl) && !"proportions" %in% names(tbl)) {
    tbl$proportions <- tbl$proportion_name
  }
  if ("centroid" %in% names(tbl) && !"centroids" %in% names(tbl)) {
    tbl$centroids <- tbl$centroid
  }
  drop <- intersect(.redundant_scenario_columns(), names(tbl))
  if (length(drop) > 0L) {
    tbl <- tbl[, setdiff(names(tbl), drop), drop = FALSE]
  }
  tbl
}

#' Unwrap a list-column `true_theta` cell
#'
#' @keywords internal
#' @noRd
.unwrap_true_theta <- function(th) {
  if (!is.list(th)) {
    return(th)
  }
  if (!is.null(th$p) && !is.null(th$mu)) {
    return(th)
  }
  if (length(th) >= 1L && is.list(th[[1L]])) {
    return(.unwrap_true_theta(th[[1L]]))
  }
  th
}

#' Guarantee an `ID` column on a scenario table
#'
#' @keywords internal
#' @noRd
.ensure_scenario_id <- function(tbl, prefix = "S") {
  if (!"ID" %in% names(tbl) || any(!nzchar(as.character(tbl$ID)))) {
    tbl$ID <- paste0(prefix, seq_len(nrow(tbl)))
  }
  tbl$ID <- as.character(tbl$ID)
  tbl
}

#' Rewrite bivariate `ID` from design columns
#'
#' @keywords internal
#' @noRd
.rewrite_bivariate_ids <- function(tbl) {
  needed <- c("scenario_idx", "variance", "proportions", "centroids")
  if (!all(needed %in% names(tbl))) {
    return(tbl)
  }
  tbl$ID <- .encode_bivariate_id(
    tbl$scenario_idx,
    tbl$variance,
    tbl$proportions,
    tbl$centroids
  )
  tbl
}

#' Rebuild descriptor metrics from `\theta`, keeping MixSim if present
#'
#' @keywords internal
#' @noRd
.descriptors_from_theta <- function(theta_tbl, old_descriptors = NULL) {
  fresh <- purrr::pmap_dfr(
    theta_tbl,
    function(ID, true_theta, ...) {
      th <- .unwrap_true_theta(true_theta)
      described <- describe_simulation_scenario(
        th,
        include_mixsim = FALSE
      )
      dplyr::bind_cols(
        tibble::tibble(ID = as.character(ID)),
        described$descriptors,
        described$supplementary
      )
    }
  )
  if (
    !is.null(old_descriptors) &&
      "mixsim_baromega" %in% names(old_descriptors) &&
      "ID" %in% names(old_descriptors)
  ) {
    old <- old_descriptors[, c("ID", "mixsim_baromega"), drop = FALSE]
    old$ID <- as.character(old$ID)
    fresh$mixsim_baromega <- NULL
    fresh <- dplyr::left_join(fresh, old, by = "ID")
  }
  fresh
}

#' Write split simulation artefacts (config, descriptors, theta, metrics)
#'
#' [run_simulation_benchmark()] still returns a single list for tests.
#' Scripts persist four RDS files keyed by `ID` so the metrics object
#' does not duplicate geometry, MixSim, or `\theta`.
#'
#' @param benchmark List from [run_simulation_benchmark()].
#' @param dir Output directory.
#' @param stem File-name stem (`bivariate`, `hybrid`, …).
#' @param config Optional scenario grid (defaults to `benchmark$config`).
#' @param rewrite_bivariate_id If `TRUE`, rebuild fig02-style IDs.
#'
#' @return Invisibly, a named list of written paths.
#' @seealso [read_simulation_artefacts()]
#' @export
write_simulation_artefacts <- function(
  benchmark,
  dir,
  stem,
  config = NULL,
  rewrite_bivariate_id = FALSE
) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(config)) {
    config <- benchmark$config
  }
  config <- tibble::as_tibble(config)
  if (isTRUE(rewrite_bivariate_id)) {
    config <- .rewrite_bivariate_ids(config)
  }
  config <- .ensure_scenario_id(config)

  theta_src <- if ("true_theta" %in% names(config)) {
    config$true_theta
  } else {
    benchmark$theta_true
  }
  theta_tbl <- tibble::tibble(
    ID = as.character(config$ID),
    true_theta = lapply(theta_src, .unwrap_true_theta)
  )

  old_desc <- benchmark$descriptors
  if (!is.null(old_desc) && nrow(old_desc) > 0L) {
    if (isTRUE(rewrite_bivariate_id)) {
      old_desc <- .rewrite_bivariate_ids(old_desc)
    }
    old_desc$ID <- as.character(old_desc$ID)
  }
  descriptors <- .descriptors_from_theta(theta_tbl, old_desc)
  keep_desc <- intersect(
    c("ID", .descriptor_metric_columns()),
    names(descriptors)
  )
  descriptors <- descriptors[, keep_desc, drop = FALSE]

  config_slim <- slim_scenario_table(config)
  if ("true_theta" %in% names(config_slim)) {
    config_slim$true_theta <- NULL
  }

  slim_metrics <- function(tbl) {
    if (is.null(tbl) || ncol(tbl) == 0L) {
      return(tbl)
    }
    tbl <- tibble::as_tibble(tbl)
    if (isTRUE(rewrite_bivariate_id)) {
      tbl <- .rewrite_bivariate_ids(tbl)
    }
    design_drop <- setdiff(
      names(config_slim),
      c("ID")
    )
    drop <- unique(
      c(.redundant_scenario_columns(), intersect(design_drop, names(tbl)))
    )
    drop <- setdiff(drop, "ID")
    if (length(drop) > 0L) {
      tbl <- tbl[, setdiff(names(tbl), drop), drop = FALSE]
    }
    tbl
  }

  metrics <- list(
    regression = list(
      global = slim_metrics(benchmark$regression$global),
      cell_type = slim_metrics(benchmark$regression$cell_type)
    ),
    monte_carlo = slim_metrics(benchmark$monte_carlo),
    optimisation = slim_metrics(benchmark$optimisation),
    call = benchmark$call
  )

  paths <- list(
    config = file.path(dir, paste0(stem, "_config.rds")),
    descriptors = file.path(dir, paste0(stem, "_descriptors.rds")),
    theta = file.path(dir, paste0(stem, "_theta.rds")),
    benchmark = file.path(dir, paste0(stem, "_benchmark.rds"))
  )
  saveRDS(config_slim, paths$config)
  saveRDS(descriptors, paths$descriptors)
  saveRDS(theta_tbl, paths$theta)
  saveRDS(metrics, paths$benchmark)
  invisible(paths)
}

#' Read split simulation artefacts and optionally reassemble a benchmark list
#'
#' @inheritParams write_simulation_artefacts
#' @param assemble If `TRUE`, join config onto metric tables for plotting
#'   helpers that still expect design columns.
#'
#' @return Named list of tibbles, or a benchmark-like list when
#'   `assemble = TRUE`.
#' @export
read_simulation_artefacts <- function(dir, stem, assemble = FALSE) {
  paths <- list(
    config = file.path(dir, paste0(stem, "_config.rds")),
    descriptors = file.path(dir, paste0(stem, "_descriptors.rds")),
    theta = file.path(dir, paste0(stem, "_theta.rds")),
    benchmark = file.path(dir, paste0(stem, "_benchmark.rds"))
  )
  out <- lapply(paths, function(p) {
    if (file.exists(p)) {
      readRDS(p)
    } else {
      NULL
    }
  })
  names(out) <- names(paths)
  if (!isTRUE(assemble)) {
    return(out)
  }
  cfg <- out$config
  metrics <- out$benchmark
  join_cfg <- function(tbl) {
    if (is.null(tbl) || is.null(cfg) || !"ID" %in% names(tbl)) {
      return(tbl)
    }
    extra <- setdiff(names(cfg), names(tbl))
    if (length(extra) == 0L) {
      return(tbl)
    }
    dplyr::left_join(tbl, cfg[, c("ID", extra), drop = FALSE], by = "ID")
  }
  theta_list <- if (!is.null(out$theta)) {
    lapply(out$theta$true_theta, .unwrap_true_theta)
  } else {
    NULL
  }
  list(
    regression = list(
      global = join_cfg(metrics$regression$global),
      cell_type = join_cfg(metrics$regression$cell_type)
    ),
    monte_carlo = join_cfg(metrics$monte_carlo),
    optimisation = join_cfg(metrics$optimisation),
    config = cfg,
    theta_true = theta_list,
    descriptors = out$descriptors,
    supplementary = if (
      !is.null(out$descriptors) && "jeffreys" %in% names(out$descriptors)
    ) {
      out$descriptors[, c("ID", "jeffreys"), drop = FALSE]
    } else {
      NULL
    },
    call = metrics$call
  )
}
