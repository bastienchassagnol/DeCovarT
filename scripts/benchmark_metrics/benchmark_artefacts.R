# Script-side helpers of the simulation pipelines (fig02, fig03, supp S1-S4)
# and of the figure books under scripts/visualisation/. Not part of the
# package namespace: they persist and reload benchmark artefacts, apply the
# fig02 / fig03 factor orders and display labels, and refresh Wald coverage
# from the expected Fisher information at the true composition.
#
# Usage: load the package first, then source this file, e.g.
#   devtools::load_all(".")
#   source(file.path("scripts", "auxiliary", "benchmark_artefacts.R"))
# The helpers call package internals (.ilr_wald_se(), .relevel_algorithm(),
# describe_simulation_scenario(), pivot_mc_estimates()).

# Unwrap a list-column `true_theta` cell
#
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

# Columns that belong on the descriptors artefact (not the design grid)
#
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

# Rebuild descriptor metrics from `\theta`, keeping MixSim if present
#
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

# Write split simulation artefacts (config, descriptors, theta, metrics)
#
# [run_simulation_benchmark()] still returns a single list for tests.
# Scripts persist four RDS files keyed by `ID` so the metrics object
# does not duplicate geometry, MixSim, or `\theta`.
#
# @param benchmark List from [run_simulation_benchmark()].
# @param dir Output directory.
# @param stem File-name stem (`bivariate`, `hybrid`, …).
# @param config Optional scenario grid (defaults to `benchmark$config`).
#   Must contain an `ID` column.
#
# @return Invisibly, a named list of written paths.
write_simulation_artefacts <- function(
  benchmark,
  dir,
  stem,
  config = NULL
) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(config)) {
    config <- benchmark$config
  }
  config <- tibble::as_tibble(config)
  if (!"ID" %in% names(config)) {
    stop("`config` must contain an `ID` column.", call. = FALSE)
  }
  config$ID <- as.character(config$ID)

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
    old_desc$ID <- as.character(old_desc$ID)
  }
  descriptors <- .descriptors_from_theta(theta_tbl, old_desc)
  keep_desc <- intersect(
    c("ID", .descriptor_metric_columns()),
    names(descriptors)
  )
  descriptors <- descriptors[, keep_desc, drop = FALSE]

  config_slim <- tibble::as_tibble(config)
  config_slim$true_theta <- NULL
  config_slim$scenario_idx <- NULL

  slim_metrics <- function(tbl) {
    if (is.null(tbl) || ncol(tbl) == 0L) {
      return(tbl)
    }
    tbl <- tibble::as_tibble(tbl)
    design_drop <- setdiff(
      names(config_slim),
      c("ID")
    )
    # Design columns repeated from the config, plus redundant grid keys.
    redundant <- c(
      "scenario_idx",
      "rho_ct1",
      "rho_ct2",
      "proportion_name",
      "centroid"
    )
    drop <- unique(c(redundant, intersect(design_drop, names(tbl))))
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

# Read split simulation artefacts and optionally reassemble a benchmark list
#
# @param assemble If `TRUE`, join config onto metric tables for plotting
#   helpers that still expect design columns.
#
# @return Named list of tibbles, or a benchmark-like list when
#   `assemble = TRUE`.
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
    out$config <- .relevel_scenario_table(out$config)
    out$descriptors <- .relevel_scenario_table(out$descriptors)
    return(out)
  }
  cfg <- .relevel_scenario_table(out$config)
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
      global = .relevel_scenario_table(join_cfg(metrics$regression$global)),
      cell_type = .relevel_scenario_table(
        join_cfg(metrics$regression$cell_type)
      )
    ),
    monte_carlo = .relevel_scenario_table(join_cfg(metrics$monte_carlo)),
    optimisation = .relevel_scenario_table(join_cfg(metrics$optimisation)),
    config = cfg,
    theta_true = theta_list,
    descriptors = .relevel_scenario_table(out$descriptors),
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

# Cell-type names from a `theta_true` element
#
.truth_cell_names <- function(th) {
  p <- th$p
  nms <- names(p)
  if (!is.null(nms)) {
    return(nms)
  }
  cn <- colnames(th$mu)
  if (!is.null(cn)) {
    return(cn)
  }
  paste0("ct", seq_along(p))
}

# Fill missing `theta_true` from a theta artefact table
#
.fill_theta_true_from_theta <- function(benchmark) {
  cfg <- benchmark$config
  truths <- benchmark$theta_true
  if (!is.null(truths) && length(truths) == nrow(cfg)) {
    return(benchmark)
  }
  th_tbl <- benchmark$theta
  if (is.null(th_tbl) || !"true_theta" %in% names(th_tbl)) {
    return(benchmark)
  }
  if ("ID" %in% names(cfg) && "ID" %in% names(th_tbl)) {
    idx <- match(as.character(cfg$ID), as.character(th_tbl$ID))
    benchmark$theta_true <- lapply(
      th_tbl$true_theta[idx],
      .unwrap_true_theta
    )
    return(benchmark)
  }
  if (nrow(th_tbl) == nrow(cfg)) {
    benchmark$theta_true <- lapply(th_tbl$true_theta, .unwrap_true_theta)
  }
  benchmark
}

# Refresh Wald coverage from expected Fisher at \eqn{p^{\star}}
#
# [expected_fisher_unconstrained()] / [vcov_ilr_delta()] depend only on
# \eqn{(p,\mu,\Sigma)}, so the Cramer--Rao SE is constant for a
# scenario (the same object as [confint.decovart_fit()]). Coverage is
# the Monte Carlo rate at which
# \eqn{\hat p_j\pm z\,\mathrm{SE}_j(p^{\star})} covers
# \eqn{p_j^{\star}}. Mean-only solvers are left unchanged.
#
.attach_expected_fisher_wald <- function(benchmark, level = 0.95) {
  z <- stats::qnorm((1 + level) / 2)
  cfg <- benchmark$config
  mc <- benchmark$monte_carlo
  if (is.null(mc) || is.null(cfg) || nrow(mc) == 0L) {
    return(benchmark)
  }
  benchmark <- .fill_theta_true_from_theta(benchmark)
  truths <- benchmark$theta_true
  if (is.null(truths) || length(truths) == 0L) {
    return(benchmark)
  }
  se_pieces <- lapply(seq_along(truths), function(i) {
    th <- .unwrap_true_theta(truths[[i]])
    nms <- .truth_cell_names(th)
    se <- .ilr_wald_se(
      stats::setNames(as.numeric(th$p), nms),
      th$mu,
      th$sigma,
      warn = FALSE
    )
    row <- tibble::tibble(cell_type = nms, theoretical_se = as.numeric(se[nms]))
    if ("ID" %in% names(cfg)) {
      row$ID <- as.character(cfg$ID[[i]])
    }
    row
  })
  se_tbl <- dplyr::bind_rows(se_pieces)
  long <- pivot_mc_estimates(benchmark)
  join_se <- intersect(c("ID", "cell_type"), names(long))
  join_se <- join_se[join_se %in% names(se_tbl)]
  if ("ID" %in% names(long)) {
    long$ID <- as.character(long$ID)
  }
  if ("ID" %in% names(se_tbl)) {
    se_tbl$ID <- as.character(se_tbl$ID)
  }
  long <- dplyr::left_join(long, se_tbl, by = join_se)
  long$is_wald <- .uses_ilr_wald_algorithm(long$algorithm)
  long$covered <- long$is_wald &
    is.finite(long$estimate) &
    is.finite(long$p_true) &
    is.finite(long$theoretical_se) &
    abs(long$estimate - long$p_true) <= z * long$theoretical_se
  method <- "wilson"
  if ("coverage_interval" %in% names(mc)) {
    methods <- unique(stats::na.omit(as.character(mc$coverage_interval)))
    if (length(methods) == 1L) {
      method <- methods[[1L]]
    }
  }
  grp <- intersect(c("ID", "algorithm", "cell_type"), names(long))
  wald_long <- long[long$is_wald, , drop = FALSE]
  if (nrow(wald_long) == 0L) {
    return(benchmark)
  }
  grouped <- dplyr::group_by(wald_long, dplyr::across(dplyr::all_of(grp)))
  upd <- dplyr::summarise(
    grouped,
    theoretical_se = mean(.data[["theoretical_se"]], na.rm = TRUE),
    covered_list = list(.data[["covered"]]),
    .groups = "drop"
  )
  intervals <- lapply(
    upd$covered_list,
    function(cv) coverage_mc_interval(cv, method = method)
  )
  upd$coverage <- vapply(intervals, `[[`, numeric(1), "coverage")
  upd$coverage_lower <- vapply(intervals, `[[`, numeric(1), "lower")
  upd$coverage_upper <- vapply(intervals, `[[`, numeric(1), "upper")
  upd$mcse_coverage <- vapply(intervals, `[[`, numeric(1), "mcse")
  upd$coverage_interval <- vapply(intervals, `[[`, character(1), "method")
  upd$mean_model_se <- upd$theoretical_se
  upd$mean_model_sd <- upd$theoretical_se
  upd$mean_interval_width <- 2 * z * upd$theoretical_se
  upd$covered_list <- NULL
  replace_cols <- c(
    "theoretical_se",
    "coverage",
    "coverage_lower",
    "coverage_upper",
    "coverage_interval",
    "mcse_coverage",
    "mean_model_se",
    "mean_model_sd",
    "mean_interval_width"
  )
  is_wald_mc <- .uses_ilr_wald_algorithm(mc$algorithm)
  mc_wald <- mc[is_wald_mc, , drop = FALSE]
  mc_rest <- mc[!is_wald_mc, , drop = FALSE]
  drop_now <- intersect(replace_cols, names(mc_wald))
  if (length(drop_now) > 0L) {
    mc_wald <- mc_wald[, setdiff(names(mc_wald), drop_now), drop = FALSE]
  }
  join_mc <- intersect(grp, names(mc_wald))
  if ("ID" %in% names(mc_wald)) {
    mc_wald$ID <- as.character(mc_wald$ID)
  }
  if ("ID" %in% names(upd)) {
    upd$ID <- as.character(upd$ID)
  }
  mc_wald <- dplyr::left_join(mc_wald, upd, by = join_mc)
  if ("empirical_sd" %in% names(mc_wald)) {
    mc_wald$se_sd_ratio <- mc_wald$theoretical_se / mc_wald$empirical_sd
  }
  if (!"theoretical_se" %in% names(mc_rest)) {
    mc_rest$theoretical_se <- NA_real_
  }
  benchmark$monte_carlo <- dplyr::bind_rows(mc_wald, mc_rest)
  benchmark$theta_true <- truths
  benchmark
}

.normalise_algorithm_key <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(as.character(x)))
}

.uses_ilr_wald_algorithm <- function(name) {
  .normalise_algorithm_key(name) %in%
    c(
      "newtonraphson",
      "newton",
      "marquardtlevenberg",
      "marquardt",
      "lbfgs",
      "lbfgsb",
      "gradient",
      "gradientdescent",
      "bfgs",
      "sa",
      "sann",
      "simulatedannealing"
    )
}

# Canonical factor orders for fig02 / fig03 scenario columns
#
# Codes are releveled first, then display labels are applied with
# [forcats::fct_relabel()] (graph families and MixSim overlap).
#
.scenario_level_orders <- function() {
  list(
    centroids = c("small_CLD", "large_CLD"),
    variance = c("homoscedastic", "heteroscedastic"),
    proportions = c(
      "balanced",
      "moderately unbalanced",
      "highly unbalanced"
    ),
    proportion_name = c(
      "balanced",
      "moderately unbalanced",
      "highly unbalanced"
    ),
    graph_ct1 = c("scale_free", "stochastic_block_model"),
    graph_ct2 = c("scale_free", "stochastic_block_model"),
    graph_ct3 = c("scale_free", "stochastic_block_model"),
    overlap_label = c("low", "moderate", "high")
  )
}

# Graph-family codes to facet / axis labels
#
.graph_display_labels <- function() {
  c(
    scale_free = "Scale-free",
    stochastic_block_model = "Cluster SBM"
  )
}

# MixSim overlap codes to facet / axis labels
#
.overlap_display_labels <- function() {
  c(
    low = "low OVL",
    moderate = "moderate OVL",
    high = "high OVL"
  )
}

# Map a labelled or raw graph / overlap string back to its code
#
.canonical_code <- function(x, labels) {
  x <- as.character(x)
  codes <- names(labels)
  labs <- unname(labels)
  out <- x
  hit_lab <- match(x, labs)
  out[!is.na(hit_lab)] <- codes[hit_lab[!is.na(hit_lab)]]
  hit_code <- x %in% codes
  out[hit_code] <- x[hit_code]
  out
}

# Relevel codes then [forcats::fct_relabel()] to display names
#
.relevel_relabel <- function(x, code_order, labels) {
  x_code <- .canonical_code(x, labels)
  x_fac <- .relevel_existing(x_code, code_order)
  lab_fun <- function(lvl) {
    mapped <- unname(labels[lvl])
    ifelse(is.na(mapped), lvl, mapped)
  }
  if (requireNamespace("forcats", quietly = TRUE)) {
    forcats::fct_relabel(x_fac, lab_fun)
  } else {
    new_levels <- lab_fun(levels(x_fac))
    factor(
      lab_fun(as.character(x_fac)),
      levels = unique(new_levels)
    )
  }
}

# Relevel a vector, keeping only levels that appear
#
.relevel_existing <- function(x, order) {
  x_chr <- as.character(x)
  present <- intersect(order, unique(x_chr))
  extra <- setdiff(unique(x_chr), present)
  levels <- c(present, extra)
  if (length(levels) == 0L) {
    return(factor(x_chr))
  }
  if (requireNamespace("forcats", quietly = TRUE) && length(present) > 0L) {
    forcats::fct_relevel(factor(x_chr, levels = levels), present)
  } else {
    factor(x_chr, levels = levels)
  }
}

# Relevel scenario and algorithm columns used in fig02 / fig03 plots
#
# Graph families and MixSim overlap are releveled as codes, then
# relabelled with [forcats::fct_relabel()] (`Scale-free`,
# `Cluster SBM`, `low OVL`, …).
#
.relevel_scenario_table <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl)) {
    return(tbl)
  }
  orders <- .scenario_level_orders()
  graph_nms <- c("graph_ct1", "graph_ct2", "graph_ct3")
  graph_labs <- .graph_display_labels()
  ovl_labs <- .overlap_display_labels()
  for (nm in names(orders)) {
    if (!nm %in% names(tbl)) {
      next
    }
    if (nm %in% graph_nms) {
      tbl[[nm]] <- .relevel_relabel(
        tbl[[nm]],
        orders[[nm]],
        graph_labs
      )
    } else if (identical(nm, "overlap_label")) {
      tbl[[nm]] <- .relevel_relabel(
        tbl[[nm]],
        orders[[nm]],
        ovl_labs
      )
    } else {
      tbl[[nm]] <- .relevel_existing(tbl[[nm]], orders[[nm]])
    }
  }
  if ("algorithm" %in% names(tbl)) {
    tbl$algorithm <- .relevel_algorithm(tbl$algorithm)
  }
  tbl
}

# Open a PDF device with cairo when it is available
#
.open_ggplot_pdf <- function(file, width, height) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      file,
      width = width,
      height = height,
      antialias = "subpixel",
      fallback_resolution = 320
    )
  } else {
    grDevices::pdf(
      file,
      width = width,
      height = height,
      useDingbats = FALSE
    )
  }
  invisible(TRUE)
}
