#' Default files for the three Monte Carlo outcomes
#'
#' @keywords internal
#' @noRd
.convergence_icon_files <- function(icon_dir = NULL) {
  dirs <- c(
    icon_dir,
    "temp_logos",
    file.path("inst", "extdata", "convergence_icons"),
    system.file(
      "extdata",
      "convergence_icons",
      package = "DeCovarT"
    )
  )
  dirs <- dirs[!is.na(dirs) & nzchar(as.character(dirs))]
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs) == 0L) {
    return(character())
  }
  names_map <- c(
    theoretical_success = "succes_converence.png",
    numerical_failure = "failed_nuermical_convergence.png",
    theoretical_failure = "failed_convergence.png"
  )
  out <- character()
  for (nm in names(names_map)) {
    for (root in dirs) {
      cand <- file.path(root, names_map[[nm]])
      if (file.exists(cand)) {
        out[[nm]] <- cand
        break
      }
    }
  }
  out
}

#' Markdown legend labels with outcome PNG icons
#'
#' @keywords internal
#' @noRd
.convergence_markdown_labels <- function(icon_dir = NULL) {
  files <- .convergence_icon_files(icon_dir)
  labs <- .convergence_outcome_labels()
  vapply(
    .convergence_outcome_levels(),
    function(nm) {
      lab <- unname(labs[[nm]])
      f <- unname(files[nm])
      if (length(f) != 1L || is.na(f) || !file.exists(f)) {
        return(lab)
      }
      f <- normalizePath(f, winslash = "/", mustWork = TRUE)
      sprintf(
        "<img src='%s' width='22' height='22'/> %s",
        f,
        lab
      )
    },
    character(1)
  )
}

#' Outcome labels and fills
#'
#' @keywords internal
#' @noRd
.convergence_outcome_levels <- function() {
  c(
    "theoretical_success",
    "numerical_failure",
    "theoretical_failure"
  )
}

#' @keywords internal
#' @noRd
.convergence_outcome_labels <- function() {
  c(
    theoretical_success = "Theoretical success",
    numerical_failure = "Numerical failure",
    theoretical_failure = "Theoretical failure"
  )
}

#' @keywords internal
#' @noRd
.convergence_outcome_fills <- function() {
  c(
    theoretical_success = "#2E7D32",
    numerical_failure = "#C62828",
    theoretical_failure = "#EDB120"
  )
}

#' Optimisation rows at the four page corners
#'
#' @keywords internal
#' @noRd
.corner_optimisation <- function(artefacts, row) {
  ids <- .corner_ids(artefacts$config, row)
  parts <- lapply(seq_along(ids), function(k) {
    if (is.na(ids[[k]])) {
      return(NULL)
    }
    sub <- .subset_benchmark_ids(artefacts, ids[[k]])
    opt <- sub$optimisation
    if (is.null(opt) || nrow(opt) == 0L) {
      return(NULL)
    }
    opt$panel <- names(ids)[[k]]
    opt
  })
  dplyr::bind_rows(parts)
}

#' Cowplot legend strip with outcome icons
#'
#' @keywords internal
#' @noRd
.outcome_icon_legend <- function(icon_dir = NULL) {
  files <- .convergence_icon_files(icon_dir)
  labs <- .convergence_outcome_labels()
  fills <- .convergence_outcome_fills()
  keys <- .convergence_outcome_levels()
  items <- lapply(keys, function(nm) {
    lab <- unname(labs[[nm]])
    fill <- unname(fills[[nm]])
    img_file <- unname(files[nm])
    g <- cowplot::ggdraw()
    if (
      length(img_file) == 1L &&
        !is.na(img_file) &&
        file.exists(img_file) &&
        requireNamespace("png", quietly = TRUE)
    ) {
      raster <- png::readPNG(img_file)
      g <- g +
        cowplot::draw_grob(
          grid::rasterGrob(raster, interpolate = TRUE),
          x = 0.04,
          y = 0.1,
          width = 0.26,
          height = 0.8
        )
    } else if (
      length(img_file) == 1L &&
        !is.na(img_file) &&
        file.exists(img_file)
    ) {
      g <- g +
        cowplot::draw_image(
          img_file,
          x = 0.04,
          y = 0.1,
          width = 0.26,
          height = 0.8
        )
    }
    g +
      cowplot::draw_label(
        lab,
        x = 0.34,
        y = 0.5,
        hjust = 0,
        size = 11,
        color = fill,
        fontface = "bold"
      )
  })
  cowplot::plot_grid(plotlist = items, nrow = 1L)
}

#' Raincloud of projected-score KKT residuals
#'
#' @param df Long optimisation table with `kkt_residual`, `algorithm`,
#'   and `panel`.
#' @param title Page title.
#'
#' @return A `ggplot`.
#' @export
plot_kkt_raincloud <- function(df, title = NULL) {
  .check_plot_dependencies(need_ggdist = TRUE)
  needed <- c("kkt_residual", "algorithm")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop("`df` must contain ", toString(missing), ".", call. = FALSE)
  }
  df <- df[is.finite(df$kkt_residual), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No finite `kkt_residual` values to plot.", call. = FALSE)
  }
  df$algorithm <- .relevel_algorithm(df$algorithm)
  df$kkt_plot <- pmax(df$kkt_residual, 1e-12)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["kkt_plot"]],
      y = .data[["algorithm"]],
      fill = .data[["algorithm"]],
      colour = .data[["algorithm"]]
    )
  ) +
    ggdist::stat_halfeye(
      ggplot2::aes(slab_colour = ggplot2::after_scale(.data[["fill"]])),
      orientation = "horizontal",
      .width = c(0.5, 0.95),
      justification = -0.08,
      point_interval = ggdist::median_qi,
      normalize = "groups",
      scale = 0.9,
      slab_linewidth = 0.6,
      slab_alpha = 0.45,
      interval_size = 2.4,
      point_size = 1.4
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::annotation_logticks(sides = "b") +
    ggplot2::labs(
      x = "KKT residual (projected score; log10)",
      y = NULL,
      fill = "Algorithm",
      colour = "Algorithm",
      title = title,
      caption = paste(
        "Projected-score residual after the estimate lies on the simplex;",
        "not a |sum(p)-1| violation before L-BFGS-B renormalisation.",
        "Inner 50% and 95% bands are Monte Carlo quantiles."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  if ("panel" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
  }
  p
}

#' Stacked bar counts of theoretical success and both failure modes
#'
#' @param df Optimisation rows with `numerical_converged`,
#'   `theoretical_converged`, `algorithm`, and `panel`.
#' @param title Page title.
#' @param icon_dir Directory containing the three outcome PNGs
#'   (`succes_converence.png`, `failed_nuermical_convergence.png`,
#'   `failed_convergence.png`). When `ggtext` is installed, those
#'   icons are rendered in the fill legend.
#'
#' @return A `ggplot` or cowplot grob.
#' @export
plot_convergence_stacked <- function(
  df,
  title = NULL,
  icon_dir = NULL
) {
  use_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  if (isTRUE(use_ggtext)) {
    .check_plot_dependencies(need_ggtext = TRUE)
  } else {
    .check_plot_dependencies(need_cowplot = TRUE)
  }
  needed <- c("numerical_converged", "theoretical_converged", "algorithm")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop("`df` must contain ", toString(missing), ".", call. = FALSE)
  }
  df$outcome <- .simulation_outcome(
    df$numerical_converged,
    df$theoretical_converged
  )
  df$outcome <- factor(
    df$outcome,
    levels = .convergence_outcome_levels()
  )
  df$algorithm <- .relevel_algorithm(df$algorithm)
  keys <- c("algorithm", "outcome")
  if ("panel" %in% names(df)) {
    keys <- c(keys, "panel")
  }
  agg <- df |>
    dplyr::count(dplyr::across(dplyr::all_of(keys)), name = "n")
  fills <- .convergence_outcome_fills()
  labs <- .convergence_outcome_labels()
  fill_labs <- if (isTRUE(use_ggtext)) {
    .convergence_markdown_labels(icon_dir)
  } else {
    labs
  }
  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = .data[["algorithm"]],
      y = .data[["n"]],
      fill = .data[["outcome"]]
    )
  ) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::scale_fill_manual(
      values = fills,
      labels = fill_labs,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Monte Carlo replicates",
      fill = NULL,
      title = title,
      caption = paste(
        "Disjoint partition: theoretical success;",
        "numerical failure (nested in theoretical failure);",
        "theoretical failure with a finite simplex.",
        "Icons: green success, red numerical failure,",
        "yellow theoretical failure."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  if (isTRUE(use_ggtext)) {
    p <- p +
      ggplot2::theme(
        legend.text = ggtext::element_markdown(size = 10, lineheight = 1.1),
        legend.key.height = grid::unit(16, "pt"),
        legend.spacing.x = grid::unit(8, "pt")
      )
  }
  if ("panel" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
  }
  if (isTRUE(use_ggtext)) {
    return(p)
  }
  cowplot::plot_grid(
    p + ggplot2::theme(legend.position = "none"),
    .outcome_icon_legend(icon_dir),
    ncol = 1L,
    rel_heights = c(1, 0.2)
  )
}

#' Expected Fisher information in ILR coordinates at p*
#'
#' @keywords internal
#' @noRd
.ilr_expected_information <- function(p, mu, Sigma) {
  info_p <- expected_fisher_unconstrained(p, mu, Sigma)
  z <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(z)
  t(jac) %*% info_p %*% jac
}

#' Whitened ILR coordinates of one Monte Carlo estimate
#'
#' \eqn{W = I_z(p*)^{1/2}(\hat z - z*)}. Boundary estimates (ILR
#' undefined) return `NULL`.
#'
#' @keywords internal
#' @noRd
.whiten_ilr_estimate <- function(p_hat, p_true, info_z) {
  floor <- 100 * .Machine$double.eps
  if (any(p_hat <= floor | p_hat >= 1 - floor)) {
    return(NULL)
  }
  z_hat <- isometric_log_ratio(p_hat)
  z_true <- isometric_log_ratio(p_true)
  dz <- as.numeric(z_hat - z_true)
  r_chol <- tryCatch(chol(info_z), error = function(e) NULL)
  if (is.null(r_chol)) {
    return(NULL)
  }
  w <- as.numeric(r_chol %*% dz)
  list(whitened = w, mahalanobis = sum(w * w))
}

#' Whitened ILR Monte Carlo table for iterative DeCovarT solvers
#'
#' Only L-BFGS, Marquardt--Levenberg and Newton--Raphson. Gradient and
#' simulated annealing are omitted.
#'
#' @keywords internal
#' @noRd
.whitened_corner_estimates <- function(artefacts, row) {
  artefacts <- .fill_theta_true_from_theta(artefacts)
  ids <- .corner_ids(artefacts$config, row)
  keep <- .convolution_algorithms()
  parts <- lapply(seq_along(ids), function(k) {
    if (is.na(ids[[k]])) {
      return(NULL)
    }
    sub <- .subset_benchmark_ids(artefacts, ids[[k]])
    opt <- sub$optimisation
    th <- sub$theta_true
    if (is.null(opt) || nrow(opt) == 0L || length(th) == 0L) {
      return(NULL)
    }
    th <- .unwrap_true_theta(th[[1L]])
    nms <- .truth_cell_names(th)
    if (!all(nms %in% names(opt))) {
      return(NULL)
    }
    p_true <- stats::setNames(as.numeric(th$p), nms)
    info_z <- tryCatch(
      .ilr_expected_information(p_true, th$mu, th$sigma),
      error = function(e) NULL
    )
    if (is.null(info_z)) {
      return(NULL)
    }
    opt <- opt[.algorithm_in(opt$algorithm, keep), , drop = FALSE]
    if (nrow(opt) == 0L) {
      return(NULL)
    }
    rows <- lapply(seq_len(nrow(opt)), function(i) {
      p_hat <- stats::setNames(as.numeric(opt[i, nms, drop = TRUE]), nms)
      w <- .whiten_ilr_estimate(p_hat, p_true, info_z)
      if (is.null(w)) {
        return(NULL)
      }
      n_ilr <- length(w$whitened)
      tibble::tibble(
        algorithm = opt$algorithm[[i]],
        sample_id = if ("sample_id" %in% names(opt)) {
          opt$sample_id[[i]]
        } else {
          paste0("sample_", i)
        },
        panel = names(ids)[[k]],
        ilr_coord = paste0("ILR", seq_len(n_ilr)),
        whitened = w$whitened,
        mahalanobis = w$mahalanobis,
        chi_df = n_ilr
      )
    })
    dplyr::bind_rows(rows)
  })
  df <- dplyr::bind_rows(parts)
  if (nrow(df) == 0L) {
    return(df)
  }
  df$algorithm <- .relevel_algorithm(df$algorithm)
  df[is.finite(df$whitened), , drop = FALSE]
}

#' One Mahalanobis row per replicate (joint \eqn{D^2})
#'
#' @keywords internal
#' @noRd
.mahalanobis_from_whitened <- function(df) {
  keys <- intersect(
    c("algorithm", "panel", "sample_id", "chi_df", "page"),
    names(df)
  )
  dplyr::distinct(
    df,
    dplyr::across(dplyr::all_of(c(keys, "mahalanobis")))
  )
}

#' Overlay a cowplot legend in one facet of a 2-by-2 page
#'
#' @keywords internal
#' @noRd
.inset_cowplot_legend <- function(
  p,
  x = 0.56,
  y = 0.06,
  width = 0.3,
  height = 0.36
) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    return(p + ggplot2::theme(legend.position = "right"))
  }
  p_leg <- p + ggplot2::theme(legend.position = "right")
  legend <- tryCatch(
    cowplot::get_legend(p_leg),
    error = function(e) {
      tryCatch(
        cowplot::get_plot_component(
          p_leg,
          "guide-box",
          return_all = TRUE
        )[[1L]],
        error = function(e2) NULL
      )
    }
  )
  if (is.null(legend)) {
    return(p_leg)
  }
  p_main <- p + ggplot2::theme(legend.position = "none")
  cowplot::ggdraw(p_main) +
    cowplot::draw_plot(legend, x, y, width, height)
}

#' Q-Q plot of whitened ILR coordinates versus N(0, 1)
#'
#' Tail-sensitive simultaneous envelope (`qqplotr` `bandType = "ts"`)
#' calibrated to the theoretical \eqn{N(0,1)} law
#' (`identity = TRUE`, `mu = 0`, `sigma = 1`). Colour and fill distinguish
#' solvers. The plot is **not** detrended.
#'
#' @param df Long table with `whitened`, `algorithm`, and `panel`.
#' @param title Page title.
#'
#' @return A `ggplot`.
#' @export
plot_mc_qq_normal <- function(df, title = NULL) {
  needed <- c("whitened", "algorithm")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop("`df` must contain ", toString(missing), ".", call. = FALSE)
  }
  df <- df[.algorithm_in(df$algorithm, .convolution_algorithms()), ]
  df$algorithm <- .relevel_algorithm(df$algorithm)
  if (requireNamespace("qqplotr", quietly = TRUE)) {
    p <- .plot_mc_qq_qqplotr(df, title)
  } else {
    p <- .plot_mc_qq_ggplot(df, title)
  }
  n_ilr <- if ("ilr_coord" %in% names(df)) {
    dplyr::n_distinct(df$ilr_coord)
  } else {
    1L
  }
  if ("panel" %in% names(df) && n_ilr > 1L) {
    p <- p + ggplot2::facet_grid(ilr_coord ~ panel)
  } else if ("panel" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
  } else if (n_ilr > 1L) {
    p <- p + ggplot2::facet_wrap(~ilr_coord, ncol = 2L)
  }
  p
}

#' @keywords internal
#' @noRd
.plot_mc_qq_qqplotr <- function(df, title) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      sample = .data[["whitened"]],
      colour = .data[["algorithm"]],
      fill = .data[["algorithm"]]
    )
  ) +
    qqplotr::stat_qq_band(
      distribution = "norm",
      dparams = list(mean = 0, sd = 1),
      bandType = "ts",
      B = 1000L,
      conf = 0.95,
      identity = TRUE,
      mu = 0,
      sigma = 1,
      detrend = FALSE,
      alpha = 0.28
    ) +
    qqplotr::stat_qq_line(
      distribution = "norm",
      dparams = list(mean = 0, sd = 1),
      identity = TRUE,
      detrend = FALSE
    ) +
    qqplotr::stat_qq_point(
      distribution = "norm",
      dparams = list(mean = 0, sd = 1),
      detrend = FALSE,
      size = 1.05,
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = "Theoretical N(0, 1) quantiles",
      y = "Whitened ILR Monte Carlo quantiles",
      colour = "Algorithm",
      fill = "Algorithm",
      title = title,
      caption = paste(
        "W = I_z(p*)^{1/2} (hat z - z*) versus N(0, 1); identity line.",
        "Tail-sensitive band (Aldor-Noiman et al.) with mu = 0, sigma = 1.",
        "LBFGS, Marquardt-Levenberg and Newton-Raphson only; not detrended."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      legend.position = "bottom"
    )
}

#' Pointwise beta Q-Q envelope without qqplotr
#'
#' @keywords internal
#' @noRd
.hazen_qq_band <- function(n, alpha = 0.05, type = c("pointwise", "ks")) {
  type <- match.arg(type)
  i <- seq_len(as.integer(n))
  p <- (i - 0.5) / n
  x <- stats::qnorm(p)
  if (identical(type, "pointwise")) {
    lo <- stats::qnorm(stats::qbeta(alpha / 2, i, n - i + 1L))
    hi <- stats::qnorm(stats::qbeta(1 - alpha / 2, i, n - i + 1L))
  } else {
    d <- sqrt(-0.5 * log(alpha / 2) / n)
    lo <- stats::qnorm(pmax(p - d, .Machine$double.eps))
    hi <- stats::qnorm(pmin(p + d, 1 - .Machine$double.eps))
  }
  data.frame(x = x, lo = lo, hi = hi, band = type)
}

#' @keywords internal
#' @noRd
.plot_mc_qq_ggplot <- function(df, title) {
  keys <- "algorithm"
  if ("ilr_coord" %in% names(df)) {
    keys <- c(keys, "ilr_coord")
  }
  if ("panel" %in% names(df)) {
    keys <- c(keys, "panel")
  }
  keys <- unique(keys[keys %in% names(df)])
  n <- max(dplyr::count(df, dplyr::across(dplyr::all_of(keys)))$n)
  if (!is.finite(n) || n < 8L) {
    n <- nrow(df)
  }
  bands <- .hazen_qq_band(n, type = "pointwise")
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      sample = .data[["whitened"]],
      colour = .data[["algorithm"]]
    )
  ) +
    ggplot2::geom_ribbon(
      data = bands,
      ggplot2::aes(
        x = .data[["x"]],
        ymin = .data[["lo"]],
        ymax = .data[["hi"]]
      ),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.35
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.4) +
    ggplot2::stat_qq(size = 1.05, alpha = 0.7) +
    ggplot2::labs(
      x = "Theoretical N(0, 1) quantiles",
      y = "Whitened ILR Monte Carlo quantiles",
      colour = "Algorithm",
      title = title,
      caption = paste(
        "W = I_z(p*)^{1/2} (hat z - z*) versus N(0, 1).",
        "qqplotr is unavailable; fallback is a Hazen pointwise beta envelope."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      legend.position = "bottom"
    )
}

#' Monte Carlo chi-square Q-Q envelope under the exact null
#'
#' @keywords internal
#' @noRd
.chisq_qq_envelope <- function(n, df, B = 2000L, conf = 0.95) {
  n <- as.integer(n)
  df <- as.numeric(df)
  i <- seq_len(n)
  p <- (i - 0.5) / n
  x <- stats::qchisq(p, df = df)
  sims <- matrix(stats::rchisq(n * B, df = df), nrow = n, ncol = B)
  sims <- apply(sims, 2L, sort)
  alpha <- (1 - conf) / 2
  data.frame(
    x = x,
    lo = apply(sims, 1L, stats::quantile, probs = alpha, names = FALSE),
    hi = apply(sims, 1L, stats::quantile, probs = 1 - alpha, names = FALSE)
  )
}

#' Q-Q plot of joint Mahalanobis D^2 versus chi-square_{J-1}
#'
#' @param df Table with `mahalanobis`, `algorithm`, `chi_df`, and `panel`.
#' @param title Page title.
#'
#' @return A `ggplot`.
#' @export
plot_mc_chi2_qq <- function(df, title = NULL) {
  needed <- c("mahalanobis", "algorithm")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop("`df` must contain ", toString(missing), ".", call. = FALSE)
  }
  df <- df[.algorithm_in(df$algorithm, .convolution_algorithms()), ]
  df$algorithm <- .relevel_algorithm(df$algorithm)
  df_chi <- if ("chi_df" %in% names(df) && any(is.finite(df$chi_df))) {
    as.integer(stats::median(df$chi_df, na.rm = TRUE))
  } else {
    1L
  }
  keys <- "algorithm"
  if ("panel" %in% names(df)) {
    keys <- c(keys, "panel")
  }
  keys <- keys[keys %in% names(df)]
  n <- max(dplyr::count(df, dplyr::across(dplyr::all_of(keys)))$n)
  if (!is.finite(n) || n < 8L) {
    n <- nrow(df)
  }
  bands <- .chisq_qq_envelope(n, df_chi)
  qfun <- function(p) stats::qchisq(p, df = df_chi)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      sample = .data[["mahalanobis"]],
      colour = .data[["algorithm"]]
    )
  ) +
    ggplot2::geom_ribbon(
      data = bands,
      ggplot2::aes(
        x = .data[["x"]],
        ymin = .data[["lo"]],
        ymax = .data[["hi"]]
      ),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.4
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.4) +
    ggplot2::stat_qq(
      distribution = qfun,
      size = 1.05,
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = bquote(Theoretical ~ chi^2[.(df_chi)] ~ quantiles),
      y = "Mahalanobis D^2 of whitened ILR",
      colour = "Algorithm",
      title = title,
      caption = paste(
        "D^2 = W^T W versus chi^2_{J-1} (joint covariance of the ILR MLE).",
        "Grey band: 95% Monte Carlo envelope under the exact chi^2 null.",
        "LBFGS, Marquardt-Levenberg and Newton-Raphson; identity reference."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      legend.position = "bottom"
    )
  if ("panel" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
  }
  p
}

#' Kolmogorov-Smirnov p-values versus N(0, 1) by algorithm
#'
#' One lollipop per scenario (panel) and solver, pooling whitened ILR
#' coordinates. Values below 0.05 reject normality at the conventional
#' 5% level.
#'
#' @param df Long whitened table.
#' @param title Page title.
#'
#' @return A `ggplot` or cowplot grob with an inset legend.
#' @export
plot_ks_normality_box <- function(df, title = NULL) {
  needed <- c("whitened", "algorithm")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop("`df` must contain ", toString(missing), ".", call. = FALSE)
  }
  df <- df[.algorithm_in(df$algorithm, .convolution_algorithms()), ]
  grp <- "algorithm"
  if ("panel" %in% names(df)) {
    grp <- c(grp, "panel")
  }
  ks_fun <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 8L || dplyr::n_distinct(round(x, 10)) < 3L) {
      return(NA_real_)
    }
    suppressWarnings(stats::ks.test(x, "pnorm")$p.value)
  }
  ks_tbl <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::summarise(
      ks_p = ks_fun(.data[["whitened"]]),
      .groups = "drop"
    )
  ks_tbl$algorithm <- .relevel_algorithm(ks_tbl$algorithm)
  ks_tbl$reject <- is.finite(ks_tbl$ks_p) & ks_tbl$ks_p < 0.05
  ks_tbl$seg_w <- ifelse(ks_tbl$reject, 1.25, 0.7)
  ks_tbl$pt_s <- ifelse(ks_tbl$reject, 3.6, 2.6)
  p <- ggplot2::ggplot(
    ks_tbl,
    ggplot2::aes(
      x = .data[["ks_p"]],
      y = .data[["algorithm"]],
      colour = .data[["algorithm"]]
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data[["ks_p"]],
        y = .data[["algorithm"]],
        yend = .data[["algorithm"]],
        linewidth = .data[["seg_w"]]
      )
    ) +
    ggplot2::geom_point(ggplot2::aes(size = .data[["pt_s"]])) +
    ggplot2::geom_vline(
      xintercept = 0.05,
      linetype = "dashed",
      colour = "#C62828"
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(
      x = "KS p-value versus N(0, 1)",
      y = NULL,
      colour = "Algorithm",
      title = title,
      caption = paste(
        "One KS score per scenario x solver (whitened ILR coordinates pooled).",
        "Dashed line: alpha = 0.05. Below 0.05 we reject normality of the",
        "Monte Carlo distribution at the conventional 5% level."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  if ("panel" %in% names(ks_tbl)) {
    p <- p + ggplot2::facet_wrap(~panel, ncol = 2L)
  }
  .inset_cowplot_legend(p)
}


#' Multi-page KKT raincloud book
#'
#' @inheritParams save_bivariate_raincloud_book
#' @export
save_bivariate_kkt_raincloud_book <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  .check_plot_dependencies(need_ggdist = TRUE, need_cowplot = TRUE)
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  grDevices::pdf(file, width = 16, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    df <- .corner_optimisation(artefacts, row)
    if (nrow(df) == 0L || !any(is.finite(df$kkt_residual))) {
      next
    }
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    p <- plot_kkt_raincloud(df, title = .page_title(row))
    print(.attach_cowplot_legend(p))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(
      dplyr::bind_rows(collected),
      data_rds,
      "kkt_raincloud"
    )
  }
  invisible(file)
}

#' Multi-page stacked-bar book of solver outcomes
#'
#' @inheritParams save_bivariate_raincloud_book
#' @param icon_dir Directory of outcome PNG icons.
#' @export
save_bivariate_convergence_book <- function(
  artefacts,
  file,
  data_rds = NULL,
  icon_dir = NULL
) {
  .check_plot_dependencies(need_cowplot = TRUE)
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  grDevices::pdf(file, width = 16, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    df <- .corner_optimisation(artefacts, row)
    if (nrow(df) == 0L) {
      next
    }
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    print(
      plot_convergence_stacked(
        df,
        title = .page_title(row),
        icon_dir = icon_dir
      )
    )
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(
      dplyr::bind_rows(collected),
      data_rds,
      "convergence_stacked"
    )
  }
  invisible(file)
}

#' Multi-page Q-Q book of whitened ILR coordinates versus N(0, 1)
#'
#' @inheritParams save_bivariate_raincloud_book
#' @export
save_bivariate_qq_book <- function(artefacts, file, data_rds = NULL) {
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  .open_ggplot_pdf(file, width = 16, height = 14)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    df <- .whitened_corner_estimates(artefacts, row)
    if (nrow(df) == 0L) {
      next
    }
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    print(plot_mc_qq_normal(df, title = .page_title(row)))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "qq_normal")
  }
  invisible(file)
}

#' Multi-page KS lollipop book of whitened ILR normality
#'
#' @inheritParams save_bivariate_raincloud_book
#' @export
save_bivariate_ks_book <- function(artefacts, file, data_rds = NULL) {
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  .open_ggplot_pdf(file, width = 14, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    df <- .whitened_corner_estimates(artefacts, row)
    if (nrow(df) == 0L) {
      next
    }
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    print(plot_ks_normality_box(df, title = .page_title(row)))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "ks_normal")
  }
  invisible(file)
}

#' Multi-page Q-Q book of joint Mahalanobis D^2 versus chi-square
#'
#' @inheritParams save_bivariate_raincloud_book
#' @export
save_bivariate_chi2_book <- function(artefacts, file, data_rds = NULL) {
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  .open_ggplot_pdf(file, width = 16, height = 14)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    df <- .whitened_corner_estimates(artefacts, row)
    if (nrow(df) == 0L) {
      next
    }
    df <- .mahalanobis_from_whitened(df)
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    print(plot_mc_chi2_qq(df, title = .page_title(row)))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "qq_chi2")
  }
  invisible(file)
}
