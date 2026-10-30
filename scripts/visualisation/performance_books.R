# Performance books: accuracy, interval coverage and width, the Monte Carlo
# law of the estimates, the Wald approximation, agreement between methods,
# convergence and cost. As in density_books.R, `.data_*` builds a table and
# `plot_*` draws its pages.

# Accuracy -----------------------------------------------------------------------

# Mean RMSE and Aitchison distance of every method in every scenario,
# against the target of each simulated bulk.
.data_errors <- function(benchmark, layout) {
  cts <- .celltype_names(benchmark)
  opt <- benchmark$optimisation
  sample_index <- as.integer(sub("sample_", "", opt$sample_id))
  estimate <- as.matrix(opt[cts])
  truth <- t(mapply(
    function(id, s) benchmark$data[[id]]$truth[cts, s],
    opt$ID,
    sample_index
  ))
  opt$rmse <- sqrt(rowMeans((estimate - truth)^2))
  opt$aitchison <- vapply(
    seq_len(nrow(opt)),
    function(i) .aitchison_distance(truth[i, ], estimate[i, ]),
    numeric(1)
  )
  opt |>
    dplyr::group_by(.data$ID, .data$algorithm) |>
    dplyr::summarise(
      rmse = mean(.data$rmse, na.rm = TRUE),
      aitchison = mean(.data$aitchison, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(algorithm = .method_factor(.data$algorithm)) |>
    .with_layout(layout)
}

.metric_names <- c(rmse = "Mean RMSE", aitchison = "Mean Aitchison distance")

# Grid cells of a page (panel_x by panel_y, for every method) that hold no
# simulated scenario: the fig05 a = 1 row outside the baseline page.
.missing_cells <- function(df) {
  if (!isTRUE(df$grid_panels[[1L]])) {
    return(df[0L, c("algorithm", "panel_x", "panel_y")])
  }
  cells <- expand.grid(
    algorithm = levels(droplevels(df$algorithm)),
    panel_x = levels(df$panel_x),
    panel_y = levels(df$panel_y)
  )
  present <- paste(df$algorithm, df$panel_x, df$panel_y)
  cells <- cells[
    !paste(cells$algorithm, cells$panel_x, cells$panel_y) %in%
      present,
  ]
  cells$algorithm <- factor(cells$algorithm, levels(droplevels(df$algorithm)))
  cells$panel_x <- factor(cells$panel_x, levels(df$panel_x))
  cells$panel_y <- factor(cells$panel_y, levels(df$panel_y))
  cells
}

# Methods left out of an accuracy book, named in its caption.
.excluded_note <- function(exclude) {
  if (length(exclude) == 0L) {
    return("")
  }
  paste0(
    "\n",
    paste(exclude, collapse = ", "),
    " not shown: its much larger error would squeeze the shared scale."
  )
}

# Drop excluded methods and their factor levels, so the facet grid
# closes up and the colour scale is computed from the methods that remain.
.drop_algorithms <- function(data, exclude) {
  data <- data[!data$algorithm %in% exclude, ]
  data$algorithm <- forcats::fct_drop(factor(data$algorithm))
  data
}

plot_metric_tiles <- function(data, metric, exclude = character()) {
  data <- .drop_algorithms(data, exclude)
  limits <- range(data[[metric]], finite = TRUE)
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ]
    df$value <- df[[metric]]
    failed <- df[!is.finite(df$value), ]
    missing <- .missing_cells(df)
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(.data$panel_x, .data$panel_y, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2)
    if (nrow(failed) > 0L) {
      p <- p +
        ggplot2::geom_tile(
          data = failed,
          ggplot2::aes(linetype = "no finite estimate in any bulk"),
          fill = "grey60",
          colour = "white",
          linewidth = 0.2
        ) +
        ggplot2::scale_linetype_manual(
          name = NULL,
          values = 1,
          guide = ggplot2::guide_legend(
            override.aes = list(fill = "grey60", colour = "grey30")
          )
        )
    }
    if (nrow(missing) > 0L) {
      p <- p +
        ggplot2::geom_text(
          data = missing,
          ggplot2::aes(.data$panel_x, .data$panel_y),
          label = "not\nsimulated",
          inherit.aes = FALSE,
          size = 2.4,
          colour = "grey45"
        )
    }
    p +
      ggplot2::facet_wrap(
        ~algorithm,
        ncol = ceiling(sqrt(dplyr::n_distinct(df$algorithm)))
      ) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_discrete(
        drop = FALSE,
        limits = if (isTRUE(df$grid_panels[[1L]])) rev
      ) +
      ggplot2::scale_fill_viridis_c(
        name = .metric_names[[metric]],
        limits = limits,
        trans = "log10",
        na.value = "grey60",
        guide = ggplot2::guide_colourbar(barwidth = grid::unit(7, "cm"))
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = df$panel_x_title[[1L]],
        y = df$panel_y_title[[1L]],
        title = paste(.metric_names[[metric]], "of every method"),
        subtitle = page,
        caption = paste0(
          "Mean over the simulated bulks of each scenario; one colour scale (log10) for every page.",
          "\nGrey tile: the method returned no finite estimate for any bulk of that scenario (solver failure).",
          .excluded_note(exclude)
        )
      ) +
      .page_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  })
}

plot_rmse_aitchison_dots <- function(data, exclude = character()) {
  data <- .drop_algorithms(data, exclude)
  rmse_limits <- range(data$rmse, finite = TRUE)
  size_limits <- range(data$aitchison, finite = TRUE)
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ]
    failed <- df[!is.finite(df$rmse), ]
    missing <- .missing_cells(df)
    p <- ggplot2::ggplot(
      df[is.finite(df$rmse), ],
      ggplot2::aes(.data$panel_x, .data$panel_y)
    ) +
      ggplot2::geom_point(ggplot2::aes(
        colour = .data$rmse,
        size = .data$aitchison
      ))
    if (nrow(failed) > 0L) {
      p <- p +
        ggplot2::geom_point(
          data = failed,
          ggplot2::aes(shape = "no finite estimate in any bulk"),
          size = 3,
          colour = "grey45"
        ) +
        ggplot2::scale_shape_manual(name = NULL, values = 4)
    }
    if (nrow(missing) > 0L) {
      p <- p +
        ggplot2::geom_text(
          data = missing,
          label = "not\nsimulated",
          size = 2.4,
          colour = "grey45"
        )
    }
    p +
      ggplot2::facet_wrap(
        ~algorithm,
        ncol = ceiling(sqrt(dplyr::n_distinct(df$algorithm)))
      ) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_discrete(
        drop = FALSE,
        limits = if (isTRUE(df$grid_panels[[1L]])) rev
      ) +
      ggplot2::scale_colour_viridis_c(
        name = "Mean RMSE",
        limits = rmse_limits,
        trans = "log10"
      ) +
      ggplot2::scale_size_continuous(
        name = "Mean Aitchison",
        limits = size_limits,
        range = c(0.8, 7)
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = df$panel_x_title[[1L]],
        y = df$panel_y_title[[1L]],
        title = "Accuracy of every method: colour RMSE, size Aitchison distance",
        subtitle = page,
        caption = paste0(
          "Means over the simulated bulks; colour and size scales shared by every page.",
          "\nCross: the method returned no finite estimate for any bulk of that scenario (solver failure).",
          .excluded_note(exclude)
        )
      ) +
      .page_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  })
}

# Interval coverage and width ----------------------------------------------------

.data_ci_forest <- function(benchmark, layout, figure) {
  intervals <- benchmark$intervals[benchmark$intervals$is_composition, ]
  intervals <- intervals[
    is.finite(intervals$mean_lower) & intervals$n_intervals > 0L,
  ]
  if (identical(figure, "fig02")) {
    intervals <- intervals[intervals$parameter == "celltype_1", ]
  }
  # Profile intervals refit from the bulk whatever the solver: keep them
  # once, under Newton-Raphson (benchmarks run before 2.5.0 repeat them).
  intervals <- intervals[
    !grepl("^profile", intervals$ci_method) |
      intervals$algorithm == "Newton-Raphson",
  ]
  intervals$algorithm <- .method_factor(intervals$algorithm)
  intervals$ci_label <- factor(
    .ci_labels[intervals$ci_method],
    unname(.ci_labels)
  )
  intervals <- .with_layout(intervals, layout)
  intervals[intervals$corner, ]
}

plot_ci_forest <- function(data) {
  data$key <- paste(data$ID, data$parameter)
  lapply(
    unique(data$key[order(data$page, data$panel, data$parameter)]),
    function(key) {
      df <- droplevels(data[data$key == key & !is.na(data$ci_label), ])
      if (nrow(df) == 0L) {
        return(NULL)
      }
      span <- range(
        c(df$mean_lower, df$mean_upper, df$mean_truth),
        finite = TRUE
      )
      # The row label carries the coverage, its Wilson Monte Carlo interval
      # and the mean width, so nothing is drawn outside the panel.
      df$row <- sprintf(
        "%s   %-4s [%4s, %4s]   w = %s",
        formatC(
          as.character(df$ci_label),
          width = -max(nchar(as.character(df$ci_label)))
        ),
        .signif_label(df$coverage, 2L),
        .signif_label(df$coverage_lower, 2L),
        .signif_label(df$coverage_upper, 2L),
        .signif_label(df$mean_width)
      )
      df$row <- stats::reorder(df$row, -as.integer(df$ci_label))
      ggplot2::ggplot(
        df,
        ggplot2::aes(y = .data$row, colour = .data$ci_label)
      ) +
        ggplot2::geom_vline(
          xintercept = df$mean_truth[[1L]],
          linetype = 2,
          colour = "grey30"
        ) +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$mean_lower, xmax = .data$mean_upper),
          width = 0.3,
          linewidth = 0.9,
          orientation = "y"
        ) +
        ggplot2::geom_point(ggplot2::aes(x = .data$mean_estimate), size = 2.4) +
        ggplot2::facet_grid(
          algorithm ~ .,
          scales = "free_y",
          space = "free_y"
        ) +
        ggplot2::scale_colour_manual(name = NULL, values = .ci_palette) +
        ggplot2::coord_cartesian(xlim = span) +
        ggplot2::labs(
          x = sprintf(
            "Proportion of %s",
            sub("celltype_", "cell type ", df$parameter[[1L]])
          ),
          y = NULL,
          title = "Confidence intervals: coverage [Wilson 95% interval] and mean width w",
          subtitle = sprintf("%s / %s", df$page[[1L]], df$panel[[1L]]),
          caption = paste(
            "Bar: mean lower to mean upper bound over the simulated bulks; point: mean estimate;",
            "dashed: mean target.\nMethods without intervals are not shown."
          )
        ) +
        .page_theme() +
        ggplot2::theme(
          strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
          axis.text.y = ggplot2::element_text(
            family = "mono",
            size = 9,
            hjust = 0
          ),
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank()
        )
    }
  )
}

# Intervals for the technical scales of the single-cell model, one page per
# scenario: mean scale c = a n, variance scale phi = a^2 n, platform scale
# a = phi / c and cells per bulk n = c^2 / phi, on log10 axes (every bound
# is positive).
.nuisance_parameters <- c(
  mean_scale = "c (mean scale, a n)",
  variance_scale = "phi (variance scale, a^2 n)",
  platform_scale_a = "a = phi / c (platform scale)",
  cell_count_n = "n = c^2 / phi (cells per bulk)"
)

.data_ci_forest_nuisance <- function(benchmark, layout) {
  intervals <- benchmark$intervals[
    !benchmark$intervals$is_composition &
      benchmark$intervals$parameter %in% names(.nuisance_parameters),
  ]
  intervals <- intervals[
    is.finite(intervals$mean_lower) &
      intervals$mean_lower > 0 &
      intervals$n_intervals > 0L,
  ]
  intervals$algorithm <- .method_factor(intervals$algorithm)
  intervals$ci_label <- factor(
    .ci_labels[intervals$ci_method],
    unname(.ci_labels)
  )
  intervals$parameter_label <- factor(
    .nuisance_parameters[intervals$parameter],
    unname(.nuisance_parameters)
  )
  .with_layout(intervals, layout)
}

plot_ci_forest_nuisance <- function(data) {
  ids <- unique(data$ID[order(data$page, data$panel)])
  lapply(ids, function(id) {
    df <- droplevels(data[data$ID == id & !is.na(data$ci_label), ])
    if (nrow(df) == 0L) {
      return(NULL)
    }
    df$row <- interaction(
      df$ci_label,
      df$algorithm,
      sep = " / ",
      lex.order = TRUE
    )
    df$row <- factor(df$row, rev(levels(droplevels(df$row))))
    df$label <- sprintf(
      "cov %s [%s, %s]",
      .signif_label(df$coverage, 2L),
      .signif_label(df$coverage_lower, 2L),
      .signif_label(df$coverage_upper, 2L)
    )
    truth <- dplyr::distinct(df, .data$parameter_label, .data$mean_truth)
    ggplot2::ggplot(df, ggplot2::aes(y = .data$row, colour = .data$ci_label)) +
      ggplot2::geom_vline(
        data = truth,
        ggplot2::aes(xintercept = .data$mean_truth),
        linetype = 2,
        colour = "grey30"
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$mean_lower, xmax = .data$mean_upper),
        width = 0.3,
        linewidth = 0.9,
        orientation = "y"
      ) +
      ggplot2::geom_point(ggplot2::aes(x = .data$mean_estimate), size = 2.2) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data$mean_estimate, label = .data$label),
        vjust = -0.9,
        size = 2.5,
        show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~parameter_label, ncol = 2L, scales = "free_x") +
      ggplot2::scale_x_log10() +
      ggplot2::scale_colour_manual(name = NULL, values = .ci_palette) +
      ggplot2::labs(
        x = "Value (log10)",
        y = NULL,
        title = "Confidence intervals for the technical scales: coverage [Wilson 95% interval]",
        subtitle = sprintf("%s / %s", df$page[[1L]], df$panel[[1L]]),
        caption = paste(
          "Bar: mean lower to mean upper bound over the simulated bulks; point: mean estimate; dashed: true value.",
          "\nIntervals formed on the log scale, so every bound is positive. n-hat is a number of cells only",
          "because the platform scale a is known to the generator; otherwise c and phi identify a n and a^2 n only."
        )
      ) +
      .page_theme() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        legend.position = "bottom",
        panel.grid.minor = ggplot2::element_blank()
      )
  })
}

# Monte Carlo law of the estimates ---------------------------------------------

.data_estimates <- function(benchmark, layout, figure) {
  long <- .estimates_long(benchmark, layout)
  long <- long[long$corner, ]
  if (identical(figure, "fig02")) {
    long <- long[long$cell_type == "celltype_1", ]
  }
  long$cell_type <- sub("celltype_", "Cell type ", long$cell_type)
  long
}

plot_estimate_laws <- function(data, style) {
  keys <- unique(data[c("page", "cell_type")])
  keys <- keys[order(keys$page, keys$cell_type), ]
  lapply(seq_len(nrow(keys)), function(k) {
    df <- data[
      data$page == keys$page[[k]] & data$cell_type == keys$cell_type[[k]],
    ]
    targets <- dplyr::summarise(
      df,
      truth = mean(.data$truth),
      .by = c("panel", "panel_x", "panel_y")
    )
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(.data$estimate, .data$algorithm, fill = .data$algorithm)
    ) +
      ggplot2::geom_vline(
        data = targets,
        ggplot2::aes(xintercept = .data$truth),
        linetype = 2,
        colour = "grey30"
      )
    p <- if (identical(style, "raincloud")) {
      p +
        ggdist::stat_slab(
          side = "top",
          scale = 0.6,
          alpha = 0.6,
          normalize = "groups",
          na.rm = TRUE
        ) +
        ggdist::stat_dots(
          side = "bottom",
          scale = 0.5,
          dotsize = 0.8,
          na.rm = TRUE
        ) +
        ggdist::stat_pointinterval(
          .width = c(0.5, 0.95),
          point_size = 1.4,
          na.rm = TRUE
        )
    } else {
      p +
        ggdist::stat_dotsinterval(
          quantiles = 50L,
          .width = c(0.5, 0.95),
          point_size = 1.4,
          na.rm = TRUE
        )
    }
    p +
      .panel_facets(df) +
      .not_simulated(df) +
      ggplot2::scale_fill_manual(values = .method_palette, guide = "none") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(
        x = paste("Estimated proportion,", tolower(keys$cell_type[[k]])),
        y = NULL,
        title = if (identical(style, "raincloud")) {
          "Monte Carlo law of the estimates: density, dots and quantile intervals"
        } else {
          "Monte Carlo law of the estimates: 50 quantile dots"
        },
        subtitle = keys$page[[k]],
        caption = paste(
          "One estimate per simulated bulk; intervals: 50% and 95% Monte Carlo quantiles;",
          "dashed: mean target. Mass at 0 or 1 is an estimate on a face."
        )
      ) +
      .page_theme()
  })
}

# Wald approximation ---------------------------------------------------------------

# Estimates of the likelihood methods in Helmert ILR coordinates, whitened
# by the information at the true composition: I_z(p*)^{1/2} (z-hat - z*).
# Population: expected information; single cell: observed profile
# information of the bulk at its realised composition. Estimates on a face
# have no ILR coordinates and are dropped.
.data_whitened <- function(benchmark, layout, figure) {
  cts <- .celltype_names(benchmark)
  opt <- benchmark$optimisation
  opt <- opt[
    opt$algorithm %in%
      .likelihood_methods &
      opt$ID %in% layout$ID[layout$corner],
  ]
  info_population <- list()
  rows <- lapply(seq_len(nrow(opt)), function(i) {
    id <- opt$ID[[i]]
    s <- as.integer(sub("sample_", "", opt$sample_id[[i]]))
    th <- benchmark$theta_true[[id]]
    p_hat <- as.numeric(opt[i, cts])
    if (!all(is.finite(p_hat)) || min(p_hat) <= 0) {
      return(NULL)
    }
    if (identical(figure, "fig05")) {
      p_true <- benchmark$data[[id]]$truth[cts, s]
      z_true <- isometric_log_ratio(p_true)
      info <- -hessian_profile_loglik_single_cell_ilr(
        z_true,
        benchmark$data[[id]]$bulk[, s],
        th$mu,
        th$sigma
      )
    } else {
      p_true <- as.numeric(th$p)
      z_true <- isometric_log_ratio(p_true)
      if (is.null(info_population[[id]])) {
        jac <- jacobian_isometric_logistic(z_true)
        info_population[[id]] <<- crossprod(
          jac,
          expected_fisher_unconstrained(p_true, th$mu, th$sigma) %*% jac
        )
      }
      info <- info_population[[id]]
    }
    root <- tryCatch(chol(info), error = function(e) NULL)
    if (is.null(root)) {
      return(NULL)
    }
    w <- drop(root %*% (isometric_log_ratio(p_hat) - z_true))
    tibble::tibble(
      ID = id,
      algorithm = opt$algorithm[[i]],
      coordinate = paste0("z", seq_along(w)),
      whitened = w,
      mahalanobis = sum(w^2),
      df = length(w),
      outcome = simulation_outcome(
        opt$numerical_converged[[i]],
        opt$theoretical_converged[[i]]
      )
    )
  })
  out <- dplyr::bind_rows(rows)
  out$algorithm <- .method_factor(out$algorithm)
  out$outcome <- factor(out$outcome, names(.outcome_labels))
  .with_layout(out, layout)
}

.outcome_labels <- c(
  theoretical_success = "Converged, local maximum",
  theoretical_failure = "Converged, not a local maximum",
  numerical_failure = "Numerical failure"
)
.outcome_palette <- c(
  theoretical_success = "#1B9E77",
  theoretical_failure = "#E6AB02",
  numerical_failure = "#D95F02"
)

plot_wald_normality_qq <- function(data) {
  limit <- max(4, stats::quantile(abs(data$whitened), 0.99, na.rm = TRUE))
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ]
    df <- df |>
      dplyr::group_by(.data$panel, .data$algorithm, .data$coordinate) |>
      dplyr::mutate(
        theoretical = stats::qnorm(stats::ppoints(dplyr::n()))[rank(
          .data$whitened,
          ties.method = "first"
        )]
      ) |>
      dplyr::ungroup()
    ggplot2::ggplot(df, ggplot2::aes(.data$theoretical, .data$whitened)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data$outcome, shape = .data$coordinate),
        size = 1.8
      ) +
      ggplot2::facet_grid(algorithm ~ panel) +
      ggplot2::scale_colour_manual(
        name = NULL,
        values = .outcome_palette,
        labels = .outcome_labels,
        drop = FALSE
      ) +
      ggplot2::scale_shape_discrete(name = "ILR coordinate") +
      ggplot2::coord_cartesian(xlim = c(-3, 3), ylim = c(-limit, limit)) +
      ggplot2::labs(
        x = "N(0, 1) quantile",
        y = "Whitened ILR estimate",
        title = "Wald approximation: whitened ILR estimates against N(0, 1)",
        subtitle = page,
        caption = paste(
          "I_z(p*)^(1/2) (z-hat - z*), information at the true composition;",
          "estimates on a face have no ILR coordinates and are not shown.",
          "\nAssesses the Wald approximation only; profile and bootstrap intervals do not rely on it."
        )
      ) +
      .page_theme()
  })
}

plot_wald_chi2_qq <- function(data) {
  data <- dplyr::distinct(
    data,
    .data$ID,
    .data$page,
    .data$panel,
    .data$algorithm,
    .data$mahalanobis,
    .data$df,
    .data$outcome
  )
  limits <- range(c(1e-3, data$mahalanobis, 30), finite = TRUE)
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ] |>
      dplyr::group_by(.data$panel, .data$algorithm) |>
      dplyr::mutate(
        n = dplyr::n(),
        theoretical = stats::qchisq(stats::ppoints(.data$n), .data$df)[rank(
          .data$mahalanobis,
          ties.method = "first"
        )]
      ) |>
      dplyr::ungroup()
    sizes <- unique(df[c("n", "df")])
    envelope <- purrr::map_dfr(seq_len(nrow(sizes)), function(k) {
      n <- sizes$n[[k]]
      sims <- withr::with_seed(1L, stats::rchisq(n * 2000L, sizes$df[[k]]))
      sims <- matrix(apply(matrix(sims, n), 2L, sort), nrow = n)
      tibble::tibble(
        n = n,
        df = sizes$df[[k]],
        theoretical = stats::qchisq(stats::ppoints(n), sizes$df[[k]]),
        lower = apply(sims, 1L, stats::quantile, 0.025),
        upper = apply(sims, 1L, stats::quantile, 0.975)
      )
    })
    bands <- dplyr::inner_join(
      dplyr::distinct(df, .data$panel, .data$algorithm, .data$n, .data$df),
      envelope,
      by = c("n", "df"),
      relationship = "many-to-many"
    )
    ggplot2::ggplot(df, ggplot2::aes(.data$theoretical, .data$mahalanobis)) +
      ggplot2::geom_ribbon(
        data = bands,
        ggplot2::aes(
          x = .data$theoretical,
          ymin = .data$lower,
          ymax = .data$upper
        ),
        inherit.aes = FALSE,
        fill = "grey80"
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
      ggplot2::geom_point(ggplot2::aes(colour = .data$outcome), size = 1.8) +
      ggplot2::facet_grid(algorithm ~ panel) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::coord_cartesian(xlim = limits, ylim = limits) +
      ggplot2::scale_colour_manual(
        name = NULL,
        values = .outcome_palette,
        labels = .outcome_labels,
        drop = FALSE
      ) +
      ggplot2::labs(
        x = "chi-square quantile, J - 1 degrees of freedom (log10)",
        y = "Squared Mahalanobis norm D^2 (log10)",
        title = "Wald approximation: D^2 of the ILR estimate against chi-square",
        subtitle = page,
        caption = paste(
          "Grey band: 95% pointwise Monte Carlo envelope of sorted chi-square samples of the same size.",
          "\nLog axes keep the bulk of the points visible next to the large values of estimates near a face."
        )
      ) +
      .page_theme()
  })
}

# Agreement between methods -------------------------------------------------------

# Mean Hellinger distance between the paired estimates of two methods over
# the bulks of a page: ||sqrt(p_a) - sqrt(p_b)|| / sqrt(2), Euclidean in
# sqrt(p) (as Ward's criterion requires) and defined at exact zeros.
.data_agreement <- function(benchmark, layout) {
  cts <- .celltype_names(benchmark)
  opt <- .with_layout(benchmark$optimisation, layout)
  root <- sqrt(pmax(as.matrix(opt[cts]), 0))
  key <- paste(opt$ID, opt$sample_id)
  methods <- levels(.method_factor(opt$algorithm))
  purrr::map_dfr(levels(droplevels(opt$page)), function(page) {
    on_page <- opt$page == page
    pairs <- expand.grid(a = methods, b = methods, stringsAsFactors = FALSE)
    pairs$hellinger <- mapply(
      function(a, b) {
        ia <- which(on_page & opt$algorithm == a)
        ib <- which(on_page & opt$algorithm == b)
        ib <- ib[match(key[ia], key[ib])]
        d <- sqrt(rowSums(
          (root[ia, , drop = FALSE] - root[ib, , drop = FALSE])^2
        )) /
          sqrt(2)
        mean(d, na.rm = TRUE)
      },
      pairs$a,
      pairs$b
    )
    tibble::tibble(page = page, pairs)
  })
}

plot_estimate_agreement <- function(data) {
  limits <- c(0, max(data$hellinger, na.rm = TRUE))
  lapply(unique(data$page), function(page) {
    df <- data[data$page == page, ]
    distance <- stats::xtabs(hellinger ~ a + b, df)
    keep <- rownames(distance)[
      is.finite(diag(distance)) & rowSums(!is.finite(distance)) == 0
    ]
    distance <- distance[keep, keep]
    hc <- stats::hclust(stats::as.dist(distance), method = "ward.D2")
    ordered <- rownames(distance)[hc$order]
    df <- df[df$a %in% keep & df$b %in% keep, ]
    df$a <- factor(df$a, ordered)
    df$b <- factor(df$b, rev(ordered))
    n <- length(ordered)
    pairwise <- distance[upper.tri(distance)]
    average <- mean(pairwise[is.finite(pairwise)])
    tree <- ggplot2::ggplot(.dendrogram_segments(hc)) +
      ggplot2::geom_hline(
        yintercept = average,
        linetype = 2,
        colour = "#C45C26"
      ) +
      ggplot2::annotate(
        "label",
        x = n + 0.45,
        y = average,
        label = paste("mean Hellinger =", .signif_label(average)),
        hjust = 1,
        vjust = -0.3,
        size = 3.6,
        colour = "#C45C26",
        fill = "white",
        linewidth = 0
      ) +
      ggplot2::geom_segment(ggplot2::aes(
        .data$x,
        .data$y,
        xend = .data$xend,
        yend = .data$yend
      )) +
      ggplot2::scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      ggplot2::labs(
        title = "Agreement of the estimates: mean Hellinger distance, Ward clustering",
        subtitle = sprintf(
          "%s / mean pairwise Hellinger distance between methods = %s",
          page,
          .signif_label(average)
        )
      ) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.subtitle = ggplot2::element_text(size = 11)
      )
    tiles <- ggplot2::ggplot(
      df,
      ggplot2::aes(.data$a, .data$b, fill = .data$hellinger)
    ) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(
        ggplot2::aes(
          label = .signif_label(.data$hellinger, 2L),
          colour = .data$hellinger > mean(limits)
        ),
        size = 3.2
      ) +
      ggplot2::scale_colour_manual(
        values = c(`FALSE` = "black", `TRUE` = "white"),
        guide = "none"
      ) +
      ggplot2::scale_fill_viridis_c(
        name = "Mean Hellinger distance",
        limits = limits,
        direction = -1,
        end = 0.9
      ) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        caption = paste(
          "Distance between the estimates of two methods on the same bulk, averaged over the bulks of the page;",
          "\n0 means identical estimates. Ward.D2 on these distances; one colour scale for every page.",
          "\nDashed line on the dendrogram: mean of the pairwise distances between methods on this page."
        )
      ) +
      .page_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        legend.key.width = grid::unit(2, "cm")
      )
    cowplot::plot_grid(
      tree,
      tiles,
      ncol = 1L,
      align = "v",
      axis = "lr",
      rel_heights = c(0.3, 1)
    ) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
  })
}

# Convergence and cost ----------------------------------------------------------

.data_optimisation <- function(benchmark, layout) {
  opt <- .with_layout(benchmark$optimisation, layout)
  opt <- opt[opt$corner, ]
  opt$algorithm <- .method_factor(opt$algorithm)
  opt$outcome <- factor(
    simulation_outcome(opt$numerical_converged, opt$theoretical_converged),
    names(.outcome_labels)
  )
  opt$memory_mib <- opt$memory_bytes / 2^20
  opt[c(
    "ID",
    "page",
    "panel",
    "panel_x",
    "panel_y",
    "grid_panels",
    "algorithm",
    "outcome",
    "kkt_residual",
    "elapsed_sec",
    "memory_mib"
  )]
}

plot_kkt_raincloud <- function(data) {
  data$kkt <- pmax(data$kkt_residual, 1e-12)
  limits <- range(data$kkt, finite = TRUE)
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page & is.finite(data$kkt), ]
    ggplot2::ggplot(
      df,
      ggplot2::aes(.data$kkt, .data$algorithm, fill = .data$algorithm)
    ) +
      ggdist::stat_slab(
        side = "top",
        scale = 0.6,
        alpha = 0.6,
        normalize = "groups"
      ) +
      ggdist::stat_dots(side = "bottom", scale = 0.5, dotsize = 0.8) +
      .panel_facets(df) +
      .not_simulated(df) +
      ggplot2::scale_x_log10(limits = limits) +
      ggplot2::scale_fill_manual(values = .method_palette, guide = "none") +
      ggplot2::labs(
        x = "KKT residual of the DeCovarT log-likelihood (log10, floored at 1e-12)",
        y = NULL,
        title = "First-order optimality of every estimate",
        subtitle = page,
        caption = "Projected score at the estimate on the simplex; 0 at a maximum of the DeCovarT likelihood."
      ) +
      .page_theme()
  })
}

plot_convergence <- function(data) {
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ]
    ggplot2::ggplot(df, ggplot2::aes(.data$algorithm, fill = .data$outcome)) +
      ggplot2::geom_bar(width = 0.75) +
      .panel_facets(df) +
      .not_simulated(df) +
      ggplot2::scale_fill_manual(
        name = NULL,
        values = .outcome_palette,
        labels = .outcome_labels,
        drop = FALSE
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Simulated bulks",
        title = "Convergence: theoretical success, theoretical and numerical failures",
        subtitle = page,
        caption = paste(
          "Numerical failure: no finite estimate on the simplex. Theoretical failure: a finite estimate",
          "that is not a local maximum of the DeCovarT likelihood."
        )
      ) +
      .page_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  })
}

plot_resource <- function(data, metric) {
  spec <- list(
    elapsed_sec = c("Run time per bulk (seconds, log10)", "Run time"),
    memory_mib = c("Peak memory per bulk (MiB, log10)", "Memory")
  )[[metric]]
  data$value <- data[[metric]]
  data <- data[is.finite(data$value) & data$value > 0, ]
  limits <- range(data$value)
  lapply(levels(droplevels(data$page)), function(page) {
    df <- data[data$page == page, ]
    ggplot2::ggplot(
      df,
      ggplot2::aes(.data$value, .data$algorithm, colour = .data$algorithm)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(height = 0.15, seed = 1L),
        size = 1,
        alpha = 0.6
      ) +
      .panel_facets(df) +
      .not_simulated(df) +
      ggplot2::scale_x_log10(limits = limits) +
      ggplot2::scale_colour_manual(values = .method_palette, guide = "none") +
      ggplot2::labs(
        x = spec[[1L]],
        y = NULL,
        title = paste(spec[[2L]], "of every method"),
        subtitle = page,
        caption = "Point estimate only; interval computations are not included. Competitors: batch time divided by the bulks."
      ) +
      .page_theme()
  })
}

# Nuisance scales (single cell) ----------------------------------------------------

# Estimated over true mean scale c = a n and variance scale phi = a^2 n, and
# phi rescaled by G / (G - J) (REML degrees of freedom).
.data_nuisance <- function(benchmark, layout) {
  opt <- benchmark$optimisation
  opt <- opt[
    opt$algorithm %in% .likelihood_methods & is.finite(opt$mean_scale),
  ]
  s <- as.integer(sub("sample_", "", opt$sample_id))
  truth <- t(mapply(
    function(id, k) benchmark$data[[id]]$nuisance[, k],
    opt$ID,
    s
  ))
  J <- length(.celltype_names(benchmark))
  opt <- .with_layout(opt, layout)
  ratio <- tibble::tibble(
    ID = opt$ID,
    algorithm = .method_factor(opt$algorithm),
    n_genes = layout$n_genes[match(opt$ID, layout$ID)],
    n_cells = layout$n_cells[match(opt$ID, layout$ID)],
    platform_scale = layout$platform_scale[match(opt$ID, layout$ID)],
    `c-hat / c` = opt$mean_scale / truth[, "mean_scale"],
    `phi-hat / phi` = opt$variance_scale / truth[, "variance_scale"],
    `phi-hat G / (G - J) / phi` = `phi-hat / phi` * n_genes / (n_genes - J)
  )
  tidyr::pivot_longer(
    ratio,
    dplyr::ends_with("phi") | dplyr::ends_with("/ c"),
    names_to = "ratio"
  )
}

plot_nuisance_ratio <- function(data) {
  data$cells <- factor(sprintf("n = %g cells", data$n_cells))
  data$genes <- factor(data$n_genes)
  data$ratio <- factor(
    data$ratio,
    c("c-hat / c", "phi-hat / phi", "phi-hat G / (G - J) / phi")
  )
  list(
    ggplot2::ggplot(
      data,
      ggplot2::aes(.data$genes, .data$value, colour = .data$algorithm)
    ) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
      ggplot2::geom_boxplot(
        outlier.shape = NA,
        position = ggplot2::position_dodge(0.8),
        width = 0.6
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = 0.8,
          seed = 1L
        ),
        size = 0.8,
        alpha = 0.5
      ) +
      ggplot2::facet_grid(ratio ~ cells, scales = "free_y") +
      ggplot2::scale_colour_manual(name = NULL, values = .method_palette) +
      ggplot2::labs(
        x = "Genes G",
        y = "Estimate / truth",
        title = "Nuisance scales: mean scale c = a n and variance scale phi = a^2 n",
        caption = paste(
          "All scenarios pooled; the ML phi-hat is biased by (G - J) / G, removed by the REML rescaling (bottom row).",
          "\nc-hat / c also absorbs the multinomial noise of the realised composition."
        )
      ) +
      .page_theme()
  )
}
