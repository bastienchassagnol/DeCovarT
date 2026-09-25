# Manuscript multi-page PDF books (not part of the installed package).
# Loaded into the DeCovarT namespace by attach_figure_books().
#
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
