#' Convolution covariance \eqn{\Sigma(p)=\sum_j p_j^2 \Sigma_j}
#'
#' @keywords internal
#' @noRd
.sigma_of_p <- function(p, Sigma) {
  .compute_global_variance(as.numeric(p), Sigma)
}

#' Draw a 2-D Gaussian sample as a data frame
#'
#' @keywords internal
#' @noRd
.rnorm_2d <- function(n, mu, sigma, component) {
  draws <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma, empirical = FALSE)
  if (is.null(dim(draws))) {
    draws <- matrix(draws, nrow = 1L)
  }
  data.frame(
    gene_1 = draws[, 1L],
    gene_2 = draws[, 2L],
    component = component
  )
}

#' Four correlation corners used on each bivariate density page
#'
#' @keywords internal
#' @noRd
.bivariate_corr_corners <- function() {
  list(
    "rho = (0, 0)" = c(0, 0),
    "rho = (-0.8, -0.8)" = c(-0.8, -0.8),
    "rho = (0.8, 0.8)" = c(0.8, 0.8),
    "rho = (-0.8, 0.8)" = c(-0.8, 0.8)
  )
}

#' Exact Gaussian probability ellipse (known mean and covariance)
#'
#' Boundary of the set
#' \eqn{(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^2_{2,\alpha}}
#' for a bivariate Gaussian with **known** \eqn{\mu} and \eqn{\Sigma}.
#' The squared Mahalanobis distance is exactly \eqn{\chi^2_2}.
#'
#' @param mu Length-2 mean.
#' @param Sigma \eqn{2 \times 2} covariance.
#' @param level Probability content (default 0.95).
#' @param n Number of boundary points.
#'
#' @return A data frame with `x` and `y`.
#' @examples
#' mu <- c(20, 22)
#' Sigma <- matrix(c(1, 0.4, 0.4, 1), 2)
#' head(gaussian_confidence_ellipse(mu, Sigma))
#' @export
gaussian_confidence_ellipse <- function(
  mu,
  Sigma,
  level = 0.95,
  n = 200L
) {
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  if (length(mu) != 2L || !all(dim(Sigma) == c(2L, 2L))) {
    stop(
      "gaussian_confidence_ellipse() is defined for G = 2.",
      call. = FALSE
    )
  }
  cutoff <- stats::qchisq(level, df = 2L)
  eig <- eigen(Sigma, symmetric = TRUE)
  ev <- pmax(eig$values, 0)
  sigma_half <- eig$vectors %*%
    diag(sqrt(ev), nrow = 2L) %*%
    t(eig$vectors)
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(theta), sin(theta))
  boundary <- matrix(mu, nrow = 2L, ncol = n) +
    sqrt(cutoff) * sigma_half %*% circle
  data.frame(x = boundary[1L, ], y = boundary[2L, ])
}

#' @keywords internal
#' @noRd
.celltype_names <- function(true_theta) {
  nms <- colnames(true_theta$mu)
  if (is.null(nms)) {
    nms <- paste0("celltype_", seq_len(ncol(true_theta$mu)))
  }
  nms
}

#' @keywords internal
#' @noRd
.celltype_colour_values <- function(cts) {
  pal <- c("#E41A1C", "#4DAF4A")
  stats::setNames(rep(pal, length.out = length(cts)), cts)
}

#' @keywords internal
#' @noRd
.celltype_shape_values <- function(cts) {
  shp <- c(16L, 17L)
  stats::setNames(rep(shp, length.out = length(cts)), cts)
}

#' Centroids and 95% ellipses for each cell type
#'
#' @keywords internal
#' @noRd
.celltype_overlay <- function(true_theta, panel = NULL) {
  cts <- .celltype_names(true_theta)
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  centroids <- data.frame(
    gene_1 = as.numeric(mu[1L, ]),
    gene_2 = as.numeric(mu[2L, ]),
    cell_type = factor(cts, levels = cts)
  )
  ellipses <- lapply(seq_along(cts), function(j) {
    ell <- gaussian_confidence_ellipse(mu[, j], Sigma[,, j])
    data.frame(
      gene_1 = ell$x,
      gene_2 = ell$y,
      cell_type = factor(cts[[j]], levels = cts)
    )
  })
  ellipses <- dplyr::bind_rows(ellipses)
  if (!is.null(panel)) {
    centroids$panel <- panel
    ellipses$panel <- panel
  }
  list(centroids = centroids, ellipses = ellipses)
}

#' @keywords internal
#' @noRd
.round_interval_labels <- function(x, digits = 2L) {
  # R TRE: use [^]] not [^\\]] (the latter never matches).
  vapply(
    as.character(x),
    function(lab) {
      m <- regexec("^\\(([^,]+),\\s*([^]]+)\\]$", lab)
      parts <- regmatches(lab, m)[[1L]]
      if (length(parts) < 3L) {
        return(lab)
      }
      lo <- format(
        round(as.numeric(parts[[2L]]), digits),
        nsmall = digits,
        trim = TRUE
      )
      hi <- format(
        round(as.numeric(parts[[3L]]), digits),
        nsmall = digits,
        trim = TRUE
      )
      paste0("(", lo, ", ", hi, "]")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' @keywords internal
#' @noRd
.scale_fill_ndensity <- function(name = "Density") {
  ggplot2::scale_fill_viridis_d(
    name = name,
    labels = function(breaks) {
      out <- .round_interval_labels(breaks, digits = 2L)
      stats::setNames(out, as.character(breaks))
    }
  )
}

#' @keywords internal
#' @noRd
.theme_density_2d <- function() {
  theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(1, 2, 1, 1),
      panel.spacing = grid::unit(0.1, "lines"),
      panel.background = ggplot2::element_rect(
        fill = "#440154",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0)
    )
}

#' @keywords internal
#' @noRd
.pad_limits <- function(x, pad = 0.06) {
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) {
    return(c(rng[[1L]] - 1, rng[[1L]] + 1))
  }
  rng + c(-1, 1) * pad * span
}

#' Shared gene-axis limits from means, covariances, and ellipses
#'
#' @keywords internal
#' @noRd
.gene_axis_limits <- function(thetas) {
  xs <- numeric(0)
  ys <- numeric(0)
  for (th in thetas) {
    if (is.null(th)) {
      next
    }
    ov <- .celltype_overlay(th)
    xs <- c(xs, ov$centroids$gene_1, ov$ellipses$gene_1)
    ys <- c(ys, ov$centroids$gene_2, ov$ellipses$gene_2)
    mu <- th$mu
    Sigma <- th$sigma
    for (j in seq_len(ncol(mu))) {
      sds <- sqrt(pmax(diag(Sigma[,, j]), 0))
      xs <- c(xs, mu[1L, j] + c(-4, 4) * sds[[1L]])
      ys <- c(ys, mu[2L, j] + c(-4, 4) * sds[[2L]])
    }
  }
  list(xlim = .pad_limits(xs), ylim = .pad_limits(ys))
}

#' @keywords internal
#' @noRd
.purified_density_draws <- function(true_theta, n, panel = NULL) {
  cts <- .celltype_names(true_theta)
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  parts <- lapply(seq_along(cts), function(j) {
    .rnorm_2d(n, mu[, j], Sigma[,, j], cts[[j]])
  })
  df <- dplyr::bind_rows(parts)
  df$cell_type <- factor(df$component, levels = cts)
  if (!is.null(panel)) {
    df$panel <- panel
  }
  df
}

#' @keywords internal
#' @noRd
.bulk_density_draws <- function(true_theta, n, panel = NULL) {
  p <- as.numeric(true_theta$p)
  mu_p <- as.numeric(true_theta$mu %*% p)
  sigma_p <- .sigma_of_p(p, true_theta$sigma)
  df <- .rnorm_2d(n, mu_p, sigma_p, "bulk")
  if (!is.null(panel)) {
    df$panel <- panel
  }
  df
}

#' @keywords internal
#' @noRd
.add_celltype_overlay_layers <- function(
  p,
  overlay,
  cts,
  ellipses = TRUE
) {
  pal <- .celltype_colour_values(cts)
  shp <- .celltype_shape_values(cts)
  if (
    isTRUE(ellipses) &&
      !is.null(overlay$ellipses) &&
      nrow(overlay$ellipses) > 0L
  ) {
    p <- p +
      ggplot2::geom_path(
        data = overlay$ellipses,
        ggplot2::aes(
          x = .data[["gene_1"]],
          y = .data[["gene_2"]],
          colour = .data[["cell_type"]],
          group = .data[["cell_type"]]
        ),
        inherit.aes = FALSE,
        linewidth = 1.35
      )
  }
  p +
    ggplot2::geom_point(
      data = overlay$centroids,
      ggplot2::aes(
        x = .data[["gene_1"]],
        y = .data[["gene_2"]],
        colour = .data[["cell_type"]],
        shape = .data[["cell_type"]]
      ),
      inherit.aes = FALSE,
      size = 3.2
    ) +
    ggplot2::scale_colour_manual(
      name = "Cell type",
      values = pal,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = "Cell type",
      values = shp,
      drop = FALSE
    )
}

#' Kernel density raster covering the full gene-axis window
#'
#' @keywords internal
#' @noRd
.kde2d_raster <- function(df, xlim, ylim, n = 80L) {
  if (nrow(df) < 2L) {
    return(data.frame(
      gene_1 = numeric(0),
      gene_2 = numeric(0),
      density = numeric(0)
    ))
  }
  kd <- MASS::kde2d(
    df$gene_1,
    df$gene_2,
    n = n,
    lims = c(xlim, ylim)
  )
  grid <- expand.grid(
    gene_1 = kd$x,
    gene_2 = kd$y,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$density <- as.vector(kd$z)
  grid
}

#' Faceted (or single) 2-D density raster over shared gene-axis limits
#'
#' @keywords internal
#' @noRd
.density_raster_table <- function(df, lims, n = 80L) {
  if (!"panel" %in% names(df)) {
    return(.kde2d_raster(df, lims$xlim, lims$ylim, n = n))
  }
  pieces <- lapply(split(df, df$panel), function(part) {
    out <- .kde2d_raster(part, lims$xlim, lims$ylim, n = n)
    if (nrow(out) == 0L) {
      return(out)
    }
    out$panel <- part$panel[[1L]]
    out
  })
  dplyr::bind_rows(pieces)
}

#' 2-D density of purified Gaussians
#'
#' A kernel-density raster fills the shared gene-axis window so panels
#' have no inner white frame. Cell-type 1 is a red circle; cell-type 2
#' is a green triangle. Outlines are exact 95% Gaussian ellipses for
#' the known means and covariances ([gaussian_confidence_ellipse()]):
#' \eqn{(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^{2}_{2,0.95}}.
#'
#' @param true_theta List with `p`, `mu`, `sigma` for \eqn{G=2}.
#' @param n Draws **per cell type**.
#'
#' @return A `ggplot`.
#' @export
plot_purified_density_2d <- function(true_theta, n = 800L) {
  cts <- .celltype_names(true_theta)
  df <- .purified_density_draws(true_theta, n)
  overlay <- .celltype_overlay(true_theta)
  lims <- .gene_axis_limits(list(true_theta))
  raster_df <- .kde2d_raster(df, lims$xlim, lims$ylim)
  p <- ggplot2::ggplot(
    raster_df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[["density"]]),
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "Density",
      na.value = "#440154"
    ) +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = "Purified Gaussians",
      caption = .purified_ellipse_caption()
    ) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  .add_celltype_overlay_layers(p, overlay, cts, ellipses = TRUE)
}

#' @keywords internal
#' @noRd
.purified_ellipse_caption <- function() {
  paste(
    "Footnote: outlines are exact 95% Gaussian regions",
    "(x - mu)^T Sigma^{-1} (x - mu) <= chi^2_{2, 0.95}",
    "(stats::qchisq(0.95, df = 2); gaussian_confidence_ellipse())."
  )
}

#' 2-D density of the bulk convolution \eqn{y\sim N(\mu p,\Sigma(p))}
#'
#' Centroids of the purified components are marked; 95% ellipses are
#' omitted so the convolution density fills the panel.
#'
#' @inheritParams plot_purified_density_2d
#' @return A `ggplot`.
#' @export
plot_bulk_convolution_density_2d <- function(true_theta, n = 1200L) {
  cts <- .celltype_names(true_theta)
  df <- .bulk_density_draws(true_theta, n)
  overlay <- .celltype_overlay(true_theta)
  lims <- .gene_axis_limits(list(true_theta))
  raster_df <- .kde2d_raster(df, lims$xlim, lims$ylim)
  p <- ggplot2::ggplot(
    raster_df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[["density"]]),
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "Density",
      na.value = "#440154"
    ) +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = "Bulk convolution"
    )
  .add_celltype_overlay_layers(p, overlay, cts, ellipses = FALSE)
}

#' @keywords internal
#' @noRd
.plot_gene_density_page <- function(
  thetas,
  n,
  which,
  title,
  share_legend = FALSE
) {
  keep <- !vapply(thetas, is.null, logical(1))
  thetas <- thetas[keep]
  panels <- factor(names(thetas), levels = names(thetas))
  cts <- .celltype_names(thetas[[1L]])
  draw_fun <- if (identical(which, "purified")) {
    .purified_density_draws
  } else {
    .bulk_density_draws
  }
  df <- dplyr::bind_rows(
    lapply(seq_along(thetas), function(i) {
      draw_fun(thetas[[i]], n, panel = as.character(panels[[i]]))
    })
  )
  df$panel <- factor(df$panel, levels = levels(panels))
  overlay <- list(
    centroids = dplyr::bind_rows(
      lapply(seq_along(thetas), function(i) {
        .celltype_overlay(
          thetas[[i]],
          panel = as.character(panels[[i]])
        )$centroids
      })
    ),
    ellipses = dplyr::bind_rows(
      lapply(seq_along(thetas), function(i) {
        .celltype_overlay(
          thetas[[i]],
          panel = as.character(panels[[i]])
        )$ellipses
      })
    )
  )
  overlay$centroids$panel <- factor(
    overlay$centroids$panel,
    levels = levels(panels)
  )
  overlay$ellipses$panel <- factor(
    overlay$ellipses$panel,
    levels = levels(panels)
  )
  lims <- .gene_axis_limits(thetas)
  raster_df <- .density_raster_table(df, lims)
  raster_df$panel <- factor(raster_df$panel, levels = levels(panels))
  overlay_ellipses <- identical(which, "purified")
  p <- ggplot2::ggplot(
    raster_df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[["density"]]),
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "Density",
      na.value = "#440154"
    ) +
    ggplot2::facet_wrap(~panel, ncol = 2L) +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = title,
      caption = if (isTRUE(overlay_ellipses)) {
        .purified_ellipse_caption()
      } else {
        NULL
      }
    )
  if (isTRUE(overlay_ellipses)) {
    p <- p +
      ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)
      )
  }
  p <- .add_celltype_overlay_layers(
    p,
    overlay,
    cts,
    ellipses = overlay_ellipses
  )
  ggplot_data <- list(
    raster = raster_df,
    centroids = overlay$centroids,
    ellipses = overlay$ellipses
  )
  if (!isTRUE(share_legend)) {
    attr(p, "ggplot_data") <- ggplot_data
    return(p)
  }
  out <- .attach_cowplot_legend(p)
  attr(out, "ggplot_data") <- ggplot_data
  out
}

#' Shared legend to the right of a faceted ggplot
#'
#' @keywords internal
#' @noRd
.attach_cowplot_legend <- function(p) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    return(p)
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
  cowplot::plot_grid(
    p_main,
    legend,
    nrow = 1L,
    rel_widths = c(4.2, 1)
  )
}

#' Log-likelihood lattice over hypothesised cell-type ratios
#'
#' Evaluates [loglik_multivariate()] of \(y=\mu p^{\star}\) on a grid
#' of \((p_1,p_2)\in(0,1)^2\) (unnormalised mixture weights).
#'
#' @keywords internal
#' @noRd
.proportion_loglik_lattice <- function(true_theta, grid = 50L, y = NULL) {
  p_true <- as.numeric(true_theta$p)
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (is.null(y)) {
    y <- drop(mu %*% p_true)
  }
  y <- as.numeric(y)
  p_seq <- seq(0.02, 0.98, length.out = grid)
  z <- matrix(NA_real_, grid, grid)
  for (i in seq_len(grid)) {
    for (j in seq_len(grid)) {
      z[i, j] <- tryCatch(
        loglik_multivariate(c(p_seq[[i]], p_seq[[j]]), y, mu, Sigma),
        error = function(e) NA_real_
      )
    }
  }
  z_true <- tryCatch(
    loglik_multivariate(p_true, y, mu, Sigma),
    error = function(e) NA_real_
  )
  list(
    p1 = p_seq,
    p2 = p_seq,
    z = z,
    p_true = p_true,
    z_true = z_true,
    y = y
  )
}

#' Contour of the bulk log-likelihood on a \eqn{(p_1,p_2)} lattice
#'
#' Axes are hypothesised cell-type ratios, not gene expression. The
#' true simulation proportions (MLE for \eqn{y=\mu p^{\star}}) are marked;
#' the dashed line is the simplex \eqn{p_1+p_2=1}.
#'
#' @inheritParams plot_purified_density_2d
#' @param grid Length of the lattice per axis.
#' @param y Optional bulk observation. Default is \eqn{\mu p^{\star}}.
#'
#' @return A `ggplot`.
#' @export
plot_bulk_loglik_surface_p <- function(
  true_theta,
  grid = 50L,
  y = NULL
) {
  lat <- .proportion_loglik_lattice(true_theta, grid = grid, y = y)
  grid_df <- expand.grid(p1 = lat$p1, p2 = lat$p2)
  grid_df$loglik <- as.vector(lat$z)
  true_df <- data.frame(
    p1 = lat$p_true[[1L]],
    p2 = lat$p_true[[2L]]
  )
  ggplot2::ggplot(
    grid_df,
    ggplot2::aes(x = .data[["p1"]], y = .data[["p2"]])
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[["loglik"]]),
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(name = "log lik.") +
    ggplot2::geom_abline(
      intercept = 1,
      slope = -1,
      linetype = 2,
      colour = "grey40",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = true_df,
      ggplot2::aes(x = .data[["p1"]], y = .data[["p2"]]),
      inherit.aes = FALSE,
      colour = "white",
      fill = "#E41A1C",
      shape = 21,
      size = 3.2,
      stroke = 0.8
    ) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(2, 4, 2, 2),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = expression(p[1]),
      y = expression(p[2]),
      title = "Bulk log-likelihood"
    )
}

#' Log-likelihood profile in the ILR coordinate \eqn{\rho\in\mathbb{R}^{J-1}}
#'
#' For the bivariate toy (\eqn{J=2}) the free coordinate is scalar. The
#' profile evaluates [loglik_multivariate_constrained()] on a grid of
#' \eqn{\rho} (Helmert ILR). The MLE for \eqn{y=\mu p^{\star}} is marked at
#' [isometric_log_ratio()]\eqn{(p^{\star})}.
#'
#' @inheritParams plot_bulk_loglik_surface_p
#' @param grid Number of \eqn{\rho} evaluation points.
#' @param rho_lim Length-2 range for \eqn{\rho}. Default spans -4 to 10,
#'   expanded if needed to include the ILR image of
#'   \eqn{p^{\star}} (unbalanced compositions sit at large positive
#'   \eqn{\rho}).
#'
#' @return A `ggplot` (likelihood versus \eqn{\rho} on a log10 y-axis).
#' @export
plot_bulk_loglik_ilr_profile <- function(
  true_theta,
  grid = 400L,
  y = NULL,
  rho_lim = NULL
) {
  p_true <- as.numeric(true_theta$p)
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (length(p_true) != 2L) {
    stop(
      "plot_bulk_loglik_ilr_profile() is defined for J = 2.",
      call. = FALSE
    )
  }
  if (is.null(y)) {
    y <- drop(mu %*% p_true)
  }
  y <- as.numeric(y)
  rho_true <- as.numeric(isometric_log_ratio(p_true))
  if (is.null(rho_lim)) {
    rho_lim <- range(-4, 10, rho_true)
  }
  rho_seq <- seq(rho_lim[[1L]], rho_lim[[2L]], length.out = grid)
  ll <- vapply(
    rho_seq,
    function(rho) {
      tryCatch(
        loglik_multivariate_constrained(rho, y, mu, Sigma),
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )
  rho_mle <- rho_true
  ll_mle <- tryCatch(
    loglik_multivariate_constrained(rho_mle, y, mu, Sigma),
    error = function(e) NA_real_
  )
  if (!any(is.finite(ll))) {
    lik <- rep(NA_real_, length(ll))
    ll_max <- NA_real_
  } else {
    ll_max <- max(ll[is.finite(ll)], na.rm = TRUE)
    lik <- exp(pmin(ll - ll_max, 0))
    lik[!is.finite(ll)] <- NA_real_
    lik <- pmax(lik, 1e-16)
  }
  lik_mle <- if (is.finite(ll_mle) && is.finite(ll_max)) {
    pmax(exp(min(ll_mle - ll_max, 0)), 1e-16)
  } else {
    NA_real_
  }
  df <- data.frame(rho = rho_seq, loglik = ll, likelihood = lik)
  mle_df <- data.frame(
    rho = rho_mle,
    loglik = ll_mle,
    likelihood = lik_mle
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["rho"]], y = .data[["likelihood"]])
  ) +
    ggplot2::geom_line(colour = "#1B4F72", linewidth = 0.8) +
    ggplot2::geom_vline(
      xintercept = rho_mle,
      linetype = 2,
      colour = "grey40",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data = mle_df,
      ggplot2::aes(x = .data[["rho"]], y = .data[["likelihood"]]),
      inherit.aes = FALSE,
      colour = "white",
      fill = "#E41A1C",
      shape = 21,
      size = 3.2,
      stroke = 0.8
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::annotation_logticks(sides = "l") +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(2, 8, 2, 2),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(xlim = rho_lim, clip = "off") +
    ggplot2::labs(
      x = expression(rho),
      y = "Relative likelihood (log10 scale)",
      title = "Bulk log-likelihood (ILR)",
      caption = paste(
        "Y-axis is L / max(L) = exp(ell - max ell) on a log10 scale",
        "(ggplot2::annotation_logticks)."
      )
    )
}

#' Draw the proportion log-likelihood surface in an open rgl device
#'
#' @keywords internal
#' @noRd
.draw_proportion_loglik_rgl <- function(lat, lab) {
  z <- lat$z
  z[!is.finite(z)] <- min(z[is.finite(z)], na.rm = TRUE)
  rgl::persp3d(
    lat$p1,
    lat$p2,
    z,
    xlab = "",
    ylab = "",
    zlab = "",
    col = "steelblue",
    alpha = 0.85,
    polygon_offset = 1,
    axes = FALSE,
    box = TRUE,
    smooth = FALSE
  )
  rgl::aspect3d(1, 1, 0.65)
  rgl::axes3d(cex = 0.55, nticks = 4L)
  rgl::title3d(
    main = lab,
    xlab = "p1",
    ylab = "p2",
    zlab = "log lik.",
    cex = 0.7,
    font = 2L,
    line = 1
  )
  p_line <- seq(0.02, 0.98, length.out = 80L)
  z_line <- vapply(
    p_line,
    function(p1) {
      i <- which.min(abs(lat$p1 - p1))
      j <- which.min(abs(lat$p2 - (1 - p1)))
      z[i, j]
    },
    numeric(1)
  )
  rgl::lines3d(
    p_line,
    1 - p_line,
    z_line,
    col = "grey30",
    lwd = 2
  )
  if (length(lat$p_true) >= 2L && is.finite(lat$z_true)) {
    zr <- diff(range(z, na.rm = TRUE))
    rgl::spheres3d(
      lat$p_true[[1L]],
      lat$p_true[[2L]],
      lat$z_true,
      radius = 0.025,
      col = "#E41A1C"
    )
    rgl::texts3d(
      lat$p_true[[1L]],
      lat$p_true[[2L]],
      lat$z_true + 0.06 * max(zr, 1),
      texts = "MLE",
      col = "#222222",
      cex = 0.55,
      adj = c(0.5, 0)
    )
  }
  invisible(lat)
}

#' Snapshot the current rgl window as a ggplot raster
#'
#' @keywords internal
#' @noRd
.rgl_window_to_ggplot <- function(title) {
  if (!requireNamespace("png", quietly = TRUE)) {
    stop(
      "plot_bulk_loglik_rgl() needs the Suggests package png.",
      call. = FALSE
    )
  }
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  rgl::snapshot3d(tmp, fmt = "png", webshot = FALSE)
  if (!file.exists(tmp)) {
    stop("rgl::snapshot3d() did not write a PNG.", call. = FALSE)
  }
  img <- png::readPNG(tmp)
  ggplot2::ggplot() +
    ggplot2::annotation_raster(
      img,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1
    ) +
    ggplot2::coord_fixed(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA)
    ) +
    ggplot2::labs(title = title)
}

#' Open a high-resolution rgl window
#'
#' @keywords internal
#' @noRd
.rgl_open_hires <- function(width = 1600L, height = 1200L) {
  rgl::open3d()
  rgl::par3d(windowRect = c(40L, 40L, 40L + width, 40L + height))
  invisible(TRUE)
}

#' rgl surface of the bulk log-likelihood on a \eqn{(p_1,p_2)} lattice
#'
#' Opens an `rgl` window, draws [rgl::persp3d()] of
#' [loglik_multivariate()] versus hypothesised ratios, and marks the
#' MLE (true simulation proportions for \eqn{y=\mu p^{\star}}) with a
#' sphere. Returns a ggplot snapshot suitable for a PDF page.
#'
#' @inheritParams plot_bulk_loglik_surface_p
#' @param title Plot title (bold).
#'
#' @return A `ggplot` raster of the `rgl` snapshot.
#' @examplesIf interactive() && requireNamespace("rgl", quietly = TRUE) && requireNamespace("png", quietly = TRUE)
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' th <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
#' plot_bulk_loglik_rgl(th, grid = 20L)
#' @export
plot_bulk_loglik_rgl <- function(
  true_theta,
  grid = 40L,
  y = NULL,
  title = "Bulk log-likelihood"
) {
  .check_suggested_package("rgl", "plot_bulk_loglik_rgl")
  lat <- .proportion_loglik_lattice(true_theta, grid = grid, y = y)
  .rgl_open_hires()
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  .draw_proportion_loglik_rgl(lat, title)
  .rgl_window_to_ggplot(title)
}

#' Tile heatmap of a bivariate metric on the (\eqn{\rho_1},\eqn{\rho_2}) plane
#'
#' One page per meta-scenario (composition \eqn{\times} variance
#' structure \eqn{\times} CLD), with one panel per solver.
#'
#' @param metrics Long table with `algorithm`, correlations, and `value`.
#' @param title Page title.
#' @return A `ggplot`.
#' @export
plot_bivariate_metric_tiles <- function(metrics, title) {
  ggplot2::ggplot(
    metrics,
    ggplot2::aes(
      x = .data[["correlation_celltype1"]],
      y = .data[["correlation_celltype2"]],
      fill = .data[["value"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(~algorithm, ncol = 4L) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_viridis_c() +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      x = expression(rho[1]),
      y = expression(rho[2]),
      fill = NULL,
      title = title
    )
}

#' @keywords internal
#' @noRd
.corner_thetas <- function(config, theta_tbl, row) {
  corners <- .bivariate_corr_corners()
  out <- lapply(names(corners), function(lab) {
    rho <- corners[[lab]]
    hit <- config$proportions == row$proportions &
      config$variance == row$variance &
      config$centroids == row$centroids &
      abs(config$correlation_celltype1 - rho[[1L]]) < 1e-8 &
      abs(config$correlation_celltype2 - rho[[2L]]) < 1e-8
    if (!any(hit)) {
      return(NULL)
    }
    id <- config$ID[which(hit)[[1L]]]
    .unwrap_true_theta(
      theta_tbl$true_theta[theta_tbl$ID == id][[1L]]
    )
  })
  names(out) <- names(corners)
  out
}

#' @keywords internal
#' @noRd
.bivariate_page_meta <- function(config) {
  config <- .relevel_scenario_table(config)
  dplyr::distinct(
    config,
    .data[["centroids"]],
    .data[["variance"]],
    .data[["proportions"]]
  ) |>
    dplyr::arrange(
      .data[["centroids"]],
      .data[["variance"]],
      .data[["proportions"]]
    )
}

#' Page title matching performance-book factor order
#'
#' @keywords internal
#' @noRd
.page_title <- function(row) {
  cents <- as.character(row$centroids)
  var <- as.character(row$variance)
  prop <- as.character(row$proportions)
  paste(cents, var, prop, sep = " / ")
}

#' Persist ggplot `data` (not rgl snapshots) for a fig02 book
#'
#' @keywords internal
#' @noRd
.write_ggplot_rds <- function(object, data_rds, stem) {
  if (is.null(data_rds) || !nzchar(as.character(data_rds)[[1L]])) {
    return(invisible(NULL))
  }
  dir.create(data_rds, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(data_rds, paste0(stem, ".rds"))
  saveRDS(object, path)
  invisible(path)
}

#' Multi-page PDF of purified densities at four correlation corners
#'
#' @param config Slim config tibble with `ID`.
#' @param theta_tbl Tibble with `ID` and `true_theta`.
#' @param file Output PDF path.
#' @param n Draws per cell type.
#' @param data_rds Optional directory; when set, writes
#'   `purified_density.rds` (ggplot raster and overlay tables).
#' @return `file`, invisibly.
#' @export
save_bivariate_purified_density_book <- function(
  config,
  theta_tbl,
  file,
  n = 600L,
  data_rds = NULL
) {
  .save_bivariate_density_book(
    config,
    theta_tbl,
    file,
    n = n,
    which = "purified",
    data_rds = data_rds
  )
}

#' Multi-page PDF of bulk convolution densities
#'
#' @inheritParams save_bivariate_purified_density_book
#' @export
save_bivariate_bulk_density_book <- function(
  config,
  theta_tbl,
  file,
  n = 800L,
  data_rds = NULL
) {
  .save_bivariate_density_book(
    config,
    theta_tbl,
    file,
    n = n,
    which = "bulk",
    data_rds = data_rds
  )
}

#' Multi-page PDF of bulk log-likelihood surfaces on \eqn{(p_1,p_2)}
#'
#' @inheritParams save_bivariate_purified_density_book
#' @export
save_bivariate_loglik_surface_p_book <- function(
  config,
  theta_tbl,
  file,
  data_rds = NULL
) {
  .save_bivariate_loglik_ggplot_book(
    config,
    theta_tbl,
    file,
    which = "surface_p",
    data_rds = data_rds
  )
}

#' Multi-page PDF of ILR log-likelihood profiles (\eqn{\rho\in\mathbb{R}})
#'
#' @inheritParams save_bivariate_purified_density_book
#' @export
save_bivariate_loglik_ilr_profile_book <- function(
  config,
  theta_tbl,
  file,
  data_rds = NULL
) {
  .save_bivariate_loglik_ggplot_book(
    config,
    theta_tbl,
    file,
    which = "ilr_profile",
    data_rds = data_rds
  )
}

#' Multi-page PDF (and optional HTML) of rgl bulk log-likelihood surfaces
#'
#' Requires `rgl` and `png` for PDF snapshots. When `html_file` is set,
#' also writes an interactive [rgl::rglwidget()] HTML (needs
#' `htmlwidgets`).
#'
#' @inheritParams save_bivariate_purified_density_book
#' @param html_file Optional path for a self-contained interactive HTML.
#' @export
save_bivariate_loglik_rgl_book <- function(
  config,
  theta_tbl,
  file,
  html_file = NULL
) {
  .check_suggested_package("rgl", "save_bivariate_loglik_rgl_book")
  .check_suggested_package("png", "save_bivariate_loglik_rgl_book")
  meta <- .bivariate_page_meta(config)
  widgets <- list()
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    thetas <- .corner_thetas(config, theta_tbl, row)
    page_title <- .page_title(row)
    .rgl_open_hires(width = 1600L, height = 1200L)
    rgl::mfrow3d(2, 2, sharedMouse = TRUE)
    labs <- names(thetas)
    for (k in seq_along(thetas)) {
      if (k > 1L) {
        rgl::next3d()
      }
      th <- thetas[[k]]
      lab <- labs[[k]]
      if (is.null(th)) {
        next
      }
      lat <- .proportion_loglik_lattice(th, grid = 35L)
      .draw_proportion_loglik_rgl(lat, lab)
    }
    print(.rgl_window_to_ggplot(page_title))
    if (
      !is.null(html_file) && requireNamespace("htmlwidgets", quietly = TRUE)
    ) {
      widgets[[length(widgets) + 1L]] <- list(
        title = page_title,
        widget = rgl::rglwidget(width = 960, height = 720)
      )
    }
    try(rgl::close3d(), silent = TRUE)
  }
  if (!is.null(html_file) && length(widgets) > 0L) {
    .save_rgl_widget_book(widgets, html_file)
  } else if (!is.null(html_file)) {
    warning(
      "Interactive HTML skipped: install htmlwidgets.",
      call. = FALSE
    )
  }
  invisible(file)
}

#' Write one interactive HTML per rgl page plus an index
#'
#' @keywords internal
#' @noRd
.save_rgl_widget_book <- function(widgets, html_file) {
  .check_suggested_package("htmlwidgets", ".save_rgl_widget_book")
  html_file <- normalizePath(html_file, mustWork = FALSE)
  dir.create(dirname(html_file), recursive = TRUE, showWarnings = FALSE)
  pages_dir <- paste0(tools::file_path_sans_ext(html_file), "_pages")
  dir.create(pages_dir, recursive = TRUE, showWarnings = FALSE)
  links <- character(length(widgets))
  for (i in seq_along(widgets)) {
    page_path <- file.path(
      pages_dir,
      sprintf("page_%02d.html", i)
    )
    htmlwidgets::saveWidget(
      widgets[[i]]$widget,
      file = page_path,
      selfcontained = TRUE,
      title = widgets[[i]]$title
    )
    rel <- file.path(
      basename(pages_dir),
      basename(page_path)
    )
    links[[i]] <- sprintf(
      "<li><a href=\"%s\">%s</a></li>",
      rel,
      gsub("([&<>\"'])", "", widgets[[i]]$title)
    )
  }
  index_body <- paste0(
    "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
    "<title>Fig02 log-likelihood (rgl)</title></head><body>",
    "<h1>Fig02 bulk log-likelihood (interactive rgl)</h1>",
    "<ol>\n",
    paste(links, collapse = "\n"),
    "\n</ol></body></html>"
  )
  writeLines(index_body, html_file)
  invisible(html_file)
}

#' @keywords internal
#' @noRd
.save_bivariate_loglik_ggplot_book <- function(
  config,
  theta_tbl,
  file,
  which = c("surface_p", "ilr_profile"),
  data_rds = NULL
) {
  which <- match.arg(which)
  config <- .relevel_scenario_table(config)
  theta_tbl <- .relevel_scenario_table(theta_tbl)
  .check_suggested_package("gridExtra", "save_bivariate_loglik_ggplot_book")
  meta <- .bivariate_page_meta(config)
  collected <- list()
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    page_title <- .page_title(row)
    thetas <- .corner_thetas(config, theta_tbl, row)
    labs <- names(thetas)
    plots <- vector("list", length(labs))
    for (k in seq_along(labs)) {
      lab <- labs[[k]]
      th <- thetas[[lab]]
      if (is.null(th)) {
        plots[[k]] <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(lab)
        next
      }
      p <- switch(
        which,
        surface_p = plot_bulk_loglik_surface_p(th),
        ilr_profile = plot_bulk_loglik_ilr_profile(th)
      )
      d <- p$data
      if (is.data.frame(d) && nrow(d) > 0L) {
        d$page <- page_title
        d$panel <- lab
        collected[[length(collected) + 1L]] <- d
      }
      plots[[k]] <- p +
        ggplot2::ggtitle(lab) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold")
        )
    }
    grob <- gridExtra::arrangeGrob(
      grobs = plots,
      ncol = 2L,
      top = grid::textGrob(
        page_title,
        gp = grid::gpar(fontface = "bold", fontsize = 14)
      )
    )
    grid::grid.draw(grob)
    if (i < nrow(meta)) {
      grid::grid.newpage()
    }
  }
  stem <- if (identical(which, "ilr_profile")) {
    "loglik_ilr_profile"
  } else {
    "loglik_surface_p"
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, stem)
  }
  invisible(file)
}

#' @keywords internal
#' @noRd
.save_bivariate_density_book <- function(
  config,
  theta_tbl,
  file,
  n,
  which,
  data_rds = NULL
) {
  config <- .relevel_scenario_table(config)
  theta_tbl <- .relevel_scenario_table(theta_tbl)
  meta <- .bivariate_page_meta(config)
  collected <- list()
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    thetas <- .corner_thetas(config, theta_tbl, row)
    p <- .plot_gene_density_page(
      thetas,
      n = n,
      which = which,
      title = .page_title(row),
      share_legend = identical(which, "bulk")
    )
    gd <- attr(p, "ggplot_data")
    if (!is.null(gd)) {
      collected[[i]] <- list(
        page = .page_title(row),
        raster = gd$raster,
        centroids = gd$centroids,
        ellipses = gd$ellipses
      )
    }
    print(p)
  }
  stem <- if (identical(which, "purified")) {
    "purified_density"
  } else {
    "bulk_density"
  }
  .write_ggplot_rds(collected, data_rds, stem)
  invisible(file)
}

#' Write RMSE / MAE / Aitchison tile PDFs for the bivariate toy
#'
#' @param artefacts List from `read_simulation_artefacts()` with
#'   `assemble = TRUE`, or the same named pieces.
#' @param dir Output directory.
#' @param data_rds Optional directory for ggplot `data` RDS files.
#' @return Named paths.
#' @export
save_bivariate_metric_heatmaps <- function(artefacts, dir, data_rds = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  agg <- .aggregate_bivariate_metrics(artefacts)

  write_metric_pdf <- function(metric, path, stem) {
    meta <- .bivariate_page_meta(agg)
    pages <- list()
    grDevices::pdf(path, width = 16, height = 10)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (i in seq_len(nrow(meta))) {
      row <- meta[i, , drop = FALSE]
      page <- agg[
        agg$proportions == row$proportions &
          agg$variance == row$variance &
          agg$centroids == row$centroids,
        ,
        drop = FALSE
      ]
      page$value <- page[[metric]]
      page$page <- .page_title(row)
      title <- paste(
        toupper(metric),
        row$proportions,
        row$variance,
        row$centroids,
        sep = " / "
      )
      print(plot_bivariate_metric_tiles(page, title))
      pages[[i]] <- page
    }
    .write_ggplot_rds(dplyr::bind_rows(pages), data_rds, stem)
    path
  }

  list(
    rmse = write_metric_pdf(
      "rmse",
      file.path(dir, "heatmap_rmse.pdf"),
      "heatmap_rmse"
    ),
    mae = write_metric_pdf(
      "mae",
      file.path(dir, "heatmap_mae.pdf"),
      "heatmap_mae"
    ),
    aitchison = write_metric_pdf(
      "aitchison",
      file.path(dir, "heatmap_aitchison.pdf"),
      "heatmap_aitchison"
    )
  )
}

#' @keywords internal
#' @noRd
.aggregate_bivariate_metrics <- function(artefacts) {
  cfg <- artefacts$config
  opt <- artefacts$optimisation
  if (is.null(opt) || !"ID" %in% names(opt)) {
    stop("`optimisation` with an `ID` column is required.", call. = FALSE)
  }
  if (!is.null(artefacts$theta) && "true_theta" %in% names(artefacts$theta)) {
    theta_map <- artefacts$theta
  } else {
    theta_map <- tibble::tibble(
      ID = cfg$ID,
      true_theta = artefacts$theta_true
    )
  }
  p_true <- purrr::pmap_dfr(
    theta_map,
    function(ID, true_theta, ...) {
      th <- .unwrap_true_theta(true_theta)
      tibble::tibble(
        ID = ID,
        p1 = as.numeric(th$p)[[1L]],
        p2 = as.numeric(th$p)[[2L]]
      )
    }
  )
  scored <- dplyr::left_join(opt, p_true, by = "ID")
  scored$mae <- 0.5 *
    (abs(scored$celltype_1 - scored$p1) + abs(scored$celltype_2 - scored$p2))
  scored$rmse_row <- sqrt(
    0.5 *
      ((scored$celltype_1 - scored$p1)^2 + (scored$celltype_2 - scored$p2)^2)
  )
  scored$aitchison <- .aitchison_pair(
    scored$p1,
    scored$p2,
    scored$celltype_1,
    scored$celltype_2
  )
  scored |>
    dplyr::group_by(.data[["ID"]], .data[["algorithm"]]) |>
    dplyr::summarise(
      rmse = mean(.data[["rmse_row"]], na.rm = TRUE),
      mae = mean(.data[["mae"]], na.rm = TRUE),
      aitchison = mean(.data[["aitchison"]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(cfg, by = "ID")
}

#' Scenario IDs at the four correlation corners of one meta-row
#'
#' @keywords internal
#' @noRd
.corner_ids <- function(config, row) {
  corners <- .bivariate_corr_corners()
  vapply(
    names(corners),
    function(lab) {
      rho <- corners[[lab]]
      hit <- config$proportions == row$proportions &
        config$variance == row$variance &
        config$centroids == row$centroids &
        abs(config$correlation_celltype1 - rho[[1L]]) < 1e-8 &
        abs(config$correlation_celltype2 - rho[[2L]]) < 1e-8
      if (!any(hit)) {
        return(NA_character_)
      }
      as.character(config$ID[which(hit)[[1L]]])
    },
    character(1)
  )
}

#' Restrict a benchmark list to selected scenario IDs
#'
#' @keywords internal
#' @noRd
.subset_benchmark_ids <- function(benchmark, ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  cfg <- benchmark$config
  cfg <- cfg[as.character(cfg$ID) %in% ids, , drop = FALSE]
  if (!is.null(benchmark$theta) && "ID" %in% names(benchmark$theta)) {
    th <- benchmark$theta
    th <- th[match(as.character(cfg$ID), as.character(th$ID)), , drop = FALSE]
    benchmark$theta <- th
    benchmark$theta_true <- lapply(th$true_theta, .unwrap_true_theta)
  } else if (
    !is.null(benchmark$theta_true) &&
      length(benchmark$theta_true) == nrow(benchmark$config)
  ) {
    keep <- as.character(benchmark$config$ID) %in% ids
    benchmark$theta_true <- benchmark$theta_true[keep]
  }
  benchmark$config <- cfg
  for (nm in c("monte_carlo", "optimisation")) {
    tbl <- benchmark[[nm]]
    if (!is.null(tbl) && "ID" %in% names(tbl)) {
      benchmark[[nm]] <- tbl[as.character(tbl$ID) %in% ids, , drop = FALSE]
    }
  }
  benchmark
}

#' @keywords internal
#' @noRd
.algorithm_in <- function(x, keep) {
  kn <- unique(vapply(keep, .normalise_algorithm_key, character(1)))
  .normalise_algorithm_key(x) %in% kn
}

#' @keywords internal
#' @noRd
.bivariate_wald_algorithms <- function() {
  c("gradient", "LBFGS", "Newton-Raphson", "Marquardt-Levenberg")
}

#' Wald forest of mean estimates at four correlation corners
#'
#' @keywords internal
#' @noRd
.plot_bivariate_wald_forest_page <- function(plot_df, title) {
  plot_df$algorithm <- .relevel_algorithm(plot_df$algorithm)
  plot_df$panel <- factor(
    plot_df$panel,
    levels = names(.bivariate_corr_corners())
  )
  plot_df$cell_type <- factor(
    plot_df$cell_type,
    levels = unique(as.character(plot_df$cell_type))
  )
  pal <- c("#E41A1C", "#4DAF4A")
  ct <- levels(plot_df$cell_type)
  names(pal) <- ct[seq_len(min(2L, length(ct)))]
  alg_lvls <- levels(droplevels(plot_df$algorithm))
  n_ct <- nlevels(plot_df$cell_type)
  dodge_width <- 0.55
  dodge_step <- dodge_width / max(n_ct, 1L)
  ct_idx <- as.numeric(plot_df$cell_type)
  offset <- (ct_idx - (n_ct + 1) / 2) * dodge_step
  plot_df$y_dodge <- as.numeric(plot_df$algorithm) + offset
  plot_df$y_arrow <- plot_df$y_dodge + 0.22
  truth <- dplyr::distinct(
    plot_df,
    .data[["panel"]],
    .data[["cell_type"]],
    .data[["p_true"]]
  )
  emp_df <- dplyr::filter(
    plot_df,
    is.finite(.data[["emp_lo"]]) & is.finite(.data[["emp_hi"]])
  )
  th_df <- dplyr::filter(
    plot_df,
    is.finite(.data[["ci_lo"]]) & is.finite(.data[["ci_hi"]])
  )
  arrow_df <- dplyr::filter(
    plot_df,
    is.finite(.data[["mean_est"]]) &
      is.finite(.data[["p_true"]]) &
      abs(.data[["mean_est"]] - .data[["p_true"]]) > 0.008
  )
  annot_df <- plot_df |>
    dplyr::group_by(.data[["panel"]], .data[["algorithm"]]) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(y_lab = as.numeric(.data[["algorithm"]]))
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["mean_est"]],
      y = .data[["y_dodge"]],
      colour = .data[["cell_type"]]
    )
  ) +
    ggplot2::geom_vline(
      data = truth,
      ggplot2::aes(
        xintercept = .data[["p_true"]],
        colour = .data[["cell_type"]]
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      data = emp_df,
      ggplot2::aes(xmin = .data[["emp_lo"]], xmax = .data[["emp_hi"]]),
      width = 0.12,
      linewidth = 0.45,
      linetype = "22",
      orientation = "y"
    ) +
    ggplot2::geom_errorbar(
      data = th_df,
      ggplot2::aes(xmin = .data[["ci_lo"]], xmax = .data[["ci_hi"]]),
      width = 0.22,
      linewidth = 0.75,
      orientation = "y"
    ) +
    ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(
        x = .data[["mean_est"]],
        xend = .data[["p_true"]],
        y = .data[["y_arrow"]],
        yend = .data[["y_arrow"]],
        colour = .data[["cell_type"]]
      ),
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed",
        angle = 25
      ),
      linewidth = 0.85,
      lineend = "butt",
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_label(
      data = annot_df,
      ggplot2::aes(
        x = Inf,
        y = .data[["y_lab"]],
        label = .data[["annot"]]
      ),
      inherit.aes = FALSE,
      hjust = 1.08,
      vjust = 0.5,
      size = 3.4,
      fontface = "bold",
      lineheight = 0.95,
      label.size = 0.25,
      label.padding = grid::unit(0.28, "lines"),
      fill = "#F4F1EA",
      colour = "grey20",
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(alg_lvls),
      labels = alg_lvls,
      expand = ggplot2::expansion(add = 0.45)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.04, 0.38))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~panel, ncol = 2L) +
    ggplot2::labs(
      x = "Mean Monte Carlo estimate",
      y = NULL,
      colour = "Cell type",
      title = title,
      caption = paste(
        "Solid whiskers: mean estimate plus or minus 1.96 times the",
        "expected-Fisher Wald SE at the true composition",
        "(confint.decovart_fit / vcov_ilr_delta).",
        "Dashed whiskers: plus or minus 1.96 times the empirical",
        "Monte Carlo SD.",
        "Horizontal arrows: bias toward the true proportion.",
        "Bold labels (right, one per solver): RMSE and coverage of",
        "those Wald intervals (identical for p1 and p2 on the unit",
        "simplex). NNLS, LSEI, and SA omitted (no convolution Wald SE)."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(4, 10, 4, 4)
    )
}

#' Build the Wald-forest table for one meta-scenario (four corners)
#'
#' @keywords internal
#' @noRd
.bivariate_wald_forest_table <- function(artefacts, ids, labels) {
  z <- stats::qnorm(0.975)
  pieces <- lapply(seq_along(ids), function(k) {
    id <- ids[[k]]
    if (is.na(id) || !nzchar(id)) {
      return(NULL)
    }
    sub <- .subset_benchmark_ids(artefacts, id)
    if (
      is.null(sub$monte_carlo) ||
        !"theoretical_se" %in% names(sub$monte_carlo) ||
        !any(is.finite(sub$monte_carlo$theoretical_se))
    ) {
      sub <- .attach_expected_fisher_wald(sub)
    }
    mc <- sub$monte_carlo
    if (is.null(mc) || nrow(mc) == 0L) {
      return(NULL)
    }
    keep <- .algorithm_in(mc$algorithm, .bivariate_wald_algorithms())
    mc <- mc[keep, , drop = FALSE]
    if (nrow(mc) == 0L) {
      return(NULL)
    }
    long <- pivot_mc_estimates(sub)
    truth <- dplyr::distinct(
      long,
      .data[["algorithm"]],
      .data[["cell_type"]],
      .data[["p_true"]]
    )
    out <- dplyr::left_join(
      mc,
      truth,
      by = intersect(c("algorithm", "cell_type", "ID"), names(truth))
    )
    if (!"p_true" %in% names(out)) {
      return(NULL)
    }
    se_th <- if ("theoretical_se" %in% names(out)) {
      out$theoretical_se
    } else {
      rep(NA_real_, nrow(out))
    }
    se_th <- ifelse(is.finite(se_th), se_th, NA_real_)
    out$mean_est <- out$p_true + out$bias
    out$ci_lo <- out$mean_est - z * se_th
    out$ci_hi <- out$mean_est + z * se_th
    out$emp_lo <- out$mean_est - z * out$empirical_sd
    out$emp_hi <- out$mean_est + z * out$empirical_sd
    cov_pct <- ifelse(
      is.finite(out$coverage),
      sprintf("%.0f%%", 100 * out$coverage),
      "NA"
    )
    rmse_lab <- ifelse(
      is.finite(out$rmse),
      sprintf("RMSE %.3f", out$rmse),
      "RMSE NA"
    )
    out$annot <- paste(rmse_lab, cov_pct, sep = "\n")
    out$panel <- labels[[k]]
    out <- out[is.finite(out$mean_est), , drop = FALSE]
    out
  })
  dplyr::bind_rows(pieces)
}

#' Monte Carlo metric dots on the 9-by-9 correlation grid
#'
#' Enrichplot-style discrete grid: colour is RMSE and size is
#' Aitchison distance (means over Monte Carlo replicates).
#'
#' @keywords internal
#' @noRd
.plot_bivariate_solver_dots_page <- function(df, title) {
  df$algorithm <- .relevel_algorithm(df$algorithm)
  rho1 <- sort(unique(df$correlation_celltype1))
  rho2 <- sort(unique(df$correlation_celltype2))
  df$rho1 <- factor(
    df$correlation_celltype1,
    levels = rho1,
    labels = format(rho1, digits = 2, trim = TRUE)
  )
  df$rho2 <- factor(
    df$correlation_celltype2,
    levels = rho2,
    labels = format(rho2, digits = 2, trim = TRUE)
  )
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["rho1"]],
      y = .data[["rho2"]],
      colour = .data[["rmse"]],
      size = .data[["aitchison"]]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~algorithm, ncol = 4L) +
    ggplot2::scale_colour_viridis_c(name = "RMSE") +
    ggplot2::scale_size_continuous(
      name = "Aitchison",
      range = c(1.5, 10)
    ) +
    ggplot2::labs(
      x = expression(rho[1]),
      y = expression(rho[2]),
      title = title,
      caption = paste(
        "Each panel is the 9-by-9 correlation factorial (81 cells).",
        "Colour: mean RMSE; size: mean Aitchison distance."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    )
}

#' Multi-page PDF of clustered algorithm-similarity heatmaps
#'
#' Twelve pages (CLD \eqn{\times} variance \eqn{\times} composition). Each
#' page is a 2-by-2 of the correlation corners, with average-linkage
#' clustering of \(1-r\) and a dendrogram to the right of the tiles,
#' with leaves flush against the heatmap.
#'
#' @inheritParams save_bivariate_metric_heatmaps
#' @param file Output PDF path.
#' @param data_rds Optional directory for ggplot `data` RDS files.
#' @return `file`, invisibly.
#' @export
save_bivariate_similarity_book <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  .check_suggested_package("gridExtra", "save_bivariate_similarity_book")
  .check_plot_dependencies(need_ggdendro = TRUE, need_cowplot = TRUE)
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  grDevices::pdf(file, width = 16, height = 14)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    ids <- .corner_ids(cfg, row)
    labs <- names(ids)
    page_title <- .page_title(row)
    panels <- vector("list", length(ids))
    for (k in seq_along(ids)) {
      if (is.na(ids[[k]])) {
        panels[[k]] <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::ggtitle(labs[[k]])
        next
      }
      sub <- .subset_benchmark_ids(artefacts, ids[[k]])
      sim <- algorithm_similarity(sub)
      sim$page <- page_title
      sim$panel <- labs[[k]]
      collected[[length(collected) + 1L]] <- sim
      p <- plot_algorithm_similarity(sub, dendrogram = TRUE)
      panels[[k]] <- cowplot::plot_grid(
        cowplot::ggdraw() +
          cowplot::draw_label(labs[[k]], fontface = "bold", size = 11),
        p,
        ncol = 1,
        rel_heights = c(0.08, 1)
      )
    }
    grob <- gridExtra::arrangeGrob(
      grobs = panels,
      ncol = 2L,
      top = grid::textGrob(
        page_title,
        gp = grid::gpar(fontface = "bold", fontsize = 14)
      )
    )
    grid::grid.draw(grob)
    if (i < nrow(meta)) {
      grid::grid.newpage()
    }
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "similarity")
  }
  invisible(file)
}

#' Multi-page PDF of Wald forests at four correlation corners
#'
#' @inheritParams save_bivariate_similarity_book
#' @export
save_bivariate_forest_book <- function(artefacts, file, data_rds = NULL) {
  .check_suggested_package("gridExtra", "save_bivariate_forest_book")
  artefacts <- .attach_expected_fisher_wald(artefacts)
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  grDevices::pdf(file, width = 18, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    ids <- .corner_ids(cfg, row)
    tbl <- .bivariate_wald_forest_table(artefacts, ids, names(ids))
    if (is.null(tbl) || nrow(tbl) == 0L) {
      next
    }
    tbl$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- tbl
    print(.plot_bivariate_wald_forest_page(tbl, .page_title(row)))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "forest")
  }
  invisible(file)
}

#' Multi-page PDF of rainclouds at four correlation corners
#'
#' @inheritParams save_bivariate_similarity_book
#' @export
save_bivariate_raincloud_book <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  .check_plot_dependencies(need_ggdist = TRUE)
  cfg <- artefacts$config
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  grDevices::pdf(file, width = 16, height = 30)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    ids <- .corner_ids(cfg, row)
    parts <- lapply(seq_along(ids), function(k) {
      if (is.na(ids[[k]])) {
        return(NULL)
      }
      sub <- .subset_benchmark_ids(artefacts, ids[[k]])
      long <- pivot_mc_estimates(sub)
      long$panel <- names(ids)[[k]]
      long
    })
    df <- dplyr::bind_rows(parts)
    if (nrow(df) == 0L) {
      next
    }
    df$panel <- factor(df$panel, levels = names(ids))
    df$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- df
    p <- plot_mc_raincloud(
      df,
      quantity = "estimate",
      include_dots = FALSE,
      dodge_width = 1,
      slab_scale = 1.05,
      slab_alpha = 0.45
    )
    p <- p +
      ggplot2::facet_wrap(~panel, ncol = 2L) +
      ggplot2::ggtitle(.page_title(row)) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
    print(p)
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "raincloud")
  }
  invisible(file)
}

#' 12-page PDF of RMSE/Aitchison solver dots on the correlation grid
#'
#' One page per meta-scenario. Each page facets solvers on the same
#' 9-by-9 \eqn{(\rho_1,\rho_2)} factorial as the metric heatmaps
#' (enrichplot-style dots: colour = RMSE, size = Aitchison).
#'
#' @inheritParams save_bivariate_similarity_book
#' @export
save_bivariate_solver_dots_book <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  agg <- .aggregate_bivariate_metrics(artefacts)
  meta <- .bivariate_page_meta(agg)
  collected <- list()
  grDevices::pdf(file, width = 16, height = 10)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    page <- agg[
      agg$proportions == row$proportions &
        agg$variance == row$variance &
        agg$centroids == row$centroids,
      ,
      drop = FALSE
    ]
    if (nrow(page) == 0L) {
      next
    }
    page$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- page
    print(.plot_bivariate_solver_dots_page(page, .page_title(row)))
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(dplyr::bind_rows(collected), data_rds, "solver_dots")
  }
  invisible(file)
}
