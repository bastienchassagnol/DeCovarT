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
    "rho = (-0.8, 0.8)" = c(-0.8, 0.8),
    "rho = (0.8, 0.8)" = c(0.8, 0.8)
  )
}

#' Exact Gaussian probability ellipse (known mean and covariance)
#'
#' Boundary of the set
#' \eqn{(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^2_{2,\alpha}}
#' for a bivariate Gaussian with **known** \(\mu\) and \(\Sigma\).
#' The squared Mahalanobis distance is exactly \(\chi^2_2\).
#'
#' @param mu Length-2 mean.
#' @param Sigma \(2\times 2\) covariance.
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
.theme_density_2d <- function() {
  theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.background = ggplot2::element_rect(
        fill = "grey12",
        colour = NA
      ),
      panel.grid.major = ggplot2::element_line(
        colour = "grey28",
        linewidth = 0.2
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
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
.add_celltype_overlay_layers <- function(p, overlay, cts) {
  pal <- .celltype_colour_values(cts)
  shp <- .celltype_shape_values(cts)
  p +
    ggplot2::geom_path(
      data = overlay$ellipses,
      ggplot2::aes(
        x = .data[["gene_1"]],
        y = .data[["gene_2"]],
        colour = .data[["cell_type"]],
        group = .data[["cell_type"]]
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
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

#' 2-D density of purified Gaussians
#'
#' Filled contours use `contour_var = "ndensity"` so peak intensity is
#' comparable across facets. Cell-type 1 is a red circle; cell-type 2
#' is a green triangle. Outlines are exact 95% Gaussian ellipses for
#' the known means and covariances ([gaussian_confidence_ellipse()]).
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
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_density_2d_filled(
      contour_var = "ndensity",
      bins = 12L,
      colour = NA
    ) +
    ggplot2::scale_fill_viridis_d(name = "Density") +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = "Purified Gaussians"
    )
  .add_celltype_overlay_layers(p, overlay, cts)
}

#' 2-D density of the bulk convolution \eqn{y\sim N(\mu p,\Sigma(p))}
#'
#' @inheritParams plot_purified_density_2d
#' @return A `ggplot`.
#' @export
plot_bulk_convolution_density_2d <- function(true_theta, n = 1200L) {
  cts <- .celltype_names(true_theta)
  df <- .bulk_density_draws(true_theta, n)
  overlay <- .celltype_overlay(true_theta)
  lims <- .gene_axis_limits(list(true_theta))
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_density_2d_filled(
      contour_var = "ndensity",
      bins = 12L,
      colour = NA
    ) +
    ggplot2::scale_fill_viridis_d(name = "Density") +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = "Bulk convolution"
    )
  .add_celltype_overlay_layers(p, overlay, cts)
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
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[["gene_1"]], y = .data[["gene_2"]])
  ) +
    ggplot2::geom_density_2d_filled(
      contour_var = "ndensity",
      bins = 12L,
      colour = NA
    ) +
    ggplot2::scale_fill_viridis_d(name = "Density") +
    ggplot2::facet_wrap(~panel, ncol = 2L) +
    ggplot2::coord_equal(
      xlim = lims$xlim,
      ylim = lims$ylim,
      expand = FALSE
    ) +
    .theme_density_2d() +
    ggplot2::labs(
      x = "Gene 1",
      y = "Gene 2",
      title = title
    )
  p <- .add_celltype_overlay_layers(p, overlay, cts)
  if (!isTRUE(share_legend)) {
    return(p)
  }
  .attach_cowplot_legend(p)
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

#' Contour of the bulk log-likelihood on a \((p_1,p_2)\) lattice
#'
#' Axes are hypothesised cell-type ratios, not gene expression. The
#' true simulation proportions are marked; the dashed line is the
#' simplex \(p_1+p_2=1\).
#'
#' @inheritParams plot_purified_density_2d
#' @param grid Length of the lattice per axis.
#' @param y Optional bulk observation. Default is \(\mu p^{\star}\).
#'
#' @return A `ggplot`.
#' @export
plot_bulk_loglik_surface <- function(
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
      colour = "grey80",
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
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = expression(p[1]),
      y = expression(p[2]),
      title = "Bulk log-likelihood"
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
    xlab = "p1",
    ylab = "p2",
    zlab = "log lik.",
    col = "steelblue",
    alpha = 0.85,
    polygon_offset = 1
  )
  rgl::aspect3d(1, 1, 0.65)
  rgl::title3d(main = lab, cex = 1.15, font = 2L)
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
    col = "grey85",
    lwd = 2
  )
  if (length(lat$p_true) >= 2L && is.finite(lat$z_true)) {
    zr <- diff(range(z, na.rm = TRUE))
    rgl::spheres3d(
      lat$p_true[[1L]],
      lat$p_true[[2L]],
      lat$z_true,
      radius = 0.03,
      col = "#E41A1C"
    )
    rgl::texts3d(
      lat$p_true[[1L]],
      lat$p_true[[2L]],
      lat$z_true + 0.08 * max(zr, 1),
      texts = "true p",
      col = "white",
      cex = 0.8,
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
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = title)
}

#' rgl surface of the bulk log-likelihood on a \((p_1,p_2)\) lattice
#'
#' Opens an `rgl` window, draws [rgl::persp3d()] of
#' [loglik_multivariate()] versus hypothesised ratios, and marks the
#' true simulation proportions with a sphere (the MLE of \(y=\mu
#' p^{\star}\) under the convolution model). Returns a ggplot snapshot
#' suitable for a PDF page.
#'
#' @inheritParams plot_bulk_loglik_surface
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
  rgl::open3d()
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
  dplyr::distinct(
    config,
    .data[["proportions"]],
    .data[["variance"]],
    .data[["centroids"]]
  )
}

#' @keywords internal
#' @noRd
.page_title <- function(row) {
  paste(row$proportions, row$variance, row$centroids, sep = " / ")
}

#' Multi-page PDF of purified densities at four correlation corners
#'
#' @param config Slim config tibble with `ID`.
#' @param theta_tbl Tibble with `ID` and `true_theta`.
#' @param file Output PDF path.
#' @param n Draws per cell type.
#' @return `file`, invisibly.
#' @export
save_bivariate_purified_density_book <- function(
  config,
  theta_tbl,
  file,
  n = 600L
) {
  .save_bivariate_density_book(
    config,
    theta_tbl,
    file,
    n = n,
    which = "purified"
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
  n = 800L
) {
  .save_bivariate_density_book(
    config,
    theta_tbl,
    file,
    n = n,
    which = "bulk"
  )
}

#' Multi-page PDF of bulk log-likelihood surfaces on \((p_1,p_2)\)
#'
#' @inheritParams save_bivariate_purified_density_book
#' @export
save_bivariate_loglik_surface_book <- function(
  config,
  theta_tbl,
  file
) {
  .save_bivariate_loglik_ggplot_book(config, theta_tbl, file)
}

#' Multi-page PDF of rgl bulk log-likelihood surfaces on \((p_1,p_2)\)
#'
#' Requires `rgl` (and `png` to embed snapshots). Does **not** fall
#' back to the ggplot contour book.
#'
#' @inheritParams save_bivariate_purified_density_book
#' @export
save_bivariate_loglik_rgl_book <- function(
  config,
  theta_tbl,
  file
) {
  .check_suggested_package("rgl", "save_bivariate_loglik_rgl_book")
  .check_suggested_package("png", "save_bivariate_loglik_rgl_book")
  meta <- .bivariate_page_meta(config)
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    thetas <- .corner_thetas(config, theta_tbl, row)
    page_title <- .page_title(row)
    rgl::open3d()
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
    try(rgl::close3d(), silent = TRUE)
  }
  invisible(file)
}

#' @keywords internal
#' @noRd
.save_bivariate_loglik_ggplot_book <- function(config, theta_tbl, file) {
  .check_suggested_package("gridExtra", "save_bivariate_loglik_surface_book")
  meta <- .bivariate_page_meta(config)
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    thetas <- .corner_thetas(config, theta_tbl, row)
    plots <- lapply(names(thetas), function(lab) {
      th <- thetas[[lab]]
      if (is.null(th)) {
        return(
          ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle(lab)
        )
      }
      plot_bulk_loglik_surface(th) +
        ggplot2::ggtitle(lab) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold")
        )
    })
    grob <- gridExtra::arrangeGrob(
      grobs = plots,
      ncol = 2L,
      top = grid::textGrob(
        .page_title(row),
        gp = grid::gpar(fontface = "bold", fontsize = 14)
      )
    )
    grid::grid.draw(grob)
    if (i < nrow(meta)) {
      grid::grid.newpage()
    }
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
  which
) {
  meta <- .bivariate_page_meta(config)
  grDevices::pdf(file, width = 14, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    thetas <- .corner_thetas(config, theta_tbl, row)
    print(
      .plot_gene_density_page(
        thetas,
        n = n,
        which = which,
        title = .page_title(row),
        share_legend = identical(which, "bulk")
      )
    )
  }
  invisible(file)
}

#' Write RMSE / MAE / Aitchison tile PDFs for the bivariate toy
#'
#' @param artefacts List from `read_simulation_artefacts()` with
#'   `assemble = TRUE`, or the same named pieces.
#' @param dir Output directory.
#' @return Named paths.
#' @export
save_bivariate_metric_heatmaps <- function(artefacts, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
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
  agg <- scored |>
    dplyr::group_by(.data[["ID"]], .data[["algorithm"]]) |>
    dplyr::summarise(
      rmse = mean(.data[["rmse_row"]], na.rm = TRUE),
      mae = mean(.data[["mae"]], na.rm = TRUE),
      aitchison = mean(.data[["aitchison"]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(cfg, by = "ID")

  write_metric_pdf <- function(metric, path) {
    meta <- dplyr::distinct(
      agg,
      .data[["proportions"]],
      .data[["variance"]],
      .data[["centroids"]]
    )
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
      title <- paste(
        toupper(metric),
        row$proportions,
        row$variance,
        row$centroids,
        sep = " / "
      )
      print(plot_bivariate_metric_tiles(page, title))
    }
    path
  }

  list(
    rmse = write_metric_pdf(
      "rmse",
      file.path(dir, "heatmap_rmse.pdf")
    ),
    mae = write_metric_pdf(
      "mae",
      file.path(dir, "heatmap_mae.pdf")
    ),
    aitchison = write_metric_pdf(
      "aitchison",
      file.path(dir, "heatmap_aitchison.pdf")
    )
  )
}
