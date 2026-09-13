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

#' Fig03 2-by-2 of CT1 / CT2 graph families (CT3 held at Scale-free)
#'
#' Labels match `.relevel_scenario_table()` display names.
#'
#' @keywords internal
#' @noRd
.hybrid_topology_panels <- function() {
  sf <- unname(.graph_display_labels()[["scale_free"]])
  sb <- unname(.graph_display_labels()[["stochastic_block_model"]])
  list(
    "CT1 Scale-free / CT2 Scale-free" = c(sf, sf),
    "CT1 Scale-free / CT2 Cluster SBM" = c(sf, sb),
    "CT1 Cluster SBM / CT2 Scale-free" = c(sb, sf),
    "CT1 Cluster SBM / CT2 Cluster SBM" = c(sb, sb)
  )
}

#' Covariance-driven (fig03) factorial: graph families plus MixSim overlap
#'
#' @keywords internal
#' @noRd
.is_hybrid_config <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  nms <- names(x)
  all(c("graph_ct1", "graph_ct2", "overlap_label") %in% nms)
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
  p <- as.numeric(true_theta$p)
  p <- p / sum(p)
  n_j <- as.integer(pmax(1L, round(n * p)))
  parts <- lapply(seq_along(cts), function(j) {
    .rnorm_2d(n_j[[j]], mu[, j], Sigma[,, j], cts[[j]])
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

#' Multivariate Gaussian density on rows of a two-column matrix
#'
#' @keywords internal
#' @noRd
.dmvnorm_rows <- function(xy, mu, sigma) {
  xy <- as.matrix(xy)
  mu <- as.numeric(mu)
  sigma <- as.matrix(sigma)
  k <- ncol(sigma)
  z <- sweep(xy, 2L, mu, "-")
  qf <- rowSums(z * t(solve(sigma, t(z))))
  log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  exp(-0.5 * (k * log(2 * pi) + log_det + qf))
}

#' Exact p-weighted mixture of purified Gaussians on a gene-axis grid
#'
#' Density \eqn{\sum_j p_j\,\varphi(x;\mu_j,\Sigma_j)}, not the
#' convolution \eqn{\mathcal{N}(\mu p,\Sigma(p))}.
#'
#' @keywords internal
#' @noRd
.purified_mixture_raster <- function(
  true_theta,
  xlim,
  ylim,
  n = 80L,
  panel = NULL
) {
  p <- as.numeric(true_theta$p)
  p <- p / sum(p)
  mu <- true_theta$mu
  sigma <- true_theta$sigma
  xs <- seq(xlim[[1L]], xlim[[2L]], length.out = n)
  ys <- seq(ylim[[1L]], ylim[[2L]], length.out = n)
  grid <- expand.grid(
    gene_1 = xs,
    gene_2 = ys,
    KEEP.OUT.ATTRS = FALSE
  )
  dens <- numeric(nrow(grid))
  xy <- cbind(grid$gene_1, grid$gene_2)
  for (j in seq_len(ncol(mu))) {
    dens <- dens + p[[j]] * .dmvnorm_rows(xy, mu[, j], sigma[,, j])
  }
  grid$density <- dens
  if (!is.null(panel)) {
    grid$panel <- panel
  }
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
#' An exact mixture raster
#' \eqn{\sum_j p_j\,\varphi(x;\boldsymbol{\mu}_j,\boldsymbol{\Sigma}_j)}
#' fills the shared gene-axis window, so rare types contribute mass in
#' proportion to the scenario's composition (not an equal-sized kernel
#' density per cell type). This is **not** the bulk convolution
#' \eqn{\mathcal{N}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))}.
#' Cell-type 1 is a red circle; cell-type 2 is a green triangle.
#' Outlines are exact 95% Gaussian ellipses for the known means and
#' covariances ([gaussian_confidence_ellipse()]):
#' \eqn{(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^{2}_{2,0.95}}.
#'
#' @param true_theta List with `p`, `mu`, `sigma` for \eqn{G=2}.
#' @param n Unused for the raster (kept for API compatibility with
#'   [plot_bulk_convolution_density_2d()]); mixing weights come from
#'   `true_theta$p`.
#'
#' @return A `ggplot`.
#' @export
plot_purified_density_2d <- function(true_theta, n = 800L) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("`n` must be a single positive number.")
  }
  cts <- .celltype_names(true_theta)
  overlay <- .celltype_overlay(true_theta)
  lims <- .gene_axis_limits(list(true_theta))
  raster_df <- .purified_mixture_raster(
    true_theta,
    lims$xlim,
    lims$ylim
  )
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
      title = "Purified Gaussians (p-weighted mixture)",
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
    "Fill is the p-weighted mixture of purified Gaussians,",
    "sum_j p_j phi(x; mu_j, Sigma_j), not the bulk convolution.",
    "Outlines are exact 95% Gaussian regions",
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
  if (identical(which, "purified")) {
    raster_df <- dplyr::bind_rows(
      lapply(seq_along(thetas), function(i) {
        .purified_mixture_raster(
          thetas[[i]],
          lims$xlim,
          lims$ylim,
          panel = as.character(panels[[i]])
        )
      })
    )
  } else {
    df <- dplyr::bind_rows(
      lapply(seq_along(thetas), function(i) {
        .bulk_density_draws(
          thetas[[i]],
          n,
          panel = as.character(panels[[i]])
        )
      })
    )
    df$panel <- factor(df$panel, levels = levels(panels))
    raster_df <- .density_raster_table(df, lims)
  }
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

#' Contour of the bulk log-likelihood on hypothesised ratios
#'
#' For \eqn{J=2} the axes are \eqn{(p_1,p_2)} on the unit square (the
#' dashed line is the simplex). For \eqn{J=3} the axes are additive
#' log-ratio coordinates
#' \eqn{\rho_1=\ln(p_1/p_3)}, \eqn{\rho_2=\ln(p_2/p_3)}
#' ([additive_log_ratio()]). The true simulation proportions (MLE for
#' \eqn{y=\mu p^{\star}}) are marked.
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
  p_true <- as.numeric(true_theta$p)
  if (length(p_true) >= 3L) {
    return(.plot_bulk_loglik_surface_alr(true_theta, grid = grid, y = y))
  }
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

#' ALR-plane bulk log-likelihood lattice for \eqn{J=3}
#'
#' @keywords internal
#' @noRd
.alr_loglik_lattice <- function(true_theta, grid = 50L, y = NULL) {
  p_true <- as.numeric(true_theta$p)
  mu <- true_theta$mu
  Sigma <- true_theta$sigma
  if (length(p_true) < 3L) {
    stop(".alr_loglik_lattice() requires J >= 3.", call. = FALSE)
  }
  if (is.null(y)) {
    y <- drop(mu %*% p_true)
  }
  y <- as.numeric(y)
  rho_true <- as.numeric(additive_log_ratio(p_true))
  span <- max(4, abs(rho_true) + 1.5)
  rho1_seq <- seq(-span, span, length.out = grid)
  rho2_seq <- seq(-span, span, length.out = grid)
  z <- matrix(NA_real_, grid, grid)
  for (i in seq_len(grid)) {
    for (j in seq_len(grid)) {
      p_hat <- additive_logistic(c(rho1_seq[[i]], rho2_seq[[j]]))
      z[i, j] <- tryCatch(
        loglik_multivariate(p_hat, y, mu, Sigma),
        error = function(e) NA_real_
      )
    }
  }
  z_true <- tryCatch(
    loglik_multivariate(p_true, y, mu, Sigma),
    error = function(e) NA_real_
  )
  list(
    rho1 = rho1_seq,
    rho2 = rho2_seq,
    z = z,
    p_true = p_true,
    rho_true = rho_true,
    z_true = z_true,
    y = y
  )
}

#' @keywords internal
#' @noRd
.plot_bulk_loglik_surface_alr <- function(
  true_theta,
  grid = 50L,
  y = NULL
) {
  lat <- .alr_loglik_lattice(true_theta, grid = grid, y = y)
  grid_df <- expand.grid(rho1 = lat$rho1, rho2 = lat$rho2)
  grid_df$loglik <- as.vector(lat$z)
  true_df <- data.frame(
    rho1 = lat$rho_true[[1L]],
    rho2 = lat$rho_true[[2L]]
  )
  ggplot2::ggplot(
    grid_df,
    ggplot2::aes(x = .data[["rho1"]], y = .data[["rho2"]])
  ) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = .data[["loglik"]]),
      interpolate = TRUE
    ) +
    ggplot2::scale_fill_viridis_c(name = "log lik.") +
    ggplot2::geom_point(
      data = true_df,
      ggplot2::aes(x = .data[["rho1"]], y = .data[["rho2"]]),
      inherit.aes = FALSE,
      colour = "white",
      fill = "#E41A1C",
      shape = 21,
      size = 3.2,
      stroke = 0.8
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(2, 4, 2, 2),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = expression(rho[1] == log(p[1] / p[3])),
      y = expression(rho[2] == log(p[2] / p[3])),
      title = "Bulk log-likelihood (ALR)"
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

#' Tile heatmap of a scenario metric on the inner 2-by-2 design
#'
#' Fig02: one page per meta-scenario (composition \eqn{\times} variance
#' \eqn{\times} CLD) with tiles on \eqn{(\rho_1,\rho_2)} and one panel
#' per solver. Fig03: tiles on the CT1 / CT2 graph families.
#'
#' @param metrics Long table with `algorithm`, design columns, and `value`.
#' @param title Page title.
#' @return A `ggplot`.
#' @export
plot_bivariate_metric_tiles <- function(metrics, title) {
  hybrid <- .is_hybrid_config(metrics)
  n_alg <- dplyr::n_distinct(metrics$algorithm)
  ncol_alg <- if (isTRUE(hybrid)) {
    min(5L, n_alg)
  } else {
    4L
  }
  if (isTRUE(hybrid)) {
    x_col <- "graph_ct1"
    y_col <- "graph_ct2"
    x_lab <- "CT1 graph"
    y_lab <- "CT2 graph"
  } else {
    x_col <- "correlation_celltype1"
    y_col <- "correlation_celltype2"
    x_lab <- expression(rho[1])
    y_lab <- expression(rho[2])
  }
  p <- ggplot2::ggplot(
    metrics,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      fill = .data[["value"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(~algorithm, ncol = ncol_alg)
  if (!isTRUE(hybrid)) {
    p <- p + ggplot2::coord_equal()
  }
  p +
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
      x = x_lab,
      y = y_lab,
      fill = NULL,
      title = title
    )
}

#' @keywords internal
#' @noRd
.corner_thetas <- function(config, theta_tbl, row) {
  ids <- .corner_ids(config, row)
  out <- lapply(ids, function(id) {
    if (is.na(id) || !nzchar(id)) {
      return(NULL)
    }
    hit <- as.character(theta_tbl$ID) == as.character(id)
    if (!any(hit)) {
      return(NULL)
    }
    .unwrap_true_theta(theta_tbl$true_theta[which(hit)[[1L]]])
  })
  names(out) <- names(ids)
  out
}

#' @keywords internal
#' @noRd
.bivariate_page_meta <- function(config) {
  config <- .relevel_scenario_table(config)
  if (.is_hybrid_config(config)) {
    dplyr::distinct(
      config,
      .data[["proportions"]],
      .data[["overlap_label"]]
    ) |>
      dplyr::arrange(
        .data[["proportions"]],
        .data[["overlap_label"]]
      )
  } else {
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
}

#' Page title matching performance-book factor order
#'
#' @keywords internal
#' @noRd
.page_title <- function(row) {
  if ("overlap_label" %in% names(row)) {
    paste(
      as.character(row$proportions[[1L]]),
      as.character(row$overlap_label[[1L]]),
      sep = " / "
    )
  } else {
    cents <- as.character(row$centroids[[1L]])
    var <- as.character(row$variance[[1L]])
    prop <- as.character(row$proportions[[1L]])
    paste(cents, var, prop, sep = " / ")
  }
}

#' Rows of a long table that belong on one performance-book page
#'
#' @keywords internal
#' @noRd
.rows_on_page <- function(df, row) {
  if ("overlap_label" %in% names(df) && "overlap_label" %in% names(row)) {
    keep <- df$proportions == row$proportions[[1L]] &
      as.character(df$overlap_label) == as.character(row$overlap_label[[1L]])
  } else {
    keep <- df$proportions == row$proportions[[1L]] &
      df$variance == row$variance[[1L]] &
      df$centroids == row$centroids[[1L]]
  }
  df[keep, , drop = FALSE]
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
#' @param n Draws for the bulk KDE. Purified panels use the exact
#'   \eqn{p}-weighted Gaussian mixture and ignore `n`.
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
  hybrid <- .is_hybrid_config(agg)

  write_metric_pdf <- function(metric, path, stem) {
    meta <- .bivariate_page_meta(agg)
    pages <- list()
    pdf_w <- if (isTRUE(hybrid)) 18 else 16
    pdf_h <- if (isTRUE(hybrid)) 6 else 10
    grDevices::pdf(path, width = pdf_w, height = pdf_h)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (i in seq_len(nrow(meta))) {
      row <- meta[i, , drop = FALSE]
      page <- .rows_on_page(agg, row)
      page$value <- page[[metric]]
      page$page <- .page_title(row)
      title <- paste(toupper(metric), .page_title(row), sep = " / ")
      print(plot_bivariate_metric_tiles(page, title))
      pages[[i]] <- page
    }
    .write_ggplot_rds(dplyr::bind_rows(pages), data_rds, stem)
    path
  }

  out <- list(
    rmse = write_metric_pdf(
      "rmse",
      file.path(dir, "heatmap_rmse.pdf"),
      "heatmap_rmse"
    )
  )
  if (!isTRUE(hybrid)) {
    out$mae <- write_metric_pdf(
      "mae",
      file.path(dir, "heatmap_mae.pdf"),
      "heatmap_mae"
    )
    out$aitchison <- write_metric_pdf(
      "aitchison",
      file.path(dir, "heatmap_aitchison.pdf"),
      "heatmap_aitchison"
    )
  }
  out
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
      tibble::tibble(ID = ID, p_true = list(as.numeric(th$p)))
    }
  )
  scored <- dplyr::left_join(opt, p_true, by = "ID")
  ct_cols <- grep("^celltype_[0-9]+$", names(scored), value = TRUE)
  if (length(ct_cols) < 2L) {
    stop("Need celltype_* columns on `optimisation`.", call. = FALSE)
  }
  p_hat <- as.matrix(scored[, ct_cols, drop = FALSE])
  p_star <- do.call(rbind, scored$p_true)
  if (ncol(p_star) != ncol(p_hat)) {
    n_use <- min(ncol(p_star), ncol(p_hat))
    p_star <- p_star[, seq_len(n_use), drop = FALSE]
    p_hat <- p_hat[, seq_len(n_use), drop = FALSE]
  }
  scored$mae <- rowMeans(abs(p_hat - p_star), na.rm = TRUE)
  scored$rmse_row <- sqrt(rowMeans((p_hat - p_star)^2, na.rm = TRUE))
  scored$aitchison <- vapply(
    seq_len(nrow(scored)),
    function(i) {
      .aitchison_distance(p_star[i, ], p_hat[i, ])
    },
    numeric(1)
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
  if (.is_hybrid_config(config)) {
    panels <- .hybrid_topology_panels()
    vapply(
      names(panels),
      function(lab) {
        pair <- panels[[lab]]
        hit <- config$proportions == row$proportions &
          as.character(config$overlap_label) ==
            as.character(row$overlap_label) &
          as.character(config$graph_ct1) == pair[[1L]] &
          as.character(config$graph_ct2) == pair[[2L]]
        if (!any(hit)) {
          return(NA_character_)
        }
        as.character(config$ID[which(hit)[[1L]]])
      },
      character(1)
    )
  } else {
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
  panel_lvls <- unique(as.character(plot_df$panel))
  hybrid_lvls <- names(.hybrid_topology_panels())
  canon <- if (any(panel_lvls %in% hybrid_lvls)) {
    hybrid_lvls
  } else {
    names(.bivariate_corr_corners())
  }
  plot_df$panel <- factor(
    plot_df$panel,
    levels = c(intersect(canon, panel_lvls), setdiff(panel_lvls, canon))
  )
  plot_df$cell_type <- factor(
    plot_df$cell_type,
    levels = unique(as.character(plot_df$cell_type))
  )
  pal <- .cell_type_colours(levels(plot_df$cell_type))
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
  if (n_ct > 2L) {
    annot_df <- plot_df
    hi <- pmax(plot_df$emp_hi, plot_df$ci_hi, plot_df$mean_est, na.rm = TRUE)
    annot_df$lab_x <- ifelse(is.finite(hi), hi, plot_df$mean_est)
    annot_df$annot <- paste0(
      sprintf("%.3f", annot_df$rmse),
      " / ",
      ifelse(
        is.finite(annot_df$coverage),
        sprintf("%.0f%%", 100 * annot_df$coverage),
        "NA"
      )
    )
    annot_df$y_lab <- annot_df$y_dodge
    annot_size <- 2.4
  } else {
    annot_df <- plot_df |>
      dplyr::group_by(.data[["panel"]], .data[["algorithm"]]) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::mutate(y_lab = as.numeric(.data[["algorithm"]]))
    annot_size <- 3.4
  }
  caption_txt <- if (n_ct > 2L) {
    paste(
      "Solid whiskers: mean estimate plus or minus 1.96 times the",
      "expected-Fisher Wald SE at the true composition when available.",
      "Dashed whiskers: plus or minus 1.96 times the empirical",
      "Monte Carlo SD.",
      "Vertical dashed lines: true proportions (one colour per type).",
      "Labels: RMSE and Wald coverage for each cell type."
    )
  } else {
    paste(
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
  }
  p <- ggplot2::ggplot(
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
    ggplot2::geom_point(size = 2.6)

  if (n_ct > 2L) {
    p <- p +
      ggplot2::geom_text(
        data = annot_df,
        ggplot2::aes(
          x = .data[["lab_x"]],
          y = .data[["y_lab"]],
          label = .data[["annot"]],
          colour = .data[["cell_type"]]
        ),
        inherit.aes = FALSE,
        hjust = -0.08,
        vjust = 0.5,
        size = annot_size,
        fontface = "bold",
        lineheight = 0.9,
        show.legend = FALSE
      )
    x_expand <- c(0.04, 0.22)
  } else {
    p <- p +
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
        size = annot_size,
        fontface = "bold",
        lineheight = 0.95,
        label.size = 0.25,
        label.padding = grid::unit(0.28, "lines"),
        fill = "#F4F1EA",
        colour = "grey20",
        show.legend = FALSE
      )
    x_expand <- c(0.04, 0.38)
  }
  p +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(alg_lvls),
      labels = alg_lvls,
      expand = ggplot2::expansion(add = 0.45)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = x_expand)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::facet_wrap(~panel, ncol = 2L) +
    ggplot2::labs(
      x = "Mean Monte Carlo estimate",
      y = NULL,
      colour = "Cell type",
      title = title,
      caption = caption_txt
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
.bivariate_wald_forest_table <- function(
  artefacts,
  ids,
  labels,
  filter_wald = TRUE
) {
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
    if (isTRUE(filter_wald)) {
      keep <- .algorithm_in(mc$algorithm, .bivariate_wald_algorithms())
      mc <- mc[keep, , drop = FALSE]
      if (nrow(mc) == 0L) {
        return(NULL)
      }
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
  hybrid <- .is_hybrid_config(df)
  if (isTRUE(hybrid)) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["graph_ct1"]],
        y = .data[["graph_ct2"]],
        colour = .data[["rmse"]],
        size = .data[["aitchison"]]
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(
        ~algorithm,
        ncol = min(5L, dplyr::n_distinct(df$algorithm))
      ) +
      ggplot2::scale_colour_viridis_c(name = "RMSE") +
      ggplot2::scale_size_continuous(
        name = "Aitchison",
        range = c(1.5, 10)
      ) +
      ggplot2::labs(
        x = "CT1 graph",
        y = "CT2 graph",
        title = title,
        caption = paste(
          "Each panel is the 2-by-2 of CT1 / CT2 graph families.",
          "Colour: mean RMSE; size: mean Aitchison distance."
        )
      )
  } else {
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
    p <- ggplot2::ggplot(
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
      )
  }
  p +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
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
  filter_wald <- !.is_hybrid_config(cfg)
  collected <- list()
  grDevices::pdf(file, width = 18, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    ids <- .corner_ids(cfg, row)
    tbl <- .bivariate_wald_forest_table(
      artefacts,
      ids,
      names(ids),
      filter_wald = filter_wald
    )
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
  hybrid <- .is_hybrid_config(cfg)
  collected <- list()
  pdf_h <- if (isTRUE(hybrid)) 22 else 30
  slab_s <- if (isTRUE(hybrid)) 1.45 else 1.05
  slab_a <- if (isTRUE(hybrid)) 0.35 else 0.45
  grDevices::pdf(file, width = 16, height = pdf_h)
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
      slab_scale = slab_s,
      slab_alpha = slab_a
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
  hybrid <- .is_hybrid_config(agg)
  pdf_w <- if (isTRUE(hybrid)) 18 else 16
  pdf_h <- if (isTRUE(hybrid)) 6 else 10
  grDevices::pdf(file, width = pdf_w, height = pdf_h)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    page <- .rows_on_page(agg, row)
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

#' Mean of SPD condition numbers of the precision matrices \eqn{\Theta_j}
#'
#' @keywords internal
#' @noRd
.mean_kappa_precision <- function(th) {
  Sigma <- th$sigma
  Theta <- th$Theta
  j_dim <- dim(Sigma)[[3L]]
  kappas <- vapply(
    seq_len(j_dim),
    function(j) {
      mat <- NULL
      if (!is.null(Theta)) {
        mat <- Theta[,, j]
      }
      if (is.null(mat)) {
        mat <- tryCatch(
          solve(Sigma[,, j]),
          error = function(e) NULL
        )
      }
      if (is.null(mat)) {
        return(NA_real_)
      }
      ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
      ev <- ev[is.finite(ev) & ev > 0]
      if (length(ev) < 1L) {
        return(NA_real_)
      }
      max(ev) / min(ev)
    },
    numeric(1)
  )
  mean(kappas, na.rm = TRUE)
}

#' Pairwise Hellinger and Jeffreys scores among cell-type Gaussians
#'
#' @keywords internal
#' @noRd
.pairwise_gaussian_scores <- function(th) {
  mu <- th$mu
  Sigma <- th$sigma
  j_dim <- ncol(mu)
  ct <- colnames(mu)
  if (is.null(ct)) {
    ct <- paste0("CT", seq_len(j_dim))
  }
  rows <- list()
  k <- 1L
  for (j in seq_len(j_dim - 1L)) {
    for (ell in seq.int(j + 1L, j_dim)) {
      rows[[k]] <- tibble::tibble(
        pair = paste(ct[[j]], ct[[ell]], sep = "--"),
        i = j,
        l = ell,
        hellinger = .hellinger_gaussian(
          mu[, j],
          mu[, ell],
          Sigma[,, j],
          Sigma[,, ell]
        ),
        jeffreys = tryCatch(
          .jeffreys_gaussian(
            mu[, j],
            mu[, ell],
            Sigma[,, j],
            Sigma[,, ell]
          ),
          error = function(e) NA_real_
        )
      )
      k <- k + 1L
    }
  }
  dplyr::bind_rows(rows)
}

#' Adjacency list (one matrix per cell type) from `true_theta`
#'
#' @keywords internal
#' @noRd
.theta_adjacency <- function(th) {
  adj <- th$adjacency
  if (is.null(adj)) {
    return(NULL)
  }
  if (is.list(adj)) {
    return(adj)
  }
  if (is.array(adj) && length(dim(adj)) == 3L) {
    lapply(seq_len(dim(adj)[[3L]]), function(j) adj[,, j])
  } else {
    NULL
  }
}

#' Scenario-level geometry metrics for the funkyheatmap
#'
#' @keywords internal
#' @noRd
.scenario_global_metrics <- function(artefacts, n_mc = 1500L) {
  cfg <- .relevel_scenario_table(artefacts$config)
  theta_tbl <- artefacts$theta
  desc <- artefacts$descriptors
  if (is.null(theta_tbl) || !"true_theta" %in% names(theta_tbl)) {
    stop("`artefacts$theta$true_theta` is required.", call. = FALSE)
  }
  desc_join <- NULL
  if (!is.null(desc) && "ID" %in% names(desc)) {
    keep <- intersect(
      c("ID", "mixsim_baromega", "network_density", "jeffreys", "h_star"),
      names(desc)
    )
    desc_join <- desc[, keep, drop = FALSE]
  }
  rows <- lapply(seq_len(nrow(cfg)), function(i) {
    id <- as.character(cfg$ID[[i]])
    th <- .unwrap_true_theta(
      theta_tbl$true_theta[as.character(theta_tbl$ID) == id][[1L]]
    )
    p <- as.numeric(th$p)
    pair <- .pairwise_gaussian_scores(th)
    ovl <- tryCatch(
      overlap_gaussian_mc(
        true_theta = th,
        n_mc = n_mc,
        seed = 20260807L + i
      ),
      error = function(e) NULL
    )
    tibble::tibble(
      ID = id,
      shannon_entropy = compute_shannon_entropy(p),
      avg_overlap = if (!is.null(ovl)) ovl$BarOmega else NA_real_,
      max_overlap = if (!is.null(ovl)) ovl$MaxOmega else NA_real_,
      max_kl = if (any(is.finite(pair$jeffreys))) {
        max(pair$jeffreys, na.rm = TRUE)
      } else {
        NA_real_
      },
      kappa_precision = .mean_kappa_precision(th),
      edge_density = NA_real_
    )
  })
  out <- dplyr::bind_rows(rows)
  if (!is.null(desc_join)) {
    out <- dplyr::left_join(out, desc_join, by = "ID")
    if ("mixsim_baromega" %in% names(out)) {
      out$avg_overlap <- ifelse(
        is.finite(out$avg_overlap),
        out$avg_overlap,
        out$mixsim_baromega
      )
    }
    if ("network_density" %in% names(out)) {
      out$edge_density <- out$network_density
    }
  }
  dplyr::left_join(out, cfg, by = "ID")
}

#' Tile heatmap of the unscaled mean signature \eqn{\mu}
#'
#' @param true_theta A `true_theta` list (`mu` is \eqn{G\times J}).
#' @param file Output PDF path.
#' @param data_rds Optional directory for the ggplot `data` RDS.
#' @return `file`, invisibly.
#' @export
save_mean_signature_heatmap <- function(
  true_theta,
  file,
  data_rds = NULL
) {
  th <- .unwrap_true_theta(true_theta)
  mu <- as.matrix(th$mu)
  genes <- rownames(mu)
  if (is.null(genes)) {
    genes <- paste0("g", seq_len(nrow(mu)))
  }
  cts <- colnames(mu)
  if (is.null(cts)) {
    cts <- paste0("celltype_", seq_len(ncol(mu)))
  }
  df <- expand.grid(
    gene = factor(genes, levels = rev(genes)),
    cell_type = factor(cts, levels = cts),
    stringsAsFactors = FALSE
  )
  df$value <- as.vector(mu)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["cell_type"]],
      y = .data[["gene"]],
      fill = .data[["value"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_viridis_c(name = "mean") +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Gene",
      title = "Mean signature (unscaled)",
      caption = paste(
        "Raw mu (G genes by J cell types); no row or column scaling."
      )
    )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file,
    plot = p,
    width = 7,
    height = 8,
    dpi = 320,
    limitsize = FALSE
  )
  .write_ggplot_rds(df, data_rds, "mean_signature")
  invisible(file)
}

#' Funky heatmap of scenario-level geometry metrics
#'
#' Rows follow the fig03 factor order (proportions, CT1/CT2 graphs,
#' MixSim overlap). No clustering.
#'
#' @inheritParams save_bivariate_similarity_book
#' @export
save_scenario_metrics_funkyheatmap <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  .check_suggested_package("funkyheatmap", "save_scenario_metrics_funkyheatmap")
  metrics <- .scenario_global_metrics(artefacts)
  metrics <- .relevel_scenario_table(metrics)
  metrics <- dplyr::arrange(
    metrics,
    .data[["proportions"]],
    .data[["graph_ct1"]],
    .data[["graph_ct2"]],
    .data[["overlap_label"]]
  )
  metrics$id <- as.character(metrics$ID)
  graph_short <- function(x) {
    dplyr::case_when(
      as.character(x) == "Scale-free" ~ "SF",
      as.character(x) == "Cluster SBM" ~ "SBM",
      TRUE ~ as.character(x)
    )
  }
  metrics$label <- paste(
    paste(
      graph_short(metrics$graph_ct1),
      graph_short(metrics$graph_ct2),
      sep = "/"
    ),
    as.character(metrics$overlap_label),
    sep = " | "
  )
  data <- tibble::tibble(
    id = metrics$id,
    scenario = metrics$label,
    shannon_entropy = metrics$shannon_entropy,
    avg_overlap = metrics$avg_overlap,
    max_overlap = metrics$max_overlap,
    max_kl = metrics$max_kl,
    kappa_precision = metrics$kappa_precision,
    edge_density = metrics$edge_density
  )
  column_info <- tibble::tibble(
    id = c(
      "scenario",
      "shannon_entropy",
      "avg_overlap",
      "max_overlap",
      "max_kl",
      "kappa_precision",
      "edge_density"
    ),
    name = c(
      "Scenario",
      "Shannon entropy",
      "Average overlap",
      "Maximum overlap",
      "Largest KL",
      "Precision kappa",
      "Mean edge density"
    ),
    geom = c("text", rep("bar", 6L)),
    group = c("id", rep("geometry", 6L)),
    palette = c("Greys", rep("Blues", 6L)),
    legend = c(FALSE, rep(TRUE, 6L)),
    width = c(4, rep(1.2, 6L))
  )
  row_info <- tibble::tibble(
    id = metrics$id,
    group = as.character(metrics$proportions)
  )
  row_groups <- tibble::tibble(
    group = unique(as.character(metrics$proportions)),
    level1 = unique(as.character(metrics$proportions))
  )
  column_groups <- tibble::tibble(
    group = c("id", "geometry"),
    palette = c("Greys", "Blues"),
    level1 = c("Id", "Geometry")
  )
  palettes <- list(Blues = "Blues", Greys = "Greys")
  p <- funkyheatmap::funky_heatmap(
    data = data,
    column_info = column_info,
    row_info = row_info,
    row_groups = row_groups,
    column_groups = column_groups,
    palettes = palettes,
    scale_column = TRUE,
    add_abc = FALSE,
    position_args = funkyheatmap::position_arguments(
      col_annot_offset = 3.5,
      expand_xmax = 2
    )
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  w <- max(if (!is.null(p$width)) p$width else 12, 14)
  h <- max(if (!is.null(p$height)) p$height else 14, 16)
  ggplot2::ggsave(
    file,
    plot = p,
    width = w,
    height = h,
    dpi = 320,
    limitsize = FALSE
  )
  .write_ggplot_rds(metrics, data_rds, "scenario_metrics")
  invisible(file)
}

#' 3-by-3 Hellinger heatmaps of cell-type pairs versus MixSim overlap
#'
#' Facets are pairwise comparisons (rows) by overlap (columns). Each
#' panel is the 2-by-2 of CT1 / CT2 graph families (CT3 is Scale-free).
#' Hellinger does not depend on \eqn{p}, so compositions are collapsed.
#'
#' @inheritParams save_bivariate_similarity_book
#' @export
save_pairwise_network_distance_heatmap <- function(
  artefacts,
  file,
  data_rds = NULL
) {
  cfg <- .relevel_scenario_table(artefacts$config)
  theta_tbl <- artefacts$theta
  pick <- cfg
  if ("proportions" %in% names(pick)) {
    bal <- as.character(pick$proportions) == "balanced"
    if (any(bal)) {
      pick <- pick[bal, , drop = FALSE]
    }
  }
  pieces <- lapply(seq_len(nrow(pick)), function(i) {
    id <- as.character(pick$ID[[i]])
    th <- .unwrap_true_theta(
      theta_tbl$true_theta[as.character(theta_tbl$ID) == id][[1L]]
    )
    pair <- .pairwise_gaussian_scores(th)
    pair$ID <- id
    pair$graph_ct1 <- pick$graph_ct1[[i]]
    pair$graph_ct2 <- pick$graph_ct2[[i]]
    pair$overlap_label <- pick$overlap_label[[i]]
    pair
  })
  df <- dplyr::bind_rows(pieces)
  pair_map <- c(
    "celltype_1--celltype_2" = "CT1--CT2",
    "celltype_1--celltype_3" = "CT1--CT3",
    "celltype_2--celltype_3" = "CT2--CT3",
    "CT1--CT2" = "CT1--CT2",
    "CT1--CT3" = "CT1--CT3",
    "CT2--CT3" = "CT2--CT3"
  )
  mapped <- unname(pair_map[as.character(df$pair)])
  df$pair <- ifelse(is.na(mapped), as.character(df$pair), mapped)
  df$pair <- factor(
    df$pair,
    levels = c("CT1--CT2", "CT1--CT3", "CT2--CT3")
  )
  df <- .relevel_scenario_table(df)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["graph_ct1"]],
      y = .data[["graph_ct2"]],
      fill = .data[["hellinger"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(pair ~ overlap_label) +
    ggplot2::scale_fill_viridis_c(name = "Hellinger") +
    ggplot2::coord_equal() +
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
      x = "CT1 graph",
      y = "CT2 graph",
      title = "Pairwise Hellinger distance of cell-type Gaussians",
      caption = paste(
        "Each facet is one cell-type pair (rows) and MixSim overlap",
        "(columns). Inner 2-by-2: CT1 / CT2 topologies; CT3 is Scale-free."
      )
    )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file,
    plot = p,
    width = 12,
    height = 11,
    dpi = 320,
    limitsize = FALSE
  )
  .write_ggplot_rds(df, data_rds, "pairwise_hellinger")
  invisible(file)
}

#' Absolute covariance on the undirected skeleton, shared colour scale
#'
#' @keywords internal
#' @noRd
.sigma_edge_weights <- function(adj, sigma) {
  adj <- as.matrix(adj)
  sigma <- as.matrix(sigma)
  w <- abs(sigma)
  w[as.matrix(adj) == 0] <- 0
  diag(w) <- 0
  w
}

#' Draw one igraph with edge widths proportional to |Sigma|
#'
#' @keywords internal
#' @noRd
.plot_weighted_skeleton <- function(
  adj,
  sigma,
  layout,
  vertex_col,
  w_max,
  main
) {
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    diag = FALSE,
    weighted = NULL
  )
  el <- igraph::as_edgelist(g, names = FALSE)
  w <- .sigma_edge_weights(adj, sigma)
  ew <- if (nrow(el) == 0L) {
    numeric(0)
  } else {
    vapply(
      seq_len(nrow(el)),
      function(k) w[el[k, 1L], el[k, 2L]],
      numeric(1)
    )
  }
  width <- if (is.finite(w_max) && w_max > 0) {
    0.35 + 6 * sqrt(pmax(ew, 0) / w_max)
  } else {
    rep(0.8, length(ew))
  }
  igraph::V(g)$color <- vertex_col
  igraph::plot.igraph(
    g,
    layout = layout,
    vertex.size = 8,
    vertex.label = NA,
    vertex.frame.color = "#2f3e4f",
    edge.color = "#4a5560",
    edge.width = width,
    main = main
  )
}

#' Three-page network book (increasing MixSim overlap)
#'
#' Each page: CT3 (Scale-free) centred at the top; 2-by-2 of CT1 / CT2
#' topologies below. Edge widths use |Sigma_ij| on the skeleton with one
#' shared scale across pages.
#'
#' @inheritParams save_bivariate_similarity_book
#' @param png_file Optional PNG copy (vignette / README).
#' @export
save_network_topology_book <- function(
  artefacts,
  file,
  data_rds = NULL,
  png_file = NULL
) {
  .check_suggested_package("igraph", "save_network_topology_book")
  cfg <- .relevel_scenario_table(artefacts$config)
  theta_tbl <- artefacts$theta
  if ("proportions" %in% names(cfg)) {
    bal <- as.character(cfg$proportions) == "balanced"
    if (any(bal)) {
      cfg <- cfg[bal, , drop = FALSE]
    }
  }
  ovl_lvls <- levels(cfg$overlap_label)
  if (is.null(ovl_lvls)) {
    ovl_lvls <- unique(as.character(cfg$overlap_label))
  }
  pal <- c(
    "Scale-free" = "#4C72B0",
    "Cluster SBM" = "#C44E52"
  )
  layouts <- list()
  w_max <- 0
  pages <- lapply(ovl_lvls, function(ovl) {
    sub <- cfg[as.character(cfg$overlap_label) == ovl, , drop = FALSE]
    panels <- .hybrid_topology_panels()
    slot_list <- list()
    for (lab in names(panels)) {
      pair <- panels[[lab]]
      hit <- as.character(sub$graph_ct1) == pair[[1L]] &
        as.character(sub$graph_ct2) == pair[[2L]]
      if (!any(hit)) {
        next
      }
      id <- as.character(sub$ID[which(hit)[[1L]]])
      th <- .unwrap_true_theta(
        theta_tbl$true_theta[as.character(theta_tbl$ID) == id][[1L]]
      )
      adj <- .theta_adjacency(th)
      slot_list[[lab]] <- list(
        th = th,
        adj = adj,
        graph_ct1 = pair[[1L]],
        graph_ct2 = pair[[2L]]
      )
      if (!is.null(adj)) {
        for (j in seq_along(adj)) {
          w_max <<- max(
            w_max,
            max(.sigma_edge_weights(adj[[j]], th$sigma[,, j]), na.rm = TRUE)
          )
        }
      }
    }
    list(overlap = ovl, slots = slot_list)
  })
  layout_key <- function(j, family) {
    paste(j, family, sep = "::")
  }
  ensure_layout <- function(j, family, adj) {
    key <- layout_key(j, family)
    if (is.null(layouts[[key]])) {
      g <- igraph::graph_from_adjacency_matrix(
        adj,
        mode = "undirected",
        diag = FALSE
      )
      layouts[[key]] <<- igraph::layout_with_fr(g)
    }
    layouts[[key]]
  }
  draw_one <- function(slot, j, family, main) {
    .plot_weighted_skeleton(
      slot$adj[[j]],
      slot$th$sigma[,, j],
      ensure_layout(j, family, slot$adj[[j]]),
      pal[[family]],
      w_max,
      main
    )
  }
  draw_page <- function(page) {
    graphics::layout(matrix(
      c(0, 1, 1, 0, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = 3L,
      byrow = TRUE
    ))
    graphics::par(mar = c(0.3, 0.2, 2.1, 0.2))
    first <- page$slots[[1L]]
    draw_one(first, 3L, "Scale-free", "CT3 Scale-free (held)")
    labs <- names(.hybrid_topology_panels())
    for (lab in labs) {
      slot <- page$slots[[lab]]
      pair <- .hybrid_topology_panels()[[lab]]
      if (is.null(slot)) {
        graphics::plot.new()
        graphics::title(paste("CT1", pair[[1L]]))
        graphics::plot.new()
        graphics::title(paste("CT2", pair[[2L]]))
        next
      }
      draw_one(slot, 1L, pair[[1L]], paste("CT1", pair[[1L]]))
      draw_one(slot, 2L, pair[[2L]], paste("CT2", pair[[2L]]))
    }
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 12, height = 11)
  for (i in seq_along(pages)) {
    graphics::par(oma = c(0.2, 0.2, 2.4, 0.2))
    draw_page(pages[[i]])
    graphics::mtext(
      pages[[i]]$overlap,
      outer = TRUE,
      cex = 1.3,
      font = 2
    )
  }
  grDevices::dev.off()
  if (!is.null(png_file) && length(pages) > 0L) {
    dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
    mid <- min(2L, length(pages))
    grDevices::png(png_file, width = 2400, height = 2200, res = 220)
    graphics::par(oma = c(0.2, 0.2, 2.4, 0.2))
    draw_page(pages[[mid]])
    graphics::mtext(
      pages[[mid]]$overlap,
      outer = TRUE,
      cex = 1.2,
      font = 2
    )
    grDevices::dev.off()
  }
  meta <- tibble::tibble(
    overlap = vapply(pages, `[[`, character(1), "overlap"),
    n_slots = vapply(pages, function(p) length(p$slots), integer(1)),
    w_max = w_max
  )
  .write_ggplot_rds(meta, data_rds, "network_topologies")
  invisible(file)
}
