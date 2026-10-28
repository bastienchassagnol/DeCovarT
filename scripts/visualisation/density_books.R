# Density books: the simulated laws (purified and bulk densities, mean
# signature, networks, latent projections, pairwise distances) and the
# log-likelihood pictures. Every `.data_*` function builds the table of one
# book from the benchmark; every `plot_*` function draws the pages of that
# book from the table alone.

# Purified and bulk densities (G = 2) ------------------------------------------

.dmvnorm_rows <- function(xy, mu, sigma) {
  z <- sweep(xy, 2L, mu, "-")
  quad <- rowSums(z * t(solve(sigma, t(z))))
  log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  exp(-0.5 * (ncol(xy) * log(2 * pi) + log_det + quad))
}

.ellipse_95 <- function(mu, sigma, n = 120L) {
  angle <- seq(0, 2 * pi, length.out = n)
  root <- t(chol(sigma)) * sqrt(stats::qchisq(0.95, df = 2L))
  xy <- t(mu + root %*% rbind(cos(angle), sin(angle)))
  data.frame(gene_1 = xy[, 1L], gene_2 = xy[, 2L])
}

# Exact densities on a grid shared by every page: the p-weighted mixture of
# the purified Gaussians, and the convolution N(mu p, sum_j p_j^2 Sigma_j)
# that DeCovarT fits; 95% ellipses and centroids of the purified types.
.data_gene_densities <- function(benchmark, layout, n_grid = 90L) {
  corners <- layout[layout$corner, ]
  thetas <- benchmark$theta_true[corners$ID]
  centres <- do.call(cbind, lapply(thetas, `[[`, "mu"))
  spread <- 4 *
    sqrt(max(unlist(lapply(thetas, function(th) {
      apply(th$sigma, 3L, diag)
    }))))
  axis <- seq(min(centres) - spread, max(centres) + spread, length.out = n_grid)
  grid <- as.matrix(expand.grid(gene_1 = axis, gene_2 = axis))
  densities <- purrr::map_dfr(seq_len(nrow(corners)), function(i) {
    th <- thetas[[i]]
    p <- as.numeric(th$p)
    purified <- Reduce(
      `+`,
      lapply(seq_along(p), function(j) {
        p[[j]] * .dmvnorm_rows(grid, th$mu[, j], th$sigma[,, j])
      })
    )
    bulk <- .dmvnorm_rows(
      grid,
      drop(th$mu %*% p),
      .compute_global_variance(p, th$sigma)
    )
    tibble::tibble(
      ID = corners$ID[[i]],
      gene_1 = grid[, 1L],
      gene_2 = grid[, 2L],
      purified = purified,
      bulk = bulk
    )
  })
  types <- purrr::map_dfr(seq_len(nrow(corners)), function(i) {
    th <- thetas[[i]]
    purrr::map_dfr(seq_len(ncol(th$mu)), function(j) {
      tibble::tibble(
        ID = corners$ID[[i]],
        cell_type = paste("Cell type", j),
        .ellipse_95(th$mu[, j], th$sigma[,, j]),
        centre_1 = th$mu[1L, j],
        centre_2 = th$mu[2L, j]
      )
    })
  })
  list(
    densities = .with_layout(densities, layout),
    types = .with_layout(types, layout),
    limits = range(axis)
  )
}

plot_gene_densities <- function(data, which) {
  title <- c(
    purified = "Purified cell-type laws, weighted by the composition",
    bulk = "Bulk convolution N(mu p, sum_j p_j^2 Sigma_j)"
  )[[which]]
  caption <- c(
    purified = paste(
      "Fill: exact density of sum_j p_j N(mu_j, Sigma_j), not the bulk law.",
      "Outlines: 95% regions of each purified type; points: centroids."
    ),
    bulk = paste(
      "Fill: exact density of the bulk law; points: centroids of the",
      "purified types. Same gene axes on every page."
    )
  )[[which]]
  lapply(levels(droplevels(data$densities$page)), function(page) {
    dens <- data$densities[data$densities$page == page, ]
    dens$density <- dens[[which]]
    types <- data$types[data$types$page == page, ]
    centres <- dplyr::distinct(
      types,
      .data$panel,
      .data$cell_type,
      .data$centre_1,
      .data$centre_2
    )
    p <- ggplot2::ggplot(dens, ggplot2::aes(.data$gene_1, .data$gene_2)) +
      ggplot2::geom_raster(
        ggplot2::aes(fill = .data$density),
        interpolate = TRUE
      ) +
      ggplot2::scale_fill_viridis_c(name = "Density") +
      ggplot2::geom_point(
        data = centres,
        ggplot2::aes(.data$centre_1, .data$centre_2, shape = .data$cell_type),
        colour = "white",
        size = 2.6
      ) +
      ggplot2::scale_shape_manual(
        name = NULL,
        values = c(16, 17),
        guide = ggplot2::guide_legend(override.aes = list(colour = "black"))
      ) +
      ggplot2::facet_wrap(~panel, ncol = 2L) +
      ggplot2::coord_equal(
        xlim = data$limits,
        ylim = data$limits,
        expand = FALSE
      ) +
      ggplot2::labs(
        x = "Gene 1",
        y = "Gene 2",
        title = title,
        subtitle = page,
        caption = caption
      ) +
      .page_theme()
    if (identical(which, "purified")) {
      p <- p +
        ggplot2::geom_path(
          data = types,
          ggplot2::aes(group = .data$cell_type, linetype = .data$cell_type),
          colour = "white",
          linewidth = 0.5
        ) +
        ggplot2::scale_linetype_manual(
          name = NULL,
          values = c(1, 2),
          guide = ggplot2::guide_legend(override.aes = list(colour = "black"))
        )
    }
    p
  })
}

# Log-likelihood lattices -------------------------------------------------------

# Expected log-likelihood E_{p*} l(p) of the convolution, without the
# constant -G/2 log(2 pi): maximised at p* exactly (information inequality).
.expected_loglik_population <- function(p, p_true, mu, Sigma) {
  model <- .sigma_p_factorisation(p, Sigma)
  truth <- .sigma_p_factorisation(p_true, Sigma)
  gap <- backsolve(model$chol, drop(mu %*% (p_true - p)), transpose = TRUE)
  -0.5 * (model$log_det + sum(gap^2) + sum(model$inverse * truth$matrix))
}

# Single-cell covariance mixture sum_j p_j Sigma_j with its Cholesky factor,
# log-determinant and inverse.
.single_cell_mixture <- function(p, Sigma) {
  mixture <- apply(sweep(Sigma, 3L, p, "*"), c(1L, 2L), sum)
  chol_factor <- chol(mixture)
  list(
    matrix = mixture,
    log_det = 2 * sum(log(diag(chol_factor))),
    inverse = chol2inv(chol_factor)
  )
}

# Expected profile log-likelihood of the single-cell model: E l(p, c, phi)
# is maximised in (c, phi) in closed form, as the observed profile is.
.expected_loglik_single_cell <- function(
  p,
  p_true,
  mu,
  Sigma,
  c_true,
  phi_true
) {
  model <- .single_cell_mixture(p, Sigma)
  truth <- .single_cell_mixture(p_true, Sigma)
  m <- drop(mu %*% p)
  m_true <- c_true * drop(mu %*% p_true)
  w_m <- drop(model$inverse %*% m)
  scale <- sum(w_m * m_true) / sum(w_m * m)
  gap <- m_true - scale * m
  phi <- (phi_true *
    sum(model$inverse * truth$matrix) +
    sum(gap * drop(model$inverse %*% gap))) /
    length(m)
  -0.5 * (length(m) * log(phi) + model$log_det + length(m))
}

# Relative log-likelihood log10(L / max L) of `f` on an n x n lattice of
# the window `xlim` x `ylim`, floored at .loglik_floor.
.relative_lattice <- function(f, xlim, ylim, n_grid) {
  grid <- expand.grid(
    x = seq(xlim[[1L]], xlim[[2L]], length.out = n_grid),
    y = seq(ylim[[1L]], ylim[[2L]], length.out = n_grid)
  )
  value <- mapply(f, grid$x, grid$y)
  grid$value <- pmax(
    (value - max(value, na.rm = TRUE)) / log(10),
    .loglik_floor
  )
  grid
}

# Zoom of the single-cell lattices: the peak narrows as the bulk sums more
# cells, so the window shrinks to the cells above the floor, plus one cell
# on each side, until that region fills the window; the last window is
# padded by 10%.
.zoomed_lattice <- function(f, n_grid, max_zoom = 6L) {
  xlim <- .balance_limits
  ylim <- .balance_limits
  for (k in seq_len(max_zoom)) {
    grid <- .relative_lattice(f, xlim, ylim, n_grid)
    above <- grid[grid$value > .loglik_floor, ]
    step <- c(diff(xlim), diff(ylim)) / (n_grid - 1L)
    box_x <- range(above$x) + c(-1, 1) * step[[1L]]
    box_y <- range(above$y) + c(-1, 1) * step[[2L]]
    if (diff(box_x) > 0.6 * diff(xlim) && diff(box_y) > 0.6 * diff(ylim)) {
      break
    }
    xlim <- box_x + c(-1, 1) * 0.1 * diff(box_x)
    ylim <- box_y + c(-1, 1) * 0.1 * diff(box_y)
  }
  grid
}

# Observed (first simulated bulk) and expected log-likelihood of one
# scenario on a lattice. fig02: unnormalised weights (p1, p2) of the
# convolution; fig03: the plane of the two balances; fig05: the plane of
# the two within-lineage balances, the between-lineage balance held at
# its true value, each surface on its own zoomed window. The Monte Carlo
# mean of the Newton-Raphson estimates is the marker of the expected
# surface.
.loglik_lattice <- function(benchmark, id, figure, n_grid) {
  th <- benchmark$theta_true[[id]]
  p_true <- as.numeric(th$p)
  y <- benchmark$data[[id]]$bulk[, 1L]
  single_cell <- identical(figure, "fig05")
  estimates <- benchmark$optimisation[
    benchmark$optimisation$ID == id &
      benchmark$optimisation$algorithm == "Newton-Raphson",
    colnames(th$mu)
  ]
  estimates <- as.matrix(estimates[stats::complete.cases(estimates), ])
  if (identical(figure, "fig02")) {
    xlim <- ylim <- c(0.01, 0.99)
    to_p <- function(x1, x2) c(x1, x2)
    truth_xy <- p_true
    estimate_xy <- estimates
  } else {
    balance <- .figure_balances(figure)
    V <- balance$basis
    z_true <- drop(crossprod(V, log(p_true)))
    xlim <- ylim <- .balance_limits
    to_p <- function(x1, x2) isometric_logistic(c(x1, x2, z_true[-(1:2)]), V)
    truth_xy <- z_true[1:2]
    interior <- estimates[apply(estimates > 0, 1L, all), , drop = FALSE]
    estimate_xy <- t(crossprod(V, t(log(interior))))[, 1:2, drop = FALSE]
  }
  if (single_cell) {
    nuisance <- benchmark$data[[id]]$nuisance
    observed <- function(p) {
      profile_loglik_single_cell_ilr(isometric_log_ratio(p), y, th$mu, th$sigma)
    }
    expected <- function(p) {
      .expected_loglik_single_cell(
        p,
        p_true,
        th$mu,
        th$sigma,
        nuisance["mean_scale", 1L],
        nuisance["variance_scale", 1L]
      )
    }
  } else {
    observed <- function(p) loglik_multivariate(p, y, th$mu, th$sigma)
    expected <- function(p) {
      .expected_loglik_population(p, p_true, th$mu, th$sigma)
    }
  }
  on_lattice <- function(f) {
    g <- function(x1, x2) {
      tryCatch(f(to_p(x1, x2)), error = function(e) NA_real_)
    }
    if (single_cell) {
      .zoomed_lattice(g, n_grid)
    } else {
      .relative_lattice(g, xlim, ylim, n_grid)
    }
  }
  surface <- dplyr::bind_rows(
    observed = on_lattice(observed),
    expected = on_lattice(expected),
    .id = "surface"
  )
  observed_grid <- surface[surface$surface == "observed", ]
  argmax <- observed_grid[which.max(observed_grid$value), c("x", "y")]
  list(
    surface = tibble::tibble(ID = id, surface),
    markers = tibble::tibble(
      ID = id,
      surface = c("observed", "observed", "expected", "expected"),
      marker = c("truth", "best", "truth", "best"),
      x = c(truth_xy[[1L]], argmax$x, truth_xy[[1L]], mean(estimate_xy[, 1L])),
      y = c(truth_xy[[2L]], argmax$y, truth_xy[[2L]], mean(estimate_xy[, 2L]))
    )
  )
}

.data_loglik_surfaces <- function(benchmark, layout, figure, n_grid) {
  ids <- layout$ID[layout$corner]
  pieces <- lapply(ids, function(id) {
    .loglik_lattice(benchmark, id, figure, n_grid)
  })
  bind <- function(name) {
    .with_layout(dplyr::bind_rows(lapply(pieces, `[[`, name)), layout)
  }
  balance <- .figure_balances(figure)
  list(
    surface = bind("surface"),
    markers = bind("markers"),
    axis_labels = if (identical(figure, "fig02")) {
      c("p1 (unnormalised weight)", "p2 (unnormalised weight)")
    } else {
      balance$labels[1:2]
    },
    simplex = identical(figure, "fig02"),
    zoomed = identical(figure, "fig05"),
    held = if (identical(figure, "fig05")) balance$labels[[3L]]
  )
}

plot_loglik_surfaces <- function(data, surface) {
  expected <- identical(surface, "expected")
  legend <- .marker_legend(c(
    "true p*",
    if (expected) {
      "Monte Carlo mean of p-hat (Newton-Raphson)"
    } else {
      "maximum on the grid"
    }
  ))
  caption <- paste0(
    if (expected) {
      "Expected log-likelihood at the true composition, maximised at p* exactly."
    } else {
      "Log-likelihood of the first simulated bulk."
    },
    " Fill: log10(L / max L), floored at ",
    .loglik_floor,
    ".",
    if (!is.null(data$held)) paste0(" ", data$held, " held at its true value."),
    if (isTRUE(data$zoomed)) {
      "\nEach panel zooms on the region above the floor (padded by 10%): the peak narrows as the number of cells grows."
    }
  )
  lapply(levels(droplevels(data$surface$page)), function(page) {
    grid <- data$surface[
      data$surface$page == page & data$surface$surface == surface,
    ]
    marks <- data$markers[
      data$markers$page == page & data$markers$surface == surface,
    ]
    p <- ggplot2::ggplot(grid, ggplot2::aes(.data$x, .data$y)) +
      ggplot2::geom_raster(
        ggplot2::aes(fill = .data$value),
        interpolate = TRUE
      ) +
      ggplot2::scale_fill_viridis_c(
        name = "log10(L / max L)",
        limits = c(.loglik_floor, 0),
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(barwidth = grid::unit(6, "cm"))
      )
    if (data$simplex) {
      p <- p +
        ggplot2::geom_abline(
          intercept = 1,
          slope = -1,
          linetype = 2,
          colour = "grey85"
        )
    }
    p <- p +
      ggplot2::geom_point(
        data = marks[marks$marker == "best", ],
        shape = 21,
        fill = "#E41A1C",
        colour = "black",
        size = 3.2
      ) +
      ggplot2::geom_point(
        data = marks[marks$marker == "truth", ],
        shape = 23,
        fill = "white",
        colour = "black",
        size = 3.4
      ) +
      .not_simulated(grid, at_origin = isTRUE(data$zoomed))
    p <- if (isTRUE(data$zoomed)) {
      p +
        .panel_facets(grid, scales = "free") +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::theme(aspect.ratio = 1)
    } else {
      p +
        .panel_facets(grid) +
        ggplot2::coord_equal(
          xlim = range(grid$x),
          ylim = range(grid$y),
          expand = FALSE
        )
    }
    p +
      ggplot2::labs(
        x = data$axis_labels[[1L]],
        y = data$axis_labels[[2L]],
        title = if (expected) {
          "Expected log-likelihood"
        } else {
          "Observed log-likelihood"
        },
        subtitle = page,
        caption = caption
      ) +
      .page_theme() +
      ggplot2::theme(legend.position = "bottom")
  }) |>
    lapply(.with_legend, legend = legend)
}

# Log-likelihood along the single ILR coordinate of fig02, on the simplex.
.data_loglik_profiles <- function(benchmark, layout, n_grid = 300L) {
  axis <- seq(-5, 5, length.out = n_grid)
  ids <- layout$ID[layout$corner]
  pieces <- lapply(ids, function(id) {
    th <- benchmark$theta_true[[id]]
    p_true <- as.numeric(th$p)
    y <- benchmark$data[[id]]$bulk[, 1L]
    values <- t(vapply(
      axis,
      function(z) {
        p <- isometric_logistic(z)
        c(
          loglik_multivariate(p, y, th$mu, th$sigma),
          .expected_loglik_population(p, p_true, th$mu, th$sigma)
        )
      },
      numeric(2)
    ))
    relative <- function(v) pmax((v - max(v)) / log(10), .loglik_floor)
    est <- benchmark$optimisation[
      benchmark$optimisation$ID == id &
        benchmark$optimisation$algorithm == "Newton-Raphson",
      colnames(th$mu)
    ]
    est <- as.matrix(est)
    est <- est[apply(est > 0, 1L, all), , drop = FALSE]
    z_hat <- (log(est[, 1L]) - log(est[, 2L])) / sqrt(2)
    observed <- relative(values[, 1L])
    list(
      curve = tibble::tibble(
        ID = id,
        z = axis,
        observed = observed,
        expected = relative(values[, 2L])
      ),
      markers = tibble::tibble(
        ID = id,
        surface = c("observed", "observed", "expected", "expected"),
        marker = c("truth", "best", "truth", "best"),
        z = c(
          isometric_log_ratio(p_true),
          axis[[which.max(observed)]],
          isometric_log_ratio(p_true),
          mean(z_hat)
        )
      )
    )
  })
  bind <- function(name) {
    .with_layout(dplyr::bind_rows(lapply(pieces, `[[`, name)), layout)
  }
  list(curve = bind("curve"), markers = bind("markers"))
}

plot_loglik_profiles <- function(data, surface) {
  expected <- identical(surface, "expected")
  legend <- .marker_legend(c(
    "true p*",
    if (expected) "Monte Carlo mean of p-hat (Newton-Raphson)" else "maximum"
  ))
  lapply(levels(droplevels(data$curve$page)), function(page) {
    curve <- data$curve[data$curve$page == page, ]
    curve$value <- curve[[surface]]
    marks <- data$markers[
      data$markers$page == page & data$markers$surface == surface,
    ]
    marks <- dplyr::left_join(
      marks,
      curve[c("ID", "z", "value")] |> dplyr::rename(grid_z = "z"),
      by = "ID",
      relationship = "many-to-many"
    ) |>
      dplyr::group_by(.data$panel, .data$marker) |>
      dplyr::slice_min(
        abs(.data$grid_z - .data$z),
        n = 1L,
        with_ties = FALSE
      ) |>
      dplyr::ungroup()
    ggplot2::ggplot(curve, ggplot2::aes(.data$z, .data$value)) +
      ggplot2::geom_line(colour = "#1D4E89", linewidth = 0.8) +
      ggplot2::geom_point(
        data = marks[marks$marker == "best", ],
        shape = 21,
        fill = "#E41A1C",
        colour = "black",
        size = 3.2
      ) +
      ggplot2::geom_point(
        data = marks[marks$marker == "truth", ],
        shape = 23,
        fill = "white",
        colour = "black",
        size = 3.4
      ) +
      ggplot2::facet_wrap(~panel, ncol = 2L) +
      ggplot2::coord_cartesian(ylim = c(.loglik_floor, 0)) +
      ggplot2::labs(
        x = "z = balance(CT1 | CT2) = (log p1 - log p2) / sqrt(2)",
        y = "log10(L / max L)",
        title = if (expected) {
          "Expected log-likelihood along the simplex"
        } else {
          "Observed log-likelihood along the simplex"
        },
        subtitle = page,
        caption = paste(
          if (expected) {
            "Expected log-likelihood at p*."
          } else {
            "First simulated bulk."
          },
          "Same axes on every panel and page."
        )
      ) +
      .page_theme()
  }) |>
    lapply(.with_legend, legend = legend)
}

# Three-dimensional relative likelihood (rgl) of the observed lattices: one
# snapshot per panel, a 2 x 2 page and the marker legend.
plot_loglik_rgl <- function(data) {
  legend <- .marker_legend(c("true p*", "maximum on the grid"))
  snapshot <- function(grid, marks, title) {
    axis <- sort(unique(grid$x))
    axis_y <- sort(unique(grid$y))
    z <- matrix(grid$value, nrow = length(axis))
    rgl::open3d(windowRect = c(40L, 40L, 1040L, 840L))
    on.exit(rgl::close3d(), add = TRUE)
    rgl::persp3d(
      axis,
      axis_y,
      z,
      col = "steelblue",
      alpha = 0.75,
      xlab = "",
      ylab = "",
      zlab = "",
      zlim = c(.loglik_floor, 0),
      axes = FALSE,
      box = TRUE
    )
    rgl::aspect3d(1, 1, 0.6)
    rgl::view3d(theta = 35, phi = 30, zoom = 0.9)
    rgl::axes3d(cex = 0.8, nticks = 4L)
    rgl::title3d(
      xlab = data$axis_labels[[1L]],
      ylab = data$axis_labels[[2L]],
      zlab = "log10(L / max L)",
      cex = 0.8
    )
    lift <- function(x, y) {
      z[which.min(abs(axis - x)), which.min(abs(axis_y - y))] + 0.6
    }
    best <- marks[marks$marker == "best", ]
    truth <- marks[marks$marker == "truth", ]
    rgl::points3d(
      best$x,
      best$y,
      lift(best$x, best$y),
      color = "#E41A1C",
      size = 18
    )
    rgl::points3d(
      truth$x,
      truth$y,
      lift(truth$x, truth$y),
      color = "white",
      size = 18
    )
    file <- tempfile(fileext = ".png")
    on.exit(unlink(file), add = TRUE)
    rgl::snapshot3d(file, webshot = FALSE)
    image <- png::readPNG(file)
    ggplot2::ggplot() +
      ggplot2::annotation_raster(image, 0, 1, 0, 1) +
      ggplot2::coord_fixed(
        ratio = dim(image)[[1L]] / dim(image)[[2L]],
        xlim = c(0, 1),
        ylim = c(0, 1),
        expand = FALSE
      ) +
      ggplot2::labs(title = title) +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
  }
  lapply(levels(droplevels(data$surface$page)), function(page) {
    grid <- data$surface[
      data$surface$page == page & data$surface$surface == "observed",
    ]
    marks <- data$markers[
      data$markers$page == page & data$markers$surface == "observed",
    ]
    panels <- lapply(levels(droplevels(grid$panel)), function(panel) {
      snapshot(
        grid[grid$panel == panel, ],
        marks[marks$panel == panel, ],
        panel
      )
    })
    body <- cowplot::plot_grid(plotlist = panels, ncol = 2L)
    header <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste("Observed relative likelihood /", page),
        fontface = "bold",
        x = 0.01,
        hjust = 0
      )
    .with_legend(
      cowplot::plot_grid(header, body, ncol = 1L, rel_heights = c(0.05, 1)),
      legend
    )
  })
}

# Mean signature, networks and pairwise distances -----------------------------

# One mean matrix per distinct reference: fig03 has one, fig05 one per
# gene count and cell-size vector.
.reference_ids <- function(layout, figure) {
  keys <- switch(
    figure,
    fig03 = c("graph_ct1", "graph_ct2"),
    fig05 = c("n_genes", "cell_size_label")
  )
  keep <- layout$proportions == "balanced"
  if (identical(figure, "fig05")) {
    keep <- keep &
      layout$platform_scale == max(layout$platform_scale) &
      layout$n_cells == max(layout$n_cells)
  }
  candidates <- layout[keep, ]
  candidates[!duplicated(candidates[c(keys, "overlap_label")]), ]
}

.data_mean_signature <- function(benchmark, layout, figure) {
  refs <- .reference_ids(layout, figure)
  refs <- if (identical(figure, "fig05")) {
    refs[!duplicated(refs[c("n_genes", "cell_size_label")]), ]
  } else {
    refs[1L, ]
  }
  purrr::map_dfr(seq_len(nrow(refs)), function(i) {
    mu <- benchmark$theta_true[[refs$ID[[i]]]]$mu
    tibble::tibble(
      reference = if (identical(figure, "fig05")) {
        sprintf(
          "G = %d, %s cell sizes",
          refs$n_genes[[i]],
          refs$cell_size_label[[i]]
        )
      } else {
        sprintf("G = %d", nrow(mu))
      },
      gene = factor(rep(rownames(mu), ncol(mu)), rev(rownames(mu))),
      cell_type = factor(rep(colnames(mu), each = nrow(mu)), colnames(mu)),
      value = as.vector(mu)
    )
  })
}

plot_mean_signature <- function(data) {
  lapply(unique(data$reference), function(ref) {
    df <- data[data$reference == ref, ]
    many <- dplyr::n_distinct(df$gene) > 40L
    ggplot2::ggplot(
      df,
      ggplot2::aes(.data$cell_type, .data$gene, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = if (many) NA else "white", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(name = "Mean expression") +
      ggplot2::labs(
        x = NULL,
        y = "Gene",
        title = "Mean signature mu",
        subtitle = ref,
        caption = "Raw means of the purified types, no row or column scaling."
      ) +
      .page_theme() +
      ggplot2::theme(
        axis.text.y = if (many) {
          ggplot2::element_blank()
        } else {
          ggplot2::element_text(size = 8)
        },
        panel.grid = ggplot2::element_blank()
      )
  })
}

# Edges of every network with the covariance they carry (width |Sigma_ij|,
# colour its sign) and its genes: the block of a stochastic block model,
# the hubs of a scale-free graph. fig03: one page per overlap target, the
# five distinct graphs (types 1 and 2 under both families, type 3
# scale-free); fig05: one page per gene count and overlap target, the four
# types. The family of each graph is `theta_true$graph_models` when stored,
# otherwise the design of the figure; SBM blocks are
# `theta_true$blocks[[j]]` when stored, otherwise contiguous thirds.
.network_families <- function(th, ref, figure) {
  if (!is.null(th$graph_models)) {
    return(unname(th$graph_models))
  }
  if (identical(figure, "fig05")) {
    return(c(
      "stochastic_block_model",
      "stochastic_block_model",
      "scale_free",
      "scale_free"
    ))
  }
  vapply(
    seq_along(th$adjacency),
    function(j) ref[[paste0("graph_ct", j)]][[1L]],
    character(1)
  )
}

.network_node_roles <- function(adj, family, blocks) {
  n_genes <- nrow(adj)
  degree <- rowSums(adj != 0)
  if (identical(family, "stochastic_block_model")) {
    if (is.null(blocks)) {
      blocks <- cut(seq_len(n_genes), 3L, labels = FALSE)
    }
    return(tibble::tibble(
      node = seq_len(n_genes),
      role = paste("block", blocks),
      hub = NA_integer_
    ))
  }
  threshold <- max(stats::quantile(degree, 0.9), 2)
  hubs <- order(degree, decreasing = TRUE)[seq_len(min(3L, n_genes))]
  hubs <- hubs[degree[hubs] >= threshold]
  tibble::tibble(
    node = seq_len(n_genes),
    role = ifelse(degree >= threshold, "hub", "other gene"),
    hub = vapply(
      seq_len(n_genes),
      function(g) {
        owner <- hubs[g == hubs | adj[g, hubs] != 0]
        if (length(owner) == 0L) NA_integer_ else owner[[1L]]
      },
      integer(1)
    )
  )
}

.data_network_topologies <- function(benchmark, layout, figure) {
  refs <- .reference_ids(layout, figure)
  if (identical(figure, "fig05")) {
    refs <- refs[refs$cell_size_label == "equal", ]
  }
  family_label <- c(scale_free = "scale-free", stochastic_block_model = "SBM")
  pieces <- purrr::map(seq_len(nrow(refs)), function(i) {
    th <- benchmark$theta_true[[refs$ID[[i]]]]
    families <- .network_families(th, refs[i, ], figure)
    page <- if (identical(figure, "fig05")) {
      sprintf(
        "G = %d / %s sibling overlap",
        refs$n_genes[[i]],
        refs$overlap_label[[i]]
      )
    } else {
      paste(refs$overlap_label[[i]], "overlap")
    }
    purrr::map(seq_along(th$adjacency), function(j) {
      graph <- if (identical(figure, "fig05")) {
        sprintf(
          "Cell type %d (lineage %d, %s)",
          j,
          (j + 1L) %/% 2L,
          family_label[[families[[j]]]]
        )
      } else {
        sprintf("Cell type %d, %s", j, family_label[[families[[j]]]])
      }
      adj <- th$adjacency[[j]]
      edges <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
      list(
        edges = tibble::tibble(
          page = page,
          graph = graph,
          n_genes = nrow(adj),
          from = edges[, 1L],
          to = edges[, 2L],
          covariance = th$sigma[,, j][edges]
        ),
        nodes = tibble::tibble(
          page = page,
          graph = graph,
          .network_node_roles(adj, families[[j]], th$blocks[[j]])
        )
      )
    })
  }) |>
    purrr::list_flatten()
  list(
    edges = dplyr::distinct(
      dplyr::bind_rows(lapply(pieces, `[[`, "edges")),
      .data$page,
      .data$graph,
      .data$from,
      .data$to,
      .keep_all = TRUE
    ),
    nodes = dplyr::distinct(
      dplyr::bind_rows(lapply(pieces, `[[`, "nodes")),
      .data$page,
      .data$graph,
      .data$node,
      .keep_all = TRUE
    )
  )
}

.node_role_palette <- c(
  "block 1" = "#1B9E77",
  "block 2" = "#7570B3",
  "block 3" = "#E7298A",
  "hub" = "#E6AB02",
  "other gene" = "grey35"
)

plot_network_topologies <- function(data) {
  width_max <- max(abs(data$edges$covariance))
  panel <- function(edges, nodes) {
    graph <- igraph::graph_from_data_frame(
      data.frame(
        from = edges$from,
        to = edges$to,
        covariance = edges$covariance
      ),
      directed = FALSE,
      vertices = data.frame(
        name = as.character(nodes$node),
        role = factor(nodes$role, names(.node_role_palette)),
        hub = nodes$hub
      )
    )
    layout <- withr::with_seed(
      1L,
      ggraph::create_layout(graph, layout = "fr")
    )
    neighbourhoods <- layout[!is.na(layout$hub), ]
    neighbourhoods$hub <- factor(neighbourhoods$hub)
    many <- nrow(nodes) > 50L
    ggraph::ggraph(layout) +
      ggforce::geom_mark_ellipse(
        data = neighbourhoods,
        ggplot2::aes(.data$x, .data$y, group = .data$hub),
        fill = "#E6AB02",
        alpha = 0.12,
        colour = NA,
        expand = grid::unit(2, "mm")
      ) +
      ggraph::geom_edge_link(ggplot2::aes(
        edge_width = abs(.data$covariance),
        edge_colour = ifelse(.data$covariance > 0, "positive", "negative")
      )) +
      ggraph::geom_node_point(
        ggplot2::aes(fill = .data$role),
        shape = 21,
        colour = "black",
        stroke = if (many) 0.15 else 0.4,
        size = if (many) 1 else 2.6,
        # Keys of roles absent from this panel (hubs in a block model,
        # blocks in a scale-free graph) keep their point.
        show.legend = TRUE
      ) +
      ggraph::scale_edge_width(
        name = latex2exp::TeX("$|\\Sigma_{ij}|$"),
        range = c(0.2, 2.4),
        limits = c(0, width_max)
      ) +
      ggraph::scale_edge_colour_manual(
        name = latex2exp::TeX("Sign of $\\Sigma_{ij}$"),
        values = c(positive = "#2B6CB0", negative = "#C53030"),
        drop = FALSE
      ) +
      ggplot2::scale_fill_manual(
        name = "Gene",
        values = .node_role_palette,
        drop = FALSE
      ) +
      ggplot2::labs(title = edges$graph[[1L]]) +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 10,
          hjust = 0.5
        ),
        legend.position = "bottom",
        legend.box = "vertical"
      )
  }
  lapply(unique(data$edges$page), function(page) {
    edges <- data$edges[data$edges$page == page, ]
    nodes <- data$nodes[data$nodes$page == page, ]
    graphs <- sort(unique(nodes$graph))
    panels <- lapply(graphs, function(g) {
      panel(edges[edges$graph == g, ], nodes[nodes$graph == g, ])
    })
    .shared_legend_page(
      panels,
      title = paste("Gene networks and the covariance on their edges /", page),
      caption = paste(
        "Edges of the precision support; width |Σij| on a scale shared by every page, colour its sign.",
        "Genes coloured by block in the stochastic block models;",
        "\nhubs (degree at or above the 90th percentile) in yellow in the scale-free graphs, with a shaded",
        "ellipse around the three largest hubs and their neighbours."
      ),
      ncol = if (length(panels) == 4L) 2L else 3L,
      legend_height = 0.2
    )
  })
}

# Pairwise discrepancies of the purified laws: Hellinger distance, MixSim
# overlap of the pair and affine-invariant Riemannian distance of the
# covariances (none depends on the composition).
.data_pairwise_distances <- function(benchmark, layout, figure) {
  refs <- .reference_ids(layout, figure)
  purrr::map_dfr(seq_len(nrow(refs)), function(i) {
    th <- benchmark$theta_true[[refs$ID[[i]]]]
    pairs <- utils::combn(ncol(th$mu), 2L)
    purrr::map_dfr(seq_len(ncol(pairs)), function(k) {
      a <- pairs[1L, k]
      b <- pairs[2L, k]
      moments <- list(mu = th$mu[, c(a, b)], sigma = th$sigma[,, c(a, b)])
      overlap <- overlap_gaussian_mc(
        true_theta = c(moments, list(p = c(0.5, 0.5))),
        n_mc = 2000L,
        seed = 20261024L + k
      )$BarOmega
      tibble::tibble(
        ID = refs$ID[[i]],
        pair = sprintf("CT%d-CT%d", a, b),
        hellinger = hellinger_gaussian(
          th$mu[, a],
          th$mu[, b],
          th$sigma[,, a],
          th$sigma[,, b]
        ),
        overlap = overlap,
        riemannian = compute_average_riemannian(moments)
      )
    })
  }) |>
    .with_layout(refs) |>
    dplyr::left_join(
      refs[c(
        "ID",
        "overlap_label",
        intersect(c("n_genes", "cell_size_label"), names(refs))
      )],
      by = "ID"
    )
}

plot_pairwise_distances <- function(data, figure) {
  specs <- list(
    hellinger = "Hellinger distance",
    overlap = "MixSim overlap of the pair",
    riemannian = "Affine-invariant Riemannian distance of Sigma_a and Sigma_b"
  )
  if (identical(figure, "fig05")) {
    data$tile_x <- factor(
      paste("G =", data$n_genes),
      paste("G =", sort(unique(data$n_genes)))
    )
    data$tile_y <- factor(data$cell_size_label)
    axes <- c("Genes", "Cell sizes")
  } else {
    data$tile_x <- data$panel_x
    data$tile_y <- data$panel_y
    axes <- c(data$panel_x_title[[1L]], data$panel_y_title[[1L]])
  }
  lapply(names(specs), function(metric) {
    data$value <- data[[metric]]
    ggplot2::ggplot(
      data,
      ggplot2::aes(.data$tile_x, .data$tile_y, fill = .data$value)
    ) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(
        ggplot2::aes(label = .signif_label(.data$value)),
        colour = "white",
        size = 3.2
      ) +
      ggplot2::facet_grid(pair ~ overlap_label) +
      ggplot2::scale_fill_viridis_c(name = specs[[metric]], end = 0.9) +
      ggplot2::labs(
        x = axes[[1L]],
        y = axes[[2L]],
        title = paste("Pairwise", tolower(specs[[metric]])),
        subtitle = "Rows: cell-type pairs; columns: overlap target",
        caption = "Purified laws only; none of these quantities depends on the composition."
      ) +
      .page_theme() +
      ggplot2::theme(legend.key.width = grid::unit(2, "cm"))
  })
}

# Latent projections ------------------------------------------------------------

# Labelled purified draws (equal counts per type) projected on the first two
# principal components, which only the means drive, and on the first two
# directions of MclustDR, which also use the covariance differences.
.data_latent_projections <- function(
  benchmark,
  layout,
  figure,
  n_per_type = 400L
) {
  refs <- .reference_ids(layout, figure)
  scree <- list()
  projections <- purrr::map_dfr(seq_len(nrow(refs)), function(i) {
    th <- benchmark$theta_true[[refs$ID[[i]]]]
    J <- ncol(th$mu)
    draws <- withr::with_seed(20261024L + i, {
      do.call(
        rbind,
        lapply(seq_len(J), function(j) {
          MASS::mvrnorm(n_per_type, th$mu[, j], th$sigma[,, j])
        })
      )
    })
    class <- factor(rep(paste("Cell type", seq_len(J)), each = n_per_type))
    pca <- stats::prcomp(draws, center = TRUE, scale. = FALSE)
    share <- pca$sdev^2 / sum(pca$sdev^2)
    # MclustDA evaluates its own calls in the caller's frame, so mclust has
    # to be attached while it runs.
    dr <- withr::with_package("mclust", {
      da <- mclust::MclustDA(
        draws,
        class,
        modelType = "EDDA",
        modelNames = "VVV",
        verbose = FALSE
      )
      mclust::MclustDR(da, lambda = 1)
    })
    page <- switch(
      figure,
      fig03 = sprintf(
        "%s overlap / %s",
        refs$overlap_label[[i]],
        refs$panel[[i]]
      ),
      fig05 = sprintf(
        "G = %d / %s sibling overlap / %s cell sizes",
        refs$n_genes[[i]],
        refs$overlap_label[[i]],
        refs$cell_size_label[[i]]
      )
    )
    # Share of the total variance of the first ten components.
    scree[[i]] <<- tibble::tibble(
      page = page,
      component = seq_len(min(10L, length(share))),
      share = share[seq_len(min(10L, length(share)))]
    )
    dplyr::bind_rows(
      tibble::tibble(
        page = page,
        method = sprintf(
          "PCA (%.0f%% + %.0f%% of the variance)",
          100 * share[[1L]],
          100 * share[[2L]]
        ),
        cell_type = class,
        dim_1 = pca$x[, 1L],
        dim_2 = pca$x[, 2L]
      ),
      tibble::tibble(
        page = page,
        method = "MclustDR (means and covariances)",
        cell_type = class,
        dim_1 = dr$dir[, 1L],
        dim_2 = dr$dir[, 2L]
      )
    )
  })
  attr(projections, "scree") <- dplyr::bind_rows(scree)
  projections
}

# Scree plot of one page: share of the total variance of the first ten
# principal components, drawn on [0, 1] so the bars sum to at most one.
.scree_inset <- function(scree) {
  ggplot2::ggplot(scree, ggplot2::aes(.data$component, .data$share)) +
    ggplot2::geom_col(fill = "grey55", colour = "grey25", linewidth = 0.1) +
    ggplot2::scale_x_continuous(breaks = c(1, 5, 10)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = "Component", y = "Share of variance", title = "Scree") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 7),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = "grey70",
        linewidth = 0.2
      )
    )
}

plot_latent_projections <- function(data) {
  scree <- attr(data, "scree")
  lapply(unique(data$page), function(page) {
    df <- data[data$page == page, ]
    methods <- unique(as.character(df$method))
    df$method <- factor(
      df$method,
      c(
        grep("^PCA", methods, value = TRUE),
        grep("^PCA", methods, value = TRUE, invert = TRUE)
      )
    )
    # Head-room above the PCA panel, where the scree inset sits.
    pca <- df[startsWith(as.character(df$method), "PCA"), ]
    room <- tibble::tibble(
      method = factor(pca$method[[1L]], levels(df$method)),
      dim_1 = mean(pca$dim_1),
      dim_2 = max(pca$dim_2) + 1.1 * diff(range(pca$dim_2))
    )
    main <- ggplot2::ggplot(
      df,
      ggplot2::aes(.data$dim_1, .data$dim_2)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data$cell_type),
        shape = 21,
        size = 0.9,
        alpha = 0.45,
        colour = "grey60",
        stroke = 0.1
      ) +
      ggplot2::stat_ellipse(
        ggplot2::aes(colour = .data$cell_type),
        level = 0.95,
        linewidth = 0.7
      ) +
      ggplot2::geom_blank(data = room) +
      ggplot2::facet_wrap(~method, scales = "free") +
      ggplot2::scale_colour_brewer(name = NULL, palette = "Dark2") +
      ggplot2::scale_fill_brewer(name = NULL, palette = "Dark2") +
      ggplot2::guides(
        colour = "none",
        fill = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
      ) +
      ggplot2::labs(
        x = "Direction 1",
        y = "Direction 2",
        title = "Purified profiles in two latent directions",
        subtitle = page,
        caption = paste0(
          "Equal numbers of draws per type; ellipses: 95% normal regions.\n",
          "PCA uses the pooled variance, MclustDR (EDDA, VVV) the class means and covariances.\n",
          "Inset: share of the total variance of each of the first ten principal components."
        )
      ) +
      .page_theme()
    if (is.null(scree)) {
      return(main)
    }
    # cowplot inset (draw_plot), over the head-room of the left panel
    cowplot::ggdraw(main) +
      cowplot::draw_plot(
        .scree_inset(scree[scree$page == page, ]),
        x = 0.085,
        y = 0.55,
        width = 0.2,
        height = 0.24
      )
  })
}
