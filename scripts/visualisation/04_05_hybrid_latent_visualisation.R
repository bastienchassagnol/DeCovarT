#' Cell-type colours for the G = 20, J = 3 latent-space book
#'
#' @keywords internal
#' @noRd
.hybrid_celltype_colours <- function(cts) {
  cts <- as.character(cts)
  pal <- c("#2B8C99", "#F4B400", "#C1443C", "#6A51A3")
  stats::setNames(rep(pal, length.out = length(cts)), cts)
}

#' First k eigenvalues as percentages of the trace
#'
#' @keywords internal
#' @noRd
.scree_percent <- function(values, k = 10L) {
  values <- as.numeric(values)
  values <- values[is.finite(values) & values >= 0]
  if (length(values) == 0L) {
    return(numeric(0))
  }
  tot <- sum(values)
  if (!is.finite(tot) || tot <= 0) {
    return(numeric(0))
  }
  k <- min(as.integer(k), length(values))
  100 * values[seq_len(k)] / tot
}

#' Scree inset (first 10 scaled contributions)
#'
#' @keywords internal
#' @noRd
.plot_scree_inset <- function(pct) {
  pct <- as.numeric(pct)
  if (length(pct) == 0L) {
    pct <- 0
  }
  df <- data.frame(
    comp = factor(seq_along(pct), levels = seq_along(pct)),
    pct = pct
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["comp"]], y = .data[["pct"]])) +
    ggplot2::geom_col(fill = "grey35", width = 0.72) +
    ggplot2::geom_line(
      ggplot2::aes(group = 1L),
      colour = "grey20",
      linewidth = 0.35
    ) +
    ggplot2::geom_point(size = 0.7, colour = "grey20") +
    ggplot2::labs(x = NULL, y = "%") +
    ggplot2::theme_classic(base_size = 7) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      axis.text = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(2, 4, 2, 2)
    )
}

#' Thomson / regression factor scores, with PCA fallback
#'
#' @keywords internal
#' @noRd
.thomson_projection <- function(X, n_factors = 2L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n_factors <- as.integer(n_factors)
  pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  ev <- pca$sdev^2
  scree <- .scree_percent(ev, k = 10L)
  fa <- tryCatch(
    stats::factanal(
      X,
      factors = n_factors,
      rotation = "none",
      scores = "regression"
    ),
    error = function(e) NULL
  )
  if (is.null(fa) || is.null(fa$scores)) {
    scores <- pca$x[, seq_len(n_factors), drop = FALSE]
    axis_pct <- .scree_percent(ev, k = n_factors)
  } else {
    scores <- as.matrix(fa$scores)
    ss <- colSums(unclass(fa$loadings)^2)
    axis_pct <- 100 * ss / ncol(X)
  }
  colnames(scores) <- paste0("Dim", seq_len(ncol(scores)))
  list(scores = scores, axis_pct = axis_pct, scree = scree)
}

#' Shared MCFA scores (`EMMIXmfa::mcfa`)
#'
#' @keywords internal
#' @noRd
.mcfa_projection <- function(
  Y,
  g = 3L,
  q = 2L,
  itmax = 150L,
  Y_score = NULL
) {
  .check_suggested_package("EMMIXmfa", ".mcfa_projection")
  Y <- as.matrix(Y)
  if (is.null(Y_score)) {
    Y_score <- Y
  } else {
    Y_score <- as.matrix(Y_score)
  }
  fit <- NULL
  utils::capture.output({
    fit <- tryCatch(
      EMMIXmfa::mcfa(
        Y,
        g = as.integer(g),
        q = as.integer(q),
        itmax = as.integer(itmax),
        nkmeans = 2L,
        nrandom = 2L,
        warn_messages = FALSE
      ),
      error = function(e) NULL
    )
  })
  if (is.null(fit)) {
    pca <- stats::prcomp(Y_score, center = TRUE, scale. = TRUE)
    scores <- pca$x[, seq_len(q), drop = FALSE]
    colnames(scores) <- paste0("Dim", seq_len(ncol(scores)))
    ev <- pca$sdev^2
    return(
      list(
        scores = scores,
        axis_pct = .scree_percent(ev, k = q),
        scree = .scree_percent(ev, k = 10L)
      )
    )
  }
  sc <- EMMIXmfa::factor_scores(fit, Y_score)
  scores <- as.matrix(sc$Umean)
  if (ncol(scores) > q) {
    scores <- scores[, seq_len(q), drop = FALSE]
  }
  colnames(scores) <- paste0("Dim", seq_len(ncol(scores)))
  omega_bar <- matrix(0, q, q)
  pi_hat <- as.numeric(fit$pivec)
  pi_hat <- pi_hat / sum(pi_hat)
  for (i in seq_len(fit$g)) {
    omega_bar <- omega_bar + pi_hat[[i]] * fit$omega[,, i]
  }
  sigma_hat <- fit$A %*% omega_bar %*% t(fit$A) + fit$D
  ev <- eigen(sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(as.numeric(ev), 0)
  list(
    scores = scores,
    axis_pct = .scree_percent(ev, k = q),
    scree = .scree_percent(ev, k = 10L)
  )
}

#' Supervised `MclustDR` scores
#'
#' @keywords internal
#' @noRd
.mclustdr_projection <- function(Y, class, q = 2L) {
  .check_suggested_package("mclust", ".mclustdr_projection")
  Y <- as.matrix(Y)
  class <- factor(class)
  da <- tryCatch(
    mclust::MclustDA(
      Y,
      class,
      modelType = "EDDA",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  dr <- NULL
  if (!is.null(da)) {
    dr <- tryCatch(
      mclust::MclustDR(da, lambda = 1),
      error = function(e) NULL
    )
  }
  if (is.null(dr) || is.null(dr$dir)) {
    pca <- stats::prcomp(Y, center = TRUE, scale. = TRUE)
    scores <- pca$x[, seq_len(q), drop = FALSE]
    colnames(scores) <- paste0("Dim", seq_len(ncol(scores)))
    ev <- pca$sdev^2
    return(
      list(
        scores = scores,
        axis_pct = .scree_percent(ev, k = q),
        scree = .scree_percent(ev, k = 10L)
      )
    )
  }
  nd <- ncol(dr$dir)
  take <- seq_len(min(q, nd))
  scores <- dr$dir[, take, drop = FALSE]
  if (ncol(scores) < q) {
    pad <- matrix(0, nrow = nrow(scores), ncol = q - ncol(scores))
    scores <- cbind(scores, pad)
  }
  colnames(scores) <- paste0("Dim", seq_len(ncol(scores)))
  ev <- as.numeric(dr$evalues)
  list(
    scores = scores,
    axis_pct = .scree_percent(ev, k = q),
    scree = .scree_percent(ev, k = 10L)
  )
}

#' Labelled draws from a Gaussian mixture (weights \eqn{p}, not \eqn{p^2})
#'
#' @keywords internal
#' @noRd
.simulate_labelled_mixture <- function(mu, Sigma, p, n) {
  mu <- as.matrix(mu)
  p <- as.numeric(p)
  p <- p / sum(p)
  j <- ncol(mu)
  g <- nrow(mu)
  z <- sample.int(j, n, replace = TRUE, prob = p)
  y <- matrix(NA_real_, n, g)
  colnames(y) <- rownames(mu)
  for (k in seq_len(j)) {
    idx <- which(z == k)
    n_k <- length(idx)
    if (n_k == 0L) {
      next
    }
    draws <- MASS::mvrnorm(
      n = n_k,
      mu = mu[, k],
      Sigma = Sigma[,, k],
      empirical = FALSE
    )
    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1L)
    }
    y[idx, ] <- draws
  }
  list(Y = y, class = z)
}

#' 2D latent scatter with cell-type ellipses and a scree inset
#'
#' @keywords internal
#' @noRd
.plot_latent_scatter <- function(
  scores_df,
  axis_pct,
  scree,
  title,
  pal
) {
  .check_suggested_package("cowplot", ".plot_latent_scatter")
  scores_df <- as.data.frame(scores_df)
  pct1 <- if (length(axis_pct) >= 1L) axis_pct[[1L]] else NA_real_
  pct2 <- if (length(axis_pct) >= 2L) axis_pct[[2L]] else NA_real_
  xlab <- sprintf("Dim. 1 (%.1f%%)", pct1)
  ylab <- sprintf("Dim. 2 (%.1f%%)", pct2)
  n_grp <- length(unique(as.character(scores_df$cell_type)))
  p <- ggplot2::ggplot(
    scores_df,
    ggplot2::aes(x = .data[["Dim1"]], y = .data[["Dim2"]])
  )
  if (n_grp == 1L) {
    p <- p +
      ggplot2::geom_density_2d_filled(
        contour_var = "ndensity",
        bins = 8L,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_viridis_d(guide = "none")
  }
  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data[["cell_type"]]),
      size = 0.55,
      alpha = 0.45
    )
  n_by_type <- table(as.character(scores_df$cell_type))
  ellipse_ok <- names(n_by_type)[n_by_type >= 4L]
  if (length(ellipse_ok) > 0L) {
    df_ell <- scores_df[
      as.character(scores_df$cell_type) %in% ellipse_ok,
      ,
      drop = FALSE
    ]
    p <- p +
      ggplot2::stat_ellipse(
        data = df_ell,
        ggplot2::aes(colour = .data[["cell_type"]]),
        level = 0.95,
        linewidth = 0.7,
        type = "norm"
      )
  }
  p <- p +
    ggplot2::scale_colour_manual(
      values = pal,
      drop = FALSE,
      name = "Cell type"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    theme_decovart_facets() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7)
    )
  inset <- .plot_scree_inset(scree)
  cowplot::ggdraw(p) +
    cowplot::draw_plot(
      inset,
      x = 0.62,
      y = 0.58,
      width = 0.32,
      height = 0.30
    )
}

#' Fig03 pages: extreme overlap \eqn{\times} extreme entropy \eqn{\times}
#' four topologies
#'
#' @keywords internal
#' @noRd
.hybrid_latent_page_meta <- function(config) {
  config <- .relevel_scenario_table(config)
  ovl_code <- .canonical_code(
    as.character(config$overlap_label),
    .overlap_display_labels()
  )
  keep_ovl <- ovl_code %in% c("low", "high")
  keep_p <- as.character(config$proportions) %in%
    c("balanced", "highly unbalanced")
  out <- config[keep_ovl & keep_p, , drop = FALSE]
  dplyr::arrange(
    out,
    .data[["proportions"]],
    .data[["overlap_label"]],
    .data[["graph_ct1"]],
    .data[["graph_ct2"]]
  )
}

#' One 3-by-2 latent-space page for a hybrid scenario
#'
#' @keywords internal
#' @noRd
.hybrid_latent_page <- function(true_theta, n, itmax) {
  th <- .unwrap_true_theta(true_theta)
  mu <- as.matrix(th$mu)
  sigma <- th$sigma
  p <- as.numeric(th$p)
  p <- p / sum(p)
  cts <- colnames(mu)
  if (is.null(cts)) {
    cts <- paste0("celltype_", seq_len(ncol(mu)))
    colnames(mu) <- cts
  }
  pal <- .hybrid_celltype_colours(cts)
  j <- ncol(mu)
  purified <- lapply(seq_len(j), function(k) {
    draws <- MASS::mvrnorm(
      n = n,
      mu = mu[, k],
      Sigma = sigma[,, k],
      empirical = FALSE
    )
    if (is.null(dim(draws))) {
      draws <- matrix(draws, nrow = 1L)
    }
    draws
  })
  names(purified) <- cts
  top <- lapply(seq_len(j), function(k) {
    proj <- .thomson_projection(purified[[k]], n_factors = 2L)
    df <- data.frame(
      Dim1 = proj$scores[, 1L],
      Dim2 = proj$scores[, 2L],
      cell_type = factor(cts[[k]], levels = cts)
    )
    .plot_latent_scatter(
      df,
      proj$axis_pct,
      proj$scree,
      title = paste0(cts[[k]], " \u00b7 Thomson FA"),
      pal = pal
    )
  })
  while (length(top) < 3L) {
    top[[length(top) + 1L]] <- cowplot::ggdraw()
  }
  sim <- simulate_bulk_mixture(mu, sigma, p = p, n = n)
  y_conv <- t(sim$Y)
  latent_rows <- lapply(seq_len(j), function(k) {
    t(sim$latent_profiles[, k, ])
  })
  y_lat <- do.call(rbind, latent_rows)
  lab_lat <- factor(rep(cts, each = n), levels = cts)
  mcfa_conv <- .mcfa_projection(
    y_conv,
    g = j,
    q = 2L,
    itmax = itmax,
    Y_score = y_lat
  )
  df_conv <- data.frame(
    Dim1 = mcfa_conv$scores[, 1L],
    Dim2 = mcfa_conv$scores[, 2L],
    cell_type = lab_lat
  )
  p_conv <- .plot_latent_scatter(
    df_conv,
    mcfa_conv$axis_pct,
    mcfa_conv$scree,
    title = "MCFA \u00b7 convolution",
    pal = pal
  )
  mix <- .simulate_labelled_mixture(mu, sigma, p, n)
  lab_mix <- factor(cts[mix$class], levels = cts)
  mcfa_mix <- .mcfa_projection(
    mix$Y,
    g = j,
    q = 2L,
    itmax = itmax
  )
  df_mix <- data.frame(
    Dim1 = mcfa_mix$scores[, 1L],
    Dim2 = mcfa_mix$scores[, 2L],
    cell_type = lab_mix
  )
  p_mix <- .plot_latent_scatter(
    df_mix,
    mcfa_mix$axis_pct,
    mcfa_mix$scree,
    title = "MCFA \u00b7 unsupervised mixture",
    pal = pal
  )
  dr <- .mclustdr_projection(mix$Y, mix$class, q = 2L)
  df_dr <- data.frame(
    Dim1 = dr$scores[, 1L],
    Dim2 = dr$scores[, 2L],
    cell_type = lab_mix
  )
  p_dr <- .plot_latent_scatter(
    df_dr,
    dr$axis_pct,
    dr$scree,
    title = "MclustDR \u00b7 supervised mixture",
    pal = pal
  )
  cowplot::plot_grid(
    top[[1L]],
    top[[2L]],
    top[[3L]],
    p_conv,
    p_mix,
    p_dr,
    ncol = 3L,
    nrow = 2L,
    labels = c("A", "B", "C", "D", "E", "F"),
    label_size = 11,
    label_fontface = "bold",
    align = "none"
  )
}

#' Six-line footnote for one latent-projection page
#'
#' @keywords internal
#' @noRd
.hybrid_latent_caption <- function() {
  paste(
    "A: Independent Thomson factanal() on purified type 1 (regression scores).",
    "B: Independent Thomson factanal() on purified type 2 (regression scores).",
    "C: Independent Thomson factanal() on purified type 3 (regression scores).",
    "D: MCFA fitted to convolution bulk draws (weights p_j^2); purified slices scored in that plane.",
    "E: MCFA fitted to unsupervised mixture draws (one type per draw, weights p_j).",
    "F: Supervised MclustDA (EDDA) then MclustDR on labelled mixture draws. Insets: first 10 eigenvalues.",
    sep = "\n"
  )
}

#' 16-page 2D latent-space book for the covariance-driven scenario
#'
#' Extreme MixSim overlap (low / high) crossed with extreme Shannon
#' composition (balanced / highly unbalanced) and the four CT1/CT2
#' graph assignments (\eqn{2 \times 2 \times 4 = 16} pages). Each page
#' is a 3-by-2 of projections: independent Thomson factor analyses
#' \insertCite{thomsonMethodsEstimatingMental1938}{DeCovarT} on the
#' purified types (top), then a shared plane from mixtures of common
#' factor analysers
#' \insertCite{MixturesFactorAnalyzers2000}{DeCovarT} on convolution
#' draws, the same MCFA on an unsupervised mixture with weights
#' \eqn{\boldsymbol{p}}, and supervised [mclust::MclustDR()]
#' \insertCite{scruccaDimensionReductionModelbased2010a,scruccaGraphicalToolsModelbased2015}{DeCovarT}
#' (bottom).
#'
#' @param artefacts Split or assembled fig03 artefacts from
#'   [read_simulation_artefacts()].
#' @param file Output PDF path.
#' @param n Draws per cell type / bulk sample (default 300).
#' @param seed Optional seed; `NULL` leaves the RNG state unchanged.
#' @param itmax EM iterations for `EMMIXmfa::mcfa()`.
#' @param data_rds Optional directory for ggplot `data` RDS files.
#'
#' @return `file`, invisibly.
#' @export
#' @seealso [simulate_bulk_mixture()], [save_hybrid_runtime_book()]
#' @references
#' \insertAllCited{}
#' @examples
#' skip <- !requireNamespace("EMMIXmfa", quietly = TRUE) ||
#'   !requireNamespace("mclust", quietly = TRUE) ||
#'   !requireNamespace("cowplot", quietly = TRUE)
#' if (!skip) {
#'   genes <- paste0("g", seq_len(6L))
#'   cts <- paste0("celltype_", 1:3)
#'   mu <- matrix(
#'     c(8, 9, 9, 8, 2, 12, 2, 3, 11, 10, 4, 5, 5, 4, 12, 3, 6, 7),
#'     nrow = 6L,
#'     dimnames = list(genes, cts)
#'   )
#'   sig <- diag(6)
#'   sigma <- array(c(sig, sig, sig), dim = c(6L, 6L, 3L))
#'   dimnames(sigma) <- list(genes, genes, cts)
#'   th <- list(p = c(0.5, 0.3, 0.2), mu = mu, sigma = sigma)
#'   artefacts <- list(
#'     config = tibble::tibble(
#'       ID = "V1",
#'       proportions = "balanced",
#'       overlap_label = "low",
#'       graph_ct1 = "scale_free",
#'       graph_ct2 = "scale_free",
#'       graph_ct3 = "scale_free"
#'     ),
#'     theta = tibble::tibble(ID = "V1", true_theta = list(th))
#'   )
#'   tf <- withr::local_tempfile(fileext = ".pdf")
#'   save_hybrid_latent_projection_book(
#'     artefacts,
#'     tf,
#'     n = 40L,
#'     seed = 1L,
#'     itmax = 20L
#'   )
#' }
save_hybrid_latent_projection_book <- function(
  artefacts,
  file,
  n = 300L,
  seed = NULL,
  itmax = 150L,
  data_rds = NULL
) {
  .check_suggested_package("cowplot", "save_hybrid_latent_projection_book")
  .check_suggested_package("EMMIXmfa", "save_hybrid_latent_projection_book")
  .check_suggested_package("mclust", "save_hybrid_latent_projection_book")
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  cfg <- .hybrid_latent_page_meta(artefacts$config)
  theta_tbl <- artefacts$theta
  if (is.null(theta_tbl) || !"true_theta" %in% names(theta_tbl)) {
    stop(
      "`artefacts$theta` with a `true_theta` list-column is required.",
      call. = FALSE
    )
  }
  theta_tbl <- .relevel_scenario_table(theta_tbl)
  collected <- list()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 20, height = 14)
  on.exit(grDevices::dev.off(), add = TRUE)
  caption <- .hybrid_latent_caption()
  for (i in seq_len(nrow(cfg))) {
    row <- cfg[i, , drop = FALSE]
    id <- as.character(row$ID[[1L]])
    hit <- as.character(theta_tbl$ID) == id
    if (!any(hit)) {
      next
    }
    th <- theta_tbl$true_theta[which(hit)[[1L]]]
    page <- .hybrid_latent_page(th, n = as.integer(n), itmax = itmax)
    title <- paste(
      as.character(row$proportions[[1L]]),
      as.character(row$overlap_label[[1L]]),
      .graph_pair_short(row$graph_ct1[[1L]], row$graph_ct2[[1L]]),
      sep = " / "
    )
    header <- cowplot::ggdraw() +
      cowplot::draw_label(title, fontface = "bold", size = 14)
    foot <- cowplot::ggdraw() +
      cowplot::draw_label(
        caption,
        size = 8,
        hjust = 0,
        vjust = 1,
        x = 0.01,
        y = 0.98,
        lineheight = 1.15
      )
    print(
      cowplot::plot_grid(
        header,
        page,
        foot,
        ncol = 1L,
        rel_heights = c(0.05, 0.78, 0.17)
      ),
      newpage = i > 1L
    )
    collected[[length(collected) + 1L]] <- data.frame(
      ID = id,
      page = title,
      stringsAsFactors = FALSE
    )
  }
  if (length(collected) > 0L) {
    .write_ggplot_rds(
      dplyr::bind_rows(collected),
      data_rds,
      "latent_projections"
    )
  }
  invisible(file)
}

#' Topology x-axis labels (SF/SF, SF/SBM, ...)
#'
#' @keywords internal
#' @noRd
.hybrid_topology_levels <- function() {
  vapply(
    .hybrid_topology_panels(),
    function(pair) {
      .graph_pair_short(pair[[1L]], pair[[2L]])
    },
    character(1)
  )
}

#' Raincloud of solver runtime or memory (one composition per page)
#'
#' @keywords internal
#' @noRd
.plot_hybrid_resource_page <- function(
  df,
  y_lab,
  title,
  x_lab = "CT1 / CT2 topology",
  facet = TRUE
) {
  dodge <- ggplot2::position_dodge(width = 0.78)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["topology"]],
      y = .data[["value"]],
      fill = .data[["algorithm"]],
      colour = .data[["algorithm"]],
      group = .data[["algorithm"]]
    )
  ) +
    ggdist::stat_halfeye(
      orientation = "vertical",
      .width = c(0.5, 0.95),
      justification = -0.06,
      point_interval = ggdist::median_qi,
      normalize = "groups",
      scale = 0.7,
      interval_size = 2.2,
      point_size = 1.3,
      slab_linewidth = 0.45,
      slab_alpha = 0.55,
      position = dodge
    ) +
    ggplot2::geom_rug(
      sides = "l",
      alpha = 0.12,
      linewidth = 0.25,
      length = ggplot2::unit(0.03, "npc"),
      inherit.aes = TRUE
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      fill = "Solver",
      colour = "Solver",
      title = title,
      caption = paste(
        "Half-eye: median and central 50% / 95% of Monte Carlo",
        "replicates. Rug: raw draws. y-axis is log10."
      )
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 9),
      legend.position = "bottom"
    )
  if (isTRUE(facet) && "overlap_label" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~overlap_label, nrow = 1L)
  }
  p
}

#' Shared worker for runtime / memory books
#'
#' @keywords internal
#' @noRd
.save_hybrid_resource_book <- function(
  artefacts,
  file,
  column,
  y_lab,
  stem,
  scale = 1,
  data_rds = NULL
) {
  .check_plot_dependencies(need_ggdist = TRUE)
  opt <- artefacts$optimisation
  if (is.null(opt) || !column %in% names(opt)) {
    stop(
      "`artefacts$optimisation` must contain `",
      column,
      "`.",
      call. = FALSE
    )
  }
  opt <- .relevel_scenario_table(opt)
  opt$topology <- .graph_pair_short(opt$graph_ct1, opt$graph_ct2)
  opt$topology <- factor(opt$topology, levels = .hybrid_topology_levels())
  opt$value <- as.numeric(opt[[column]]) * scale
  opt$value[!is.finite(opt$value) | opt$value <= 0] <- NA_real_
  opt <- opt[is.finite(opt$value), , drop = FALSE]
  meta <- dplyr::distinct(opt, .data[["proportions"]])
  meta <- dplyr::arrange(meta, .data[["proportions"]])
  collected <- list()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 16, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    prop <- as.character(meta$proportions[[i]])
    page <- opt[as.character(opt$proportions) == prop, , drop = FALSE]
    if (nrow(page) == 0L) {
      next
    }
    page$page <- prop
    collected[[length(collected) + 1L]] <- page
    print(
      .plot_hybrid_resource_page(page, y_lab, prop),
      newpage = i > 1L
    )
  }
  if (length(collected) > 0L) {
    keep <- c(
      "ID",
      "sample_id",
      "algorithm",
      "topology",
      "overlap_label",
      "proportions",
      "value",
      "page"
    )
    keep <- intersect(keep, names(dplyr::bind_rows(collected)))
    .write_ggplot_rds(
      dplyr::bind_rows(collected)[, keep, drop = FALSE],
      data_rds,
      stem
    )
  }
  invisible(file)
}

#' Solver wall-clock time (one page per composition)
#'
#' Facets MixSim overlap; topology is on the x-axis; colour / fill /
#' grouping are the five fig03 solvers. y is elapsed seconds on a log10
#' scale, with a left-side rug of the raw Monte Carlo draws.
#'
#' @inheritParams save_hybrid_latent_projection_book
#' @return `file`, invisibly.
#' @export
#' @seealso [save_hybrid_memory_book()], [plot_mc_raincloud()]
#' @examples
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   opt <- tibble::tibble(
#'     ID = rep("V1", 20),
#'     sample_id = paste0("s", 1:20),
#'     algorithm = rep(c("lsei", "LBFGS"), each = 10),
#'     elapsed_sec = runif(20, 0.01, 0.2),
#'     graph_ct1 = "scale_free",
#'     graph_ct2 = "scale_free",
#'     overlap_label = "low",
#'     proportions = "balanced"
#'   )
#'   artefacts <- list(optimisation = opt)
#'   tf <- withr::local_tempfile(fileext = ".pdf")
#'   save_hybrid_runtime_book(artefacts, tf)
#' }
save_hybrid_runtime_book <- function(artefacts, file, data_rds = NULL) {
  .save_hybrid_resource_book(
    artefacts,
    file,
    column = "elapsed_sec",
    y_lab = "Elapsed time (s)",
    stem = "runtime",
    scale = 1,
    data_rds = data_rds
  )
}

#' Solver peak memory (one page per composition)
#'
#' Same layout as [save_hybrid_runtime_book()]. Memory is plotted in
#' mebibytes (`memory_bytes / 2^20`) on a log10 y-axis.
#'
#' @inheritParams save_hybrid_latent_projection_book
#' @return `file`, invisibly.
#' @export
#' @seealso [save_hybrid_runtime_book()]
#' @examples
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   opt <- tibble::tibble(
#'     ID = rep("V1", 20),
#'     sample_id = paste0("s", 1:20),
#'     algorithm = rep(c("lsei", "LBFGS"), each = 10),
#'     memory_bytes = runif(20, 3.8e8, 4.2e8),
#'     graph_ct1 = "scale_free",
#'     graph_ct2 = "scale_free",
#'     overlap_label = "low",
#'     proportions = "balanced"
#'   )
#'   artefacts <- list(optimisation = opt)
#'   tf <- withr::local_tempfile(fileext = ".pdf")
#'   save_hybrid_memory_book(artefacts, tf)
#' }
save_hybrid_memory_book <- function(artefacts, file, data_rds = NULL) {
  .save_hybrid_resource_book(
    artefacts,
    file,
    column = "memory_bytes",
    y_lab = "Peak memory (MiB)",
    stem = "memory",
    scale = 1 / (1024^2),
    data_rds = data_rds
  )
}

#' Fig02 runtime / memory worker (one page per CLD / variance / composition)
#'
#' @keywords internal
#' @noRd
.save_bivariate_resource_book <- function(
  artefacts,
  file,
  column,
  y_lab,
  stem,
  scale = 1,
  data_rds = NULL
) {
  .check_plot_dependencies(need_ggdist = TRUE)
  cfg <- artefacts$config
  opt <- artefacts$optimisation
  if (is.null(opt) || !column %in% names(opt)) {
    stop(
      "`artefacts$optimisation` must contain `",
      column,
      "`.",
      call. = FALSE
    )
  }
  cfg <- .relevel_scenario_table(cfg)
  opt <- .relevel_scenario_table(opt)
  meta <- .bivariate_page_meta(cfg)
  collected <- list()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = 16, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(meta))) {
    row <- meta[i, , drop = FALSE]
    ids <- .corner_ids(cfg, row)
    parts <- lapply(seq_along(ids), function(k) {
      id <- ids[[k]]
      if (is.na(id) || !nzchar(id)) {
        return(NULL)
      }
      piece <- opt[as.character(opt$ID) == id, , drop = FALSE]
      if (nrow(piece) == 0L) {
        return(NULL)
      }
      piece$topology <- names(ids)[[k]]
      piece
    })
    page <- dplyr::bind_rows(parts)
    if (nrow(page) == 0L) {
      next
    }
    page$value <- as.numeric(page[[column]]) * scale
    page$value[!is.finite(page$value) | page$value <= 0] <- NA_real_
    page <- page[is.finite(page$value), , drop = FALSE]
    if (nrow(page) == 0L) {
      next
    }
    page$topology <- factor(page$topology, levels = names(ids))
    page$page <- .page_title(row)
    collected[[length(collected) + 1L]] <- page
    print(
      .plot_hybrid_resource_page(
        page,
        y_lab,
        page$page[[1L]],
        x_lab = "Correlation corner",
        facet = FALSE
      ),
      newpage = i > 1L
    )
  }
  if (length(collected) > 0L) {
    keep <- intersect(
      c(
        "ID",
        "sample_id",
        "algorithm",
        "topology",
        "proportions",
        "value",
        "page"
      ),
      names(dplyr::bind_rows(collected))
    )
    .write_ggplot_rds(
      dplyr::bind_rows(collected)[, keep, drop = FALSE],
      data_rds,
      stem
    )
  }
  invisible(file)
}

#' Solver wall-clock time for the bivariate toy (12 pages)
#'
#' One page per CLD / variance / composition. The x-axis is the four
#' correlation corners; colour, fill, and grouping are solvers. y is
#' elapsed seconds on a log10 scale, with a left-side rug.
#'
#' @inheritParams save_hybrid_runtime_book
#' @return `file`, invisibly.
#' @export
#' @seealso [save_hybrid_runtime_book()], [save_bivariate_memory_book()]
#' @examples
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   opt <- tibble::tibble(
#'     ID = rep(c("A", "B", "C", "D"), each = 8),
#'     sample_id = paste0("s", 1:32),
#'     algorithm = rep(c("lsei", "LBFGS"), 16),
#'     elapsed_sec = runif(32, 0.01, 0.2)
#'   )
#'   cfg <- tibble::tibble(
#'     ID = c("A", "B", "C", "D"),
#'     centroids = "small_CLD",
#'     variance = "homoscedastic",
#'     proportions = "balanced",
#'     correlation_celltype1 = c(0, -0.8, 0.8, -0.8),
#'     correlation_celltype2 = c(0, -0.8, 0.8, 0.8)
#'   )
#'   artefacts <- list(config = cfg, optimisation = opt)
#'   tf <- withr::local_tempfile(fileext = ".pdf")
#'   save_bivariate_runtime_book(artefacts, tf)
#' }
save_bivariate_runtime_book <- function(artefacts, file, data_rds = NULL) {
  .save_bivariate_resource_book(
    artefacts,
    file,
    column = "elapsed_sec",
    y_lab = "Elapsed time (s)",
    stem = "runtime",
    scale = 1,
    data_rds = data_rds
  )
}

#' Solver peak memory for the bivariate toy (12 pages)
#'
#' Same layout as [save_bivariate_runtime_book()]. Memory is plotted in
#' mebibytes (`memory_bytes / 2^20`) on a log10 y-axis.
#'
#' @inheritParams save_hybrid_runtime_book
#' @return `file`, invisibly.
#' @export
#' @seealso [save_bivariate_runtime_book()]
#' @examples
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   opt <- tibble::tibble(
#'     ID = rep(c("A", "B", "C", "D"), each = 8),
#'     sample_id = paste0("s", 1:32),
#'     algorithm = rep(c("lsei", "LBFGS"), 16),
#'     memory_bytes = runif(32, 3.8e8, 4.2e8)
#'   )
#'   cfg <- tibble::tibble(
#'     ID = c("A", "B", "C", "D"),
#'     centroids = "small_CLD",
#'     variance = "homoscedastic",
#'     proportions = "balanced",
#'     correlation_celltype1 = c(0, -0.8, 0.8, -0.8),
#'     correlation_celltype2 = c(0, -0.8, 0.8, 0.8)
#'   )
#'   artefacts <- list(config = cfg, optimisation = opt)
#'   tf <- withr::local_tempfile(fileext = ".pdf")
#'   save_bivariate_memory_book(artefacts, tf)
#' }
save_bivariate_memory_book <- function(artefacts, file, data_rds = NULL) {
  .save_bivariate_resource_book(
    artefacts,
    file,
    column = "memory_bytes",
    y_lab = "Peak memory (MiB)",
    stem = "memory",
    scale = 1 / (1024^2),
    data_rds = data_rds
  )
}
