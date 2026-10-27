# Shared pieces of the figure books: palettes, factor orders and labels,
# the page layout of each figure, the balance basis of the log-likelihood
# pictures, shared legends and multi-page saving.

# Palettes ---------------------------------------------------------------------

.method_levels <- c(
  "lsei",
  "NNLS",
  "RLM",
  "CIBERSORT",
  "Newton-Raphson",
  "L-BFGS-B",
  "Firth",
  "MEAD",
  "RNA-Sieve",
  "ReDeconv"
)

.method_palette <- c(
  "lsei" = "#999999",
  "NNLS" = "#5F5F5F",
  "RLM" = "#B8A07E",
  "CIBERSORT" = "#E69F00",
  "Newton-Raphson" = "#1D4E89",
  "L-BFGS-B" = "#56B4E9",
  "Firth" = "#009E73",
  "MEAD" = "#CC79A7",
  "RNA-Sieve" = "#D55E00",
  "ReDeconv" = "#8C6D31"
)

.ci_levels <- c(
  "wald",
  "ilr",
  "ilr_reml",
  "profile",
  "profile_reml",
  "bootstrap",
  "oracle",
  "native",
  "native (transcript fractions)",
  "fisher",
  "godambe",
  "native (released code)"
)

.ci_labels <- c(
  "wald" = "Wald (delta)",
  "ilr" = "ILR, mapped back",
  "ilr_reml" = "ILR, REML",
  "profile" = "Profile likelihood",
  "profile_reml" = "Profile, REML",
  "bootstrap" = "Parametric bootstrap",
  "oracle" = "Oracle (Fisher at p*)",
  "native" = "Method's own interval",
  "native (transcript fractions)" = "Own interval, transcript fractions",
  "fisher" = "Fisher information (Wald)",
  "godambe" = "Godambe sandwich",
  "native (released code)" = "Released code, as published"
)

.ci_palette <- c(
  "Wald (delta)" = "#1D4E89",
  "ILR, mapped back" = "#56B4E9",
  "ILR, REML" = "#7FB8DA",
  "Profile likelihood" = "#009E73",
  "Profile, REML" = "#66C2A5",
  "Parametric bootstrap" = "#E69F00",
  "Oracle (Fisher at p*)" = "#6B6B6B",
  "Method's own interval" = "#CC79A7",
  "Own interval, transcript fractions" = "#CC79A7",
  "Fisher information (Wald)" = "#B2182B",
  "Godambe sandwich" = "#D55E00",
  "Released code, as published" = "#BDBDBD"
)

.method_factor <- function(x) {
  present <- unique(as.character(x))
  factor(
    x,
    levels = c(
      intersect(.method_levels, present),
      setdiff(present[!is.na(present)], .method_levels)
    )
  )
}

.likelihood_methods <- c("Newton-Raphson", "L-BFGS-B", "Firth")

# Theme of every page: the package theme on a white background, with the
# type sizes of a full-page figure.
.page_theme <- function(base_size = 11) {
  theme_decovart_facets(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.subtitle = ggplot2::element_text(size = 11, colour = "grey25"),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 9,
        colour = "grey25",
        lineheight = 1.1
      ),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "bottom"
    )
}

# Pages and panels -------------------------------------------------------------

# The four correlation pairs of the fig02 density, forest and diagnostic
# pages: no correlation, both negative, both positive, and opposite signs.
.corner_rho <- c("0_0", "-0.8_-0.8", "0.8_0.8", "-0.8_0.8")

# Scenario table with the page and panel of every scenario. A page groups
# the scenarios that share every factor except the two shown inside it:
# fig02 the two within-type correlations, fig03 the graph families of
# types 1 and 2, fig05 the cell sizes and the number of cells per bulk.
.scenario_layout <- function(config, figure) {
  config <- tibble::as_tibble(config)
  if (identical(figure, "fig02")) {
    config$proportions <- factor(
      config$proportions,
      c("balanced", "moderately unbalanced", "highly unbalanced")
    )
    config$centroids <- factor(config$centroids, c("small_CLD", "large_CLD"))
    config$variance <- factor(
      config$variance,
      c("homoscedastic", "heteroscedastic")
    )
    config <- dplyr::arrange(
      config,
      .data$centroids,
      .data$variance,
      .data$proportions
    )
    config$page <- paste(
      sub("_CLD", " centroid gap", config$centroids),
      config$variance,
      config$proportions,
      sep = " / "
    )
    config$panel_x <- factor(config$correlation_celltype1)
    config$panel_y <- factor(config$correlation_celltype2)
    config$panel_x_title <- "Correlation, cell type 1"
    config$panel_y_title <- "Correlation, cell type 2"
    config$panel <- sprintf(
      "rho = (%+.1f, %+.1f)",
      config$correlation_celltype1,
      config$correlation_celltype2
    )
    config$corner <- paste(
      config$correlation_celltype1,
      config$correlation_celltype2,
      sep = "_"
    ) %in%
      .corner_rho
  } else if (identical(figure, "fig03")) {
    graph <- c(scale_free = "Scale-free", stochastic_block_model = "SBM")
    config$proportions <- factor(
      config$proportions,
      c("balanced", "moderately unbalanced", "highly unbalanced")
    )
    config$overlap_label <- factor(
      config$overlap_label,
      c("low", "moderate", "high")
    )
    config <- dplyr::arrange(config, .data$proportions, .data$overlap_label)
    config$page <- paste(
      config$proportions,
      paste(config$overlap_label, "overlap"),
      sep = " / "
    )
    config$panel_x <- factor(unname(graph[config$graph_ct1]), graph)
    config$panel_y <- factor(unname(graph[config$graph_ct2]), graph)
    config$panel_x_title <- "Graph of cell type 1"
    config$panel_y_title <- "Graph of cell type 2"
    config$panel <- paste0(
      "CT1 ",
      config$panel_x,
      " / CT2 ",
      config$panel_y
    )
    config$corner <- TRUE
  } else {
    # Pages: genes, then cells per bulk, composition and sibling overlap;
    # inside a page, platform scale a in rows and cell sizes in columns.
    config$proportions <- forcats::fct_relevel(
      factor(config$proportions),
      "balanced",
      "unbalanced"
    )
    config$overlap_label <- forcats::fct_relevel(
      factor(config$overlap_label),
      "low",
      "high"
    )
    config <- dplyr::arrange(
      config,
      .data$n_genes,
      .data$n_cells,
      .data$proportions,
      .data$overlap_label,
      .data$platform_scale,
      .data$cell_size_label
    )
    config$page <- sprintf(
      "G = %d / %s cells per bulk / %s / %s sibling overlap",
      config$n_genes,
      formatC(config$n_cells, format = "g"),
      config$proportions,
      config$overlap_label
    )
    config$panel_x <- forcats::fct_relevel(
      factor(ifelse(
        config$cell_size_label == "equal",
        "equal sizes",
        "sizes 1:1:5:5"
      )),
      "equal sizes",
      "sizes 1:1:5:5"
    )
    config$panel_y <- factor(
      paste("a =", config$platform_scale),
      paste("a =", sort(unique(config$platform_scale)))
    )
    config$panel_x_title <- "Cell sizes"
    config$panel_y_title <- "Platform scale"
    config$panel <- paste(config$panel_y, config$panel_x, sep = " / ")
    config$corner <- TRUE
  }
  config$grid_panels <- identical(figure, "fig05")
  config$page <- factor(config$page, unique(config$page))
  config$panel <- factor(
    config$panel,
    unique(config$panel[order(config$panel_y, config$panel_x)])
  )
  config
}

# Panels of a page: a 2 x 2 grid (rows panel_y, columns panel_x, as two
# strip lines) when the figure has one, every cell shown even when it was
# not simulated; otherwise one panel per scenario. facet_wrap() keeps free
# scales independent per panel, which the zoomed likelihood pages need.
.panel_facets <- function(df, scales = "fixed") {
  if (isTRUE(df$grid_panels[[1L]])) {
    ggplot2::facet_wrap(
      ggplot2::vars(.data$panel_y, .data$panel_x),
      ncol = nlevels(df$panel_x),
      drop = FALSE,
      scales = scales
    )
  } else {
    ggplot2::facet_wrap(~panel, ncol = 2L, scales = scales)
  }
}

# Text in the grid cells of a page that hold no simulated scenario, in the
# top-left corner; a panel with free scales has no range of its own, so the
# zoomed pages place it at (x, y) = (0, 0) instead.
.not_simulated <- function(df, at_origin = FALSE) {
  if (!isTRUE(df$grid_panels[[1L]])) {
    return(NULL)
  }
  cells <- expand.grid(
    panel_y = levels(df$panel_y),
    panel_x = levels(df$panel_x)
  )
  present <- paste(df$panel_y, df$panel_x)
  cells <- cells[!paste(cells$panel_y, cells$panel_x) %in% present, ]
  if (nrow(cells) == 0L) {
    return(NULL)
  }
  cells$panel_y <- factor(cells$panel_y, levels(df$panel_y))
  cells$panel_x <- factor(cells$panel_x, levels(df$panel_x))
  ggplot2::geom_text(
    data = cells,
    ggplot2::aes(
      x = if (at_origin) 0 else -Inf,
      y = if (at_origin) 0 else Inf
    ),
    label = "not simulated",
    hjust = if (at_origin) 0.5 else -0.15,
    vjust = if (at_origin) 0.5 else 1.8,
    colour = "grey45",
    size = 3.4,
    inherit.aes = FALSE
  )
}

# Attach the page and panel columns of `layout` to a result table.
.with_layout <- function(tbl, layout) {
  keep <- c(
    "ID",
    "page",
    "panel",
    "panel_x",
    "panel_y",
    "panel_x_title",
    "panel_y_title",
    "corner",
    "grid_panels"
  )
  tbl <- tbl[setdiff(names(tbl), setdiff(keep, "ID"))]
  dplyr::inner_join(tbl, layout[keep], by = "ID")
}

.celltype_names <- function(benchmark) {
  colnames(benchmark$theta_true[[1L]]$mu)
}

# Long table of the Monte Carlo estimates: one row per sample, method and
# cell type, with the target of that sample.
.estimates_long <- function(benchmark, layout) {
  cts <- .celltype_names(benchmark)
  opt <- benchmark$optimisation
  opt$sample_index <- as.integer(sub("sample_", "", opt$sample_id))
  long <- tidyr::pivot_longer(
    opt[c("ID", "algorithm", "sample_index", "numerical_converged", cts)],
    tidyr::all_of(cts),
    names_to = "cell_type",
    values_to = "estimate"
  )
  long$truth <- mapply(
    function(id, s, ct) benchmark$data[[id]]$truth[ct, s],
    long$ID,
    long$sample_index,
    long$cell_type
  )
  long$algorithm <- .method_factor(long$algorithm)
  .with_layout(long, layout)
}

# Balance basis ----------------------------------------------------------------

# ILR basis of a sequential binary partition: each column contrasts the
# parts in `plus` with those in `minus`, with coefficients
# sqrt(s / (r (r + s))) and -sqrt(r / (s (r + s))) (r = |plus|, s = |minus|).
.balance_basis <- function(n_parts, partitions) {
  V <- vapply(
    partitions,
    function(partition) {
      r <- length(partition$plus)
      s <- length(partition$minus)
      v <- numeric(n_parts)
      v[partition$plus] <- sqrt(s / (r * (r + s)))
      v[partition$minus] <- -sqrt(r / (s * (r + s)))
      v
    },
    numeric(n_parts)
  )
  V <- matrix(V, nrow = n_parts)
  stopifnot(
    max(abs(crossprod(V) - diag(ncol(V)))) < 1e-12,
    max(abs(colSums(V))) < 1e-12
  )
  V
}

# The balances follow the lineages: fig03 contrasts the two close types,
# then both against the third; fig05 contrasts the siblings of each lineage,
# then the two lineages (held at its true value in the pictures).
.figure_balances <- function(figure) {
  switch(
    figure,
    fig02 = list(
      basis = .balance_basis(2L, list(list(plus = 1L, minus = 2L))),
      labels = "z = balance(CT1 | CT2)"
    ),
    fig03 = list(
      basis = .balance_basis(
        3L,
        list(list(plus = 1L, minus = 2L), list(plus = 1:2, minus = 3L))
      ),
      labels = c("z1 = balance(CT1 | CT2)", "z2 = balance(CT1, CT2 | CT3)")
    ),
    fig05 = list(
      basis = .balance_basis(
        4L,
        list(
          list(plus = 1L, minus = 2L),
          list(plus = 3L, minus = 4L),
          list(plus = 1:2, minus = 3:4)
        )
      ),
      labels = c(
        "z1 = balance(CT1 | CT2)",
        "z2 = balance(CT3 | CT4)",
        "z3 = balance(CT1, CT2 | CT3, CT4)"
      )
    )
  )
}

.balance_limits <- c(-4, 4)
.loglik_floor <- -10

# Legends and saving -----------------------------------------------------------

# A page of panels sharing one legend: the legend is taken once from the
# first panel and placed under the grid.
.shared_legend_page <- function(
  panels,
  title,
  caption = NULL,
  ncol = 2L,
  legend_height = 0.09
) {
  legend <- cowplot::get_legend(
    panels[[1L]] + ggplot2::theme(legend.position = "bottom")
  )
  stripped <- lapply(panels, function(p) {
    p + ggplot2::theme(legend.position = "none")
  })
  grid <- cowplot::plot_grid(plotlist = stripped, ncol = ncol)
  page <- cowplot::plot_grid(
    grid,
    legend,
    ncol = 1L,
    rel_heights = c(1, legend_height)
  )
  page <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        title,
        fontface = "bold",
        size = 15,
        x = 0.01,
        hjust = 0
      ),
    page,
    cowplot::ggdraw() +
      cowplot::draw_label(
        if (is.null(caption)) "" else caption,
        size = 9,
        x = 0.01,
        hjust = 0,
        colour = "grey25"
      ),
    ncol = 1L,
    rel_heights = c(0.05, 1, if (is.null(caption)) 0.01 else 0.08)
  )
  page +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Write a list of pages (ggplots or cowplot drawings) to one PDF, one page
# per element; cairo keeps the Greek letters and the vector graphics.
.save_pages <- function(pages, file, width, height) {
  pages <- Filter(Negate(is.null), pages)
  if (length(pages) == 0L) {
    return(invisible(NULL))
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  book <- gridExtra::marrangeGrob(
    lapply(pages, function(page) {
      if (inherits(page, "patchwork")) {
        patchwork::patchworkGrob(page)
      } else {
        ggplot2::ggplotGrob(page)
      }
    }),
    nrow = 1L,
    ncol = 1L,
    top = NULL
  )
  ggplot2::ggsave(
    file,
    book,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    limitsize = FALSE
  )
  invisible(file)
}

.signif_label <- function(x, digits = 3L) {
  ifelse(
    is.finite(x),
    formatC(signif(x, digits), format = "g", digits = digits),
    "NA"
  )
}

# Marker legend of the likelihood pictures: white diamond for the true
# composition, red circle for the maximum (observed surface) or the Monte
# Carlo mean of the estimates (expected surface), grey dots for estimates.
.marker_legend <- function(labels) {
  shapes <- c(23, 21, 21)[seq_along(labels)]
  fills <- c("white", "#E41A1C", "grey55")[seq_along(labels)]
  sizes <- c(4, 3.6, 2)[seq_along(labels)]
  df <- data.frame(
    marker = factor(labels, labels),
    x = seq_along(labels),
    y = 1
  )
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$x, .data$y, shape = .data$marker)
  ) +
    ggplot2::geom_point(fill = "white") +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = stats::setNames(shapes, labels),
      guide = ggplot2::guide_legend(
        nrow = 1L,
        override.aes = list(fill = fills, size = sizes, colour = "black")
      )
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
  cowplot::get_legend(p)
}

# Stack a page, a legend grob and an optional caption line.
.with_legend <- function(page, legend, height = 0.06) {
  cowplot::plot_grid(page, legend, ncol = 1L, rel_heights = c(1, height)) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Segments of a dendrogram drawn above leaves placed at 1, ..., n in the
# order of `hc$order`.
.dendrogram_segments <- function(hc) {
  n <- length(hc$order)
  leaf_x <- numeric(n)
  leaf_x[hc$order] <- seq_len(n)
  node_x <- numeric(n - 1L)
  node_h <- hc$height
  child <- function(k) {
    if (k < 0) {
      list(x = leaf_x[[-k]], h = 0)
    } else {
      list(x = node_x[[k]], h = node_h[[k]])
    }
  }
  segs <- vector("list", n - 1L)
  for (i in seq_len(n - 1L)) {
    a <- child(hc$merge[i, 1L])
    b <- child(hc$merge[i, 2L])
    node_x[[i]] <- (a$x + b$x) / 2
    segs[[i]] <- data.frame(
      x = c(a$x, b$x, a$x),
      xend = c(a$x, b$x, b$x),
      y = c(a$h, b$h, node_h[[i]]),
      yend = c(node_h[[i]], node_h[[i]], node_h[[i]])
    )
  }
  do.call(rbind, segs)
}
