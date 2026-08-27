#!/usr/bin/env Rscript
# Illustration of the four strata of the 3-simplex (J = 4 cell types) and
# the tangent cone that governs the null distribution of the DeCovarT
# likelihood-ratio statistic on each of them.
#
# Usage (from the repository root):
#   Rscript scripts/auxiliary/generate_simplex_boundary_figure.R
#
# Output: vignettes/figures/simplex_boundary_tangent_cones.png

library(ggplot2)

output_path <- file.path(
  "vignettes",
  "figures",
  "simplex_boundary_tangent_cones.png"
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Vertex positions for a tetrahedron seen from outside: ct1, ct2 and ct4
# form the visible hull and ct3 is the far vertex, so the three edges
# incident to it read as receding. The barycentric embedding
# p -> sum_j p_j v_j is affine, hence maps the simplex onto that hull.
tetrahedron <- rbind(
  ct1 = c(-1.00, -0.58),
  ct2 = c(1.00, -0.58),
  ct3 = c(0.10, 0.02),
  ct4 = c(0.00, 1.12)
)

embed <- function(p) {
  coords <- drop(matrix(p, nrow = 1L) %*% tetrahedron)
  data.frame(x = coords[[1L]], y = coords[[2L]])
}

vertices <- data.frame(
  x = tetrahedron[, 1L],
  y = tetrahedron[, 2L],
  celltype = rownames(tetrahedron),
  # Nudge labels away from the hull centroid.
  label_x = tetrahedron[, 1L] * 1.13 + c(0, 0, 0.28, 0),
  label_y = tetrahedron[, 2L] * 1.13 + c(-0.12, -0.12, 0.14, 0.12),
  stringsAsFactors = FALSE
)

edge_pairs <- utils::combn(4L, 2L)
edges <- do.call(
  rbind,
  lapply(seq_len(ncol(edge_pairs)), function(k) {
    from <- vertices[edge_pairs[1L, k], c("x", "y")]
    to <- vertices[edge_pairs[2L, k], c("x", "y")]
    data.frame(
      x = from$x,
      y = from$y,
      xend = to$x,
      yend = to$y,
      # Edges touching the far vertex ct3 are drawn as hidden lines.
      hidden = 3L %in% edge_pairs[, k]
    )
  })
)

# One panel per stratum: q counts the constraints active on the boundary.
strata <- list(
  list(
    label = "Interior: every cell type present",
    p = c(0.40, 0.25, 0.20, 0.15),
    q = 0L,
    highlight = NULL
  ),
  list(
    label = "Face: one cell type absent",
    p = c(0.50, 0.30, 0.20, 0.00),
    q = 1L,
    highlight = c(1L, 2L, 3L)
  ),
  list(
    label = "Edge: two cell types absent",
    p = c(0.60, 0.40, 0.00, 0.00),
    q = 2L,
    highlight = c(1L, 2L)
  ),
  list(
    label = "Vertex: a single pure cell type",
    p = c(1.00, 0.00, 0.00, 0.00),
    q = 3L,
    highlight = 1L
  )
)

mixture_label <- function(q) {
  if (q == 0L) {
    return("chi[1]^2~~(Wilks)")
  }
  weights <- stats::dbinom(0:q, size = q, prob = 0.5)
  terms <- sprintf(
    "%s*chi[%d]^2",
    MASS::fractions(weights),
    0:q
  )
  paste(terms, collapse = " + ")
}

panel_levels <- vapply(strata, function(s) s$label, character(1))

panel_frame <- function(component) {
  do.call(
    rbind,
    lapply(strata, function(s) {
      value <- component(s)
      if (is.null(value) || nrow(value) == 0L) {
        return(NULL)
      }
      value$panel <- factor(s$label, levels = panel_levels)
      value
    })
  )
}

edges_all <- panel_frame(function(s) edges)
vertices_all <- panel_frame(function(s) vertices)

points_all <- panel_frame(function(s) {
  out <- embed(s$p)
  out$q <- s$q
  out
})

# Shade the active face (q = 1) or thicken the active edge / vertex.
faces_all <- panel_frame(function(s) {
  if (is.null(s$highlight) || length(s$highlight) != 3L) {
    return(NULL)
  }
  vertices[s$highlight, c("x", "y")]
})

active_edges_all <- panel_frame(function(s) {
  if (is.null(s$highlight) || length(s$highlight) != 2L) {
    return(NULL)
  }
  from <- vertices[s$highlight[[1L]], c("x", "y")]
  to <- vertices[s$highlight[[2L]], c("x", "y")]
  data.frame(x = from$x, y = from$y, xend = to$x, yend = to$y)
})

captions_all <- panel_frame(function(s) {
  data.frame(
    x = 0,
    y = min(vertices$y) - 0.42,
    text = sprintf(
      "q = %d active constraint%s",
      s$q,
      if (s$q == 1L) "" else "s"
    )
  )
})

mixtures_all <- panel_frame(function(s) {
  data.frame(
    x = 0,
    y = min(vertices$y) - 0.72,
    text = mixture_label(s$q)
  )
})

highlight_colour <- "#B2182B"
context_colour <- "#4D4D4D"

plot_object <- ggplot() +
  geom_polygon(
    data = faces_all,
    aes(x = .data$x, y = .data$y),
    fill = highlight_colour,
    alpha = 0.16
  ) +
  geom_segment(
    data = edges_all,
    aes(
      x = .data$x,
      y = .data$y,
      xend = .data$xend,
      yend = .data$yend,
      linetype = .data$hidden
    ),
    colour = context_colour,
    linewidth = 0.35,
    alpha = 0.6
  ) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22")) +
  geom_segment(
    data = active_edges_all,
    aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    colour = highlight_colour,
    linewidth = 1.5
  ) +
  geom_point(
    data = vertices_all,
    aes(x = .data$x, y = .data$y),
    colour = context_colour,
    size = 1.6
  ) +
  geom_text(
    data = vertices_all,
    aes(x = .data$label_x, y = .data$label_y, label = .data$celltype),
    colour = context_colour,
    size = 3.1
  ) +
  geom_point(
    data = points_all,
    aes(x = .data$x, y = .data$y),
    colour = highlight_colour,
    size = 3.4
  ) +
  geom_text(
    data = captions_all,
    aes(x = .data$x, y = .data$y, label = .data$text),
    size = 3.2,
    colour = context_colour
  ) +
  geom_text(
    data = mixtures_all,
    aes(x = .data$x, y = .data$y, label = .data$text),
    parse = TRUE,
    size = 3.4,
    colour = highlight_colour
  ) +
  facet_wrap(~panel, nrow = 1L) +
  coord_equal(clip = "off") +
  labs(
    title = paste(
      "Where the estimate sits on the simplex decides the null",
      "distribution of the likelihood-ratio statistic"
    ),
    subtitle = paste(
      "Cellular compositions of J = 4 cell types live on the tetrahedron",
      "\u0394\u00b3. Each stratum is approximated by a different tangent",
      "cone,\nso the boundary-corrected calibration is a binomial mixture",
      "of chi-square laws (Chernoff 1954; Self and Liang 1987)."
    )
  ) +
  theme_void(base_size = 11) +
  guides(linetype = "none") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0, size = 9, colour = context_colour),
    strip.text = element_text(size = 10, face = "bold", margin = margin(b = 6)),
    plot.margin = margin(12, 22, 16, 22)
  )

ggsave(
  filename = output_path,
  plot = plot_object,
  width = 12.4,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

message("Wrote ", normalizePath(output_path))
