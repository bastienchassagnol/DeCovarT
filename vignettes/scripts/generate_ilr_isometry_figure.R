#!/usr/bin/env Rscript
# Isometric log-ratio (ILR) chart of the unit simplex, and the Gram-Schmidt
# construction of the Helmert basis used by DeCovarT.
#
#   A  J = 2: the simplex is a segment; the ILR chart is a line.
#   B  J = 3: an ILR grid and two segments of equal Aitchison length, drawn
#      on the ternary triangle.
#   C  The same objects in the ILR plane, where they are straight and of
#      equal length.
#   D  Gram-Schmidt on the successive differences e_k - e_{k+1}, inside the
#      plane 1^perp, returns the Helmert columns v_1, v_2.
#
# Usage (from the repository root):
#   Rscript vignettes/scripts/generate_ilr_isometry_figure.R
#
# Output: vignettes/figures/ilr_isometry_gram_schmidt.png

library(ggplot2)

output_path <- file.path(
  "vignettes",
  "figures",
  "ilr_isometry_gram_schmidt.png"
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

colour_p <- c(p1 = "#1d4e89", p2 = "#c45c26", p3 = "#7b3294")
colour_ink <- "#1c2b24"
colour_grid <- "#b8b2a7"

theme_slide <- function() {
  ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e1d8"),
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = "#0e4d3a",
        size = 15
      ),
      plot.subtitle = ggplot2::element_text(colour = "#3d4a42", size = 11),
      axis.title = ggplot2::element_text(colour = "#1c2b24", size = 12),
      axis.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(face = "bold", colour = "#0e4d3a"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Same construction as DeCovarT::helmert_basis(): column k contrasts the
# first k parts with part k + 1.
helmert_basis <- function(n_parts) {
  v <- matrix(0, nrow = n_parts, ncol = n_parts - 1L)
  for (k in seq_len(n_parts - 1L)) {
    v[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
    v[k + 1L, k] <- -k / sqrt(k * (k + 1))
  }
  v
}

ilr_inverse <- function(rho, V) {
  eta <- drop(V %*% rho)
  w <- exp(eta - max(eta))
  w / sum(w)
}

ilr <- function(p, V) drop(crossprod(V, log(p)))

# Ternary embedding: p1 at (0, 0), p2 at (1, 0), p3 at the apex.
ternary <- function(p) {
  data.frame(x = p[[2L]] + p[[3L]] / 2, y = p[[3L]] * sqrt(3) / 2)
}

# --- Panel A: J = 2 ----------------------------------------------------------

V2 <- helmert_basis(2L)
rho_line <- seq(-3, 3, by = 0.5)
p1_line <- vapply(rho_line, function(r) ilr_inverse(r, V2)[[1L]], numeric(1))

track_chart <- 0
track_simplex <- 1
panel_a_points <- rbind(
  data.frame(
    x = (rho_line + 3) / 6,
    y = track_chart,
    rho = rho_line
  ),
  data.frame(x = p1_line, y = track_simplex, rho = rho_line)
)
panel_a_links <- data.frame(
  x = (rho_line + 3) / 6,
  xend = p1_line,
  y = track_chart,
  yend = track_simplex
)
rho_ticks <- data.frame(x = (c(-3, 0, 3) + 3) / 6, label = c("-3", "0", "3"))

panel_a <- ggplot() +
  annotate(
    "segment",
    x = 0,
    xend = 1,
    y = c(track_chart, track_simplex),
    yend = c(track_chart, track_simplex),
    colour = colour_ink,
    linewidth = 0.6
  ) +
  geom_segment(
    data = panel_a_links,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = colour_grid,
    linewidth = 0.35
  ) +
  geom_point(
    data = panel_a_points,
    aes(x = x, y = y, colour = rho),
    size = 2.4
  ) +
  annotate(
    "text",
    x = c(0, 1),
    y = track_simplex + 0.16,
    label = c("p[2] == 1", "p[1] == 1"),
    parse = TRUE,
    colour = colour_p[c("p2", "p1")],
    size = 3.6
  ) +
  geom_text(
    data = rho_ticks,
    aes(x = x, y = track_chart - 0.16, label = label),
    colour = colour_ink,
    size = 3.4
  ) +
  scale_colour_gradient2(
    low = colour_p[["p2"]],
    mid = "#d9d3c7",
    high = colour_p[["p1"]],
    midpoint = 0,
    name = expression(rho)
  ) +
  scale_y_continuous(
    breaks = c(track_chart, track_simplex),
    labels = c(
      expression("ILR chart " * rho %in% R),
      expression("simplex " * Delta^1)
    ),
    limits = c(-0.35, 1.3)
  ) +
  scale_x_continuous(breaks = NULL) +
  labs(
    title = "A. J = 2: segment onto a line",
    subtitle = "Equal steps in rho crowd towards the vertices",
    x = NULL,
    y = NULL
  ) +
  theme_slide() +
  theme(
    panel.grid.major = element_blank(),
    legend.key.width = grid::unit(1.2, "cm"),
    legend.key.height = grid::unit(0.3, "cm")
  )

# --- Panels B and C: J = 3 ---------------------------------------------------

V3 <- helmert_basis(3L)
grid_levels <- seq(-3, 3, by = 1)
grid_path <- seq(-4.5, 4.5, length.out = 200)

grid_lines <- do.call(
  rbind,
  c(
    lapply(grid_levels, function(level) {
      data.frame(
        rho1 = level,
        rho2 = grid_path,
        family = "rho1",
        line = paste0("rho1_", level)
      )
    }),
    lapply(grid_levels, function(level) {
      data.frame(
        rho1 = grid_path,
        rho2 = level,
        family = "rho2",
        line = paste0("rho2_", level)
      )
    })
  )
)
grid_ternary <- do.call(
  rbind,
  lapply(seq_len(nrow(grid_lines)), function(i) {
    ternary(ilr_inverse(c(grid_lines$rho1[[i]], grid_lines$rho2[[i]]), V3))
  })
)
grid_lines <- cbind(grid_lines, grid_ternary)

# Two segments of Aitchison length 1: one at the barycentre, one near the
# p3 vertex. Both are the same translation in rho.
segment_step <- c(1, 0)
segment_starts <- list(centre = c(-0.5, 0), vertex = c(-0.5, -2.4))
segments_rho <- do.call(
  rbind,
  lapply(names(segment_starts), function(name) {
    start <- segment_starts[[name]]
    t <- seq(0, 1, length.out = 60)
    data.frame(
      segment = name,
      rho1 = start[[1L]] + t * segment_step[[1L]],
      rho2 = start[[2L]] + t * segment_step[[2L]]
    )
  })
)
segments_ternary <- cbind(
  segments_rho,
  do.call(
    rbind,
    lapply(seq_len(nrow(segments_rho)), function(i) {
      ternary(ilr_inverse(
        c(segments_rho$rho1[[i]], segments_rho$rho2[[i]]),
        V3
      ))
    })
  )
)

p_example <- c(0.5, 0.3, 0.2)
rho_example <- ilr(p_example, V3)
example_ternary <- ternary(p_example)

triangle <- data.frame(
  x = c(0, 1, 0.5, 0),
  y = c(0, 0, sqrt(3) / 2, 0)
)
vertex_labels <- data.frame(
  x = c(-0.05, 1.05, 0.5),
  y = c(-0.05, -0.05, sqrt(3) / 2 + 0.06),
  label = c("p[1] == 1", "p[2] == 1", "p[3] == 1"),
  colour = unname(colour_p)
)
segment_colours <- c(centre = "#0e4d3a", vertex = "#b8860b")
segment_labels <- c(
  centre = "near the barycentre",
  vertex = "near the vertex p3"
)

panel_b <- ggplot() +
  geom_path(
    data = grid_lines,
    aes(x = x, y = y, group = line, linetype = family),
    colour = colour_grid,
    linewidth = 0.35
  ) +
  geom_path(data = triangle, aes(x = x, y = y), colour = colour_ink) +
  geom_path(
    data = segments_ternary,
    aes(x = x, y = y, colour = segment, group = segment),
    linewidth = 1.6,
    lineend = "round"
  ) +
  geom_point(
    data = example_ternary,
    aes(x = x, y = y),
    shape = 21,
    fill = "white",
    colour = colour_ink,
    size = 3,
    stroke = 1
  ) +
  annotate(
    "text",
    x = example_ternary$x + 0.04,
    y = example_ternary$y + 0.03,
    label = "(0.5, 0.3, 0.2)",
    hjust = 0,
    size = 3.4,
    colour = colour_ink
  ) +
  geom_text(
    data = vertex_labels,
    aes(x = x, y = y, label = label),
    colour = vertex_labels$colour,
    parse = TRUE,
    size = 3.8
  ) +
  scale_colour_manual(
    values = segment_colours,
    labels = segment_labels,
    name = "Aitchison length 1:"
  ) +
  scale_linetype_manual(
    values = c(rho1 = "solid", rho2 = "22"),
    labels = c(
      rho1 = expression("constant " * rho[1]),
      rho2 = expression("constant " * rho[2])
    ),
    name = NULL
  ) +
  coord_equal(xlim = c(-0.1, 1.1), ylim = c(-0.1, 0.95), expand = FALSE) +
  labs(
    title = "B. J = 3: the simplex",
    subtitle = "ILR grid lines are curves; equal lengths look unequal",
    x = NULL,
    y = NULL
  ) +
  theme_slide() +
  theme(
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0, "pt")
  )

panel_c <- ggplot() +
  geom_path(
    data = grid_lines,
    aes(x = rho1, y = rho2, group = line, linetype = family),
    colour = colour_grid,
    linewidth = 0.35
  ) +
  geom_path(
    data = segments_ternary,
    aes(x = rho1, y = rho2, colour = segment, group = segment),
    linewidth = 1.6,
    lineend = "round"
  ) +
  geom_point(
    aes(x = rho_example[[1L]], y = rho_example[[2L]]),
    shape = 21,
    fill = "white",
    colour = colour_ink,
    size = 3,
    stroke = 1
  ) +
  annotate(
    "text",
    x = rho_example[[1L]] + 0.2,
    y = rho_example[[2L]] + 0.25,
    label = sprintf("(%.3f, %.3f)", rho_example[[1L]], rho_example[[2L]]),
    hjust = 0,
    size = 3.4,
    colour = colour_ink
  ) +
  scale_colour_manual(
    values = segment_colours,
    labels = segment_labels,
    name = "Aitchison length 1:"
  ) +
  scale_linetype_manual(
    values = c(rho1 = "solid", rho2 = "22"),
    labels = c(
      rho1 = expression("constant " * rho[1]),
      rho2 = expression("constant " * rho[2])
    ),
    name = NULL
  ) +
  coord_equal(xlim = c(-3.4, 3.4), ylim = c(-3.4, 3.4), expand = FALSE) +
  labs(
    title = "C. J = 3: the ILR chart",
    subtitle = "Orthonormal grid; Euclidean length = Aitchison length",
    x = expression(rho[1] == sqrt(1 / 2) * log(p[1] / p[2])),
    y = expression(rho[2] == sqrt(2 / 3) * log(sqrt(p[1] * p[2]) / p[3]))
  ) +
  theme_slide() +
  theme(
    panel.grid.major = element_blank(),
    legend.box = "vertical",
    legend.spacing.y = grid::unit(0, "pt")
  )

# --- Panel D: Gram-Schmidt in the plane 1^perp ------------------------------

# Vectors of R^3 orthogonal to 1 are drawn in the plane they span, with
# v1 as horizontal and v2 as vertical axis (an isometric drawing chart).
plane_coords <- function(w) drop(crossprod(V3, w))

a1 <- c(1, -1, 0)
a2 <- c(0, 1, -1)
v1 <- a1 / sqrt(sum(a1^2))
projection <- sum(a2 * v1) * v1
residual <- a2 - projection
v2 <- residual / sqrt(sum(residual^2))
stopifnot(max(abs(cbind(v1, v2) - V3)) < 1e-12)

centred_axes <- t(vapply(
  1:3,
  function(j) plane_coords(diag(3)[, j] - 1 / 3),
  numeric(2)
))

arrow_df <- function(w, role) {
  xy <- plane_coords(w)
  data.frame(x = 0, y = 0, xend = xy[[1L]], yend = xy[[2L]], role = role)
}
gs_arrows <- rbind(
  arrow_df(a1, "a1"),
  arrow_df(a2, "a2"),
  arrow_df(v1, "v1"),
  arrow_df(v2, "v2")
)
proj_xy <- plane_coords(projection)
a2_xy <- plane_coords(a2)

gs_colours <- c(a1 = "#8a8f98", a2 = "#8a8f98", v1 = "#0e4d3a", v2 = "#b8860b")
gs_labels <- data.frame(
  x = c(gs_arrows$xend * 1.12),
  y = c(gs_arrows$yend * 1.12) + c(-0.12, 0.08, 0.12, 0),
  label = c(
    "a[1] == e[1] - e[2]",
    "a[2] == e[2] - e[3]",
    "v[1] == a[1] / '||' * a[1] * '||'",
    "v[2] %prop% a[2] - (a[2] %.% v[1]) * v[1]"
  ),
  colour = gs_colours[gs_arrows$role],
  hjust = c(0.5, 0.5, 0.5, -0.05)
)
axis_labels <- data.frame(
  x = centred_axes[, 1L] * 1.18,
  y = centred_axes[, 2L] * 1.18,
  label = c(
    "clr(e[1])",
    "clr(e[2])",
    "clr(e[3])"
  ),
  colour = unname(colour_p)
)

panel_d <- ggplot() +
  geom_segment(
    data = data.frame(
      xend = centred_axes[, 1L],
      yend = centred_axes[, 2L],
      colour = unname(colour_p)
    ),
    aes(x = 0, y = 0, xend = xend, yend = yend),
    colour = unname(colour_p),
    linewidth = 0.5,
    linetype = "22"
  ) +
  geom_text(
    data = axis_labels,
    aes(x = x, y = y, label = label),
    colour = axis_labels$colour,
    parse = TRUE,
    size = 3.4
  ) +
  annotate(
    "segment",
    x = a2_xy[[1L]],
    xend = proj_xy[[1L]],
    y = a2_xy[[2L]],
    yend = proj_xy[[2L]],
    colour = colour_ink,
    linewidth = 0.4,
    linetype = "dotted"
  ) +
  annotate(
    "segment",
    x = 0,
    xend = proj_xy[[1L]],
    y = 0,
    yend = proj_xy[[2L]],
    colour = "#c45c26",
    linewidth = 1.4
  ) +
  annotate(
    "text",
    x = proj_xy[[1L]],
    y = -0.16,
    label = "(a[2] %.% v[1]) * v[1]",
    parse = TRUE,
    colour = "#c45c26",
    size = 3.2
  ) +
  geom_segment(
    data = gs_arrows,
    aes(x = x, y = y, xend = xend, yend = yend, colour = role),
    linewidth = 1,
    arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed")
  ) +
  geom_text(
    data = gs_labels,
    aes(x = x, y = y, label = label),
    colour = gs_labels$colour,
    hjust = gs_labels$hjust,
    parse = TRUE,
    size = 3.2
  ) +
  scale_colour_manual(values = gs_colours, guide = "none") +
  coord_equal(xlim = c(-1.6, 2.1), ylim = c(-1.25, 1.55), expand = FALSE) +
  labs(
    title = "D. Gram-Schmidt builds the Helmert basis",
    subtitle = expression(
      "In the plane orthogonal to " * bold(1) * ", orthogonalise " * a[k] ==
        e[k] - e[k + 1]
    ),
    x = expression(v[1] * " axis"),
    y = expression(v[2] * " axis")
  ) +
  theme_slide() +
  theme(axis.text = element_blank(), panel.grid.major = element_blank())

# --- Assemble ---------------------------------------------------------------

ragg::agg_png(
  output_path,
  width = 13,
  height = 12,
  units = "in",
  res = 200,
  background = "white"
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2L, 2L)))
panels <- list(panel_a, panel_d, panel_b, panel_c)
positions <- list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
for (i in seq_along(panels)) {
  print(
    panels[[i]],
    vp = grid::viewport(
      layout.pos.row = positions[[i]][[1L]],
      layout.pos.col = positions[[i]][[2L]]
    )
  )
}
invisible(grDevices::dev.off())
message("Wrote ", output_path)
