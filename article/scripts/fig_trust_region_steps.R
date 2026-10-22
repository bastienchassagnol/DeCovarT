#!/usr/bin/env Rscript
# Trust-region step control of stats::nlminb on the single-cell profile
# log-likelihood (Remark rem:step-size of the appendix).
#
#   A  ILR plane of a J = 3 composition, contours of the profile
#      log-likelihood. Every trial point that nlminb evaluated, joined to
#      the iterate it started from; the circle around that iterate has the
#      length of the trial step, which equals the trust radius when the
#      step is on the boundary of the ball.
#   B  Zoom on the last iterations.
#   C  Length of every trial step, on a log scale.
#   D  Actual over predicted increase of every trial, with the PORT
#      thresholds 1e-4 (reject below) and 0.1 (accept, cut the radius).
#
# Nothing is re-implemented: the trials are the objective evaluations of
# nlminb, the iterates are the points where it asked for the gradient,
# and the predicted increase is the quadratic model at the iterate.
#
# Usage (from the repository root):
#   Rscript article/scripts/fig_trust_region_steps.R
#
# Output: article/figs/fig_trust_region_steps.svg (vector, the source) and
# article/figs/fig_trust_region_steps.pdf (Inkscape export of that SVG,
# included by pdflatex).

devtools::load_all(".", quiet = TRUE)
library(ggplot2)

output_svg <- file.path("article", "figs", "fig_trust_region_steps.svg")
output_pdf <- sub("\\.svg$", ".pdf", output_svg)

colour_accept <- "#1d4e89"
colour_reject <- "#c45c26"
colour_superseded <- "#7b3294"
outcome_colours <- c(
  accepted = colour_accept,
  superseded = colour_superseded,
  rejected = colour_reject
)
colour_ink <- "#1c2b24"

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
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# --- One single-cell bulk sample, J = 3 ---------------------------------

moments <- withr::with_seed(3L, {
  simulate_hierarchical_grn_moments(
    n_genes = 10L,
    n_celltypes = 3L,
    mean_scale = 10,
    target_cosine = 0.3,
    nonnegative = TRUE,
    precision_shift = 0.5,
    precision_scale = 0.3,
    graph_model = "erdos_renyi"
  )
})
ref <- rescale_single_cell_moments(moments, c(1, 2, 5))
sim <- simulate_bulk_from_cells(
  ref$generative_mean,
  ref$generative_sigma,
  p = c(0.6, 0.3, 0.1),
  n_cells = 1e4,
  platform_scale = 1,
  n_samples = 1L,
  seed = 4L
)
y <- sim$Y[, 1L]
mu <- ref$mean_signature_matrix
sigma <- ref$Sigma

loglik <- function(z) profile_loglik_single_cell_ilr(z, y, mu, sigma)
score <- function(z) gradient_profile_loglik_single_cell_ilr(z, y, mu, sigma)
hessian <- function(z) hessian_profile_loglik_single_cell_ilr(z, y, mu, sigma)

# --- Run nlminb from a far start and record what it evaluates ------------

# A start near the corner p = (0.02, 0.02, 0.96), far from the maximum, so
# that the radius has to grow, overshoot and shrink.
z_start <- isometric_log_ratio(c(0.02, 0.02, 0.96))
events <- list()
record <- function(type, z) {
  events[[length(events) + 1L]] <<- list(type = type, z = z)
}
fit <- stats::nlminb(
  z_start,
  objective = function(z) {
    record("f", z)
    -loglik(z)
  },
  gradient = function(z) {
    record("g", z)
    -score(z)
  },
  hessian = function(z) -hessian(z),
  control = list(rel.tol = 1e-8, x.tol = 1e-8, xf.tol = 1e-8, abs.tol = 0)
)
stopifnot(fit$convergence == 0L)

# Every objective call after the first is a trial step from the latest
# point where the gradient was requested (the current iterate). A trial is
# accepted when the next gradient request is at that trial. A trial that
# increased the log-likelihood but was not accepted was superseded: PORT
# enlarged the radius, tried a longer step from the same iterate, and
# restored the best step it had saved.
is_gradient <- vapply(events, function(e) e$type == "g", logical(1))
steps <- list()
for (i in which(!is_gradient)[-1L]) {
  base_event <- max(which(is_gradient[seq_len(i - 1L)]))
  z0 <- events[[base_event]]$z
  z1 <- events[[i]]$z
  delta <- z1 - z0
  if (sqrt(sum(delta^2)) < 1e-10) {
    next
  }
  later <- which(is_gradient & seq_along(events) > i)
  accepted <- length(later) > 0L &&
    max(abs(events[[later[[1L]]]]$z - z1)) < 1e-12
  g <- score(z0)
  h <- hessian(z0)
  predicted <- sum(g * delta) + 0.5 * drop(crossprod(delta, h %*% delta))
  actual <- loglik(z1) - loglik(z0)
  outcome <- if (accepted) {
    "accepted"
  } else if (actual > 0) {
    "superseded"
  } else {
    "rejected"
  }
  steps[[length(steps) + 1L]] <- data.frame(
    trial = length(steps) + 1L,
    base_1 = z0[[1L]],
    base_2 = z0[[2L]],
    rho_1 = z1[[1L]],
    rho_2 = z1[[2L]],
    length = sqrt(sum(delta^2)),
    ratio = actual / predicted,
    outcome = outcome
  )
}
steps <- do.call(rbind, steps)
outcome_levels <- c("accepted", "superseded", "rejected")
steps$outcome <- factor(steps$outcome, outcome_levels)
iterates <- unique(do.call(
  rbind,
  lapply(events[is_gradient], function(e) e$z)
))
iterate_df <- data.frame(rho_1 = iterates[, 1L], rho_2 = iterates[, 2L])
message(
  nrow(steps),
  " trials, ",
  sum(steps$outcome == "accepted"),
  " accepted; nlminb iterations: ",
  fit$iterations
)

# --- Contours of the profile log-likelihood ------------------------------

contour_grid <- function(x_range, y_range, n = 110L) {
  grid <- expand.grid(
    rho_1 = seq(x_range[[1L]], x_range[[2L]], length.out = n),
    rho_2 = seq(y_range[[1L]], y_range[[2L]], length.out = n)
  )
  grid$loglik <- vapply(
    seq_len(nrow(grid)),
    function(i) loglik(c(grid$rho_1[[i]], grid$rho_2[[i]])),
    numeric(1)
  )
  grid$loglik[!is.finite(grid$loglik)] <- NA_real_
  grid
}

circle_df <- function(steps_subset) {
  angle <- seq(0, 2 * pi, length.out = 120L)
  do.call(
    rbind,
    lapply(seq_len(nrow(steps_subset)), function(i) {
      s <- steps_subset[i, ]
      data.frame(
        trial = s$trial,
        outcome = s$outcome,
        rho_1 = s$base_1 + s$length * cos(angle),
        rho_2 = s$base_2 + s$length * sin(angle)
      )
    })
  )
}

path_panel <- function(
  x_range,
  y_range,
  steps_subset,
  title,
  label_size,
  legend = "bottom"
) {
  grid <- contour_grid(x_range, y_range)
  maximum <- max(grid$loglik, na.rm = TRUE)
  breaks <- maximum - c(0.5, 2, 5, 10, 20, 40, 80, 160, 320)
  ggplot() +
    geom_contour(
      data = grid,
      aes(rho_1, rho_2, z = loglik),
      breaks = breaks,
      colour = "#b8b2a7",
      linewidth = 0.35
    ) +
    geom_path(
      data = circle_df(steps_subset),
      aes(rho_1, rho_2, group = trial, colour = outcome, linetype = outcome),
      linewidth = 0.45
    ) +
    geom_segment(
      data = steps_subset,
      aes(
        x = base_1,
        y = base_2,
        xend = rho_1,
        yend = rho_2,
        colour = outcome
      ),
      linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")
    ) +
    geom_point(
      data = iterate_df,
      aes(rho_1, rho_2),
      colour = colour_ink,
      size = 1.8
    ) +
    geom_text(
      data = steps_subset[steps_subset$length > 0.02, ],
      aes(rho_1, rho_2, label = trial, colour = outcome),
      size = label_size,
      vjust = -0.7,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = outcome_colours,
      name = "trial step",
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = c(accepted = "solid", superseded = "dotdash", rejected = "dashed"),
      name = "trial step",
      drop = FALSE
    ) +
    coord_equal(xlim = x_range, ylim = y_range, expand = FALSE) +
    labs(
      title = title,
      x = "ILR coordinate 1",
      y = "ILR coordinate 2"
    ) +
    theme_slide() +
    theme(legend.position = legend)
}

full_range_1 <- range(c(steps$rho_1, steps$base_1)) + c(-1, 1)
full_range_2 <- range(c(steps$rho_2, steps$base_2)) + c(-1, 1)
panel_a <- path_panel(
  full_range_1,
  full_range_2,
  steps,
  "A  Trial steps and their balls",
  3.2
)

late <- steps[steps$length < 1.5 & steps$trial > 3L, ]
optimum <- fit$par
panel_b <- path_panel(
  optimum[[1L]] + c(-0.45, 0.45),
  optimum[[2L]] + c(-0.75, 0.6),
  late,
  "B  Zoom near the maximum",
  3.2,
  legend = "none"
)

panel_c <- ggplot(steps, aes(trial, length, colour = outcome)) +
  geom_line(aes(group = 1L), colour = "#b8b2a7", linewidth = 0.4) +
  geom_point(size = 2.4) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = colour_ink) +
  annotate(
    "text",
    x = max(steps$trial),
    y = 1.6,
    label = "first radius = 1",
    hjust = 1,
    size = 3.2,
    colour = colour_ink
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = outcome_colours,
    guide = "none"
  ) +
  labs(
    title = "C  Step length",
    x = "trial",
    y = "step length (Euclidean, ILR)"
  ) +
  theme_slide()

ratio_floor <- -1
steps$ratio_plot <- pmax(steps$ratio, ratio_floor)
panel_d <- ggplot(steps, aes(trial, ratio_plot, colour = outcome)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = ratio_floor,
    ymax = 1e-4,
    fill = colour_reject,
    alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 1e-4,
    ymax = 0.1,
    fill = "#d9d3c7",
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = c(1e-4, 0.1, 1),
    linetype = c("dashed", "dashed", "dotted"),
    colour = colour_ink
  ) +
  annotate(
    "text",
    x = max(steps$trial) + 0.4,
    y = c(-0.5, 0.35, 1.25),
    label = c("reject, shrink", "accept, shrink", "model exact"),
    hjust = 1,
    size = 3.2,
    colour = colour_ink
  ) +
  geom_point(size = 2.4) +
  scale_colour_manual(
    values = outcome_colours,
    guide = "none"
  ) +
  coord_cartesian(ylim = c(ratio_floor, max(steps$ratio, na.rm = TRUE) + 0.3)) +
  labs(
    title = "D  Actual / predicted increase",
    subtitle = "values below -1 drawn at -1",
    x = "trial",
    y = "ratio"
  ) +
  theme_slide()

# --- Assemble with grid layout and write the SVG -------------------------

draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      3L,
      2L,
      widths = grid::unit(c(1, 1.15), "null"),
      heights = grid::unit(c(1.2, 1, 1), "null")
    )
  ))
  place <- function(plot, row, col) {
    print(plot, vp = grid::viewport(layout.pos.row = row, layout.pos.col = col))
  }
  place(panel_a, 1:3, 1L)
  place(panel_b, 1L, 2L)
  place(panel_c, 2L, 2L)
  place(panel_d, 3L, 2L)
}

dir.create(dirname(output_svg), recursive = TRUE, showWarnings = FALSE)
svglite::svglite(output_svg, width = 12, height = 12, bg = "white")
draw_figure()
grDevices::dev.off()

inkscape <- Sys.which("inkscape")
if (nzchar(inkscape)) {
  system2(
    inkscape,
    c(shQuote(output_svg), "--export-type=pdf", paste0("--export-filename=", shQuote(output_pdf)))
  )
} else {
  message("Inkscape not found: export ", output_svg, " to PDF by hand.")
}
print(steps[, c("trial", "length", "ratio", "outcome")], digits = 3)
sessioninfo::session_info(pkgs = "attached")
