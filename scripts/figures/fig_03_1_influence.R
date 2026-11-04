###############################################################################
# Figures for vignette 03-1, #sec-fs-leverage ("Leverage is not influence"):
# regression diagnostics for outlier cells (figure 1) and outlier genes
# (figure 2), after the plot.lm diagnostics which = 1, 2, 5 and 6.
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_influence.R
# Writes
#   vignettes/figures/fig_03_1_influence_cells.{svg,png}
#   vignettes/figures/fig_03_1_influence_genes.{svg,png}
# Seed: 20261101 (one stream for the simulation and the randomised quantile
# residuals; ggrepel label placement uses the same seed).
#
# Figure 1, cells (eq-fs-nb-main of 03-1). One gene, J = 4 types with
# C_j = 2000, 300, 60, 25 cells, size factors s_c ~ log-normal(0, 0.5),
#   x_c ~ NB(m_c, psi),  log m_c = beta_0 + beta_{j(c)} + log s_c,
#   Var(x_c) = m_c + m_c^2 / psi,  psi = 5,
# type means per unit depth exp(beta_0 + beta_j) = 1, 2, 0.5, 30 (type 4 is
# a rare type carrying a high-expression marker). Injected cells: two
# doublet-like cells of type 4 with counts 4 x their fitted mean, three cells
# of type 1 (size factor between 1 and 2) with counts 8 x their fitted mean.
# "Fitted" means the M1 fit that includes them: the counts are set by a
# fixed-point iteration (count = round(factor x fitted mean), refit).
# Fit: MASS::glm.nb, offset log s_c, psi estimated by maximum likelihood.
#   IRLS weight        w_c  = m_c^2 / V(m_c)
#   leverage           h_c  = w_c / sum_{c': j(c') = j(c)} w_c'
#   standardised Pearson r^SP_c = (x_c - m_c) / sqrt(phi V(m_c) (1 - h_c))
#   Cook's distance    D_c  = r_c^2 h_c / (J phi (1 - h_c)^2),
#                           r_c the Pearson residual,
# with phi = 1, the convention of cooks.distance() and rstandard() on a
# negbin fit (summary.negbin uses dispersion = 1); the Pearson estimate
# phi-hat = X^2 / (C - J) is printed for comparison.
# Randomised quantile residuals (Dunn and Smyth 1996):
#   u_c ~ U(F(x_c - 1), F(x_c)),  r^Q_c = Phi^{-1}(u_c),
# F the fitted NB cdf, computed through the upper tail for stability.
#
# Figure 2, genes. Bulk per cell regressed on the reference signature across
# G = 500 genes, J = 4 types, signature means mu_gj ~ log-normal(0, 1.2),
# true composition p = (0.5, 0.3, 0.15, 0.05):
#   y_g = b_g sum_j p_j mu_gj + eps_g,   eps_g ~ N(0, tau_g^2),
#   tau_g^2 = kappa * { S2_g / n + sum_j p_j^2 V(mu_gj) / C_j },
#   S2_g    = sum_j p_j V(mu_gj) + sum_j p_j mu_gj^2 - (sum_j p_j mu_gj)^2,
#   V(mu)   = mu + mu^2 / psi,  psi = 5,  n = 1e4 cells,  C_j as above.
# S2_g / n is the variance of the per-cell mean of n cells drawn from the NB
# cell mixture; the second term is the sampling variance of the reference
# means estimated from C_j cells. kappa = 100 stands for the between-sample
# variance missing from both terms; it sets the noise level relative to the
# discordance and does not change p-hat or h_g (it cancels in the weights).
# b_g = 4 for 10 platform-discordant genes, 1 otherwise. Two consistent genes
# have a mean over types 50 x the median: E1 with an ordinary profile, E2
# with a rare-type (type 4) profile. Fit: weighted least squares, weights
# 1 / tau_g^2 (true p, oracle weights), no intercept, no constraint.
#   h_g = w_g x_g' (X' W X)^{-1} x_g,
#   r_g = sqrt(w_g) (y_g - x_g' p-hat) / (sigma-hat sqrt(1 - h_g)),
#   D_g = e_g^2 w_g h_g / (J sigma-hat^2 (1 - h_g)^2),
#   sigma-hat^2 = sum_g w_g e_g^2 / (G - J),
# checked against lm(..., weights =), hatvalues(), rstandard(),
# cooks.distance().
# Because tau_g^2 grows like mu^2, scaling a gene's profile by 50 scales
# w_g by about 1/2500: gene leverage depends on the direction of the profile
# (how specific it is to a rare type), not on its level. The script prints
# the OLS leverage of E1 and E2 for comparison.
###############################################################################

library(ggplot2)
library(patchwork)
library(DeCovarT) # theme_decovart_facets(), the repository theme

seed <- 20261101
set.seed(seed)

out_dir <- file.path("vignettes", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_series <- "#1d4e89"
col_second <- "#c45c26"
col_guide <- "#1c2b24"
col_grey <- "#a7a9ac"
col_green <- "#5b8c5a"
r_levels <- c(1, 2, 3, 6) # lines of equal |standardised residual|, panel d

# ---- shared drawing helpers -------------------------------------------------
theme_fig <- theme_decovart_facets(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.tag = element_text(face = "bold", size = 12),
    plot.margin = margin(4, 8, 4, 4)
  )

lab_plain <- function(x) {
  vapply(
    x,
    function(v) {
      if (is.na(v)) {
        return(NA_character_)
      }
      format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
    },
    character(1)
  )
}

# randomised quantile residual for a discrete cdf, through the upper tail:
# 1 - u = S(y - 1) - v (S(y - 1) - S(y)), S(k) = P(Y > k)
rq_resid_nb <- function(y, mu, size) {
  s_lo <- stats::pnbinom(y - 1, mu = mu, size = size, lower.tail = FALSE)
  s_hi <- stats::pnbinom(y, mu = mu, size = size, lower.tail = FALSE)
  v <- stats::runif(length(y))
  stats::qnorm(s_lo - v * (s_lo - s_hi), lower.tail = FALSE)
}

label_layer <- function(xvar, yvar, size = 3) {
  ggrepel::geom_text_repel(
    data = function(d) d[d$label != "", ],
    aes(.data[[xvar]], .data[[yvar]], label = label),
    colour = col_guide,
    size = size,
    seed = seed,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.colour = col_guide,
    segment.size = 0.25,
    max.overlaps = Inf
  )
}

# (a) residual against fitted value, log x axis
panel_resid_fitted <- function(d, pts, x_lab, y_lab) {
  ggplot(d) +
    geom_hline(yintercept = 0, colour = col_guide, linewidth = 0.4) +
    geom_hline(
      yintercept = c(-3, 3),
      colour = col_guide,
      linewidth = 0.3,
      linetype = "22"
    ) +
    pts("fit", "res") +
    label_layer("fit", "res") +
    scale_x_log10(labels = lab_plain) +
    labs(x = x_lab, y = y_lab) +
    theme_fig
}

# (b) normal QQ plot; y = x is the N(0, 1) reference, the quartile line is the
# plot.lm / qqline reference
panel_qq <- function(d, pts, y_lab, quartile_line = FALSE) {
  p <- ggplot(d) +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = col_guide,
      linewidth = 0.4,
      linetype = if (quartile_line) "13" else "solid"
    )
  if (quartile_line) {
    qy <- stats::quantile(d$res, c(0.25, 0.75), names = FALSE)
    qx <- stats::qnorm(c(0.25, 0.75))
    slope <- diff(qy) / diff(qx)
    p <- p +
      geom_abline(
        intercept = qy[1] - slope * qx[1],
        slope = slope,
        colour = col_guide,
        linewidth = 0.4
      ) +
      annotate(
        "text",
        x = max(d$q_theo),
        y = min(d$res),
        label = "solid: quartile line\ndotted: y = x",
        hjust = 1,
        vjust = 0,
        lineheight = 0.95,
        size = 2.8,
        colour = col_guide
      )
  }
  p +
    pts("q_theo", "res") +
    label_layer("q_theo", "res") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
    labs(x = "N(0, 1) quantile", y = y_lab) +
    theme_fig
}

# (c) standardised residual against leverage, Cook contours
#   r = +- sqrt(D J (1 - h) / h),  D = 0.5 (dashed) and 1 (solid)
panel_resid_leverage <- function(d, pts, J, h_lim, r_lim, y_lab) {
  h <- 10^seq(log10(h_lim[1]), log10(h_lim[2]), length.out = 300)
  cont <- do.call(
    rbind,
    lapply(c(0.5, 1), function(D) {
      r <- sqrt(D * J * (1 - h) / h)
      rbind(
        data.frame(h = h, r = r, D = D, side = "upper"),
        data.frame(h = h, r = -r, D = D, side = "lower")
      )
    })
  )
  # labels in the top band: Cook 0.5 left of where its curve leaves the
  # panel, Cook 1 right of it (h* solves r(h*) = upper y limit)
  r_top <- r_lim[2] - 0.05 * diff(r_lim)
  h_star <- c(0.5, 1) * J / (r_top^2 + c(0.5, 1) * J)
  lab <- data.frame(
    h = h_star * c(1 / 1.15, 1.25),
    r = r_top,
    D = c(0.5, 1),
    hjust = c(1, 0)
  )
  ggplot(d) +
    geom_hline(yintercept = 0, colour = col_guide, linewidth = 0.4) +
    geom_line(
      data = cont,
      aes(h, r, group = interaction(D, side), linetype = factor(D)),
      colour = col_guide,
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    geom_text(
      data = lab,
      aes(h, r, label = paste("Cook", D), hjust = hjust),
      colour = col_guide,
      size = 2.8
    ) +
    pts("h", "r_sp") +
    label_layer("h", "r_sp") +
    scale_linetype_manual(values = c("0.5" = "22", "1" = "solid")) +
    scale_x_log10(labels = lab_plain) +
    coord_cartesian(xlim = h_lim, ylim = r_lim, expand = FALSE) +
    labs(x = "leverage h (log scale)", y = y_lab) +
    theme_fig
}

# (d) Cook's distance against h / (1 - h), log-log; dotted lines of equal
# |standardised residual|: D = r^2 x / J; dashed: the 4/n cut. Cook values
# below the floor are drawn on the lower edge.
panel_cook_leverage <- function(d, pts, J, x_lim, cook_lim, n_obs, n_lab) {
  x <- 10^seq(log10(x_lim[1]), log10(x_lim[2]), length.out = 100)
  cont <- do.call(
    rbind,
    lapply(r_levels, function(r) data.frame(x = x, D = r^2 * x / J, r = r))
  )
  x_lab <- x_lim[2] / 1.08
  lab <- data.frame(x = x_lab, D = r_levels^2 * x_lab / J, r = r_levels)
  lab <- lab[lab$D < cook_lim[2], ]
  lab$text <- as.character(lab$r)
  lab$text[nrow(lab)] <- paste0("|r| = ", lab$r[nrow(lab)])
  d$cook_plot <- pmax(d$cook, cook_lim[1])
  ggplot(d) +
    geom_line(
      data = cont,
      aes(x, D, group = r),
      colour = col_guide,
      linewidth = 0.35,
      linetype = "13"
    ) +
    geom_text(
      data = lab,
      aes(x, D, label = text),
      colour = col_guide,
      size = 2.7,
      hjust = 1,
      vjust = -0.4
    ) +
    geom_hline(
      yintercept = 4 / n_obs,
      colour = col_guide,
      linewidth = 0.3,
      linetype = "42"
    ) +
    annotate(
      "text",
      x = x_lim[1] * 1.1,
      y = 4 / n_obs,
      label = n_lab,
      hjust = 0,
      vjust = -0.4,
      size = 2.8,
      colour = col_guide
    ) +
    pts("h_ratio", "cook_plot") +
    label_layer("h_ratio", "cook_plot") +
    scale_x_log10(labels = lab_plain) +
    scale_y_log10(labels = scales::label_log()) +
    coord_cartesian(xlim = x_lim, ylim = cook_lim, expand = FALSE) +
    labs(x = "h / (1 - h) (log scale)", y = "Cook's distance (log scale)") +
    theme_fig
}

assemble <- function(pa, pb, pc, pd) {
  (pa + pb) /
    (pc + pd) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")
}

save_figure <- function(p, name, width = 8, height = 7.6) {
  svglite::svglite(
    file.path(out_dir, paste0(name, ".svg")),
    width = width,
    height = height,
    bg = "white"
  )
  print(p)
  grDevices::dev.off()
  ragg::agg_png(
    file.path(out_dir, paste0(name, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 200,
    background = "white"
  )
  print(p)
  grDevices::dev.off()
  invisible(NULL)
}

###############################################################################
# Figure 1: cells
###############################################################################

C_j <- c(2000L, 300L, 60L, 25L)
J <- length(C_j)
n_cells <- sum(C_j)
psi <- 5
mean_unit <- c(1, 2, 0.5, 30)
type_lab <- sprintf("type %d (%d cells)", seq_len(J), C_j)
cell_type <- factor(rep(type_lab, C_j), levels = type_lab)
type_id <- as.integer(cell_type)
s_c <- exp(stats::rnorm(n_cells, mean = 0, sd = 0.5))
m_true <- mean_unit[type_id] * s_c
x <- stats::rnbinom(n_cells, mu = m_true, size = psi)

# injected cells
id_doublet <- sort(sample(which(type_id == 4L), 2L))
id_outlier <- sort(sample(which(type_id == 1L & s_c > 1 & s_c < 2), 3L))
id_inj <- c(id_doublet, id_outlier)
factor_inj <- c(rep(4, 2L), rep(8, 3L))
lab_inj <- c("D1", "D2", "O1", "O2", "O3")

cells <- data.frame(x = x, type = cell_type, s = s_c)
fit_m1 <- function(d) {
  MASS::glm.nb(x ~ type + offset(log(s)), data = d)
}
theta_clean <- fit_m1(cells)$theta # before injection, for comparison
cells$x[id_inj] <- round(factor_inj * m_true[id_inj])
for (it in seq_len(50L)) {
  fit <- fit_m1(cells)
  x_new <- round(factor_inj * stats::fitted(fit)[id_inj])
  if (all(x_new == cells$x[id_inj])) {
    break
  }
  cells$x[id_inj] <- x_new
}
fit <- fit_m1(cells)
stopifnot(all(
  cells$x[id_inj] == round(factor_inj * stats::fitted(fit)[id_inj])
))

theta_hat <- fit$theta
m_hat <- unname(stats::fitted(fit))
V_hat <- m_hat + m_hat^2 / theta_hat
w <- m_hat^2 / V_hat
h <- w / stats::ave(w, type_id, FUN = sum)
r_p <- (cells$x - m_hat) / sqrt(V_hat)
phi_one <- 1
phi_pearson <- sum(r_p^2) / (n_cells - J)
r_sp <- r_p / sqrt(phi_one * (1 - h))
cook <- r_p^2 * h / (J * phi_one * (1 - h)^2)
r_q <- rq_resid_nb(cells$x, m_hat, theta_hat)

# checks against the fitted object. fit$weights are the working weights of
# the last IRLS step, so they match m-hat^2 / V(m-hat) only up to the
# convergence tolerance (relative error about 1e-7); with fit$weights the
# block formula for h is exact to machine precision.
w_irls <- unname(fit$weights)
stopifnot(
  isTRUE(all.equal(w_irls, w, tolerance = 1e-6)),
  isTRUE(all.equal(
    unname(stats::hatvalues(fit)),
    w_irls / stats::ave(w_irls, type_id, FUN = sum),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(unname(stats::hatvalues(fit)), h, tolerance = 1e-6)),
  isTRUE(all.equal(
    unname(stats::residuals(fit, type = "pearson")),
    r_p,
    tolerance = 1e-6
  )),
  isTRUE(all.equal(
    unname(stats::rstandard(fit, type = "pearson")),
    r_sp,
    tolerance = 1e-6
  )),
  isTRUE(all.equal(
    unname(stats::cooks.distance(fit)),
    cook,
    tolerance = 1e-6
  )),
  summary(fit)$dispersion == 1
)
# the Pearson-dispersion version differs from cooks.distance() by 1 / phi-hat
stopifnot(
  !isTRUE(all.equal(
    unname(stats::cooks.distance(fit)),
    cook / phi_pearson,
    tolerance = 1e-6
  ))
)

cells$fit <- m_hat
cells$res <- r_q
cells$h <- h
cells$h_ratio <- h / (1 - h)
cells$r_sp <- r_sp
cells$cook <- cook
cells$label <- ""
cells$label[id_inj] <- lab_inj
cells$injected <- "consistent cell"
cells$injected[id_doublet] <- "doublet-like, type 4 (D)"
cells$injected[id_outlier] <- "outlier, type 1 (O)"
cells$injected <- factor(
  cells$injected,
  levels = c(
    "consistent cell",
    "doublet-like, type 4 (D)",
    "outlier, type 1 (O)"
  )
)
cells$q_theo <- NA_real_
cells$q_theo[order(cells$res)] <- stats::qnorm(stats::ppoints(n_cells))
# drawing order: large types first, injected cells last
cells <- cells[order(type_id), ]
cells <- rbind(cells[cells$label == "", ], cells[cells$label != "", ])

type_cols <- stats::setNames(
  c(col_grey, col_green, col_series, col_second),
  type_lab
)
cell_points <- function(xvar, yvar) {
  list(
    geom_point(
      data = function(d) d[d$label == "", ],
      aes(.data[[xvar]], .data[[yvar]], colour = type),
      shape = 16,
      size = 0.9,
      alpha = 0.55
    ),
    geom_point(
      data = function(d) d[d$label != "", ],
      aes(.data[[xvar]], .data[[yvar]], fill = type, shape = injected),
      colour = col_guide,
      size = 2.6,
      stroke = 0.5
    ),
    scale_colour_manual(values = type_cols, drop = FALSE),
    scale_fill_manual(values = type_cols, guide = "none"),
    scale_shape_manual(
      values = c("doublet-like, type 4 (D)" = 24, "outlier, type 1 (O)" = 22),
      drop = TRUE
    ),
    guides(
      colour = guide_legend(
        nrow = 2,
        override.aes = list(size = 2.6, alpha = 1),
        order = 1
      ),
      shape = guide_legend(
        nrow = 2,
        override.aes = list(fill = c(col_second, col_grey), size = 2.6),
        order = 2
      )
    )
  )
}

fig_cells <- assemble(
  panel_resid_fitted(
    cells,
    cell_points,
    x_lab = "fitted mean (log scale)",
    y_lab = "randomised quantile residual"
  ),
  panel_qq(cells, cell_points, y_lab = "randomised quantile residual"),
  panel_resid_leverage(
    cells,
    cell_points,
    J,
    h_lim = c(min(cells$h) / 1.3, 0.2),
    r_lim = c(-2.5, 9.5),
    y_lab = "standardised Pearson residual"
  ),
  panel_cook_leverage(
    cells,
    cell_points,
    J,
    x_lim = c(min(cells$h_ratio) / 1.3, 0.25),
    cook_lim = c(1e-6, 6),
    n_obs = n_cells,
    n_lab = "4/C"
  )
)
save_figure(fig_cells, "fig_03_1_influence_cells")

###############################################################################
# Figure 2: genes
###############################################################################

G <- 500L
n_bulk <- 1e4
kappa <- 100
p_true <- c(0.5, 0.3, 0.15, 0.05)
V_nb <- function(m) m + m^2 / psi
mu <- matrix(stats::rlnorm(G * J, meanlog = 0, sdlog = 1.2), G, J)
med_mean <- stats::median(rowMeans(mu))

# two consistent genes with a mean over types 50 x the median:
# E1 keeps an ordinary profile, E2 is a type-4 marker
id_ext <- 1:2
prof_e1 <- mu[1, ] / mean(mu[1, ])
prof_e2 <- c(0.02, 0.02, 0.02, 1) / mean(c(0.02, 0.02, 0.02, 1))
mu[1, ] <- 50 * med_mean * prof_e1
mu[2, ] <- 50 * med_mean * prof_e2
id_disc <- sort(sample(setdiff(seq_len(G), id_ext), 10L))
b <- rep(1, G)
b[id_disc] <- 4

S2 <- drop(V_nb(mu) %*% p_true) +
  drop(mu^2 %*% p_true) -
  drop(mu %*% p_true)^2
tau2 <- kappa * (S2 / n_bulk + drop(V_nb(mu) %*% (p_true^2 / C_j)))
y_bulk <- b * drop(mu %*% p_true) + stats::rnorm(G, sd = sqrt(tau2))

w_g <- 1 / tau2
Xw <- mu * sqrt(w_g)
XtWX_inv <- solve(crossprod(Xw))
p_hat <- drop(XtWX_inv %*% crossprod(Xw, y_bulk * sqrt(w_g)))
# the unconstrained fit is interior (every p-hat_j > 0), so the non-negative
# least-squares solution coincides with it
stopifnot(all(p_hat > 0))
y_hat <- drop(mu %*% p_hat)
e_g <- y_bulk - y_hat
h_g <- rowSums((Xw %*% XtWX_inv) * Xw)
sigma2_hat <- sum(w_g * e_g^2) / (G - J)
r_std_g <- sqrt(w_g) * e_g / sqrt(sigma2_hat * (1 - h_g))
cook_g <- e_g^2 * w_g * h_g / (J * sigma2_hat * (1 - h_g)^2)

fit_lm <- stats::lm(y_bulk ~ 0 + mu, weights = w_g)
stopifnot(
  isTRUE(all.equal(unname(stats::coef(fit_lm)), p_hat, tolerance = 1e-8)),
  isTRUE(all.equal(unname(stats::hatvalues(fit_lm)), h_g, tolerance = 1e-8)),
  isTRUE(all.equal(
    unname(stats::rstandard(fit_lm)),
    r_std_g,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(
    unname(stats::cooks.distance(fit_lm)),
    cook_g,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(sum(h_g), J))
)
# leverage the extreme genes would have without weights (ordinary least squares)
h_ols <- rowSums((mu %*% solve(crossprod(mu))) * mu)

share4 <- mu[, 4] / rowSums(mu)
class_lev <- c(
  "consistent gene",
  "consistent, type-4 marker (share \u2265 0.8)",
  "platform-discordant, b = 4 (d)",
  "consistent, mean 50 \u00d7 median (E)"
)
gene_class <- rep(class_lev[1], G)
gene_class[share4 >= 0.8] <- class_lev[2]
gene_class[id_disc] <- class_lev[3]
gene_class[id_ext] <- class_lev[4]
genes <- data.frame(
  fit = y_hat,
  res = r_std_g,
  h = h_g,
  h_ratio = h_g / (1 - h_g),
  r_sp = r_std_g,
  cook = cook_g,
  class = factor(gene_class, levels = class_lev),
  label = ""
)
genes$label[id_ext] <- c("E1", "E2")
genes$label[id_disc] <- paste0("d", seq_along(id_disc))
genes$q_theo <- NA_real_
genes$q_theo[order(genes$res)] <- stats::qnorm(stats::ppoints(G))
genes <- genes[order(match(genes$class, class_lev)), ]

class_cols <- stats::setNames(
  c(col_grey, col_green, col_second, col_series),
  class_lev
)
class_shapes <- stats::setNames(c(16, 16, 17, 15), class_lev)
class_sizes <- stats::setNames(c(1.2, 1.7, 2.3, 2.6), class_lev)
gene_points <- function(xvar, yvar) {
  list(
    geom_point(
      aes(
        .data[[xvar]],
        .data[[yvar]],
        colour = class,
        shape = class,
        size = class
      ),
      alpha = 0.85
    ),
    scale_colour_manual(values = class_cols, drop = FALSE),
    scale_shape_manual(values = class_shapes, drop = FALSE),
    scale_size_manual(values = class_sizes, drop = FALSE),
    guides(colour = guide_legend(nrow = 2))
  )
}

fig_genes <- assemble(
  panel_resid_fitted(
    genes,
    gene_points,
    x_lab = "fitted bulk value (log scale)",
    y_lab = "standardised residual"
  ),
  panel_qq(
    genes,
    gene_points,
    y_lab = "standardised residual",
    quartile_line = TRUE
  ),
  panel_resid_leverage(
    genes,
    gene_points,
    J,
    h_lim = c(min(genes$h) / 1.3, 0.12),
    r_lim = c(-2.5, 12),
    y_lab = "standardised residual"
  ),
  panel_cook_leverage(
    genes,
    gene_points,
    J,
    x_lim = c(min(genes$h_ratio) / 1.3, 0.14),
    cook_lim = c(1e-7, 3),
    n_obs = G,
    n_lab = "4/G"
  )
)
save_figure(fig_genes, "fig_03_1_influence_genes")

###############################################################################
# Numbers for the caption
###############################################################################

cat("\n==== Figure 1: cells ====\n")
cat("fixed-point iterations:", it, "\n")
cat(
  "theta-hat:",
  round(theta_hat, 3),
  "(before injection:",
  round(theta_clean, 3),
  "); Pearson phi-hat:",
  round(phi_pearson, 3),
  "\n"
)
cat(
  "cooks.distance() matches phi = 1 (summary(fit)$dispersion =",
  summary(fit)$dispersion,
  "); with the Pearson phi-hat every Cook is divided by",
  round(phi_pearson, 3),
  "\n"
)
print(stats::aggregate(
  h ~ type,
  data = cells,
  FUN = function(z) {
    signif(c(min = min(z), median = stats::median(z), max = max(z)), 3)
  }
))
idx <- match(lab_inj, cells$label)
inj_tab <- data.frame(
  label = lab_inj,
  type = as.character(cells$type[idx]),
  count = cells$x[idx],
  fitted = round(cells$fit[idx], 2),
  true_mean = round(m_true[id_inj], 2),
  h = signif(cells$h[idx], 3),
  r_sp = round(cells$r_sp[idx], 2),
  r_q = round(cells$res[idx], 2),
  cook = signif(cells$cook[idx], 3),
  cook_rank = rank(-cells$cook)[idx]
)
print(inj_tab)
cons <- cells$label == ""
t4 <- cons & as.integer(cells$type) == 4L
cat(
  "consistent type-4 cells: max Cook",
  signif(max(cells$cook[t4]), 3),
  "; max |r_sp|",
  round(max(abs(cells$r_sp[t4])), 2),
  "\n"
)
cat(
  "max Cook, consistent cells of types 1-3:",
  signif(max(cells$cook[cons & !t4]), 3),
  "\n"
)
tid <- as.integer(cells$type)
cat(
  "cells with Cook > 0.5:",
  sum(cells$cook > 0.5),
  "; > 1:",
  sum(cells$cook > 1),
  "; > 4/C =",
  signif(4 / n_cells, 3),
  ":",
  sum(cells$cook > 4 / n_cells),
  "(by type:",
  paste(tabulate(tid[cells$cook > 4 / n_cells], J), collapse = ", "),
  "; injected:",
  sum(cells$cook[!cons] > 4 / n_cells),
  ")\n"
)
cat("cells drawn on the Cook floor of panel d:", sum(cells$cook < 1e-6), "\n")
cat(
  "|r_q| > 3:",
  sum(abs(cells$res) > 3),
  "(injected:",
  sum(abs(cells$res[!cons]) > 3),
  ")\n"
)

cat("\n==== Figure 2: genes ====\n")
cat(
  "p-hat:",
  round(p_hat, 4),
  "; normalised:",
  round(p_hat / sum(p_hat), 4),
  "\n"
)
cat("sigma-hat^2:", round(sigma2_hat, 3), "\n")
cat(
  "median CV tau_g / (sum_j p_j mu_gj):",
  round(stats::median(sqrt(tau2) / drop(mu %*% p_true)), 3),
  "\n"
)
cat(
  "leverage: median",
  signif(stats::median(h_g), 3),
  "; max",
  signif(max(h_g), 3),
  "; mean J/G =",
  J / G,
  "\n"
)
mk <- gene_class == class_lev[2]
cat(
  "type-4 markers (share >= 0.8, consistent):",
  sum(mk),
  "; h range",
  signif(range(h_g[mk]), 3),
  "; max Cook",
  signif(max(cook_g[mk]), 3),
  "\n"
)
print(data.frame(
  label = c("E1", "E2"),
  h_wls = signif(h_g[id_ext], 3),
  h_ols = signif(h_ols[id_ext], 3),
  r = round(r_std_g[id_ext], 2),
  cook = signif(cook_g[id_ext], 3),
  share4 = round(share4[id_ext], 3)
))
print(data.frame(
  label = paste0("d", seq_along(id_disc)),
  h = signif(h_g[id_disc], 3),
  r = round(r_std_g[id_disc], 2),
  cook = signif(cook_g[id_disc], 3),
  share4 = round(share4[id_disc], 3),
  cook_rank = rank(-cook_g)[id_disc]
))
cat(
  "genes with Cook > 0.5:",
  sum(cook_g > 0.5),
  "; > 1:",
  sum(cook_g > 1),
  "; > 4/G =",
  4 / G,
  ":",
  sum(cook_g > 4 / G),
  "(discordant among them:",
  sum(cook_g[id_disc] > 4 / G),
  ")\n"
)
cat(
  "|r| > 3:",
  sum(abs(r_std_g) > 3),
  "(discordant:",
  sum(abs(r_std_g[id_disc]) > 3),
  ")\n"
)
cat("genes drawn on the Cook floor of panel d:", sum(cook_g < 1e-7), "\n")

utils::sessionInfo()
