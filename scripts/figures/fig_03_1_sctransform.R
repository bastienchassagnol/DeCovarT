###############################################################################
# Figure for vignette 03-1, #fig-sc-transform: the SCTransform (v2) workflow
# on simulated counts. Four ggplot2 panels, assembled with the notation in
# vignettes/figures/src/sctransform/sctransform.tex (LuaLaTeX, TikZ).
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_sctransform.R
#   cd vignettes/figures/src/sctransform && latexmk -lualatex sctransform.tex
#   pdftocairo -svg sctransform.pdf ../../fig_03_1_sctransform.svg
#
# Model (eq-fs-nb-null of 03-1): x_gc ~ NB(m_gc, psi_g),
# log m_gc = beta_0g + log s_c, library size as offset, no cell type in the
# null. Two cell types are simulated so that marker genes keep their
# between-type signal in the Pearson residuals. Seed: 20261101.
###############################################################################

library(ggplot2)
library(DeCovarT) # theme_decovart_facets(), the repository theme

set.seed(20261101)

out_dir <- file.path("vignettes", "figures", "src", "sctransform")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_series <- "#1d4e89"
col_second <- "#c45c26"
col_guide <- "#1c2b24"
col_bars <- "#b9b2a5"

# ---- simulate ---------------------------------------------------------------
n_genes <- 2000L
n_cells <- 2000L
cell_type <- factor(rep(c("type 1", "type 2"), each = n_cells / 2))
s_c <- exp(stats::rnorm(n_cells, mean = 0, sd = 0.5)) # size factors
beta0 <- stats::rnorm(n_genes, mean = -1.5, sd = 1.3) # log mean per unit depth
lfc <- numeric(n_genes)
markers <- sample.int(n_genes, 200L)
lfc[markers] <- sample(c(-1.2, 1.2), 200L, replace = TRUE)
# true NB dispersion grows with the mean (psi -> Poisson limit is psi = Inf)
psi_true <- 1 + 12 * (1 - exp(-exp(beta0)))

type2 <- as.integer(cell_type == "type 2")
log_m <- outer(beta0, rep(1, n_cells)) +
  outer(lfc, type2) +
  outer(rep(1, n_genes), log(s_c))
m_true <- exp(log_m)
counts <- matrix(
  stats::rnbinom(length(m_true), mu = m_true, size = rep(psi_true, n_cells)),
  nrow = n_genes
)
keep <- rowSums(counts > 0) >= 20 # detection filter before fitting
counts <- counts[keep, ]
is_marker <- seq_len(n_genes)[keep] %in% markers

# ---- fit M0 gene by gene: offset log s_c, intercept only ---------------------
# MLE of beta_0g with an offset and no covariate: exp(beta_0g) = sum x / sum s.
mu_hat_unit <- rowSums(counts) / sum(s_c)
m_hat <- outer(mu_hat_unit, s_c)
theta_hat <- vapply(
  seq_len(nrow(counts)),
  function(g) {
    th <- tryCatch(
      suppressWarnings(MASS::theta.ml(counts[g, ], m_hat[g, ], limit = 50)),
      error = function(e) NA_real_
    )
    min(as.numeric(th), 1e4)
  },
  numeric(1)
)
gene_mean <- rowMeans(counts)
gene_var <- apply(counts, 1, stats::var)

# empirical-Bayes regularisation across genes: robust smooth of log10 theta
# against log10 mean, fitted on genes whose theta is estimable (below the cap;
# capped genes are indistinguishable from Poisson and, with failed fits
# (theta = 0), are left out of the fit)
ok <- is.finite(theta_hat) & theta_hat > 0 & theta_hat < 1e4
x_log <- log10(gene_mean[ok])
y_log <- log10(theta_hat[ok])
smooth <- stats::lowess(x_log, y_log, f = 0.3)
smooth_fun <- function(x) {
  stats::approx(smooth$x, smooth$y, xout = x, rule = 2, ties = mean)$y
}
theta_reg <- 10^smooth_fun(log10(gene_mean))

# Pearson residuals with the regularised theta, clipped at +- sqrt(C/30)
clip <- sqrt(n_cells / 30)
resid <- (counts - m_hat) / sqrt(m_hat + m_hat^2 / theta_reg)
resid <- pmin(pmax(resid, -clip), clip)

# points drawn: a random subset keeps the vector file small; every fit above
# uses all genes and cells
genes_draw <- sort(sample.int(nrow(counts), 800L))
cells_draw <- sort(sample.int(n_cells, 600L))

# one marker gene, moderately expressed, for panels b and d
cand <- which(is_marker & gene_mean > 1 & gene_mean < 4)
g_show <- cand[1]

theme_fig <- theme_decovart_facets(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  )

save_panel <- function(p, name, width = 3.1, height = 2.6) {
  ggsave(
    filename = file.path(out_dir, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}

# ---- (a) mean-variance --------------------------------------------------------
grid_m <- 10^seq(log10(min(gene_mean)), log10(max(gene_mean)), length.out = 200)
theta_curve <- 10^smooth_fun(log10(grid_m))
curves <- rbind(
  data.frame(m = grid_m, v = grid_m, law = "Poisson"),
  data.frame(
    m = grid_m,
    v = grid_m + grid_m^2 / theta_curve,
    law = "NB, smoothed dispersion"
  )
)
pa <- ggplot(data.frame(m = gene_mean, v = gene_var)[genes_draw, ], aes(m, v)) +
  geom_point(colour = col_bars, size = 0.5, alpha = 0.7) +
  geom_line(
    data = curves,
    aes(m, v, linetype = law),
    colour = col_series,
    linewidth = 0.7
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_linetype_manual(values = c("solid", "22")) +
  labs(x = "gene mean", y = "gene variance") +
  theme_fig
save_panel(pa, "panel_a_mean_variance")

# ---- (b) depth regression with offset ------------------------------------------
df_b <- data.frame(s = s_c, x = counts[g_show, ], type = cell_type)[
  cells_draw,
]
line_b <- data.frame(s = sort(s_c), m = mu_hat_unit[g_show] * sort(s_c))
pb <- ggplot(df_b, aes(s, x)) +
  geom_jitter(
    aes(colour = type),
    width = 0,
    height = 0.15,
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(data = line_b, aes(s, m), colour = col_guide, linewidth = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = c(col_series, col_second)) +
  labs(x = "size factor", y = "UMI count") +
  theme_fig
save_panel(pb, "panel_b_depth")

# ---- (c) regularisation of the dispersion -------------------------------------
df_c <- data.frame(m = gene_mean, theta = theta_hat)[
  intersect(genes_draw, which(ok)),
]
line_c <- data.frame(m = 10^smooth$x, theta = 10^smooth$y)
line_c <- line_c[!duplicated(line_c$m), ]
line_c <- line_c[order(line_c$m), ]
pc <- ggplot(df_c, aes(m, theta)) +
  geom_point(colour = col_bars, size = 0.5, alpha = 0.7) +
  geom_line(data = line_c, colour = col_second, linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "gene mean", y = "dispersion estimate") +
  theme_fig
save_panel(pc, "panel_c_regularisation")

# ---- (d) Pearson residuals: no depth trend, type signal kept -------------------
df_d <- data.frame(s = s_c, r = resid[g_show, ], type = cell_type)[cells_draw, ]
pd <- ggplot(df_d, aes(s, r)) +
  geom_hline(yintercept = 0, colour = col_guide, linewidth = 0.4) +
  geom_point(aes(colour = type), size = 0.5, alpha = 0.5) +
  scale_x_log10() +
  scale_colour_manual(values = c(col_series, col_second)) +
  labs(x = "size factor", y = "Pearson residual") +
  theme_fig
save_panel(pd, "panel_d_residuals")

cat(
  "genes kept:",
  nrow(counts),
  "; shown gene mean:",
  round(gene_mean[g_show], 2),
  "; residual clip:",
  round(clip, 2),
  "\n"
)
utils::sessionInfo()
