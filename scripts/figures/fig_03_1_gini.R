###############################################################################
# Figure for vignette 03-1, #fig-gini-derivation (#sec-gini-entropy): the Gini
# coefficient of one gene across cell types, from the mean absolute
# difference of the type means to the share form of eq-gini. Five ggplot2
# panels, one per step, assembled with the step boxes in
# vignettes/figures/src/gini/gini.tex (LuaLaTeX, TikZ, awesomebox).
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_gini.R
#   cd vignettes/figures/src/gini && latexmk -lualatex gini.tex
#   pdftocairo -svg gini.pdf ../../fig_03_1_gini.svg
#   pdftocairo -png -r 200 -singlefile gini.pdf ../../fig_03_1_gini
#   latexmk -c gini.tex
#
# Notation (03-1, eq-gini): gene g fixed, J types j, l; type means (counts)
# mu_gj, mean over types mu_bar_g, shares u_gj = mu_gj / sum_l mu_gl, ranks k
# with mu_g(1) <= ... <= mu_g(J). Steps:
#   A  Delta_g = J^-2 sum_j sum_l |mu_gj - mu_gl|        (count units)
#   B  Gini_g = Delta_g / (2 mu_bar_g)                     (scale-free)
#   C  Gini_g = (J^2 mu_bar_g)^-1 sum_j (2 R_gj - J - 1) mu_gj, R_gj = rank of mu_gj
#   D  max Gini_g = (J - 1) / J; J Gini_g / (J - 1) in [0, 1]
#   E  Gini_g = (2J)^-1 sum_j sum_l |u_gj - u_gl| = J^-1 sum_j (2 R_gj - J - 1) u_gj
# Two genes over J = 8 types, integer type means, every mean >= 1 (so both
# pass the presence gate). Deterministic: no random numbers are drawn (the
# seed is set for completeness). Every identity is checked with stopifnot().
###############################################################################

library(ggplot2)
library(DeCovarT) # theme_decovart_facets(), the repository theme

set.seed(20261101)

out_dir <- file.path("vignettes", "figures", "src", "gini")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_low <- "#1d4e89"
col_high <- "#c45c26"
col_guide <- "#1c2b24"
col_ref <- "grey45"
gene_cols <- c(low = col_low, high = col_high)
txt_size <- 8 / .pt # 8 pt annotations

# ---- data and the five forms of the coefficient -------------------------------
J <- 8L
genes <- list(
  low = c(12, 15, 9, 14, 11, 13, 10, 16),
  high = c(1, 2, 1, 3, 1, 2, 4, 86)
)
stopifnot(
  all(lengths(genes) == J),
  all(unlist(genes) >= 1),
  all(unlist(genes) == round(unlist(genes)))
)

gini_steps <- function(mu) {
  J <- length(mu)
  k <- seq_len(J)
  w <- 2 * k - J - 1
  u <- mu / sum(mu)
  mu_bar <- mean(mu)
  delta <- sum(abs(outer(mu, mu, "-"))) / J^2 # A
  list(
    delta = delta,
    mu_bar = mu_bar,
    gini_b = delta / (2 * mu_bar), # B
    weighted_sum = sum(w * sort(mu)),
    gini_c = sum(w * sort(mu)) / (J^2 * mu_bar), # C
    gini_e_pairs = sum(abs(outer(u, u, "-"))) / (2 * J), # E, pairs
    gini_e_ranks = sum(w * sort(u)) / J, # E, ranks
    corrected = J * (delta / (2 * mu_bar)) / (J - 1) # D
  )
}

res <- lapply(genes, gini_steps)

for (g in names(genes)) {
  r <- res[[g]]
  r10 <- gini_steps(10 * genes[[g]])
  stopifnot(
    # B = C = E (both forms): one coefficient, four formulas
    isTRUE(all.equal(r$gini_b, r$gini_c)),
    isTRUE(all.equal(r$gini_b, r$gini_e_pairs)),
    isTRUE(all.equal(r$gini_b, r$gini_e_ranks)),
    # A depends on scale, B does not
    isTRUE(all.equal(r10$delta, 10 * r$delta)),
    isTRUE(all.equal(r10$gini_b, r$gini_b)),
    # shares are the mean-relative values divided by J
    isTRUE(all.equal(
      genes[[g]] / sum(genes[[g]]),
      (genes[[g]] / r$mu_bar) / J
    )),
    # D: corrected value in [0, 1], below the maximum before correction
    r$gini_b <= (J - 1) / J,
    r$corrected >= 0,
    r$corrected <= 1
  )
}

# D: the maximum (J - 1) / J is reached when one type carries all the mass
one_type <- c(rep(0, J - 1), 1)
r_one <- gini_steps(one_type)
stopifnot(
  isTRUE(all.equal(r_one$gini_b, (J - 1) / J)),
  isTRUE(all.equal(r_one$gini_c, (J - 1) / J)),
  isTRUE(all.equal(r_one$corrected, 1))
)

# printed values (exact rationals for these counts)
stopifnot(
  isTRUE(all.equal(res$low$mu_bar, 12.5)),
  isTRUE(all.equal(res$high$mu_bar, 12.5)),
  isTRUE(all.equal(res$low$delta, 21 / 8)),
  isTRUE(all.equal(res$high$delta, 77 / 4)),
  res$low$weighted_sum == 84,
  res$high$weighted_sum == 616,
  isTRUE(all.equal(res$low$gini_b, 21 / 200)),
  isTRUE(all.equal(res$high$gini_b, 77 / 100)),
  isTRUE(all.equal(res$low$corrected, 3 / 25)),
  isTRUE(all.equal(res$high$corrected, 22 / 25))
)

summary_tab <- data.frame(
  gene = names(res),
  mu_bar = vapply(res, `[[`, numeric(1), "mu_bar"),
  delta = vapply(res, `[[`, numeric(1), "delta"),
  weighted_sum = vapply(res, `[[`, numeric(1), "weighted_sum"),
  gini = vapply(res, `[[`, numeric(1), "gini_b"),
  corrected = vapply(res, `[[`, numeric(1), "corrected")
)
print(summary_tab, row.names = FALSE)

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# ---- shared theme and output --------------------------------------------------
theme_fig <- theme_decovart_facets(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7.5, colour = col_guide),
    axis.title = element_text(size = 8.5, colour = col_guide),
    plot.margin = margin(3, 5, 2, 3)
  )

panel_w <- 3.4 # in, every panel
panel_h <- 1.4

save_panel <- function(p, name, width = panel_w, height = panel_h) {
  ggsave(
    filename = file.path(out_dir, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}

lab_type_j <- expression("type " * italic(j))
lab_type_l <- expression("type " * italic("\u2113"))

df_types <- do.call(
  rbind,
  lapply(names(genes), function(g) {
    mu <- genes[[g]]
    data.frame(
      gene = g,
      j = seq_len(J),
      mu = mu,
      ratio = mu / mean(mu),
      u = mu / sum(mu)
    )
  })
)
df_types$gene <- factor(df_types$gene, levels = c("low", "high"))

# ---- A: pairwise absolute differences ----------------------------------------
df_a <- do.call(
  rbind,
  lapply(names(genes), function(g) {
    mu <- genes[[g]]
    d <- expand.grid(j = seq_len(J), l = seq_len(J))
    d$diff <- abs(mu[d$j] - mu[d$l])
    d$gene <- g
    d
  })
)
strip_a <- c(
  low = sprintf(
    "'low:  '*Delta[g]=='%s'",
    format(res$low$delta)
  ),
  high = sprintf(
    "'high:  '*Delta[g]=='%s'",
    format(res$high$delta)
  )
)
df_a$strip <- factor(
  strip_a[df_a$gene],
  levels = strip_a[c("low", "high")]
)
pa <- ggplot(df_a, aes(j, l, fill = gene, alpha = diff)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_wrap(~strip, labeller = label_parsed) +
  scale_fill_manual(values = gene_cols, guide = "none") +
  scale_alpha_continuous(
    range = c(0.04, 1),
    limits = c(0, 86),
    transform = "sqrt",
    breaks = c(1, 5, 20, 50, 85),
    name = expression(abs(mu[gj] - mu[g * "\u2113"]))
  ) +
  scale_x_continuous(breaks = c(1, 4, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1, 4, 8), expand = c(0, 0)) +
  labs(x = lab_type_j, y = lab_type_l) +
  theme_fig +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3, "mm"),
    strip.text = element_text(size = 8, margin = margin(2, 2, 2, 2)),
    panel.spacing = unit(0.9, "lines"),
    panel.grid = element_blank()
  )
save_panel(pa, "panel_a_pairwise")

# ---- B: means relative to the mean of the gene ---------------------------------
ann_b <- data.frame(
  x = 0.55,
  y = c(6.2, 4.9),
  gene = c("low", "high"),
  label = c(
    sprintf(
      "'low:  '*Delta[g]/(2*bar(mu)[g])=='%s / %s = %s'",
      fmt(res$low$delta, 3),
      fmt(2 * res$low$mu_bar, 0),
      fmt(res$low$gini_b)
    ),
    sprintf(
      "'high:  '*Delta[g]/(2*bar(mu)[g])=='%s / %s = %s'",
      fmt(res$high$delta, 2),
      fmt(2 * res$high$mu_bar, 0),
      fmt(res$high$gini_b)
    )
  )
)
pb <- ggplot(df_types, aes(factor(j), ratio, fill = gene)) +
  geom_col(
    position = position_dodge(width = 0.82),
    width = 0.78,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 1,
    colour = col_guide,
    linetype = "22",
    linewidth = 0.45
  ) +
  geom_text(
    data = ann_b,
    aes(x, y, label = label, colour = gene),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    size = txt_size
  ) +
  scale_fill_manual(values = gene_cols) +
  scale_colour_manual(values = gene_cols) +
  scale_y_continuous(breaks = c(0, 1, 3, 5, 7), limits = c(0, 7.2)) +
  labs(
    x = lab_type_j,
    y = expression(mu[gj] / bar(mu)[g])
  ) +
  theme_fig +
  theme(panel.grid.major.x = element_blank())
save_panel(pb, "panel_b_scale_free")

# ---- C: sorted means against their rank, with the weights ---------------------
df_c <- do.call(
  rbind,
  lapply(names(genes), function(g) {
    data.frame(gene = g, k = seq_len(J), mu_k = sort(genes[[g]]))
  })
)
df_c$gene <- factor(df_c$gene, levels = c("low", "high"))
ann_c <- data.frame(
  x = 0.85,
  y = c(78, 58),
  gene = c("low", "high"),
  label = c(
    sprintf(
      "'low:  %d / %d = %s'",
      res$low$weighted_sum,
      J^2 * res$low$mu_bar,
      fmt(res$low$gini_c)
    ),
    sprintf(
      "'high:  %d / %d = %s'",
      res$high$weighted_sum,
      J^2 * res$high$mu_bar,
      fmt(res$high$gini_c)
    )
  )
)
pc <- ggplot(df_c, aes(k, mu_k, colour = gene)) +
  geom_vline(
    xintercept = (J + 1) / 2,
    colour = col_guide,
    linetype = "22",
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  geom_text(
    data = ann_c,
    aes(x, y, label = label, colour = gene),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    size = txt_size
  ) +
  scale_colour_manual(values = gene_cols) +
  scale_x_continuous(
    breaks = seq_len(J),
    limits = c(0.7, J + 0.3),
    expand = c(0, 0),
    sec.axis = sec_axis(
      transform = ~.,
      breaks = seq_len(J),
      labels = 2 * seq_len(J) - J - 1,
      name = expression("weight " * 2 * italic(R)[italic(gj)] - italic(J) - 1)
    )
  ) +
  scale_y_continuous(breaks = c(0, 40, 80), limits = c(0, 90)) +
  labs(
    x = expression("rank " * italic(R)[italic(gj)]),
    y = expression("order statistic " * mu[g(R)])
  ) +
  theme_fig
save_panel(pc, "panel_c_ordering")

# ---- D: raw and corrected coefficient, and the maximum ------------------------
df_d <- data.frame(
  row = c("low", "high", "one"),
  raw = c(res$low$gini_b, res$high$gini_b, r_one$gini_b),
  corrected = c(res$low$corrected, res$high$corrected, r_one$corrected)
)
df_d$row <- factor(df_d$row, levels = c("one", "high", "low"))
df_d$label <- sprintf(
  "%s \u2192 %s",
  fmt(df_d$raw),
  c(fmt(df_d$corrected[1:2]), "1")
)
df_d$label_x <- c(df_d$corrected[1] + 0.035, df_d$raw[2:3] - 0.035)
df_d$label_hjust <- c(0, 1, 1)
row_cols <- c(low = col_low, high = col_high, one = col_ref)
pd <- ggplot(df_d, aes(y = row, colour = row)) +
  geom_vline(
    xintercept = (J - 1) / J,
    colour = col_guide,
    linetype = "22",
    linewidth = 0.45
  ) +
  geom_vline(xintercept = 1, colour = col_guide, linewidth = 0.45) +
  annotate(
    "text",
    x = (J - 1) / J - 0.015,
    y = 3.55,
    label = "max~(J-1)/J=='0.875'",
    parse = TRUE,
    hjust = 1,
    size = txt_size,
    colour = col_guide
  ) +
  geom_segment(
    aes(x = raw, xend = corrected, yend = row),
    linewidth = 0.55,
    arrow = arrow(length = unit(1.4, "mm"), type = "closed")
  ) +
  geom_point(
    aes(x = raw),
    shape = 21,
    fill = "white",
    size = 1.9,
    stroke = 0.8
  ) +
  geom_point(aes(x = corrected), size = 1.9) +
  geom_text(
    aes(x = label_x, label = label, hjust = label_hjust),
    size = txt_size
  ) +
  scale_colour_manual(values = row_cols) +
  scale_y_discrete(
    labels = c(low = "low-Gini", high = "high-Gini", one = "one type only")
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1.04),
    expand = c(0, 0.01)
  ) +
  coord_cartesian(ylim = c(0.6, 3.75), clip = "off") +
  labs(
    x = expression(
      Gini[g] * " (open)  " %->% "  " * J ~ Gini[g] / (J - 1) * " (filled)"
    ),
    y = NULL
  ) +
  theme_fig +
  theme(panel.grid.major.y = element_blank())
save_panel(pd, "panel_d_correction")

# ---- E: type shares ------------------------------------------------------------
ann_e <- data.frame(
  x = 0.55,
  y = c(0.79, 0.63),
  gene = c("low", "high"),
  label = c(
    sprintf(
      "'low:  '*Gini[g]=='%s'*',  corrected '*'%s'",
      fmt(res$low$gini_e_pairs),
      fmt(res$low$corrected)
    ),
    sprintf(
      "'high:  '*Gini[g]=='%s'*',  corrected '*'%s'",
      fmt(res$high$gini_e_pairs),
      fmt(res$high$corrected)
    )
  )
)
pe <- ggplot(df_types, aes(factor(j), u, fill = gene)) +
  geom_col(
    position = position_dodge(width = 0.82),
    width = 0.78,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = 1 / J,
    colour = col_guide,
    linetype = "22",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 7.45,
    y = 1 / J + 0.085,
    hjust = 1,
    label = "1/J=='0.125'",
    parse = TRUE,
    size = txt_size,
    colour = col_guide
  ) +
  geom_text(
    data = ann_e,
    aes(x, y, label = label, colour = gene),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    size = txt_size
  ) +
  scale_fill_manual(values = gene_cols) +
  scale_colour_manual(values = gene_cols) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75),
    labels = c("0", "0.25", "0.5", "0.75"),
    limits = c(0, 0.9)
  ) +
  labs(x = lab_type_j, y = expression(u[gj])) +
  theme_fig +
  theme(panel.grid.major.x = element_blank())
save_panel(pe, "panel_e_shares", height = 1.6) # taller box E (two-line formula)

utils::sessionInfo()
