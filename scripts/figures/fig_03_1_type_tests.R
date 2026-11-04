###############################################################################
# Figure for vignette 03-1, #fig-type-tests: four tests of the cell-type
# effect on single-cell counts, compared on simulated genes whose null is
# either homoscedastic or heteroscedastic across types.
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_type_tests.R
# Output: vignettes/figures/fig_03_1_type_tests.{svg,png}. Run time: about
# 2.2 min on one core (reference BLAS), printed at the end.
#
# Simulation. J = 4 types with C_j = 2000, 300, 60, 25 cells (C = 2385),
# size factors s_c ~ log-normal(0, 0.5^2), G = 1000 genes with baseline
# log mean per unit depth beta_0g ~ N(-1, 1.2^2), and
#   x_gc ~ NB(m_gc, psi), Var(x_gc) = m_gc + m_gc^2 / psi,
#   log m_gc = beta_0g + lfc_{g, j(c)} + log s_c.
# - genes 1-600, homoscedastic null: lfc = 0, one psi_g, log-uniform on
#   [2, 20];
# - genes 601-900, heteroscedastic null: lfc = 0, psi_gj = 20 in types 1 and
#   2 and psi_gj = 1 in types 3 and 4: the type changes the shape of the
#   count distribution, not its mean, and breaks the variance function that
#   the parametric tests assume;
# - genes 901-1000, alternatives: |lfc| ~ U(0.5, 1.5), random sign, in one
#   or two randomly chosen types; psi_g as for the homoscedastic null.
# No presence gate (eq-fs-presence) is applied: every gene is tested; the
# number of genes the gate would remove is printed.
#
# Models (03-1): M0 (eq-fs-nb-null), log m_gc = beta_0g + log s_c;
# M1 (eq-fs-nb-main), M0 + type main effect, J - 1 = 3 df.
# Tests, one per gene:
# 1. LR, analysis of deviance: M1 by MASS::glm.nb (offset log s_c), M0
#    refitted with the same psi-hat; Lambda = D0 - D1 against chi^2_{J-1}.
# 2. One-way ANOVA F on the shifted logarithm log(x / s + y0),
#    y0 = 1 / (4 alpha), alpha = 0.05 fixed (Ahlmann-Eltze & Huber 2023),
#    against F_{J-1, C-J}.
# 3. Quasi-likelihood F (eq-fs-qlf, deviance numerator as in edgeR and
#    glmGamPoi::test_de()): F = (Lambda / (J - 1)) / phi-hat,
#    phi-hat = X^2_Pearson(M1) / (C - J), against F_{J-1, C-J}. This is the
#    classical QL F: glmGamPoi additionally shrinks phi-hat across genes by
#    empirical Bayes (and moves the denominator df), which is not done here.
# 4. Permutation, Freedman-Lane (Winkler et al. 2014): Pearson residuals of
#    M0 (label-free, psi refitted under M0), Welch's G statistic (Winkler
#    et al. 2014, Eq. 6-7; for a one-way design Welch's 1951 v^2), B = 2000
#    label permutations, the same permutation for every gene,
#    p = (1 + #{G* >= G}) / (B + 1). The plain F on the same residuals is
#    permuted alongside for comparison.
#
# Seed: 20261101 (simulation, then the permutations).
###############################################################################

library(ggplot2)
library(grid)
library(DeCovarT) # theme_decovart_facets(), the repository theme

t_start <- proc.time()
set.seed(20261101)

out_dir <- file.path("vignettes", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_series <- "#1d4e89"
col_second <- "#c45c26"
col_guide <- "#1c2b24"

# ---- simulate ---------------------------------------------------------------
n_type <- c(2000L, 300L, 60L, 25L)
n_types <- length(n_type)
n_cells <- sum(n_type)
type <- factor(rep(seq_len(n_types), n_type))
s_c <- exp(stats::rnorm(n_cells, mean = 0, sd = 0.5))
log_s <- log(s_c)

n_null_homo <- 600L
n_null_het <- 300L
n_alt <- 100L
n_genes <- n_null_homo + n_null_het + n_alt
gene_class <- factor(
  rep(c("homo", "het", "alt"), c(n_null_homo, n_null_het, n_alt)),
  levels = c("homo", "het", "alt")
)
beta0 <- stats::rnorm(n_genes, mean = -1, sd = 1.2)

# dispersion per gene and type (genes x types)
psi_g <- exp(stats::runif(n_genes, log(2), log(20)))
psi_gj <- matrix(psi_g, n_genes, n_types)
het <- gene_class == "het"
psi_gj[het, ] <- matrix(c(20, 20, 1, 1), sum(het), n_types, byrow = TRUE)

# log fold changes per gene and type (genes x types)
lfc <- matrix(0, n_genes, n_types)
for (g in which(gene_class == "alt")) {
  hit <- sample.int(n_types, sample.int(2L, 1L))
  lfc[g, hit] <- sample(c(-1, 1), length(hit), replace = TRUE) *
    stats::runif(length(hit), 0.5, 1.5)
}

j_c <- as.integer(type)
m_true <- exp(beta0 + lfc[, j_c] + rep(log_s, each = n_genes))
counts <- matrix(
  stats::rnbinom(length(m_true), mu = m_true, size = psi_gj[, j_c]),
  nrow = n_genes
)
stopifnot(all(rowSums(counts) > 0))

# presence gate of eq-fs-presence (detected in >= 5% of the cells of every
# type and non-zero variance in every type): reported, not applied
det_frac <- sapply(seq_len(n_types), function(j) {
  rowMeans(counts[, j_c == j] > 0)
})
var_j <- sapply(seq_len(n_types), function(j) {
  apply(counts[, j_c == j], 1, stats::var)
})
pass_gate <- apply(det_frac, 1, min) >= 0.05 & apply(var_j, 1, min) > 0
cat("genes failing the presence gate:", sum(!pass_gate), "of", n_genes, "\n")
print(table(gene_class, pass_gate))

# ---- tests 1 and 3: NB GLM per gene (LR and QL F) ---------------------------
X0 <- matrix(1, n_cells, 1)
X1 <- stats::model.matrix(~type)
fit_gene <- function(x) {
  d <- data.frame(x = x, type = type, log_s = log_s)
  f1 <- suppressWarnings(
    MASS::glm.nb(x ~ type + offset(log_s), data = d)
  )
  theta <- f1$theta
  fam <- MASS::negative.binomial(theta)
  # both models refitted at the final theta, so that D0 - D1 compares two
  # nested fits with the same variance function
  g1 <- suppressWarnings(
    stats::glm.fit(X1, x, offset = log_s, family = fam)
  )
  g0 <- suppressWarnings(
    stats::glm.fit(X0, x, offset = log_s, family = fam)
  )
  m1 <- g1$fitted.values
  pearson1 <- sum((x - m1)^2 / (m1 + m1^2 / theta))
  # M0 with its own theta: the label-free fit whose residuals are permuted
  f0 <- suppressWarnings(MASS::glm.nb(x ~ offset(log_s), data = d))
  m0 <- f0$fitted.values
  r0 <- (x - m0) / sqrt(m0 + m0^2 / f0$theta)
  list(
    stat = c(
      theta1 = theta,
      theta0 = f0$theta,
      D0 = g0$deviance,
      D1 = g1$deviance,
      pearson1 = pearson1
    ),
    r0 = r0
  )
}
t_glm <- system.time(
  fits <- lapply(seq_len(n_genes), function(g) fit_gene(counts[g, ]))
)
cat("NB GLM fits:", round(t_glm[["elapsed"]], 1), "s\n")
glm_stat <- do.call(rbind, lapply(fits, `[[`, "stat"))
resid0 <- do.call(rbind, lapply(fits, `[[`, "r0")) # genes x cells
rm(fits)

df1 <- n_types - 1
df2 <- n_cells - n_types
lambda <- pmax(glm_stat[, "D0"] - glm_stat[, "D1"], 0)
phi_hat <- glm_stat[, "pearson1"] / df2
f_ql <- (lambda / df1) / phi_hat
logp_lr <- stats::pchisq(lambda, df1, lower.tail = FALSE, log.p = TRUE)
logp_ql <- stats::pf(f_ql, df1, df2, lower.tail = FALSE, log.p = TRUE)

# ---- test 2: one-way ANOVA F on the shifted logarithm -----------------------
alpha_vst <- 0.05
y0 <- 1 / (4 * alpha_vst)
Z <- stats::model.matrix(~ type - 1) # cells x types indicator
oneway_f <- function(Y, Z) {
  n_j <- colSums(Z)
  Yc <- Y - rowMeans(Y) # centred: grand mean 0, better conditioned sums
  sums <- Yc %*% Z
  ss_b <- rowSums(sweep(sums^2, 2, n_j, "/"))
  ss_w <- rowSums(Yc^2) - ss_b
  (ss_b / (ncol(Z) - 1)) / (ss_w / (nrow(Z) - ncol(Z)))
}
y_vst <- log(sweep(counts, 2, s_c, "/") + y0)
f_log <- oneway_f(y_vst, Z)
logp_log <- stats::pf(f_log, df1, df2, lower.tail = FALSE, log.p = TRUE)

# ---- test 4: permutation of M0 Pearson residuals, Welch's G -----------------
# Freedman-Lane permutes the residuals of the reduced model and adds them back
# to its fit. Here the reduced model is M0 (intercept and offset only), the
# statistic is computed on the residuals themselves, and the fitted part
# carries no type information, so permuting the M0 residuals against fixed
# labels is the same as permuting the labels against fixed residuals. One
# permutation of the C labels is applied to all genes at once.
welch_and_f <- function(S1, S2, n_j, tot, tot2) {
  # S1, S2: genes x J sums of r and r^2 per type; tot, tot2: their row totals.
  # Returns Welch's G, the one-way F on the same values, and Welch's
  # denominator df nu (used only for the analytic diagnostic below).
  J <- length(n_j)
  n <- sum(n_j)
  mean_j <- sweep(S1, 2, n_j, "/")
  var_j <- sweep(S2 - sweep(S1^2, 2, n_j, "/"), 2, n_j - 1, "/")
  w <- sweep(1 / var_j, 2, n_j, "*")
  w_tot <- rowSums(w)
  r_w <- rowSums(w * mean_j) / w_tot
  num <- rowSums(w * (mean_j - r_w)^2) / (J - 1)
  a <- rowSums(sweep((1 - w / w_tot)^2, 2, n_j - 1, "/"))
  den <- 1 + 2 * (J - 2) / (J^2 - 1) * a
  ss_b <- rowSums(sweep(S1^2, 2, n_j, "/")) - tot^2 / n
  ss_w <- tot2 - rowSums(sweep(S1^2, 2, n_j, "/"))
  cbind(
    G = num / den,
    F = (ss_b / (J - 1)) / (ss_w / (n - J)),
    nu = (J^2 - 1) / (3 * a)
  )
}
resid0_sq <- resid0^2
tot <- rowSums(resid0)
tot2 <- rowSums(resid0_sq)
obs <- welch_and_f(resid0 %*% Z, resid0_sq %*% Z, n_type, tot, tot2)
# check the hand-coded statistics against stats::oneway.test on three genes
for (g in c(1L, 650L, 950L)) {
  d_g <- data.frame(r = resid0[g, ], type = type)
  stopifnot(
    isTRUE(all.equal(
      unname(obs[g, "G"]),
      unname(stats::oneway.test(r ~ type, d_g, var.equal = FALSE)$statistic)
    )),
    isTRUE(all.equal(
      unname(obs[g, "F"]),
      unname(stats::oneway.test(r ~ type, d_g, var.equal = TRUE)$statistic)
    ))
  )
}

n_perm <- 2000L
perm_lab <- replicate(n_perm, sample(j_c)) # cells x B, shared by all genes
exceed <- matrix(0L, n_genes, 2, dimnames = list(NULL, c("G", "F")))
chunk <- 250L
t_perm <- system.time({
  for (b0 in seq(1L, n_perm, by = chunk)) {
    bb <- b0:min(b0 + chunk - 1L, n_perm)
    # stacked indicator matrices: column (b - 1) * J + j flags type j in
    # permutation b, so one matrix product gives the type sums of all genes
    # for a block of permutations
    Zp <- matrix(0, n_cells, length(bb) * n_types)
    idx <- cbind(
      rep(seq_len(n_cells), length(bb)),
      as.vector(sweep(
        perm_lab[, bb, drop = FALSE],
        2,
        (seq_along(bb) - 1L) * n_types,
        "+"
      ))
    )
    Zp[idx] <- 1
    S1 <- resid0 %*% Zp
    S2 <- resid0_sq %*% Zp
    for (k in seq_along(bb)) {
      cols <- (k - 1L) * n_types + seq_len(n_types)
      st <- welch_and_f(S1[, cols], S2[, cols], n_type, tot, tot2)
      exceed <- exceed + (st[, c("G", "F")] >= obs[, c("G", "F")])
    }
  }
})
cat("permutations:", round(t_perm[["elapsed"]], 1), "s\n")
p_perm <- (1 + exceed) / (n_perm + 1)

# ---- p-values, BH, error rates ----------------------------------------------
# permF, the plain F permuted on the same residuals, separates the statistic
# from the reference distribution; it is drawn in panels a and c only.
tests <- c("LR", "QLF", "logF", "Welch", "permF")
test_lab <- c(
  LR = "LR (deviance)",
  QLF = "QL F",
  logF = "ANOVA F, shifted log",
  Welch = "Welch G, permuted",
  permF = "F, permuted"
)
logp <- cbind(
  LR = logp_lr,
  QLF = logp_ql,
  logF = logp_log,
  Welch = log(p_perm[, "G"]),
  permF = log(p_perm[, "F"])
)
pval <- exp(logp)
# diagnostic, not drawn: Welch's G against its analytic F_{J-1, nu} tail
p_welch_an <- stats::pf(obs[, "G"], df1, obs[, "nu"], lower.tail = FALSE)
p_all <- cbind(pval, Welch_analytic = p_welch_an)

type1 <- sapply(c("homo", "het"), function(cl) {
  colMeans(p_all[gene_class == cl, ] < 0.05)
})
type1_gate <- sapply(c("homo", "het"), function(cl) {
  colMeans(p_all[gene_class == cl & pass_gate, ] < 0.05)
})
colnames(type1_gate) <- paste0(colnames(type1_gate), "_gated")
cat("\ntype I error at 0.05 (rows: tests; columns: null class; gated =\n")
cat("only the genes that pass the presence gate)\n")
print(round(cbind(type1, type1_gate), 3))
cat(
  "Monte Carlo SE at 0.05:",
  round(sqrt(0.05 * 0.95 / c(n_null_homo, n_null_het)), 4),
  "\n"
)

bh <- apply(p_all, 2, stats::p.adjust, method = "BH") < 0.05
bh_tab <- data.frame(
  test = colnames(bh),
  discoveries = colSums(bh),
  false_homo = colSums(bh & gene_class == "homo"),
  false_het = colSums(bh & gene_class == "het"),
  true = colSums(bh & gene_class == "alt")
)
bh_tab$fdr <- with(
  bh_tab,
  ifelse(discoveries > 0, (false_homo + false_het) / discoveries, 0)
)
bh_tab$power <- bh_tab$true / n_alt
bh_tab$power_gated <- colSums(bh & gene_class == "alt" & pass_gate) /
  sum(gene_class == "alt" & pass_gate)
cat("\nBH at 0.05 (power_gated: alternatives that pass the presence gate)\n")
print(bh_tab, row.names = FALSE, digits = 3)
cat(
  "unadjusted power at 0.05:",
  round(colMeans(p_all[gene_class == "alt", ] < 0.05), 2),
  "\n"
)

# ---- Kendall's tau between rankings -----------------------------------------
# tau-b between the p-value rankings (log scale for the analytic tests, so no
# tie at underflow), for the four tests of the comparison. The permutation
# p-value has a floor of 1 / (B + 1), and the ties there are kept; breaking
# them by the observed G is printed as a check.
tau_tests <- c("LR", "QLF", "logF", "Welch")
tau_mat <- function(keys) {
  tt <- stats::cor(keys[, tau_tests], method = "kendall")
  dimnames(tt) <- list(tau_tests, tau_tests)
  tt
}
tau_all <- tau_mat(logp)
tau_alt <- tau_mat(logp[gene_class == "alt", ])
rank_key <- logp
rank_key[, "Welch"] <- order(order(p_perm[, "G"], -obs[, "G"]))
cat("\nKendall tau-b, all genes\n")
print(round(tau_all, 2))
cat("Kendall tau-b, alternatives\n")
print(round(tau_alt, 2))
cat("check, Welch ties broken by G: all genes / alternatives\n")
print(round(tau_mat(rank_key)["Welch", ], 2))
print(round(tau_mat(rank_key[gene_class == "alt", ])["Welch", ], 2))
cat(
  "alternatives at the permutation floor 1/(B+1):",
  sum(p_perm[gene_class == "alt", "G"] == 1 / (n_perm + 1)),
  "\n"
)
cat(
  "median psi-hat (M1), heteroscedastic null:",
  round(stats::median(glm_stat[het, "theta1"]), 2),
  "; median phi-hat, heteroscedastic / homoscedastic null:",
  round(stats::median(phi_hat[het]), 3),
  "/",
  round(stats::median(phi_hat[gene_class == "homo"]), 3),
  "\n"
)

# ---- figure -----------------------------------------------------------------
theme_fig <- theme_decovart_facets(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    plot.title = element_text(face = "bold", size = 12, colour = col_guide),
    plot.title.position = "plot",
    plot.margin = margin(6, 10, 6, 6)
  )
# blue: deviance-based (LR, QL F); orange: shifted-log ANOVA; near-black:
# permutation of the M0 residuals (Welch G solid, plain F dashed)
test_col <- c(col_series, col_series, col_second, col_guide, col_guide)
test_lty <- c("solid", "42", "solid", "solid", "42")
names(test_col) <- names(test_lty) <- test_lab

# (a) QQ plot of null p-values, with a pointwise 95% band for uniform order
# statistics (Beta(i, n - i + 1))
class_lab <- c(
  homo = sprintf("homoscedastic null (%d genes)", n_null_homo),
  het = sprintf("heteroscedastic null (%d genes)", n_null_het)
)
qq <- do.call(
  rbind,
  lapply(c("homo", "het"), function(cl) {
    sel <- gene_class == cl
    n <- sum(sel)
    do.call(
      rbind,
      lapply(tests, function(tt) {
        data.frame(
          class = class_lab[[cl]],
          test = test_lab[[tt]],
          expected = -log10(stats::ppoints(n)),
          observed = sort(-logp[sel, tt] / log(10), decreasing = TRUE)
        )
      })
    )
  })
)
band <- do.call(
  rbind,
  lapply(c("homo", "het"), function(cl) {
    n <- sum(gene_class == cl)
    i <- seq_len(n)
    data.frame(
      class = class_lab[[cl]],
      expected = -log10(stats::ppoints(n)),
      lo = -log10(stats::qbeta(0.975, i, n - i + 1)),
      hi = -log10(stats::qbeta(0.025, i, n - i + 1))
    )
  })
)
qq$class <- factor(qq$class, levels = class_lab)
band$class <- factor(band$class, levels = class_lab)
qq$test <- factor(qq$test, levels = test_lab)
y_cap <- min(8, max(qq$observed)) # values above 8 would be clipped
pa <- ggplot(qq, aes(expected, pmin(observed, y_cap))) +
  geom_ribbon(
    data = band,
    aes(x = expected, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "grey88"
  ) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.4) +
  geom_line(aes(colour = test, linetype = test), linewidth = 0.8) +
  facet_wrap(~class, nrow = 1) +
  scale_colour_manual(values = test_col) +
  scale_linetype_manual(values = test_lty) +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(
    title = "a   Null p-values against the uniform",
    x = expression("expected " * -log[10] * italic(p)),
    y = expression("observed " * -log[10] * italic(p))
  ) +
  theme_fig +
  guides(
    colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1)),
    linetype = guide_legend(nrow = 1)
  ) +
  theme(legend.key.width = unit(2, "lines"))

# (b) Kendall's tau between rankings
short_lab <- c(LR = "LR", QLF = "QL F", logF = "ANOVA", Welch = "Welch")
n_tau <- length(tau_tests)
tau_df <- do.call(
  rbind,
  lapply(
    list(
      list(tau_all, sprintf("all %d genes", n_genes)),
      list(tau_alt, sprintf("%d alternatives", n_alt))
    ),
    function(z) {
      data.frame(
        set = z[[2]],
        row = factor(
          short_lab[rep(tau_tests, n_tau)],
          levels = rev(short_lab)
        ),
        col = factor(
          short_lab[rep(tau_tests, each = n_tau)],
          levels = short_lab
        ),
        tau = as.vector(z[[1]])
      )
    }
  )
)
tau_df$set <- factor(tau_df$set, levels = unique(tau_df$set))
pb <- ggplot(tau_df, aes(col, row, fill = tau)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.2f", tau), colour = tau > 0.6),
    size = 3.6
  ) +
  facet_wrap(~set, nrow = 1) +
  scale_fill_gradient(
    low = "white",
    high = col_series,
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(`TRUE` = "white", `FALSE` = col_guide),
    guide = "none"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = expression(
      bold("b   Kendall's ") *
        bolditalic(tau) *
        bold(" between p-value rankings")
    ),
    x = NULL,
    y = NULL
  ) +
  theme_fig +
  theme(panel.grid = element_blank())

# (c) BH at 0.05: discoveries, empirical FDR, power; axes start at zero
metric_lab <- c(
  discoveries = "BH discoveries",
  fdr = "empirical FDR",
  power = "power"
)
bh_main <- bh_tab[match(tests, bh_tab$test), ]
dots <- do.call(
  rbind,
  lapply(names(metric_lab), function(mm) {
    data.frame(
      test = factor(test_lab[bh_main$test], levels = rev(test_lab)),
      metric = factor(metric_lab[[mm]], levels = metric_lab),
      value = bh_main[[mm]]
    )
  })
)
refs <- data.frame(
  metric = factor(metric_lab[c("discoveries", "fdr")], levels = metric_lab),
  value = c(n_alt, 0.05)
)
x_range <- data.frame(
  metric = factor(rep(metric_lab, each = 2), levels = metric_lab),
  value = c(
    0,
    1.1 * max(bh_main$discoveries),
    0,
    max(0.1, 1.1 * max(bh_main$fdr)),
    0,
    1
  ),
  test = test_lab[["LR"]]
)
pc <- ggplot(dots, aes(value, test)) +
  geom_blank(data = x_range) +
  geom_vline(
    data = refs,
    aes(xintercept = value),
    colour = col_guide,
    linewidth = 0.4,
    linetype = "22"
  ) +
  geom_point(colour = col_series, size = 2.6) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_y_discrete(limits = unname(rev(test_lab))) +
  scale_x_continuous(
    breaks = function(lims) pretty(lims, n = 3),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(title = "c   Benjamini-Hochberg at 0.05", x = NULL, y = NULL) +
  theme_fig +
  theme(panel.grid.minor = element_blank())

# ---- assemble with grid.layout and write --------------------------------------
draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(
    layout = grid.layout(
      nrow = 2,
      ncol = 2,
      heights = unit(c(1, 0.85), "null"),
      widths = unit(c(1, 1.05), "null")
    )
  ))
  print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(pb, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pc, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  popViewport()
}
fig_w <- 10
fig_h <- 7.4
svglite::svglite(
  file.path(out_dir, "fig_03_1_type_tests.svg"),
  width = fig_w,
  height = fig_h,
  bg = "white"
)
draw_figure()
invisible(grDevices::dev.off())
ragg::agg_png(
  file.path(out_dir, "fig_03_1_type_tests.png"),
  width = fig_w,
  height = fig_h,
  units = "in",
  res = 200,
  background = "white"
)
draw_figure()
invisible(grDevices::dev.off())

t_total <- proc.time() - t_start
cat("\ntotal run time:", round(t_total[["elapsed"]] / 60, 2), "min\n")
utils::sessionInfo()
