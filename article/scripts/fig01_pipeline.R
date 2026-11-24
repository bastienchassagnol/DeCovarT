#!/usr/bin/env Rscript
# Figure 1 of article/main_layout.tex: the DeCovarT pipeline applied to the
# Suppinger gastruloids. This script draws the ggplot2 panels only; the
# layout, the schematic shapes and the graphical model are drawn by
# article/scripts/fig01_pipeline.tex, which includes these panels.
#
#   a1  presence gate          symbolic (seeded toy counts)
#   a2  permuted Welch G       symbolic (seeded toy null)
#   a3  Benjamini-Hochberg     symbolic (seeded toy p-values)
#   a4  gene funnel            real counts, tier (i) and tier (ii) tables
#   b1  shrinkage of R_j       symbolic (seeded toy correlation)
#   b2  bulk log-likelihood    symbolic curve, no values on the y axis
#   d   96 h compositions      real DeCovarT estimates, standard and early Wnt
#
# Inputs, read from the reproducibility repository (sibling folder):
#   output/feature_selection/tier_01/gene_filter_summary.csv
#   output/feature_selection/tier_02/selected_<time>.csv
#   output/deconvolution/proportions_long.csv
#   output/deconvolution/observed_proportions.csv
#   data/dictionaries/suppinger_celltype_colours.csv
#
# Build, from the DeCovarT repository root:
#   Rscript --no-save --no-restore article/scripts/fig01_pipeline.R
#   cd article/scripts && latexmk -lualatex -outdir=../figs fig01_pipeline.tex
#   latexmk -c -outdir=../figs fig01_pipeline.tex
#   cd ../figs && pdftocairo -svg fig01_pipeline.pdf fig01_pipeline.svg
#
# Outputs: article/figs/fig01_pipeline/panel_*.{svg,pdf} (svglite, then
# Inkscape for the PDF that LuaLaTeX includes) and
# article/figs/fig01_pipeline/values.tex (gene counts and cell-type colours
# for TikZ).
# Each panel is written at the size of its slot in the .tex, in mm, so
# nothing is rescaled and the 7 pt type stays 7 pt on the A4 page.

library(ggplot2)

repro <- Sys.getenv(
  "DECOVART_REPRO",
  file.path("..", "DeCovarT_reproducibility")
)
if (!dir.exists(repro)) {
  stop(
    "Reproducibility repository not found at ", repro,
    "; set DECOVART_REPRO to its path."
  )
}
read_repro <- function(...) {
  utils::read.csv(
    file.path(repro, ...),
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8",
    check.names = FALSE
  )
}

out_dir <- file.path("article", "figs", "fig01_pipeline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- House style -----------------------------------------------------------

colour_series <- "#1d4e89"
colour_mark <- "#c45c26"
colour_ink <- "#1c2b24"
colour_grid <- "#e6e1d8"
colour_fail <- "#b9b4aa"
font_family <- "Noto Sans"
pt <- 7

theme_slide <- function(base_size = pt) {
  ggplot2::theme_minimal(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = colour_grid, linewidth = 0.25),
      axis.title = ggplot2::element_text(colour = colour_ink, size = base_size),
      axis.text = ggplot2::element_text(colour = colour_ink, size = base_size),
      strip.text = ggplot2::element_text(
        face = "bold", colour = "#0e4d3a", size = base_size, margin = ggplot2::margin(0, 0, 1, 0)
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(1, 1, 1, 1),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA)
    )
}
mm_text <- pt / ggplot2::.pt # geom_text size giving 7 pt

save_panel <- function(plot, name, width_mm, height_mm) {
  svg <- file.path(out_dir, paste0("panel_", name, ".svg"))
  pdf <- sub("\\.svg$", ".pdf", svg)
  ggplot2::ggsave(
    svg, plot,
    device = svglite::svglite,
    width = width_mm, height = height_mm, units = "mm",
    bg = "transparent"
  )
  inkscape <- Sys.which("inkscape")
  if (!nzchar(inkscape)) stop("Inkscape is needed to turn ", svg, " into PDF.")
  status <- system2(
    inkscape,
    c(shQuote(svg), "--export-type=pdf", paste0("--export-filename=", shQuote(pdf))),
    stdout = FALSE, stderr = FALSE
  )
  if (status != 0L || !file.exists(pdf)) stop("Inkscape failed on ", svg)
  invisible(pdf)
}

times <- c("48h", "72h", "96h")
time_colours <- c("48h" = "#9db8d9", "72h" = "#5a86b8", "96h" = colour_series)

# --- a1: presence gate (symbolic) -----------------------------------------
# Eight genes in four cell types. Gene 3 is a perfect marker (zero outside
# one type) and gene 6 is null everywhere: both fail the gate.

a1_data <- withr::with_seed(20261123L, {
  grid <- expand.grid(type = 1:4, gene = 1:8)
  grid$value <- stats::rgamma(nrow(grid), shape = 2, rate = 1.2)
  grid$value[grid$gene == 3 & grid$type != 2] <- 0
  grid$value[grid$gene == 3 & grid$type == 2] <- 4
  grid$value[grid$gene == 6] <- 0
  grid$fails <- grid$gene %in% c(3, 6)
  grid
})
a1_crosses <- data.frame(gene = c(3, 6), label = c("marker", "null"))

panel_a1 <- ggplot(a1_data, aes(gene, type)) +
  geom_tile(aes(fill = value), colour = "white", linewidth = 0.4) +
  geom_tile(
    data = a1_data[a1_data$fails, ], fill = NA,
    colour = colour_mark, linewidth = 0.25, linetype = "22"
  ) +
  annotate(
    "segment",
    x = a1_crosses$gene - 0.45, xend = a1_crosses$gene + 0.45,
    y = 0.5, yend = 4.5, colour = colour_mark, linewidth = 0.5
  ) +
  geom_text(
    data = a1_crosses, aes(x = gene, y = 4.95, label = label),
    inherit.aes = FALSE, size = mm_text, colour = colour_mark,
    family = font_family, vjust = 0
  ) +
  scale_fill_gradient(low = "#f4f1ea", high = colour_series) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = c(0, 0.9))) +
  labs(x = "genes", y = "cell types") +
  theme_slide() +
  theme(axis.text = element_blank(), panel.grid.major = element_blank())

# --- a2: permuted Welch G (symbolic) --------------------------------------

a2_null <- withr::with_seed(20261124L, data.frame(g = stats::rchisq(1999L, df = 4)))
g_observed <- 13.5
panel_a2 <- ggplot(a2_null, aes(g)) +
  geom_histogram(
    aes(fill = g >= g_observed),
    binwidth = 0.75, boundary = 0, colour = "white", linewidth = 0.15
  ) +
  geom_vline(xintercept = g_observed, colour = colour_mark, linewidth = 0.6) +
  annotate(
    "text", x = g_observed + 0.4, y = Inf, vjust = 1.3, hjust = 0,
    label = "observed", size = mm_text, colour = colour_mark, family = font_family
  ) +
  scale_fill_manual(values = c(`FALSE` = "#d9d3c7", `TRUE` = colour_mark)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Welch G, permuted labels", y = NULL) +
  theme_slide() +
  theme(axis.text = element_blank(), panel.grid.major = element_blank())

# --- a3: Benjamini-Hochberg (symbolic) ------------------------------------

a3_data <- withr::with_seed(20261125L, {
  p <- sort(c(stats::rbeta(12L, 0.25, 30), stats::runif(18L)))
  m <- length(p)
  q <- 0.05
  k <- seq_len(m)
  passes <- k <= max(c(0L, which(p <= q * k / m)))
  data.frame(k = k, p = p, passes = passes, m = m, q = q)
})
a3_zoom <- a3_data[a3_data$k <= 16L, ]
panel_a3 <- ggplot(a3_zoom, aes(k, p)) +
  geom_abline(
    intercept = 0, slope = a3_data$q[1] / a3_data$m[1],
    colour = colour_ink, linetype = "dashed", linewidth = 0.4
  ) +
  geom_point(aes(colour = passes), size = 1.1) +
  annotate(
    "text", x = 16, y = 0.05 * 16 / a3_data$m[1], hjust = 1, vjust = -0.6,
    label = "0.05 k / m", size = mm_text, colour = colour_ink, family = font_family
  ) +
  scale_colour_manual(values = c(`FALSE` = colour_fail, `TRUE` = colour_mark)) +
  coord_cartesian(ylim = c(0, 0.06)) +
  labs(x = "rank k of the p-value", y = NULL) +
  theme_slide() +
  theme(axis.text = element_blank())

# --- a4: gene funnel (real counts) -----------------------------------------

filters <- read_repro("output", "feature_selection", "tier_01", "gene_filter_summary.csv")
n_deposited <- unique(filters$genes[filters$rule == "deposited genes"])
stopifnot(length(n_deposited) == 1L)
n_presence <- vapply(times, function(tm) {
  filters$genes[filters$time_point == tm & filters$rule == "enters tier (ii)"]
}, numeric(1))
selected <- lapply(times, function(tm) {
  read_repro("output", "feature_selection", "tier_02", paste0("selected_", tm, ".csv"))$gene
})
names(selected) <- times
n_welch <- vapply(selected, length, integer(1))
n_union <- length(unique(unlist(selected)))

stage_levels <- c(
  "deposited symbols",
  "(1) present in every type\n48 / 72 / 96 h",
  "(2)+(3) Welch G, BH 0.05\n48 / 72 / 96 h",
  "panel, union over times"
)
funnel <- rbind(
  data.frame(stage = stage_levels[1], time = "all", genes = n_deposited, slot = 0),
  data.frame(stage = stage_levels[2], time = times, genes = n_presence, slot = c(-1, 0, 1)),
  data.frame(stage = stage_levels[3], time = times, genes = n_welch, slot = c(-1, 0, 1)),
  data.frame(stage = stage_levels[4], time = "union", genes = n_union, slot = 0)
)
funnel$row <- match(funnel$stage, rev(stage_levels))
funnel$y <- funnel$row + 0.27 * -funnel$slot
funnel$height <- ifelse(funnel$time %in% times, 0.24, 0.55)
funnel$fill <- ifelse(
  funnel$time %in% times, time_colours[funnel$time],
  ifelse(funnel$time == "union", colour_mark, "#8a8f8b")
)
# One label per stage; the three dodged bars are 48, 72 and 96 h, top down.
stage_label <- function(rows) {
  if (all(rows$time %in% times)) {
    paste(format(rows$genes, big.mark = ",", trim = TRUE), collapse = " / ")
  } else {
    format(rows$genes, big.mark = ",", trim = TRUE)
  }
}
funnel_labels <- do.call(rbind, lapply(split(funnel, funnel$row), function(rows) {
  data.frame(row = rows$row[1], x = max(rows$genes) / 2 + 400, label = stage_label(rows))
}))
half_width <- n_deposited / 2

panel_a4 <- ggplot(funnel) +
  geom_rect(
    aes(
      xmin = -genes / 2, xmax = genes / 2,
      ymin = y - height / 2, ymax = y + height / 2, fill = fill
    ),
    colour = NA
  ) +
  geom_text(
    data = funnel_labels, aes(x = x, y = row, label = label),
    hjust = 0, size = mm_text, colour = colour_ink, family = font_family
  ) +
  scale_fill_identity() +
  scale_y_continuous(
    breaks = seq_along(stage_levels), labels = rev(stage_levels),
    expand = expansion(add = 0.3)
  ) +
  scale_x_continuous(limits = c(-half_width, half_width + 8000), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_slide() +
  theme(axis.text.x = element_blank(), panel.grid.major = element_blank())

# --- b1: shrinkage of a correlation matrix (symbolic) ----------------------

b1_r <- withr::with_seed(20261126L, {
  g <- 8L
  loadings <- matrix(stats::rnorm(g * 2L), g, 2L)
  s <- tcrossprod(loadings) + diag(stats::runif(g, 0.3, 0.8))
  stats::cov2cor(s)
})
lambdas <- c(0, 0.5, 1)
b1_data <- do.call(rbind, lapply(lambdas, function(lam) {
  r <- (1 - lam) * b1_r + lam * diag(nrow(b1_r))
  d <- expand.grid(row = seq_len(nrow(r)), col = seq_len(ncol(r)))
  d$value <- r[cbind(d$row, d$col)]
  d$lambda <- lam
  d
}))
b1_data$lambda <- factor(
  b1_data$lambda, levels = lambdas,
  labels = paste0("lambda == ", lambdas)
)
panel_b1 <- ggplot(b1_data, aes(col, -row, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~lambda, nrow = 1L, labeller = label_parsed) +
  scale_fill_gradient2(
    low = colour_series, mid = "#f7f5f0", high = colour_mark,
    midpoint = 0, limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_slide() +
  theme(axis.text = element_blank(), panel.grid.major = element_blank(), panel.spacing = unit(2, "mm"))

# --- b2: profile log-likelihood of the bulk against lambda (symbolic) ------
# Shape only: a concave curve with its maximum inside (0, 1). The grid of
# points mimics the swept levels; the y axis carries no values.

b2_grid <- c(0, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 1)
b2_fun <- function(lam) -(1.6 * (lam - 0.5)^2 + 2.2 * exp(-18 * lam))
b2_curve <- data.frame(lambda = seq(0, 1, length.out = 201L))
b2_curve$ll <- b2_fun(b2_curve$lambda)
b2_points <- data.frame(lambda = b2_grid, ll = b2_fun(b2_grid))
lambda_hat <- b2_points$lambda[which.max(b2_points$ll)]

panel_b2 <- ggplot(b2_curve, aes(lambda, ll)) +
  geom_vline(xintercept = lambda_hat, colour = colour_mark, linetype = "dashed", linewidth = 0.4) +
  geom_line(colour = colour_series, linewidth = 0.6) +
  geom_point(data = b2_points, colour = colour_series, size = 1) +
  geom_point(
    data = b2_points[b2_points$lambda == lambda_hat, ],
    colour = colour_mark, size = 1.8
  ) +
  annotate(
    "text", x = lambda_hat + 0.03, y = min(b2_curve$ll), hjust = 0, vjust = 0,
    label = "hat(lambda)", parse = TRUE, size = mm_text * 1.15,
    colour = colour_mark, family = font_family
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(add = 0.02)) +
  labs(
    x = expression(paste("shrinkage level ", lambda)),
    y = "profile log-likelihood\nof the bulk replicates"
  ) +
  theme_slide() +
  theme(axis.text.y = element_blank())

# --- d: compositions at 96 h (real estimates) ------------------------------

colours_tab <- read_repro("data", "dictionaries", "suppinger_celltype_colours.csv")
type_colours <- stats::setNames(colours_tab$colour, colours_tab$celltypeannotation)

method_d <- "DeCovarT shrinkage, bulk likelihood | welch"
props <- read_repro("output", "deconvolution", "proportions_long.csv")
props <- props[props$method == method_d & props$time == "96h", ]
if (nrow(props) == 0L) stop("No 96 h rows for ", method_d)
observed <- read_repro("output", "deconvolution", "observed_proportions.csv")
observed <- observed[observed$time == "96h", ]

types_96 <- colours_tab$celltypeannotation[colours_tab$celltypeannotation %in% props$cell_type]
arm_labels <- c(
  standard = "standard Wnt (pulse at 48 h)",
  early = "early Wnt (pulse at 24 h)"
)
props$arm <- factor(arm_labels[props$arm], levels = arm_labels)
props$cell_type <- factor(props$cell_type, levels = rev(types_96))
means <- stats::aggregate(proportion ~ arm + cell_type, data = props, FUN = mean)
observed$arm <- factor(arm_labels[["standard"]], levels = arm_labels)
observed$cell_type <- factor(observed$cell_type, levels = rev(types_96))

panel_d <- ggplot(props, aes(proportion, cell_type, colour = cell_type)) +
  geom_errorbar(
    data = means, aes(xmin = proportion, xmax = proportion),
    width = 0.7, linewidth = 0.9, orientation = "y"
  ) +
  geom_point(size = 1.4, alpha = 0.85) +
  geom_point(
    data = observed, aes(x = observed), shape = 5, size = 1.8,
    colour = colour_ink, stroke = 0.4
  ) +
  facet_wrap(~arm, nrow = 1L) +
  scale_colour_manual(values = type_colours) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  labs(x = "estimated proportion (dots: 3 bulk replicates; bar: mean)", y = NULL) +
  theme_slide() +
  theme(panel.spacing = unit(5, "mm"))

# --- Values and cell-type colours for TikZ ---------------------------------

tikz_names <- c(
  "Ectopic pluripotency" = "ctEctopic",
  "Naïve pluripotency" = "ctNaive",
  "Epiblast" = "ctEpiblast",
  "Primitive streak" = "ctStreak",
  "Caudal mesoderm" = "ctCaudal",
  "Pre-somitic mesoderm" = "ctPSM",
  "Somite differentiation front" = "ctSomite",
  "Gut" = "ctGut"
)
stopifnot(all(names(tikz_names) %in% names(type_colours)))
big <- function(x) format(x, big.mark = ",", trim = TRUE)
writeLines(
  c(
    "% Written by article/scripts/fig01_pipeline.R; do not edit by hand.",
    "% Colours: data/dictionaries/suppinger_celltype_colours.csv (reproducibility repository).",
    sprintf(
      "\\definecolor{%s}{HTML}{%s}",
      tikz_names, sub("^#", "", type_colours[names(tikz_names)])
    ),
    sprintf("\\newcommand{\\nDeposited}{%s}", big(n_deposited)),
    sprintf("\\newcommand{\\nPresence}{%s}", paste(big(n_presence), collapse = " / ")),
    sprintf("\\newcommand{\\nWelch}{%s}", paste(big(n_welch), collapse = " / ")),
    sprintf("\\newcommand{\\nUnion}{%s}", big(n_union)),
    sprintf("\\newcommand{\\nGenesD}{%s}", big(unique(props$genes)))
  ),
  file.path(out_dir, "values.tex"),
  useBytes = TRUE
)

# --- Write the panels at their slot sizes (mm) ------------------------------

save_panel(panel_a1, "a1_presence", 34, 23)
save_panel(panel_a2, "a2_welch", 34, 23)
save_panel(panel_a3, "a3_bh", 34, 23)
save_panel(panel_a4, "a4_funnel", 112, 28)
save_panel(panel_b1, "b1_shrinkage", 60, 24)
save_panel(panel_b2, "b2_loglik", 50, 34)
save_panel(panel_d, "d_composition", 112, 30)

cat(
  "Genes: deposited", n_deposited,
  "| presence", paste(n_presence, collapse = "/"),
  "| Welch G + BH", paste(n_welch, collapse = "/"),
  "| union", n_union, "\n"
)
print(means, digits = 2)
sessioninfo::session_info(pkgs = "attached")
