#!/usr/bin/env Rscript
# Figure 4 of article/main_layout.tex: DeCovarT on the Suppinger gastruloids.
# This script draws the ggplot2 panels c-j and arranges them with
# cowplot::plot_grid() into one block; article/scripts/fig04_suppinger.tex
# draws the schematic panels a (pipeline) and b (design) in TikZ and
# includes the block below them, on one portrait A4 page.
#
#   c  genes removed per tier (i) rule      fs_tier1_gene_filters.rds
#   d  gate and rank (permuted Welch G)     fs_tier2_gate_rank.rds
#   e  UpSet of the confirmed genes (A)     fs_deg_superiority_upset.rds
#   f  curated enrichment themes            fs_enrichment_lollipop.rds
#   g  pooled UMAP, one facet per time      fs_tier1_umap_pooled.rds
#   h  Kendall tau of the mean profiles     profile_correlation.csv (Welch)
#   i  benchmark on the ladder of models    benchmark_metrics.csv (full run)
#   j  early against standard Wnt (CLR)     proportions_long.csv (final model)
#
# Every panel is drawn from the outputs of the reproducibility repository
# (sibling folder, or DECOVART_REPRO); nothing is typed by hand.
#
# Build, from the DeCovarT repository root:
#   Rscript --no-save --no-restore article/scripts/fig04_suppinger.R
#   cd article/scripts && latexmk -lualatex -outdir=../figs fig04_suppinger.tex
#   latexmk -c -outdir=../figs fig04_suppinger.tex
#   cd ../figs && pdftocairo -svg fig04_suppinger.pdf fig04_suppinger.svg
#
# Outputs: article/figs/fig04_suppinger/panel_block.{svg,pdf} (svglite, then
# Inkscape for the PDF that LuaLaTeX includes) and da_contrasts.csv (the
# numbers of panel j). The block is written at the size of its slot in the
# .tex, in mm, so the 7 pt type stays 7 pt on the page.

library(ggplot2)

# cowplot measures the grobs on a null device; the default pdf() one has no
# Noto Sans, so use ragg, which finds fonts through systemfonts like svglite.
cowplot::set_null_device("agg")

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
rds_repro <- function(...) readRDS(file.path(repro, ...))

out_dir <- file.path("article", "figs", "fig04_suppinger")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- House style (as article/scripts/fig01_pipeline.R) ---------------------

colour_series <- "#1d4e89"
colour_mark <- "#c45c26"
colour_ink <- "#1c2b24"
colour_grid <- "#e6e1d8"
colour_context <- "#9a9286"
font_family <- "Noto Sans"
pt <- 7
small <- 6

theme_slide <- function(base_size = pt) {
  ggplot2::theme_minimal(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = colour_grid, linewidth = 0.25),
      axis.title = ggplot2::element_text(colour = colour_ink, size = base_size),
      axis.text = ggplot2::element_text(colour = colour_ink, size = small),
      strip.text = ggplot2::element_text(
        face = "bold", colour = "#0e4d3a", size = base_size, margin = ggplot2::margin(0, 0, 1, 0)
      ),
      legend.position = "none",
      legend.text = ggplot2::element_text(size = small),
      legend.key.size = ggplot2::unit(2.5, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.spacing = ggplot2::unit(1.5, "mm"),
      plot.margin = ggplot2::margin(4, 2, 1, 2),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA)
    )
}
mm_text <- function(size_pt) size_pt / ggplot2::.pt # geom_text size in pt

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
pretty_time <- function(x) sub("h$", " h", x)

# Scored cell types, numbered 1-9 in the dictionary (developmental) order;
# the numbers are the key shared by panels g and h.
dict <- read_repro("data", "dictionaries", "suppinger_celltype_colours.csv")
scored <- c(
  "Naïve pluripotency", "Exiting naïve pluripotency", "Ectopic pluripotency",
  "Epiblast", "Primitive streak", "Caudal mesoderm", "Pre-somitic mesoderm",
  "Somite differentiation front", "Gut"
)
scored <- dict$celltypeannotation[dict$celltypeannotation %in% scored]
stopifnot(length(scored) == 9L)
type_colours <- stats::setNames(dict$colour, dict$celltypeannotation)[scored]
type_number <- stats::setNames(seq_along(scored), scored)

# --- c: genes removed per tier (i) rule -------------------------------------

filters <- rds_repro(
  "output", "feature_selection", "tier_01", "figures", "fs_tier1_gene_filters.rds"
)$gene_summary
rules <- c(
  "zero in every cell of the slot" = "Zero in the slot",
  "not shared one-to-one with the bulk" = "Not one-to-one\nwith the bulk",
  "fails the presence gate (>= 20 cells, >= 2 in every type, variance > 0)" =
    "Fails the\npresence gate"
)
c_data <- filters[filters$rule %in% names(rules), , drop = FALSE]
stopifnot(nrow(c_data) == 9L)
c_data$label <- unname(rules[c_data$rule])
c_data$label <- factor(
  c_data$label,
  levels = names(sort(tapply(c_data$genes, c_data$label, sum)))
)
c_data$time <- pretty_time(c_data$time_point)
c_data$inside <- c_data$genes > 8000
n_deposited <- unique(filters$genes[filters$rule == "deposited genes"])
stopifnot(length(n_deposited) == 1L)

panel_c <- ggplot(c_data, aes(genes, label)) +
  geom_col(fill = colour_series, width = 0.7) +
  geom_text(
    aes(label = scales::comma(genes), hjust = ifelse(inside, 1.08, -0.12),
        colour = inside),
    size = mm_text(small), family = font_family
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = colour_ink)) +
  scale_x_continuous(
    breaks = c(0, 10000, 20000),
    labels = c("0", "10k", "20k"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  facet_wrap(~time, nrow = 1) +
  labs(x = paste0("Genes removed (of ", scales::comma(n_deposited), ")"), y = NULL) +
  theme_slide() +
  theme(panel.grid.major.y = element_blank())

# --- d: gate and rank ---------------------------------------------------------

n_perm <- 2000L # B of the permuted Welch G (docs/02-feature-selection.qmd)
gate_colours <- c(
  "Fails BH" = colour_context,
  "Passes BH 0.05" = colour_mark,
  "At the floor 1/(B+1)" = colour_series
)
ranked <- rds_repro(
  "output", "feature_selection", "tier_02", "figures", "fs_tier2_gate_rank.rds"
)
stopifnot(identical(names(ranked), times))
d_data <- do.call(rbind, lapply(times, function(tm) {
  r <- ranked[[tm]]
  data.frame(
    time = pretty_time(tm),
    mlog_p = -log10(r$pval),
    # 4 genes have R2_D of -1e-4 to -1e-5 (fit tolerance); drawn at 0, not dropped.
    r2 = pmax(r$deviance_r2, 0),
    status = factor(as.character(r$gate_class), levels = names(gate_colours))
  )
}))
stopifnot(!anyNA(d_data$status))
d_cuts <- data.frame(
  time = pretty_time(times),
  cut = vapply(times, function(tm) {
    r <- ranked[[tm]]
    -log10(max(r$pval[r$pass_fdr %in% TRUE]))
  }, numeric(1))
)
floor_x <- -log10(1 / (n_perm + 1))
# Larger class drawn first, so the smaller ones stay visible (#55, #61).
d_data <- d_data[order(-table(d_data$status)[d_data$status]), , drop = FALSE]
d_data$x <- withr::with_seed(20261123L, jitter(d_data$mlog_p, amount = 0.04))
d_top <- max(d_data$r2)
d_counts <- stats::aggregate(r2 ~ time + status, d_data, length)
d_counts <- d_counts[order(d_counts$time, d_counts$status), ]
d_counts$y <- d_top * c(1, 0.86, 0.72)[as.integer(d_counts$status)]^2
short_status <- c("Fails BH" = "Fails", "Passes BH 0.05" = "Passes", "At the floor 1/(B+1)" = "At floor")
d_counts$label <- paste0(short_status[as.character(d_counts$status)], ": ", scales::comma(d_counts$r2))
d_lines <- data.frame(
  time = pretty_time(times[1L]),
  x = c(d_cuts$cut[1L], floor_x),
  y = d_top * 0.45^2,
  label = c("BH", "floor")
)

panel_d <- ggplot(d_data, aes(x, r2, colour = status)) +
  geom_point(size = 0.15, alpha = 0.35, stroke = 0) +
  geom_vline(data = d_cuts, aes(xintercept = cut), linetype = "dashed",
             colour = colour_mark, linewidth = 0.3) +
  geom_vline(xintercept = floor_x, linetype = "dotted", colour = colour_series,
             linewidth = 0.3) +
  geom_text(
    data = d_lines, aes(x, y, label = label), inherit.aes = FALSE,
    hjust = 1.1, size = mm_text(small), family = font_family, colour = colour_ink
  ) +
  geom_text(
    data = d_counts, aes(0, y, label = label), hjust = 0, vjust = 1,
    size = mm_text(small), family = font_family
  ) +
  scale_colour_manual(values = gate_colours) +
  scale_y_sqrt(breaks = c(0, 0.05, 0.2, 0.4, 0.6)) +
  facet_wrap(~time, nrow = 1) +
  labs(
    x = expression("\u2212" * log[10] ~ italic(p) ~ "(permuted Welch" ~ italic(G) * ")"),
    y = expression("Deviance share" ~ italic(R)[D]^2)
  ) +
  theme_slide()

# --- e: UpSet of the confirmed genes (hypothesis A) -------------------------

upset <- rds_repro(
  "output", "feature_selection", "deg", "figures", "fs_deg_superiority_upset.rds"
)
sets <- setdiff(names(upset), "gene")
stopifnot(length(sets) == 11L, all(vapply(upset[sets], is.logical, TRUE)))
names(upset)[match(sets, names(upset))] <- sub("^(\\d+)h: ", "\\1 h  ", sets)
sets <- sub("^(\\d+)h: ", "\\1 h  ", sets)
# As scripts/02_04_vis_deg.R:86 (min_size 5, 25 largest intersections).
panel_e <- ComplexUpset::upset(
  upset,
  intersect = sets,
  min_size = 5,
  n_intersections = 25,
  width_ratio = 0.25,
  name = NULL,
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      fill = colour_series, text = list(size = 0), text_colors = c(on_background = NA, on_bar = NA)
    )
  ),
  set_sizes = ComplexUpset::upset_set_size(
    geom = geom_bar(fill = colour_context)
  ) + scale_y_continuous(breaks = c(0, 2000), labels = c("0", "2k")),
  stripes = "white",
  matrix = ComplexUpset::intersection_matrix(geom = geom_point(size = 1.1)),
  themes = ComplexUpset::upset_modify_themes(list(
    default = theme_slide(),
    "Intersection size" = theme_slide() + theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      panel.grid.major.x = element_blank()
    ),
    intersections_matrix = theme_slide() + theme(
      axis.text.x = element_blank(), axis.title = element_blank(),
      axis.text.y = element_text(size = small), panel.grid = element_blank()
    ),
    overall_sizes = theme_slide() + theme(
      axis.text.y = element_blank(), axis.title.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
  ))
)

# --- f: curated enrichment themes -------------------------------------------

lolli <- rds_repro(
  "output", "feature_selection", "enrichment", "figures", "fs_enrichment_lollipop.rds"
)
keys <- c("time_point", "cell_type", "method", "database", "id")
f_data <- merge(lolli$themes[, c("theme", keys)], lolli$terms[, c(keys, "description", "padj")], by = keys)
stopifnot(nrow(f_data) == nrow(lolli$themes), all(f_data$cell_type %in% scored))
short_terms <- c(
  "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)" =
    "IGF transport and uptake by IGFBPs",
  "cellular response to leukemia inhibitory factor" = "cellular response to LIF"
)
f_data$term <- ifelse(
  f_data$description %in% names(short_terms),
  short_terms[f_data$description],
  sub("^HALLMARK_", "", f_data$description)
)
f_data$term <- ifelse(grepl("^[A-Z0-9_]+$", f_data$term), gsub("_", " ", f_data$term), f_data$term)
f_data$mlog <- -log10(f_data$padj)
f_data$label <- paste0(
  f_data$term, " (", type_number[f_data$cell_type], ", ", pretty_time(f_data$time_point), ")"
)
f_data <- f_data[order(f_data$theme, f_data$mlog), , drop = FALSE]
f_data$label <- factor(f_data$label, levels = unique(f_data$label))
f_data$theme <- factor(
  vapply(f_data$theme, function(s) paste(strwrap(s, 16), collapse = "\n"), ""),
  levels = vapply(unique(lolli$themes$theme), function(s) paste(strwrap(s, 16), collapse = "\n"), "")
)
db_colours <- c("GO BP" = colour_series, "Reactome" = colour_mark, "Hallmark" = "#5b8c5a")

panel_f <- ggplot(f_data, aes(mlog, label, colour = database)) +
  geom_segment(aes(x = 0, xend = mlog, yend = label), linewidth = 0.4) +
  geom_point(aes(shape = method), size = 1.3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = colour_ink, linewidth = 0.3) +
  scale_colour_manual(values = db_colours, name = NULL) +
  scale_shape_manual(values = c(ORA = 16, GSEA = 17), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = expression("\u2212" * log[10] ~ "BH-adjusted" ~ italic(p)), y = NULL) +
  theme_slide() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = small, face = "bold", colour = "#0e4d3a"),
    axis.text.y = element_text(size = small - 0.5),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.8, "mm"),
    legend.position = "inside",
    legend.position.inside = c(1, 0),
    legend.justification = c(1, 0),
    legend.box = "vertical",
    legend.spacing.y = unit(0, "mm")
  )

# --- g: pooled UMAP, one facet per time -------------------------------------

umap <- rds_repro(
  "output", "feature_selection", "tier_01", "figures", "fs_tier1_umap_pooled.rds"
)
scores <- umap$scores
stopifnot(all(scores$cell_type %in% scored))
# Grey context (all times) drawn from a seeded sample of 5,000 cells, the
# same in every facet; the coloured cells of each time are all drawn. Without
# the sample the SVG holds 3 x 20,437 grey points and grows past 40 MB.
bg_cells <- withr::with_seed(20261123L, sample.int(nrow(scores), 5000L))
g_bg <- do.call(rbind, lapply(times, function(tm) {
  data.frame(umap1 = scores$umap1[bg_cells], umap2 = scores$umap2[bg_cells], time = pretty_time(tm))
}))
g_fg <- data.frame(
  umap1 = scores$umap1, umap2 = scores$umap2,
  time = pretty_time(scores$time_point),
  cell_type = factor(scores$cell_type, levels = scored)
)
g_fg <- g_fg[order(-table(g_fg$cell_type)[g_fg$cell_type]), , drop = FALSE]
g_med <- stats::aggregate(cbind(umap1, umap2) ~ cell_type + time, g_fg, stats::median)
g_med$number <- type_number[as.character(g_med$cell_type)]
key_labels <- paste(type_number, scored)

panel_g <- ggplot(g_fg, aes(umap1, umap2)) +
  geom_point(data = g_bg, colour = "#d9d5cc", size = 0.3, stroke = 0) +
  geom_point(aes(colour = cell_type), size = 0.18, alpha = 0.6, stroke = 0) +
  geom_label(
    data = g_med, aes(label = number), size = mm_text(small), family = font_family,
    label.padding = unit(0.6, "pt"), linewidth = 0.15, colour = colour_ink,
    fill = "white", alpha = 0.85
  ) +
  scale_colour_manual(values = type_colours, labels = key_labels, drop = FALSE, name = NULL) +
  facet_wrap(~time, nrow = 1) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(ncol = 3, override.aes = list(size = 1.6, alpha = 1))) +
  theme_slide() +
  theme(
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = colour_grid, fill = NA, linewidth = 0.3),
    legend.position = "bottom",
    legend.key.spacing.y = unit(0, "mm"),
    legend.box.margin = margin(-3, 0, 0, 0)
  )

# --- h: Kendall tau between the mean profiles, Welch panel ------------------

cors <- read_repro("output", "deconvolution", "profile_correlation.csv")
cors <- cors[cors$gene_set == "welch", , drop = FALSE]
stopifnot(all(cors$type_a %in% scored), all(times %in% cors$time))
# Order of the book (R/04_01_vis_deconvolution.R:176): average linkage on
# 1 - Pearson, within each time.
h_order <- lapply(split(cors, cors$time), function(x) {
  types <- unique(x$type_a)
  m <- matrix(1, length(types), length(types), dimnames = list(types, types))
  m[cbind(x$type_a, x$type_b)] <- x$pearson
  types[stats::hclust(stats::as.dist(1 - m), method = "average")$order]
})
h_key <- function(type, time) paste(type, time, sep = "___")
h_levels <- unlist(lapply(names(h_order), function(tm) h_key(h_order[[tm]], tm)))
h_data <- data.frame(
  time = pretty_time(cors$time),
  a = factor(h_key(cors$type_a, cors$time), levels = h_levels),
  b = factor(h_key(cors$type_b, cors$time), levels = rev(h_levels)),
  kendall = cors$kendall
)
h_low <- floor(min(h_data$kendall) * 10) / 10
h_number <- function(x) as.character(type_number[sub("___.*$", "", x)])

panel_h <- ggplot(h_data, aes(a, b, fill = kendall)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sub("^0", "", formatC(kendall, format = "f", digits = 2)),
        colour = kendall > (h_low + 1) / 2),
    size = mm_text(5), family = font_family
  ) +
  scale_fill_gradient(low = "#f3f1ec", high = colour_series, limits = c(h_low, 1),
                      breaks = c(h_low, 1), name = expression("Kendall" ~ tau),
                      guide = guide_colourbar(display = "rectangles", nbin = 40)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = colour_ink), guide = "none") +
  scale_x_discrete(labels = h_number, expand = c(0, 0)) +
  scale_y_discrete(labels = h_number, expand = c(0, 0)) +
  facet_wrap(~time, nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme_slide() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = small),
    legend.title.position = "top",
    legend.key.width = unit(8, "mm"),
    legend.key.height = unit(1.6, "mm"),
    legend.box.margin = margin(-3, 0, 0, 0)
  )

# --- i: benchmark on the ladder of noise models -------------------------------

bench <- read_repro("output", "deconvolution", "benchmark_metrics.csv")
bench <- bench[bench$gene_set == "welch" & bench$family != "DeCovarT sweep", , drop = FALSE]
# Rungs of the ladder (vignettes/05-1-review-distribution-aware.qmd,
# fig-ladder; drawn in Figure 1c). MEAD and DECALS: rung 6, as on the ladder,
# although their point estimates use working weights.
comparator_rung <- c(
  NNLS = 1, RLM = 1, DWLS = 2, SCDC = 2, DeMixSC = 2,
  `RNA-Sieve` = 4, ReDeconv = 4, MEAD = 6, DECALS = 6, BayesPrism = 7
)
decovart_rows <- c(
  "DeCovarT shrinkage, bulk likelihood | welch" = "DeCovarT, raw, free scales (final)",
  "DeCovarT size-normalised shrinkage, free scales, bulk likelihood | welch" =
    "DeCovarT, size-normalised",
  "DeCovarT cell CV: factor, per type | welch" = "DeCovarT, factor, held-out cells"
)
comp <- bench[bench$family == "comparator", , drop = FALSE]
comp$base <- sub(" .*$", "", comp$algorithm)
stopifnot(setequal(unique(comp$base), names(comparator_rung)))
# Several versions of one method (DeMixSC): keep the one with the lowest
# mean RMSE over the nine standard-Wnt replicates.
pooled <- tapply(comp$rmse, comp$method, mean)
best <- vapply(split(comp$method, comp$base), function(m) {
  m <- unique(m)
  m[which.min(pooled[m])]
}, "")
comp <- comp[comp$method %in% best, , drop = FALSE]
comp$label <- ifelse(comp$base == "DeMixSC", "DeMixSC (weights)",
  ifelse(comp$base == "BayesPrism", "BayesPrism (Bayesian)",
    ifelse(comp$base == "MEAD", "MEAD", comp$algorithm)))
comp$rung <- comparator_rung[comp$base]
deco <- bench[bench$method %in% names(decovart_rows), , drop = FALSE]
deco$label <- unname(decovart_rows[deco$method])
deco$rung <- 4
i_data <- rbind(
  comp[, c("label", "rung", "time", "replicate", "rmse")],
  deco[, c("label", "rung", "time", "replicate", "rmse")]
)
stopifnot(all(table(i_data$label) == 9L), length(unique(i_data$label)) == 13L)
rung_names <- c(
  `1` = "1 OLS", `2` = "2 WLS", `4` = "4 QL", `6` = "6 GEE", `7` = "7 GLMM"
)
i_data$rung <- factor(rung_names[as.character(i_data$rung)], levels = rev(rung_names))
i_data$time <- pretty_time(i_data$time)
i_order <- names(sort(tapply(i_data$rmse, i_data$label, mean), decreasing = TRUE))
i_data$label <- factor(i_data$label, levels = i_order)
i_data$bayes <- grepl("Bayesian", i_data$label)
i_mean <- stats::aggregate(rmse ~ label + rung + time + bayes, i_data, mean)
i_bands <- data.frame(
  rung = factor(rung_names, levels = rev(rung_names)),
  alpha = c(0.05, 0.10, 0.17, 0.24, 0.31)
)

panel_i <- ggplot(i_data, aes(rmse, label, colour = bayes)) +
  geom_rect(
    data = i_bands, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = alpha),
    inherit.aes = FALSE, fill = "#5a86b8"
  ) +
  geom_point(size = 0.7, alpha = 0.6, stroke = 0) +
  geom_point(data = i_mean, size = 1.6) +
  scale_alpha_identity() +
  scale_colour_manual(values = c(`TRUE` = colour_mark, `FALSE` = colour_ink)) +
  scale_x_continuous(limits = c(0, NA), breaks = c(0, 0.2, 0.4), labels = c("0", "0.2", "0.4"),
                     expand = expansion(mult = c(0.02, 0.06))) +
  facet_grid(rung ~ time, scales = "free_y", space = "free_y") +
  labs(x = "RMSE of the proportions, standard-Wnt replicates", y = NULL) +
  theme_slide() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, size = small),
    panel.spacing.y = unit(0.6, "mm")
  )

# --- j: early against standard Wnt, DeCovarT final model --------------------

final_model <- "DeCovarT shrinkage, bulk likelihood | welch"
props <- read_repro("output", "deconvolution", "proportions_long.csv")
props <- props[props$method == final_model, , drop = FALSE]
stopifnot(
  all(props$cell_type %in% scored),
  all(table(props$time, props$arm, props$cell_type)[table(props$time, props$arm, props$cell_type) > 0] == 3L),
  all(abs(tapply(props$proportion, props$sample_id, sum) - 1) < 1e-6)
)
# Analysis fixed before any estimate was seen, in
# DeCovarT_reproducibility/scripts/06_01_plot_early_wnt_compositions.R
# (lines 67-69 and 247-275): centred log-ratio within each sample over the
# types of its reference, pseudo-count 1e-3; estimate = mean early CLR minus
# mean standard CLR; 95 % percentile bootstrap that resamples replicates
# within each arm (2,000 draws); range of the nine pairwise differences.
seed <- 20261101L
n_boot <- 2000L
pseudo_count <- 1e-3
log_p <- log(props$proportion + pseudo_count)
props$clr <- log_p - stats::ave(log_p, props$sample_id, FUN = mean)
resample <- function(v) v[sample.int(length(v), replace = TRUE)]
contrast_one <- function(clr_early, clr_standard) {
  boot <- replicate(n_boot, mean(resample(clr_early)) - mean(resample(clr_standard)))
  pairs <- outer(clr_early, clr_standard, "-")
  c(
    estimate = mean(clr_early) - mean(clr_standard),
    boot_lower = stats::quantile(boot, 0.025, names = FALSE),
    boot_upper = stats::quantile(boot, 0.975, names = FALSE),
    pair_min = min(pairs),
    pair_max = max(pairs)
  )
}
props$cell_type <- factor(props$cell_type, levels = scored)
set.seed(seed)
j_data <- do.call(rbind, lapply(split(props, list(props$time, props$cell_type), drop = TRUE), function(g) {
  data.frame(
    time = g$time[1L], cell_type = g$cell_type[1L],
    t(contrast_one(g$clr[g$arm == "early"], g$clr[g$arm == "standard"]))
  )
}))
rownames(j_data) <- NULL
utils::write.csv(
  j_data[order(j_data$time, j_data$cell_type), ],
  file.path(out_dir, "da_contrasts.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
log2_cols <- c("estimate", "boot_lower", "boot_upper", "pair_min", "pair_max")
j_data[log2_cols] <- j_data[log2_cols] / log(2) # natural-log CLR to log2 units
j_data$time <- pretty_time(j_data$time)
j_data$cell_type <- factor(j_data$cell_type, levels = rev(scored))
j_lim <- max(abs(unlist(j_data[log2_cols]))) * 1.05
j_arrows <- data.frame(
  time = pretty_time(times[1L]),
  x = c(-0.12, 0.12) * j_lim, xend = c(-0.55, 0.55) * j_lim,
  label = c("higher in standard", "higher in early"),
  hjust = c(1, 0)
)

panel_j <- ggplot(j_data, aes(y = cell_type, colour = cell_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = colour_ink, linewidth = 0.3) +
  geom_linerange(aes(xmin = pair_min, xmax = pair_max), colour = "#d9d3c7", linewidth = 2.2) +
  geom_errorbar(aes(xmin = boot_lower, xmax = boot_upper), width = 0.35,
                orientation = "y", linewidth = 0.45) +
  geom_point(aes(x = estimate), size = 1.6) +
  geom_segment(
    data = j_arrows, aes(x = x, xend = xend, y = Inf, yend = Inf), inherit.aes = FALSE,
    arrow = grid::arrow(length = unit(1.2, "mm"), type = "closed"),
    colour = colour_ink, linewidth = 0.3
  ) +
  scale_colour_manual(values = type_colours) +
  scale_x_continuous(limits = c(-j_lim, j_lim)) +
  facet_grid(time ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = expression(log[2] ~ "CLR ratio, early / standard Wnt"), y = NULL) +
  theme_slide() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    plot.margin = margin(9, 2, 1, 2)
  )
panel_j <- panel_j +
  geom_text(
    data = j_arrows, aes(x = xend, y = Inf, label = label, hjust = hjust),
    inherit.aes = FALSE, vjust = -1.2, size = mm_text(small), family = font_family,
    colour = colour_ink
  )

# --- Arrange the block (cowplot) and write it --------------------------------

block_width <- 194
block_height <- 195
labelled <- function(..., labels, rel_widths) {
  cowplot::plot_grid(
    ..., nrow = 1, labels = labels, rel_widths = rel_widths,
    label_fontfamily = font_family, label_size = 8, label_fontface = "bold",
    label_x = 0, label_y = 1, hjust = 0, vjust = 1, align = "h", axis = "tb"
  )
}
row_1 <- labelled(panel_c, panel_d, labels = c("c", "d"), rel_widths = c(0.40, 0.60))
row_2 <- labelled(panel_e, panel_f, labels = c("e", "f"), rel_widths = c(0.50, 0.50))
row_3 <- labelled(panel_g, panel_h, labels = c("g", "h"), rel_widths = c(0.60, 0.40))
row_4 <- labelled(panel_i, panel_j, labels = c("i", "j"), rel_widths = c(0.55, 0.45))
block <- cowplot::plot_grid(row_1, row_2, row_3, row_4, ncol = 1,
                            rel_heights = c(33, 66, 43, 53))
save_panel(block, "block", block_width, block_height)

print(j_data[, c("time", "cell_type", log2_cols)], digits = 2)
print(tapply(i_data$rmse, i_data$label, mean), digits = 3)
sessioninfo::session_info(pkgs = "attached")
