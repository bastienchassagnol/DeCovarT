# Two-panel schematic of shared-plus-private mean signatures (no toy
# PCA panel). Run from the package root:
#   Rscript scripts/auxiliary/generate_mean_signature_schematic.R

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 to regenerate the mean-signature schematic.")
}

n_genes <- 12L
n_celltypes <- 3L
rho <- 0.35
mean_scale <- 10
shared_col <- "#9AA0A6"

gene_names <- paste0("g", seq_len(n_genes))
type_names <- paste0("Type ", seq_len(n_celltypes))

mu <- DeCovarT::generate_mean_signature_matrix(
  n_genes = n_genes,
  n_celltypes = n_celltypes,
  mean_scale = mean_scale,
  target_cosine = rho,
  gene_names = gene_names,
  celltype_names = type_names
)

block_size <- n_genes %/% n_celltypes
block_id <- rep(seq_len(n_celltypes), each = block_size)
heat_df <- data.frame(
  gene = factor(rep(rev(gene_names), n_celltypes), levels = rev(gene_names)),
  celltype = factor(
    rep(type_names, each = n_genes),
    levels = type_names
  ),
  value = as.vector(mu),
  block = factor(
    rep(block_id, n_celltypes),
    levels = seq_len(n_celltypes)
  )
)

panel_a <- ggplot2::ggplot(
  heat_df,
  ggplot2::aes(x = celltype, y = gene, fill = value)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::scale_fill_gradient(
    low = "#F4F1EA",
    high = "#1B4F72",
    name = "Mean"
  ) +
  ggplot2::labs(
    title = "A  Gene-by-cell-type mean signature",
    x = NULL,
    y = "Genes (private blocks)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

u <- rep(1 / sqrt(n_genes), n_genes)
v <- matrix(0, n_genes, n_celltypes)
for (j in seq_len(n_celltypes)) {
  idx <- which(block_id == j)
  v[idx, j] <- 1 / sqrt(length(idx))
}

vec_long <- rbind(
  data.frame(
    component = "shared u",
    gene = factor(gene_names, levels = gene_names),
    type = "shared",
    value = u
  ),
  do.call(
    rbind,
    lapply(seq_len(n_celltypes), function(j) {
      data.frame(
        component = paste0("private v", j),
        gene = factor(gene_names, levels = gene_names),
        type = type_names[[j]],
        value = v[, j]
      )
    })
  )
)
vec_long$component <- factor(
  vec_long$component,
  levels = c("shared u", paste0("private v", seq_len(n_celltypes)))
)

blend_df <- do.call(
  rbind,
  lapply(seq_len(n_celltypes), function(j) {
    data.frame(
      celltype = type_names[[j]],
      part = factor(
        c("shared", "private"),
        levels = c("shared", "private")
      ),
      weight = c(sqrt(rho), sqrt(1 - rho))
    )
  })
)

panel_b_top <- ggplot2::ggplot(
  vec_long,
  ggplot2::aes(x = gene, y = component, fill = value)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::scale_fill_gradient(
    low = "#F4F1EA",
    high = "#4A4A4A",
    name = "Loading"
  ) +
  ggplot2::labs(
    title = "B  Shared direction and orthogonal private blocks",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

panel_b_bot <- ggplot2::ggplot(
  blend_df,
  ggplot2::aes(x = celltype, y = weight, fill = part)
) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::scale_fill_manual(
    values = c(shared = shared_col, private = "#1B4F72"),
    labels = c(
      shared = expression(sqrt(rho) ~ u),
      private = expression(sqrt(1 - rho) ~ v[j])
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Blend weight"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.major.x = ggplot2::element_blank()
  )

out_path <- file.path(
  "vignettes",
  "figures",
  "fig_mean_signature_shared_private.png"
)

grDevices::png(
  out_path,
  width = 11.2,
  height = 5.2,
  units = "in",
  res = 180
)
grid::grid.newpage()
layout_mat <- grid::grid.layout(
  nrow = 2,
  ncol = 2,
  widths = grid::unit(c(1.15, 1), "null"),
  heights = grid::unit(c(1.35, 0.85), "null")
)
grid::pushViewport(grid::viewport(layout = layout_mat))
print(
  panel_a,
  vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1)
)
print(panel_b_top, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(panel_b_bot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
grDevices::dev.off()

message("Wrote: ", out_path)
