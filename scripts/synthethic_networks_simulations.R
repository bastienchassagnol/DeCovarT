library(ComplexHeatmap)
library(circlize)
library(visNetwork)
library(ggplot2)
library(tidyr)

# devtools::load_all()

set.seed(42)
synthethic_networks_structure <- simulate_hierarchical_grn_moments(
  n_expressed_genes = 100L,
  mean_lower_expressed = 5,
  mean_upper_expressed = 20,
  mean_lower_background = 1,
  mean_upper_background = 3,
  library_size = 2,
  alpha = 2,
  precision_shift = 0.1,
  precision_scale = 0.3,
  child_perturbation_sd = 0.1,
  graph_model = "power_law",
  graph_params = list(power = 1, edges_per_node = 2)
)

# ---- 1. ComplexHeatmap: parent vs child mean profiles side by side --------

parent_means <- synthethic_networks_structure$parent_parameters$mean_profiles
child_means <- synthethic_networks_structure$child_parameters$mean_profiles

shared_col_fun <- circlize::colorRamp2(
  seq(0, max(parent_means, child_means), length.out = 9),
  rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
)

ht_parents <- ComplexHeatmap::Heatmap(
  t(parent_means),
  name = "Gene Expression",
  col = shared_col_fun,
  column_title = "Parent cell types",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_title = "Genes",
  column_names_rot = 45,
  width = unit(3, "cm")
)

ht_children <- Heatmap(
  t(child_means),
  name = "Expression ",
  col = shared_col_fun,
  column_title = "Child cell types",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 45,
  width = unit(6, "cm")
)

draw(ht_parents + ht_children, column_title = "Mean expression profiles")


# ---- 2. visNetwork: adjacency & precision-weighted graphs ----------------

adj_mat <- synthethic_networks_structure$graph_structure$adjacency_matrix
prec_mat <- synthethic_networks_structure$graph_structure$normalised_precision

edges_idx <- which(adj_mat == 1 & upper.tri(adj_mat), arr.ind = TRUE)
parent_assignment <- ifelse(
  seq_len(nrow(adj_mat)) <= ncol(parent_means) / 2,
  "parent_1", "parent_2"
)
nodes <- data.frame(
  id = seq_len(nrow(adj_mat)),
  label = rownames(adj_mat),
  group = parent_assignment,
  title = paste0(
    rownames(adj_mat), "<br>Block: ", parent_assignment
  ),
  stringsAsFactors = FALSE
)

## ---- i. plot unweighted graph ----
edges_unweighted <- data.frame(
  from = edges_idx[, "row"],
  to = edges_idx[, "col"],
  stringsAsFactors = FALSE
)
visNetwork(nodes, edges_unweighted, main = "Adjacency graph") |>
  visGroups(groupname = "parent_1", color = "#E41A1C") |>
  visGroups(groupname = "parent_2", color = "#377EB8") |>
  visEdges(color = list(color = "#AAAAAA", highlight = "#333333")) |>
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
  visLegend()


## ---- i. plot weighted graph ----
prec_weights <- abs(prec_mat[edges_idx])
scaled_weights <- 1 + 9 * (prec_weights - min(prec_weights)) /
  (max(prec_weights) - min(prec_weights) + .Machine$double.eps)

edges_weighted <- data.frame(
  from = edges_idx[, "row"],
  to = edges_idx[, "col"],
  value = scaled_weights,
  title = paste0("precision: ", round(prec_weights, 4)),
  stringsAsFactors = FALSE
)

visNetwork(nodes, edges_weighted, main = "Precision-weighted graph") |>
  visGroups(groupname = "parent_1", color = "#E41A1C") |>
  visGroups(groupname = "parent_2", color = "#377EB8") |>
  visEdges(
    color = list(color = "#AAAAAA", highlight = "#333333"),
    scaling = list(min = 1, max = 8)
  ) |>
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
  visLegend()


# ---- 3. Mean-variance facet plot for child marginal variances ------------

child_var_data <- do.call(rbind, lapply(
  seq_len(nrow(child_means)),
  function(k) {
    mu_k <- child_means[k, ]
    cov_k <- synthethic_networks_structure$child_parameters$covariance_matrices[,, k]
    data.frame(
      child = rownames(child_means)[k],
      mean = mu_k,
      variance = diag(cov_k)
    )
  }
))

ggplot(child_var_data, aes(x = mean, y = variance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~child, nrow = 2, ncol = 2) +
  labs(x = expression(mu[g]), y = expression(sigma[g]^2)) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

