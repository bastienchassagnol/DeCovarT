synthethic_networks_structure <- simulate_hierarchical_grn_moments(
  n_expressed_genes = 5L,
  mean_lower_expressed = 5,
  mean_upper_expressed = 20,
  mean_lower_background = 1,
  mean_upper_background = 3,
  library_size = 0,
  alpha = 2,
  precision_shift = 0.1,
  precision_scale = 0.3,
  child_perturbation_sd = 0.1,
  graph_model = "power_law",
  graph_params = list(power = 1, edges_per_node = 2)
)
