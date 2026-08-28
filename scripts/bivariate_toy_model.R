###########################################################################
###########################################################################
###                                                                     ###
###                   BIVARIATE TOY MODEL FIGURES                       ###
###                                                                     ###
###########################################################################
###########################################################################

# Post-process benchmark RDS artefacts (tables and JOBIM heatmaps).
# To regenerate simulations, run scripts/run_bivariate_toy_benchmark.R first.
#
# Scenario configuration: scripts/configure_bivariate_toy_scenarios.R
# Wrapper API: run_simulation_benchmark()

library(dplyr)
library(tinytable)

results_dir <- file.path("simulations", "results")
bivariate_simulation <- readRDS(
  file.path(results_dir, "bivariate_scenario.rds")
)
bivariate_configuration <- readRDS(
  file.path(results_dir, "complete_bivariate_configuration.rds")
)

reduced_bivariate_configuration <- bivariate_configuration |>
  dplyr::mutate(
    ID = factor(
      ID,
      levels = unique(bivariate_configuration$ID),
      ordered = TRUE
    )
  ) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(
    Entropy = entropy,
    OVL = mean(overlap),
    Proportions = purrr::map_chr(
      true_parameters,
      ~ paste(.x$p, collapse = " / ")
    ),
    Means = purrr::map_chr(
      true_parameters,
      ~ paste0(
        "(",
        paste0(.x$mu[, 1], collapse = ","),
        ");(",
        paste0(.x$mu[, 2], collapse = ","),
        ")"
      )
    ),
    Variance = purrr::map_chr(
      true_parameters,
      ~ paste(
        c(.x$sigma[1, 1, 1], .x$sigma[2, 2, 1]),
        collapse = " / "
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::distinct()

saveRDS(
  reduced_bivariate_configuration,
  file.path(results_dir, "reduced_bivariate_configuration.rds")
)

tinytable::tt(
  reduced_bivariate_configuration,
  caption = paste(
    "The 8 general scenarios tested to compare the performance of DeCovarT",
    "vs standard linear deconvolution model"
  )
) |>
  tinytable::style_tt(j = 2:6, align = "c")

reduced_bivariate_simulation <- bivariate_simulation |>
  dplyr::mutate(
    algorithm = factor(
      algorithm,
      levels = c(
        "nnls",
        "lsei",
        "gradient",
        "hessian",
        "DeCoVarT",
        "optim",
        "barrier",
        "SA",
        "LBFGS",
        "Newton-Raphson",
        "Marquardt-Levenberg"
      )
    )
  ) |>
  dplyr::mutate(
    algorithm = forcats::fct_recode(algorithm, Levenberg = "DeCoVarT") |>
      forcats::fct_relevel()
  ) |>
  dplyr::mutate(dplyr::across(where(is.numeric), signif, digits = 4)) |>
  dplyr::select(
    -dplyr::any_of(c(
      "proportions",
      "true_parameters",
      "variance",
      "overlap",
      "model_coef_determination",
      "model_coef_determination_adjusted",
      "model_cor",
      "entropy",
      "centroids"
    ))
  ) |>
  dplyr::relocate(ID, .before = correlation_celltype1) |>
  dplyr::rename_with(
    ~ gsub("celltype_", "p", .x),
    .cols = dplyr::matches("celltype_")
  )

data_dir <- file.path("data", "bivariate")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  reduced_bivariate_simulation,
  file.path(data_dir, "bivariate_parameters.rds")
)

if (file.exists(file.path(results_dir, "high_overlap_version2.rds"))) {
  highly_overlapping_simulation <- readRDS(
    file.path(results_dir, "high_overlap_version2.rds")
  )

  data <- highly_overlapping_simulation |>
    dplyr::filter(proportions == "balanced" & variance == "homoscedastic") |>
    plot_correlation_Heatmap(score_variable = "model_mse")
  data1 <- data$lsei@matrix
  data2 <- data$`tricky optim`@matrix

  common_min <- min(c(data1, data2))
  common_max <- max(c(data1, data2))
  col_fun <- circlize::colorRamp2(c(common_min, common_max), c("blue", "red"))

  global_heatmap <- ComplexHeatmap::Heatmap(
    data1,
    col = col_fun,
    heatmap_legend_param = list(title = "MSE"),
    row_title = "Corr cell type 1",
    cluster_rows = FALSE,
    row_names_gp = grid::gpar(fontsize = 8),
    row_labels = colnames(data1),
    row_title_gp = grid::gpar(fontsize = 10),
    column_names_rot = 0,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 8),
    column_labels = colnames(data1),
    width = grid::unit(8, "cm"),
    height = grid::unit(8, "cm"),
    column_title_gp = grid::gpar(fontsize = 10),
    column_title = "Corr cell type 2"
  ) +
    ComplexHeatmap::Heatmap(
      data2,
      col = col_fun,
      show_heatmap_legend = FALSE,
      heatmap_legend_param = list(title = "MSE"),
      row_title = "Corr cell type 1",
      cluster_rows = FALSE,
      row_names_gp = grid::gpar(fontsize = 8),
      row_labels = colnames(data2),
      row_title_gp = grid::gpar(fontsize = 10),
      column_names_rot = 0,
      cluster_columns = FALSE,
      column_names_gp = grid::gpar(fontsize = 8),
      column_labels = colnames(data2),
      width = grid::unit(8, "cm"),
      height = grid::unit(8, "cm"),
      column_title_gp = grid::gpar(fontsize = 10),
      column_title = "Corr cell type 2"
    )
}
