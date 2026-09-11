# FIG02_POSTPROCESS_ONLY=1 redraws density + performance.
# This helper redraws only performance books from saved artefacts.

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- "output/fig02"
PERF_DIR <- file.path(OUT_DIR, "performance_visualisations")
GGPLOT_RDS_DIR <- file.path(OUT_DIR, "ggplot_rds")
dir.create(PERF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GGPLOT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

artefacts <- read_simulation_artefacts(OUT_DIR, "bivariate", assemble = TRUE)
artefacts$theta <- readRDS(file.path(OUT_DIR, "bivariate_theta.rds"))
artefacts <- DeCovarT:::.attach_expected_fisher_wald(artefacts)
write_simulation_artefacts(
  artefacts,
  OUT_DIR,
  "bivariate",
  config = artefacts$config
)

skip_heavy <- identical(Sys.getenv("FIG02_SKIP_HEAVY", "0"), "1")
if (!isTRUE(skip_heavy)) {
  message("heatmaps…")
  save_bivariate_metric_heatmaps(
    artefacts,
    PERF_DIR,
    data_rds = GGPLOT_RDS_DIR
  )
}

message("raincloud…")
save_bivariate_raincloud_book(
  artefacts,
  file.path(PERF_DIR, "raincloud.pdf"),
  data_rds = GGPLOT_RDS_DIR
)

message("forest…")
save_bivariate_forest_book(
  artefacts,
  file.path(PERF_DIR, "forest.pdf"),
  data_rds = GGPLOT_RDS_DIR
)

message("similarity…")
save_bivariate_similarity_book(
  artefacts,
  file.path(PERF_DIR, "similarity.pdf"),
  data_rds = GGPLOT_RDS_DIR
)

if (!isTRUE(skip_heavy)) {
  message("solver dots…")
  save_bivariate_solver_dots_book(
    artefacts,
    file.path(PERF_DIR, "solver_dots.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
}

message("done")
print(list.files(PERF_DIR))
