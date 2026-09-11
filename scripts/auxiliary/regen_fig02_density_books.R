# One-off: regenerate fig02 density / loglik books from saved artefacts.
# Rscript --no-save --no-restore scripts/auxiliary/regen_fig02_density_books.R
# Needs a working OpenGL display for rgl PDF snapshots (not RGL_USE_NULL).

devtools::load_all(".", quiet = TRUE)

OUT_DIR <- "output/fig02"
DENSITY_DIR <- file.path(OUT_DIR, "density_visualisations")
GGPLOT_RDS_DIR <- file.path(OUT_DIR, "ggplot_rds")
dir.create(DENSITY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GGPLOT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

cfg <- DeCovarT:::.relevel_scenario_table(
  readRDS(file.path(OUT_DIR, "bivariate_config.rds"))
)
theta_tbl <- DeCovarT:::.relevel_scenario_table(
  readRDS(file.path(OUT_DIR, "bivariate_theta.rds"))
)
message(
  "pages: ",
  nrow(dplyr::distinct(cfg, proportions, variance, centroids))
)

old <- file.path(DENSITY_DIR, "loglik_surface.pdf")
if (file.exists(old)) {
  unlink(old)
  message("removed stale loglik_surface.pdf")
}

message("purified…")
save_bivariate_purified_density_book(
  cfg,
  theta_tbl,
  file.path(DENSITY_DIR, "purified_density.pdf"),
  n = 600L,
  data_rds = GGPLOT_RDS_DIR
)
message("bulk…")
save_bivariate_bulk_density_book(
  cfg,
  theta_tbl,
  file.path(DENSITY_DIR, "bulk_density.pdf"),
  n = 800L,
  data_rds = GGPLOT_RDS_DIR
)
message("surface_p…")
save_bivariate_loglik_surface_p_book(
  cfg,
  theta_tbl,
  file.path(DENSITY_DIR, "loglik_surface_p.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
message("ilr profile…")
save_bivariate_loglik_ilr_profile_book(
  cfg,
  theta_tbl,
  file.path(DENSITY_DIR, "loglik_ilr_profile.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
skip_rgl <- identical(Sys.getenv("FIG02_SKIP_RGL", "0"), "1")
if (isTRUE(skip_rgl)) {
  message("skipping rgl (FIG02_SKIP_RGL=1)")
} else {
  message("rgl pdf+html…")
  save_bivariate_loglik_rgl_book(
    cfg,
    theta_tbl,
    file.path(DENSITY_DIR, "loglik_rgl.pdf"),
    html_file = file.path(DENSITY_DIR, "loglik_rgl.html")
  )
}
message("done")
print(list.files(DENSITY_DIR))
