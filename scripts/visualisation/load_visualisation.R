# Load the figure-book helpers of the simulation pipelines (fig02, fig03).
#
# The four 05_0x scripts below call one another and the auxiliary helpers,
# so they are always sourced together, after the package itself:
#   devtools::load_all(".")
#   source(file.path("scripts", "visualisation", "load_visualisation.R"))

for (helper in c(
  file.path("scripts", "auxiliary", "benchmark_artefacts.R"),
  file.path("scripts", "auxiliary", "alr_coordinates.R"),
  file.path("scripts", "auxiliary", "display_metrics.R"),
  file.path("scripts", "visualisation", "05_03_bivariate_visualisation.R"),
  file.path("scripts", "visualisation", "05_04_hybrid_latent_visualisation.R"),
  file.path("scripts", "visualisation", "05_05_convergence_visualisation.R"),
  file.path("scripts", "visualisation", "05_06_convergence_books.R")
)) {
  source(helper)
}
rm(helper)
