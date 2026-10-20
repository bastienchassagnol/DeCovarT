# Additive log-ratio maps, kept for figure scripts only.
# The package estimators use the isometric log-ratio chart.

additive_logistic <- function(rho) {
  scores <- c(rho, 0)
  weights <- exp(scores - max(scores))
  weights / sum(weights)
}

additive_log_ratio <- function(p) {
  num_cells <- length(p)
  log(p[seq_len(num_cells - 1L)] / p[[num_cells]])
}
