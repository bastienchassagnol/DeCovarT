# Numbers quoted in the Results of vignettes/R3-bivariate-single-cell.qmd.
# Usage: Rscript scripts/fig05_summary_numbers.R [output/fig05/fig05_benchmark.rds]
#
# Monte Carlo standard errors: RMSE by the delta method,
# sd(e^2) / (2 RMSE sqrt(N)); proportions (coverage, rejection) by the
# Wilson interval of coverage_mc_interval().
suppressMessages({
  library(dplyr)
  library(tidyr)
})
options(width = 200, pillar.sigfig = 3, dplyr.summarise.inform = FALSE)
args <- commandArgs(trailingOnly = TRUE)
file <- if (length(args)) args[[1]] else "output/fig05/fig05_benchmark.rds"
b <- readRDS(file)
cfg <- b$config
cat("meta:", b$meta$date, b$meta$git_commit, "N =", b$meta$n_replicates, "\n")
show <- function(x, digits = 3) print(as.data.frame(x), digits = digits)
contrast_bin <- function(x) factor(round(x, 1))
p_true <- setNames(vapply(b$theta_true, function(th) th$p[[1L]], 0), names(b$theta_true))

opt <- b$optimisation
opt$p1_true <- p_true[opt$ID]
opt$err <- opt$celltype_1 - opt$p1_true
opt$contrast <- abs(opt$correlation_celltype1 - opt$correlation_celltype2)

cat("\n== Fisher information of p1 at p*: covariance share, single cell and population ==\n")
show(
  cfg |>
    mutate(
      f_cov_sc = fisher_cov_sc / (fisher_mean_sc + fisher_cov_sc),
      f_cov_pop = fisher_cov_pop / (fisher_mean_pop + fisher_cov_pop),
      equal = correlation_celltype1 == correlation_celltype2
    ) |>
    summarise(
      median_f_cov_sc = median(f_cov_sc),
      median_f_cov_pop = median(f_cov_pop),
      f_cov_sc_equal_R = max(f_cov_sc[equal]),
      f_cov_pop_equal_R = median(f_cov_pop[equal]),
      f_cov_sc_opposite = f_cov_sc[correlation_celltype1 == -0.8 & correlation_celltype2 == 0.8],
      f_cov_pop_opposite = f_cov_pop[correlation_celltype1 == -0.8 & correlation_celltype2 == 0.8],
      I_sc_rho00 = (fisher_mean_sc + fisher_cov_sc)[correlation_celltype1 == 0 & correlation_celltype2 == 0],
      I_pop_rho00 = (fisher_mean_pop + fisher_cov_pop)[correlation_celltype1 == 0 & correlation_celltype2 == 0],
      .by = c(centroids, proportions)
    )
)

cat("\n== RMSE of p1 per scenario: median over the 81 correlation pairs [median MC SE] ==\n")
rm <- opt |>
  summarise(
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    mcse = sd(err^2, na.rm = TRUE) / (2 * sqrt(mean(err^2, na.rm = TRUE)) * sqrt(sum(is.finite(err)))),
    bias = mean(err, na.rm = TRUE),
    fail = mean(!numerical_converged),
    face = mean(n_zero > 0, na.rm = TRUE),
    .by = c(ID, centroids, proportions, algorithm, contrast, overlap,
            correlation_celltype1, correlation_celltype2)
  )
show(
  rm |>
    summarise(
      value = sprintf("%.3f [%.4f]", median(rmse), median(mcse)),
      .by = c(centroids, proportions, algorithm)
    ) |>
    pivot_wider(names_from = algorithm, values_from = value)
)
cat("\nnumerical failures (max over scenarios):\n")
show(rm |> summarise(fail = max(fail), .by = algorithm))
cat("\nfraction of bulks with an estimate on a face (median over scenarios):\n")
show(
  rm |>
    summarise(face = median(face), .by = c(centroids, proportions, algorithm)) |>
    pivot_wider(names_from = algorithm, values_from = face)
)
cat("\nbias of p1 (median over scenarios):\n")
show(
  rm |>
    summarise(bias = median(bias), .by = c(centroids, proportions, algorithm)) |>
    pivot_wider(names_from = algorithm, values_from = bias)
)

cat("\n== RMSE ratio to lsei on the same bulks ==\n")
ratio <- rm |>
  left_join(
    rm |> filter(algorithm == "lsei") |> select(ID, rmse_lsei = rmse),
    by = "ID"
  ) |>
  mutate(ratio = rmse / rmse_lsei) |>
  filter(algorithm != "lsei")
cat("median ratio by |rho1 - rho2|, Scalar MLE:\n")
show(
  ratio |>
    filter(algorithm == "Scalar MLE") |>
    summarise(ratio = median(ratio), .by = c(centroids, proportions, contrast)) |>
    mutate(contrast = contrast_bin(contrast)) |>
    summarise(ratio = median(ratio), .by = c(centroids, proportions, contrast)) |>
    pivot_wider(names_from = contrast, values_from = ratio)
)
cat("\nshare of scenarios with ratio < 0.95 (gain) and > 1.05 (loss):\n")
show(
  ratio |>
    summarise(
      gain = mean(ratio < 0.95),
      loss = mean(ratio > 1.05),
      min = min(ratio),
      max = max(ratio),
      .by = c(centroids, proportions, algorithm)
    ) |>
    arrange(algorithm, centroids, proportions)
)
cat("\nmedian ratio (Scalar MLE) by overlap tercile x |rho1 - rho2| tercile, within centroids:\n")
show(
  ratio |>
    filter(algorithm == "Scalar MLE") |>
    mutate(
      overlap_t = ntile(overlap, 3),
      contrast_t = ntile(contrast, 3),
      .by = centroids
    ) |>
    summarise(ratio = median(ratio), .by = c(centroids, overlap_t, contrast_t)) |>
    pivot_wider(names_from = contrast_t, values_from = ratio, names_prefix = "contrast_t") |>
    arrange(centroids, overlap_t)
)

cat("\n== Package solvers against the scalar reference MLE ==\n")
ref <- opt |>
  filter(algorithm == "Scalar MLE") |>
  select(ID, sample_id, p1_ref = celltype_1, regret_ref = loglik_regret)
agree <- opt |>
  filter(algorithm %in% c("Newton-Raphson", "L-BFGS-B", "Marquardt-Levenberg")) |>
  inner_join(ref, by = c("ID", "sample_id")) |>
  mutate(
    gap = abs(celltype_1 - p1_ref),
    # loglik(reference) - loglik(solver) = regret(solver) - regret(reference)
    shortfall = loglik_regret - regret_ref
  )
show(
  agree |>
    summarise(
      median_gap = median(gap, na.rm = TRUE),
      gap_gt_0.01 = mean(gap > 0.01, na.rm = TRUE),
      shortfall_gt_1e3 = mean(shortfall > 1e-3, na.rm = TRUE),
      shortfall_lt_m1e6 = mean(shortfall < -1e-6, na.rm = TRUE),
      .by = c(centroids, proportions, algorithm)
    ) |>
    arrange(algorithm, centroids, proportions)
)
cat("\nmedian time per bulk (ms):\n")
show(opt |> summarise(ms = 1000 * median(elapsed_sec), .by = algorithm))

cat("\n== Coverage of nominal 95% intervals for p1, mean over scenarios [width] ==\n")
iv <- b$intervals |>
  filter(parameter == "celltype_1", n_intervals > 0)
show(
  iv |>
    filter(algorithm == "Newton-Raphson", ci_method != "bootstrap") |>
    summarise(
      value = sprintf("%.3f [%.3f]", mean(coverage, na.rm = TRUE), mean(mean_width, na.rm = TRUE)),
      .by = c(centroids, proportions, ci_method)
    ) |>
    pivot_wider(names_from = ci_method, values_from = value)
)
cat("\nWald and ILR under the other solvers (mean coverage):\n")
show(
  iv |>
    filter(algorithm != "Newton-Raphson") |>
    summarise(coverage = mean(coverage, na.rm = TRUE), .by = c(centroids, proportions, algorithm, ci_method)) |>
    pivot_wider(names_from = c(algorithm, ci_method), values_from = coverage)
)
cat("\nintervals at the boundary or undefined (share of bulks without an interval, Newton):\n")
show(
  iv |>
    filter(algorithm %in% c("Newton-Raphson", "L-BFGS-B")) |>
    summarise(
      missing = mean(1 - n_intervals / b$meta$n_replicates),
      .by = c(proportions, algorithm, ci_method)
    ) |>
    pivot_wider(names_from = c(algorithm, ci_method), values_from = missing)
)
cat("\nbootstrap scenarios (Newton-Raphson): coverage [lower, upper] and width per method\n")
boot_ids <- unique(iv$ID[iv$ci_method == "bootstrap"])
show(
  iv |>
    filter(ID %in% boot_ids, algorithm == "Newton-Raphson") |>
    mutate(rho = paste0("(", correlation_celltype1, ", ", correlation_celltype2, ")")) |>
    transmute(
      proportions,
      rho,
      ci_method,
      value = sprintf("%.3f [%.3f, %.3f] w=%.3f", coverage, coverage_lower, coverage_upper, mean_width)
    ) |>
    pivot_wider(names_from = ci_method, values_from = value) |>
    arrange(proportions, rho)
)

cat("\n== Paired comparison: single-cell law against population law, same draws ==\n")
paired <- b$paired |>
  summarise(
    rmse_sc = sqrt(mean((p1_single_cell - p1_true)^2)),
    rmse_pop = sqrt(mean((p1_population - p1_true)^2)),
    # MC SE of the log ratio from the paired squared errors (delta method)
    .by = c(ID, method)
  ) |>
  left_join(cfg[c("ID", "centroids", "proportions", "correlation_celltype1", "correlation_celltype2", "covariance_contrast")], by = "ID") |>
  mutate(log2_ratio = log2(rmse_sc / rmse_pop))
show(
  paired |>
    summarise(
      rmse_sc = median(rmse_sc),
      rmse_pop = median(rmse_pop),
      median_log2_ratio = median(log2_ratio),
      share_sc_better = mean(log2_ratio < 0),
      log2_ratio_equal_R = median(log2_ratio[correlation_celltype1 == correlation_celltype2]),
      log2_ratio_opposite = log2_ratio[correlation_celltype1 == -0.8 & correlation_celltype2 == 0.8],
      .by = c(centroids, proportions, method)
    ) |>
    arrange(method, centroids, proportions)
)
cat("\nmedian log2 ratio (Scalar MLE) by |rho1 - rho2|:\n")
show(
  paired |>
    filter(method == "Scalar MLE") |>
    mutate(contrast = contrast_bin(covariance_contrast)) |>
    summarise(log2_ratio = median(log2_ratio), .by = c(centroids, proportions, contrast)) |>
    pivot_wider(names_from = contrast, values_from = log2_ratio)
)

cat("\n== Boundary likelihood-ratio test of p2 = 0 ==\n")
wilson <- function(x) {
  ci <- DeCovarT::coverage_mc_interval(x)
  sprintf("%.4f [%.4f, %.4f]", ci$coverage, ci$lower, ci$upper)
}
cat("power at the design composition (chi-bar-square, 5%), pooled over 81 pairs:\n")
show(
  b$boundary_lrt$power |>
    left_join(cfg[c("ID", "centroids", "proportions")], by = "ID") |>
    summarise(rate = wilson(p_value < 0.05), .by = c(centroids, proportions))
)
cat("\nsize at p = (1, 0), N =", b$boundary_lrt$settings$n_null, ", B =", b$boundary_lrt$settings$b_null, ":\n")
size <- b$boundary_lrt$size |>
  left_join(b$boundary_lrt$null_config[c("ID", "centroids", "correlation_celltype1", "correlation_celltype2")], by = "ID")
show(
  size |>
    summarise(
      chibar = wilson(p_chibar < 0.05),
      monte_carlo = wilson(p_monte_carlo < 0.05),
      share_T_zero = mean(statistic < 1e-8),
      .by = c(ID, centroids, correlation_celltype1, correlation_celltype2)
    )
)
cat("pooled: chi-bar", wilson(size$p_chibar < 0.05), " Monte Carlo", wilson(size$p_monte_carlo < 0.05), "\n")
check <- size |> filter(is.finite(package_statistic)) |> mutate(diff = statistic - package_statistic)
cat(
  "\npackage lrt_single_cell() against the scalar statistic on", nrow(check), "null bulks:",
  "\n  |diff| < 1e-4:", round(mean(abs(check$diff) < 1e-4), 4),
  " scalar larger by > 1e-4:", round(mean(check$diff > 1e-4), 4),
  " package larger by > 1e-4:", round(mean(check$diff < -1e-4), 4),
  " max diff:", signif(max(abs(check$diff)), 3), "\n"
)
