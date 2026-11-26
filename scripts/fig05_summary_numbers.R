# Numbers quoted in the Results of vignettes/R3-bivariate-single-cell.qmd.
# Usage: Rscript scripts/fig05_summary_numbers.R [output/fig05/fig05_benchmark.rds]
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
opt <- b$optimisation
cts <- c("celltype_1", "celltype_2")
s <- as.integer(sub("sample_", "", opt$sample_id))
truth <- t(mapply(function(id, k) b$data[[id]]$truth[cts, k], opt$ID, s))
opt$err <- opt$celltype_1 - truth[, 1]
opt$p1_true <- truth[, 1]
opt$arm <- ifelse(
  grepl("oracle n", opt$algorithm),
  "oracle",
  ifelse(grepl("nominal n", opt$algorithm), "nominal", "none")
)
opt$method <- sub(" \\((oracle|nominal) n\\)$", "", opt$algorithm)
opt$corner <- paste(opt$correlation_celltype1, opt$correlation_celltype2)

cat(
  "\n== RMSE of p1 (sqrt mean sq error over bulks), median over 81 scenarios ==\n"
)
rm <- opt |>
  summarise(
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    bias = mean(err, na.rm = TRUE),
    fail = mean(!numerical_converged),
    .by = c(ID, proportions, n_cells, method, arm, corner)
  )
print(
  as.data.frame(
    rm |>
      summarise(
        rmse = median(rmse),
        .by = c(proportions, n_cells, method, arm)
      ) |>
      pivot_wider(names_from = c(method), values_from = rmse) |>
      arrange(proportions, n_cells, arm)
  ),
  digits = 3
)

cat("\n== bias of p1 (median over scenarios) ==\n")
print(
  as.data.frame(
    rm |>
      summarise(
        bias = median(bias),
        .by = c(proportions, n_cells, method, arm)
      ) |>
      pivot_wider(names_from = c(method), values_from = bias) |>
      arrange(proportions, n_cells, arm)
  ),
  digits = 2
)

cat("\n== corners: RMSE ==\n")
print(
  as.data.frame(
    rm |>
      filter(
        corner %in% c("0 0", "-0.8 0.8", "0.8 0.8", "-0.8 -0.8"),
        arm != "nominal"
      ) |>
      select(proportions, n_cells, corner, method, rmse) |>
      pivot_wider(names_from = method, values_from = rmse) |>
      arrange(proportions, n_cells, corner)
  ),
  digits = 3
)

cat("\n== failures ==\n")
print(
  as.data.frame(
    opt |>
      summarise(
        fail = mean(!numerical_converged),
        local_max = mean(local_maximum, na.rm = TRUE),
        n_zero = mean(n_zero > 0, na.rm = TRUE),
        kkt_med = median(kkt_residual, na.rm = TRUE),
        kkt_q99 = quantile(kkt_residual, 0.99, na.rm = TRUE),
        sec = median(elapsed_sec),
        .by = c(method, arm)
      )
  ),
  digits = 3
)

cat("\n== optimiser agreement (oracle arm): |NR - LBFGS|, |NR - ML| ==\n")
w <- opt |>
  filter(
    arm == "oracle",
    method %in% c("Newton-Raphson", "L-BFGS-B", "Marquardt-Levenberg")
  ) |>
  select(ID, sample_id, proportions, n_cells, method, celltype_1) |>
  pivot_wider(names_from = method, values_from = celltype_1)
print(
  as.data.frame(
    w |>
      summarise(
        lb_med = median(abs(`Newton-Raphson` - `L-BFGS-B`), na.rm = TRUE),
        lb_q95 = quantile(
          abs(`Newton-Raphson` - `L-BFGS-B`),
          0.95,
          na.rm = TRUE
        ),
        lb_gt01 = mean(abs(`Newton-Raphson` - `L-BFGS-B`) > 0.01, na.rm = TRUE),
        ml_med = median(
          abs(`Newton-Raphson` - `Marquardt-Levenberg`),
          na.rm = TRUE
        ),
        ml_q95 = quantile(
          abs(`Newton-Raphson` - `Marquardt-Levenberg`),
          0.95,
          na.rm = TRUE
        ),
        ml_gt01 = mean(
          abs(`Newton-Raphson` - `Marquardt-Levenberg`) > 0.01,
          na.rm = TRUE
        ),
        lb_face = mean(`L-BFGS-B` %in% c(0, 1)),
        .by = c(proportions, n_cells)
      )
  ),
  digits = 3
)

cat("\n== loglik regret by method (oracle): median, share > 1e-3 ==\n")
print(
  as.data.frame(
    opt |>
      filter(arm != "nominal") |>
      summarise(
        reg_med = median(loglik_regret, na.rm = TRUE),
        reg_pos = mean(loglik_regret > 1e-3, na.rm = TRUE),
        .by = c(method, n_cells)
      )
  ),
  digits = 3
)

iv <- b$intervals |> filter(parameter == "celltype_1")
iv$arm <- ifelse(
  grepl("oracle n", iv$algorithm),
  "oracle",
  ifelse(grepl("nominal n", iv$algorithm), "nominal", "none")
)
iv$method <- sub(" \\((oracle|nominal) n\\)$", "", iv$algorithm)
iv <- left_join(
  iv,
  cfg[c(
    "ID",
    "proportions",
    "n_cells",
    "correlation_celltype1",
    "correlation_celltype2"
  )] |>
    rename_with(~ paste0(.x, "_c"), -ID),
  by = "ID"
)
cat("\n== coverage of p1, mean over scenarios ==\n")
print(
  as.data.frame(
    iv |>
      filter(!ci_method %in% c("bootstrap", "bootstrap_nprior")) |>
      summarise(
        cov = mean(coverage, na.rm = TRUE),
        .by = c(proportions_c, n_cells_c, method, arm, ci_method)
      ) |>
      pivot_wider(names_from = c(arm), values_from = cov) |>
      arrange(method, ci_method, proportions_c, n_cells_c)
  ),
  digits = 3
)
cat("\n== width ==\n")
print(
  as.data.frame(
    iv |>
      filter(
        !ci_method %in% c("bootstrap", "bootstrap_nprior"),
        method %in% c("Newton-Raphson", "RNA-Sieve", "DeCovarT diagonal")
      ) |>
      summarise(
        w = median(mean_width, na.rm = TRUE),
        .by = c(proportions_c, n_cells_c, method, arm, ci_method)
      ) |>
      pivot_wider(names_from = c(arm), values_from = w) |>
      arrange(method, ci_method, proportions_c, n_cells_c)
  ),
  digits = 3
)
cat("\n== bootstrap scenarios ==\n")
boot_ids <- unique(iv$ID[iv$ci_method %in% c("bootstrap", "bootstrap_nprior")])
print(
  as.data.frame(
    iv |>
      filter(ID %in% boot_ids, method == "Newton-Raphson") |>
      select(
        ID,
        arm,
        ci_method,
        coverage,
        coverage_lower,
        coverage_upper,
        mean_width
      ) |>
      arrange(ID, arm, ci_method)
  ),
  digits = 3
)

cat("\n== boundary LRT ==\n")
l <- b$boundary_lrt |>
  left_join(cfg[c("ID", "proportions", "n_cells")], by = "ID")
print(
  as.data.frame(
    l |>
      filter(is.finite(p_value)) |>
      mutate(absent = p2_realised == 0) |>
      summarise(
        rate = mean(p_value < 0.05),
        n = n(),
        .by = c(proportions, n_cells, arm, absent)
      ) |>
      arrange(absent, proportions, n_cells, arm)
  ),
  digits = 3
)

cat("\n== Newton / lsei RMSE ratio (oracle), median by |rho1 - rho2| ==\n")
d <- rm |>
  filter(arm == "oracle", method %in% c("Newton-Raphson", "lsei")) |>
  select(ID, proportions, n_cells, corner, method, rmse) |>
  pivot_wider(names_from = method, values_from = rmse) |>
  mutate(
    gap = round(
      abs(
        as.numeric(sub(" .*", "", corner)) - as.numeric(sub(".* ", "", corner))
      ),
      1
    )
  )
print(
  as.data.frame(
    d |>
      summarise(
        ratio = median(`Newton-Raphson` / lsei),
        .by = c(proportions, n_cells, gap)
      ) |>
      pivot_wider(names_from = gap, values_from = ratio)
  ),
  digits = 3
)
print(
  as.data.frame(
    d |>
      summarise(
        better = mean(`Newton-Raphson` < 0.95 * lsei),
        worse = mean(`Newton-Raphson` > 1.05 * lsei),
        .by = n_cells
      )
  ),
  digits = 3
)

cat("\n== nominal / oracle RMSE ratio ==\n")
o <- rm |>
  filter(
    method %in% c("Newton-Raphson", "lsei", "DeCovarT diagonal"),
    arm != "none"
  ) |>
  select(ID, proportions, n_cells, corner, method, arm, rmse) |>
  pivot_wider(names_from = arm, values_from = rmse)
print(
  as.data.frame(
    o |>
      summarise(
        ratio = median(nominal / oracle),
        q90 = quantile(nominal / oracle, 0.9),
        .by = c(method, proportions, n_cells)
      ) |>
      arrange(method, proportions, n_cells)
  ),
  digits = 3
)
print(
  as.data.frame(
    o |>
      filter(method == "Newton-Raphson", corner %in% c("0 0", "-0.8 0.8")) |>
      select(proportions, n_cells, corner, oracle, nominal)
  ),
  digits = 3
)

cat("\n== mean estimate, SB90_Ba_n1e3 ==\n")
print(
  as.data.frame(
    opt |>
      filter(ID == "SB90_Ba_n1e3") |>
      summarise(
        estimate = mean(celltype_1),
        target = mean(p1_true),
        .by = algorithm
      )
  ),
  digits = 3
)
