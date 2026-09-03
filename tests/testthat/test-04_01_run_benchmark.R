test_that("compute_benchmark_metrics scores against true ratios", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  true_p <- c(0.4, 0.6)
  y <- drop(mu %*% true_p)
  estimated_p <- c(0.45, 0.55)

  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = true_p
  )

  expect_named(scores, c("regression", "monte_carlo", "optimisation"))
  expect_named(scores$regression, c("global", "cell_type"))
  expect_equal(scores$regression$global$tv, .tv(true_p, estimated_p))
  expect_equal(scores$regression$global$rmse, .rmse(true_p, estimated_p))
  expect_equal(scores$regression$global$sdid, .sdid(true_p, estimated_p))
  expect_equal(scores$regression$global$maxae, .maxae(true_p, estimated_p))
  expect_identical(nrow(scores$regression$cell_type), 2L)
  expect_identical(nrow(scores$monte_carlo), 2L)
  expect_gt(scores$regression$global$rmse, 0)
})

test_that("compute_benchmark_metrics scores reconstituted bulk without truth", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  estimated_p <- c(0.5, 0.5)
  y <- drop(mu %*% c(0.4, 0.6))

  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = NULL
  )

  expect_named(scores, c("regression", "monte_carlo", "optimisation"))
  expect_true(is.na(scores$regression$global$tv))
  expect_true(is.finite(scores$regression$global$rmse))
  expect_identical(nrow(scores$regression$cell_type), 0L)
  expect_identical(nrow(scores$monte_carlo), 0L)
})

test_that("simplex metrics attain the vertex-swap bounds", {
  p <- c(1, 0)
  p_hat <- c(0, 1)
  expect_equal(.tv(p, p_hat), 1)
  expect_equal(.rmse(p, p_hat), 1)
  expect_equal(.maxae(p, p_hat), 1)
  expect_equal(.angular_distance(p, p_hat), 1)
  expect_equal(.sdid(p, p_hat), 1)
  expect_equal(.mae(p, p_hat), 1)
  expect_equal(.tv(p, p), 0)
})

test_that("KKT residual vanishes for a constant ambient score", {
  p <- c(0.2, 0.3, 0.5)
  grad <- c(4, 4, 4)
  expect_equal(.kkt_residual(p, grad), 0, tolerance = 1e-10)
})

test_that("cell-type Pearson is computed across samples, not within a sample", {
  mu <- matrix(
    c(20, 22, 24, 22, 20, 18),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:3))
  )
  true_p <- rbind(
    c(0.5, 0.3, 0.1),
    c(0.3, 0.5, 0.2),
    c(0.2, 0.2, 0.7)
  )
  p_hat <- true_p
  y <- mu %*% true_p
  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = p_hat,
    true_ratios = true_p
  )
  expect_equal(scores$regression$cell_type$pearson, c(1, 1, 1))
  expect_equal(unname(scores$regression$global$tv), c(0, 0, 0))
})

test_that("Monte Carlo block reports ADEMP columns separately from global scores", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  true_p <- cbind(c(0.4, 0.6), c(0.4, 0.6), c(0.4, 0.6))
  p_hat <- cbind(c(0.42, 0.58), c(0.38, 0.62), c(0.41, 0.59))
  se <- matrix(0.02, nrow = 2, ncol = 3)
  y <- mu %*% true_p
  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = p_hat,
    true_ratios = true_p,
    se = se
  )
  expect_named(
    scores$monte_carlo,
    c(
      "algorithm",
      "cell_type",
      "bias",
      "empirical_sd",
      "mean_model_sd",
      "mean_model_se",
      "se_sd_ratio",
      "rmse",
      "coverage",
      "coverage_lower",
      "coverage_upper",
      "coverage_interval",
      "mean_interval_width",
      "mcse_coverage"
    )
  )
  expect_named(
    scores$regression$global,
    c(
      "sample_id",
      "algorithm",
      "tv",
      "rmse",
      "angular",
      "sdid",
      "maxae",
      "reconstitution_mae",
      "reconstitution_cor"
    )
  )
  expect_named(
    scores$regression$cell_type,
    c(
      "algorithm",
      "cell_type",
      "pearson",
      "presence_f1",
      "false_positive_mass"
    )
  )
  expect_equal(unname(scores$monte_carlo$bias), c(0.01, -0.01) / 3)
  expect_true(all(scores$monte_carlo$coverage == 1))
})

test_that("presence F1 and false-positive mass detect spillover onto a null type", {
  mu <- matrix(
    c(20, 22, 24, 22, 20, 18),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:3))
  )
  true_p <- cbind(c(0.5, 0.5, 0), c(0.6, 0.4, 0), c(0.7, 0.3, 0))
  p_hat <- cbind(c(0.45, 0.45, 0.1), c(0.55, 0.35, 0.1), c(0.65, 0.25, 0.1))
  y <- mu %*% true_p
  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = p_hat,
    true_ratios = true_p
  )
  ct3 <- scores$regression$cell_type$cell_type == "ct3"
  expect_equal(scores$regression$cell_type$false_positive_mass[ct3], 0.1)
  expect_equal(scores$regression$cell_type$presence_f1[ct3], 0)
})

test_that("Monte Carlo coverage MCSE matches the binomial formula", {
  expect_equal(.mcse_coverage(c(1, 1, 1, 0)), sqrt(0.75 * 0.25 / 4))
})

test_that("coverage_mc_interval defaults to Wilson and can use Wald", {
  covered <- c(TRUE, TRUE, TRUE, FALSE)
  wilson <- coverage_mc_interval(covered, method = "wilson")
  wald <- coverage_mc_interval(covered, method = "wald")
  ac <- coverage_mc_interval(covered, method = "agresti_coull")
  expect_equal(wilson$coverage, 0.75)
  expect_equal(wilson$mcse, sqrt(0.75 * 0.25 / 4))
  expect_equal(wald$lower, max(0, 0.75 - stats::qnorm(0.975) * wald$mcse))
  expect_gt(wilson$lower, 0)
  expect_lt(wilson$upper, 1)
  expect_true(ac$lower >= 0 && ac$upper <= 1)
  expect_identical(wilson$method, "wilson")
})

test_that("deconvolute_ratios returns a three-block metric list for nnls", {
  skip_if_not_installed("nnls")
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(genes, cts)
  )
  bulk <- cbind(
    drop(mu %*% c(0.5, 0.5)),
    drop(mu %*% c(0.3, 0.7))
  )
  dimnames(bulk) <- list(genes, paste0("sample_", 1:2))

  out <- deconvolute_ratios(
    signature_matrix = mu,
    bulk_expression = bulk,
    true_ratios = cbind(c(0.5, 0.5), c(0.3, 0.7)),
    deconvolution_functions = list(
      "nnls" = list(FUN = deconvolute_ratios_nnls)
    ),
    cores = 1
  )

  expect_named(out, c("regression", "monte_carlo", "optimisation"))
  expect_identical(nrow(out$optimisation), 2L)
  expect_true(all(c("ct1", "ct2") %in% names(out$optimisation)))
  expect_true("elapsed_sec" %in% names(out$optimisation))
  expect_true("kkt_residual" %in% names(out$optimisation))
  expect_identical(nrow(out$regression$global), 2L)
  expect_identical(nrow(out$regression$cell_type), 2L)
})

test_that("deconvolute_ratios errors on missing values (G2.13)", {
  skip_if_not_installed("nnls")
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(genes, cts)
  )
  bulk <- cbind(drop(mu %*% c(0.5, 0.5)))
  dimnames(bulk) <- list(genes, "sample_1")
  fns <- list("nnls" = list(FUN = deconvolute_ratios_nnls))

  bulk_na <- bulk
  bulk_na[1, 1] <- NA_real_
  expect_error(
    deconvolute_ratios(
      signature_matrix = mu,
      bulk_expression = bulk_na,
      deconvolution_functions = fns,
      cores = 1
    ),
    "missing"
  )

  mu_na <- mu
  mu_na[1, 1] <- NA_real_
  expect_error(
    deconvolute_ratios(
      signature_matrix = mu_na,
      bulk_expression = bulk,
      deconvolution_functions = fns,
      cores = 1
    ),
    "missing"
  )
})

test_that(".match_arg_case_insensitive is case-insensitive (G2.3b)", {
  expect_identical(
    .match_arg_case_insensitive("SIGMA", c("either", "sigma", "Theta")),
    "sigma"
  )
  expect_identical(
    .match_arg_case_insensitive("theta", c("either", "sigma", "Theta")),
    "Theta"
  )
  expect_error(
    .match_arg_case_insensitive("not_a_choice", c("either", "sigma", "Theta")),
    "must be one of"
  )
})

test_that("deconvolute_ratios uses the bundled toy fixture (G5.1)", {
  skip_if_not_installed("nnls")
  toy <- .toy_deconvolution()
  out <- deconvolute_ratios(
    signature_matrix = toy$signature_matrix,
    bulk_expression = toy$bulk_expression,
    true_ratios = toy$true_ratios,
    Sigma = toy$Sigma,
    deconvolution_functions = list(
      "nnls" = list(FUN = deconvolute_ratios_nnls)
    ),
    cores = 1
  )
  expect_named(out, c("regression", "monte_carlo", "optimisation"))
  expect_identical(nrow(out$optimisation), ncol(toy$bulk_expression))
})

test_that("deconvolute_ratios errors on Inf, negatives, and J > G", {
  skip_if_not_installed("nnls")
  toy <- .toy_deconvolution()
  fns <- list("nnls" = list(FUN = deconvolute_ratios_nnls))

  bulk_inf <- toy$bulk_expression
  bulk_inf[1, 1] <- Inf
  expect_error(
    deconvolute_ratios(
      signature_matrix = toy$signature_matrix,
      bulk_expression = bulk_inf,
      deconvolution_functions = fns,
      cores = 1
    ),
    "Inf"
  )

  bulk_neg <- toy$bulk_expression
  bulk_neg[1, 1] <- -1
  expect_error(
    deconvolute_ratios(
      signature_matrix = toy$signature_matrix,
      bulk_expression = bulk_neg,
      deconvolution_functions = fns,
      cores = 1
    ),
    "non-negative"
  )

  expect_error(
    deconvolute_ratios(
      signature_matrix = as.data.frame(toy$signature_matrix),
      bulk_expression = toy$bulk_expression,
      deconvolution_functions = fns,
      cores = 1
    ),
    "data.frame"
  )

  wide_mu <- cbind(toy$signature_matrix, extra = c(21, 19))
  expect_error(
    deconvolute_ratios(
      signature_matrix = wide_mu,
      bulk_expression = toy$bulk_expression,
      deconvolution_functions = fns,
      cores = 1
    ),
    "Undetermined"
  )
})

test_that("deconvolute_ratios rejects log2 mixing (RE2.3)", {
  skip_if_not_installed("nnls")
  toy <- .toy_deconvolution()
  expect_error(
    deconvolute_ratios(
      signature_matrix = toy$signature_matrix,
      bulk_expression = toy$bulk_expression,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls)
      ),
      scaled = TRUE,
      cores = 1
    ),
    "log2"
  )
})

test_that(".ensure_file_suffix adds or rejects extensions (G4.0)", {
  expect_identical(.ensure_file_suffix("scores", "rds"), "scores.rds")
  expect_identical(.ensure_file_suffix("scores.RDS", "rds"), "scores.RDS")
  expect_error(.ensure_file_suffix("scores.csv", "rds"), "must use suffix")
})

test_that(".write_artefact writes RDS under with_tempfile (G4.0)", {
  toy <- .toy_deconvolution()
  withr::with_tempfile("tf", fileext = ".rds", {
    path <- .write_artefact(toy, tf, kind = "rds")
    expect_identical(tools::file_ext(path), "rds")
    expect_identical(readRDS(path)$true_ratios, toy$true_ratios)
  })
})
