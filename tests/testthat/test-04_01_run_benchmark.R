# Benchmark driver: deconvolute_ratios(), compute_benchmark_metrics() and
# the per-sample diagnostics.

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

  expect_named(
    scores,
    c("regression", "monte_carlo", "optimisation", "intervals")
  )
  expect_named(scores$regression, c("global", "cell_type"))
  expect_equal(scores$regression$global$tv, 0.05)
  expect_equal(scores$regression$global$rmse, 0.05)
  expect_equal(scores$regression$global$maxae, 0.05)
  expect_true(is.na(scores$regression$global$reconstitution_mae))
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

  expect_named(
    scores,
    c("regression", "monte_carlo", "optimisation", "intervals")
  )
  expect_true(is.na(scores$regression$global$tv))
  expect_true(is.finite(scores$regression$global$rmse))
  expect_identical(nrow(scores$regression$cell_type), 0L)
  expect_identical(nrow(scores$monte_carlo), 0L)
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
      "theoretical_se",
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

  expect_named(
    out,
    c("regression", "monte_carlo", "optimisation", "intervals")
  )
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
  expect_named(
    out,
    c("regression", "monte_carlo", "optimisation", "intervals")
  )
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
  expect_warning(
    deconvolute_ratios(
      signature_matrix = toy$signature_matrix,
      bulk_expression = bulk_neg,
      deconvolution_functions = fns,
      cores = 1
    ),
    "negative"
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

test_that("Newton-Raphson attaches ILR Wald SEs; nnls does not", {
  skip_if_not_installed("nnls")
  toy <- .toy_deconvolution()
  out <- withr::with_seed(
    11L,
    deconvolute_ratios(
      signature_matrix = toy$signature_matrix,
      bulk_expression = toy$bulk_expression,
      true_ratios = toy$true_ratios,
      Sigma = toy$Sigma,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls),
        "Newton-Raphson" = list(
          FUN = deconvolute_ratios_Newton_Raphson,
          additional_parameters = list(itmax = 40L, epsilon = 1e-4)
        )
      ),
      cores = 1L
    )
  )
  mc <- out$monte_carlo
  newton <- mc[mc$algorithm == "Newton-Raphson", , drop = FALSE]
  nnls <- mc[mc$algorithm == "nnls", , drop = FALSE]
  expect_true(all(is.finite(newton$mean_model_se)))
  expect_true(all(is.finite(newton$se_sd_ratio)))
  expect_true(all(is.finite(newton$coverage)))
  expect_true(all(is.na(nnls$mean_model_se)))
  expect_true(all(is.na(nnls$coverage)))
})

test_that("solvers store one convergence record, closed forms store NA", {
  mu <- matrix(
    c(20, 40, 15, 40, 20, 25),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("ct", 1:2))
  )
  sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
  y <- drop(mu %*% c(0.6, 0.4))
  newton <- deconvolute_ratios_Newton_Raphson(
    y,
    mu,
    sigma,
    return_model = TRUE
  )$convergence
  marquardt <- deconvolute_ratios_Marquardt_Levenberg(
    y,
    mu,
    sigma,
    return_model = TRUE
  )$convergence
  expect_named(newton, c("code", "iterations", "message", "rdm"))
  expect_named(marquardt, c("code", "iterations", "message", "rdm"))
  expect_type(newton$iterations, "integer")
  expect_true(is.na(newton$rdm))
  expect_true(is.finite(marquardt$rdm))
  out <- suppressWarnings(deconvolute_ratios(
    signature_matrix = mu,
    bulk_expression = cbind(y, y),
    Sigma = sigma,
    deconvolution_functions = list(
      nnls = list(FUN = deconvolute_ratios_nnls),
      newton = list(FUN = deconvolute_ratios_Newton_Raphson)
    ),
    cores = 1L
  ))
  conv <- out$optimisation$converged
  algorithm <- out$optimisation$algorithm
  expect_true(all(is.na(vapply(
    conv[algorithm == "nnls"],
    function(x) x$iterations,
    integer(1)
  ))))
  expect_true(all(
    vapply(conv[algorithm == "newton"], function(x) x$code, integer(1)) == 0L
  ))
})

test_that("optimisation table stores Hessian extrema and converged lists", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  scores <- compute_benchmark_metrics(
    y = drop(mu %*% c(0.4, 0.6)),
    mean_signature_matrix = mu,
    estimated_p = c(0.45, 0.55),
    true_ratios = c(0.4, 0.6)
  )
  opt <- scores$optimisation
  expect_true("converged" %in% names(opt))
  expect_true("local_maximum" %in% names(opt))
  expect_true("min_eigenvalue" %in% names(opt))
  expect_type(opt$converged, "list")
  expect_true(is.na(opt$converged[[1L]]$iterations))
})

test_that(".process_memory_bytes is a finite NA or a single byte count", {
  mem <- DeCovarT:::.process_memory_bytes()
  expect_type(mem, "double")
  expect_length(mem, 1L)
  expect_true(is.na(mem) || is.finite(mem))
})
