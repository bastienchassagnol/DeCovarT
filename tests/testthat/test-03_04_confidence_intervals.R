# ILR back-transformation, Firth and the single-cell Schur intervals.

test_that("two cell types map the Wald endpoints of rho", {
  mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
  colnames(mu) <- c("ct1", "ct2")
  sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
  y <- drop(mu %*% c(0.6, 0.4))

  fit <- confint_decovart(y, mu, sigma)
  expect_identical(fit$method, "ilr")
  expect_identical(fit$mapping, "endpoint")
  expect_true(all(fit$interval > 0))
  expect_true(all(fit$interval < 1))
  expect_true(all(fit$interval[, 1] <= fit$estimate))
  expect_true(all(fit$interval[, 2] >= fit$estimate))

  manual <- rbind(
    isometric_logistic(fit$rho_interval[[1]]),
    isometric_logistic(fit$rho_interval[[2]])
  )
  expect_equal(
    fit$interval[, 1],
    pmin(manual[1, ], manual[2, ]),
    ignore_attr = TRUE
  )
})

test_that("replicate columns divide the ILR covariance", {
  setup <- .inference_setup()
  p <- c(ct1 = 0.5, ct2 = 0.3, ct3 = 0.2)
  one <- .population_chart_covariance(
    p,
    setup$signature,
    setup$Sigma,
    n_obs = 1L
  )
  two <- .population_chart_covariance(
    p,
    setup$signature,
    setup$Sigma,
    n_obs = 2L
  )
  expect_equal(two$vcov_p, one$vcov_p / 2)
  expect_lt(max(abs(rowSums(one$vcov_p))), 1e-8)
})

test_that("three cell types use the Monte Carlo ILR map", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.5, 0.3, 0.2))
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  fit <- confint_decovart(
    y,
    setup$signature,
    setup$Sigma,
    n_draws = 80L,
    seed = 11
  )
  restored <- if (is.null(old)) {
    !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), old)
  }
  expect_true(restored)
  expect_identical(fit$mapping, "monte_carlo")
  expect_true(all(fit$interval > 0))
  expect_true(all(fit$interval < 1))
  expect_identical(fit$n_draws, 80L)
})

test_that("Firth raises the penalised population objective", {
  mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
  colnames(mu) <- c("ct1", "ct2")
  sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
  y <- as.matrix(drop(mu %*% c(0.6, 0.4)))
  start <- restricted_mle_decovart(y, mu, sigma)$coefficients
  start <- repair_simplex(start, open = TRUE, nms = colnames(mu))
  base <- .pooled_population_loglik(start, y, mu, sigma) +
    .jeffreys_penalty_population(isometric_log_ratio(start), mu, sigma)
  fit <- fit_firth_decovart(y, mu, sigma, itmax = 40L)
  expect_gte(fit$penalised_loglik + 1e-4, base)
  expect_equal(sum(fit$coefficients), 1)
  expect_true(all(fit$coefficients > 0))
  wald <- confint_decovart(y, mu, sigma, method = "firth", itmax = 40L)
  expect_identical(wald$method, "firth")
  expect_identical(dim(wald$interval_ilr), c(2L, 2L))
})

test_that("the reference dispatcher resamples cells when donors are absent", {
  mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
  colnames(mu) <- c("ct1", "ct2")
  sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
  refs <- withr::with_seed(
    1L,
    lapply(seq_len(2), function(j) {
      draws <- MASS::mvrnorm(n = 4, mu = mu[, j], Sigma = sigma[,, j])
      t(draws)
    })
  )
  names(refs) <- colnames(mu)
  y <- drop(mu %*% c(0.6, 0.4))
  fit <- confint_decovart(
    y,
    method = "reference",
    reference_profiles = refs,
    n_boot = 3L,
    seed = 2
  )
  expect_identical(fit$resample, "cells")
  expect_true(fit$regenerate_bulk)
  expect_identical(dim(fit$interval), c(2L, 2L))

  donor_ids <- lapply(refs, function(x) rep(c("a", "b"), length.out = ncol(x)))
  donors <- confint_decovart(
    y,
    method = "reference",
    reference_profiles = refs,
    donor_ids = donor_ids,
    n_boot = 2L,
    seed = 3
  )
  expect_identical(donors$resample, "donors")
  expect_false(donors$regenerate_bulk)
})

test_that("single-cell Schur intervals stay on the simplex and on the scales", {
  setup <- .single_cell_interval_setup()
  fit <- confint_single_cell(setup$y, setup$mu, setup$sigma)
  expect_identical(fit$method, "ilr")
  expect_identical(fit$mapping, "endpoint")
  expect_true(all(fit$interval > 0))
  expect_true(all(fit$interval < 1))
  expect_true(all(fit$interval[, 1] <= fit$estimate))
  expect_true(all(fit$estimate <= fit$interval[, 2]))
  expect_lt(max(abs(rowSums(fit$vcov))), 1e-6)
  expect_identical(dim(fit$nuisance$vcov), c(2L, 2L))
  expect_true(all(fit$nuisance$interval[, 1] < fit$nuisance$estimate))
  expect_true(all(fit$nuisance$interval[, 2] > fit$nuisance$estimate))
  expect_true(all(fit$nuisance$interval > 0))

  direct <- deconvolute_single_cell(setup$y, setup$mu, setup$sigma)
  expect_equal(direct$vcov, fit$vcov)
  expect_identical(dim(direct$vcov_nuisance), c(2L, 2L))
})

test_that("single-cell profile and boundary tests are calibrated by stratum", {
  setup <- .single_cell_interval_setup()
  profile <- confint_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    method = "profile"
  )
  expect_true(all(profile$interval[, 1] <= profile$estimate + 1e-6))
  expect_true(all(profile$interval[, 2] >= profile$estimate - 1e-6))
  expect_true(all(profile$interval >= 0))
  expect_true(all(profile$interval <= 1))

  interior <- lrt_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    null_value = c(ct2 = 0.3)
  )
  far <- lrt_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    null_value = c(ct2 = 0.05)
  )
  expect_identical(interior$calibration, "chi-square (Wilks)")
  expect_lt(interior$statistic, far$statistic)

  boundary <- lrt_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    null_value = c(ct2 = 0)
  )
  expect_identical(boundary$calibration, "chi-bar-square (Self-Liang)")
  expect_gte(boundary$statistic, 0)
  expect_true(boundary$p_value > 0)
  expect_true(boundary$p_value <= 1)
})

test_that("single-cell Firth and bootstrap intervals have the right shape", {
  setup <- .single_cell_interval_setup()
  start <- deconvolute_single_cell(setup$y, setup$mu, setup$sigma)
  base <- start$loglik +
    .jeffreys_penalty_single_cell(start$rho, setup$y, setup$mu, setup$sigma)
  fit <- fit_firth_single_cell(setup$y, setup$mu, setup$sigma, itmax = 30L)
  expect_gte(fit$penalised_loglik + 1e-3, base)
  expect_equal(sum(fit$coefficients), 1)

  boot <- bootstrap_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    n_boot = 3L,
    itmax = 25L
  )
  expect_identical(dim(boot$interval), c(2L, 2L))
  expect_identical(dim(boot$scale_interval), c(4L, 2L))
  expect_true(all(boot$scale_interval > 0, na.rm = TRUE))
  expect_equal(
    boot$scales["platform_scale_a", ] * boot$scales["cell_count_n", ],
    boot$scales["mean_scale", ]
  )

  boot_lbfgsb <- bootstrap_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    n_boot = 3L,
    itmax = 25L,
    optimiser = "lbfgsb"
  )
  expect_identical(dim(boot_lbfgsb$interval), c(2L, 2L))
  expect_true(all(boot_lbfgsb$interval >= 0 & boot_lbfgsb$interval <= 1))
})

test_that("single-cell reference bootstrap regenerates without donors", {
  setup <- .single_cell_interval_setup()
  refs <- withr::with_seed(
    4L,
    lapply(seq_len(2), function(j) {
      draws <- MASS::mvrnorm(
        n = 6,
        mu = setup$mu[, j],
        Sigma = setup$sigma[,, j]
      )
      t(draws)
    })
  )
  names(refs) <- colnames(setup$mu)
  fit <- reference_bootstrap_single_cell(
    setup$y,
    refs,
    n_boot = 2L,
    itmax = 20L
  )
  expect_identical(fit$resample, "cells")
  expect_true(fit$regenerate_bulk)
  expect_identical(dim(fit$interval), c(2L, 2L))
})

# Numerical checks ---------------------------------------------------------

.pooled_population_setup <- function() {
  withr::with_seed(11L, {
    moments <- simulate_hierarchical_grn_moments(
      n_genes = 10L,
      n_celltypes = 2L,
      target_cosine = 0.3,
      precision_shift = 0.3,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
    y <- simulate_bulk_mixture(
      moments$mean_profiles,
      moments$covariance_matrices,
      p = c(0.7, 0.3),
      n = 5L
    )$Y
  })
  list(
    y = y,
    mu = moments$mean_profiles,
    sigma = moments$covariance_matrices
  )
}

test_that("the pooled restricted fit reaches the maximum, not a vertex", {
  # Regression: BFGS from the barycentre took one step to a vertex, where
  # the likelihood is flat, and stopped there.
  setup <- .pooled_population_setup()
  fit <- restricted_mle_decovart(setup$y, setup$mu, setup$sigma)
  pooled <- function(a) {
    .pooled_population_loglik(c(a, 1 - a), setup$y, setup$mu, setup$sigma)
  }
  best <- stats::optimize(pooled, c(1e-6, 1 - 1e-6), maximum = TRUE)
  expect_equal(unname(fit$coefficients[[1L]]), best$maximum, tolerance = 1e-4)
  expect_gte(fit$loglik, best$objective - 1e-8)
  expect_named(fit$convergence, c("code", "iterations", "message", "rdm"))
})

test_that("profile endpoints sit where the deviance meets the cut-off", {
  setup <- .pooled_population_setup()
  interval <- confint_profile_decovart(setup$y, setup$mu, setup$sigma)
  full <- restricted_mle_decovart(setup$y, setup$mu, setup$sigma)
  cutoff <- stats::qchisq(0.95, df = 1L)
  for (bound in c("lower", "upper")) {
    value <- interval["celltype_1", bound]
    restricted <- restricted_mle_decovart(
      setup$y,
      setup$mu,
      setup$sigma,
      fixed = c(celltype_1 = value)
    )
    expect_equal(
      2 * (full$loglik - restricted$loglik),
      cutoff,
      tolerance = 1e-5
    )
  }
})

test_that("profile roots return the face or NA, never the estimate", {
  expect_identical(
    .profile_root(function(x) -1, estimate = 0.4, face = 1),
    1
  )
  expect_identical(
    .profile_root(function(x) -1, estimate = 0.4, face = 0),
    0
  )
  expect_warning(
    root <- .profile_root(function(x) 1, estimate = 0.4, face = 1),
    "did not reach the maximum"
  )
  expect_true(is.na(root))
  expect_equal(
    .profile_root(function(x) x - 0.7, estimate = 0.4, face = 1),
    0.7,
    tolerance = 1e-6
  )
})

test_that("a null at a vertex counts J - 1 boundary constraints", {
  expect_identical(
    .null_constraint_counts(c(ct1 = 1), n_celltypes = 3L, tol = 1e-6),
    list(n_constraints = 2L, n_boundary = 2L, df_interior = 0L)
  )
  expect_identical(
    .null_constraint_counts(c(ct3 = 0), n_celltypes = 3L, tol = 1e-6),
    list(n_constraints = 1L, n_boundary = 1L, df_interior = 0L)
  )
  expect_identical(
    .null_constraint_counts(
      c(ct1 = 0.5, ct2 = 0.5),
      n_celltypes = 3L,
      tol = 1e-6
    ),
    list(n_constraints = 2L, n_boundary = 1L, df_interior = 1L)
  )
  expect_identical(
    .null_constraint_counts(c(ct1 = 0.4), n_celltypes = 3L, tol = 1e-6),
    list(n_constraints = 1L, n_boundary = 0L, df_interior = 1L)
  )
})

.single_cell_three_types <- function() {
  moments <- withr::with_seed(31L, {
    simulate_hierarchical_grn_moments(
      n_genes = 10L,
      n_celltypes = 3L,
      mean_scale = 10,
      target_cosine = 0,
      nonnegative = TRUE,
      precision_shift = 0.5,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  })
  ref <- rescale_single_cell_moments(moments, c(1, 2, 3))
  sim <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = c(0.5, 0.3, 0.2),
    n_cells = 1e4,
    platform_scale = 1,
    n_samples = 1L,
    seed = 32L
  )
  list(y = sim$Y[, 1L], mu = ref$mean_signature_matrix, sigma = ref$Sigma)
}

test_that("single-cell covariances invert the joint observed information", {
  setup <- .single_cell_three_types()
  fit <- deconvolute_single_cell(setup$y, setup$mu, setup$sigma)
  theta <- c(fit$rho, log(fit$mean_scale), log(fit$variance_scale))
  full <- function(theta) {
    loglik_single_cell(
      isometric_logistic(theta[1:2]),
      setup$y,
      setup$mu,
      setup$sigma,
      exp(theta[[3L]]),
      exp(theta[[4L]])
    )
  }
  joint <- solve(-numDeriv::hessian(full, theta))
  scale_jac <- diag(c(fit$mean_scale, fit$variance_scale))
  expect_equal(unname(fit$vcov_rho), joint[1:2, 1:2], tolerance = 1e-5)
  expect_equal(
    unname(fit$vcov_nuisance),
    scale_jac %*% joint[3:4, 3:4] %*% scale_jac,
    tolerance = 1e-5
  )
})

test_that("fixed scales use the fixed-scale information", {
  setup <- .single_cell_three_types()
  profile <- deconvolute_single_cell(setup$y, setup$mu, setup$sigma)
  fixed <- deconvolute_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    mean_scale = profile$mean_scale,
    variance_scale = profile$variance_scale
  )
  h_fixed <- hessian_loglik_single_cell_ilr(
    fixed$rho,
    setup$y,
    setup$mu,
    setup$sigma,
    profile$mean_scale,
    profile$variance_scale
  )
  expect_equal(fixed$vcov_rho, solve(-.symmetrise(h_fixed)), tolerance = 1e-8)
  expect_true(all(is.na(fixed$vcov_nuisance)))
  # Estimating the scales can only widen the composition variance.
  expect_true(all(diag(profile$vcov_rho) >= diag(fixed$vcov_rho) - 1e-12))
})

test_that("bootstraps report failed refits", {
  setup <- .single_cell_interval_setup()
  boot <- withr::with_seed(
    3L,
    bootstrap_single_cell(setup$y, setup$mu, setup$sigma, n_boot = 3L)
  )
  expect_identical(boot$n_failed, sum(is.na(boot$estimates[1L, ])))
})

test_that("the REML calibration widens single-cell intervals by t and G / df", {
  setup <- .single_cell_three_types()
  ml <- confint_single_cell(setup$y, setup$mu, setup$sigma, method = "wald")
  reml <- confint_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    method = "wald",
    variance = "reml"
  )
  n_genes <- length(setup$y)
  df <- n_genes - 3L
  expect_identical(reml$df, df)
  expect_equal(reml$vcov, ml$vcov * n_genes / df)
  half_ml <- (ml$interval[, 2L] - ml$interval[, 1L]) / 2
  half_reml <- (reml$interval[, 2L] - reml$interval[, 1L]) / 2
  expect_equal(
    half_reml / half_ml,
    rep(sqrt(n_genes / df) * stats::qt(0.975, df) / stats::qnorm(0.975), 3L),
    ignore_attr = TRUE
  )
  phi <- reml$nuisance$estimate[["variance_scale"]]
  expect_equal(
    unname(reml$nuisance$interval["variance_scale", ]),
    n_genes * phi / stats::qchisq(c(0.975, 0.025), df)
  )
})

test_that("the REML profile cut-off is the exact linear-model one", {
  setup <- .single_cell_interval_setup()
  ml <- confint_profile_single_cell(setup$y, setup$mu, setup$sigma)
  reml <- confint_profile_single_cell(
    setup$y,
    setup$mu,
    setup$sigma,
    variance = "reml"
  )
  expect_true(all(reml[, "lower"] <= ml[, "lower"] + 1e-8))
  expect_true(all(reml[, "upper"] >= ml[, "upper"] - 1e-8))
})

test_that("platform scale and cell count intervals are positive and consistent", {
  setup <- .single_cell_interval_setup()
  fit <- confint_single_cell(setup$y, setup$mu, setup$sigma, method = "wald")
  nuisance <- fit$nuisance
  expect_identical(
    rownames(nuisance$interval),
    c("mean_scale", "variance_scale", "platform_scale_a", "cell_count_n")
  )
  expect_true(all(nuisance$interval > 0))
  expect_true(all(nuisance$interval[, 1] < nuisance$estimate))
  expect_true(all(nuisance$interval[, 2] > nuisance$estimate))
  expect_equal(
    nuisance$estimate[["platform_scale_a"]] *
      nuisance$estimate[["cell_count_n"]],
    nuisance$estimate[["mean_scale"]]
  )
})
