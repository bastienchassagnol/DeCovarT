.single_cell_deriv_setup <- function() {
  withr::with_seed(7L, {
    mean_signature_matrix <- matrix(
      c(4, 1, 0.4, 0.2, 0.3, 3, 0.5, 2, 0.2, 0.4, 5, 1),
      nrow = 4,
      ncol = 3
    )
    sigma <- array(0, dim = c(4, 4, 3))
    for (j in seq_len(3L)) {
      draw <- matrix(stats::rnorm(16), 4, 4)
      sigma[,, j] <- crossprod(draw) / 4 + diag(0.4, 4)
    }
    p <- c(0.2, 0.5, 0.3)
    y <- 15 * drop(mean_signature_matrix %*% p) + stats::rnorm(4)
    list(
      mean_signature_matrix = mean_signature_matrix,
      Sigma = sigma,
      p = p,
      z = isometric_log_ratio(p),
      y = y,
      mean_scale = 15,
      variance_scale = 6
    )
  })
}

test_that("single-cell score and Hessian match numDeriv", {
  setup <- .single_cell_deriv_setup()
  args <- list(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma,
    mean_scale = setup$mean_scale,
    variance_scale = setup$variance_scale
  )
  gradient_numerical <- do.call(
    numDeriv::grad,
    c(
      list(
        func = loglik_single_cell,
        x = setup$p,
        method = "Richardson",
        method.args = list(eps = 1e-4, r = 6)
      ),
      args
    )
  )
  expect_equal(
    gradient_loglik_single_cell(
      setup$p,
      setup$y,
      setup$mean_signature_matrix,
      setup$Sigma,
      setup$mean_scale,
      setup$variance_scale
    ),
    gradient_numerical,
    tolerance = 1e-4
  )

  hessian_numerical <- do.call(
    numDeriv::hessian,
    c(
      list(
        func = loglik_single_cell,
        x = setup$p,
        method = "Richardson",
        method.args = list(eps = 1e-4, r = 6)
      ),
      args
    )
  )
  hessian_theoretical <- hessian_loglik_single_cell(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma,
    setup$mean_scale,
    setup$variance_scale
  )
  expect_true(isSymmetric(hessian_theoretical, tol = 1e-8))
  expect_equal(hessian_theoretical, hessian_numerical, tolerance = 1e-4)
})

test_that("ILR score and Hessian of the single-cell model match numDeriv", {
  setup <- .single_cell_deriv_setup()
  args <- list(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma,
    mean_scale = setup$mean_scale,
    variance_scale = setup$variance_scale
  )
  gradient_numerical <- do.call(
    numDeriv::grad,
    c(
      list(
        func = loglik_single_cell_ilr,
        x = setup$z,
        method = "Richardson",
        method.args = list(eps = 1e-4, r = 6)
      ),
      args
    )
  )
  expect_equal(
    gradient_loglik_single_cell_ilr(
      setup$z,
      setup$y,
      setup$mean_signature_matrix,
      setup$Sigma,
      setup$mean_scale,
      setup$variance_scale
    ),
    as.numeric(gradient_numerical),
    tolerance = 1e-4
  )
  hessian_numerical <- do.call(
    numDeriv::hessian,
    c(
      list(
        func = loglik_single_cell_ilr,
        x = setup$z,
        method = "Richardson",
        method.args = list(eps = 1e-4, r = 6)
      ),
      args
    )
  )
  hessian_theoretical <- hessian_loglik_single_cell_ilr(
    setup$z,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma,
    setup$mean_scale,
    setup$variance_scale
  )
  expect_equal(hessian_theoretical, hessian_numerical, tolerance = 1e-4)
})

test_that("observed Hessian stays symmetric for four cell types", {
  withr::local_seed(3)
  n_genes <- 6L
  n_types <- 4L
  mu <- matrix(stats::rexp(n_genes * n_types, 0.2), n_genes, n_types)
  sigma <- array(0, dim = c(n_genes, n_genes, n_types))
  for (j in seq_len(n_types)) {
    draw <- matrix(stats::rnorm(n_genes^2), n_genes)
    sigma[,, j] <- crossprod(draw) / n_genes + diag(0.3, n_genes)
  }
  p <- c(0.4, 0.3, 0.2, 0.1)
  y <- 10 * drop(mu %*% p) + stats::rnorm(n_genes)
  hessian <- hessian_loglik_single_cell(p, y, mu, sigma, 10, 4)
  numerical <- numDeriv::hessian(
    loglik_single_cell,
    p,
    y = y,
    mean_signature_matrix = mu,
    Sigma = sigma,
    mean_scale = 10,
    variance_scale = 4
  )
  expect_true(isSymmetric(hessian, tol = 1e-8))
  expect_equal(hessian, numerical, tolerance = 1e-5)
})

test_that("stage i recovers a noise-free composition and its GLS scale", {
  mu <- matrix(
    c(4, 0.2, 0.3, 1, 0.2, 3, 0.4, 2),
    nrow = 4
  )
  sigma <- array(rep(diag(0.2, 4), 2), dim = c(4, 4, 2))
  p <- c(0.7, 0.3)
  y <- 25 * drop(mu %*% p)
  y <- y + c(0.4, -0.2, 0.1, 0.3)
  stage <- single_cell_nuisance_scales(y, mu, sigma)
  expect_equal(stage$p, p, tolerance = 0.05)
  expect_gt(stage$mean_scale, 0)
  expect_gt(stage$variance_scale, 0)
  expect_equal(
    stage$n_cells,
    stage$mean_scale^2 / stage$variance_scale
  )
  expect_equal(stage$platform_scale, stage$variance_scale / stage$mean_scale)
})

.single_cell_grn <- function() {
  withr::with_seed(20261022L, {
    simulate_hierarchical_grn_moments(
      n_genes = 10L,
      n_celltypes = 2L,
      mean_scale = 10,
      target_cosine = 0,
      nonnegative = TRUE,
      precision_shift = 0.5,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  })
}

test_that("cellular fits recover p and c on the four scenarios", {
  moments <- .single_cell_grn()
  design <- single_cell_scenario_design()
  n_samples <- 4L
  for (scenario in design$scenarios) {
    ref <- rescale_single_cell_moments(
      moments,
      scenario$cell_size,
      reference = "cellular_raw"
    )
    for (scale in design$platform_scale) {
      sim <- simulate_bulk_from_cells(
        ref$generative_mean,
        ref$generative_sigma,
        p = scenario$proportions,
        n_cells = design$n_cells,
        platform_scale = scale,
        n_samples = n_samples,
        seed = 1000L + as.integer(10 * scale)
      )
      p_error <- numeric(n_samples)
      c_error <- numeric(n_samples)
      phi_hat <- numeric(n_samples)
      for (i in seq_len(n_samples)) {
        fit <- deconvolute_single_cell(
          sim$Y[, i],
          ref$mean_signature_matrix,
          ref$Sigma,
          estimand = "cellular"
        )
        expect_identical(fit$estimand, "cellular")
        expect_identical(fit$convergence, 0L)
        p_error[[i]] <- max(abs(fit$p - sim$p_realised[, i]))
        c_error[[i]] <- abs(fit$mean_scale - sim$mean_scale[[i]]) /
          sim$mean_scale[[i]]
        phi_hat[[i]] <- fit$variance_scale
      }
      expect_lt(max(p_error), 0.05)
      expect_lt(max(c_error), 0.10)
      # phi is a scale with G = 10 degrees of freedom, so one sample
      # has relative sd sqrt(2/G). The median absolute relative error
      # should stay inside that sampling range.
      phi_true <- sim$variance_scale[[1L]]
      phi_relative <- abs(phi_hat - phi_true) / phi_true
      expect_lt(median(phi_relative), 0.75)
    }
  }
})

test_that("Mahalanobis scale is unbiased when the composition is known", {
  moments <- .single_cell_grn()
  ref <- rescale_single_cell_moments(
    moments,
    c(1, 5),
    reference = "cellular_raw"
  )
  sim <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = c(0.8, 0.2),
    n_cells = 10000,
    platform_scale = 2,
    n_samples = 30L,
    seed = 19
  )
  phi_hat <- numeric(30L)
  c_hat <- numeric(30L)
  for (i in seq_len(30L)) {
    stage <- single_cell_nuisance_scales(
      sim$Y[, i],
      ref$mean_signature_matrix,
      ref$Sigma,
      p = sim$p_realised[, i]
    )
    phi_hat[[i]] <- stage$variance_scale
    c_hat[[i]] <- stage$mean_scale
  }
  expect_lt(
    abs(mean(c_hat) - sim$mean_scale[[1L]]) / sim$mean_scale[[1L]],
    0.02
  )
  expect_lt(
    abs(mean(phi_hat) - sim$variance_scale[[1L]]) /
      sim$variance_scale[[1L]],
    0.2
  )
})

test_that("a library-normalised reference targets transcript fractions", {
  moments <- .single_cell_grn()
  cell_size <- c(1, 5)
  generative <- rescale_single_cell_moments(
    moments,
    cell_size,
    reference = "cellular_raw"
  )
  transcript <- rescale_single_cell_moments(
    moments,
    cell_size,
    reference = "transcript_raw"
  )
  sim <- simulate_bulk_from_cells(
    generative$generative_mean,
    generative$generative_sigma,
    p = c(0.5, 0.5),
    n_cells = 10000,
    platform_scale = 1,
    n_samples = 4L,
    seed = 17
  )
  q_error <- numeric(4L)
  for (i in seq_len(4L)) {
    fit <- deconvolute_single_cell(
      sim$Y[, i],
      transcript$mean_signature_matrix,
      transcript$Sigma,
      estimand = "transcript",
      cell_size = cell_size
    )
    expect_identical(fit$convergence, 0L)
    q_error[[i]] <- max(abs(fit$q - sim$q[, i]))
    size_total <- sum(sim$p_realised[, i] * cell_size)
    expect_equal(
      fit$mean_scale,
      sim$mean_scale[[i]] * size_total,
      tolerance = 0.05
    )
  }
  expect_lt(max(q_error), 0.03)
})

test_that("profile score and Hessian match numDeriv in p and ILR", {
  # G - J residual degrees of freedom must be several, or phi_hat is tiny
  # and finite differences lose accuracy on the profile.
  setup <- withr::with_seed(8L, {
    n_genes <- 8L
    mu <- matrix(stats::rexp(n_genes * 3L, 0.3), n_genes, 3L)
    sigma <- array(0, dim = c(n_genes, n_genes, 3L))
    for (j in seq_len(3L)) {
      draw <- matrix(stats::rnorm(n_genes^2), n_genes)
      sigma[,, j] <- crossprod(draw) / n_genes + diag(0.4, n_genes)
    }
    p <- c(0.2, 0.5, 0.3)
    list(
      mean_signature_matrix = mu,
      Sigma = sigma,
      p = p,
      z = isometric_log_ratio(p),
      y = 15 * drop(mu %*% p) + 4 * stats::rnorm(n_genes)
    )
  })
  profile_p <- function(p) {
    scales <- .single_cell_scales_at_p(
      setup$y,
      setup$mean_signature_matrix,
      setup$Sigma,
      p
    )
    loglik_single_cell(
      p,
      setup$y,
      setup$mean_signature_matrix,
      setup$Sigma,
      scales$mean_scale,
      scales$variance_scale
    )
  }
  scales <- .single_cell_scales_at_p(
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma,
    setup$p
  )
  analytic <- .profile_hessian_single_cell(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma,
    scales$mean_scale,
    scales$variance_scale
  )
  expect_equal(
    analytic,
    numDeriv::hessian(profile_p, setup$p),
    tolerance = 1e-5
  )
  args <- list(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  expect_equal(
    do.call(gradient_profile_loglik_single_cell_ilr, c(list(setup$z), args)),
    numDeriv::grad(
      function(z) do.call(profile_loglik_single_cell_ilr, c(list(z), args)),
      setup$z
    ),
    tolerance = 1e-6
  )
  expect_equal(
    do.call(hessian_profile_loglik_single_cell_ilr, c(list(setup$z), args)),
    numDeriv::hessian(
      function(z) do.call(profile_loglik_single_cell_ilr, c(list(z), args)),
      setup$z
    ),
    tolerance = 1e-5
  )
})

test_that("profile Hessian matches numDeriv for four cell types", {
  withr::local_seed(5)
  n_genes <- 7L
  n_types <- 4L
  mu <- matrix(stats::rexp(n_genes * n_types, 0.2), n_genes, n_types)
  sigma <- array(0, dim = c(n_genes, n_genes, n_types))
  for (j in seq_len(n_types)) {
    draw <- matrix(stats::rnorm(n_genes^2), n_genes)
    sigma[,, j] <- crossprod(draw) / n_genes + diag(0.3, n_genes)
  }
  z <- isometric_log_ratio(c(0.4, 0.3, 0.2, 0.1))
  y <- 10 * drop(mu %*% c(0.4, 0.3, 0.2, 0.1)) + stats::rnorm(n_genes)
  analytic <- hessian_profile_loglik_single_cell_ilr(z, y, mu, sigma)
  expect_true(isSymmetric(analytic, tol = 1e-8))
  expect_equal(
    analytic,
    numDeriv::hessian(
      profile_loglik_single_cell_ilr,
      z,
      y = y,
      mean_signature_matrix = mu,
      Sigma = sigma
    ),
    tolerance = 1e-5
  )
})

test_that("profile Newton reaches the joint maximum on the four scenarios", {
  moments <- .single_cell_grn()
  design <- single_cell_scenario_design()
  for (scenario in design$scenarios) {
    ref <- rescale_single_cell_moments(moments, scenario$cell_size)
    for (scale in design$platform_scale) {
      sim <- simulate_bulk_from_cells(
        ref$generative_mean,
        ref$generative_sigma,
        p = scenario$proportions,
        n_cells = design$n_cells,
        platform_scale = scale,
        n_samples = 3L,
        seed = 2000L + as.integer(10 * scale)
      )
      for (i in seq_len(3L)) {
        y <- sim$Y[, i]
        fit <- deconvolute_single_cell(
          y,
          ref$mean_signature_matrix,
          ref$Sigma
        )
        expect_identical(fit$fit_type, "profile")
        expect_identical(fit$convergence, 0L)
        # Independent reference: BFGS on the closed-form
        # profile, started at the fit, must not find a better point.
        reference <- stats::optim(
          fit$rho,
          function(z) {
            -profile_loglik_single_cell_ilr(
              z,
              y,
              ref$mean_signature_matrix,
              ref$Sigma
            )
          },
          method = "BFGS",
          control = list(reltol = 1e-14, maxit = 500L)
        )
        expect_lt(
          max(abs(isometric_logistic(reference$par) - fit$p)),
          1e-6
        )
        expect_lt(-reference$value - fit$loglik, 1e-8)
      }
    }
  }
})

test_that("two-stage and fixed-scale fits are labelled and differ", {
  moments <- .single_cell_grn()
  ref <- rescale_single_cell_moments(moments, c(1, 5))
  sim <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = c(0.5, 0.5),
    n_cells = 10000,
    platform_scale = 1,
    n_samples = 1L,
    seed = 11
  )
  y <- sim$Y[, 1]
  joint <- deconvolute_single_cell(y, ref$mean_signature_matrix, ref$Sigma)
  two_stage <- deconvolute_single_cell(
    y,
    ref$mean_signature_matrix,
    ref$Sigma,
    profile = FALSE
  )
  fixed <- deconvolute_single_cell(
    y,
    ref$mean_signature_matrix,
    ref$Sigma,
    mean_scale = sim$mean_scale[[1L]],
    variance_scale = sim$variance_scale[[1L]]
  )
  expect_identical(two_stage$fit_type, "two_stage")
  expect_identical(fixed$fit_type, "fixed_scales")
  expect_identical(fixed$mean_scale, sim$mean_scale[[1L]])
  expect_gte(joint$loglik, two_stage$loglik - 1e-10)
})

test_that("the fit needs more genes than cell types", {
  mu <- diag(3)
  sigma <- array(rep(diag(3), 3), dim = c(3, 3, 3))
  expect_error(
    deconvolute_single_cell(c(1, 2, 3), mu, sigma),
    "G > J"
  )
})

test_that("both single-cell optimisers are equivariant to the bulk scale", {
  moments <- withr::with_seed(
    20261022L,
    simulate_hierarchical_grn_moments(
      n_genes = 12L,
      n_celltypes = 3L,
      mean_scale = 10,
      target_cosine = 0.3,
      nonnegative = TRUE,
      precision_shift = 0.5,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  )
  ref <- rescale_single_cell_moments(moments, cell_size = c(1, 1, 5))
  sim <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = c(0.5, 0.3, 0.2),
    n_cells = 1000,
    platform_scale = 5,
    seed = 3L
  )
  y <- sim$Y[, 1]
  for (optimiser in c("newton", "lbfgsb")) {
    fit <- deconvolute_single_cell(
      y,
      ref$mean_signature_matrix,
      ref$Sigma,
      optimiser = optimiser,
      epsilon = 1e-10
    )
    scaled <- deconvolute_single_cell(
      7 * y,
      ref$mean_signature_matrix,
      ref$Sigma,
      optimiser = optimiser,
      epsilon = 1e-10
    )
    expect_equal(scaled$p, fit$p, tolerance = 1e-5)
    expect_equal(scaled$mean_scale, 7 * fit$mean_scale, tolerance = 1e-5)
    expect_equal(
      scaled$variance_scale,
      49 * fit$variance_scale,
      tolerance = 1e-5
    )
  }
  newton <- deconvolute_single_cell(y, ref$mean_signature_matrix, ref$Sigma)
  lbfgsb <- deconvolute_single_cell(
    y,
    ref$mean_signature_matrix,
    ref$Sigma,
    optimiser = "lbfgsb",
    epsilon = 1e-12
  )
  expect_equal(lbfgsb$p, newton$p, tolerance = 1e-4)
})
