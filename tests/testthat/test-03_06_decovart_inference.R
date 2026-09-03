.inference_setup <- function() {
  signature <- matrix(
    c(20, 40, 15, 40, 20, 25, 25, 30, 35),
    nrow = 3,
    dimnames = list(paste0("gene_", 1:3), paste0("ct", 1:3))
  )
  Sigma <- array(
    0,
    dim = c(3, 3, 3),
    dimnames = list(
      rownames(signature),
      rownames(signature),
      colnames(signature)
    )
  )
  for (j in 1:3) {
    Sigma[,, j] <- diag(3) * c(1, 1.5, 1.2)[j]
  }
  list(signature = signature, Sigma = Sigma)
}


test_that("loglik_multivariate is the Gaussian log-density up to a constant", {
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  y <- drop(setup$signature %*% p) + c(0.4, -0.2, 0.1)

  covariance <- .compute_global_variance(p, setup$Sigma)
  residual <- y - drop(setup$signature %*% p)
  reference <- -0.5 *
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) -
    0.5 * drop(t(residual) %*% solve(covariance) %*% residual)

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma),
    reference,
    tolerance = 1e-10
  )
})


test_that("expected Fisher information equals the mean of -Hessian", {
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  information <- expected_fisher_unconstrained(p, setup$signature, setup$Sigma)

  draws <- withr::with_seed(
    42L,
    MASS::mvrnorm(
      n = 4000L,
      mu = drop(setup$signature %*% p),
      Sigma = .compute_global_variance(p, setup$Sigma)
    )
  )
  monte_carlo <- Reduce(
    `+`,
    lapply(
      seq_len(nrow(draws)),
      function(b) {
        -hessian_loglik_unconstrained(
          p,
          draws[b, ],
          setup$signature,
          setup$Sigma
        )
      }
    )
  ) /
    nrow(draws)

  expect_equal(monte_carlo, information, tolerance = 0.02)
})


test_that("additive_logistic stays finite for extreme ALR coordinates", {
  p <- additive_logistic(c(800, 1200))

  expect_false(anyNA(p))
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_identical(which.max(p), 2L)
})


test_that("isometric_logistic stays finite for extreme ILR coordinates", {
  p <- isometric_logistic(c(800, 1200))

  expect_false(anyNA(p))
  expect_equal(sum(p), 1, tolerance = 1e-12)
})


test_that("restricted_mle_decovart honours the fixed coordinates", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.5, 0.5, 0))

  restricted <- restricted_mle_decovart(
    y,
    setup$signature,
    setup$Sigma,
    fixed = c(ct3 = 0)
  )

  expect_identical(restricted$coefficients[["ct3"]], 0)
  expect_equal(sum(restricted$coefficients), 1, tolerance = 1e-8)
  expect_lte(
    restricted$loglik,
    loglik_multivariate(
      restricted$coefficients,
      y,
      setup$signature,
      setup$Sigma
    ) +
      1e-8
  )
})


test_that("profile log-likelihood peaks at the unrestricted MLE", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.5, 0.3, 0.2))

  profile <- profile_loglik_decovart(
    y,
    setup$signature,
    setup$Sigma,
    celltype = "ct3",
    value = c(0.05, 0.2, 0.5)
  )

  expect_length(profile, 3L)
  expect_identical(unname(which.max(profile)), 2L)
})


test_that("chi_bar_square_pvalue matches the closed-form mixtures", {
  # One active constraint: 1/2 chi2_0 + 1/2 chi2_1.
  expect_equal(
    chi_bar_square_pvalue(2.7055, n_boundary = 1L),
    0.5 * stats::pchisq(2.7055, df = 1L, lower.tail = FALSE),
    tolerance = 1e-8
  )
  # chi2_0 is an atom, so the strict tail caps the p-value at 1 - w_0.
  expect_equal(chi_bar_square_pvalue(0, n_boundary = 1L), 0.5)
  # No boundary constraint recovers Wilks.
  expect_equal(
    chi_bar_square_pvalue(3.8415, n_boundary = 0L, df_interior = 1L),
    stats::pchisq(3.8415, df = 1L, lower.tail = FALSE),
    tolerance = 1e-8
  )
  expect_error(chi_bar_square_pvalue(1, 0L, 0L), "least one constraint")
  expect_error(chi_bar_square_pvalue(1, 1L, weights = 1), "length")
})


test_that("lrt_decovart separates a true from a false boundary null", {
  setup <- .inference_setup()
  replicates <- withr::with_seed(
    11L,
    simulate_bulk_mixture(
      signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      p = c(0.5, 0.5, 0),
      n = 6L
    )$Y
  )

  true_null <- lrt_decovart(
    bulk_expression = replicates,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma,
    null_value = c(ct3 = 0)
  )
  expect_identical(true_null$n_boundary, 1L)
  expect_identical(true_null$calibration, "chi-bar-square (Self-Liang)")
  expect_gte(true_null$statistic, 0)
  expect_gt(true_null$p_value, 0.05)

  false_null <- lrt_decovart(
    bulk_expression = replicates,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma,
    null_value = c(ct1 = 0.05)
  )
  expect_identical(false_null$n_boundary, 0L)
  expect_identical(false_null$calibration, "chi-square (Wilks)")
  expect_lt(false_null$p_value, 0.01)
})


test_that("profile intervals bracket the estimate and stay on the simplex", {
  setup <- .inference_setup()
  replicates <- withr::with_seed(
    7L,
    simulate_bulk_mixture(
      signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      p = c(0.45, 0.35, 0.20),
      n = 6L
    )$Y
  )

  interval <- confint_profile_decovart(
    bulk_expression = replicates,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma
  )

  expect_identical(dim(interval), c(3L, 3L))
  expect_named(interval[1, ], c("estimate", "lower", "upper"))
  expect_true(all(interval[, "lower"] <= interval[, "estimate"] + 1e-6))
  expect_true(all(interval[, "upper"] >= interval[, "estimate"] - 1e-6))
  expect_true(all(interval >= 0))
  expect_true(all(interval <= 1))
})


test_that("bootstrap_decovart reproduces the atom at zero under H0", {
  setup <- .inference_setup()
  replicates <- withr::with_seed(
    11L,
    simulate_bulk_mixture(
      signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      p = c(0.5, 0.5, 0),
      n = 4L
    )$Y
  )

  calibration <- withr::with_seed(
    5L,
    bootstrap_decovart(
      bulk_expression = replicates,
      mean_signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      null_value = c(ct3 = 0),
      n_boot = 60L
    )
  )

  expect_length(calibration$null_statistics, 60L)
  expect_true(all(calibration$null_statistics >= 0))
  # Chernoff / Self-Liang: about half the replicates land on the face.
  expect_gt(mean(calibration$null_statistics < 1e-8), 0.3)
  expect_lt(mean(calibration$null_statistics < 1e-8), 0.7)
  expect_gt(calibration$p_value, 0.05)
})


test_that("boundary_diagnostics flags faces and certifies interior optima", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.45, 0.35, 0.20))

  on_face <- boundary_diagnostics(
    c(0.5, 0.5, 1e-12),
    y,
    setup$signature,
    setup$Sigma
  )
  expect_true(on_face$near_boundary)
  expect_false(on_face$local_maximum)

  mle <- suppressWarnings(deconvolute_ratios_Marquardt_Levenberg(
    y = y,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma,
    itmax = 200
  ))
  interior <- boundary_diagnostics(
    as.numeric(mle),
    y,
    setup$signature,
    setup$Sigma
  )
  expect_false(interior$near_boundary)
  expect_lt(interior$max_eigenvalue, 0)
})


test_that("multistart_decovart explores several starts", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.45, 0.35, 0.20))

  restarts <- withr::with_seed(
    3L,
    multistart_decovart(
      y = y,
      mean_signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      n_starts = 3L,
      itmax = 100
    )
  )

  expect_identical(dim(restarts$starts), c(3L, 4L))
  expect_length(restarts$logliks, 4L)
  expect_type(restarts$multimodal, "logical")
  expect_equal(sum(restarts$coefficients), 1, tolerance = 1e-6)
})


test_that("fit_decovart stores boundary diagnostics and multi-start spread", {
  setup <- .inference_setup()
  replicates <- withr::with_seed(
    7L,
    simulate_bulk_mixture(
      signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      p = c(0.45, 0.35, 0.20),
      n = 2L
    )$Y
  )

  fit <- withr::with_seed(
    2L,
    fit_decovart(
      signature_matrix = setup$signature,
      bulk_expression = replicates,
      Sigma = setup$Sigma,
      itmax = 100,
      n_starts = 2L
    )
  )

  expect_length(fit$diagnostics, 2L)
  expect_named(
    fit$diagnostics[[1L]],
    c(
      "boundary_distance",
      "near_boundary",
      "score_norm",
      "max_eigenvalue",
      "local_maximum"
    )
  )
  expect_false(is.null(fit$convergence[[1L]]$multimodal))
})


test_that("loglik_multivariate matches dmvnorm up to the 2 pi constant", {
  skip_if_not_installed("mvtnorm")
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  y <- drop(setup$signature %*% p) + c(0.4, -0.2, 0.1)
  covariance <- .compute_global_variance(p, setup$Sigma)
  n_genes <- length(y)

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma) -
      0.5 * n_genes * log(2 * pi),
    as.numeric(mvtnorm::dmvnorm(
      y,
      mean = drop(setup$signature %*% p),
      sigma = covariance,
      log = TRUE
    )),
    tolerance = 1e-10
  )
})


test_that("Cholesky backsolve matches the explicit-inverse Mahalanobis term", {
  setup <- .inference_setup()
  p <- c(0.45, 0.35, 0.20)
  y <- drop(setup$signature %*% p) + c(0.1, -0.3, 0.2)
  residual <- as.numeric(y - drop(setup$signature %*% p))
  sigma_p <- .sigma_p_factorisation(p, setup$Sigma)
  z <- backsolve(sigma_p$chol, residual, transpose = TRUE)

  expect_equal(
    sum(z * z),
    .inner_product(residual, sigma_p$inverse),
    tolerance = 1e-10
  )
})


test_that("cell-type permutation leaves the Gaussian log-likelihood invariant", {
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  perm <- c(3L, 1L, 2L)
  y <- drop(setup$signature %*% p)
  mu_perm <- setup$signature[, perm, drop = FALSE]
  colnames(mu_perm) <- colnames(setup$signature)
  sigma_perm <- setup$Sigma[,, perm, drop = FALSE]
  if (!is.null(dimnames(setup$Sigma))) {
    dimnames(sigma_perm) <- dimnames(setup$Sigma)
  }

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma),
    loglik_multivariate(p[perm], y, mu_perm, sigma_perm),
    tolerance = 1e-10
  )

  gene_perm <- c(3L, 1L, 2L)
  y_gene <- y[gene_perm]
  expect_false(
    isTRUE(all.equal(
      loglik_multivariate(p, y, setup$signature, setup$Sigma),
      loglik_multivariate(p, y_gene, setup$signature, setup$Sigma),
      tolerance = 1e-8
    ))
  )
})


test_that("equivariance_check_decovart recovers the relabelled MLE", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.7, 0.2, 0.1))
  perm <- c(3L, 1L, 2L)

  check <- equivariance_check_decovart(
    y,
    setup$signature,
    setup$Sigma,
    perm = perm,
    itmax = 80
  )
  expect_equal(check$p_star, check$p_expected, tolerance = 1e-5)
  expect_lt(abs(check$loglik_diff), 1e-6)
  expect_identical(check$perm, perm)
})


test_that("reference_bootstrap_decovart resamples cells within type", {
  setup <- .inference_setup()
  p <- c(0.55, 0.30, 0.15)
  refs <- lapply(seq_len(3), function(j) {
    draws <- MASS::mvrnorm(
      n = 10L,
      mu = setup$signature[, j],
      Sigma = setup$Sigma[,, j]
    )
    out <- t(draws)
    rownames(out) <- rownames(setup$signature)
    out
  })
  names(refs) <- colnames(setup$signature)
  y <- drop(setup$signature %*% p)

  expect_error(
    reference_bootstrap_decovart(y, refs, n_boot = 2L, itmax = 40),
    "donor_ids"
  )

  boot <- withr::with_seed(
    13L,
    reference_bootstrap_decovart(
      y,
      refs,
      method = "cells",
      n_boot = 12L,
      itmax = 80
    )
  )
  expect_identical(boot$method, "cells")
  expect_identical(dim(boot$estimates), c(3L, 12L))
  expect_named(boot$p_hat, colnames(setup$signature))
  expect_equal(colSums(boot$estimates), rep(1, 12L), tolerance = 1e-6)
  expect_identical(dim(boot$interval), c(3L, 2L))
})


test_that("reference_bootstrap_decovart resamples donors within type", {
  setup <- .inference_setup()
  p <- c(0.55, 0.30, 0.15)
  refs <- lapply(seq_len(3), function(j) {
    draws <- MASS::mvrnorm(
      n = 12L,
      mu = setup$signature[, j],
      Sigma = setup$Sigma[,, j]
    )
    out <- t(draws)
    rownames(out) <- rownames(setup$signature)
    out
  })
  names(refs) <- colnames(setup$signature)
  donor_ids <- lapply(refs, function(x) {
    rep(paste0("d", 1:3), length.out = ncol(x))
  })
  y <- drop(setup$signature %*% p)

  boot <- withr::with_seed(
    14L,
    reference_bootstrap_decovart(
      y,
      refs,
      donor_ids = donor_ids,
      method = "donors",
      n_boot = 8L,
      itmax = 80
    )
  )
  expect_identical(boot$method, "donors")
  expect_identical(dim(boot$estimates), c(3L, 8L))
  expect_equal(colSums(boot$estimates), rep(1, 8L), tolerance = 1e-6)
})


test_that("reference_bootstrap_decovart draws Dirichlet compositions", {
  setup <- .inference_setup()
  refs <- lapply(seq_len(3), function(j) {
    draws <- MASS::mvrnorm(
      n = 8L,
      mu = setup$signature[, j],
      Sigma = setup$Sigma[,, j]
    )
    out <- t(draws)
    rownames(out) <- rownames(setup$signature)
    out
  })
  names(refs) <- colnames(setup$signature)
  y <- drop(setup$signature %*% c(0.5, 0.3, 0.2))

  boot <- withr::with_seed(
    15L,
    reference_bootstrap_decovart(
      y,
      refs,
      method = "dirichlet",
      n_boot = 8L,
      itmax = 80
    )
  )
  expect_identical(boot$method, "dirichlet")
  expect_identical(dim(boot$p_simulated), c(3L, 8L))
  expect_equal(colSums(boot$p_simulated), rep(1, 8L), tolerance = 1e-10)
  expect_equal(colSums(boot$estimates), rep(1, 8L), tolerance = 1e-6)
})
