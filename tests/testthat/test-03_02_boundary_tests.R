# Boundary tests: restricted MLE on simplex faces, chi-bar-square
# calibration, likelihood-ratio tests and boundary diagnostics.
#
#' @srrstats {RE4.3} Interior and boundary LRTs (`lrt_decovart()`) are
#'   exercised below.
#' @srrstats {G5.8} Zero-proportion (active-face) nulls use the
#'   chi-bar-square mixture rather than a Wald / Wilks \(\chi^2_1\).

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
    y = replicates,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma,
    null_value = c(ct3 = 0)
  )
  expect_identical(true_null$n_boundary, 1L)
  expect_identical(true_null$calibration, "chi-bar-square (Self-Liang)")
  expect_gte(true_null$statistic, 0)
  expect_gt(true_null$p_value, 0.05)

  false_null <- lrt_decovart(
    y = replicates,
    mean_signature_matrix = setup$signature,
    Sigma = setup$Sigma,
    null_value = c(ct1 = 0.05)
  )
  expect_identical(false_null$n_boundary, 0L)
  expect_identical(false_null$calibration, "chi-square (Wilks)")
  expect_lt(false_null$p_value, 0.01)
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
  expect_lt(interior$min_eigenvalue, 0)
})
