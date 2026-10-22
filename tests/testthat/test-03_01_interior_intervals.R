# Interior intervals: expected Fisher information (Wald) and
# profile-likelihood intervals.
#
#' @srrstats {RE4.3} Profile intervals (`confint_profile_decovart()`) are
#'   exercised below; Wald intervals in test-02_03_population_decovart_fit.R.

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

test_that("restricted fits on p_j = c peak at the unrestricted MLE", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.5, 0.3, 0.2))

  profile <- vapply(
    c(0.05, 0.2, 0.5),
    function(value) {
      restricted_mle_decovart(
        y,
        setup$signature,
        setup$Sigma,
        fixed = c(ct3 = value)
      )$loglik
    },
    numeric(1)
  )

  expect_length(profile, 3L)
  expect_identical(which.max(profile), 2L)
  expect_error(
    restricted_mle_decovart(y, setup$signature, setup$Sigma, fixed = c(0.2)),
    "named by columns"
  )
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
    y = replicates,
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
