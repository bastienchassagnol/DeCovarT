# Shared scenario (3 genes, 3 cell types) for the derivative checks below.
# A three-component simplex yields a non-trivial (2 x 2) constrained Hessian.
.decovart_deriv_setup <- function() {
  withr::with_seed(7L, {
    mean_signature_matrix <- matrix(
      c(20, 40, 15, 40, 20, 25, 25, 30, 35),
      nrow = 3,
      dimnames = list(paste0("gene_", 1:3), paste0("celltype_", 1:3))
    )
    Sigma <- array(0, dim = c(3, 3, 3))
    Sigma[,, 1] <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3)
    Sigma[,, 2] <- matrix(c(2, -0.2, 0, -0.2, 1.5, 0.3, 0, 0.3, 1), 3, 3)
    Sigma[,, 3] <- matrix(
      c(1.2, 0.1, 0.4, 0.1, 0.8, -0.1, 0.4, -0.1, 1.1),
      3,
      3
    )
    p <- c(0.2, 0.3, 0.5)
    y <- mean_signature_matrix %*% p + rnorm(3)
    list(
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      p = p,
      rho = additive_log_ratio(p),
      y = y
    )
  })
}

# Richardson extrapolation is the most stable finite-difference scheme for
# validating the closed-form derivatives (see ?numDeriv::grad).
.richardson_args <- list(eps = 1e-4, r = 6)


test_that("Unconstrained gradient and Hessian of the objective match numDeriv", {
  setup <- .decovart_deriv_setup()

  gradient_numerical <- numDeriv::grad(
    loglik_multivariate,
    setup$p,
    method = "Richardson",
    method.args = .richardson_args,
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  gradient_theoretical <- gradient_loglik_unconstrained(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    gradient_numerical,
    gradient_theoretical,
    tolerance = 1e-4
  )

  hessian_numerical <- numDeriv::hessian(
    loglik_multivariate,
    setup$p,
    method = "Richardson",
    method.args = .richardson_args,
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  hessian_theoretical <- hessian_loglik_unconstrained(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    hessian_numerical,
    hessian_theoretical,
    tolerance = 1e-4
  )
})


test_that("Constrained gradient and Hessian of the objective match numDeriv", {
  setup <- .decovart_deriv_setup()

  gradient_numerical <- numDeriv::grad(
    loglik_multivariate_constrained,
    setup$rho,
    method = "Richardson",
    method.args = .richardson_args,
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  gradient_theoretical <- as.numeric(gradient_loglik_constrained(
    setup$rho,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  ))
  testthat::expect_equal(
    gradient_numerical,
    gradient_theoretical,
    tolerance = 1e-4
  )

  hessian_numerical <- numDeriv::hessian(
    loglik_multivariate_constrained,
    setup$rho,
    method = "Richardson",
    method.args = .richardson_args,
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  hessian_theoretical <- hessian_loglik_constrained(
    setup$rho,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    hessian_numerical,
    hessian_theoretical,
    tolerance = 1e-4
  )
})


test_that("Additive logistic Jacobian and Hessian match numDeriv", {
  setup <- .decovart_deriv_setup()

  jacobian_numerical <- numDeriv::jacobian(additive_logistic, setup$rho)
  jacobian_theoretical <- jacobian_additive_logistic(setup$rho)
  testthat::expect_equal(
    jacobian_numerical,
    jacobian_theoretical,
    tolerance = 1e-4
  )

  # Hessian tensor is stored slice-wise as [ , , i] = d^2 p_i / (d rho d rho^T)
  hessian_theoretical <- hessian_additive_logistic(setup$rho)
  for (i in seq_along(setup$p)) {
    hessian_numerical <- numDeriv::hessian(
      function(rho) additive_logistic(rho)[i],
      setup$rho
    )
    testthat::expect_equal(
      hessian_numerical,
      hessian_theoretical[,, i],
      tolerance = 1e-4
    )
  }
})


test_that("Every deconvolute_ratios_* solver returns a valid simplex on the challenging 3-component scenario", {
  setup <- .decovart_deriv_setup()

  # All frequentist DeCovarT solvers exported from
  # R/03_03_DeCovarT_estimate_ratios_frequentist.R.
  deconvolution_solvers <- list(
    Marquardt_Levenberg = deconvolute_ratios_Marquardt_Levenberg,
    simulated_annealing = deconvolute_ratios_simulated_annealing,
    L_BFGS_B = deconvolute_ratios_L_BFGS_B,
    Newton_Raphson = deconvolute_ratios_Newton_Raphson,
    gradient_descent = deconvolute_ratios_gradient_descent
  )
  celltype_names <- colnames(setup$mean_signature_matrix)

  for (solver_name in names(deconvolution_solvers)) {
    inferred <- withr::with_seed(
      3L,
      suppressWarnings(
        deconvolution_solvers[[solver_name]](
          y = setup$y,
          mean_signature_matrix = setup$mean_signature_matrix,
          Sigma = setup$Sigma,
          true_ratios = setup$p,
          epsilon = 10^-4,
          itmax = 200
        )
      )
    )

    estimated_ratios <- unlist(inferred[celltype_names], use.names = FALSE)

    testthat::expect_length(estimated_ratios, length(setup$p))
    testthat::expect_false(anyNA(estimated_ratios), info = solver_name)
    testthat::expect_true(
      all(estimated_ratios >= 0 & estimated_ratios <= 1),
      info = solver_name
    )
    testthat::expect_equal(sum(estimated_ratios), 1, tolerance = 1e-6)
  }
})

# test_that("Benchmark standard deconvolution algorithms against DeCovarT", {
#   deconvolution_functions <- list(
#     "lm" = list(FUN = deconvolute_ratios_abbas),
#     "nnls" = list(FUN = deconvolute_ratios_nnls),
#     "lsei" = list(FUN = deconvolute_ratios_deconRNASeq),
#     "LBFGS" = list(
#       FUN = deconvolute_ratios_L_BFGS_B,
#       additional_parameters = list(epsilon = 10^-3, itmax = 200)
#     ),
#     "gradient" = list(
#       FUN = deconvolute_ratios_first_order,
#       additional_parameters = list(epsilon = 10^-3, itmax = 200)
#     ),
#     "Newton-Raphson" = list(
#       FUN = deconvolute_ratios_second_order,
#       additional_parameters = list(epsilon = 10^-3, itmax = 200)
#     ),
#     "Marquardt-Levenberg" = list(
#       FUN = deconvolute_ratios_Marquardt_Levenberg,
#       additional_parameters = list(epsilon = 10^-3, itmax = 200)
#     ),
#     "SA" = list(
#       FUN = deconvolute_ratios_simulated_annealing,
#       additional_parameters = list(epsilon = 10^-3, itmax = 200)
#     )
#   )

#   bivariate_scenario <- withr::with_seed(
#     seed = 3L,
#     benchmark_bivariate_gaussian_convolutions(
#       proportions = list("balanced" = c(0.50, 0.50)),
#       n = 2,
#       corr_sequence = c(-0.75, 0.75),
#       signature_matrices = list(
#         "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
#       )
#     ),
#     .rng_kind = "L'Ecuyer-CMRG"
#   )
#   bivariate_configuraton <- readRDS(testthat::test_path(
#     "fixtures",
#     "bivariate_configuraton.rds"
#   ))
#   bivariate_estimation <- readRDS(testthat::test_path(
#     "fixtures",
#     "bivariate_estimation.rds"
#   ))

#   expect_equal(bivariate_configuraton, bivariate_scenario$config)
#   expect_equal(bivariate_estimation, bivariate_scenario$simulations)
# })
