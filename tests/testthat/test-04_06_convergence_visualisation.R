.toy_theta_2d <- function() {
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  list(p = c(ct1 = 0.5, ct2 = 0.5), mu = mu, sigma = Sigma)
}

test_that("whitened ILR is zero at the true composition", {
  th <- .toy_theta_2d()
  info_z <- DeCovarT:::.ilr_expected_information(th$p, th$mu, th$sigma)
  w <- DeCovarT:::.whiten_ilr_estimate(th$p, th$p, info_z)
  expect_type(w, "list")
  expect_equal(w$whitened, 0, tolerance = 1e-10)
  expect_equal(w$mahalanobis, 0, tolerance = 1e-10)
})

test_that("plot_mc_qq_normal and plot_mc_chi2_qq build ggplots", {
  skip_if_not_installed("ggplot2")
  df <- tibble::tibble(
    algorithm = rep(c("LBFGS", "Marquardt-Levenberg"), each = 40L),
    panel = "rho=(0,0)",
    sample_id = rep(paste0("s", seq_len(40L)), times = 2L),
    ilr_coord = "ILR1",
    whitened = stats::rnorm(80L)
  )
  df$mahalanobis <- df$whitened^2
  df$chi_df <- 1L
  p_qq <- plot_mc_qq_normal(df, title = "test")
  expect_s3_class(p_qq, "ggplot")
  p_chi <- plot_mc_chi2_qq(
    DeCovarT:::.mahalanobis_from_whitened(df),
    title = "test"
  )
  expect_s3_class(p_chi, "ggplot")
})

test_that("plot_ks_normality_box returns a grob with an inset legend", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  df <- tibble::tibble(
    algorithm = rep(c("LBFGS", "Newton-Raphson"), each = 40L),
    panel = "rho=(0,0)",
    whitened = c(stats::rnorm(40L), stats::rnorm(40L, mean = 0.8))
  )
  p <- plot_ks_normality_box(df, title = "test")
  expect_true(inherits(p, "gg") || inherits(p, "gtable"))
})

test_that("expected log-likelihood at p* is below the mean-as-y loglik", {
  th <- .toy_theta_2d()
  y_mean <- drop(th$mu %*% th$p)
  ll_mean <- loglik_multivariate(th$p, y_mean, th$mu, th$sigma)
  ell_exp <- DeCovarT:::.expected_loglik_multivariate(
    th$p,
    th$p,
    th$mu,
    th$sigma
  )
  g <- nrow(th$mu)
  expect_equal(ell_exp, ll_mean - g / 2, tolerance = 1e-8)
  ell_far <- DeCovarT:::.expected_loglik_multivariate(
    c(0.9, 0.1),
    th$p,
    th$mu,
    th$sigma
  )
  expect_lt(ell_far, ell_exp)
})

test_that("plot_convergence_stacked builds a legend grob", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  df <- tibble::tibble(
    algorithm = rep(c("LBFGS", "ols"), each = 20L),
    panel = "rho=(0,0)",
    numerical_converged = c(rep(TRUE, 18L), rep(FALSE, 22L)),
    theoretical_converged = c(rep(TRUE, 15L), rep(FALSE, 25L))
  )
  p <- plot_convergence_stacked(
    df,
    title = "test",
    icon_dir = file.path("inst", "extdata", "convergence_icons")
  )
  expect_true(
    inherits(p, "gg") || inherits(p, "gtable") || inherits(p, "grob")
  )
})
