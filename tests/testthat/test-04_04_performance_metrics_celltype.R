# Cell-type performance metrics across Monte Carlo replicates.

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

test_that("presence F1 treats non-finite estimates as absent", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  true_p <- cbind(c(0.5, 0.5), c(0.6, 0.4), c(0.7, 0.3))
  p_hat <- cbind(c(0.45, 0.55), c(NA_real_, NA_real_), c(0.65, 0.35))
  y <- mu %*% true_p
  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = p_hat,
    true_ratios = true_p
  )
  expect_true(all(is.finite(scores$regression$cell_type$presence_f1)))
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
  expect_equal(
    coverage_mc_interval(c(1, 1, 1, 0), method = "wald")$mcse,
    sqrt(0.75 * 0.25 / 4)
  )
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
