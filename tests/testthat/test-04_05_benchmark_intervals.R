test_that("interval table treats a missing parameter row as no interval", {
  p_hat <- matrix(c(0.4, 0.6, 0.5, 0.5), 2, dimnames = list(c("a", "b"), NULL))
  full <- matrix(
    c(0.3, 0.5, 0.7, 0.9, 1, 2),
    3,
    dimnames = list(c("a", "b", "c"), c("lower", "upper"))
  )
  partial <- full["a", , drop = FALSE]
  intervals <- list(list(wald = full), list(wald = partial))
  out <- .interval_table(
    intervals,
    p_hat = p_hat,
    p_true = p_hat,
    nuisance_hat = NULL,
    nuisance_true = NULL
  )
  expect_identical(sort(out$parameter), c("a", "b"))
  expect_identical(out$n_intervals[out$parameter == "b"], 1L)
})
