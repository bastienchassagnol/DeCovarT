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
  expect_identical(out$n_numerical_failures, c(0L, 0L))
})

test_that("a non-finite estimate is a numerical failure, not a missing interval", {
  p_hat <- matrix(
    c(NA_real_, NA_real_, 0.4, 0.6),
    2,
    dimnames = list(c("a", "b"), NULL)
  )
  full <- matrix(
    c(0.2, 0.5, 0.6, 0.9),
    2,
    dimnames = list(c("a", "b"), c("lower", "upper"))
  )
  intervals <- list(NULL, list(wald = full))
  out <- .interval_table(
    intervals,
    p_hat = p_hat,
    p_true = matrix(
      c(0.4, 0.6, 0.4, 0.6),
      2,
      dimnames = list(c("a", "b"), NULL)
    )
  )
  row <- out[out$parameter == "a", ]
  expect_identical(row$n_numerical_failures, 1L)
  expect_identical(row$n_intervals, 1L)
})
