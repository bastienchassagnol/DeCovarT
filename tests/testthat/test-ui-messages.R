test_that("cli UI aborts keep test-facing substrings", {
  expect_error(
    .match_arg_case_insensitive("not_a_choice", c("either", "sigma")),
    "must be one of"
  )
  expect_error(
    .assert_no_missing(c(1, NA_real_), "bulk_expression"),
    "missing"
  )
  expect_error(
    .assert_no_missing(c(1, Inf), "bulk_expression"),
    "Inf"
  )
})

test_that("scenario labels prefer ID then metadata", {
  row <- tibble::tibble(
    ID = "B1_Ho",
    true_theta = list(list(p = 1)),
    proportion_name = "balanced"
  )
  expect_identical(.scenario_label(row), "B1_Ho")
  row$ID <- NA_character_
  expect_match(.scenario_label(row), "proportion_name=balanced")
})
