# G3.0: never assert exact equality of floating-point values.
.tol_srr <- 100 * .Machine$double.eps

skip_if_not_extended <- function() {
  testthat::skip_if(
    Sys.getenv("DECOVART_EXTENDED_TESTS") != "true",
    "Set DECOVART_EXTENDED_TESTS=true to run extended tests."
  )
}

.toy_deconvolution <- function() {
  path <- system.file(
    "extdata",
    "toy_deconvolution.rds",
    package = "DeCovarT"
  )
  if (!nzchar(path) || !file.exists(path)) {
    path <- testthat::test_path("fixtures", "toy_deconvolution.rds")
  }
  readRDS(path)
}
