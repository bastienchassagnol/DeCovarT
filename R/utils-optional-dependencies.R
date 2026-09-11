#' Check that an optional Suggests package is installed
#'
#' @param pkg Package name.
#' @param fn Calling function name (for the error message).
#' @keywords internal
#' @noRd
.check_suggested_package <- function(pkg, fn) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .ui_abort(c(
      "{.fn {fn}} requires the optional package {.pkg {pkg}}.",
      "i" = "Install it with install.packages(\"{pkg}\")."
    ))
  }
  invisible(TRUE)
}

#' Exported function from an optional package (no `pkg::` in source)
#'
#' `R CMD check` flags `pkg::fun` when `pkg` is not in Imports. Optional
#' helpers that must not pull Suggests onto CI use this instead.
#'
#' @keywords internal
#' @noRd
.suggested_export <- function(pkg, name) {
  getExportedValue(pkg, name)
}
