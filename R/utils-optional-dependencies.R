#' Check that an optional Suggests package is installed
#'
#' @param pkg Package name.
#' @param fn Calling function name (for the error message).
#' @keywords internal
#' @noRd
.check_suggested_package <- function(pkg, fn) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      fn,
      "() requires the optional package '",
      pkg,
      "'. Install with install.packages(\"",
      pkg,
      "\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
