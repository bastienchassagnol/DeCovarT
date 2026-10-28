# nocov start
# Package load / attach hooks -----------------------------------------------

#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  attached <- paste0("package:", search())
  lines <- c(
    cli::format_inline("DeCovarT {ver}: covariance-aware bulk deconvolution."),
    cli::format_inline(
      "Main entry points: {.fn fit_decovart} and {.fn deconvolute_ratios}."
    )
  )
  if ("package:MASS" %in% attached && "package:dplyr" %in% attached) {
    lines <- c(
      lines,
      cli::format_inline(
        "Both {.pkg MASS} and {.pkg dplyr} are attached; ",
        "prefer {.fn dplyr::select}."
      )
    )
  }
  packageStartupMessage(paste(lines, collapse = "\n"))
}
# nocov end
