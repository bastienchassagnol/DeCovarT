# Package load / attach hooks -----------------------------------------------

#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Reserved for package-level options (none required yet).
  invisible()
}

#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  attached <- paste0("package:", search())
  mass_dplyr <- "package:MASS" %in% attached && "package:dplyr" %in% attached
  if (.has_cli()) {
    lines <- c(
      "DeCovarT {ver}: covariance-aware bulk deconvolution.",
      "Main entry points: {.fn fit_decovart} and {.fn deconvolute_ratios}."
    )
    if (mass_dplyr) {
      lines <- c(
        lines,
        "Both {.pkg MASS} and {.pkg dplyr} are attached; prefer {.fn dplyr::select}."
      )
    }
    envir <- environment()
    formatted <- vapply(
      lines,
      function(x) cli::format_inline(x, .envir = envir),
      character(1L)
    )
    packageStartupMessage(paste(formatted, collapse = "\n"))
    return(invisible())
  }
  msg <- c(
    paste0("DeCovarT ", ver, ": covariance-aware bulk deconvolution."),
    "Main entry points: fit_decovart() and deconvolute_ratios()."
  )
  if (mass_dplyr) {
    msg <- c(
      msg,
      "Note: both MASS and dplyr are attached; dplyr::select() may be masked",
      "by MASS::select(). Prefer dplyr::select() / MASS::select() explicitly."
    )
  }
  packageStartupMessage(paste(msg, collapse = "\n"))
}
