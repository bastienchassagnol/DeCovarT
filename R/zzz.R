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
  msg <- c(
    paste0("DeCovarT ", ver, ": covariance-aware bulk deconvolution."),
    "Main entry points: fit_decovart() and deconvolute_ratios()."
  )

  # Highlight a common masking risk when MASS is attached after tidyselect APIs.
  attached <- paste0("package:", search())
  if ("package:MASS" %in% attached && "package:dplyr" %in% attached) {
    msg <- c(
      msg,
      "Note: both MASS and dplyr are attached; dplyr::select() may be masked",
      "by MASS::select(). Prefer dplyr::select() / MASS::select() explicitly."
    )
  }

  packageStartupMessage(paste(msg, collapse = "\n"))
}
