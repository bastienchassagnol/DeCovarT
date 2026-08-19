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
  imports <- tryCatch(
    utils::packageDescription(pkgname, fields = "Imports"),
    error = function(e) NA_character_
  )
  import_pkgs <- character()
  if (!is.na(imports) && nzchar(imports)) {
    import_pkgs <- trimws(strsplit(imports, ",", fixed = TRUE)[[1L]])
    import_pkgs <- sub("\\s*\\(.*$", "", import_pkgs)
  }

  msg <- c(
    paste0("DeCovarT ", ver, " - covariance-aware bulk deconvolution."),
    "Estimate cellular proportions from bulk RNA-seq using Gaussian",
    "convolutions of purified means and covariances (ALR + MLE).",
    "Main entry points: fit_decovart() and deconvolute_ratios().",
    "Website: https://bastienchassagnol.github.io/DeCovarT/"
  )
  if (length(import_pkgs)) {
    msg <- c(
      msg,
      paste0(
        "Imports (loaded with the package): ",
        toString(import_pkgs),
        "."
      )
    )
  }

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
