#' Compare lazy-data compression sizes (gzip / bzip2 / xz)
#'
#' Installs a built package three times with different
#' `--data-compress=` options and reports compressed Rdata sizes in KiB.
#'
#' @param pkg Path to a source package tarball or directory.
#' @return Named numeric vector of ceiling sizes in KiB.
#' @keywords internal
CheckLazyDataCompression <- function(pkg) {
  pkg_name <- sub("_.*", "", basename(pkg))
  lib <- tempfile()
  dir.create(lib)
  zs <- c("gzip", "bzip2", "xz")
  res <- integer(3)
  names(res) <- zs
  for (z in zs) {
    opts <- c(
      paste0("--data-compress=", z),
      "--no-libs",
      "--no-help",
      "--no-demo",
      "--no-exec",
      "--no-test-load"
    )
    utils::install.packages(
      pkg,
      lib,
      INSTALL_opts = opts,
      repos = NULL,
      quiet = TRUE
    )
    res[z] <- file.size(file.path(lib, pkg_name, "data", "Rdata.rdb"))
  }
  ceiling(res / 1024)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
  print(CheckLazyDataCompression(args[[1L]]))
}
