CheckLazyDataCompression <- function(pkg) {
  pkg_name <- sub("_.*", "", pkg)
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
