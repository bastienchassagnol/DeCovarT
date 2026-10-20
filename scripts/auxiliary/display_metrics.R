# Display helpers used by figure scripts and not by the package API.

simulation_outcome <- function(numerical_converged, theoretical_converged) {
  num <- as.logical(numerical_converged)
  th <- as.logical(theoretical_converged)
  n <- length(num)
  if (length(th) != n) {
    th <- rep_len(th, n)
  }
  out <- rep("numerical_failure", n)
  ok <- !is.na(num) & num
  success <- ok & !is.na(th) & th
  out[success] <- "theoretical_success"
  out[ok & !success] <- "theoretical_failure"
  out
}

relative_log10_likelihood <- function(z, z_ref = NULL) {
  d <- dim(z)
  z <- as.numeric(z)
  if (is.null(z_ref) || !is.finite(z_ref)) {
    finite <- z[is.finite(z)]
    z_ref <- if (length(finite) == 0L) {
      0
    } else {
      max(finite)
    }
  }
  rel <- exp(z - z_ref)
  rel[!is.finite(rel) | rel <= 0] <- .Machine$double.xmin
  out <- log10(pmax(rel, .Machine$double.xmin))
  if (is.null(d)) {
    return(out)
  }
  array(out, dim = d)
}
