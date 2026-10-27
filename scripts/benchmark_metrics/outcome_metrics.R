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

# Aitchison (clr Euclidean) distance on the simplex
#
.aitchison_distance <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  if (length(p) != length(p_hat) || length(p) < 2L) {
    return(NA_real_)
  }
  if (any(!is.finite(p)) || any(!is.finite(p_hat))) {
    return(NA_real_)
  }
  eps <- .Machine$double.eps
  p <- pmax(p, eps)
  p_hat <- pmax(p_hat, eps)
  p <- p / sum(p)
  p_hat <- p_hat / sum(p_hat)
  clr <- function(x) {
    log(x) - mean(log(x))
  }
  sqrt(sum((clr(p) - clr(p_hat))^2))
}
