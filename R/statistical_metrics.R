#' Repair a numeric vector onto the unit simplex
#'
#' @description
#' Clips numerical under/overflow and renormalises so that
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1} and
#' \eqn{\boldsymbol{p}\ge\mathbf{0}}. This is a **repair /
#' renormalisation** step for estimated proportions, not a Euclidean
#' projection onto the simplex and not a statistical-identifiability
#' constraint.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param tolerance Non-negative tolerance for treating entries as zero
#'   (default `100 * .Machine$double.eps`).
#'
#' @return Numeric vector on the simplex \eqn{\Delta^{J-1}}.
#'
#' @seealso [compositions::clo()] for compositional closure.
#' @export
repair_simplex <- function(p, tolerance = 100 * .Machine$double.eps) {
  if (!is.numeric(p) || length(p) == 0L || anyNA(p)) {
    stop("`p` must be a non-empty numeric vector without missing values.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0) {
    stop("`tolerance` must be a single non-negative number.")
  }

  p[abs(p) <= tolerance] <- 0

  if (any(p < 0)) {
    stop("`p` contains negative values beyond numerical tolerance.")
  }

  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    stop("`p` must have a positive, finite sum.")
  }

  p <- p / total
  p[abs(p - 1) <= tolerance] <- 1
  p
}

#' Normalised Shannon entropy of a discrete distribution
#'
#' @description
#' For a probability vector \eqn{\boldsymbol{p}\in\Delta^{J-1}} over
#' \eqn{J} classes (here: cell types), the Shannon entropy is
#' \deqn{
#'   H(\boldsymbol{p})
#'   =
#'   -\sum_{j=1}^{J} p_j \log p_j,
#' }
#' with the convention \eqn{0\log 0 = 0}. Dividing by the maximum
#' entropy \eqn{\log J} (uniform over all \eqn{J} classes) yields
#' **Pielou's evenness**
#' \deqn{
#'   H^{\star}(\boldsymbol{p})
#'   =
#'   \frac{H(\boldsymbol{p})}{\log J}
#'   \in[0,1],
#' }
#' so \eqn{H^{\star}=0} for a Dirac mass on one type and
#' \eqn{H^{\star}=1} for the uniform distribution over the \eqn{J}
#' cell types. Zero masses are dropped only inside the sum; the
#' normaliser still uses the original class count \eqn{J}.
#'
#' This is preferable to changing the logarithm base to the number of
#' *positive* masses \eqn{J'}, which would renormalise only on the
#' support and hide sparsity relative to the full panel of cell types.
#'
#' @section Illustration:
#' Panel B of the figure contrasts a peaked (type-specific, low entropy)
#' share vector with a flat (multi-type, high entropy) one for \eqn{J=5}
#' classes.
#'
#' \if{html}{\figure{gini_vs_entropy_specificity.png}{options: width="100\%"
#'   alt="Gini versus Shannon entropy"}}
#' \if{latex}{\figure{gini_vs_entropy_specificity.png}{options: width=5.5in}}
#'
#' @param ratios Numeric vector of cellular proportions
#'   \eqn{\boldsymbol{p}\in\Delta^{J-1}} (length \eqn{J}).
#' @return Scalar normalised entropy \eqn{H^{\star}\in[0,1]}.
#' @export
#' @examples
#' compute_shannon_entropy(c(1, 0, 0)) # Dirac -> 0
#' compute_shannon_entropy(rep(1 / 3, 3)) # uniform -> 1
compute_shannon_entropy <- function(ratios) {
  if (!is.numeric(ratios) || length(ratios) == 0L || anyNA(ratios)) {
    stop("`ratios` must be a non-empty numeric vector without missing values.")
  }
  if (min(ratios) < 0 || max(ratios) > 1) {
    stop("Probabilities must be strictly included between 0 and 1.")
  }

  # J = number of cell types in the full panel (before dropping zeros)
  J <- length(ratios)
  if (J == 1L) {
    return(0)
  }

  ratios <- ratios[ratios > 0]
  ratios <- ratios / sum(ratios)
  # Pielou evenness: H / log(J)
  -sum(ratios * log(ratios)) / log(J)
}

#' Average pairwise overlap of a Gaussian mixture
#'
#' @description
#' Returns MixSim's **BarOmega**: the unweighted mean of pairwise
#' overlaps
#' \deqn{
#'   \overline{\omega}
#'   =
#'   \frac{2}{J(J-1)}
#'   \sum_{1\le j<\ell\le J}
#'   \bigl(\Omega_{j\ell}+\Omega_{\ell j}\bigr)
#'   \in[0,1]
#' }
#' (up to the MixSim numerical convention), where
#' \eqn{\Omega_{j\ell}=\Pr_{X\sim f_j}(X\text{ classified as }\ell)}
#' already uses the mixture weights \eqn{\boldsymbol{p}} inside the
#' Bayes / MAP rule of [MixSim::overlap()]. Do **not** multiply the
#' directional masses by \eqn{p_j} again.
#'
#' @param true_theta List with `p` (length \eqn{J}), `mu` (\eqn{G\times J}
#'   mean matrix) and `sigma` (\eqn{G\times G\times J} covariance array),
#'   as in a GMM
#'   \eqn{(\boldsymbol{p},\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})}.
#' @param J Number of cell types (components). Defaults to
#'   `length(true_theta$p)`.
#' @return Scalar average pairwise overlap (MixSim `BarOmega`).
#' @export
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_overlap(theta)
compute_average_overlap <- function(true_theta, J = length(true_theta$p)) {
  stopifnot(is.list(true_theta))
  p <- true_theta$p
  mu <- as.matrix(true_theta$mu)
  sigma <- true_theta$sigma

  if (length(p) != J) {
    stop("`J` must equal length(true_theta$p).")
  }
  if (is.null(dim(sigma)) || length(dim(sigma)) != 3L) {
    stop("`true_theta$sigma` must be a G x G x J array.")
  }

  G <- dim(sigma)[[1L]]
  n_celltypes <- dim(sigma)[[3L]]
  if (dim(sigma)[[2L]] != G || n_celltypes != length(p)) {
    stop("`sigma` dims must be G x G x J with J = length(p).")
  }

  # DeCovarT stores mu as G x J; MixSim::overlap expects Mu as J x G
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_mixsim <- t(mu)
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_mixsim <- mu
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  MixSim::overlap(Pi = p, Mu = mu_mixsim, S = sigma)$BarOmega
}
