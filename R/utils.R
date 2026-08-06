#' Project estimated proportions onto the unit simplex
#'
#' @description
#' Clips numerical under/overflow and renormalises so that
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1} and
#' \eqn{\boldsymbol{p}\ge\mathbf{0}}.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#'
#' @return Numeric vector on the simplex \eqn{\Delta^{J-1}}.
#' @export
enforce_identifiability <- function(p) {
  machine_limit <- .Machine$double.eps
  p[p < 100 * machine_limit] <- 0 # remove negative ratios
  p <- p / sum(p)
  p[p > 1 - 100 * machine_limit] <- 1 # ensure unit simplex constraint
  return(p)
}

is_positive_definite <- function(expression, tol = 1e-6) {
  eigen_values <- eigen(expression, symmetric = TRUE)$values
  return(all(eigen_values >= -tol * abs(eigen_values[1])))
}

#' Normalised Shannon entropy of a discrete distribution
#'
#' @description
#' For \eqn{\boldsymbol{p}} on the simplex (after dropping zeros),
#' \deqn{
#'   H(\boldsymbol{p})
#'   =
#'   -\sum_{j}p_j\log_{J'}p_j
#'   \in[0,1],
#' }
#' with \eqn{J'} the number of positive masses (1 = uniform).
#'
#' @param ratios Numeric vector \eqn{\boldsymbol{p}}.
#' @return Scalar entropy in \eqn{[0,1]}.
#' @export
compute_shannon_entropy <- function(ratios) {
  if (min(ratios) < 0 || max(ratios) > 1) {
    stop("Probabilities must be stricly included between 0 and 1")
  }

  # normalization process + remove NULL components, as information fisher is not modified by empty classes
  ratios <- ratios[ratios != 0]
  ratios <- ratios / sum(ratios)

  # entropy included between 0 (one component storing all information) and 1 (uniform distribution, balanced classes)
  return(-sum(ratios * logb(ratios, base = length(ratios))))
}

#' Average pairwise overlap of a Gaussian mixture
#'
#' @description
#' Approximates mixture overlap from pairwise misclassification probabilities
#' returned by [MixSim::overlap()], weighted by
#' \eqn{p_i\Omega_{ij}+p_j\Omega_{ji}} and averaged over pairs.
#'
#' @param true_theta List with `p`, `mu`, `sigma` as in a GMM
#'   \eqn{(\boldsymbol{p},\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})}.
#' @param k Number of components \eqn{J} (default `length(true_theta$p)`).
#' @return Scalar average overlap.
#' @export
compute_average_overlap <- function(true_theta, k = length(true_theta$p)) {
  # generate relevant values for the computation of the overlap
  misclassif_mat <- MixSim::overlap(
    Pi = true_theta$p,
    Mu = as.matrix(true_theta$mu),
    S = as.matrix(true_theta$sigma)
  )$OmegaMap
  pairwise_overlap <- c()
  p <- true_theta$p

  # generate the average of pairwise overlaps
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      pairwise_overlap <- c(
        pairwise_overlap,
        misclassif_mat[i, j] * p[i] + misclassif_mat[j, i] * p[j]
      )
    }
  }
  return(mean(pairwise_overlap))
}
