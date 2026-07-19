#' Check whether the estimation has been trapped in the boundary space
#'
#' * Function `enforce_identifiability` normalises the ratios to sum to 1,
#' and set ratios close to zero or 1, to the limits of the boundary space
#'
#' @param p the estimated ratios
#'
#' @return a numeric vector of size \eqn{J}, but assuring that negative
#' ratios were set to 0, and the unit simplex constraint is endorsed

#' @export

enforce_identifiability <- function(p) {
  machine_limit <- .Machine$double.eps
  p[p < 100 * machine_limit] <- 0 # remove negative ratios
  p <- p / sum(p)
  p[p > 1 - 100 * machine_limit] <- 1 # ensure unit simplex constraint
  return(p)
}


#' Control parameters output
#'
#' This step ensures that the estimates returned are uniquely ordered by
#' partial ordering on the means, and that the sum-o-one constraint, that may
#' be violated by numerical artefacts, is enforced
#'
#'
#' @param theta estimation of the parameters returned either by an initialisation
#' algorithm or by an EM algorithm on an univariate or multivariate GMM MLE estimation
#' * The proportions `p`: \eqn{p} of each component (must be included between 0 and 1, and sum to one overall)
#' * The mean matrix `mu`: \eqn{\mathrm{\mu}=(\mu_{i,j}) \in \mathbb{R}^{n \times k}}, with each column
#' giving the mean values of the variables within a given component
#' * The 3-dimensional covariance matrix array `sigma`: \eqn{\mathrm{\sigma}=(\sigma_{i,j,l}) \in \mathbb{R}^{n \times n \times k}}, with each matrix
#' \eqn{\sigma_{..l}, l \in \{ 1, \ldots, k\}} storing the covariance matrix of a given component,
#' whose diagonal terms correspond to the variance of each variable, and off-terms diagonal elements return the covariance matrix
#'
#' @return a list of the estimates, uniquely identified, by ranking each component
#' based on the ordering of their means

enforce_parameter_identifiability <- function(theta) {
  k <- length(theta$p)
  # first component = one whose means are smaller on the first dimension
  # if equality, redirect to second column, etc...
  ordered_components <- do.call(order, t(theta$mu) |> as.data.frame())
  ordered_theta <- list(
    p = theta$p[ordered_components],
    mu = theta$mu[, ordered_components],
    sigma = theta$sigma[,, ordered_components]
  )

  # enforce sum-to-one constraint
  ordered_theta$p <- ordered_theta$p / sum(ordered_theta$p) # ordered_theta$p[k] <- 1 - sum(ordered_theta$p[-k])
  return(ordered_theta)
}


is_positive_definite <- function(expression, tol = 1e-6) {
  eigen_values <- eigen(expression, symmetric = TRUE)$values # already sorted by decreasing order
  return(all(eigen_values >= -tol * abs(eigen_values[1])))
}

#' Compute the shannon entropy of a discrete distribution, normalised from 0 to 1 (equibalanced classes)
#'
#' @author Bastien CHASSAGNOL
#'
#' @param ratios vector of the proportions of the mixture
#' @return the entropy score
#' @export

compute_shannon_entropy <- function(ratios) {
  if (min(ratios) < 0 | max(ratios) > 1) {
    stop("Probabilities must be stricly included between 0 and 1")
  }

  # normalization process + remove NULL components, as information fisher is not modified by empty classes
  ratios <- ratios[ratios != 0]
  ratios <- ratios / sum(ratios)

  # entropy included between 0 (one component storing all information) and 1 (uniform distribution, balanced classes)
  return(-sum(ratios * logb(ratios, base = length(ratios))))
}

#' Compute average overlap between components
#'
#' Internally, it is the function MixSim::overlap which is used to generate an approximate pairwise
#' probability to wrongfully assign one component to another. Unfortunately, this function does not
#' the global overlap that we approximate there by averaging pairwise overlaps + compute the overlap
#' between two components accounting for their respective proportions in the mixture
#'
#' @author Bastien CHASSAGNOL
#'
#' @param true_theta the parameters of the GMM
#' @param k the number of components
#' @return the average overlap
#'
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
