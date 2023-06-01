#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

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
  p <- p/sum(p); p[p > 1 - 100 * machine_limit] <- 1 # ensure unit simplex constraint
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
  ordered_components <- do.call(order, t(theta$mu) %>% as.data.frame())
  ordered_theta <- list(
    p = theta$p[ordered_components],
    mu = theta$mu[, ordered_components],
    sigma = theta$sigma[, , ordered_components]  )
  
  # enforce sum-to-one constraint
  ordered_theta$p <- ordered_theta$p / sum(ordered_theta$p) # ordered_theta$p[k] <- 1 - sum(ordered_theta$p[-k])
  return(ordered_theta)
}


#' Check whether the estimation has been trapped in the boundary space
#'
#' * Function `check_parameters` asserts at each step of the maximisation,
#' we do not fall in a degenerate case or a non invertible. This especially occurs when one of 
#' the ratios converge to 0 or 1, implying to decrease by a factor 1 the dimension.
#'
#' @param p the ratios estimated
#' @return a boolean, set to FALSE whenever one of the ratios is numerically NULL

#' @export

check_parameters <- function(p) {
  machine_limit <- .Machine$double.eps
  if (any(p < 100 * machine_limit | p > 1 - 100 * machine_limit)) {
    warning(paste0("Cell ratio with index ", which(p < 100 * machine_limit), " is missing in the sample."))
    return(FALSE)
  }
  else {return(TRUE)}
}



isConstant <- function(x) {
  return(length(unique(x))==1)
}

is_continuous <- function(var) {
  return (any(c(is(var, "numeric"), is(var, "integer"), is(var, "Date"),
            is(var, "POSIXct"), is(var, "POSIXt"))))
}

# compute_interval <- function(x) {
#   return(purrr::map_dbl(strsplit(x, split = " to "), ~as.numeric(.x) %>% mean() %>% round()))
# }

is_positive_definite <- function(expression, tol=1e-6) {
  eigen_values <- eigen(expression, symmetric = TRUE)$values # already sorted by decreasing order
  return(all(eigen_values >= -tol * abs(eigen_values[1])))
}

# trace operator
tr <- function(A) {
  return(sum(diag(A)))
}

# matrix raised to power, R operator
`%^%` <- function(A, int_pow) {
  if(!is.integer(int_pow)) {
    stop("Power must be an integer")
  }
  else {
    for (i in 1:(int_pow-1)) {
      A <- A %*% A
    }
  }
  return(A)
}




# get the true stand deviation of a sample
theoric_sd<- function(x) {
  n <- length(x)
  standard_deviation <- ifelse(n==1, 0, sd(x) * sqrt((n-1)/n))
  return(standard_deviation)
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


#' Compute the Hellinger distance
#' @param mu1,Sigma1 mean and (co)variance of the first Gaussian distribution (can be either univariate or multivariate)
#' @param mu2,Sigma2 mean and (co)variance of the second Gaussian distribution
#' @return the Hellinger distance between two multivariate Gaussian distributions
hellinger <- function(mu1, Sigma1, mu2, Sigma2) {
  p <- length(mu1)
  d <- mu1 - mu2
  # in univariate dimension
  if (p == 1) {
    vars <- Sigma1^2 + Sigma2^2
    bc <- sqrt(2 * Sigma1 * Sigma2 / vars) * exp((-1 / 4) * d^2 / vars) # Bhattacharyya coefficient
    return(sqrt(abs(1 - bc)))
  }
  # in multivariate dimension
  else {
    vars <- (Sigma1 + Sigma2) / 2
    hell_dist <- det(Sigma1)^(1 / 4) * det(Sigma2)^(1 / 4) / det(vars)^(1 / 2) *
      exp((-1 / 8) * .maha_distance(d, vars))
    return(sqrt(1 - hell_dist) %>% as.numeric())
  }
}

hellinger_average <- function(p, signature_matrix, cov_matrix) {
  pairwise_hellinger <- c()
  for (i in 1:(ncol(signature_matrix)-1)) {
    for (j in (i+1):ncol(signature_matrix)) {
      # compute Hellinger distance between component i and component j, weighted by their respective proportion
      hellinger_value <- hellinger(mu1 = p[i] * signature_matrix[, i], Sigma1 = p[i]^2 * cov_matrix[,,i],
                                   mu2 = p[j] * signature_matrix[, j], Sigma2 = p[j]^2 * cov_matrix[,,j])
      pairwise_hellinger <- c(pairwise_hellinger, hellinger_value)
    }
  }
  return(mean(pairwise_hellinger))
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
      pairwise_overlap <- c(pairwise_overlap, misclassif_mat[i, j] * p[i] + misclassif_mat[j, i] * p[j])
    }
  }
  return(mean(pairwise_overlap))
}






