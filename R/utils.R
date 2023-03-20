isConstant <- function(x) {
  return(length(unique(x))==1)
}

is_continuous <- function(var) {
  return (any(c(is(var, "numeric"), is(var, "integer"), is(var, "Date"),
            is(var, "POSIXct"), is(var, "POSIXt")))) 
}

compute_interval <- function(x) {
  return(purrr::map_dbl(strsplit(x, split = " to "), ~as.numeric(.x) %>% mean() %>% round()))
}

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


# equals 0 when two distances are equal, and 1 when they have distinct support
hellinger <- function (mu1, Sigma1, mu2, Sigma2) {
  Sigma1 <- p1^2*Sigma1; Sigma2 <- p2^2*Sigma2
  p <- length(mu1);   d <- mu1 - mu2
  vars <- (Sigma1 + Sigma2)/2
  # in univariate dimension
  if (p == 1) {
    if (abs(Sigma1) < .Machine$double.eps | abs(Sigma2) < 
        .Machine$double.eps) {
      stop("At least one variance is zero")
    }
    d <- sqrt(Sigma1 * Sigma2/ vars) * 
      exp((-1/4) * d^2/(2*vars))
    return(sqrt(1 - d))
  }
  # in multivariate dimension
  else {
    if (abs(det(Sigma1)) < .Machine$double.eps | abs(det(Sigma2)) < 
        .Machine$double.eps) {
      stop("One of the sample variances is degenerate")
    }
    hell_dist <- det(Sigma1)^(1/4) * det(Sigma2)^(1/4) / det(vars)^(1/2) * 
      exp((-1/8) * maha(d, vars))  
    return(sqrt(1 - hell_dist) %>% as.numeric())
  }
}

# same result as Hellinger function, but borrowed from gaussDiff package
hellinger_v2 <- function(mu1,Sigma1,mu2,Sigma2,s=0.5){
  N <- nrow(Sigma1)
  I <- diag(1,N)
  inv.Sigma1 <- solve(Sigma1)
  inv.Sigma2 <- solve(Sigma2)
  sig1inv.sig2 <- inv.Sigma1%*%Sigma2
  sig2inv.sig1 <- inv.Sigma2%*%Sigma1
  d <- det(s*I+(1-s)*sig1inv.sig2)**(-s/2)*
    det((1-s)*I+s*sig2inv.sig1)**(-(1-s)/2)*
    exp(0.5*(maha(s*inv.Sigma2%*%mu2+(1-s)*inv.Sigma1%*%mu1,s*inv.Sigma2+(1-s)*inv.Sigma1)-
               s*maha(mu2,Sigma2)-(1-s)*maha(mu1,Sigma1)))
  return(as.vector(1-d) %>% sqrt())
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






