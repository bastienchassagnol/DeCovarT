###########################################################################
###########################################################################
###                                                                     ###
###                     OWN DECONVOLUTION ALGORITHM                     ###
###                                                                     ###
###########################################################################
###########################################################################

#################################################################
##                      utility functions                      ##
#################################################################

# mapping function, from theta to ratios p
mapping_function <- function(theta) {
  p <- c(exp(theta[1:length(theta)]),1)
  return (p/sum(p))
}

mapping_function_univariate <- function(theta, index=1) {
  return(mapping_function(theta)[index])
}

# reciprocal of the mapping function
inverse_mapping_function <- function(p) {
  num_cells <- length(p)
  return(log(p[1:num_cells-1]/p[num_cells]))
}

# compute the mahalabonis distance
.maha_distance <- function(x,A){
  d <- t(x)%*%solve(A)%*%x # solve A returns the reverted function
  return(d %>% as.numeric()) # supposed to be a scalar
}

# compute the dot product distance, quite similar to maha,
# but without reverting the matrix, and with possibly a different second vector
.dot_product <- function(x,A, y=x){
  d <- t(x)%*%A%*%y # solve A returns the reverted function
  return(d %>% as.numeric()) # supposed to be a scalar
}

.compute_global_variance <- function(p, Sigma) {
  ###  Sigma and TensorA packages  
  # global_cov <- matrix(0, nrow = dim(Sigma)[1], ncol=dim(Sigma)[2])
  # for (j in 1:length(p)) {
  #   global_cov <- global_cov + p[j]^2*Sigma[,,j]
  # }
  # 
  # # convert to tensors
  # num_genes <- dim(Sigma)[1]; num_celltypes <- length(p)
  # p <- to.tensor(p,c("j"=num_celltypes)); Sigma <- to.tensor(c(Sigma),c("i"=num_genes,"k"=num_genes,"j"=num_celltypes)) 
  
  # Einstein summation, then conversion back to matrix
  # global_cov <- p^2 %e% Sigma %>% to.matrix.tensor(j="i") %>% as.matrix()
  
  global_cov <- tensor::tensor(p^2, Sigma, alongA = 1, alongB = 3)
  
  # check the validity of the covariance computed
  if (!is_positive_definite(global_cov) | det(global_cov) ==0) 
    stop("Global covariance matrix is either not invertible or not positive definite")
  # with tensor
  return(global_cov)
}


#################################################################
##                   first order derivatives                   ##
#################################################################

# log-likelihood multivariate function
loglik_multivariate <- function(p, y, X, Sigma) {
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  log_lik <- -log(det(global_cov_matrix)) - 1/2 * .maha_distance(y - X %*% p, global_cov_matrix) %>% as.numeric()
  return (log_lik)
}


loglik_multivariate_constrained <- function(theta, y, X, Sigma) {
  # switch from variable
  p <- mapping_function(theta); log_lik <- loglik_multivariate(p, y, X, Sigma)
  return (log_lik)
}


# Jacobian mapping function 
jacobian_mapping_function <- function(theta) {
  denominator <- (sum(exp(theta)) + 1)^2; size_var <- length(theta)
  jacobian_matrix <- matrix(0, nrow = size_var, ncol = size_var)
  for (i in 1:size_var) {
    for (j in i:size_var) {
      # diagonal elements 
      if(i==j) {
        jacobian_matrix[i, j] <- (exp(theta[i])*(sum(exp(theta[-i]))+1))
      }
      else {
        jacobian_matrix[i, j] <- -exp(theta[i])*exp(theta[j])
      }
    }
  }
  # ensure symmetry
  jacobian_matrix[lower.tri(jacobian_matrix)] <- jacobian_matrix[upper.tri(jacobian_matrix)]
  jacobian_matrix <- rbind(jacobian_matrix, -exp(theta))
  return(jacobian_matrix/denominator)
}

gradient_loglik_unconstrained <- function(p, y, X, Sigma) {
  # compute general covariance and its reverse
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  global_precision_matrix <- solve(global_cov_matrix)
  
  # compute the gradient itself
  gradient_unconstrained <- c()
  for (j in 1:length(p)) {
    gradient_unconstrained <- c(gradient_unconstrained,
                                -2*p[j]*tr(global_precision_matrix %*% Sigma[,,j]) +
                                  .dot_product(y - X %*% p, global_precision_matrix, X[,j]) +
                                  p[j] * .dot_product(y - X %*% p, global_precision_matrix %*% Sigma[,,j] %*% global_precision_matrix)) 
    
  }
  return(gradient_unconstrained)
}

gradient_loglik_constrained <- function (theta, y, X, Sigma) {
  p <- mapping_function(theta)
  gradient_constrained <- gradient_loglik_unconstrained(p, y, X, Sigma) %*% jacobian_mapping_function(theta)
  return(gradient_constrained)
}


##################################################################
##                   second order derivatives                   ##
##################################################################
# Hessian mapping function 
hessian_mapping_function <- function(theta) {
  A <- sum(exp(theta)) + 1; denominator <- A^3
  size_var <- length(theta); J <- size_var + 1
  hessian_array <- array(0, dim = c(size_var, size_var, J))
  for (i in 1:size_var) {
    B <- sum(exp(theta[-i]))+1
    # other p_j, j< J Hessian derivation
    for (j in 1:size_var) {
      for (k in j:size_var) {
        # diagonal elements 
        if(j==k) {
          if (i==j) {
            hessian_array[i, i, i] <- B * exp(theta[i])*(B - exp(theta[i])) # condition d)
          }
          else {
            hessian_array[j, j, i] <- exp(theta[i]) * exp(theta[j]) * (-A + 2 * exp(theta[j])) # condition c)
          }
        } 
        # off diagonal terms
        else {
          if (!length(intersect(i, c(j, k)))) { # all indexes are different, situation b) 
            # alternative condition setting: (i!=j)!=k
            hessian_array[j, k, i] <- 2 * exp(theta[i]) * exp(theta[j])*exp(theta[k])
          }
          else {
            l <- setdiff(c(j, k), i) # situation a), l is the operator, either k or j, different from i
            hessian_array[j, k, i] <- exp(theta[i]) * exp(theta[l])*(-B + exp(theta[i]))
          }
        }
      }
    }
    # ensure symmetry
    hessian_array[,,i][lower.tri(hessian_array[,,i])] <- hessian_array[,,i][upper.tri(hessian_array[,,i])]
  }
  # last p_J Hessian component
  for (j in 1:size_var) {
    for (k in j:size_var) {
      # diagonal elements 
      if(j==k) {
        hessian_array[j, j, J] <- exp(theta[j])*(-A+2*exp(theta[j])) # condition e)
      }
      else {
        hessian_array[j, k, J] <- 2 * exp(theta[j])*exp(theta[k]) # condition f)
      }
    }
  }
  hessian_array[,,J][lower.tri(hessian_array[,,J])] <- hessian_array[,,J][upper.tri(hessian_array[,,J])]
  return(hessian_array/denominator)
}

hessian_loglik_unconstrained <- function(p, y, X, Sigma) {
  num_celltypes <- length(p); hessian_unconstrained <- matrix(0, nrow = num_celltypes, ncol = num_celltypes)
  global_precision_matrix <- .compute_global_variance (p, Sigma) %>% solve()
  for (i in 1:num_celltypes) {
    for (j in i:num_celltypes) {
      hessian_unconstrained[i, j] <- 4 * p[i] * p[j] * 
        tr(global_precision_matrix %*% Sigma[,,i] %*% global_precision_matrix %*% Sigma[,,j]) - 
        .dot_product(X[,i], global_precision_matrix, X[,j]) -
        2 * p[i] * .dot_product(y-X %*% p, 
                                global_precision_matrix %*% Sigma[,,i] %*% global_precision_matrix, X[,j]) - 
        2 * p[j] * .dot_product(y-X %*% p, 
                                global_precision_matrix %*% Sigma[,,j] %*% global_precision_matrix, X[,i]) - 
        4 * p[i] * p[j] * .dot_product(y-X %*% p, 
                                       global_precision_matrix %*% Sigma[,,j] %*% 
                                         global_precision_matrix %*% Sigma[,,i] %*%
                                         global_precision_matrix)
      if (i==j) { # add diagonal terms
        hessian_unconstrained[i,i] <- hessian_unconstrained[i,i] - 
          2 * tr(global_precision_matrix %*% Sigma[,,i]) + 
          .dot_product(y-X %*% p, global_precision_matrix %*% Sigma[,,i] %*% global_precision_matrix)
      }
    }
  }
  # enforce symmetry
  hessian_unconstrained[lower.tri(hessian_unconstrained)] <- hessian_unconstrained[upper.tri(hessian_unconstrained)]
  return(hessian_unconstrained)
}

hessian_loglik_constrained <- function (theta, y, X, Sigma) {
  p <- mapping_function(theta)
  # t(J_psi) X H_log X J_psi + sum over number of ratios of grad_log X H_psi
  hessian_constrained <- t(jacobian_mapping_function(theta)) %*% hessian_loglik_unconstrained(p, y, X, Sigma) %*% jacobian_mapping_function(theta) +
    tensor::tensor(A = gradient_loglik_unconstrained(p, y, X, Sigma), B = hessian_mapping_function(theta), alongA = 1, alongB = 3) %>% as.matrix()
  return(hessian_constrained)
}


#################################################################
##                    iterated descent algorithms              ##
#################################################################


# deconvolution, using general optimisation function
deconvolute_ratios_corr_decon <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_values <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  estimated_ratios <- optim(par=initial_values,fn=loglik_multivariate, gr=NULL,y=y, X=X, Sigma=Sigma,
                            control=list(fnscale=-1),method="BFGS")$par %>% 
    mapping_function() %>% stats::setNames(colnames(X)) # ensure non-negativity constraint
  
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
  return(metrics_scores)
}
