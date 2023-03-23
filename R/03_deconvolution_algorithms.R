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
#' Title
#'
#' @param theta 
#'
#' @return
#' @export
#'

mapping_function <- function(theta) {
  p <- c(exp(theta[1:length(theta)]),1)
  return (p/sum(p))
}

#' Title
#'
#' @param theta 
#' @param index 
#'
#' @return
#' @export
#'

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
  # if (!is_positive_definite(global_cov) | det(global_cov) ==0)
  #   stop("Global covariance matrix is either not invertible or not positive definite")
  # with tensor
  return(global_cov)
}


#################################################################
##                   first order derivatives                   ##
#################################################################

# log-likelihood multivariate function
loglik_multivariate <- function(p, y, X, Sigma) {
  # print(paste("estimated value of p is ", p))
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  # print(global_cov_matrix); print("\n\n")
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

# deconvolution, using simulated annealing function
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_simulated_annealing <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # gr is not used in the simulated annealing approach
  estimated_theta <- optim(par=initial_theta,fn=loglik_multivariate_constrained, y=y, X=X, Sigma=Sigma,
                            control=list(fnscale=-1, maxit=10000),method="SANN")$par
  estimated_p <- mapping_function(estimated_theta) %>% stats::setNames(colnames(X))

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


# deconvolution, using L-BFGS-B function
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_LBFGS <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  estimated_p <- optim(par=initial_p,fn=loglik_multivariate, gr=gradient_loglik_unconstrained,
                       y=y, X=X, Sigma=Sigma,
                       control=list(fnscale=-1,maxit=200, lmm=10, factr=1e-3),method="L-BFGS-B",
                       lower=rep(0, length(initial_p)), upper=rep(1, length(initial_p)))$par %>%
    stats::setNames(colnames(X)) # ensure non-negativity constraint + set fnscale to -1, since maximization is searched for

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

# deconvolution, using constrained barrier
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param epsilon 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_constrOptim <- function(y, X, Sigma, true_ratios=NULL, epsilon=10^-8, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  ui <- diag(nrow=length(initial_p)); ci <- rep(0, length(initial_p)) # encode the non-negativity constraint

  # encode the sum-to-one constraint, converting inequality into equality
  # ui <- diag(nrow=length(initial_p)-1); ci <- rep(0, length(initial_p))
  # ui <- rbind(ui, rep(1, length(initial_p)), rep(-1, length(initial_p))); ci <- c(ci,1,-1)
  ui <- rbind(ui, rep(1, length(initial_p))); ci <- c(ci,1) # can not make an exact equality


  # add some perturbation to start the iteration
  estimated_p <- stats::constrOptim(theta=initial_p + epsilon,f=loglik_multivariate, grad=gradient_loglik_unconstrained,
                                    ui=ui, ci=ci, control=list(fnscale=-1,maxit=200, lmm=10, factr=1e-3),
                                    method="BFGS", outer.iterations=200, outer.eps=1e-4,  y=y, X=X, Sigma=Sigma)$par %>%
    stats::setNames(colnames(X))

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


# deconvolution, using standard constrained Newton-Raphson approach
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_nlm <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)


  # with nlm method (use minus to maximise, instead of minimise)
  # f_list <- function(p, y, X, Sigma) {- loglik_multivariate_constrained(p, y, X, Sigma)}  # function(x) -loglik_multivariate_constrained(x)
  # attr(f_list, "gradient") <- function(p, y, X, Sigma) {- gradient_loglik_constrained(p, y, X, Sigma)}
  # attr(f_list, "hessian") <- function(p, y, X, Sigma) {-hessian_loglik_constrained(p, y, X, Sigma)}
  # estimated_theta <- stats::nlm(f=f_list, p=initial_theta, y=y, X=X, Sigma=Sigma, stepmax=1,
  #                                           ndigit=8, gradtol=1e-4, steptol=1e-4, iterlim=200)$estimate

  # with nlmimb package method (outdated, but works well for our scenario)
  estimated_theta <- stats::nlminb(start=initial_theta, objective=function(p, y, X, Sigma) {- loglik_multivariate_constrained(p, y, X, Sigma)},
                                   gradient=function(p, y, X, Sigma) {- gradient_loglik_constrained(p, y, X, Sigma)},
                                   hessian=function(p, y, X, Sigma) {-hessian_loglik_constrained(p, y, X, Sigma)},
                                   y=y, X=X, Sigma=Sigma,
                                   control=list(eval.max=10,iter.max=200, rel.tol=1e-4, x.tol=1e-4, xf.tol=1e-4, abs.tol=1e-4))$par

  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) # ensure non-negativity constraint

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)}



# deconvolution, using constrained LM parametrisation
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_DeCoVarT <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # set minimize to false; partialH=2
  estimated_theta <- marqLevAlg::marqLevAlg(b=initial_theta, fn=loglik_multivariate_constrained,
                                            gr=gradient_loglik_constrained, hess = hessian_loglik_constrained,
                                            epsa=1e-04, epsb = 1e-04, epsd = 1e-04, minimize=FALSE, multipleTry = 10,
                                            y=y, X=X, Sigma=Sigma)$b
  if (is.na(estimated_theta)) {
    output_lm <- capture.output(marqLevAlg::marqLevAlg(b=initial_theta, fn=loglik_multivariate_constrained,
                                                       gr=gradient_loglik_constrained, hess = hessian_loglik_constrained,
                                                       epsa=1e-04, epsb = 1e-04, epsd = 1e-04, minimize=FALSE, multipleTry = 10,
                                                       y=y, X=X, Sigma=Sigma)) # add partialH and blinding?
    
    estimated_theta <- output_lm[3] %>% stringr::str_match_all("[0-9,\\.]+") %>% 
      unlist %>% as.numeric # retrieve last estimae then
  }
  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>% 
    check_parameters() # ensure non-negativity constraint and remove numerical underflow

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


# deconvolution, using the most basic optimisation approach
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_basic_optim <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  # set minimize to false
  estimated_p <- stats::optim(par=initial_p, fn=loglik_multivariate, gr=NULL,
                                            y=y, X=X, Sigma=Sigma, method="BFGS",
                                  control=list(fnscale=-1,maxit=200, reltol=1e-4,lmm=10, factr=1e-3))$par %>%
    stats::setNames(colnames(X))

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


# deconvolution, using the constrained gradient approach
#' Title
#'
#' @param y 
#' @param X 
#' @param Sigma 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_constrained_optim <- function(y, X, Sigma, true_ratios=NULL, ...) {
  initial_p <- rep(1/ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  
  initial_theta <- inverse_mapping_function(initial_p)
  
  estimated_theta <- stats::optim(par=initial_theta,fn=loglik_multivariate_constrained, gr = NULL, y=y, X=X, Sigma=Sigma,
                           control=list(fnscale=-1, maxit=200, reltol=1e-4,lmm=10, factr=1e-3, maxit=200),method="BFGS")$par
  estimated_p <- mapping_function(estimated_theta) %>% stats::setNames(colnames(X))
  

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}







############################################################################
############################################################################
###                                                                      ###
###                    OTHER DECONVOLUTION ALGORITHMS                    ###
###                                                                      ###
############################################################################
############################################################################

# Core algorithm (where svm regression is performed for each mixture)
#' Title
#'
#' @param y 
#' @param X 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_CIBERSORT <- function(y, X, true_ratios=NULL, ...){
     #the set of nu values tested for best performance
    range.nu<-seq(0.2, 0.8, 0.3)
    model<-e1071::best.svm(X,y, type="nu-regression",kenel="linear",scale=F,
                      nu = range.nu, tunecontrol=e1071::tune.control(sampling = "fix", fix=0.75)) #determine best fitted model

    e1071::tune.svm(X,y,type="nu-regression",kernel="linear",scale=F,
                    nu = 0.25, tunecontrol=e1071::tune.control(sampling = "fix", fix=0.75))

    model <- e1071::svm(X, y=y)

    #get and normalize coefficients (sum to 1 and no negative coefficients)
    estimated_p <- t(model$coefs) %*% model$SV
    estimated_p[which(estimated_p<0)]<-0
    estimated_p <- estimated_p/sum(estimated_p)

    metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
      dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
    return(metrics_scores)
}

# linear regression
#' Title
#'
#' @param y 
#' @param X 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_abbas <- function(y, X, true_ratios=NULL, ...) {
  estimated_p <- stats::lsfit(X, y, intercept = F)$coefficients

  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p[which(estimated_p<0)]<-0; estimated_p <- estimated_p/sum(estimated_p)
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

# robust linear regression
#' Title
#'
#' @param y 
#' @param X 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'

deconvolute_ratios_monaco <- function(y, X, true_ratios=NULL, ...) {
  estimated_p <- MASS::rlm(y ~ X+ 0, method = c("M"))$coefficients; names(estimated_p) <- colnames(X)

  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p[which(estimated_p<0)]<-0; estimated_p <- estimated_p/sum(estimated_p)
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


# quadratic programming (lsei function, other possibility is to use lsqlin function)
#' Title
#'
#' @param y 
#' @param X 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
deconvolute_ratios_deconRNASeq <- function(y, X, true_ratios=NULL, ...) {

EE <- rep(1, ncol(X)); FF <- 1 # encode the sum-to-one constraint
GG <- diag(nrow=ncol(X)); HH <- rep(0, ncol(X)) # encode the non-negativity constraint

estimated_p <- limSolve::lsei(X, y, EE, FF, GG, HH)$X
metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
  dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
return(metrics_scores)
}

# nnls scoring
#' Title
#'
#' @param y 
#' @param X 
#' @param true_ratios 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
deconvolute_ratios_nnls <- function(y, X, true_ratios=NULL, ...) {
  estimated_p <- nnls::nnls(X, y)$x

  # normalize coefficients (sum to 1, as non-negativity is already enforced)
  estimated_p <- estimated_p/sum(estimated_p); names(estimated_p) <- colnames(X)
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}














# compute summary scores
#' Title
#'
#' @param y 
#' @param X 
#' @param estimated_p 
#' @param true_ratios 
#'
#' @return
#' @export
#'

compute_benchmark_metrics <- function(y, X, estimated_p, true_ratios=NULL) {
  n <- nrow(X); k <- ncol(X)
  df_res <- n - k + 1 # number of free parameters: only k - 1 parameters must be learnt, with sum-to-one constraint
  df_tot <- n - 1 # no intercept for the moment, in the model (so n-1, or n?)
  if(!is.null(true_ratios)) { # when the true parameters are known
    scores <- tibble::tibble(model_mse = Metrics::mse(true_ratios, estimated_p),
                             model_rmse = Metrics::rmse(true_ratios, estimated_p),
                             model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_p)),
                             model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
                             model_mae = Metrics::mae(true_ratios, estimated_p),
                             model_cor = suppressWarnings(cor(true_ratios, estimated_p, method = "pearson")))

  }
  else { # when they are unknown
    predicted_values <- as.vector(X%*% estimated_p)
    scores <- tibble::tibble(model_mse = Metrics::mse(y, predicted_values),
                             model_rmse = Metrics::rmse(y, predicted_values),
                             model_mae = Metrics::mae(y, predicted_values),
                             model_cor = suppressWarnings(cor(y, predicted_values, method = "pearson")))

    # model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_p)),
    # model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
    # # model_ccc= epiR::epi.ccc(y, predicted_values)
  }
  return(scores)
}




# main function, launching all deconvolution algorithms
#' Title
#'
#' @param signature_matrix 
#' @param bulk_expression 
#' @param scaled 
#' @param true_ratios 
#' @param Sigma 
#' @param deconvolution_functions 
#'
#' @return
#' @export
#'

deconvolute_ratios <- function(signature_matrix, bulk_expression, scaled=F, true_ratios=NULL, Sigma=NULL,
                               cores = getOption("mc.cores", parallel::detectCores()), verbose=FALSE,
                               deconvolution_functions=list("QP" =list(FUN=deconvolute_ratios_deconRNASeq,additionnal_parameters=NULL),
                                                            "lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL),
                                                            "rlm"=list(FUN=deconvolute_ratios_monaco,additionnal_parameters=NULL),
                                                            "nnls" = list(FUN=deconvolute_ratios_nnls,additionnal_parameters=NULL))){
  #read in data
  if (!is.matrix (signature_matrix) | is.null(row.names (signature_matrix)))
    stop ("required format for signature is expression matrix, with rownames as genes")
  X <- tibble::as_tibble(signature_matrix, rownames="GENE_SYMBOL")
  if (!is.matrix (bulk_expression) | is.null(row.names (bulk_expression)))
    stop ("required format for mixture is expression matrix, with rownames as genes")
  Y <- tibble::as_tibble(bulk_expression, rownames="GENE_SYMBOL")
  
  # remove potential missing data from Y database
  Y <- Y %>% tidyr::drop_na()
  
  # intersect genes (we only keep genes that are common to both data bases)
  common_genes<-intersect(X$GENE_SYMBOL, Y$GENE_SYMBOL)
  
  if (length(common_genes)/dim(X)[1]<0.5) {
    stop (paste("Only", length(common_genes)/dim(X)[1],"fraction of genes are used in the signature matrix\n.
                  Half of common genes are required at least"))}
  X <- X %>% dplyr::filter(GENE_SYMBOL %in% common_genes) %>%
    dplyr::arrange(GENE_SYMBOL) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix()
  Y <- Y %>% dplyr::filter(GENE_SYMBOL %in% common_genes) %>%
    dplyr::arrange(GENE_SYMBOL) %>%
    dplyr::select(where(is.numeric))
  
  # estimation itself
  deconvolution_estimates <- purrr::imap_dfr(deconvolution_functions, function(deconvolution_function, deconvolution_name) {
    additionnal_parameters <- deconvolution_function$additionnal_parameters
    metric_scores <- parallel::mclapply(1:ncol(Y), function(i) {
      # metric_scores <- tibble::tibble(); for (i in 1:ncol(Y)) {
      success_estimation <- tryCatch({
        estimated_p <- do.call(deconvolution_function$FUN, c(list("y"=Y[,i] %>% as.matrix(), "X"=X, "Sigma"=Sigma, "true_ratios"=true_ratios),
                                                             additionnal_parameters))
      },
      error = function(e) {
       # dir.create("/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs", showWarnings = F)
        if (verbose) {
          print(e)
          saveRDS(list("y"=Y[,i] %>% as.matrix(), "X"=X, "Sigma"=Sigma, "error"=e,
                       "function_name" = deconvolution_name, "function_role"=deconvolution_function),
                  file = paste0("/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs/erreur_",
                                i,"_function_", deconvolution_name, ".rds"))
        }
        
        return(e)
      })
      if (!inherits(success_estimation, "error")) {
        return(estimated_p)
      }
    }, mc.cores = cores) %>% dplyr::bind_rows()  # one estimation with a given deconvolution algorithm terminated
    
    if (nrow(metric_scores)!=0) {
      metric_scores <- metric_scores %>% dplyr::mutate(
        OMIC_ID=paste0("sample_", 1:nrow(metric_scores)), deconvolution_name=deconvolution_name)
    }
    return(metric_scores)})
  return (deconvolution_estimates)
}




