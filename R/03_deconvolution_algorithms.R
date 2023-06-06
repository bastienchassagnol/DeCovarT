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


#' @title Mapping function
#'
#' @description The mapping function from unconstrained parameter \eqn{theta},
#' living in \eqn{\mathbb{R}^{J-1}} to parameter vector of the cellular ratios
#' \eqn{p}, subjected to the unit simplex constraint.
#'
#' @param theta The unconstrained parameter, living in \eqn{\mathbb{R}^{J-1}}
#'
#' @return The numeric vector of size \eqn{J},
#' storing the constrained ratios.
#' @export

mapping_function <- function(theta) {
  p <- c(exp(theta[1:length(theta)]), 1)
  return(p / sum(p))
}

# reciprocal of the mapping function
inverse_mapping_function <- function(p) {
  num_cells <- length(p)
  return(log(p[1:num_cells - 1] / p[num_cells]))
}

# compute the mahalabonis distance
.maha_distance <- function(x, A) {
  d <- t(x) %*% solve(A) %*% x # solve A returns the reverted function
  return(d %>% as.numeric()) # supposed to be a scalar
}

# compute the dot product distance, quite similar to maha,
# but without reverting the matrix, and with possibly a different second vector
.dot_product <- function(x, A, y = x) {
  d <- t(x) %*% A %*% y # solve A returns the reverted function
  return(d %>% as.numeric()) # supposed to be a scalar
}

.compute_global_variance <- function(p, Sigma) {
  ###  Sigma and TensorA packages
  # global_cov <- matrix(0, nrow = dim(Sigma)[1], ncol=dim(Sigma)[2])
  # for (j in 1:length(p)) {
  #   global_cov <- global_cov + p[j]^2*Sigma[,,j]
  # }
  #
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
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  log_lik <- -log(det(global_cov_matrix)) - 1 / 2 * .maha_distance(y - X %*% p, global_cov_matrix) %>% as.numeric()
  return(log_lik)
}


loglik_multivariate_constrained <- function(theta, y, X, Sigma) {
  # switch from variable
  p <- mapping_function(theta)
  log_lik <- loglik_multivariate(p, y, X, Sigma)
  if (!check_parameters(p)) { # if the ratios returned present numerical underflows
    warning(paste("Thee ratios are given by", paste(signif(p, digits = 5), collapse = "//"), "and loglik is: ", log_lik))
    # if (length(theta) > 1) {
    #   # more than two components, then we can drop the minimal index
    #   irrelevant_ratios <- which(p < 100 * .Machine$double.eps)
    #   p <- p[-irrelevant_ratios]; X <- X[, -irrelevant_ratios]; Sigma <- Sigma[,,-irrelevant_ratios]
    #   log_lik <- loglik_multivariate(p, y, X, Sigma) # update with modified parameters
    # }
  }
  return(log_lik)
}


# Jacobian mapping function
jacobian_mapping_function <- function(theta) {
  denominator <- (sum(exp(theta)) + 1)^2
  size_var <- length(theta)
  jacobian_matrix <- matrix(0, nrow = size_var, ncol = size_var)
  for (i in 1:size_var) {
    for (j in i:size_var) {
      # diagonal elements
      if (i == j) {
        jacobian_matrix[i, j] <- (exp(theta[i]) * (sum(exp(theta[-i])) + 1))
      } else {
        jacobian_matrix[i, j] <- -exp(theta[i]) * exp(theta[j])
      }
    }
  }
  # ensure symmetry
  jacobian_matrix[lower.tri(jacobian_matrix)] <- jacobian_matrix[upper.tri(jacobian_matrix)]
  jacobian_matrix <- rbind(jacobian_matrix, -exp(theta))
  return(jacobian_matrix / denominator)
}

gradient_loglik_unconstrained <- function(p, y, X, Sigma) {
  # compute general covariance and its reverse
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  global_precision_matrix <- solve(global_cov_matrix)

  # compute the gradient itself
  gradient_unconstrained <- c()
  for (j in 1:length(p)) {
    gradient_unconstrained <- c(
      gradient_unconstrained,
      -2 * p[j] * tr(global_precision_matrix %*% Sigma[, , j]) +
        .dot_product(y - X %*% p, global_precision_matrix, X[, j]) +
        p[j] * .dot_product(y - X %*% p, global_precision_matrix %*% Sigma[, , j] %*% global_precision_matrix)
    )
  }
  return(gradient_unconstrained)
}

gradient_loglik_constrained <- function(theta, y, X, Sigma) {
  p <- mapping_function(theta)
  gradient_constrained <- gradient_loglik_unconstrained(p, y, X, Sigma) %*% jacobian_mapping_function(theta)
  return(gradient_constrained)
}


##################################################################
##                   second order derivatives                   ##
##################################################################
# Hessian mapping function
hessian_mapping_function <- function(theta) {
  A <- sum(exp(theta)) + 1
  denominator <- A^3
  size_var <- length(theta)
  J <- size_var + 1
  hessian_array <- array(0, dim = c(size_var, size_var, J))
  for (i in 1:size_var) {
    B <- sum(exp(theta[-i])) + 1
    # other p_j, j< J Hessian derivation
    for (j in 1:size_var) {
      for (k in j:size_var) {
        # diagonal elements
        if (j == k) {
          if (i == j) {
            hessian_array[i, i, i] <- B * exp(theta[i]) * (B - exp(theta[i])) # condition d)
          } else {
            hessian_array[j, j, i] <- exp(theta[i]) * exp(theta[j]) * (-A + 2 * exp(theta[j])) # condition c)
          }
        }
        # off diagonal terms
        else {
          if (!length(intersect(i, c(j, k)))) { # all indexes are different, situation b)
            # alternative condition setting: (i!=j)!=k
            hessian_array[j, k, i] <- 2 * exp(theta[i]) * exp(theta[j]) * exp(theta[k])
          } else {
            l <- setdiff(c(j, k), i) # situation a), l is the operator, either k or j, different from i
            hessian_array[j, k, i] <- exp(theta[i]) * exp(theta[l]) * (-B + exp(theta[i]))
          }
        }
      }
    }
    # ensure symmetry
    hessian_array[, , i][lower.tri(hessian_array[, , i])] <- hessian_array[, , i][upper.tri(hessian_array[, , i])]
  }
  # last p_J Hessian component
  for (j in 1:size_var) {
    for (k in j:size_var) {
      # diagonal elements
      if (j == k) {
        hessian_array[j, j, J] <- exp(theta[j]) * (-A + 2 * exp(theta[j])) # condition e)
      } else {
        hessian_array[j, k, J] <- 2 * exp(theta[j]) * exp(theta[k]) # condition f)
      }
    }
  }
  hessian_array[, , J][lower.tri(hessian_array[, , J])] <- hessian_array[, , J][upper.tri(hessian_array[, , J])]
  return(hessian_array / denominator)
}

hessian_loglik_unconstrained <- function(p, y, X, Sigma) {
  num_celltypes <- length(p)
  hessian_unconstrained <- matrix(0, nrow = num_celltypes, ncol = num_celltypes)
  global_precision_matrix <- .compute_global_variance(p, Sigma) %>% solve()
  for (i in 1:num_celltypes) {
    for (j in i:num_celltypes) {
      hessian_unconstrained[i, j] <- 4 * p[i] * p[j] *
        tr(global_precision_matrix %*% Sigma[, , i] %*% global_precision_matrix %*% Sigma[, , j]) -
        .dot_product(X[, i], global_precision_matrix, X[, j]) -
        2 * p[i] * .dot_product(
          y - X %*% p,
          global_precision_matrix %*% Sigma[, , i] %*% global_precision_matrix, X[, j]
        ) -
        2 * p[j] * .dot_product(
          y - X %*% p,
          global_precision_matrix %*% Sigma[, , j] %*% global_precision_matrix, X[, i]
        ) -
        4 * p[i] * p[j] * .dot_product(
          y - X %*% p,
          global_precision_matrix %*% Sigma[, , j] %*%
            global_precision_matrix %*% Sigma[, , i] %*%
            global_precision_matrix
        )
      if (i == j) { # add diagonal terms
        hessian_unconstrained[i, i] <- hessian_unconstrained[i, i] -
          2 * tr(global_precision_matrix %*% Sigma[, , i]) +
          .dot_product(y - X %*% p, global_precision_matrix %*% Sigma[, , i] %*% global_precision_matrix)
      }
    }
  }
  # enforce symmetry
  hessian_unconstrained[lower.tri(hessian_unconstrained)] <- hessian_unconstrained[upper.tri(hessian_unconstrained)]
  return(hessian_unconstrained)
}

hessian_loglik_constrained <- function(theta, y, X, Sigma) {
  p <- mapping_function(theta)
  # t(J_psi) X H_log X J_psi + sum over number of ratios of grad_log X H_psi
  hessian_constrained <- t(jacobian_mapping_function(theta)) %*% hessian_loglik_unconstrained(p, y, X, Sigma) %*% jacobian_mapping_function(theta) +
    tensor::tensor(A = gradient_loglik_unconstrained(p, y, X, Sigma), B = hessian_mapping_function(theta), alongA = 1, alongB = 3) %>% as.matrix()
  return(hessian_constrained)
}


#################################################################
##                    iterated descent algorithms              ##
#################################################################


#' Deconvolution algorithm itself, for a given sample.
#'
#' @param y Parameter `y`: \eqn{\boldsymbol{y}=(y_{g}) \in \mathbb{R}^{G}},
#' storing the measured expression of the `G` genes in the heterogeneous sample
#' @param X Parameter `mu`: \eqn{\boldsymbol{\mu}=(\mu_{g,j}) \in \mathbb{R}^{G \times J}},
#' storing in each each column the averaged expression of the `G` genes of the `J` cell populations.
#' @param Sigma Optional, the 3-dimensional covariance matrix array:
#'  \eqn{\mathrm{\Sigma}=(\Sigma_{l, k, j}) \in \mathbb{R}^{G \times G \times J}}, with each matrix
#' \eqn{\Sigma_{..j}, j \in \{ 1, \ldots, J\}} storing the covariance matrix
#' describing the covariance transcriptomic structure of a given cell population \eqn{j}
#' @param true_ratios Optional,  vector of size \eqn{J}, storing the normalised proportions
#' of the cell populations supposed present in the sample. If provided, summary metrics
#' will then be computed against the ones returned by the deconvolution algorithms provided.
#' @param epsilon,itmax Stopping criterion of the deconvolution algorithm,
#' respectively measuring the absolute convergence of the log-likelihood and
#' constraining the maximal number of iterations that the deconvolution algorithm performs
#'
#' @inherit compute_benchmark_metrics return
#' @export


deconvolute_ratios_DeCoVarT <- function(y, X, Sigma, true_ratios = NULL,
                                        epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # set minimize to false; partialH=2
  invisible(utils::capture.output(estimated_theta <- marqLevAlg::marqLevAlg(
    b = initial_theta, fn = loglik_multivariate_constrained,
    gr = gradient_loglik_constrained, hess = hessian_loglik_constrained,
    epsa = epsilon, epsb = epsilon, epsd = epsilon, minimize = FALSE, multipleTry = 1,
    y = y, X = X, Sigma = Sigma, maxiter = itmax
  )$b))
  if (all(is.na(estimated_theta))) { # retrieve the last non missing estimate
    output_lm <- utils::capture.output(marqLevAlg::marqLevAlg(
      b = initial_theta, fn = loglik_multivariate_constrained,
      gr = gradient_loglik_constrained, hess = hessian_loglik_constrained,
      epsa = epsilon, epsb = epsilon, epsd = epsilon, minimize = FALSE, multipleTry = 1,
      y = y, X = X, Sigma = Sigma, maxiter = itmax
    )) # add partialH and blinding?

    estimated_theta <- output_lm[grep("b : ", output_lm, value = F)] %>%
      stringr::str_match_all("[0-9,\\.]+") %>%
      unlist() %>%
      as.numeric() # retrieve last estimate before failure
  }
  estimated_p <- mapping_function(estimated_theta) %>%
    enforce_identifiability() %>% # ensure non-negativity constraint and remove numerical underflow
    stats::setNames(colnames(X))

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT Uses SA (for simulated annealing)  to infer the simulated
#' ratios, see also [stats::optim()] with `method="SANN"` for more details

deconvolute_ratios_simulated_annealing <- function(y, X, Sigma, true_ratios = NULL,
                                                   epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # gr is not used in the simulated annealing approach
  # In SANNN, maxit is the total number of point evaluations, and not the maximum number of iterations
  estimated_theta <- stats::optim(
    par = initial_theta, fn = loglik_multivariate_constrained, y = y, X = X, Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax), method = "SANN"
  )$par
  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()


  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT A variant of the standard BFGS quasi-Newton method,
#' that allows addtional box constraints (here we impose the ratios to be strictly included between 0 and 1)
#' when inferring the simulated ratios, see also [stats::optim()] with `method="L-BFGS-B"` for more details

deconvolute_ratios_LBFGS <- function(y, X, Sigma, true_ratios = NULL,
                                     epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  estimated_p <- stats::optim(
    par = initial_p, fn = loglik_multivariate, gr = gradient_loglik_unconstrained,
    y = y, X = X, Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax, lmm = 1, factr = epsilon * 10), method = "L-BFGS-B",
    lower = rep(0, length(initial_p)), upper = rep(1, length(initial_p))
  )$par %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT An adaptive barrier algorithm enforcing linear inequality constraints.
#' See also [stats::constrOptim()] for more details. Unfortunately, strict equality constraints coupled with inequality boxes
#' are not possible in this method, so we just impose that that the ratios are included between 0 and 1, and
#' that the sum should be inferior to the actual observed global bulk expression.

deconvolute_ratios_constrOptim <- function(y, X, Sigma, true_ratios = NULL,
                                           epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  ui <- diag(nrow = length(initial_p))
  ci <- rep(0, length(initial_p)) # encode the non-negativity constraint

  # encode the sum-to-one constraint, converting inequality into equality
  # ui <- diag(nrow=length(initial_p)-1); ci <- rep(0, length(initial_p))
  # ui <- rbind(ui, rep(1, length(initial_p)), rep(-1, length(initial_p))); ci <- c(ci,1,-1)
  ui <- rbind(ui, rep(1, length(initial_p)))
  ci <- c(ci, 1) # can not make an exact equality


  # add some perturbation to start the iteration
  estimated_p <- stats::constrOptim(
    theta = initial_p + epsilon, f = loglik_multivariate, grad = gradient_loglik_unconstrained,
    ui = ui, ci = ci, control = list(fnscale = -1, maxit = itmax, reltol = epsilon, abstol = epsilon),
    method = "BFGS", outer.iterations = itmax, outer.eps = epsilon, y = y, X = X, Sigma = Sigma
  )$par %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT A standard second-order descent based algorithm,
#' which reveals equivalent to perform a Newton's Raphson algorithm to retrieve the roots of the
#' gradient.  See also [stats::nlminb] for more details.

deconvolute_ratios_second_order <- function(y, X, Sigma, true_ratios = NULL,
                                            epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)


  # with nlm method (use minus to maximise, instead of minimise)
  # f_list <- function(p, y, X, Sigma) {- loglik_multivariate_constrained(p, y, X, Sigma)}  # function(x) -loglik_multivariate_constrained(x)
  # attr(f_list, "gradient") <- function(p, y, X, Sigma) {- gradient_loglik_constrained(p, y, X, Sigma)}
  # attr(f_list, "hessian") <- function(p, y, X, Sigma) {-hessian_loglik_constrained(p, y, X, Sigma)}
  # estimated_theta <- stats::nlm(f=f_list, p=initial_theta, y=y, X=X, Sigma=Sigma, stepmax=1,
  #                                           ndigit=8, gradtol=1e-4, steptol=1e-4, iterlim=200)$estimate

  # with nlmimb package method (outdated, but works well for our scenario)
  estimated_theta <- stats::nlminb(
    start = initial_theta, objective = function(p, y, X, Sigma) {
      -loglik_multivariate_constrained(p, y, X, Sigma)
    },
    gradient = function(p, y, X, Sigma) {
      -gradient_loglik_constrained(p, y, X, Sigma)
    },
    hessian = function(p, y, X, Sigma) {
      -hessian_loglik_constrained(p, y, X, Sigma)
    },
    y = y, X = X, Sigma = Sigma,
    control = list(eval.max = 1, iter.max = itmax, rel.tol = epsilon, x.tol = epsilon, xf.tol = epsilon, abs.tol = epsilon)
  )$par

  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT The most basic optimisation approach possible,
#' implemented as a standard of the expected worst case, and to check whether computing and reparametrising the
#' log-likelihood function was indeed worthy.
#' To do so, we perform simply a basic BFGS descent on the original, non-constrained log-likelihood function,
#' without even providing an explicit formula of the gradient or the hessian. However, we map back the returned estimates
#' to the unit simplex constraint.

deconvolute_ratios_basic_optim <- function(y, X, Sigma, true_ratios = NULL,
                                           epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  # set minimize to false
  estimated_p <- stats::optim(
    par = initial_p, fn = loglik_multivariate, gr = NULL,
    y = y, X = X, Sigma = Sigma, method = "BFGS",
    control = list(fnscale = -1, maxit = itmax, reltol = epsilon, abstol = epsilon)
  )$par %>%
    stats::setNames(colnames(X))
  estimated_p <- enforce_identifiability(estimated_p) # normalize coefficients

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT A standard first-order descent based algorithm,
#' using the BFGS algorithm, see also [stats::optim] with option `method="BFGS`. We provide an explicit formula
#' of the reparametrised log-likelihood function, as well as its gradient.

deconvolute_ratios_first_order <- function(y, X, Sigma, true_ratios = NULL,
                                           epsilon = 10^-4, itmax = 200) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  initial_theta <- inverse_mapping_function(initial_p)

  estimated_theta <- stats::optim(
    par = initial_theta, fn = loglik_multivariate_constrained, gr = gradient_loglik_constrained, y = y, X = X, Sigma = Sigma,
    control = list(fnscale = -1, reltol = epsilon, abstol = epsilon, maxit = itmax), method = "BFGS"
  )$par
  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()


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

#' @describeIn deconvolute_ratios_DeCoVarT From this algorithm, providing
#' any explicit prior of the covariance matrix is pointless, since we are,
#' to our knowledge, the first ones to account for it explicitly. Here,
#' a custom implementation of the CIBERSORT algorithm.

deconvolute_ratios_CIBERSORT <- function(y, X, true_ratios = NULL) {
  # the set of nu values tested for best performance
  range.nu <- seq(0.2, 0.8, 0.3)
  model <- e1071::best.svm(X, y,
    type = "nu-regression", kernel = "linear", scale = F,
    nu = range.nu, tunecontrol = e1071::tune.control(sampling = "fix", fix = 0.75)
  ) # determine best fitted model

  e1071::tune.svm(X, y,
    type = "nu-regression", kernel = "linear", scale = F,
    nu = 0.25, tunecontrol = e1071::tune.control(sampling = "fix", fix = 0.75)
  )

  model <- e1071::svm(X, y = y)

  # get and normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- t(model$coefs) %*% model$SV %>% enforce_identifiability()
  names(estimated_p) <- colnames(X)

  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @importFrom Rdpack reprompt
#' @describeIn deconvolute_ratios_DeCoVarT  Here,
#' a standard linear approach, as performed in \insertCite{abbas_etal09;textual}{DeCovarT} as it can be computed
#' with function [stats::lsfit()]. Nevertheless, similar to any other deconvolution methods,
#' inferred ratios are normalised back to the unit simplex space.

deconvolute_ratios_abbas <- function(y, X, true_ratios = NULL) {
  estimated_p <- stats::lsfit(X, y, intercept = F)$coefficients

  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- estimated_p %>% enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT  Here,
#' a robust linear approach, as performed in \insertCite{monaco_etal19;textual}{DeCovarT}
#' as it can be computed with function [MASS::rlm()].

deconvolute_ratios_monaco <- function(y, X, true_ratios = NULL) {
  estimated_p <- MASS::rlm(y ~ X + 0, method = c("M"))$coefficients
  names(estimated_p) <- colnames(X)

  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- estimated_p %>% enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT Here, function [nnls::nnls()]
#' is used as an interface to the Lawson-Hanson NNLS implementation, with the additional
#' constraints that the raturned ratios can not be negative. It is thus less restrictive than function
#' [deconvolute_ratios_deconRNASeq].

deconvolute_ratios_nnls <- function(y, X, true_ratios = NULL) {
  estimated_p <- nnls::nnls(X, y)$x
  names(estimated_p) <- colnames(X)

  # normalize coefficients (sum to 1, as non-negativity is already enforced)
  estimated_p <- estimated_p %>% enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT  Here,
#' standard linear least squares optimisation, but accounting for
#' the two equality and inequality constraints of the unit simplex. Similar to the implementation
#' of the `deconRNASeq` algorithm, see also [limSolve::lsei()]
#' or the `lsqlin` function in Matlab for additional details.
#' @references
#' \insertAllCited{}

deconvolute_ratios_deconRNASeq <- function(y, X, true_ratios = NULL) {
  EE <- rep(1, ncol(X))
  FF <- 1 # encode the sum-to-one constraint
  GG <- diag(nrow = ncol(X))
  HH <- rep(0, ncol(X)) # encode the non-negativity constraint

  estimated_p <- limSolve::lsei(X, y, EE, FF, GG, HH)$X %>% enforce_identifiability()
  names(estimated_p) <- colnames(X)
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_p, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}



#' @title Compute summary metrics evaluating the quality of the estimate
#'
#' @description Compute metrics, either comparing th estimated ratios with a gold standard
#' or the divergence between the reconstituted virtual mixture,
#' using deterministic rule \eqn{\boldsymbol{\hat{y}}=\boldsymbol{X} \times \hat{\boldsymbol{p}}}
#' and the actual measured one
#'
#' @inheritParams deconvolute_ratios_DeCoVarT
#' @param estimated_p The ratios estimated by your favourite deconvolution algorithm
#'
#' @return A `tibble`, with the following scores:
#' * mse and rmse, for respectively \emph{mean} and \emph{root mean squared error}. See also the [Metrics::mse()] function.
#' * mae, for \emph{mean absolute error}. See also the [Metrics::mae()] function.
#' * \eqn{R^2} and adjusted \eqn{R^2}, corresponding to the percentage of variance
#' captured by the linear regression model. See also the [Metrics::rse()] function.
#' * cor, for the Pearson correlation between the estimated and true cellular ratios
#' giving the mean values of the variables within a given component. See also the [stats::cor()] function.
#' @export

compute_benchmark_metrics <- function(y, X, estimated_p, true_ratios = NULL) {
  n <- nrow(X)
  k <- ncol(X)
  df_res <- n - k + 1 # number of free parameters: only k - 1 parameters must be learnt, with sum-to-one constraint
  df_tot <- n - 1 # no intercept for the moment, in the model (so n-1, or n?)
  if (!is.null(true_ratios)) { # when the true parameters are known
    scores <- tibble::tibble(
      model_mse = Metrics::mse(true_ratios, estimated_p),
      model_rmse = Metrics::rmse(true_ratios, estimated_p),
      model_mae = Metrics::mae(true_ratios, estimated_p),
      model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_p)),
      model_coef_determination_adjusted = max(0, 1 - (1 - .data$model_coef_determination) * df_tot / df_res),
      model_cor = suppressWarnings(stats::cor(true_ratios, estimated_p, method = "pearson"))
    )
  } else { # when they are unknown
    predicted_values <- as.vector(X %*% estimated_p)
    scores <- tibble::tibble(
      model_mse = Metrics::mse(y, predicted_values),
      model_rmse = Metrics::rmse(y, predicted_values),
      model_mae = Metrics::mae(y, predicted_values),
      model_cor = suppressWarnings(stats::cor(y, predicted_values, method = "pearson"))
    )

    # invisible(utils::capture.output
    # model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_p)),
    # model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
    # # model_ccc= epiR::epi.ccc(y, predicted_values)
  }
  return(scores)
}

#' Main function of the package: deconvolute in parallel mixture samples
#'
#' @author Bastien CHASSAGNOL
#'
#' @param signature_matrix Parameter `mu`: \eqn{\boldsymbol{\mu}=(\mu_{g,j}) \in \mathbb{R}^{G \times J}},
#' storing in each each column the averaged expression of the `G` genes used to deconvolve all cell populations.
#' Name your `colnames` as the cell populations, and provide `rownames` argument for the name of the genes.
#' By convention, we use HGNC symbols.
#' @param bulk_expression Parameter `y`: \eqn{\boldsymbol{y}=(\mu_{g,i}) \in \mathbb{R}^{G \times I}},
#' storing in each each column the measured expression of the `G` genes in a heterogeneous sample, using any RNASeq or microarray technology.
#' Provide the sample ID for each of your samples in column, and the name of your genes in `rownames`.
#' @param true_ratios If available (for instance, in the context of a virtual benchmark, or if some standard cytometry techniques provide them),
#' vector of size \eqn{J}, storing the normalised proportions of the cell populations supposed present in the sample. Summary metrics
#' will then be computed against the ones returned by the deconvolution algorithms provided.
#' @param Sigma Only relevant for deconvolution algorithms which require a prior estimate
#' of the transcriptomic  covariance for each of the purified cell populations.
#' A 3-dimensional covariance matrix array is expected:
#'  \eqn{\mathrm{\Sigma}=(\Sigma_{l, k, j}) \in \mathbb{R}^{G \times G \times J}} to
#' parametrise the covariance transcriptomic structure of the \eqn{J} cell populations estimated.
#' @param deconvolution_functions The deconvolution functions themselves, a list
#' with for each item two attributes to be filled with:
#' * `FUN`: the function itself (not a string, but indeed any deconvolution function
#' integrating the default parameters listed in )
#' `additional_parameters` by default, set to NULL. If your deconvolution function
#' integrates any specific, additional parameter.
#' @param scaled Whether we should scale or not the dataset. By default, we consider that the provided dataset is in its original raw space,
#' and we do not scale the dataset, since our deconvolution algorithm assumes a multivariate Gaussian distribution on the raw counts themselves.
#' @param cores For a parallel estimation of ratios in a series of bulk samples,
#' assign a number of cores strictly inferior to the number of cores avalaible
#' on your machine. By default, the maximum, minus in Unix systems, and for
#' OS compatibility, only one on Windows machines
#' @return A `tibble` storing for each row the measured cell proportions, as well as some summary metrics.
#' We ensure for each deconvolution algorithm that the returned estimates respect the unit simplex constraint,
#' with function [enforce_identifiability()].
#' @export
#'
#' @seealso [deconvolute_ratios_DeCoVarT()], to deconvolve a single, already normalised sample

deconvolute_ratios <- function(signature_matrix, bulk_expression,
                               true_ratios = NULL, Sigma = NULL,
                               deconvolution_functions = NULL, scaled = FALSE,
                               cores = ifelse(.Platform$OS.type == "unix", getOption("mc.cores", parallel::detectCores()), 1)) {
  # read in data
  if (!is.matrix(signature_matrix) | is.null(row.names(signature_matrix))) {
    stop("required format for signature is expression matrix, with rownames as genes")
  }
  X <- tibble::as_tibble(signature_matrix, rownames = "GENE_SYMBOL")
  if (!is.matrix(bulk_expression) | is.null(row.names(bulk_expression))) {
    stop("required format for mixture is expression matrix, with rownames as genes")
  }
  Y <- tibble::as_tibble(bulk_expression, rownames = "GENE_SYMBOL")

  # remove potential missing data from Y database
  Y <- Y %>% tidyr::drop_na()

  # intersect genes (we only keep genes that are common to both data bases)
  common_genes <- intersect(X$GENE_SYMBOL, Y$GENE_SYMBOL)

  if (length(common_genes) / dim(X)[1] < 0.5) {
    stop(paste("Only", length(common_genes) / dim(X)[1], "fraction of genes are used in the signature matrix\n.
                  Half of common genes are required at least"))
  }
  X <- X %>%
    dplyr::filter(.data$GENE_SYMBOL %in% common_genes) %>%
    dplyr::arrange(.data$GENE_SYMBOL) %>%
    dplyr::select(dplyr::where(is.numeric)) %>%
    as.matrix()
  Y <- Y %>%
    dplyr::filter(.data$GENE_SYMBOL %in% common_genes) %>%
    dplyr::arrange(.data$GENE_SYMBOL) %>%
    dplyr::select(dplyr::where(is.numeric))

  if (scaled) { # log-2 normalise
    Y <- log2(Y)
  }
  X <- log2(X)

  # estimation itself
  deconvolution_estimates <- purrr::imap_dfr(deconvolution_functions, function(deconvolution_function, algorithm) {
    additional_parameters <- deconvolution_function$additional_parameters
    metric_scores <- parallel::mclapply(1:ncol(Y), function(i) {
      # metric_scores <- tibble::tibble(); for (i in 1:ncol(Y)) {
      success_estimation <- tryCatch(
        {
          list_arguments <- c(list("y" = Y[, i] %>% as.matrix(), "X" = X, "Sigma" = Sigma, "true_ratios" = true_ratios), additional_parameters)
          # we only keep the arguments needed by the required function
          estimated_p <- do.call(deconvolution_function$FUN, list_arguments[methods::formalArgs(deconvolution_function$FUN)])
        },
        error = function(e) {
          # dir.create("/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs", showWarnings = F)
          warning(paste(e, "\n"))
          dir.create("./simulations/erreurs", showWarnings = F, recursive = TRUE)
          # "y"=Y[,i] %>% as.matrix(), "X"=X, "Sigma"=Sigma,
          saveRDS(
            c(
              list_arguments[methods::formalArgs(deconvolution_function$FUN)],
              list(
                "CALL" = call(algorithm, list_arguments[methods::formalArgs(deconvolution_function$FUN)]),
                "error" = e, "function_name" = algorithm
              )
            ),
            file = paste0(
              "/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs/erreur_",
              i, "_function_", algorithm, ".rds"
            )
          )
          return(e)
        }
      )
      if (!inherits(success_estimation, "error")) {
        return(estimated_p)
      }
    }, mc.cores = cores) %>% dplyr::bind_rows() # one estimation with a given deconvolution algorithm terminated

    if (nrow(metric_scores) != 0) {
      metric_scores <- metric_scores %>% dplyr::mutate(
        OMIC_ID = paste0("sample_", 1:nrow(metric_scores)), algorithm = algorithm
      )
    }
    return(metric_scores)
  })
  return(deconvolution_estimates)
}
