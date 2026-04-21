
# Algebraic operations ------------------------------------------------------

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

# compute a dot product
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

  return(global_cov)
}


# First-order derivatives -------------------------------------------------



# log-likelihood multivariate function
loglik_multivariate <- function(p, y, X, Sigma) {
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  log_lik <- -log(det(global_cov_matrix)) -
    1 / 2 * .maha_distance(y - X %*% p, global_cov_matrix) %>% as.numeric()
  return(log_lik)
}


loglik_multivariate_constrained <- function(theta, y, X, Sigma) {
  # switch from variable
  p <- mapping_function(theta)
  log_lik <- loglik_multivariate(p, y, X, Sigma)
  if (!check_parameters(p)) {
    # if the ratios returned present numerical underflows
    warning(paste(
      "Thee ratios are given by",
      paste(signif(p, digits = 5), collapse = "//"),
      "and loglik is: ",
      log_lik
    ))
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
  jacobian_matrix[lower.tri(jacobian_matrix)] <- jacobian_matrix[upper.tri(
    jacobian_matrix
  )]
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
      -2 *
        p[j] *
        tr(global_precision_matrix %*% Sigma[,, j]) +
        .dot_product(y - X %*% p, global_precision_matrix, X[, j]) +
        p[j] *
          .dot_product(
            y - X %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix
          )
    )
  }
  return(gradient_unconstrained)
}


## Actual first-order likelihood cancelled out -------------------------------------------

gradient_loglik_constrained <- function(theta, y, X, Sigma) {
  p <- mapping_function(theta)
  gradient_constrained <- gradient_loglik_unconstrained(p, y, X, Sigma) %*%
    jacobian_mapping_function(theta)
  return(gradient_constrained)
}



# Second-order derivatives ------------------------------------------------


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
            hessian_array[j, j, i] <- exp(theta[i]) *
              exp(theta[j]) *
              (-A + 2 * exp(theta[j])) # condition c)
          }
        } else {
          # off diagonal terms
          if (!length(intersect(i, c(j, k)))) {
            # all indexes are different, situation b)
            # alternative condition setting: (i!=j)!=k
            hessian_array[j, k, i] <- 2 *
              exp(theta[i]) *
              exp(theta[j]) *
              exp(theta[k])
          } else {
            l <- setdiff(c(j, k), i) # situation a), l is the operator, either k or j, different from i
            hessian_array[j, k, i] <- exp(theta[i]) *
              exp(theta[l]) *
              (-B + exp(theta[i]))
          }
        }
      }
    }
    # ensure symmetry
    hessian_array[,, i][lower.tri(hessian_array[,, i])] <- hessian_array[,,
      i
    ][upper.tri(hessian_array[,, i])]
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
  hessian_array[,, J][lower.tri(hessian_array[,, J])] <- hessian_array[,,
    J
  ][upper.tri(hessian_array[,, J])]
  return(hessian_array / denominator)
}

hessian_loglik_unconstrained <- function(p, y, X, Sigma) {
  num_celltypes <- length(p)
  hessian_unconstrained <- matrix(0, nrow = num_celltypes, ncol = num_celltypes)
  global_precision_matrix <- .compute_global_variance(p, Sigma) %>% solve()
  for (i in 1:num_celltypes) {
    for (j in i:num_celltypes) {
      hessian_unconstrained[i, j] <- 4 *
        p[i] *
        p[j] *
        tr(
          global_precision_matrix %*%
            Sigma[,, i] %*%
            global_precision_matrix %*%
            Sigma[,, j]
        ) -
        .dot_product(X[, i], global_precision_matrix, X[, j]) -
        2 *
          p[i] *
          .dot_product(
            y - X %*% p,
            global_precision_matrix %*% Sigma[,, i] %*% global_precision_matrix,
            X[, j]
          ) -
        2 *
          p[j] *
          .dot_product(
            y - X %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix,
            X[, i]
          ) -
        4 *
          p[i] *
          p[j] *
          .dot_product(
            y - X %*% p,
            global_precision_matrix %*%
              Sigma[,, j] %*%
              global_precision_matrix %*%
              Sigma[,, i] %*%
              global_precision_matrix
          )
      if (i == j) {
        # add diagonal terms
        hessian_unconstrained[i, i] <- hessian_unconstrained[i, i] -
          2 * tr(global_precision_matrix %*% Sigma[,, i]) +
          .dot_product(
            y - X %*% p,
            global_precision_matrix %*% Sigma[,, i] %*% global_precision_matrix
          )
      }
    }
  }
  # enforce symmetry
  hessian_unconstrained[lower.tri(
    hessian_unconstrained
  )] <- hessian_unconstrained[upper.tri(hessian_unconstrained)]
  return(hessian_unconstrained)
}


## Actual second-order function optimised ------------------------------------------------

hessian_loglik_constrained <- function(theta, y, X, Sigma) {
  p <- mapping_function(theta)
  # t(J_psi) X H_log X J_psi + sum over number of ratios of grad_log X H_psi
  hessian_constrained <- t(jacobian_mapping_function(theta)) %*%
    hessian_loglik_unconstrained(p, y, X, Sigma) %*%
    jacobian_mapping_function(theta) +
    tensor::tensor(
      A = gradient_loglik_unconstrained(p, y, X, Sigma),
      B = hessian_mapping_function(theta),
      alongA = 1,
      alongB = 3
    ) %>%
      as.matrix()
  return(hessian_constrained)
}



# DeCovarT core optimisation algorithms -----------------------------------


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

deconvolute_ratios_DeCoVarT <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # set minimize to false; partialH=2
  invisible(utils::capture.output(
    estimated_theta <- marqLevAlg::marqLevAlg(
      b = initial_theta,
      fn = loglik_multivariate_constrained,
      gr = gradient_loglik_constrained,
      hess = hessian_loglik_constrained,
      epsa = epsilon,
      epsb = epsilon,
      epsd = epsilon,
      minimize = FALSE,
      multipleTry = 1,
      y = y,
      X = X,
      Sigma = Sigma,
      maxiter = itmax
    )$b
  ))
  if (all(is.na(estimated_theta))) {
    # retrieve the last non missing estimate
    output_lm <- utils::capture.output(marqLevAlg::marqLevAlg(
      b = initial_theta,
      fn = loglik_multivariate_constrained,
      gr = gradient_loglik_constrained,
      hess = hessian_loglik_constrained,
      epsa = epsilon,
      epsb = epsilon,
      epsd = epsilon,
      minimize = FALSE,
      multipleTry = 1,
      y = y,
      X = X,
      Sigma = Sigma,
      maxiter = itmax
    )) # add partialH and blinding?

    estimated_theta <- output_lm[grep("b : ", output_lm, value = F)] %>%
      stringr::str_match_all("[0-9,\\.]+") %>%
      unlist() %>%
      as.numeric() # retrieve last estimate before failure
  }
  estimated_p <- mapping_function(estimated_theta) %>%
    enforce_identifiability() %>% # ensure non-negativity constraint and remove numerical underflow
    stats::setNames(colnames(X))

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT Uses SA (for simulated annealing)  to infer the simulated
#' ratios, see also [stats::optim()] with `method="SANN"` for more details

deconvolute_ratios_simulated_annealing <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # gr is not used in the simulated annealing approach
  # In SANNN, maxit is the total number of point evaluations, and not the maximum number of iterations
  estimated_theta <- stats::optim(
    par = initial_theta,
    fn = loglik_multivariate_constrained,
    y = y,
    X = X,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax),
    method = "SANN"
  )$par
  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT A variant of the standard BFGS quasi-Newton method,
#' that allows addtional box constraints (here we impose the ratios to be strictly included between 0 and 1)
#' when inferring the simulated ratios, see also [stats::optim()] with `method="L-BFGS-B"` for more details

deconvolute_ratios_LBFGS <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  estimated_p <- stats::optim(
    par = initial_p,
    fn = loglik_multivariate,
    gr = gradient_loglik_unconstrained,
    y = y,
    X = X,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax, lmm = 1, factr = epsilon * 10),
    method = "L-BFGS-B",
    lower = rep(0, length(initial_p)),
    upper = rep(1, length(initial_p))
  )$par %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT An adaptive barrier algorithm enforcing linear inequality constraints.
#' See also [stats::constrOptim()] for more details. Unfortunately, strict equality constraints coupled with inequality boxes
#' are not possible in this method, so we just impose that that the ratios are included between 0 and 1, and
#' that the sum should be inferior to the actual observed global bulk expression.

deconvolute_ratios_constrOptim <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
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
    theta = initial_p + epsilon,
    f = loglik_multivariate,
    grad = gradient_loglik_unconstrained,
    ui = ui,
    ci = ci,
    control = list(
      fnscale = -1,
      maxit = itmax,
      reltol = epsilon,
      abstol = epsilon
    ),
    method = "BFGS",
    outer.iterations = itmax,
    outer.eps = epsilon,
    y = y,
    X = X,
    Sigma = Sigma
  )$par %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_DeCoVarT A standard second-order descent based algorithm,
#' which reveals equivalent to perform a Newton's Raphson algorithm to retrieve the roots of the
#' gradient.  See also [stats::nlminb] for more details.

deconvolute_ratios_Newton_Raphson <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)

  # with nlmimb package method (outdated, but works well for our scenario)
  estimated_theta <- stats::nlminb(
    start = initial_theta,
    objective = function(p, y, X, Sigma) {
      -loglik_multivariate_constrained(p, y, X, Sigma)
    },
    gradient = function(p, y, X, Sigma) {
      -gradient_loglik_constrained(p, y, X, Sigma)
    },
    hessian = function(p, y, X, Sigma) {
      -hessian_loglik_constrained(p, y, X, Sigma)
    },
    y = y,
    X = X,
    Sigma = Sigma,
    control = list(
      eval.max = 1,
      iter.max = itmax,
      rel.tol = epsilon,
      x.tol = epsilon,
      xf.tol = epsilon,
      abs.tol = epsilon
    )
  )$par

  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_DeCoVarT A standard first-order descent based algorithm,
#' using the BFGS algorithm, see also [stats::optim] with option `method="BFGS`. We provide an explicit formula
#' of the reparametrised log-likelihood function, as well as its gradient.

deconvolute_ratios_first_order <- function(
  y,
  X,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(X), ncol(X)) # consider by hypothesis equi-balanced proportions between cell populations

  initial_theta <- inverse_mapping_function(initial_p)

  estimated_theta <- stats::optim(
    par = initial_theta,
    fn = loglik_multivariate_constrained,
    gr = gradient_loglik_constrained,
    y = y,
    X = X,
    Sigma = Sigma,
    control = list(
      fnscale = -1,
      reltol = epsilon,
      abstol = epsilon,
      maxit = itmax
    ),
    method = "BFGS"
  )$par
  estimated_p <- mapping_function(estimated_theta) %>%
    stats::setNames(colnames(X)) %>%
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    X,
    estimated_p,
    true_ratios
  ) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

