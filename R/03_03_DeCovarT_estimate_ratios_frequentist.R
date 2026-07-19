# Algebraic operations ------------------------------------------------------

#' Softmax mapping from unconstrained parameters to the simplex
#'
#' @description
#' Implements the reparametrisation
#' \eqn{\boldsymbol{\psi}:\boldsymbol{\theta}\mapsto\boldsymbol{p}} used in the
#' article, sending unconstrained coordinates
#' \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}} to cellular proportions
#' \eqn{\boldsymbol{p}\in\Delta^{J-1}}.
#'
#' @details
#' With \eqn{A=\sum_{k=1}^{J-1}\mathrm{e}^{\theta_k}+1},
#' \deqn{
#'   p_j=\frac{\mathrm{e}^{\theta_j}}{A}\quad(j<J),\qquad
#'   p_J=\frac{1}{A}.
#' }
#' Equivalently,
#' \eqn{\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\theta})} with
#' \eqn{\boldsymbol{\psi}(\boldsymbol{\theta})\propto
#' (\mathrm{e}^{\theta_1},\ldots,\mathrm{e}^{\theta_{J-1}},1)^{\mathsf{T}}}.
#'
#' @param theta Numeric vector \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @return Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}} on the unit
#'   simplex (\eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1},
#'   \eqn{\boldsymbol{p}\ge\mathbf{0}}).
#'
#' @seealso The inverse map is documented as `inverse_mapping_function()`
#'   on this help page.
#' @export
mapping_function <- function(theta) {
  p <- c(exp(theta[1:length(theta)]), 1)
  return(p / sum(p))
}

#' Inverse simplex mapping \eqn{\boldsymbol{p}\mapsto\boldsymbol{\theta}}
#'
#' Recovers \eqn{\theta_j=\ln(p_j/p_J)} for \eqn{j=1,\ldots,J-1}.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}} on the open simplex.
#'
#' @return Numeric vector \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @rdname mapping_function
#' @keywords internal
inverse_mapping_function <- function(p) {
  num_cells <- length(p)
  return(log(p[1:num_cells - 1] / p[num_cells]))
}

#' @keywords internal
#' @noRd
.maha_distance <- function(mean_signature_matrix, A) {
  d <- t(mean_signature_matrix) %*% solve(A) %*% mean_signature_matrix # solve A returns the reverted function
  return(d |> as.numeric()) # supposed to be a scalar
}

#' @keywords internal
#' @noRd
.dot_product <- function(mean_signature_matrix, A, y = mean_signature_matrix) {
  d <- t(mean_signature_matrix) %*% A %*% y # solve A returns the reverted function
  return(d |> as.numeric()) # supposed to be a scalar
}

#' Bulk mixture covariance \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @description
#' Assembles the conditional covariance of the Gaussian convolution
#' \deqn{
#'   \boldsymbol{\Sigma}(\boldsymbol{p})
#'   =\sum_{j=1}^{J}p_j^{2}\,\boldsymbol{\Sigma}_j,
#' }
#' stored as slices of the array `Sigma`.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param Sigma Array in \eqn{\mathcal{M}_{G\times G\times J}} whose slice
#'   \eqn{\boldsymbol{\Sigma}_j=} `Sigma[,, j]` is the cell-type covariance.
#'
#' @return Symmetric matrix
#'   \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})\in\mathcal{M}_{G\times G}}.
#'
#' @keywords internal
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

# Log-likelihood function, aka the objective function ----------------------

#' Unconstrained DeCovarT log-likelihood
#'
#' @description
#' Evaluates the conditional log-likelihood
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} of a bulk
#' profile under the Gaussian convolution model of the article,
#' \deqn{
#'   \boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})
#'   \sim\mathcal{N}_{G}\!\bigl(\boldsymbol{\mu}\boldsymbol{p},\,
#'   \boldsymbol{\Sigma}(\boldsymbol{p})\bigr),
#' }
#' with plug-in parameters
#' \eqn{\boldsymbol{\zeta}=(\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\}_{j=1}^{J})}
#' and mixture covariance
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_{j}p_j^{2}\boldsymbol{\Sigma}_j}.
#'
#' @details
#' Up to an additive constant independent of \eqn{\boldsymbol{p}},
#' \deqn{
#'   \ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})
#'   =
#'   -\log\det\boldsymbol{\Sigma}(\boldsymbol{p})
#'   -\tfrac{1}{2}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\mathsf{T}}
#'   \boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}).
#' }
#' Argument `mean_signature_matrix` stores \eqn{\boldsymbol{\mu}}.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param y Numeric vector (or one-column matrix)
#'   \eqn{\boldsymbol{y}\in\mathbb{R}^{G}}.
#' @param mean_signature_matrix Numeric matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}.
#' @param Sigma Array of cell-type covariances in
#'   \eqn{\mathcal{M}_{G\times G\times J}}.
#'
#' @return Scalar log-likelihood value.
#'
#' @keywords internal
#' @seealso [gradient_loglik_unconstrained()], [mapping_function()]
loglik_multivariate <- function(p, y, mean_signature_matrix, Sigma) {
  global_cov_matrix <- .compute_global_variance(p, Sigma)
  log_lik <- -log(det(global_cov_matrix)) -
    1 /
      2 *
      .maha_distance(y - mean_signature_matrix %*% p, global_cov_matrix) |>
        as.numeric()
  return(log_lik)
}


#' Constrained log-likelihood
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{\psi}(\boldsymbol{\theta}))}
#'
#' @description
#' Composes [loglik_multivariate()] with [mapping_function()], so that
#' optimisation may be performed over
#' \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @inheritParams loglik_multivariate
#' @param theta Numeric vector \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @return Scalar log-likelihood on the constrained manifold.
#'
#' @keywords internal
loglik_multivariate_constrained <- function(
  theta,
  y,
  mean_signature_matrix,
  Sigma
) {
  # switch from variable
  p <- mapping_function(theta)
  log_lik <- loglik_multivariate(p, y, mean_signature_matrix, Sigma)
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    # if the ratios returned present numerical underflows
    warning(paste(
      "Thee ratios are given by",
      paste(signif(p, digits = 5), collapse = "//"),
      "and loglik is: ",
      log_lik
    ))
  }
  return(log_lik)
}


# First-order derivatives -------------------------------------------------

#' Jacobian \eqn{\mathbf{J}_{\boldsymbol{\psi}}} of the simplex mapping
#'
#' @description
#' Returns
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}\in\mathcal{M}_{J\times(J-1)}} with
#' entries
#' \eqn{(\mathbf{J}_{\boldsymbol{\psi}})_{i,j}=\partial p_i/\partial\theta_j}.
#'
#' @param theta Numeric vector \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @return Numeric matrix of size \eqn{J\times(J-1)}.
#'
#' @keywords internal
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

#' Gradient \eqn{\nabla_{\boldsymbol{p}}\ell} of the unconstrained log-likelihood
#'
#' @description
#' Analytic gradient of [loglik_multivariate()] with respect to
#' \eqn{\boldsymbol{p}}. Writing
#' \eqn{\boldsymbol{\Theta}=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}} and
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}, the
#' \eqn{j}-th coordinate is
#' \deqn{
#'   \frac{\partial\ell}{\partial p_j}
#'   =
#'   -2p_j\,\mathrm{Tr}\!\bigl(\boldsymbol{\Theta}\boldsymbol{\Sigma}_j\bigr)
#'   +\boldsymbol{r}^{\mathsf{T}}\boldsymbol{\Theta}\boldsymbol{\mu}_{\cdot j}
#'   +p_j\,\boldsymbol{r}^{\mathsf{T}}
#'   \boldsymbol{\Theta}\boldsymbol{\Sigma}_j\boldsymbol{\Theta}\boldsymbol{r}.
#' }
#'
#' @details
#' Unit tests compare this analytic gradient to a numerical reference from
#' [numDeriv::grad()] applied to [loglik_multivariate()]. For that check the
#' Richardson method is preferred; main `method.args` knobs:
#' \describe{
#'   \item{\code{eps}}{Initial finite-difference step (default `1e-4`).}
#'   \item{\code{r}}{Number of Richardson extrapolations (default `4`; tests
#'   use `6`). Raising `r` usually improves accuracy more safely than
#'   shrinking `eps` alone.}
#'   \item{\code{d}, \code{v}}{Relative step factor and geometric reduction
#'   between extrapolations (default `v = 2`).}
#'   \item{\code{zero.tol}, \code{show.details}}{See `?numDeriv::grad`.}
#' }
#' Alternative `method` values: `"simple"` and `"complex"`.
#'
#' @inheritParams loglik_multivariate
#'
#' @return Numeric vector
#'   \eqn{\nabla_{\boldsymbol{p}}\ell\in\mathbb{R}^{J}}.
#'
#' @keywords internal
#' @seealso [numDeriv::grad()], [hessian_loglik_unconstrained()]
gradient_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
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
        sum(diag(global_precision_matrix %*% Sigma[,, j])) +
        .dot_product(
          y - mean_signature_matrix %*% p,
          global_precision_matrix,
          mean_signature_matrix[, j]
        ) +
        p[j] *
          .dot_product(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix
          )
    )
  }
  return(gradient_unconstrained)
}


#' Constrained gradient via the chain rule
#'
#' @description
#' Returns
#' \deqn{
#'   \nabla_{\boldsymbol{\theta}}\ell
#'   =
#'   \bigl(\nabla_{\boldsymbol{p}}\ell\bigr)^{\mathsf{T}}
#'   \mathbf{J}_{\boldsymbol{\psi}}(\boldsymbol{\theta}),
#' }
#' i.e. first-order chain rule for
#' \eqn{\ell\circ\boldsymbol{\psi}}.
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Numeric vector in \eqn{\mathbb{R}^{J-1}}.
#'
#' @keywords internal
gradient_loglik_constrained <- function(
  theta,
  y,
  mean_signature_matrix,
  Sigma
) {
  p <- mapping_function(theta)
  gradient_constrained <- gradient_loglik_unconstrained(
    p,
    y,
    mean_signature_matrix,
    Sigma
  ) %*%
    jacobian_mapping_function(theta)
  return(gradient_constrained)
}


# Second-order derivatives ------------------------------------------------

#' Second derivatives of \eqn{\boldsymbol{\psi}}
#'
#' @description
#' Tensor of mixed partials
#' \eqn{\partial^{2}p_i/(\partial\theta_k\partial\theta_j)}, stored as an
#' array of size \eqn{(J-1)\times(J-1)\times J}.
#'
#' @param theta Numeric vector \eqn{\boldsymbol{\theta}\in\mathbb{R}^{J-1}}.
#'
#' @return Numeric array used in the constrained Hessian chain rule.
#'
#' @keywords internal
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

#' Hessian \eqn{\mathbf{H}} of the unconstrained log-likelihood
#'
#' @description
#' Analytic Hessian
#' \eqn{\mathbf{H}\in\mathcal{M}_{J\times J}} with entries
#' \eqn{\mathbf{H}_{i,j}=\partial^{2}\ell/(\partial p_i\partial p_j)},
#' matching the matricial formulae of the article (quadratic forms in
#' \eqn{\boldsymbol{\Theta}}, \eqn{\boldsymbol{\Sigma}_i},
#' \eqn{\boldsymbol{\mu}_{\cdot i}} and residual
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}).
#'
#' @inheritParams loglik_multivariate
#'
#' @return Symmetric numeric matrix \eqn{\mathbf{H}}.
#'
#' @keywords internal
hessian_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
  num_celltypes <- length(p)
  hessian_unconstrained <- matrix(0, nrow = num_celltypes, ncol = num_celltypes)
  global_precision_matrix <- .compute_global_variance(p, Sigma) |> solve()
  for (i in 1:num_celltypes) {
    for (j in i:num_celltypes) {
      hessian_unconstrained[i, j] <- 4 *
        p[i] *
        p[j] *
        sum(diag(
          global_precision_matrix %*%
            Sigma[,, i] %*%
            global_precision_matrix %*%
            Sigma[,, j]
        )) -
        .dot_product(
          mean_signature_matrix[, i],
          global_precision_matrix,
          mean_signature_matrix[, j]
        ) -
        2 *
          p[i] *
          .dot_product(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, i] %*% global_precision_matrix,
            mean_signature_matrix[, j]
          ) -
        2 *
          p[j] *
          .dot_product(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix,
            mean_signature_matrix[, i]
          ) -
        4 *
          p[i] *
          p[j] *
          .dot_product(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*%
              Sigma[,, j] %*%
              global_precision_matrix %*%
              Sigma[,, i] %*%
              global_precision_matrix
          )
      if (i == j) {
        # add diagonal terms
        hessian_unconstrained[i, i] <- hessian_unconstrained[i, i] -
          2 * sum(diag(global_precision_matrix %*% Sigma[,, i])) +
          .dot_product(
            y - mean_signature_matrix %*% p,
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


#' Constrained Hessian of \eqn{\ell\circ\boldsymbol{\psi}}
#'
#' @description
#' Second-order chain rule
#' \deqn{
#'   \mathbf{H}_{\boldsymbol{\theta}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   \mathbf{H}_{\boldsymbol{p}}
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   +\sum_{i=1}^{J}
#'   \frac{\partial\ell}{\partial p_i}\,
#'   \frac{\partial^{2}p_i}{\partial\boldsymbol{\theta}\partial\boldsymbol{\theta}^{\mathsf{T}}}.
#' }
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Symmetric matrix in \eqn{\mathcal{M}_{(J-1)\times(J-1)}}.
#'
#' @keywords internal
hessian_loglik_constrained <- function(theta, y, mean_signature_matrix, Sigma) {
  p <- mapping_function(theta)
  # t(J_psi) mean_signature_matrix H_log mean_signature_matrix J_psi + sum over number of ratios of grad_log mean_signature_matrix H_psi
  hessian_constrained <- t(jacobian_mapping_function(theta)) %*%
    hessian_loglik_unconstrained(p, y, mean_signature_matrix, Sigma) %*%
    jacobian_mapping_function(theta) +
    tensor::tensor(
      A = gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
      B = hessian_mapping_function(theta),
      alongA = 1,
      alongB = 3
    ) |>
      as.matrix()
  return(hessian_constrained)
}


# DeCovarT core optimisation algorithms -----------------------------------

#' DeCovarT MLE of cellular proportions for one bulk sample
#'
#' @description
#' Estimates
#' \eqn{\hat{\boldsymbol{p}}=\arg\max_{\boldsymbol{p}}
#' \ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} under the
#' Gaussian convolution model
#' \deqn{
#'   \boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})
#'   \sim\mathcal{N}_{G}\!\bigl(\boldsymbol{\mu}\boldsymbol{p},\,
#'   \boldsymbol{\Sigma}(\boldsymbol{p})\bigr),
#'   \qquad
#'   \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_{j=1}^{J}p_j^{2}\boldsymbol{\Sigma}_j,
#' }
#' subject to the simplex constraint
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1}, \eqn{\boldsymbol{p}\ge\mathbf{0}}.
#' Optimisation is performed in unconstrained coordinates
#' \eqn{\boldsymbol{\theta}} via \eqn{\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\theta})}
#' (Marquardt–Levenberg default; see other methods below).
#'
#' @param y Bulk expression vector
#'   \eqn{\boldsymbol{y}\in\mathbb{R}^{G}} (one heterogeneous sample).
#' @param mean_signature_matrix Mean signature
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (columns = cell types).
#' @param Sigma Array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j=1}^{J}\in\mathcal{M}_{G\times G\times J}}
#'   of cell-type covariances.
#' @param true_ratios Optional ground-truth
#'   \eqn{\boldsymbol{p}^{\star}\in\mathbb{R}^{J}} used only for
#'   benchmark scores.
#' @param epsilon,itmax Absolute convergence tolerance and maximum number of
#'   iterations for the optimiser.
#'
#' @inherit compute_benchmark_metrics return
#' @export
#' @seealso [deconvolute_ratios()], [mapping_function()]
deconvolute_ratios_Marquardt_Levenberg <- function(
  y,
  mean_signature_matrix,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(mean_signature_matrix), ncol(mean_signature_matrix)) # consider by hypothesis equi-balanced proportions between cell populations
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
      mean_signature_matrix = mean_signature_matrix,
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
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      maxiter = itmax
    )) # add partialH and blinding?

    estimated_theta <- output_lm[grep("b : ", output_lm, value = F)] |>
      stringr::str_match_all("[0-9,\\.]+") |>
      unlist() |>
      as.numeric() # retrieve last estimate before failure
  }
  estimated_p <- mapping_function(estimated_theta) |>
    enforce_identifiability() |> # ensure non-negativity constraint and remove numerical underflow
    stats::setNames(colnames(mean_signature_matrix))

  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_Marquardt_Levenberg Simulated annealing on
#'   \eqn{\boldsymbol{\theta}} ([stats::optim()] with `method = "SANN"`).

deconvolute_ratios_simulated_annealing <- function(
  y,
  mean_signature_matrix,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(mean_signature_matrix), ncol(mean_signature_matrix)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)
  # gr is not used in the simulated annealing approach
  # In SANNN, maxit is the total number of point evaluations, and not the maximum number of iterations
  estimated_theta <- stats::optim(
    par = initial_theta,
    fn = loglik_multivariate_constrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax),
    method = "SANN"
  )$par
  estimated_p <- mapping_function(estimated_theta) |>
    stats::setNames(colnames(mean_signature_matrix)) |>
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Box-constrained L-BFGS-B
#'   directly in \eqn{\boldsymbol{p}} ([stats::optim()] `method = "L-BFGS-B"`).

deconvolute_ratios_L_BFGS_B <- function(
  y,
  mean_signature_matrix,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(mean_signature_matrix), ncol(mean_signature_matrix)) # consider by hypothesis equi-balanced proportions between cell populations
  estimated_p <- stats::optim(
    par = initial_p,
    fn = loglik_multivariate,
    gr = gradient_loglik_unconstrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax, lmm = 1, factr = epsilon * 10),
    method = "L-BFGS-B",
    lower = rep(0, length(initial_p)),
    upper = rep(1, length(initial_p))
  )$par |>
    stats::setNames(colnames(mean_signature_matrix)) |>
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Newton–Raphson /
#'   `nlminb` on \eqn{\boldsymbol{\theta}} using analytic gradient and Hessian
#'   ([stats::nlminb()]).

deconvolute_ratios_Newton_Raphson <- function(
  y,
  mean_signature_matrix,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(mean_signature_matrix), ncol(mean_signature_matrix)) # consider by hypothesis equi-balanced proportions between cell populations
  initial_theta <- inverse_mapping_function(initial_p)

  # with nlmimb package method (outdated, but works well for our scenario)
  estimated_theta <- stats::nlminb(
    start = initial_theta,
    objective = function(p, y, mean_signature_matrix, Sigma) {
      -loglik_multivariate_constrained(p, y, mean_signature_matrix, Sigma)
    },
    gradient = function(p, y, mean_signature_matrix, Sigma) {
      -gradient_loglik_constrained(p, y, mean_signature_matrix, Sigma)
    },
    hessian = function(p, y, mean_signature_matrix, Sigma) {
      -hessian_loglik_constrained(p, y, mean_signature_matrix, Sigma)
    },
    y = y,
    mean_signature_matrix = mean_signature_matrix,
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

  estimated_p <- mapping_function(estimated_theta) |>
    stats::setNames(colnames(mean_signature_matrix)) |>
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg BFGS quasi-Newton ascent
#'   on \eqn{\boldsymbol{\theta}} ([stats::optim()] `method = "BFGS"`).

deconvolute_ratios_gradient_descent <- function(
  y,
  mean_signature_matrix,
  Sigma,
  true_ratios = NULL,
  epsilon = 10^-4,
  itmax = 200
) {
  initial_p <- rep(1 / ncol(mean_signature_matrix), ncol(mean_signature_matrix)) # consider by hypothesis equi-balanced proportions between cell populations

  initial_theta <- inverse_mapping_function(initial_p)

  estimated_theta <- stats::optim(
    par = initial_theta,
    fn = loglik_multivariate_constrained,
    gr = gradient_loglik_constrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(
      fnscale = -1,
      reltol = epsilon,
      abstol = epsilon,
      maxit = itmax
    ),
    method = "BFGS"
  )$par
  estimated_p <- mapping_function(estimated_theta) |>
    stats::setNames(colnames(mean_signature_matrix)) |>
    enforce_identifiability()

  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}
