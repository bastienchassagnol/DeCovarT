# Population-level convolution: Sigma(p) = sum_j p_j^2 Sigma_j.
# The single-cell model is linear in p_j and lives in
# R/02_02_single_cell_decovart_estimate_ratios.R.

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
#' @inheritParams loglik_multivariate
#'
#' @return Symmetric matrix
#'   \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})\in\mathcal{M}_{G\times G}}.
#'
#' @keywords internal
#' @examples
#' p <- c(0.6, 0.4)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' .compute_global_variance(p, Sigma)
#' @export
.compute_global_variance <- function(p, Sigma) {
  .mixture_covariance(p^2, Sigma)
}

#' Cached Cholesky factorisation of \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @description
#' Assembles [.compute_global_variance()] and factorises it once via
#' [base::chol()], returning the matrix, its upper Cholesky factor, the
#' log-determinant (from the factor's diagonal) and the precision
#' \eqn{\boldsymbol{\Omega}(\boldsymbol{p})} (from the same factor, via
#' [base::chol2inv()]). The covariance is treated as dense: no sparsity or
#' low-rank structure is assumed.
#'
#' @details
#' Newton-type solvers evaluate the log-likelihood, score and Hessian at the
#' same trial \eqn{\boldsymbol{p}} within one iteration. The shared
#' single-slot cache (keyed on \eqn{p_j^2} and `Sigma`) factorises
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} once (\eqn{O(G^{3})}) and
#' returns the cached result to the later calls.
#'
#' @inheritParams loglik_multivariate
#'
#' @return A list with `matrix` (\eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}),
#'   `chol` (upper-triangular Cholesky factor), `log_det`
#'   (\eqn{\log\det\boldsymbol{\Sigma}(\boldsymbol{p})}) and `inverse`
#'   (\eqn{\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}}).
#'
#' @keywords internal
#' @examples
#' p <- c(0.6, 0.4)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' names(.sigma_p_factorisation(p, Sigma))
#' @export
.sigma_p_factorisation <- function(p, Sigma) {
  .mixture_factorisation(p^2, Sigma)
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
#'   -\tfrac{1}{2}\log\det\boldsymbol{\Sigma}(\boldsymbol{p})
#'   -\tfrac{1}{2}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\mathsf{T}}
#'   \boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}).
#' }
#'
#' Computationally this is the same Cholesky-and-backsolve evaluation as
#' `mvtnorm::dmvnorm(..., log = TRUE)` (Genz and Bretz), omitting only the
#' additive \eqn{-\tfrac{G}{2}\log(2\pi)} that does not depend on
#' \eqn{\boldsymbol{p}}. The cached factor from
#' [.sigma_p_factorisation()] supplies the upper-triangular
#' \eqn{\boldsymbol{R}} with
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\boldsymbol{R}^{\mathsf{T}}\boldsymbol{R}};
#' the Mahalanobis term is then
#' \eqn{\lVert\boldsymbol{R}^{-\mathsf{T}}
#' (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})\rVert^{2}},
#' obtained by [base::backsolve()] without forming the explicit inverse.
#' The inverse is still cached because the analytic score and Hessian
#' need \eqn{\boldsymbol{\Omega}(\boldsymbol{p})
#' =\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}}. A QR factorisation of
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} would be a more expensive route
#' to the same SPD quantities; Cholesky is the natural factorisation.
#'
#' Argument `mean_signature_matrix` stores the plug-in mean signature
#' \eqn{\boldsymbol{\mu}}. Latent sample-specific profiles
#' \eqn{\boldsymbol{x}_{\cdot j}} are **not** observed; the frequentist
#' likelihood treats \eqn{\boldsymbol{\mu}} as a fixed proxy. Estimating those
#' latents jointly with \eqn{\boldsymbol{p}} requires a Bayesian / MAP step
#' (see `.map_gaussian_convolution()`).
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param y Numeric vector (or one-column matrix)
#'   \eqn{\boldsymbol{y}\in\mathbb{R}^{G}}.
#' @param mean_signature_matrix Numeric matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (plug-in means).
#' @param Sigma Array of cell-type covariances in
#'   \eqn{\mathcal{M}_{G\times G\times J}}.
#'
#' @return Scalar log-likelihood value.
#'
#' @keywords internal
#' @export
#' @examples
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' loglik_multivariate(p, y, mu, Sigma)
#' @seealso [gradient_loglik_unconstrained()], [isometric_logistic()]
loglik_multivariate <- function(p, y, mean_signature_matrix, Sigma) {
  sigma_p <- .sigma_p_factorisation(p, Sigma)
  residual <- as.numeric(y - drop(mean_signature_matrix %*% p))
  # Mahalanobis distance via the cached Cholesky factor, matching
  # mvtnorm::dmvnorm: solve R^T z = r, then r^T Sigma^{-1} r = ||z||^2.
  z <- backsolve(sigma_p$chol, residual, transpose = TRUE)
  -1 / 2 * sigma_p$log_det - 1 / 2 * sum(z * z)
}


#' Constrained log-likelihood
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{\psi}(\boldsymbol{z}))}
#'
#' @description
#' Composes [loglik_multivariate()] with [isometric_logistic()], so that
#' optimisation may be performed over unconstrained ILR coordinates
#' \eqn{\boldsymbol{z}\in\mathbb{R}^{J-1}}.
#'
#' @inheritParams loglik_multivariate
#' @param z Numeric vector \eqn{\boldsymbol{z}\in\mathbb{R}^{J-1}}.
#' @param V Optional ILR basis; see [isometric_logistic()].
#'
#' @return Scalar log-likelihood on the constrained manifold.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' loglik_multivariate_constrained(isometric_log_ratio(p), y, mu, Sigma)
#' @export
loglik_multivariate_constrained <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  log_lik <- loglik_multivariate(p, y, mean_signature_matrix, Sigma)
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    warning(
      "The ratios are given by ",
      paste(signif(p, digits = 5), collapse = "//"),
      " and loglik is: ",
      log_lik,
      call. = FALSE
    )
  }
  return(log_lik)
}


# First-order derivatives -------------------------------------------------

#' Gradient \eqn{\nabla_{\boldsymbol{p}}\ell} of the unconstrained log-likelihood
#'
#' @description
#' Analytic gradient of [loglik_multivariate()] with respect to
#' \eqn{\boldsymbol{p}}. Writing
#' \eqn{\boldsymbol{\Omega}=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}} and
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}, the
#' \eqn{j}-th coordinate is
#' \deqn{
#'   \frac{\partial\ell}{\partial p_j}
#'   =
#'   -p_j\,\mathrm{Tr}\!\bigl(\boldsymbol{\Omega}\boldsymbol{\Sigma}_j\bigr)
#'   +\boldsymbol{r}^{\mathsf{T}}\boldsymbol{\Omega}\boldsymbol{\mu}_{\cdot j}
#'   +p_j\,\boldsymbol{r}^{\mathsf{T}}
#'   \boldsymbol{\Omega}\boldsymbol{\Sigma}_j\boldsymbol{\Omega}\boldsymbol{r}.
#' }
#' The determinant score is \eqn{-p_j\mathrm{Tr}(\boldsymbol{\Omega}
#' \boldsymbol{\Sigma}_j)} because
#' \eqn{\partial\boldsymbol{\Sigma}/\partial p_j=2p_j\boldsymbol{\Sigma}_j}
#' enters \eqn{-\tfrac{1}{2}\log\det\boldsymbol{\Sigma}(\boldsymbol{p})}; the
#' two residual terms are unaffected by that factor.
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
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' gradient_loglik_unconstrained(p, y, mu, Sigma)
#' @export
#' @seealso [numDeriv::grad()], [hessian_loglik_unconstrained()]
gradient_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
  omega <- .sigma_p_factorisation(p, Sigma)$inverse
  # w = Omega r: every term of the score is a product with it.
  w <- drop(omega %*% (y - mean_signature_matrix %*% p))
  # Sigma_j and Omega are symmetric, so Tr(Omega Sigma_j) = sum(Omega * Sigma_j).
  trace_term <- vapply(
    seq_along(p),
    function(j) sum(omega * Sigma[,, j]),
    numeric(1)
  )
  residual_term <- vapply(
    seq_along(p),
    function(j) sum(w * (Sigma[,, j] %*% w)),
    numeric(1)
  )
  as.numeric(
    -p * trace_term + crossprod(mean_signature_matrix, w) + p * residual_term
  )
}


#' Constrained gradient via the chain rule
#'
#' @description
#' Returns
#' \deqn{
#'   \nabla_{\boldsymbol{z}}\ell
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   \nabla_{\boldsymbol{p}}\ell
#'   =
#'   \mathbf{V}^{\mathsf{T}}
#'   \mathbf{S}(\boldsymbol{p})
#'   \nabla_{\boldsymbol{p}}\ell,
#' }
#' i.e. the first-order chain rule for
#' \eqn{\ell\circ\boldsymbol{\psi}} in ILR coordinates.
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Numeric vector in \eqn{\mathbb{R}^{J-1}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' gradient_loglik_constrained(isometric_log_ratio(p), y, mu, Sigma)
#' @export
gradient_loglik_constrained <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  .ilr_gradient(
    z,
    gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
    V
  )
}


# Second-order derivatives ------------------------------------------------

#' Hessian \eqn{\mathbf{H}} of the unconstrained log-likelihood
#'
#' @description
#' Analytic Hessian
#' \eqn{\mathbf{H}\in\mathcal{M}_{J\times J}} with entries
#' \eqn{\mathbf{H}_{i,j}=\partial^{2}\ell/(\partial p_i\partial p_j)},
#' matching the matrix formulae of the article (quadratic forms in
#' \eqn{\boldsymbol{\Omega}}, \eqn{\boldsymbol{\Sigma}_i},
#' \eqn{\boldsymbol{\mu}_{\cdot i}} and residual
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}).
#'
#' @details
#' The log-determinant contributes
#' \eqn{-\delta_{ij}\mathrm{Tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}_j)
#' +2p_ip_j\mathrm{Tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}_i
#' \boldsymbol{\Omega}\boldsymbol{\Sigma}_j)}, i.e. half the coefficients of
#' the pre-2.3.0 objective, which used
#' \eqn{-\log\det\boldsymbol{\Sigma}(\boldsymbol{p})}. Residual terms are
#' unchanged. Taking expectations under
#' \eqn{\boldsymbol{r}\sim\mathcal{N}_G(\boldsymbol{0},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))} cancels the determinant and
#' residual traces, leaving exactly
#' \eqn{\mathbb{E}[-\mathbf{H}]=I(\boldsymbol{p})} of
#' [expected_fisher_unconstrained()].
#'
#' @inheritParams loglik_multivariate
#'
#' @return Symmetric numeric matrix \eqn{\mathbf{H}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' hessian_loglik_unconstrained(p, y, mu, Sigma)
#' @export
hessian_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
  blocks <- .information_blocks(p, mean_signature_matrix, Sigma)
  w <- drop(blocks$precision %*% (y - mean_signature_matrix %*% p))
  # Column j of s is Sigma_j Omega r, column j of u is Omega Sigma_j Omega r.
  s <- matrix(
    vapply(
      seq_along(p),
      function(j) drop(Sigma[,, j] %*% w),
      numeric(length(w))
    ),
    nrow = length(w)
  )
  u <- blocks$precision %*% s
  # cross[i, j] = p_i r' Omega Sigma_i Omega mu_j
  cross <- p * crossprod(u, mean_signature_matrix)
  hessian <- blocks$cov -
    blocks$mean -
    2 * (cross + t(cross)) -
    4 * tcrossprod(p) * crossprod(s, u)
  diag(hessian) <- diag(hessian) - blocks$traces + colSums(s * w)
  unname(hessian)
}


#' Mean and covariance blocks of the expected information
#'
#' Shared by [hessian_loglik_unconstrained()],
#' [expected_fisher_unconstrained()] and the scenario descriptors:
#' `mean` is \eqn{\boldsymbol{\mu}^{\mathsf{T}}\boldsymbol{\Omega}
#' \boldsymbol{\mu}}, `cov` is \eqn{2p_jp_k\,\mathrm{Tr}(\boldsymbol{\Omega}
#' \boldsymbol{\Sigma}_j\boldsymbol{\Omega}\boldsymbol{\Sigma}_k)},
#' `traces` is \eqn{\mathrm{Tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}_j)}
#' and `precision` is \eqn{\boldsymbol{\Omega}(\boldsymbol{p})}.
#'
#' @noRd
.information_blocks <- function(p, mean_signature_matrix, Sigma) {
  omega <- .sigma_p_factorisation(p, Sigma)$inverse
  n_celltypes <- length(p)
  omega_sigma <- lapply(seq_len(n_celltypes), function(j) omega %*% Sigma[,, j])
  trace_products <- matrix(0, n_celltypes, n_celltypes)
  for (j in seq_len(n_celltypes)) {
    for (k in j:n_celltypes) {
      # Tr(A B) = sum(A * t(B))
      trace_products[j, k] <- sum(omega_sigma[[j]] * t(omega_sigma[[k]]))
      trace_products[k, j] <- trace_products[j, k]
    }
  }
  list(
    precision = omega,
    traces = vapply(omega_sigma, function(a) sum(diag(a)), numeric(1)),
    mean = unname(crossprod(
      mean_signature_matrix,
      omega %*% mean_signature_matrix
    )),
    cov = 2 * tcrossprod(p) * trace_products
  )
}


#' Constrained Hessian of \eqn{\ell\circ\boldsymbol{\psi}}
#'
#' @description
#' Second-order chain rule in ILR coordinates
#' \deqn{
#'   \mathbf{H}_{\boldsymbol{z}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   \mathbf{H}_{\boldsymbol{p}}
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   +\sum_{i=1}^{J}
#'   \frac{\partial\ell}{\partial p_i}\,
#'   \mathbf{H}_{\psi_i}.
#' }
#' The second summand cannot be dropped away from stationarity. At an
#' interior KKT point \eqn{\nabla_p\ell=\lambda\mathbf{1}} it vanishes
#' because \eqn{\sum_i\mathbf{H}_{\psi_i}=\mathbf{0}}.
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Symmetric matrix in \eqn{\mathcal{M}_{(J-1)\times(J-1)}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' hessian_loglik_constrained(isometric_log_ratio(p), y, mu, Sigma)
#' @export
hessian_loglik_constrained <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  .ilr_hessian(
    z,
    gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
    hessian_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
    V
  )
}


# DeCovarT core optimisation algorithms -----------------------------------

#' Open-simplex start for ILR solvers
#'
#' The convolution log-likelihood is not globally concave, so the start
#' can change which basin the solver enters. Draw several independent
#' Dirichlet starts (different RNG streams or
#' [multistart_decovart()]) rather than a single draw when probing
#' multimodality.
#'
#' @param n_celltypes Integer \eqn{J}.
#' @param initial_p One of:
#'   * `NULL` or `"barycentre"`: equi-balanced
#'     \eqn{(1/J,\ldots,1/J)};
#'   * a numeric vector of length \eqn{J}: used as-is after
#'     [repair_simplex()];
#'   * `"dirichlet"`: one Dirichlet\eqn{(\alpha,\ldots,\alpha)} draw
#'     (see `dirichlet_alpha`);
#'   * `"qp"` (aliases `"deconrnaseq"`, `"lsei"`): mean-only simplex QP
#'     from [deconvolute_ratios_deconrnaseq()] (`y` and
#'     `mean_signature_matrix` required).
#' @param nms Optional names (cell-type colnames).
#' @param y,mean_signature_matrix Bulk and signature, required for a QP
#'   start.
#' @param dirichlet_alpha Positive concentration, recycled to length
#'   \eqn{J}. The default `1` is uniform on the simplex.
#'   \eqn{\alpha>1} concentrates mass near the barycentre;
#'   \eqn{\alpha<1} puts extra mass near faces (boundary-biased
#'   restarts).
#'
#' @return A length-\eqn{J} open-simplex vector.
#'
#' @examples
#' starting_simplex(3L)
#' set.seed(1)
#' starting_simplex(3L, "dirichlet")
#' @export
#' @seealso [multistart_decovart()], [deconvolute_ratios_deconrnaseq()]
starting_simplex <- function(
  n_celltypes,
  initial_p = NULL,
  nms = NULL,
  y = NULL,
  mean_signature_matrix = NULL,
  dirichlet_alpha = 1
) {
  n_celltypes <- as.integer(n_celltypes)
  kind <- if (is.null(initial_p)) {
    "barycentre"
  } else if (is.character(initial_p) && length(initial_p) == 1L) {
    key <- tolower(initial_p)
    if (key %in% c("barycentre", "balanced", "equi")) {
      "barycentre"
    } else if (identical(key, "dirichlet")) {
      "dirichlet"
    } else if (key %in% c("qp", "deconrnaseq", "lsei")) {
      "qp"
    } else {
      stop(
        "`initial_p` character start must be 'barycentre', ",
        "'dirichlet', or 'qp'.",
        call. = FALSE
      )
    }
  } else if (is.numeric(initial_p) && length(initial_p) == n_celltypes) {
    "numeric"
  } else {
    stop(
      "`initial_p` must be NULL, a length-",
      n_celltypes,
      " numeric vector, or 'dirichlet' / 'qp'.",
      call. = FALSE
    )
  }
  if (identical(kind, "barycentre")) {
    p <- rep(1 / n_celltypes, n_celltypes)
  } else if (identical(kind, "numeric")) {
    p <- as.numeric(initial_p)
  } else if (identical(kind, "dirichlet")) {
    alpha <- rep_len(as.numeric(dirichlet_alpha), n_celltypes)
    if (
      !length(dirichlet_alpha) %in% c(1L, n_celltypes) ||
        any(!is.finite(alpha) | alpha <= 0)
    ) {
      stop(
        "`dirichlet_alpha` must be positive and finite, of length 1 or J.",
        call. = FALSE
      )
    }
    # One Dirichlet draw: normalised independent Gamma(alpha_j, 1) variates.
    p <- stats::rgamma(n_celltypes, shape = alpha, rate = 1)
    p <- p / sum(p)
  } else {
    if (is.null(y) || is.null(mean_signature_matrix)) {
      stop(
        "A QP start requires `y` and `mean_signature_matrix`.",
        call. = FALSE
      )
    }
    p <- deconvolute_ratios_deconrnaseq(y, mean_signature_matrix)
  }
  repair_simplex(p, open = TRUE, nms = nms)
}

#' Refit DeCovarT from several random starts to probe multimodality
#'
#' @description
#' The realised DeCovarT log-likelihood is not globally concave, so a
#' successful convergence code, a small relative-distance-to-minimum and a
#' negative-definite local Hessian together establish a **local** maximum
#' only. Refitting from independent Dirichlet\eqn{(1,\ldots,1)} starts and
#' comparing the attained log-likelihoods is the cheapest direct probe of
#' a second, equally good or better mode.
#'
#' @inheritParams loglik_multivariate
#' @param n_starts Number of Dirichlet restarts (in addition to the
#'   equi-balanced start).
#' @param dirichlet_alpha Concentration for those random starts (default
#'   `1`). See [starting_simplex()].
#' @param solver Solver accepting `initial_p` and `return_model`;
#'   defaults to [deconvolute_ratios_Newton_Raphson()].
#' @param loglik_tol Two starts count as reaching the same mode when their
#'   log-likelihoods differ by less than this amount.
#' @param ... Passed to `solver` (e.g. `epsilon`, `itmax`).
#'
#' @return A list with `coefficients` (best fit), `loglik`,
#'   `loglik_range`, `starts` (matrix of per-start estimates),
#'   `logliks`, and `multimodal`.
#'
#' @seealso [boundary_diagnostics()], [fit_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' set.seed(1)
#' multistart_decovart(y, mu, Sigma, n_starts = 3L)$multimodal
multistart_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  n_starts = 5L,
  solver = deconvolute_ratios_Newton_Raphson,
  loglik_tol = 1e-4,
  dirichlet_alpha = 1,
  ...
) {
  n_celltypes <- ncol(mean_signature_matrix)
  n_starts <- as.integer(n_starts)
  starts <- vector("list", n_starts + 1L)
  starts[[1L]] <- starting_simplex(
    n_celltypes,
    "barycentre",
    colnames(mean_signature_matrix)
  )
  for (s in seq_len(n_starts)) {
    starts[[s + 1L]] <- starting_simplex(
      n_celltypes,
      "dirichlet",
      colnames(mean_signature_matrix),
      dirichlet_alpha = dirichlet_alpha
    )
  }
  estimates <- matrix(
    NA_real_,
    nrow = n_celltypes,
    ncol = length(starts),
    dimnames = list(colnames(mean_signature_matrix), NULL)
  )
  logliks <- rep(NA_real_, length(starts))
  for (s in seq_along(starts)) {
    fit <- tryCatch(
      suppressWarnings(solver(
        y = y,
        mean_signature_matrix = mean_signature_matrix,
        Sigma = Sigma,
        return_model = TRUE,
        initial_p = starts[[s]],
        ...
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    estimates[, s] <- fit$coefficients
    logliks[[s]] <- fit$loglik
  }
  if (!any(is.finite(logliks))) {
    stop("Every random start failed to converge.", call. = FALSE)
  }
  best <- which.max(logliks)
  loglik_range <- diff(range(logliks, na.rm = TRUE))
  list(
    coefficients = estimates[, best],
    loglik = logliks[[best]],
    loglik_range = loglik_range,
    starts = estimates,
    logliks = logliks,
    multimodal = loglik_range > loglik_tol
  )
}

#' Shared return value of the DeCovarT solvers
#'
#' Names the estimate by cell type and, when `return_model = TRUE`, adds
#' the ILR coordinates, the attained log-likelihood and one convergence
#' record with the same fields for every solver: `code` (the solver's
#' own code, `istop` for Marquardt--Levenberg), `iterations`, `message`
#' and `rdm` (relative distance to the maximum, Marquardt--Levenberg
#' only).
#'
#' @noRd
.solver_result <- function(
  estimated_p,
  estimated_z,
  y,
  mean_signature_matrix,
  Sigma,
  return_model,
  code,
  iterations,
  message = NULL,
  rdm = NA_real_
) {
  names(estimated_p) <- colnames(mean_signature_matrix)
  if (!isTRUE(return_model)) {
    return(estimated_p)
  }
  message <- as.character(message)
  list(
    coefficients = estimated_p,
    z = estimated_z,
    loglik = loglik_multivariate(
      estimated_p,
      y,
      mean_signature_matrix,
      Sigma
    ),
    convergence = list(
      code = as.integer(code),
      iterations = as.integer(iterations),
      message = if (length(message) == 1L && nzchar(message)) {
        message
      } else {
        NA_character_
      },
      rdm = as.numeric(rdm)
    )
  )
}

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
#' Optimisation is performed in unconstrained ILR coordinates
#' \eqn{\boldsymbol{z}\in\mathbb{R}^{J-1}} via
#' \eqn{\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{z})}
#' (Marquardt–Levenberg default; see other methods below and
#' the
#' [population generative model](https://bastienchassagnol.github.io/DeCovarT/articles/01-1-population-generative-optimisation.html)).
#'
#' @details
#' **Plug-in signature.** Argument `mean_signature_matrix` is the mean
#' \eqn{\boldsymbol{\mu}}, used as a proxy for the unobserved cell-type
#' profiles \eqn{\boldsymbol{x}_{\cdot j}}. This is the frequentist plug-in;
#' recovering sample-specific latents is a Bayesian / MAP problem
#' (`.map_gaussian_convolution()`).
#'
#' @inheritParams loglik_multivariate
#' @param epsilon,itmax Absolute convergence tolerance and maximum number of
#'   iterations for the optimiser (same roles as `reltol` / `maxit` in
#'   [stats::optim()]).
#' @param return_model If `TRUE`, return a named list with coefficients,
#'   ILR coordinates, log-likelihood and optimiser diagnostics instead of
#'   the proportion vector.
#' @param initial_p Optional start: `NULL` / `"barycentre"` (default
#'   equi-balanced), a length-\eqn{J} numeric vector, `"dirichlet"`, or
#'   `"qp"` (mean-only simplex QP). See [starting_simplex()]. Starts
#'   on a simplex face are nudged into the interior so the ILR map is
#'   defined (ILR methods) and so L-BFGS-B does not start on a
#'   degenerate \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}.
#' @param dirichlet_alpha Dirichlet concentration when
#'   `initial_p = "dirichlet"` (default `1`: uniform on the simplex;
#'   \eqn{\alpha>1} centre-biased; \eqn{\alpha<1} face-biased).
#'   Independent draws are independent restarts; use
#'   [multistart_decovart()] to sequence several.
#'
#' @return Named numeric vector \eqn{\hat{\boldsymbol{p}}} on the simplex
#'   (ILR methods), or that list when `return_model = TRUE`.
#'   Benchmark metrics are computed by [deconvolute_ratios()].
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' y <- drop(mu %*% c(0.6, 0.4) + rnorm(2, sd = 0.1))
#' deconvolute_ratios_Marquardt_Levenberg(y, mu, Sigma, itmax = 50)
#' @export
#' @seealso [deconvolute_ratios()], [isometric_logistic()],
#'   [starting_simplex()]
deconvolute_ratios_Marquardt_Levenberg <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
) {
  initial_p <- starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix),
    y,
    mean_signature_matrix,
    dirichlet_alpha
  )
  initial_z <- isometric_log_ratio(initial_p)
  # marqLevAlg() only negates fn/gr internally when minimize = FALSE; hess()
  # is passed through unchanged regardless of `minimize`, because its
  # Cholesky-type step-and-RDM routines (dsinv/dchole) always assume the
  # "minimisation" convention (positive-definite at the optimum). Since we
  # maximise the log-likelihood, hessian_loglik_constrained() is negative-
  # definite at the optimum and must be negated by hand, or dsinv() fails
  # (ier = -1) on every iteration: the relative-distance-to-minimum (RDM)
  # criterion of Commenges et al. (2006) then never leaves its "unevaluated"
  # sentinel (rdm = epsd + 1), so the algorithm always exhausts `maxiter`
  # (istop = 2) rather than stopping once genuinely close to the maximum.
  # Send the console output of marqLevAlg to the null device.
  fit <- withr::with_output_sink(
    nullfile(),
    marqLevAlg::marqLevAlg(
      b = initial_z,
      fn = loglik_multivariate_constrained,
      gr = gradient_loglik_constrained,
      hess = function(z, y, mean_signature_matrix, Sigma) {
        -hessian_loglik_constrained(z, y, mean_signature_matrix, Sigma)
      },
      epsa = epsilon,
      epsb = epsilon,
      epsd = epsilon,
      minimize = FALSE,
      multipleTry = 1,
      y = y,
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      maxiter = itmax
    )
  )
  if (anyNA(fit$b)) {
    warning(
      "marqLevAlg::marqLevAlg() returned no usable estimate (istop = ",
      fit$istop,
      "); falling back to the equi-balanced initial guess.",
      call. = FALSE
    )
    estimated_z <- initial_z
  } else {
    if (fit$istop != 1) {
      warning(
        "marqLevAlg::marqLevAlg() stopped without jointly satisfying the ",
        "three Commenges et al. (2006) convergence criteria (istop = ",
        fit$istop,
        ", rdm = ",
        signif(fit$rdm, 3),
        "); returning the last iterate.",
        call. = FALSE
      )
    }
    estimated_z <- fit$b
  }
  .solver_result(
    isometric_logistic(estimated_z),
    estimated_z,
    y,
    mean_signature_matrix,
    Sigma,
    return_model,
    code = fit$istop,
    iterations = fit$ni,
    rdm = fit$rdm
  )
}


#' @describeIn deconvolute_ratios_Marquardt_Levenberg Simulated annealing on
#'   ILR coordinates \eqn{\boldsymbol{z}} ([stats::optim()] with
#'   `method = "SANN"`).
#' @export
deconvolute_ratios_simulated_annealing <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
) {
  initial_p <- starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix),
    y,
    mean_signature_matrix,
    dirichlet_alpha
  )
  initial_z <- isometric_log_ratio(initial_p)
  # gr is not used in the simulated annealing approach
  # In SANN, maxit is the total number of point evaluations, not iterations
  fit <- stats::optim(
    par = initial_z,
    fn = loglik_multivariate_constrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax),
    method = "SANN"
  )
  .solver_result(
    repair_simplex(isometric_logistic(fit$par)),
    fit$par,
    y,
    mean_signature_matrix,
    Sigma,
    return_model,
    code = fit$convergence,
    iterations = fit$counts[["function"]],
    message = fit$message
  )
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Box-constrained L-BFGS-B
#'   ([stats::optim()] `method = "L-BFGS-B"`) on the closure
#'   \eqn{\boldsymbol{p}=\boldsymbol{w}/\mathbf{1}^{\mathsf{T}}\boldsymbol{w}}
#'   of the box \eqn{\boldsymbol{w}\ge\mathbf{0}}. The box covers the closed
#'   simplex, so the solver maximises the same likelihood as Newton--Raphson
#'   and can return exact zeros, and a stationary point of the box is one of
#'   the simplex. `epsilon` is the relative reduction of the log-likelihood
#'   at which the solver stops (`factr = epsilon / .Machine$double.eps`).
#'   Wald standard errors use the ILR expected-Fisher map
#'   [vcov_ilr_delta()] and are undefined on a face.
#' @export
deconvolute_ratios_L_BFGS_B <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
) {
  initial_p <- starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix),
    y,
    mean_signature_matrix,
    dirichlet_alpha
  )
  fit <- .lbfgsb_closure(
    as.numeric(initial_p),
    loglik = function(p) {
      loglik_multivariate(p, y, mean_signature_matrix, Sigma)
    },
    score = function(p) {
      gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma)
    },
    epsilon = epsilon,
    itmax = itmax
  )
  .solver_result(
    fit$par,
    if (all(fit$par > 0)) isometric_log_ratio(fit$par) else NULL,
    y,
    mean_signature_matrix,
    Sigma,
    return_model,
    code = fit$convergence,
    iterations = fit$iterations,
    message = fit$message
  )
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Newton–Raphson on ILR
#'   coordinates \eqn{\boldsymbol{z}} with the analytic score and observed
#'   Hessian, safeguarded by the PORT trust region of [stats::nlminb()].
#' @export
deconvolute_ratios_Newton_Raphson <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
) {
  initial_p <- starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix),
    y,
    mean_signature_matrix,
    dirichlet_alpha
  )
  initial_z <- isometric_log_ratio(initial_p)

  fit <- .newton_ilr(
    initial_z,
    objective = function(z) {
      loglik_multivariate_constrained(z, y, mean_signature_matrix, Sigma)
    },
    gradient = function(z) {
      gradient_loglik_constrained(z, y, mean_signature_matrix, Sigma)
    },
    hessian = function(z) {
      hessian_loglik_constrained(z, y, mean_signature_matrix, Sigma)
    },
    epsilon = epsilon,
    itmax = itmax
  )
  .solver_result(
    isometric_logistic(fit$par),
    fit$par,
    y,
    mean_signature_matrix,
    Sigma,
    return_model,
    code = fit$convergence,
    iterations = fit$iterations,
    message = fit$message
  )
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg BFGS quasi-Newton ascent
#'   on ILR coordinates \eqn{\boldsymbol{z}} ([stats::optim()]
#'   `method = "BFGS"`).
#' @export
deconvolute_ratios_gradient_descent <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
) {
  initial_p <- starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix),
    y,
    mean_signature_matrix,
    dirichlet_alpha
  )
  initial_z <- isometric_log_ratio(initial_p)

  fit <- stats::optim(
    par = initial_z,
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
  )
  .solver_result(
    repair_simplex(isometric_logistic(fit$par)),
    fit$par,
    y,
    mean_signature_matrix,
    Sigma,
    return_model,
    code = fit$convergence,
    iterations = fit$counts[["function"]],
    message = fit$message
  )
}
