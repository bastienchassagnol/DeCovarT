# ILR chart of the simplex, shared by the population convolution
# (R/02_02_population_decovart_estimate_ratios.R) and the single-cell
# mixture (R/02_02_single_cell_decovart_estimate_ratios.R).

#' Repair a numeric vector onto the unit simplex
#'
#' @description
#' Clips numerical under/overflow and renormalises so that
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1} and
#' \eqn{\boldsymbol{p}\ge\mathbf{0}}. This is a **repair /
#' renormalisation** step for estimated proportions, not a Euclidean
#' projection onto the simplex and not a statistical-identifiability
#' constraint.
#'
#' When `open = TRUE`, entries on the closed-simplex boundary are nudged
#' into the relative interior so that [isometric_log_ratio()] is defined
#' (ILR / Marquardt / Newton starts). That open-simplex path is the
#' former role of [starting_simplex()].
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param tolerance Non-negative tolerance for treating entries as zero
#'   (default `100 * .Machine$double.eps`).
#' @param open Logical; if `TRUE`, push the repaired vector off the
#'   boundary into the open simplex.
#' @param nms Optional names attached to the returned vector (e.g.
#'   cell-type colnames).
#'
#' @return Numeric vector on the simplex \eqn{\Delta^{J-1}} (open when
#'   `open = TRUE`).
#'
#' @seealso [compositions::clo()] for compositional closure.
#'
#' @srrstats {G3.0} Default `tolerance` is `100 * .Machine$double.eps`;
#'   callers must not test floating-point equality with `==`.
#'
#' @examples
#' repair_simplex(c(0.2, 0.3, 0.5 + 1e-12))
#' repair_simplex(c(1, 0, 0), open = TRUE)
#' @export
repair_simplex <- function(
  p,
  tolerance = 100 * .Machine$double.eps,
  open = FALSE,
  nms = NULL
) {
  if (!is.numeric(p) || length(p) == 0L || anyNA(p)) {
    stop("`p` must be a non-empty numeric vector without missing values.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0) {
    stop("`tolerance` must be a single non-negative number.")
  }

  p[abs(p) <= tolerance] <- 0

  if (any(p < 0)) {
    stop("`p` contains negative values beyond numerical tolerance.")
  }

  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    stop("`p` must have a positive, finite sum.")
  }

  p <- p / total
  p[abs(p - 1) <= tolerance] <- 1

  if (isTRUE(open)) {
    floor <- tolerance
    if (any(p < floor) || any(p > 1 - floor)) {
      p <- pmax(p, floor)
      p <- p / sum(p)
    }
  }

  if (!is.null(nms)) {
    names(p) <- nms
  }
  p
}

#' Helmert contrast matrix for isometric log-ratio coordinates
#'
#' @description
#' Returns the \eqn{J\times(J-1)} Helmert sub-matrix \eqn{\mathbf{V}}
#' with
#' \deqn{
#'   \mathbf{V}^{\mathsf{T}}\mathbf{V}=\mathbf{I}_{J-1},\qquad
#'   \mathbf{V}^{\mathsf{T}}\mathbf{1}=\mathbf{0}.
#' }
#' Column \eqn{k} contrasts the first \eqn{k} parts against part
#' \eqn{k+1}. This is the standard ILR basis used with
#' [isometric_log_ratio()] / [isometric_logistic()]; any other valid ILR
#' basis is \eqn{\mathbf{V}Q} for an orthogonal \eqn{Q}, which rotates
#' coordinates but leaves simplex geometry, eigenvalues of quadratic
#' forms, and \eqn{\mathrm{Var}(\hat{\boldsymbol{p}})} unchanged.
#'
#' @param n_parts Integer \eqn{J\ge 2} (number of simplex parts).
#'
#' @return Numeric matrix \eqn{J\times(J-1)}.
#'
#' @examples
#' V <- helmert_basis(3L)
#' crossprod(V)
#' drop(crossprod(V, rep(1, 3)))
#' @export
helmert_basis <- function(n_parts) {
  j <- as.integer(n_parts)
  if (length(j) != 1L || is.na(j) || j < 2L) {
    stop("`n_parts` must be an integer of at least 2.", call. = FALSE)
  }
  v <- matrix(0, nrow = j, ncol = j - 1L)
  for (k in seq_len(j - 1L)) {
    scale_k <- sqrt(k * (k + 1))
    v[seq_len(k), k] <- 1 / scale_k
    v[k + 1L, k] <- -k / scale_k
  }
  v
}

#' Isometric logistic transform (ILR coordinates to the simplex)
#'
#' @description
#' Inverse isometric log-ratio map
#' \eqn{\boldsymbol{\psi}:\boldsymbol{z}\mapsto\boldsymbol{p}} with a
#' Helmert basis \eqn{\mathbf{V}},
#' \deqn{
#'   \boldsymbol{p}
#'   =
#'   \operatorname{softmax}(\mathbf{V}\boldsymbol{z})
#'   =
#'   \mathcal{C}\bigl(\exp(\mathbf{V}\boldsymbol{z})\bigr).
#' }
#' No cell type is pinned as a reference. Evaluated with a log-sum-exp
#' shift so large
#' \eqn{\|\boldsymbol{z}\|} cannot overflow.
#'
#' @param z Numeric vector \eqn{\boldsymbol{z}\in\mathbb{R}^{J-1}}.
#' @param V Optional Helmert (or other ILR) basis with \eqn{J} rows and
#'   \eqn{J-1} columns. Defaults to [helmert_basis()] for
#'   \eqn{J=\mathrm{length}(z)+1}.
#'
#' @return Numeric vector \eqn{\boldsymbol{p}} on the unit simplex.
#'
#' @examples
#' z <- c(0.2, -0.5)
#' p <- isometric_logistic(z)
#' sum(p)
#' isometric_log_ratio(p)
#' isometric_logistic(c(800, 1200))
#' @seealso [isometric_log_ratio()], [helmert_basis()]
#' @export
isometric_logistic <- function(z, V = NULL) {
  z <- as.numeric(z)
  if (is.null(V)) {
    V <- helmert_basis(length(z) + 1L)
  }
  eta <- drop(V %*% z)
  weights <- exp(eta - max(eta))
  weights / sum(weights)
}

#' Isometric log-ratio transform \eqn{\boldsymbol{p}\mapsto\boldsymbol{z}}
#'
#' \deqn{
#'   \boldsymbol{z}
#'   =
#'   \mathbf{V}^{\mathsf{T}}\log\boldsymbol{p}
#'   =
#'   \mathbf{V}^{\mathsf{T}}\operatorname{clr}(\boldsymbol{p}),
#' }
#' which is well-defined because \eqn{\mathbf{V}^{\mathsf{T}}\mathbf{1}=0}.
#' Requires a strictly positive composition (open simplex).
#'
#' @param p Numeric vector on the open simplex.
#' @param V Optional ILR basis; see [isometric_logistic()].
#'
#' @return Numeric vector \eqn{\boldsymbol{z}\in\mathbb{R}^{J-1}}.
#'
#' @rdname isometric_logistic
#' @export
isometric_log_ratio <- function(p, V = NULL) {
  p <- as.numeric(p)
  if (is.null(V)) {
    V <- helmert_basis(length(p))
  }
  drop(crossprod(V, log(p)))
}

#' Jacobian of the isometric logistic map
#'
#' \deqn{
#'   \mathbf{J}_{\boldsymbol{\psi}}(\boldsymbol{z})
#'   =
#'   \frac{\partial\boldsymbol{p}}{\partial\boldsymbol{z}^{\mathsf{T}}}
#'   =
#'   \mathbf{S}(\boldsymbol{p})\mathbf{V},
#' }
#' with \eqn{\mathbf{S}(\boldsymbol{p})=\operatorname{diag}(\boldsymbol{p})
#' -\boldsymbol{p}\boldsymbol{p}^{\mathsf{T}}}.
#'
#' @inheritParams isometric_logistic
#'
#' @return Numeric matrix \eqn{J\times(J-1)}.
#'
#' @examples
#' jacobian_isometric_logistic(c(0.2, -0.5))
#' @export
jacobian_isometric_logistic <- function(z, V = NULL) {
  z <- as.numeric(z)
  if (is.null(V)) {
    V <- helmert_basis(length(z) + 1L)
  }
  p <- isometric_logistic(z, V)
  s <- diag(p) - tcrossprod(p)
  s %*% V
}

#' Hessian tensor of the isometric logistic map
#'
#' Slice \eqn{i} is
#' \deqn{
#'   \mathbf{H}_{\psi_i}
#'   =
#'   p_i\bigl(\boldsymbol{q}_i\boldsymbol{q}_i^{\mathsf{T}}
#'   -\mathbf{C}_V\bigr),
#' }
#' where \eqn{\bar{\boldsymbol{v}}=\mathbf{V}^{\mathsf{T}}\boldsymbol{p}},
#' \eqn{\boldsymbol{q}_i=\mathbf{V}_{i\cdot}-\bar{\boldsymbol{v}}} and
#' \eqn{\mathbf{C}_V=\mathbf{V}^{\mathsf{T}}\mathbf{S}(\boldsymbol{p})
#' \mathbf{V}}. Array shape \eqn{(J-1)\times(J-1)\times J}.
#'
#' @inheritParams isometric_logistic
#'
#' @return Numeric array used in the constrained Hessian chain rule.
#'
#' @examples
#' hessian_isometric_logistic(c(0.2, -0.5))
#' @export
hessian_isometric_logistic <- function(z, V = NULL) {
  z <- as.numeric(z)
  if (is.null(V)) {
    V <- helmert_basis(length(z) + 1L)
  }
  p <- isometric_logistic(z, V)
  j <- length(p)
  d <- ncol(V)
  s <- diag(p) - tcrossprod(p)
  c_v <- crossprod(V, s %*% V)
  vbar <- drop(crossprod(p, V))
  out <- array(0, dim = c(d, d, j))
  for (i in seq_len(j)) {
    q <- V[i, ] - vbar
    out[,, i] <- p[i] * (tcrossprod(q) - c_v)
  }
  out
}

#' Score in ILR coordinates from the score in the composition
#'
#' First-order chain rule \eqn{\mathbf{J}_{\psi}^{\top}U}.
#'
#' @noRd
.ilr_gradient <- function(z, score_p, V = NULL) {
  drop(crossprod(jacobian_isometric_logistic(z, V), score_p))
}

#' Hessian in ILR coordinates from the score and Hessian in the composition
#'
#' Second-order chain rule
#' \eqn{\mathbf{J}_{\psi}^{\top}\mathbf{H}\mathbf{J}_{\psi}
#' +\sum_{i}U_i\mathbf{H}_{\psi_i}}. The second term vanishes only at an
#' interior stationary point.
#'
#' @noRd
.ilr_hessian <- function(z, score_p, hessian_p, V = NULL) {
  jac <- jacobian_isometric_logistic(z, V)
  curvature <- as.matrix(tensor::tensor(
    A = score_p,
    B = hessian_isometric_logistic(z, V),
    alongA = 1L,
    alongB = 3L
  ))
  out <- crossprod(jac, hessian_p %*% jac) + curvature
  (out + t(out)) / 2
}
