# Boundary-aware inference for DeCovarT: maximum likelihood on simplex
# faces, likelihood-ratio tests calibrated by chi-bar-square mixtures when a
# null sits on a face, edge or vertex of the simplex, and the diagnostics
# that tell an interior optimum from a boundary one.
# See vignette("theory-population-inference").

# Maximum likelihood on a face ---------------------------------------------

#' Maximum likelihood on the closed simplex, optionally on a face
#'
#' @description
#' Maximises the DeCovarT log-likelihood, summed over the replicate bulk
#' columns of `y`, over the closed simplex, while holding a (possibly empty)
#' subset of coordinates at prescribed values. It is the fitting engine of
#' [lrt_decovart()], [confint_profile_decovart()], [bootstrap_decovart()]
#' and [reference_bootstrap_decovart()]. Writing \eqn{A} for the
#' constrained index set and \eqn{s=\sum_{j\in A}c_j}, the free block is
#' reparametrised as
#' \deqn{
#'   \boldsymbol{p}_{A^{c}}
#'   =
#'   (1-s)\,\boldsymbol{\psi}(\boldsymbol{z}),
#'   \qquad
#'   \boldsymbol{z}\in\mathbb{R}^{|A^{c}|-1},
#' }
#' with \eqn{\boldsymbol{\psi}} the isometric logistic map
#' ([isometric_logistic()]) on a Helmert basis of the free face.
#' The constrained coordinates are *substituted*
#' rather than pushed through a logarithm, so a null such as
#' \eqn{p_j=0} is representable exactly: this is what makes boundary
#' likelihood-ratio tests computable (see [lrt_decovart()]).
#'
#' @details
#' Zero mass on \eqn{A} keeps
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^{2}
#' \boldsymbol{\Sigma}_j} positive definite as long as at least one
#' coordinate remains strictly positive. The free face is fitted by
#' Newton-Raphson with the PORT trust region of [stats::nlminb()], from
#' its barycentre, with the pooled analytic score and Hessian chained
#' through \eqn{(1-s)\mathbf{J}_{\boldsymbol{\psi}}}. The bounded first
#' step matters: with several columns the pooled score is large, and an
#' unbounded quasi-Newton step can land on a vertex, where the likelihood
#' is flat, and stop there.
#'
#' Without `fixed`, the interior ILR optimum is compared with the restricted
#' fit on the face where it stops (coordinates below `1e-3` set to zero),
#' and the larger likelihood is kept: the ILR chart cannot represent an
#' exact zero, so this is how a boundary MLE is reached.
#'
#' @inheritParams loglik_multivariate
#' @param y Numeric vector of length \eqn{G}, or \eqn{G\times N} matrix of
#'   replicate bulk profiles that share one composition; their
#'   log-likelihoods are summed (the pooled, population-level fit).
#' @param fixed Numeric vector of constrained ratios named by cell types,
#'   e.g. `c(celltype_1 = 0)`; the names must be columns of
#'   `mean_signature_matrix`. `NULL` (default) fits the unrestricted MLE.
#' @param epsilon,itmax Relative convergence tolerance and iteration
#'   budget passed to [stats::nlminb()].
#'
#' @return A list with `coefficients` (the fitted \eqn{J} vector),
#'   `loglik` (summed over the columns of `y`), `fixed`, and `convergence`.
#'
#' @seealso [lrt_decovart()], [confint_profile_decovart()],
#'   [loglik_multivariate()]
#' @family decovart_inference
#' @keywords internal
#' @export
#' @examples
#' mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
#' colnames(mu) <- paste0("ct", 1:3)
#' Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
#' y <- drop(mu %*% c(0.5, 0.5, 0))
#' restricted_mle_decovart(y, mu, Sigma, fixed = c(ct3 = 0))$coefficients
restricted_mle_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  fixed = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  y <- as.matrix(y)
  if (nrow(y) != nrow(mean_signature_matrix)) {
    stop(
      "`y` must have ",
      nrow(mean_signature_matrix),
      " genes (rows), got ",
      nrow(y),
      ".",
      call. = FALSE
    )
  }
  .assert_no_missing(y, "bulk expression")
  nms <- colnames(mean_signature_matrix)
  interior <- .fit_on_face(
    nms,
    fixed,
    loglik = function(p) {
      .pooled_population_loglik(p, y, mean_signature_matrix, Sigma)
    },
    score = function(p) {
      .pooled_population_score(p, y, mean_signature_matrix, Sigma)
    },
    hessian = function(p) {
      .pooled_population_hessian(p, y, mean_signature_matrix, Sigma)
    },
    epsilon = epsilon,
    itmax = itmax
  )
  n_celltypes <- length(nms)
  if (length(fixed) == 0L) {
    face <- which(interior$coefficients < 1e-3)
    if (length(face) > 0L && length(face) < n_celltypes) {
      on_face <- restricted_mle_decovart(
        y,
        mean_signature_matrix,
        Sigma,
        fixed = stats::setNames(rep(0, length(face)), nms[face]),
        epsilon = epsilon,
        itmax = itmax
      )
      if (is.finite(on_face$loglik) && on_face$loglik >= interior$loglik) {
        return(on_face)
      }
    }
  }
  interior
}

# Boundary-aware null distribution -----------------------------------------

#' Chi-bar-square tail probability for boundary likelihood-ratio tests
#'
#' @description
#' Upper-tail probability of the mixture
#' \deqn{
#'   \bar{\chi}^{2}
#'   =
#'   \sum_{i=0}^{q}
#'   w_i\,\chi^{2}_{s+i},
#'   \qquad
#'   w_i=\binom{q}{i}2^{-q},
#' }
#' with \eqn{q} constraints active on the boundary of the parameter space
#' and \eqn{s} interior constraints. This is Case 9 of
#' \insertCite{selfAsymptoticPropertiesMaximum1987}{DeCovarT}, which
#' generalises the \eqn{\tfrac{1}{2}\chi^{2}_{0}+\tfrac{1}{2}\chi^{2}_{1}}
#' result of \insertCite{hermanDistributionLikelihoodRatio1954}{DeCovarT}.
#' For \eqn{q=0} the mixture collapses to the ordinary
#' \eqn{\chi^{2}_{s}} of \insertCite{wilksLargeSampleDistributionLikelihood1938}{DeCovarT}.
#'
#' @details
#' The binomial weights hold when the block of the Fisher information
#' corresponding to the active constraints is diagonal after
#' orthogonalisation against the nuisance block. Under correlated
#' constraints the mixing weights depend on the geometry of the
#' approximating tangent cone and the binomial choice is only an
#' approximation; Self and Liang also exhibit a nuisance-parameter
#' configuration whose limit is **not** a chi-square mixture. Supply
#' `weights` explicitly, or calibrate by simulation with
#' [bootstrap_decovart()], whenever \eqn{q>1} and the constrained
#' cell types share signal.
#'
#' @param statistic Numeric vector of observed likelihood-ratio
#'   statistics \eqn{D=2(\hat{\ell}_1-\hat{\ell}_0)}.
#' @param n_boundary Integer \eqn{q}: number of constraints active on the
#'   boundary (e.g. null ratios set to zero).
#' @param df_interior Integer \eqn{s}: number of interior (two-sided)
#'   constraints tested simultaneously.
#' @param weights Optional numeric vector of length `n_boundary + 1`
#'   giving the mixing weights \eqn{w_0,\ldots,w_q}. Defaults to the
#'   binomial weights.
#'
#' @return Numeric vector of upper-tail probabilities.
#'
#' @section Atom at zero:
#' \eqn{\chi^{2}_{0}} is a point mass at the origin, so its strict upper
#' tail \eqn{\Pr(\chi^{2}_{0}>D)} is zero for every \eqn{D\ge 0}. The
#' component therefore only removes weight: with one active constraint and
#' no interior constraint the p-value is
#' \eqn{\tfrac{1}{2}\Pr(\chi^{2}_{1}>D)}, capped at \eqn{0.5} when
#' \eqn{D=0}. Note [stats::pchisq()] switches to the *closed* tail at
#' exactly `q = 0` for `df = 0`; the strict convention is used here so the
#' p-value stays continuous in `statistic`.
#'
#' @references
#' \insertRef{selfAsymptoticPropertiesMaximum1987}{DeCovarT}
#'
#' \insertRef{hermanDistributionLikelihoodRatio1954}{DeCovarT}
#'
#' \insertRef{wilksLargeSampleDistributionLikelihood1938}{DeCovarT}
#'
#' @seealso [lrt_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' # One active boundary constraint: 50:50 mixture of chi2_0 and chi2_1.
#' chi_bar_square_pvalue(2.71, n_boundary = 1L)
#' # Interior (Wilks) calibration.
#' chi_bar_square_pvalue(3.84, n_boundary = 0L, df_interior = 1L)
chi_bar_square_pvalue <- function(
  statistic,
  n_boundary = 1L,
  df_interior = 0L,
  weights = NULL
) {
  n_boundary <- as.integer(n_boundary)
  df_interior <- as.integer(df_interior)
  if (n_boundary < 0L || df_interior < 0L) {
    stop("`n_boundary` and `df_interior` must be non-negative.", call. = FALSE)
  }
  if (n_boundary == 0L && df_interior == 0L) {
    stop("At least one constraint must be tested.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- stats::dbinom(0:n_boundary, size = n_boundary, prob = 0.5)
  }
  if (length(weights) != n_boundary + 1L) {
    stop("`weights` must have length `n_boundary + 1`.", call. = FALSE)
  }
  weights <- weights / sum(weights)
  degrees <- df_interior + 0:n_boundary
  tails <- vapply(
    degrees,
    function(df) {
      if (df == 0L) {
        # Atom at the origin: strict upper tail is identically zero.
        return(rep(0, length(statistic)))
      }
      stats::pchisq(statistic, df = df, lower.tail = FALSE)
    },
    numeric(length(statistic))
  )
  if (is.null(dim(tails))) {
    tails <- matrix(tails, nrow = 1L)
  }
  drop(tails %*% weights)
}

# Likelihood-ratio test ----------------------------------------------------

#' Number of constraints a null imposes, and how many lie on a face
#'
#' Fixing \eqn{k} coordinates whose sum is below one imposes \eqn{k}
#' constraints, one per fixed zero lying on a face. When the fixed values
#' sum to one, every unfixed coordinate is forced to zero as well, so the
#' null fixes the whole composition: \eqn{J-1} constraints, and every zero
#' coordinate, fixed or implied, counts as a boundary constraint. A null
#' \eqn{p_j=1} therefore has \eqn{J-1} boundary constraints.
#'
#' @noRd
.null_constraint_counts <- function(null_value, n_celltypes, tol) {
  null_value <- as.numeric(null_value)
  if (sum(null_value) >= 1 - tol) {
    n_constraints <- n_celltypes - 1L
    implied_zero <- n_celltypes - length(null_value)
    n_boundary <- sum(null_value <= tol) + implied_zero
  } else {
    n_constraints <- length(null_value)
    n_boundary <- sum(null_value <= tol)
  }
  n_boundary <- as.integer(min(n_boundary, n_constraints))
  list(
    n_constraints = as.integer(n_constraints),
    n_boundary = n_boundary,
    df_interior = as.integer(n_constraints) - n_boundary
  )
}

#' Likelihood-ratio test for cellular ratios
#'
#' @description
#' Tests \eqn{H_0:p_j=c_j,\,j\in A} against the unrestricted simplex with
#' the profile likelihood-ratio statistic
#' \deqn{
#'   D
#'   =
#'   2\bigl\{
#'     \ell(\hat{\boldsymbol{p}})
#'     -\ell(\hat{\boldsymbol{p}}^{(0)})
#'   \bigr\},
#' }
#' where \eqn{\hat{\boldsymbol{p}}^{(0)}} is the restricted MLE of
#' [restricted_mle_decovart()]. Unlike the Wald interval of
#' [confint.decovart_fit()], \eqn{D} is invariant under the additive
#' log-ratio reparametrisation, because a smooth bijection leaves
#' likelihood ratios unchanged.
#'
#' @details
#' **Interior null.** If every tested ratio lies strictly inside
#' \eqn{(0,1)}, \eqn{D} is asymptotically \eqn{\chi^{2}_{|A|}}
#' \insertCite{wilksLargeSampleDistributionLikelihood1938}{DeCovarT}.
#'
#' **Boundary null.** A null such as \eqn{p_j=0} lets the parameter move
#' in one direction only, so the local parameter space is a half-line
#' rather than a vector space and Wilks' theorem fails. The limit is the
#' chi-bar-square mixture of [chi_bar_square_pvalue()]: for a single
#' active constraint,
#' \eqn{D\rightsquigarrow\tfrac{1}{2}\chi^{2}_{0}+\tfrac{1}{2}\chi^{2}_{1}}
#' \insertCite{hermanDistributionLikelihoodRatio1954;textual}{DeCovarT},
#' generalised to several active constraints by
#' \insertCite{selfAsymptoticPropertiesMaximum1987;textual}{DeCovarT}.
#' Boundary calibration is selected automatically when any tested value is
#' within `boundary_tol` of 0, or when the tested values sum to one, which
#' forces every other ratio to zero (a null \eqn{p_j=1} fixes \eqn{J-1}
#' coordinates on faces); override with `n_boundary`. The binomial weights
#' of [chi_bar_square_pvalue()] are exact for orthogonal constraints; for
#' several correlated faces they approximate the weights of the tangent
#' cone, and [bootstrap_decovart()] with `null_value` gives a calibration
#' that does not rely on them.
#'
#' **Replication.** Both calibrations are *population* statements: they
#' need replicate bulk samples that genuinely share one composition
#' \eqn{\boldsymbol{p}}. A single \eqn{G}-vector supplies no replication
#' (genes are dependent through \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})},
#' and adding bulk samples with distinct compositions adds parameters as
#' fast as observations). Pass several replicate columns in `y`, or
#' calibrate by [bootstrap_decovart()].
#'
#' @inheritParams restricted_mle_decovart
#' @param null_value Numeric vector of hypothesised ratios named by cell
#'   types, e.g. `c(celltype_2 = 0)`.
#' @param n_boundary Optional integer override for the number of active
#'   boundary constraints; `NULL` infers it from `null_value`, counting
#'   the zeros implied when the fixed values sum to one.
#' @param boundary_tol Distance from 0 or 1 below which a null value
#'   counts as being on the boundary.
#'
#' @return A one-row data frame with columns `statistic`,
#'   `n_boundary`, `df_interior`, `p_value`, `loglik_full`,
#'   `loglik_null`, and `calibration`.
#'
#' @references
#' \insertRef{wilksLargeSampleDistributionLikelihood1938}{DeCovarT}
#'
#' \insertRef{selfAsymptoticPropertiesMaximum1987}{DeCovarT}
#'
#' @seealso [chi_bar_square_pvalue()], [confint_profile_decovart()],
#'   [bootstrap_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
#' colnames(mu) <- paste0("ct", 1:3)
#' Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
#' y <- drop(mu %*% c(0.5, 0.5, 0))
#' lrt_decovart(y, mu, Sigma, null_value = c(ct3 = 0))
lrt_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  null_value,
  n_boundary = NULL,
  boundary_tol = 1e-6,
  epsilon = 10^-8,
  itmax = 500L
) {
  full <- restricted_mle_decovart(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  null <- restricted_mle_decovart(
    y,
    mean_signature_matrix,
    Sigma,
    fixed = null_value,
    epsilon = epsilon,
    itmax = itmax
  )
  # The alternative is the CLOSED simplex, so its supremum is at least the
  # restricted maximum. An interior ILR optimiser only approaches a face,
  # hence a small negative gap is the signature of a genuine boundary MLE
  # rather than of a failed optimisation.
  gap <- null$loglik - full$loglik
  if (gap > 1e-3) {
    warning(
      "The restricted fit beat the unrestricted fit by ",
      signif(gap, 3),
      " log-likelihood units, which points to a failed unrestricted ",
      "optimisation rather than a boundary optimum. Retry with a larger ",
      "`itmax` or screen starts with multistart_decovart().",
      call. = FALSE
    )
  }
  loglik_full <- max(full$loglik, null$loglik)
  statistic <- 2 * (loglik_full - null$loglik)
  counts <- .null_constraint_counts(
    null_value,
    ncol(mean_signature_matrix),
    boundary_tol
  )
  if (!is.null(n_boundary)) {
    counts$df_interior <- counts$n_constraints - as.integer(n_boundary)
    counts$n_boundary <- as.integer(n_boundary)
  }
  n_boundary <- counts$n_boundary
  df_interior <- counts$df_interior
  data.frame(
    statistic = statistic,
    n_boundary = n_boundary,
    df_interior = df_interior,
    p_value = chi_bar_square_pvalue(
      statistic,
      n_boundary = n_boundary,
      df_interior = df_interior
    ),
    loglik_full = loglik_full,
    loglik_null = null$loglik,
    calibration = if (n_boundary > 0L) {
      "chi-bar-square (Self-Liang)"
    } else {
      "chi-square (Wilks)"
    },
    row.names = toString(names(null_value)),
    stringsAsFactors = FALSE
  )
}

# Boundary diagnostics -----------------------------------------------------

#' Boundary and stationarity diagnostics for one DeCovarT fit
#'
#' @description
#' Separates three claims that a bare optimiser return code conflates:
#' the solver stopped, the iterate is genuinely stationary with the right
#' local curvature, and the estimate sits close to a simplex face. None of
#' them implies global uniqueness, which the DeCovarT log-likelihood does
#' not have for a single observed sample (see
#' `vignette("theory-population-inference")`).
#'
#' @details
#' Reported fields are the ILR score norm
#' \eqn{\lVert\nabla_{\boldsymbol{z}}\ell\rVert}, the extreme eigenvalues
#' \eqn{\lambda_{\min}(\mathbf{H}_{\boldsymbol{z}})} and
#' \eqn{\lambda_{\max}(\mathbf{H}_{\boldsymbol{z}})} of the ILR Hessian
#' (a local maximum requires the Hessian to be negative definite, so
#' both must be strictly negative; \eqn{\lambda_{\min}<0<\lambda_{\max}}
#' is a saddle), `boundary_distance` \eqn{=\min_j\hat{p}_j}, and the
#' flags `near_boundary` and `local_maximum`.
#'
#' `boundary_tol` is a **statistical** warning threshold for Wald / ILR
#' linearisation, deliberately much larger than the machine-precision guard
#' that decides whether a logarithm is representable. A fit with
#' \eqn{\min_j\hat{p}_j\ll 1} is not evidence of optimiser failure: the
#' solver may be correctly approaching a genuine boundary optimum. It *is*
#' evidence that interior Wald intervals should be replaced by the profile
#' or boundary-calibrated tests of [lrt_decovart()] and
#' [confint_profile_decovart()].
#'
#' @inheritParams loglik_multivariate
#' @param p Estimated proportions of length \eqn{J}.
#' @param boundary_tol Threshold on \eqn{\min_j\hat{p}_j} below which the
#'   estimate is flagged as near-boundary.
#' @param score_tol Threshold on the ILR score norm below which the
#'   iterate counts as stationary.
#'
#' @return A one-row data frame of diagnostics.
#'
#' @seealso [fit_decovart()], [lrt_decovart()]
#' @family decovart_inference
#' @keywords internal
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' boundary_diagnostics(c(0.6, 0.4), y, mu, Sigma)
boundary_diagnostics <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  boundary_tol = 1e-8,
  score_tol = 1e-4
) {
  boundary_distance <- min(p)
  near_boundary <- boundary_distance < boundary_tol
  score_norm <- NA_real_
  min_eigenvalue <- NA_real_
  max_eigenvalue <- NA_real_
  if (!near_boundary && length(p) > 1L) {
    z <- isometric_log_ratio(p)
    score_norm <- sqrt(sum(
      gradient_loglik_constrained(z, y, mean_signature_matrix, Sigma)^2
    ))
    hessian_z <- hessian_loglik_constrained(
      z,
      y,
      mean_signature_matrix,
      Sigma
    )
    ev <- eigen(hessian_z, symmetric = TRUE, only.values = TRUE)$values
    min_eigenvalue <- min(ev)
    max_eigenvalue <- max(ev)
  }
  negative_definite <- isTRUE(max_eigenvalue < 0) &&
    isTRUE(min_eigenvalue < 0)
  data.frame(
    boundary_distance = boundary_distance,
    near_boundary = near_boundary,
    score_norm = score_norm,
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = max_eigenvalue,
    local_maximum = isTRUE(score_norm < score_tol) &&
      isTRUE(negative_definite),
    stringsAsFactors = FALSE
  )
}
