# Asymptotic inference for DeCovarT: profile likelihood, likelihood-ratio
# tests with boundary (chi-bar-square) calibration, and the parametric
# bootstrap. See vignette("DeCovarT-MLE-properties").

# Restricted maximisation ------------------------------------------------

#' Restricted maximum likelihood with fixed cellular ratios
#'
#' @description
#' Maximises the DeCovarT log-likelihood over the simplex while holding a
#' subset of coordinates at prescribed values. Writing \eqn{A} for the
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
#' coordinate remains strictly positive. Optimisation uses
#' [stats::optim()] (`"BFGS"`) with the analytic score obtained by
#' chaining [gradient_loglik_unconstrained()] through
#' \eqn{(1-s)\mathbf{J}_{\boldsymbol{\psi}}} on the free-face ILR chart.
#'
#' @inheritParams loglik_multivariate
#' @param fixed Named or integer-indexed numeric vector of constrained
#'   ratios, e.g. `c(celltype_1 = 0)`. Names are matched against
#'   `colnames(mean_signature_matrix)`.
#' @param epsilon,itmax Relative convergence tolerance and iteration
#'   budget passed to [stats::optim()].
#'
#' @return A list with `coefficients` (the restricted \eqn{J} vector),
#'   `loglik`, `fixed`, and `convergence`.
#'
#' @seealso [lrt_decovart()], [profile_loglik_decovart()],
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
  fixed,
  epsilon = 10^-8,
  itmax = 500L
) {
  nms <- colnames(mean_signature_matrix)
  n_celltypes <- ncol(mean_signature_matrix)
  idx <- .resolve_celltype_index(names(fixed), nms, n_celltypes)
  fixed_value <- as.numeric(fixed)
  if (anyNA(fixed_value) || any(fixed_value < 0) || sum(fixed_value) > 1) {
    stop(
      "`fixed` must hold non-negative ratios summing to at most 1.",
      call. = FALSE
    )
  }
  free <- setdiff(seq_len(n_celltypes), idx)
  remaining <- 1 - sum(fixed_value)
  assemble <- function(z) {
    p <- numeric(n_celltypes)
    p[idx] <- fixed_value
    n_free <- length(free)
    if (n_free == 1L) {
      p[free] <- remaining
    } else if (n_free > 1L) {
      if (remaining <= 0) {
        p[free] <- 0
      } else {
        p[free] <- remaining * isometric_logistic(z)
      }
    }
    names(p) <- nms
    p
  }
  if (length(free) <= 1L || remaining <= 0) {
    p_hat <- assemble(numeric(0))
    return(list(
      coefficients = p_hat,
      loglik = loglik_multivariate(p_hat, y, mean_signature_matrix, Sigma),
      fixed = stats::setNames(fixed_value, nms[idx]),
      convergence = list(code = 0L, iterations = 0L, message = "closed form")
    ))
  }
  objective <- function(z) {
    loglik_multivariate(assemble(z), y, mean_signature_matrix, Sigma)
  }
  score <- function(z) {
    p <- assemble(z)
    grad_p <- gradient_loglik_unconstrained(
      p,
      y,
      mean_signature_matrix,
      Sigma
    )
    drop(grad_p[free] %*% (remaining * jacobian_isometric_logistic(z)))
  }
  fit <- stats::optim(
    par = rep(0, length(free) - 1L),
    fn = objective,
    gr = score,
    method = "BFGS",
    control = list(fnscale = -1, reltol = epsilon, maxit = itmax)
  )
  list(
    coefficients = assemble(fit$par),
    loglik = fit$value,
    fixed = stats::setNames(fixed_value, nms[idx]),
    convergence = list(
      code = fit$convergence,
      iterations = unname(fit$counts[["function"]]),
      message = fit$message
    )
  )
}

#' Profile log-likelihood of one cellular ratio
#'
#' @description
#' Returns \eqn{\ell^{\mathrm{prof}}_j(c)=\max\{\ell(\boldsymbol{p}):
#' \boldsymbol{p}\in\Delta^{J-1},\,p_j=c\}}, i.e. the restricted maximum of
#' [restricted_mle_decovart()] for a single coordinate. Profiling
#' concentrates out the remaining \eqn{J-2} free proportions instead of
#' relying on a quadratic (Wald) approximation, so it is reliable even
#' when \eqn{\hat{p}_j} sits close to a simplex face.
#'
#' @inheritParams restricted_mle_decovart
#' @param celltype Cell-type name or integer index \eqn{j}.
#' @param value Numeric vector of candidate ratios \eqn{c\in[0,1]}.
#'
#' @return Numeric vector of profile log-likelihoods, named by `value`.
#'
#' @seealso [lrt_decovart()], [confint_profile_decovart()]
#' @family decovart_inference
#' @keywords internal
#' @export
#' @examples
#' mu <- matrix(c(20, 22, 18, 22, 20, 24), nrow = 2)
#' colnames(mu) <- paste0("ct", 1:3)
#' Sigma <- array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
#' y <- drop(mu %*% c(0.5, 0.3, 0.2))
#' profile_loglik_decovart(y, mu, Sigma, "ct3", c(0, 0.2, 0.4))
profile_loglik_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  celltype,
  value,
  epsilon = 10^-8,
  itmax = 500L
) {
  nms <- colnames(mean_signature_matrix)
  idx <- .resolve_celltype_index(
    celltype,
    nms,
    ncol(mean_signature_matrix)
  )
  if (length(idx) != 1L) {
    stop("`celltype` must identify exactly one cell type.", call. = FALSE)
  }
  out <- vapply(
    value,
    function(v) {
      fixed <- stats::setNames(v, nms[idx])
      restricted_mle_decovart(
        y,
        mean_signature_matrix,
        Sigma,
        fixed = fixed,
        epsilon = epsilon,
        itmax = itmax
      )$loglik
    },
    numeric(1)
  )
  stats::setNames(out, format(value, trim = TRUE))
}

# Boundary-aware null distribution ---------------------------------------

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

# Likelihood-ratio test --------------------------------------------------

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
#' within `boundary_tol` of 0 or 1; override with `n_boundary`.
#'
#' **Replication.** Both calibrations are *population* statements: they
#' need replicate bulk samples that genuinely share one composition
#' \eqn{\boldsymbol{p}}. A single \eqn{G}-vector supplies no replication
#' (genes are dependent through \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})},
#' and adding bulk samples with distinct compositions adds parameters as
#' fast as observations). Use `bulk_expression` with several replicate
#' columns, or calibrate by [bootstrap_decovart()].
#'
#' @inheritParams restricted_mle_decovart
#' @param bulk_expression Numeric matrix \eqn{\boldsymbol{Y}\in
#'   \mathcal{M}_{G\times N}} of replicate bulk profiles sharing one
#'   composition, used instead of `y` for the pooled (population) test.
#' @param null_value Named numeric vector of hypothesised ratios, e.g.
#'   `c(celltype_2 = 0)`.
#' @param n_boundary Optional integer override for the number of active
#'   boundary constraints; `NULL` infers it from `null_value`.
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
  y = NULL,
  mean_signature_matrix,
  Sigma,
  null_value,
  bulk_expression = NULL,
  n_boundary = NULL,
  boundary_tol = 1e-6,
  epsilon = 10^-8,
  itmax = 500L
) {
  y_mat <- .as_replicate_matrix(y, bulk_expression, nrow(mean_signature_matrix))
  nms <- colnames(mean_signature_matrix)
  idx <- .resolve_celltype_index(
    names(null_value),
    nms,
    ncol(mean_signature_matrix)
  )
  full <- .pooled_mle(
    y_mat,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  null <- .pooled_mle(
    y_mat,
    mean_signature_matrix,
    Sigma,
    fixed = stats::setNames(as.numeric(null_value), nms[idx]),
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
  on_boundary <- as.numeric(null_value) <= boundary_tol |
    as.numeric(null_value) >= 1 - boundary_tol
  if (is.null(n_boundary)) {
    n_boundary <- sum(on_boundary)
  }
  n_boundary <- as.integer(n_boundary)
  df_interior <- length(null_value) - n_boundary
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
    row.names = toString(nms[idx]),
    stringsAsFactors = FALSE
  )
}

#' Profile-likelihood confidence intervals for cellular ratios
#'
#' @description
#' Inverts the likelihood-ratio test of [lrt_decovart()] coordinate by
#' coordinate: the interval collects every \eqn{c} with
#' \eqn{2\{\ell(\hat{\boldsymbol{p}})-\ell^{\mathrm{prof}}_j(c)\}
#' \le\chi^{2}_{1,1-\alpha}}. Because likelihood ratios are invariant
#' under reparametrisation, the interval respects
#' \eqn{0\le p_j\le 1} without the delta-method linearisation that makes
#' Wald intervals unreliable near a simplex face.
#'
#' @inheritParams lrt_decovart
#' @param level Confidence level \eqn{1-\alpha}.
#' @param celltypes Cell types to profile; defaults to all.
#'
#' @return A matrix with one row per cell type and columns `estimate`,
#'   `lower`, `upper`.
#'
#' @seealso [confint.decovart_fit()] for the Wald analogue,
#'   [lrt_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' confint_profile_decovart(y, mu, Sigma)
confint_profile_decovart <- function(
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  level = 0.95,
  celltypes = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  y_mat <- .as_replicate_matrix(y, bulk_expression, nrow(mean_signature_matrix))
  nms <- colnames(mean_signature_matrix)
  n_celltypes <- ncol(mean_signature_matrix)
  idx <- if (is.null(celltypes)) {
    seq_len(n_celltypes)
  } else {
    .resolve_celltype_index(celltypes, nms, n_celltypes)
  }
  full <- .pooled_mle(
    y_mat,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  cutoff <- stats::qchisq(level, df = 1L)
  out <- matrix(
    NA_real_,
    nrow = length(idx),
    ncol = 3L,
    dimnames = list(nms[idx], c("estimate", "lower", "upper"))
  )
  for (k in seq_along(idx)) {
    j <- idx[[k]]
    p_hat <- full$coefficients[[j]]
    deviance <- function(v) {
      restricted <- .pooled_mle(
        y_mat,
        mean_signature_matrix,
        Sigma,
        fixed = stats::setNames(v, nms[j]),
        epsilon = epsilon,
        itmax = itmax
      )
      2 * (full$loglik - restricted$loglik) - cutoff
    }
    out[k, ] <- c(
      p_hat,
      .profile_root(deviance, lower = 0, upper = p_hat),
      .profile_root(deviance, lower = p_hat, upper = 1)
    )
  }
  out
}

# Parametric bootstrap ---------------------------------------------------

#' Parametric bootstrap for DeCovarT proportions and boundary tests
#'
#' @description
#' The conditional law of the bulk profile is fully specified, so the
#' asymptotic calibrations of [lrt_decovart()] can be replaced by
#' simulation. Given plug-in moments and a composition
#' \eqn{\boldsymbol{p}}, the function draws
#' \eqn{\boldsymbol{Y}^{(b)}\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))}, refits DeCovarT, and returns the
#' empirical distribution of \eqn{\hat{\boldsymbol{p}}^{(b)}} together
#' with percentile intervals.
#'
#' @details
#' Supplying `null_value` switches to a restricted bootstrap: replicates
#' are generated under the restricted MLE and the likelihood-ratio
#' statistic is recomputed on each, giving a Monte Carlo p-value that does
#' not rely on the \eqn{50{:}50} asymptotic weights. For a composite null
#' the remaining nuisance proportions are profiled at their restricted
#' estimate, so the procedure is a plug-in (approximate) rather than an
#' exact bootstrap.
#'
#' Reference uncertainty is **not** propagated here: \eqn{\boldsymbol{\mu}}
#' and \eqn{\boldsymbol{\Sigma}_j} are treated as known, exactly as in the
#' frequentist likelihood. To resample the purified observations that
#' produced those plug-in moments, use
#' [reference_bootstrap_decovart()].
#'
#' @inheritParams lrt_decovart
#' @param p Composition to simulate from. Defaults to the MLE computed
#'   from `y` / `bulk_expression`.
#' @param n_boot Number of bootstrap replicates.
#' @param level Confidence level for percentile intervals.
#' @param n_replicates Number of bulk columns simulated per replicate;
#'   defaults to the number of observed columns.
#'
#' @return A list with `estimates` (\eqn{J\times} `n_boot` matrix),
#'   `interval` (percentile bounds), `p_simulated`, and, when
#'   `null_value` is supplied, `statistic` (observed \eqn{D}),
#'   `null_statistics`, and `p_value`.
#'
#' @seealso [lrt_decovart()], [chi_bar_square_pvalue()],
#'   [reference_bootstrap_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' set.seed(1)
#' bootstrap_decovart(y, mu, Sigma, n_boot = 20)$interval
bootstrap_decovart <- function(
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  p = NULL,
  null_value = NULL,
  n_boot = 500L,
  level = 0.95,
  n_replicates = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  y_mat <- .as_replicate_matrix(y, bulk_expression, nrow(mean_signature_matrix))
  nms <- colnames(mean_signature_matrix)
  if (is.null(n_replicates)) {
    n_replicates <- ncol(y_mat)
  }
  observed <- .pooled_mle(
    y_mat,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  statistic <- NULL
  if (is.null(null_value)) {
    p_simulated <- if (is.null(p)) observed$coefficients else p
  } else {
    idx <- .resolve_celltype_index(
      names(null_value),
      nms,
      ncol(mean_signature_matrix)
    )
    restricted <- .pooled_mle(
      y_mat,
      mean_signature_matrix,
      Sigma,
      fixed = stats::setNames(as.numeric(null_value), nms[idx]),
      epsilon = epsilon,
      itmax = itmax
    )
    p_simulated <- if (is.null(p)) restricted$coefficients else p
    statistic <- max(2 * (observed$loglik - restricted$loglik), 0)
  }
  estimates <- matrix(
    NA_real_,
    nrow = length(p_simulated),
    ncol = n_boot,
    dimnames = list(nms, NULL)
  )
  null_statistics <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    y_boot <- .simulate_bulk_replicates(
      p_simulated,
      mean_signature_matrix,
      Sigma,
      n_replicates
    )
    boot_full <- .pooled_mle(
      y_boot,
      mean_signature_matrix,
      Sigma,
      epsilon = epsilon,
      itmax = itmax
    )
    estimates[, b] <- boot_full$coefficients
    if (!is.null(statistic)) {
      boot_null <- .pooled_mle(
        y_boot,
        mean_signature_matrix,
        Sigma,
        fixed = stats::setNames(as.numeric(null_value), nms[idx]),
        epsilon = epsilon,
        itmax = itmax
      )
      null_statistics[[b]] <- 2 *
        (max(boot_full$loglik, boot_null$loglik) - boot_null$loglik)
    }
  }
  alpha <- (1 - level) / 2
  interval <- t(apply(
    estimates,
    1L,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE
  ))
  out <- list(
    estimates = estimates,
    interval = interval,
    p_simulated = p_simulated
  )
  if (!is.null(statistic)) {
    out$statistic <- statistic
    out$null_statistics <- null_statistics
    # +1 smoothing keeps the Monte Carlo p-value strictly positive.
    out$p_value <- (1 + sum(null_statistics >= statistic, na.rm = TRUE)) /
      (1 + sum(is.finite(null_statistics)))
  }
  out
}

# Label-equivariance sanity check ----------------------------------------

#' Permutation-equivariance check for a labelled reference
#'
#' @description
#' Software check, not a bootstrap and not a test of the
#' bulk-to-signature match. In a *reference-based* model the cell-type
#' names are anchored by the signature columns, so shuffling those names
#' (or the gene axis) is not a valid uncertainty procedure. What *is*
#' required is that the estimator be equivariant: if
#' \eqn{\boldsymbol{\mu}^{\star}=\boldsymbol{\mu}Q} for a permutation
#' matrix \eqn{Q}, the fitted composition should satisfy
#' \eqn{\hat{\boldsymbol{p}}^{\star}\approx Q^{\top}\hat{\boldsymbol{p}}}
#' and the reconstructed bulk should be unchanged.
#'
#' @details
#' Columns of \eqn{\boldsymbol{\mu}} and matching slices of
#' \eqn{\boldsymbol{\Sigma}_j} are reordered by `perm` while the *names*
#' stay in the original order, matching the convention of
#' `vignette("DeCovarT-MLE-properties")`. The labelled MLE on the
#' relabelled reference is then compared with
#' \eqn{\hat{\boldsymbol{p}}_{\mathrm{perm}}}.
#'
#' Use [reference_bootstrap_decovart()] for percentile intervals that
#' resample the experimental units of the reference (donors or cells)
#' or redraw compositions from a Dirichlet law.
#'
#' @inheritParams bootstrap_decovart
#' @param perm Integer permutation of `1:J`. Defaults to reversing the
#'   column order so the check is deterministic.
#'
#' @return A list with `perm`, `p_hat`, `p_star` (fit on the relabelled
#'   reference), `p_expected` (\eqn{\hat{\boldsymbol{p}}} reordered by
#'   `perm`), `max_abs_diff`, and `loglik_diff`.
#'
#' @seealso [reference_bootstrap_decovart()], [loglik_multivariate()]
#' @family decovart_inference
#' @keywords internal
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.7, 0.3))
#' equivariance_check_decovart(y, mu, Sigma)$max_abs_diff
equivariance_check_decovart <- function(
  y = NULL,
  mean_signature_matrix,
  Sigma,
  bulk_expression = NULL,
  perm = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  y_mat <- .as_replicate_matrix(y, bulk_expression, nrow(mean_signature_matrix))
  n_celltypes <- ncol(mean_signature_matrix)
  if (is.null(perm)) {
    perm <- rev(seq_len(n_celltypes))
  }
  perm <- as.integer(perm)
  if (
    length(perm) != n_celltypes ||
      !identical(sort(perm), seq_len(n_celltypes))
  ) {
    stop("`perm` must be a permutation of 1:J.", call. = FALSE)
  }
  observed <- .pooled_mle(
    y_mat,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  relabelled <- .permute_celltype_reference(
    mean_signature_matrix,
    Sigma,
    perm
  )
  fitted <- .pooled_mle(
    y_mat,
    relabelled$mean_signature_matrix,
    relabelled$Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  p_hat <- observed$coefficients
  p_expected <- p_hat[perm]
  names(p_expected) <- names(p_hat)
  list(
    perm = perm,
    p_hat = p_hat,
    p_star = fitted$coefficients,
    p_expected = p_expected,
    max_abs_diff = max(abs(fitted$coefficients - p_expected)),
    loglik_diff = fitted$loglik - observed$loglik
  )
}

# Reference-sample bootstrap ---------------------------------------------

#' Reference-based bootstrap for signature and composition uncertainty
#'
#' @description
#' The frequentist DeCovarT likelihood treats
#' \eqn{\boldsymbol{\mu}} and each \eqn{\boldsymbol{\Sigma}_j} as known.
#' In a labelled, reference-based model those moments come from purified
#' profiles whose *cell-type names stay attached*. This bootstrap
#' resamples the experimental units that generated the reference (donors,
#' by default), or redraws compositions from a Dirichlet law, then
#' re-estimates the plug-in moments and refits
#' \insertCite{efronBootstrapMethodsAnother1979}{DeCovarT}.
#'
#' @details
#' Permuting gene order, or permuting cell-type labels of an already
#' averaged signature, is **not** a bootstrap: the first destroys the
#' gene-wise pairing of bulk and reference, and the second is algebraic
#' equivariance rather than a source of uncertainty (see
#' [equivariance_check_decovart()] and
#' `vignette("DeCovarT-MLE-properties")`). The three `method` options
#' instead resample units that actually vary.
#'
#' \describe{
#'   \item{`donors`}{Cluster bootstrap: resample donor identifiers
#'     *within each cell type*, take every purified column of the
#'     sampled donors (a donor drawn twice is included twice), rebuild
#'     \eqn{(\boldsymbol{\mu},\boldsymbol{\Sigma}_j)}, and refit.
#'     This is the default: it targets biological / reference-population
#'     uncertainty rather than technical cell-level noise alone.}
#'   \item{`cells`}{Resample purified columns independently within each
#'     cell type (with replacement) and rebuild the plug-in moments.
#'     Use this when donor labels are unavailable.}
#'   \item{`dirichlet`}{Hold the original plug-in moments fixed, draw
#'     \eqn{\boldsymbol{p}^{(b)}\sim\mathrm{Dirichlet}(\boldsymbol{\alpha})},
#'     simulate bulk profiles from the Gaussian convolution, and refit.
#'     This is a composition-sweep diagnostic, not a confidence interval
#'     for one observed bulk.}
#' }
#'
#' For `donors` and `cells` the observed bulk is kept fixed unless
#' `regenerate_bulk = TRUE`, in which case each replicate also draws
#' \eqn{\boldsymbol{Y}^{(b)}} from the bootstrapped moments at `p`
#' (the observed MLE by default). That is the pipeline
#' reference units \eqn{\to S^{(b)}\to Y^{(b)}\to\hat{\boldsymbol{p}}^{(b)}}.
#'
#' When a cell type has fewer purified columns than genes, the sample
#' covariance is singular. A ridge `ridge * mean(diag(S))` is then added
#' on the diagonal so that \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} remains
#' factorisable. Set `ridge = 0` to disable that guard (and accept
#' possible Cholesky failures).
#'
#' @inheritParams bootstrap_decovart
#' @param reference_profiles Named list of \eqn{G\times n_j} matrices,
#'   one purified sample matrix per cell type. Names must match the
#'   intended column order of the signature. Alternatively, a single
#'   \eqn{G\times n_{\mathrm{ref}}} matrix together with
#'   `cell_type_labels`.
#' @param cell_type_labels Optional character or factor vector of length
#'   `ncol(reference_profiles)` when `reference_profiles` is a matrix.
#' @param donor_ids Donor (or other clustering) labels aligned with the
#'   purified columns. A named list of vectors matching
#'   `reference_profiles`, or a vector of length
#'   `ncol(reference_profiles)` when the reference is a single matrix.
#'   Required for `method = "donors"`.
#' @param method One of `"donors"` (default), `"cells"`, or
#'   `"dirichlet"`.
#' @param dirichlet_alpha Positive Dirichlet concentration: a scalar
#'   recycled to \eqn{J}, or a length-\eqn{J} vector. Ignored unless
#'   `method = "dirichlet"`. The default `1` is uniform on the simplex.
#' @param regenerate_bulk If `TRUE` (`donors` / `cells` only), simulate
#'   a new bulk from the bootstrapped moments instead of keeping the
#'   observed \eqn{\boldsymbol{Y}} fixed.
#' @param ridge Non-negative ridge multiple for the sample covariance of
#'   each cell type. The default `1e-6` is relative to the mean diagonal.
#'
#' @return A list with `method`, `estimates` (\eqn{J\times} `n_boot`),
#'   `interval` (percentile bounds), `mean_signature_matrix` / `Sigma`
#'   from the original reference, and `p_hat` from the original plug-in
#'   fit when a bulk is supplied. For `method = "dirichlet"`,
#'   `p_simulated` stores the drawn compositions.
#'
#' @seealso [bootstrap_decovart()], [equivariance_check_decovart()]
#' @references
#' \insertRef{efronBootstrapMethodsAnother1979}{DeCovarT}
#' @family decovart_inference
#' @export
#' @examples
#' set.seed(1)
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' refs <- lapply(seq_len(2), function(j) {
#'   t(MASS::mvrnorm(n = 8, mu = mu[, j], Sigma = Sigma[, , j]))
#' })
#' names(refs) <- colnames(mu)
#' donor_ids <- lapply(refs, function(x) {
#'   rep(paste0("d", 1:2), length.out = ncol(x))
#' })
#' y <- drop(mu %*% c(0.6, 0.4))
#' reference_bootstrap_decovart(
#'   y,
#'   refs,
#'   donor_ids = donor_ids,
#'   n_boot = 15
#' )$interval
reference_bootstrap_decovart <- function(
  y = NULL,
  reference_profiles,
  bulk_expression = NULL,
  cell_type_labels = NULL,
  donor_ids = NULL,
  method = c("donors", "cells", "dirichlet"),
  n_boot = 199L,
  level = 0.95,
  dirichlet_alpha = 1,
  regenerate_bulk = FALSE,
  p = NULL,
  n_replicates = NULL,
  ridge = 1e-6,
  epsilon = 10^-8,
  itmax = 500L
) {
  method <- match.arg(method)
  bundle <- .as_reference_bundle(
    reference_profiles,
    cell_type_labels,
    donor_ids
  )
  refs <- bundle$refs
  donor_list <- bundle$donor_ids
  moments0 <- .moments_from_profiles(refs, ridge = ridge)
  nms <- colnames(moments0$mean_signature_matrix)
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  if (identical(method, "donors")) {
    if (is.null(donor_list)) {
      stop(
        "Supply `donor_ids` for method = \"donors\", ",
        "or use method = \"cells\".",
        call. = FALSE
      )
    }
    .check_donors_per_type(donor_list)
  }

  has_bulk <- !is.null(y) || !is.null(bulk_expression)
  if (!identical(method, "dirichlet") && !has_bulk) {
    stop(
      "Supply `y` or `bulk_expression` for method = \"",
      method,
      "\".",
      call. = FALSE
    )
  }

  y_mat <- NULL
  observed_coef <- NULL
  if (has_bulk) {
    y_mat <- .as_replicate_matrix(
      y,
      bulk_expression,
      nrow(moments0$mean_signature_matrix)
    )
    observed <- .pooled_mle(
      y_mat,
      moments0$mean_signature_matrix,
      moments0$Sigma,
      epsilon = epsilon,
      itmax = itmax
    )
    observed_coef <- observed$coefficients
  }
  if (is.null(n_replicates)) {
    n_replicates <- if (is.null(y_mat)) 1L else ncol(y_mat)
  }
  p_sim_base <- if (is.null(p)) observed_coef else p
  if (
    identical(method, "dirichlet") ||
      isTRUE(regenerate_bulk)
  ) {
    if (is.null(p_sim_base) && !identical(method, "dirichlet")) {
      stop(
        "`p` is required to regenerate the bulk when no observed ",
        "fit is available.",
        call. = FALSE
      )
    }
  }

  alpha_dir <- NULL
  if (identical(method, "dirichlet")) {
    alpha_dir <- .parse_dirichlet_alpha(dirichlet_alpha, nms)
  }

  estimates <- matrix(
    NA_real_,
    nrow = length(nms),
    ncol = n_boot,
    dimnames = list(nms, NULL)
  )
  p_simulated <- NULL
  if (identical(method, "dirichlet")) {
    p_simulated <- matrix(
      NA_real_,
      nrow = length(nms),
      ncol = n_boot,
      dimnames = list(nms, NULL)
    )
  }

  for (b in seq_len(n_boot)) {
    if (identical(method, "dirichlet")) {
      p_b <- .rdirichlet_one(alpha_dir)
      p_simulated[, b] <- p_b
      y_fit <- .simulate_bulk_replicates(
        p_b,
        moments0$mean_signature_matrix,
        moments0$Sigma,
        n_replicates
      )
      mu_fit <- moments0$mean_signature_matrix
      sigma_fit <- moments0$Sigma
    } else {
      refs_b <- if (identical(method, "donors")) {
        .resample_donors_within_type(refs, donor_list)
      } else {
        .resample_cells_within_type(refs)
      }
      moments_b <- .moments_from_profiles(refs_b, ridge = ridge)
      mu_fit <- moments_b$mean_signature_matrix
      sigma_fit <- moments_b$Sigma
      if (isTRUE(regenerate_bulk)) {
        y_fit <- .simulate_bulk_replicates(
          p_sim_base,
          mu_fit,
          sigma_fit,
          n_replicates
        )
      } else {
        y_fit <- y_mat
      }
    }
    fit_b <- .pooled_mle(
      y_fit,
      mu_fit,
      sigma_fit,
      epsilon = epsilon,
      itmax = itmax
    )
    estimates[, b] <- fit_b$coefficients
  }
  alpha <- (1 - level) / 2
  interval <- t(apply(
    estimates,
    1L,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE
  ))
  out <- list(
    method = method,
    estimates = estimates,
    interval = interval,
    p_hat = observed_coef,
    mean_signature_matrix = moments0$mean_signature_matrix,
    Sigma = moments0$Sigma
  )
  if (!is.null(p_simulated)) {
    out$p_simulated <- p_simulated
  }
  out
}

#' Coerce purified profiles and optional donor labels
#'
#' @noRd
.as_reference_bundle <- function(
  reference_profiles,
  cell_type_labels,
  donor_ids
) {
  if (is.list(reference_profiles) && !is.data.frame(reference_profiles)) {
    refs <- .as_reference_list(reference_profiles, cell_type_labels = NULL)
    return(list(
      refs = refs,
      donor_ids = .align_donor_ids(refs, donor_ids)
    ))
  }
  mat <- as.matrix(reference_profiles)
  if (is.null(cell_type_labels)) {
    stop(
      "Supply `cell_type_labels` when `reference_profiles` is a matrix.",
      call. = FALSE
    )
  }
  if (length(cell_type_labels) != ncol(mat)) {
    stop(
      "`cell_type_labels` must have length ncol(reference_profiles).",
      call. = FALSE
    )
  }
  labs <- as.character(cell_type_labels)
  split_idx <- split(seq_len(ncol(mat)), factor(labs, levels = unique(labs)))
  refs <- lapply(split_idx, function(idx) mat[, idx, drop = FALSE])
  .check_reference_list(refs)
  donors <- NULL
  if (!is.null(donor_ids)) {
    if (length(donor_ids) != ncol(mat)) {
      stop(
        "`donor_ids` must have length ncol(reference_profiles).",
        call. = FALSE
      )
    }
    d <- as.character(donor_ids)
    donors <- lapply(split_idx, function(idx) d[idx])
  }
  list(refs = refs, donor_ids = donors)
}

#' Coerce purified profiles to a named list of gene-by-sample matrices
#'
#' @noRd
.as_reference_list <- function(reference_profiles, cell_type_labels) {
  if (is.list(reference_profiles) && !is.data.frame(reference_profiles)) {
    refs <- lapply(reference_profiles, as.matrix)
    if (is.null(names(refs)) || any(!nzchar(names(refs)))) {
      stop("`reference_profiles` must be a named list.", call. = FALSE)
    }
    n_genes <- vapply(refs, nrow, integer(1))
    if (length(unique(n_genes)) != 1L) {
      stop(
        "Every reference matrix must have the same number of genes.",
        call. = FALSE
      )
    }
    .check_reference_list(refs)
    return(refs)
  }
  mat <- as.matrix(reference_profiles)
  if (is.null(cell_type_labels)) {
    stop(
      "Supply `cell_type_labels` when `reference_profiles` is a matrix.",
      call. = FALSE
    )
  }
  if (length(cell_type_labels) != ncol(mat)) {
    stop(
      "`cell_type_labels` must have length ncol(reference_profiles).",
      call. = FALSE
    )
  }
  labs <- as.character(cell_type_labels)
  split_idx <- split(seq_len(ncol(mat)), factor(labs, levels = unique(labs)))
  refs <- lapply(split_idx, function(idx) mat[, idx, drop = FALSE])
  .check_reference_list(refs)
  refs
}

#' @noRd
.check_reference_list <- function(refs) {
  for (nm in names(refs)) {
    if (ncol(refs[[nm]]) < 2L) {
      stop(
        "Cell type ",
        nm,
        " needs at least two purified columns.",
        call. = FALSE
      )
    }
  }
  invisible(refs)
}

#' @noRd
.align_donor_ids <- function(refs, donor_ids) {
  if (is.null(donor_ids)) {
    return(NULL)
  }
  if (!(is.list(donor_ids) && !is.data.frame(donor_ids))) {
    stop(
      "`donor_ids` must be a named list matching `reference_profiles`.",
      call. = FALSE
    )
  }
  if (!setequal(names(donor_ids), names(refs))) {
    stop(
      "`donor_ids` names must match `reference_profiles`.",
      call. = FALSE
    )
  }
  out <- lapply(names(refs), function(nm) {
    d <- donor_ids[[nm]]
    if (length(d) != ncol(refs[[nm]])) {
      stop(
        "Donor labels for cell type ",
        nm,
        " must have length ncol of that matrix.",
        call. = FALSE
      )
    }
    as.character(d)
  })
  names(out) <- names(refs)
  out
}

#' @noRd
.check_donors_per_type <- function(donor_ids) {
  n_donors <- vapply(
    donor_ids,
    function(d) length(unique(d)),
    integer(1)
  )
  bad <- names(donor_ids)[n_donors < 2L]
  if (length(bad) > 0L) {
    stop(
      "method = \"donors\" needs at least two donors in every cell type (",
      toString(bad),
      ").",
      call. = FALSE
    )
  }
  invisible(donor_ids)
}

#' @noRd
.resample_cells_within_type <- function(refs) {
  lapply(refs, function(x) {
    x[, sample.int(ncol(x), replace = TRUE), drop = FALSE]
  })
}

#' @noRd
.resample_donors_within_type <- function(refs, donor_ids) {
  out <- Map(
    function(x, d) {
      donors <- unique(d)
      boot <- sample(donors, size = length(donors), replace = TRUE)
      cols <- unlist(
        lapply(boot, function(id) which(d == id)),
        use.names = FALSE
      )
      x[, cols, drop = FALSE]
    },
    refs,
    donor_ids[names(refs)]
  )
  names(out) <- names(refs)
  out
}

#' @noRd
.parse_dirichlet_alpha <- function(alpha, nms) {
  n_celltypes <- length(nms)
  if (length(alpha) == 1L) {
    alpha <- rep(alpha, n_celltypes)
  }
  if (length(alpha) != n_celltypes) {
    stop(
      "`dirichlet_alpha` must be length 1 or J.",
      call. = FALSE
    )
  }
  if (any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop("`dirichlet_alpha` must be positive and finite.", call. = FALSE)
  }
  stats::setNames(as.numeric(alpha), nms)
}

#' @noRd
.rdirichlet_one <- function(alpha) {
  draws <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  draws / sum(draws)
}

#' Plug-in mean signature and covariances from purified profiles
#'
#' @noRd
.moments_from_profiles <- function(reference_profiles, ridge = 1e-6) {
  nms <- names(reference_profiles)
  n_genes <- nrow(reference_profiles[[1L]])
  n_celltypes <- length(reference_profiles)
  gene_names <- rownames(reference_profiles[[1L]])
  mu <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = n_celltypes,
    dimnames = list(gene_names, nms)
  )
  Sigma <- array(
    NA_real_,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, nms)
  )
  for (j in seq_len(n_celltypes)) {
    x <- reference_profiles[[j]]
    mu[, j] <- rowMeans(x)
    s <- stats::cov(t(x))
    if (ridge > 0) {
      s <- s + diag(ridge * mean(diag(s)), n_genes)
    }
    Sigma[,, j] <- s
  }
  list(mean_signature_matrix = mu, Sigma = Sigma)
}

#' Relabel columns of the signature and matching covariance slices
#'
#' Used to demonstrate that the log-likelihood is equivariant under a
#' joint permutation of cell-type labels. Names are kept in the original
#' order so the *labelled* composition is the permuted vector.
#'
#' @noRd
.permute_celltype_reference <- function(
  mean_signature_matrix,
  Sigma,
  perm
) {
  nms <- colnames(mean_signature_matrix)
  mu <- mean_signature_matrix[, perm, drop = FALSE]
  colnames(mu) <- nms
  sig <- Sigma[,, perm, drop = FALSE]
  if (!is.null(dimnames(Sigma)[[3]])) {
    dimnames(sig)[[3]] <- nms
  }
  list(mean_signature_matrix = mu, Sigma = sig)
}

# Convergence and boundary diagnostics -----------------------------------

#' Boundary and stationarity diagnostics for one DeCovarT fit
#'
#' @description
#' Separates three claims that a bare optimiser return code conflates:
#' the solver stopped, the iterate is genuinely stationary with the right
#' local curvature, and the estimate sits close to a simplex face. None of
#' them implies global uniqueness, which the DeCovarT log-likelihood does
#' not have for a single observed sample (see
#' `vignette("DeCovarT-MLE-properties")`).
#'
#' @details
#' Reported fields are the ILR score norm
#' \eqn{\lVert\nabla_{\boldsymbol{z}}\ell\rVert}, the largest
#' eigenvalue \eqn{\lambda_{\max}(\mathbf{H}_{\boldsymbol{z}})} (negative
#' at a local maximum), `boundary_distance` \eqn{=\min_j\hat{p}_j}, and the
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
  curvature <- NA_real_
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
    curvature <- max(
      eigen(hessian_z, symmetric = TRUE, only.values = TRUE)$values
    )
  }
  data.frame(
    boundary_distance = boundary_distance,
    near_boundary = near_boundary,
    score_norm = score_norm,
    max_eigenvalue = curvature,
    local_maximum = isTRUE(score_norm < score_tol) && isTRUE(curvature < 0),
    stringsAsFactors = FALSE
  )
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
#' @param n_starts Number of random starts (in addition to the
#'   equi-balanced start).
#' @param solver Solver accepting `initial_p` and `return_model`;
#'   defaults to [deconvolute_ratios_Marquardt_Levenberg()].
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
  solver = deconvolute_ratios_Marquardt_Levenberg,
  loglik_tol = 1e-4,
  ...
) {
  n_celltypes <- ncol(mean_signature_matrix)
  n_starts <- as.integer(n_starts)
  starts <- vector("list", n_starts + 1L)
  starts[[1L]] <- rep(1 / n_celltypes, n_celltypes)
  for (s in seq_len(n_starts)) {
    gammas <- stats::rgamma(n_celltypes, shape = 1, rate = 1)
    starts[[s + 1L]] <- gammas / sum(gammas)
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

# Internal helpers -------------------------------------------------------

#' Resolve cell-type names or indices against a signature matrix
#'
#' @param x Character names or numeric indices.
#' @param nms Column names of the signature matrix.
#' @param n_celltypes Number of cell types.
#'
#' @noRd
.resolve_celltype_index <- function(x, nms, n_celltypes) {
  if (is.null(x) || length(x) == 0L) {
    stop(
      "Cell types must be named or indexed; got an empty selector.",
      call. = FALSE
    )
  }
  if (is.character(x)) {
    idx <- match(x, nms)
    if (anyNA(idx)) {
      stop(
        "Unknown cell type(s): ",
        toString(x[is.na(idx)]),
        ".",
        call. = FALSE
      )
    }
    return(idx)
  }
  idx <- as.integer(x)
  if (anyNA(idx) || any(idx < 1L) || any(idx > n_celltypes)) {
    stop("Cell-type indices must lie in 1:", n_celltypes, ".", call. = FALSE)
  }
  idx
}

#' Coerce one bulk vector or a replicate matrix to a \eqn{G\times N} matrix
#'
#' @param y Optional numeric vector.
#' @param bulk_expression Optional numeric matrix.
#' @param n_genes Expected number of genes.
#'
#' @noRd
.as_replicate_matrix <- function(y, bulk_expression, n_genes) {
  if (!is.null(bulk_expression)) {
    mat <- as.matrix(bulk_expression)
  } else if (!is.null(y)) {
    mat <- matrix(as.numeric(y), ncol = 1L)
  } else {
    stop("Supply either `y` or `bulk_expression`.", call. = FALSE)
  }
  if (nrow(mat) != n_genes) {
    stop(
      "Bulk profiles must have ",
      n_genes,
      " genes, got ",
      nrow(mat),
      ".",
      call. = FALSE
    )
  }
  .assert_no_missing(mat, "bulk expression")
  mat
}

#' Pooled MLE across replicate bulk profiles sharing one composition
#'
#' Sums the per-column log-likelihoods, so the asymptotics in the number
#' of replicates are the ones that justify Wilks / chi-bar-square
#' calibration. `fixed` constrains a subset of coordinates.
#'
#' @noRd
.pooled_mle <- function(
  y_mat,
  mean_signature_matrix,
  Sigma,
  fixed = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  nms <- colnames(mean_signature_matrix)
  n_celltypes <- ncol(mean_signature_matrix)
  pooled_loglik <- function(p) {
    sum(vapply(
      seq_len(ncol(y_mat)),
      function(i) {
        loglik_multivariate(
          p,
          y_mat[, i, drop = TRUE],
          mean_signature_matrix,
          Sigma
        )
      },
      numeric(1)
    ))
  }
  if (is.null(fixed)) {
    idx <- integer(0)
    fixed_value <- numeric(0)
    free <- seq_len(n_celltypes)
    remaining <- 1
  } else {
    idx <- .resolve_celltype_index(names(fixed), nms, n_celltypes)
    fixed_value <- as.numeric(fixed)
    free <- setdiff(seq_len(n_celltypes), idx)
    remaining <- 1 - sum(fixed_value)
  }
  assemble <- function(z) {
    p <- numeric(n_celltypes)
    p[idx] <- fixed_value
    n_free <- length(free)
    if (n_free == 1L) {
      p[free] <- remaining
    } else if (n_free > 1L) {
      if (remaining <= 0) {
        p[free] <- 0
      } else {
        p[free] <- remaining * isometric_logistic(z)
      }
    }
    names(p) <- nms
    p
  }
  if (length(free) <= 1L || remaining <= 0) {
    p_hat <- assemble(numeric(0))
    return(list(
      coefficients = p_hat,
      loglik = pooled_loglik(p_hat),
      convergence = list(code = 0L, iterations = 0L, message = "closed form")
    ))
  }
  score <- function(z) {
    p <- assemble(z)
    jac <- remaining * jacobian_isometric_logistic(z)
    grad_p <- rowSums(vapply(
      seq_len(ncol(y_mat)),
      function(i) {
        gradient_loglik_unconstrained(
          p,
          y_mat[, i, drop = TRUE],
          mean_signature_matrix,
          Sigma
        )
      },
      numeric(n_celltypes)
    ))
    drop(grad_p[free] %*% jac)
  }
  fit <- stats::optim(
    par = rep(0, length(free) - 1L),
    fn = function(z) pooled_loglik(assemble(z)),
    gr = score,
    method = "BFGS",
    control = list(fnscale = -1, reltol = epsilon, maxit = itmax)
  )
  interior <- list(
    coefficients = assemble(fit$par),
    loglik = fit$value,
    convergence = list(
      code = fit$convergence,
      iterations = unname(fit$counts[["function"]]),
      message = fit$message
    )
  )
  # ILR cannot represent an exact zero. If the interior chart stops on a
  # face, compare the matching restricted fit and keep the better
  # likelihood (closed-simplex supremum).
  if (is.null(fixed)) {
    face_idx <- which(interior$coefficients < 1e-3)
    if (length(face_idx) > 0L && length(face_idx) < n_celltypes) {
      face <- .pooled_mle(
        y_mat,
        mean_signature_matrix,
        Sigma,
        fixed = stats::setNames(
          rep(0, length(face_idx)),
          nms[face_idx]
        ),
        epsilon = epsilon,
        itmax = itmax
      )
      if (is.finite(face$loglik) && face$loglik >= interior$loglik) {
        return(face)
      }
    }
  }
  interior
}

#' Simulate replicate bulk profiles from the convolution model
#'
#' @noRd
.simulate_bulk_replicates <- function(
  p,
  mean_signature_matrix,
  Sigma,
  n_replicates
) {
  mean_bulk <- drop(mean_signature_matrix %*% p)
  covariance <- .compute_global_variance(p, Sigma)
  draws <- MASS::mvrnorm(
    n = n_replicates,
    mu = mean_bulk,
    Sigma = covariance
  )
  out <- matrix(
    as.numeric(t(draws)),
    nrow = length(mean_bulk),
    ncol = n_replicates
  )
  rownames(out) <- rownames(mean_signature_matrix)
  out
}

#' Bracket-and-solve one profile-likelihood endpoint
#'
#' Returns the boundary itself when the deviance never crosses the
#' cut-off inside the bracket (a genuinely one-sided interval).
#'
#' @noRd
.profile_root <- function(deviance, lower, upper) {
  if (!is.finite(lower) || !is.finite(upper) || upper - lower < 1e-10) {
    return(if (upper - lower < 1e-10) lower else NA_real_)
  }
  f_lower <- tryCatch(deviance(lower), error = function(e) NA_real_)
  f_upper <- tryCatch(deviance(upper), error = function(e) NA_real_)
  if (anyNA(c(f_lower, f_upper))) {
    return(NA_real_)
  }
  if (f_lower * f_upper > 0) {
    # Deviance stays below the cut-off across the bracket: the profile
    # interval runs into the simplex face.
    return(if (f_lower < 0) lower else upper)
  }
  tryCatch(
    stats::uniroot(deviance, lower = lower, upper = upper, tol = 1e-8)$root,
    error = function(e) NA_real_
  )
}
