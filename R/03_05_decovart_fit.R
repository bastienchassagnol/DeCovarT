#' Fit the DeCovarT Gaussian-convolution model
#'
#' @description
#' Maximum-likelihood estimator of cellular proportions
#' \eqn{\boldsymbol{p}_{\cdot i}} under the multivariate Gaussian convolution
#' \deqn{
#'   \boldsymbol{y}_{\cdot i}\mid(\boldsymbol{\mu},\boldsymbol{\Sigma}_j,
#'   \boldsymbol{p}_{\cdot i})
#'   \sim
#'   \mathcal{N}_{G}\bigl(
#'     \boldsymbol{\mu}\boldsymbol{p}_{\cdot i},\;
#'     \boldsymbol{\Sigma}(\boldsymbol{p}_{\cdot i})
#'   \bigr),
#'   \qquad
#'   \boldsymbol{\Sigma}(\boldsymbol{p})
#'   =\sum_{j=1}^{J}p_j^{2}\boldsymbol{\Sigma}_j.
#' }
#' This is a **variance / likelihood specification**, not ordinary least
#' squares: cell-type covariances enter the residual law and cannot be
#' written as extra columns of \eqn{\boldsymbol{\mu}}. There is therefore no
#' `formula` / `lm` interface, and no `predict()` method for forecasting
#' bulk expression.
#'
#' The returned object is an S3 class `decovart_fit` with accessors
#' [coef()], [fitted()], [residuals()], [vcov()], [nobs()], [confint()],
#' [print()], [summary()] and [plot()]. Residuals are
#' \eqn{\boldsymbol{Y}-\boldsymbol{\mu}\hat{\boldsymbol{P}}}; they are **not** OLS
#' residuals. Goodness of fit is the convolution log-likelihood, not
#' residual sum of squares.
#'
#' @param signature_matrix Mean signature
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (rownames = genes,
#'   colnames = cell types).
#' @param bulk_expression Bulk matrix
#'   \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}} (rownames = genes,
#'   colnames = samples).
#' @param true_ratios Optional ground-truth proportions (\eqn{J} vector or
#'   \eqn{J\times N} matrix). Stored on the fit; not used by the MLE.
#' @param Sigma Array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}} of
#'   cell-type covariances.
#' @param method Optimiser; one of `"Marquardt-Levenberg"`,
#'   `"L-BFGS-B"`, `"Newton-Raphson"` (case-insensitive). These three
#'   maps already land on the simplex (ILR or \eqn{p/\sum p}); they do
#'   **not** call [repair_simplex()].
#' @param epsilon,itmax Absolute convergence tolerance and iteration
#'   budget, in the same roles as `reltol` / `maxit` in
#'   [stats::optim()].
#' @param standardise If `TRUE`, apply a **gene-wise** affine z-score
#'   computed once from \eqn{\boldsymbol{\mu}} to bulk, means and covariances
#'   (see Details). Cell-type-wise or sample-wise transforms are not
#'   supported.
#' @param scaled Deprecated. `TRUE` (log2 mixing) always errors.
#' @param n_starts Number of additional random Dirichlet starts per
#'   sample. With `n_starts > 0` the best log-likelihood is kept and
#'   [multistart_decovart()] records the spread of attained optima, the
#'   only direct probe of multimodality (the realised log-likelihood is
#'   not globally concave).
#' @param boundary_tol Threshold on \eqn{\min_j\hat{p}_j} below which a
#'   sample is flagged `near_boundary` by [boundary_diagnostics()]. Wald
#'   intervals are unreliable there; prefer
#'   [confint_profile_decovart()] or [lrt_decovart()].
#'
#' @details
#' **Standardisation.** CIBERSORT requires non-negative expression, no
#' missing values, and a non-log linear scale
#' \insertCite{newmanRobustEnumerationCell2015}{DeCovarT}. A logarithm is
#' concave, so Jensen's inequality shifts first and second moments and
#' \eqn{\log(\boldsymbol{\mu}\boldsymbol{p})\neq(\log\boldsymbol{\mu})\boldsymbol{p}}.
#' A gene-wise affine map
#' \eqn{\boldsymbol{x}\mapsto D^{-1}(\boldsymbol{x}-\boldsymbol{m})} with the
#' **same** \eqn{D,\boldsymbol{m}} on \eqn{\boldsymbol{Y}}, \eqn{\boldsymbol{\mu}} and
#' \eqn{\boldsymbol{\Sigma}_j^\star=D^{-1}\boldsymbol{\Sigma}_j D^{-1}} leaves the
#' MLE of \eqn{\boldsymbol{p}} unchanged (equivariance).
#'
#' **Wald covariance.** Let
#' \eqn{\boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}}.
#' The expected Fisher information of the unconstrained mean--covariance
#' map (multivariate normal; see e.g. the Wikipedia entry *Fisher
#' information*, multivariate normal) is
#' \deqn{
#'   I(\boldsymbol{p})_{jk}
#'   =
#'   \boldsymbol{\mu}_{\cdot j}^{\top}
#'   \boldsymbol{\Theta}(\boldsymbol{p})
#'   \boldsymbol{\mu}_{\cdot k}
#'   +
#'   2 p_j p_k\,
#'   \mathrm{tr}\bigl(
#'     \boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}_j
#'     \boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}_k
#'   \bigr).
#' }
#' Cramer--Rao gives
#' \eqn{\mathrm{Var}(\hat{\boldsymbol{z}})\succeq I_{\boldsymbol{z}}^{-1}}
#' in ILR coordinates, with
#' \eqn{I_{\boldsymbol{z}}
#' =\mathbf{J}_{\boldsymbol{\psi}}^{\top} I(\boldsymbol{p})
#' \mathbf{J}_{\boldsymbol{\psi}}} and
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}
#' =\mathbf{S}(\boldsymbol{p})\mathbf{V}}
#' ([jacobian_isometric_logistic()]). The delta method maps the bound
#' back to the simplex:
#' \deqn{
#'   \mathrm{Var}(\hat{\boldsymbol{p}})
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   I_{\boldsymbol{z}}^{-1}
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}.
#' }
#'
#' @return An object of class `decovart_fit`.
#'
#' @examples
#' toy <- readRDS(system.file(
#'   "extdata", "toy_deconvolution.rds",
#'   package = "DeCovarT"
#' ))
#' fit <- fit_decovart(
#'   signature_matrix = toy$signature_matrix,
#'   bulk_expression = toy$bulk_expression,
#'   Sigma = toy$Sigma,
#'   itmax = 40
#' )
#' coef(fit)
#' nobs(fit)
#' @srrstats {RE1.2} Input classes are numeric matrices / arrays.
#' @srrstats {RE1.3} Gene rownames and sample colnames of `Y` are kept on
#'   fitted values and residuals; cell-type colnames of `mu` label
#'   `coef()`.
#' @srrstats {RE1.4} Gaussian convolution, PD Sigma, simplex `p`.
#' @srrstats {RE2.0} ILR reparametrisation for Marquardt and Newton.
#' @srrstats {RE2.1} Missing values error in `.prepare_deconvolution_inputs()`.
#' @srrstats {RE2.3} Gene-wise affine `standardise`; log2 mixing errors.
#' @srrstats {RE2.4} Collinear signature columns warn (RE2.4a, RE2.4b).
#' @srrstats {RE3.0} Non-convergence of Marquardt--Levenberg warns.
#' @srrstats {RE3.2} Defaults `epsilon = 1e-4`, `itmax = 200`.
#' @srrstats {RE3.3} Both knobs are arguments (optim-style).
#' @srrstats {RE4.0} `decovart_fit` stores coefficients, log-likelihood,
#'   Fisher `vcov`, and convergence diagnostics.
#' @srrstats {RE4.2} `coef()` returns \(\hat{\boldsymbol{P}}\) (\(J\times N\)).
#' @srrstats {RE4.3} `confint()` uses the ILR delta-method Wald bound.
#' @srrstats {RE4.5} `nobs()` is \(N\), with attributes `n_genes` and
#'   `n_celltypes`.
#' @srrstats {RE4.6} `vcov()` is the Cramer--Rao / delta-method matrix.
#' @srrstats {RE4.7} Convergence codes and iteration counts are stored.
#' @srrstats {RE4.8} Observed bulk `Y` is stored as `bulk_expression`.
#' @srrstats {RE4.9} `fitted()` is \(\boldsymbol{\mu}\hat{\boldsymbol{P}}\).
#' @srrstats {RE4.10} `residuals()` is \(\boldsymbol{Y}-\hat{\boldsymbol{Y}}\);
#'   these are convolution residuals, not OLS residuals. Goodness of fit
#'   is the MLE log-likelihood (RE4.11), not residual sum of squares.
#' @srrstats {RE4.11} `summary()` reports log-likelihood (and AIC).
#' @srrstats {RE4.12} ILR maps are [isometric_logistic()] /
#'   [isometric_log_ratio()] with [helmert_basis()].
#' @srrstats {RE4.13} Signature and `Sigma` are stored on the fit.
#' @srrstats {RE4.17} `print.decovart_fit()` shows proportions and
#'   log-likelihood.
#' @srrstats {RE4.18} `summary.decovart_fit()` adds Wald SEs and
#'   convergence.
#' @srrstats {RE6.0} `plot.decovart_fit()` compares observed and fitted
#'   bulk expression (not estimated vs true proportions).
#' @srrstats {RE6.1} Dispatch is `plot.decovart_fit`.
#' @srrstats {RE6.2} Default plot is observed vs fitted bulk profiles.
#' @srrstats {RE7.2} Dimnames of `Y` and `mu` are retained.
#' @srrstats {RE7.3} Tests exercise `coef`, `fitted`, `residuals`,
#'   `vcov`, `nobs`, `print`, `summary`, `plot`.
#' @references
#' \insertRef{newmanRobustEnumerationCell2015}{DeCovarT}
#' @importFrom stats nobs coef vcov fitted residuals confint
#' @export
#' @family decovart_fit
#' @seealso [deconvolute_ratios()],
#'   [deconvolute_ratios_Marquardt_Levenberg()],
#'   [expected_fisher_unconstrained()], [vcov_ilr_delta()],
#'   [coef.decovart_fit()], [vcov.decovart_fit()],
#'   [confint.decovart_fit()]
fit_decovart <- function(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  method = c(
    "Marquardt-Levenberg",
    "L-BFGS-B",
    "Newton-Raphson"
  ),
  epsilon = 10^-4,
  itmax = 200,
  standardise = FALSE,
  scaled = FALSE,
  n_starts = 0L,
  boundary_tol = 1e-8
) {
  method <- .match_arg_case_insensitive(
    method,
    c("Marquardt-Levenberg", "L-BFGS-B", "Newton-Raphson")
  )
  if (is.null(Sigma)) {
    stop(
      "`Sigma` is required: DeCovarT is a variance model, not OLS.",
      call. = FALSE
    )
  }
  aligned <- .prepare_deconvolution_inputs(
    signature_matrix = signature_matrix,
    bulk_expression = bulk_expression,
    true_ratios = true_ratios,
    Sigma = Sigma,
    standardise = standardise,
    scaled = scaled
  )
  mu <- aligned$signature_matrix
  y_mat <- aligned$bulk_expression
  sigma_arr <- aligned$Sigma
  solver <- switch(
    method,
    "Marquardt-Levenberg" = deconvolute_ratios_Marquardt_Levenberg,
    "L-BFGS-B" = deconvolute_ratios_L_BFGS_B,
    "Newton-Raphson" = deconvolute_ratios_Newton_Raphson
  )
  n_samples <- ncol(y_mat)
  n_celltypes <- ncol(mu)
  coef_mat <- matrix(
    NA_real_,
    nrow = n_celltypes,
    ncol = n_samples,
    dimnames = list(colnames(mu), colnames(y_mat))
  )
  loglik <- stats::setNames(rep(NA_real_, n_samples), colnames(y_mat))
  vcov_list <- vector("list", n_samples)
  names(vcov_list) <- colnames(y_mat)
  convergence <- vector("list", n_samples)
  names(convergence) <- colnames(y_mat)
  diagnostics <- vector("list", n_samples)
  names(diagnostics) <- colnames(y_mat)
  n_starts <- as.integer(n_starts)
  for (i in seq_len(n_samples)) {
    y_i <- y_mat[, i, drop = TRUE]
    fit_i <- solver(
      y = y_i,
      mean_signature_matrix = mu,
      Sigma = sigma_arr,
      epsilon = epsilon,
      itmax = itmax,
      return_model = TRUE
    )
    convergence[[i]] <- fit_i$convergence
    if (n_starts > 0L) {
      restarts <- multistart_decovart(
        y = y_i,
        mean_signature_matrix = mu,
        Sigma = sigma_arr,
        n_starts = n_starts,
        solver = solver,
        epsilon = epsilon,
        itmax = itmax
      )
      if (restarts$loglik > fit_i$loglik) {
        fit_i$coefficients <- restarts$coefficients
        fit_i$loglik <- restarts$loglik
      }
      convergence[[i]]$loglik_range <- restarts$loglik_range
      convergence[[i]]$multimodal <- restarts$multimodal
    }
    coef_mat[, i] <- fit_i$coefficients
    loglik[[i]] <- fit_i$loglik
    vcov_list[[i]] <- vcov_ilr_delta(
      fit_i$coefficients,
      mu,
      sigma_arr
    )
    diagnostics[[i]] <- boundary_diagnostics(
      fit_i$coefficients,
      y_i,
      mu,
      sigma_arr,
      boundary_tol = boundary_tol
    )
  }
  near_boundary <- vapply(
    diagnostics,
    function(d) isTRUE(d$near_boundary),
    logical(1)
  )
  if (any(near_boundary)) {
    warning(
      "Estimated proportions reach a simplex face in sample(s) ",
      toString(names(diagnostics)[near_boundary]),
      ". This may be a genuine boundary optimum, but ILR Wald intervals ",
      "are then undefined; use confint_profile_decovart() or ",
      "lrt_decovart().",
      call. = FALSE
    )
  }
  fitted_mat <- mu %*% coef_mat
  dimnames(fitted_mat) <- dimnames(y_mat)
  structure(
    list(
      coefficients = coef_mat,
      fitted.values = fitted_mat,
      residuals = y_mat - fitted_mat,
      vcov = vcov_list,
      loglik = loglik,
      convergence = convergence,
      diagnostics = diagnostics,
      method = method,
      epsilon = epsilon,
      itmax = itmax,
      standardise = isTRUE(standardise),
      centre = aligned$centre,
      scale = aligned$scale,
      n_genes = nrow(mu),
      n_celltypes = n_celltypes,
      n_samples = n_samples,
      signature_matrix = mu,
      bulk_expression = y_mat,
      Sigma = sigma_arr,
      true_ratios = aligned$true_ratios
    ),
    class = "decovart_fit"
  )
}

#' Expected Fisher information of unconstrained \eqn{\boldsymbol{p}}
#'
#' @description
#' For the multivariate-normal mean--covariance map of the DeCovarT
#' convolution,
#' \eqn{\boldsymbol{y}\sim\mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))} with
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^{2}\boldsymbol{\Sigma}_j}
#' and precision
#' \eqn{\boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}},
#' the expected Fisher information has entries
#' \deqn{
#'   I(\boldsymbol{p})_{jk}
#'   =
#'   \boldsymbol{\mu}_{\cdot j}^{\top}
#'   \boldsymbol{\Theta}(\boldsymbol{p})
#'   \boldsymbol{\mu}_{\cdot k}
#'   +
#'   2 p_j p_k\,
#'   \mathrm{tr}\bigl(
#'     \boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}_j
#'     \boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}_k
#'   \bigr).
#' }
#' The first summand is the mean contribution
#' (an \eqn{\boldsymbol{\Theta}}-inner product of signature columns); the
#' second is the covariance contribution of the quadratic map
#' \eqn{\boldsymbol{p}\mapsto\boldsymbol{\Sigma}(\boldsymbol{p})}. See the
#' multivariate-normal formula on
#' <https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution>.
#'
#' @param p Numeric proportions on the open simplex.
#' @param mean_signature_matrix Mean signature \eqn{\boldsymbol{\mu}}
#'   (\eqn{G\times J}).
#' @param Sigma Cell-type covariances \eqn{G\times G\times J}.
#'
#' @return Symmetric \eqn{J\times J} expected Fisher information matrix
#'   \eqn{I(\boldsymbol{p})}.
#'
#' @seealso [vcov_ilr_delta()], [vcov.decovart_fit()],
#'   [confint.decovart_fit()], [.inner_product()]
#'
#' @keywords internal
#' @export
expected_fisher_unconstrained <- function(
  p,
  mean_signature_matrix,
  Sigma
) {
  n_celltypes <- length(p)
  theta <- .sigma_p_factorisation(p, Sigma)$inverse
  info <- matrix(0, n_celltypes, n_celltypes)
  for (j in seq_len(n_celltypes)) {
    mu_j <- mean_signature_matrix[, j, drop = TRUE]
    sigma_j <- Sigma[,, j]
    for (k in seq_len(n_celltypes)) {
      mu_k <- mean_signature_matrix[, k, drop = TRUE]
      sigma_k <- Sigma[,, k]
      mean_term <- .inner_product(mu_j, theta, mu_k)
      cov_term <- 2 *
        p[[j]] *
        p[[k]] *
        sum(diag(theta %*% sigma_j %*% theta %*% sigma_k))
      info[j, k] <- mean_term + cov_term
    }
  }
  info
}

#' Cramer--Rao / ILR delta-method covariance of \eqn{\hat{\boldsymbol{p}}}
#'
#' @description
#' Maps the expected Fisher information of unconstrained proportions
#' through the isometric log-ratio (ILR) chart and back to the simplex
#' via the delta method.
#'
#' Let \eqn{\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{z})}
#' with Jacobian
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}=\mathbf{S}(\boldsymbol{p})
#' \mathbf{V}}
#' ([jacobian_isometric_logistic()]). Fisher information transforms as
#' \deqn{
#'   I_{\boldsymbol{z}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   I(\boldsymbol{p})
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   =
#'   \mathbf{V}^{\mathsf{T}}
#'   \mathbf{S}(\boldsymbol{p})
#'   I(\boldsymbol{p})
#'   \mathbf{S}(\boldsymbol{p})
#'   \mathbf{V}.
#' }
#' The first-order delta method then yields
#' \deqn{
#'   \mathrm{Var}(\hat{\boldsymbol{p}})
#'   \approx
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   I_{\boldsymbol{z}}^{-1}
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}.
#' }
#' This simplex covariance is invariant to orthogonal rotations of
#' \eqn{\mathbf{V}}. The construction is undefined on the simplex
#' boundary (the log-ratio chart blows up); the function then returns
#' `NA` with a warning. [vcov_alr_delta()] is the ALR-chart analogue
#' used only for reference-invariance checks.
#'
#' @inheritParams expected_fisher_unconstrained
#'
#' @return Symmetric \eqn{J\times J} asymptotic covariance of
#'   \eqn{\hat{\boldsymbol{p}}}, or a matrix of `NA` if the bound is
#'   undefined / singular.
#'
#' @seealso [expected_fisher_unconstrained()], [vcov.decovart_fit()],
#'   [confint.decovart_fit()], [jacobian_isometric_logistic()],
#'   [vcov_alr_delta()]
#'
#' @keywords internal
#' @export
vcov_ilr_delta <- function(p, mean_signature_matrix, Sigma) {
  nms <- names(p)
  n_celltypes <- length(p)
  out <- matrix(
    NA_real_,
    n_celltypes,
    n_celltypes,
    dimnames = list(nms, nms)
  )
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    warning(
      "Proportions on the simplex boundary; Wald vcov is undefined.",
      call. = FALSE
    )
    return(out)
  }
  info_p <- expected_fisher_unconstrained(p, mean_signature_matrix, Sigma)
  z <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(z)
  info_z <- t(jac) %*% info_p %*% jac
  vcov_z <- tryCatch(
    solve(info_z),
    error = function(e) {
      warning(
        "Expected Fisher information in ILR coordinates is singular.",
        call. = FALSE
      )
      NULL
    }
  )
  if (is.null(vcov_z)) {
    return(out)
  }
  vcov_p <- jac %*% vcov_z %*% t(jac)
  dimnames(vcov_p) <- list(nms, nms)
  vcov_p
}

#' Cramer--Rao / ALR delta-method covariance of \eqn{\hat{\boldsymbol{p}}}
#'
#' @description
#' Maps the expected Fisher information of unconstrained proportions
#' through the additive log-ratio (ALR) chart and back to the simplex
#' via the delta method.
#'
#' Let \eqn{\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})} with
#' Jacobian
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}
#' =\partial\boldsymbol{\psi}/\partial\boldsymbol{\rho}^{\top}}
#' ([jacobian_additive_logistic()]). Fisher information transforms as
#' the covariant quadratic form
#' \deqn{
#'   I_{\boldsymbol{\rho}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\top}
#'   I(\boldsymbol{p})
#'   \mathbf{J}_{\boldsymbol{\psi}}.
#' }
#' Under a regular large-sample regime the MLE in ALR coordinates is
#' asymptotically normal,
#' \eqn{\hat{\boldsymbol{\rho}}
#' \overset{a}{\sim}
#' \mathcal{N}(\boldsymbol{\rho}_{0}, I_{\boldsymbol{\rho}}^{-1})}.
#' The first-order delta method then yields the same simplex covariance
#' as [vcov_ilr_delta()] when both charts are transformed correctly:
#' \deqn{
#'   \mathrm{Var}(\hat{\boldsymbol{p}})
#'   \approx
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   I_{\boldsymbol{\rho}}^{-1}
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\top}.
#' }
#' This helper is kept for reference-invariance checks; [vcov.decovart_fit()]
#' uses [vcov_ilr_delta()]. The construction is undefined on the simplex
#' boundary (the ALR chart blows up); the function then returns `NA` with
#' a warning. See also
#' <https://en.wikipedia.org/wiki/Delta_method> and
#' <https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution>.
#'
#' @inheritParams expected_fisher_unconstrained
#'
#' @return Symmetric \eqn{J\times J} asymptotic covariance of
#'   \eqn{\hat{\boldsymbol{p}}}, or a matrix of `NA` if the bound is
#'   undefined / singular.
#'
#' @seealso [vcov_ilr_delta()], [jacobian_additive_logistic()]
#'
#' @keywords internal
#' @export
vcov_alr_delta <- function(p, mean_signature_matrix, Sigma) {
  nms <- names(p)
  n_celltypes <- length(p)
  out <- matrix(
    NA_real_,
    n_celltypes,
    n_celltypes,
    dimnames = list(nms, nms)
  )
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    warning(
      "Proportions on the simplex boundary; Wald vcov is undefined.",
      call. = FALSE
    )
    return(out)
  }
  info_p <- expected_fisher_unconstrained(p, mean_signature_matrix, Sigma)
  rho <- additive_log_ratio(p)
  jac <- jacobian_additive_logistic(rho)
  info_rho <- t(jac) %*% info_p %*% jac
  vcov_rho <- tryCatch(
    solve(info_rho),
    error = function(e) {
      warning(
        "Expected Fisher information in ALR coordinates is singular.",
        call. = FALSE
      )
      NULL
    }
  )
  if (is.null(vcov_rho)) {
    return(out)
  }
  vcov_p <- jac %*% vcov_rho %*% t(jac)
  dimnames(vcov_p) <- list(nms, nms)
  vcov_p
}

#' @rdname fit_decovart
#' @param x,object A `decovart_fit`.
#' @param ... Passed to [graphics::plot()] (plot method) or ignored.
#' @export
#' @method print decovart_fit
print.decovart_fit <- function(x, ...) {
  cat(
    "DeCovarT fit (",
    x$method,
    ")\n",
    "  genes: ",
    x$n_genes,
    ", cell types: ",
    x$n_celltypes,
    ", samples: ",
    x$n_samples,
    "\n",
    sep = ""
  )
  cat("Proportions:\n")
  print(x$coefficients, ...)
  cat("Log-likelihood:", toString(signif(x$loglik, 6)), "\n")
  invisible(x)
}

#' @rdname fit_decovart
#' @export
#' @method summary decovart_fit
summary.decovart_fit <- function(object, ...) {
  se <- vapply(
    object$vcov,
    function(v) {
      d <- diag(v)
      d[!is.finite(d)] <- NA_real_
      sqrt(pmax(d, 0))
    },
    numeric(object$n_celltypes)
  )
  if (is.null(dim(se))) {
    se <- matrix(
      se,
      ncol = 1L,
      dimnames = list(rownames(object$coefficients), NULL)
    )
  } else {
    dimnames(se) <- dimnames(object$coefficients)
  }
  n_par <- object$n_celltypes - 1L
  aic <- -2 * object$loglik + 2 * n_par
  structure(
    list(
      coefficients = object$coefficients,
      se = se,
      loglik = object$loglik,
      aic = aic,
      method = object$method,
      convergence = object$convergence,
      n_genes = object$n_genes,
      n_celltypes = object$n_celltypes,
      n_samples = object$n_samples
    ),
    class = "summary.decovart_fit"
  )
}

#' @rdname fit_decovart
#' @export
#' @method print summary.decovart_fit
print.summary.decovart_fit <- function(x, ...) {
  cat(
    "DeCovarT summary (",
    x$method,
    ")\nGoodness of fit is the convolution MLE ",
    "log-likelihood, not least squares.\n",
    sep = ""
  )
  for (i in seq_len(x$n_samples)) {
    cat("\nSample ", colnames(x$coefficients)[[i]], "\n", sep = "")
    tab <- cbind(
      Estimate = x$coefficients[, i],
      `Std. Error` = x$se[, i]
    )
    print(tab, ...)
    cat(
      "  logLik = ",
      signif(x$loglik[[i]], 6),
      ", AIC = ",
      signif(x$aic[[i]], 6),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method coef decovart_fit
coef.decovart_fit <- function(object, ...) {
  object$coefficients
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method fitted decovart_fit
fitted.decovart_fit <- function(object, ...) {
  object$fitted.values
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method residuals decovart_fit
residuals.decovart_fit <- function(object, ...) {
  object$residuals
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method vcov decovart_fit
vcov.decovart_fit <- function(object, ...) {
  if (object$n_samples == 1L) {
    return(object$vcov[[1L]])
  }
  object$vcov
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method nobs decovart_fit
nobs.decovart_fit <- function(object, ...) {
  n <- object$n_samples
  attr(n, "n_genes") <- object$n_genes
  attr(n, "n_celltypes") <- object$n_celltypes
  attr(n, "n_samples") <- object$n_samples
  n
}

#' @rdname fit_decovart
#' @param parm Kept for compatibility with [stats::confint()]; all
#'   simplex coordinates are always returned (subsetting by `parm` is
#'   not implemented). Removing this formal would break the S3 method
#'   contract with the generic.
#' @param level Confidence level \eqn{1-\alpha} (default `0.95`). Wald
#'   intervals use asymptotic normality of the MLE with standard errors
#'   from [vcov_ilr_delta()] /
#'   [expected_fisher_unconstrained()] (see Details of [fit_decovart()]
#'   and of [vcov_ilr_delta()]):
#'   \eqn{\hat{p}_j \pm z_{1-\alpha/2}\,\mathrm{SE}_j}.
#' @export
#' @family decovart_fit
#' @method confint decovart_fit
confint.decovart_fit <- function(object, parm, level = 0.95, ...) {
  # `parm` is required by stats::confint; subsetting is not implemented.
  if (!missing(parm)) {
    invisible(parm)
  }
  alpha <- (1 - level) / 2
  z <- stats::qnorm(c(alpha, 1 - alpha))
  cf <- object$coefficients
  out <- vector("list", object$n_samples)
  names(out) <- colnames(cf)
  for (i in seq_len(object$n_samples)) {
    se <- sqrt(pmax(diag(object$vcov[[i]]), 0))
    interval <- cbind(cf[, i] + z[[1L]] * se, cf[, i] + z[[2L]] * se)
    colnames(interval) <- paste0(
      format(100 * c(alpha, 1 - alpha), trim = TRUE),
      "%"
    )
    rownames(interval) <- rownames(cf)
    out[[i]] <- interval
  }
  if (object$n_samples == 1L) {
    return(out[[1L]])
  }
  out
}

#' @rdname fit_decovart
#' @export
#' @family decovart_fit
#' @method plot decovart_fit
plot.decovart_fit <- function(x, ...) {
  y_obs <- as.vector(x$bulk_expression)
  y_hat <- as.vector(x$fitted.values)
  graphics::plot(
    y_obs,
    y_hat,
    xlab = "Observed bulk expression",
    ylab = "Fitted bulk expression",
    ...
  )
  graphics::abline(0, 1, lty = 2)
  invisible(x)
}
