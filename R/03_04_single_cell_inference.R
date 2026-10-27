# Intervals for the single-cell mixture. The composition interval uses the
# profile log-likelihood (the Schur complement that removes c and phi).
# The two scales use the other Schur complement. The default map from
# the ILR coordinate to the simplex is the same inverse chart as
# confint_decovart(). Higher-order and sandwich corrections for a
# plug-in covariance are not implemented.
# See vignette("theory-population-inference") and
# vignette("theory-single-cell-generative-optimisation") and
# vignette("theory-single-cell-inference").

#' @noRd
.single_cell_profile_value <- function(p, y, mean_signature_matrix, Sigma) {
  scales <- .single_cell_scales_at_p(
    y,
    mean_signature_matrix,
    Sigma,
    p,
    strict = FALSE
  )
  if (is.null(scales)) {
    return(-Inf)
  }
  loglik_single_cell(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    scales$mean_scale,
    scales$variance_scale
  )
}

#' Maximum likelihood on a face of the single-cell model
#'
#' @description
#' Profiles the two technical scales out. Without `fixed`, this is
#' [deconvolute_single_cell()]. With `fixed`, the free block is an ILR
#' chart of the remaining face, so a coordinate can be held at zero, and
#' it is fitted by Newton-Raphson on the profile log-likelihood with the
#' PORT trust region of [stats::nlminb()], as [restricted_mle_decovart()].
#'
#' @inheritParams deconvolute_single_cell
#' @param fixed Named ratios held fixed, as in [restricted_mle_decovart()].
#'
#' @return A list with `coefficients` (the fitted composition), `loglik`,
#'   `fixed` and `convergence`.
#'
#' @keywords internal
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' restricted_mle_single_cell(y, mu, sigma, fixed = c(ct2 = 0.3))$loglik
restricted_mle_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  fixed = NULL,
  epsilon = 1e-8,
  itmax = 100L
) {
  checked <- .check_single_cell_inputs(y, mean_signature_matrix, Sigma)
  y <- checked$y
  mean_signature_matrix <- checked$mean
  nms <- colnames(mean_signature_matrix)
  if (is.null(nms)) {
    nms <- paste0("ct", seq_len(checked$n_types))
    colnames(mean_signature_matrix) <- nms
  }
  if (length(fixed) == 0L) {
    fit <- deconvolute_single_cell(
      y,
      mean_signature_matrix,
      Sigma,
      epsilon = epsilon,
      itmax = itmax,
      profile = TRUE
    )
    return(list(
      coefficients = stats::setNames(fit$composition, nms),
      loglik = fit$loglik,
      fixed = stats::setNames(numeric(0), character(0)),
      convergence = list(
        code = fit$convergence,
        iterations = fit$iterations,
        message = fit$message,
        rdm = NA_real_
      )
    ))
  }
  at_scales <- function(p, fun) {
    scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
    fun(
      p,
      y,
      mean_signature_matrix,
      Sigma,
      scales$mean_scale,
      scales$variance_scale
    )
  }
  .fit_on_face(
    nms,
    fixed,
    loglik = function(p) {
      .single_cell_profile_value(p, y, mean_signature_matrix, Sigma)
    },
    score = function(p) at_scales(p, gradient_loglik_single_cell),
    hessian = function(p) at_scales(p, .profile_hessian_single_cell),
    epsilon = epsilon,
    itmax = itmax
  )
}

#' Small-sample calibration of the single-cell intervals
#'
#' `"ml"` keeps the maximum-likelihood variance scale and normal
#' quantiles. `"reml"` counts the \eqn{J} mean parameters (the scale
#' \eqn{c} and \eqn{J-1} composition coordinates): the residual degrees
#' of freedom are \eqn{\nu=G-J}, the covariances are inflated by
#' \eqn{G/\nu} (the REML divisor of \eqn{\hat\phi}), and quantiles come from
#' \eqn{t_{\nu}}, as in a linear model with unknown variance.
#'
#' @noRd
.single_cell_calibration <- function(variance, n_genes, n_types) {
  if (identical(variance, "reml")) {
    df <- n_genes - n_types
    list(variance = "reml", df = df, inflation = n_genes / df)
  } else {
    list(variance = "ml", df = Inf, inflation = 1)
  }
}

#' Platform scale and number of cells from the two fitted scales
#'
#' \eqn{c=a\,n} and \eqn{\phi=a^{2}n} give \eqn{a=\phi/c} and
#' \eqn{n=c^{2}/\phi}; \eqn{n} is a number of cells only when \eqn{a} is
#' known. Rows `platform_scale_a` and `cell_count_n` are appended to a
#' vector or to a matrix with rows `mean_scale` and `variance_scale`.
#' @noRd
.derived_scales <- function(scales) {
  if (is.null(dim(scales))) {
    return(c(
      scales,
      platform_scale_a = scales[["variance_scale"]] / scales[["mean_scale"]],
      cell_count_n = scales[["mean_scale"]]^2 / scales[["variance_scale"]]
    ))
  }
  rbind(
    scales,
    platform_scale_a = scales["variance_scale", ] / scales["mean_scale", ],
    cell_count_n = scales["mean_scale", ]^2 / scales["variance_scale", ]
  )
}

#' @noRd
.nuisance_intervals <- function(uncertainty, level, calibration, n_genes) {
  estimate <- c(
    mean_scale = uncertainty$mean_scale,
    variance_scale = uncertainty$variance_scale
  )
  vcov <- as.matrix(uncertainty$vcov_nuisance) * calibration$inflation
  se <- sqrt(pmax(diag(vcov), 0))
  names(se) <- names(estimate)
  df <- calibration$df
  interval_wald <- .wald_bounds(estimate, se, level, df = df)
  log_estimate <- log(estimate)
  names(log_estimate) <- names(estimate)
  interval <- exp(.wald_bounds(log_estimate, se / estimate, level, df = df))
  if (is.finite(df)) {
    # G phi_hat / phi ~ chi^2_df, the exact interval of a linear model.
    alpha <- (1 - level) / 2
    interval["variance_scale", ] <- n_genes *
      estimate[["variance_scale"]] /
      stats::qchisq(c(1 - alpha, alpha), df = df)
  }
  rownames(interval) <- names(estimate)
  # log a = log phi - log c and log n = 2 log c - log phi: delta method on
  # the log scale, so both bounds stay positive.
  log_vcov <- vcov / tcrossprod(estimate)
  contrasts <- rbind(
    platform_scale_a = c(-1, 1),
    cell_count_n = c(2, -1)
  )
  derived <- .derived_scales(estimate)[rownames(contrasts)]
  derived_se <- sqrt(pmax(
    rowSums((contrasts %*% log_vcov) * contrasts),
    0
  ))
  interval <- rbind(
    interval,
    exp(.wald_bounds(log(derived), derived_se, level, df = df))
  )
  estimate <- c(estimate, derived)
  list(
    estimate = estimate,
    vcov = vcov,
    interval = interval,
    interval_wald = interval_wald
  )
}

#' @noRd
.composition_intervals <- function(
  composition,
  y,
  mean_signature_matrix,
  Sigma,
  level,
  how,
  n_draws,
  seed,
  calibration
) {
  uncertainty <- vcov_single_cell(
    composition,
    y,
    mean_signature_matrix,
    Sigma,
    warn = FALSE
  )
  uncertainty$vcov_rho <- uncertainty$vcov_rho * calibration$inflation
  uncertainty$vcov <- uncertainty$vcov * calibration$inflation
  nuisance <- .nuisance_intervals(
    uncertainty,
    level,
    calibration,
    n_genes = length(y)
  )
  if (anyNA(uncertainty$vcov_rho)) {
    interval <- .empty_proportion_interval(composition, level)
    mapped <- list(
      interval = interval,
      rho = uncertainty$rho,
      rho_interval = NULL,
      mapping = NA_character_,
      n_draws = NA_integer_
    )
  } else if (identical(how, "ilr")) {
    mapped <- .map_wald_rho(
      rho = uncertainty$rho,
      vcov_rho = uncertainty$vcov_rho,
      p_names = names(composition),
      level = level,
      n_draws = n_draws,
      seed = seed,
      df = calibration$df
    )
  } else {
    se <- sqrt(pmax(diag(uncertainty$vcov), 0))
    names(se) <- names(composition)
    mapped <- list(
      interval = .wald_bounds(composition, se, level, df = calibration$df),
      rho = uncertainty$rho,
      rho_interval = NULL,
      mapping = "delta",
      n_draws = NA_integer_
    )
  }
  list(
    estimate = composition,
    interval = mapped$interval,
    rho = mapped$rho,
    rho_interval = mapped$rho_interval,
    vcov = uncertainty$vcov,
    vcov_rho = uncertainty$vcov_rho,
    mapping = mapped$mapping,
    n_draws = mapped$n_draws,
    variance = calibration$variance,
    df = calibration$df,
    nuisance = nuisance
  )
}

# Profile intervals and likelihood-ratio tests ----------------------------

#' Profile-likelihood intervals for one single-cell bulk sample
#'
#' @description
#' Inverts the likelihood-ratio test for one coordinate of the profile
#' log-likelihood. The two scales are concentrated out at every trial
#' value. An endpoint that lands on a face is not a Wilks interval;
#' use [lrt_single_cell()] for a boundary null.
#'
#' @inheritParams restricted_mle_single_cell
#' @param level Confidence level.
#' @param celltypes Cell types to profile. Defaults to every column.
#' @param variance `"ml"` (default) uses the \eqn{\chi^2_1} cut-off.
#'   `"reml"` uses \eqn{G\log\{1+F_{1,\nu;1-\alpha}/\nu\}} with
#'   \eqn{\nu=G-J}, the exact cut-off of the profile likelihood ratio
#'   when a variance scale is profiled out of a linear model.
#'
#' @return A matrix with columns `estimate`, `lower` and `upper`.
#'
#' @seealso [confint_single_cell()], [lrt_single_cell()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' confint_profile_single_cell(y, mu, sigma)
confint_profile_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  level = 0.95,
  celltypes = colnames(mean_signature_matrix),
  variance = c("ml", "reml"),
  epsilon = 1e-8,
  itmax = 100L
) {
  level <- .check_level(level)
  variance <- match.arg(variance)
  if (is.null(celltypes)) {
    stop("`mean_signature_matrix` needs column names.", call. = FALSE)
  }
  full <- restricted_mle_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  n_genes <- length(y)
  df <- n_genes - ncol(mean_signature_matrix)
  # With phi profiled out, 2 * log LR = G log(1 + F / df) for an F(1, df)
  # statistic in a linear model; "reml" uses that exact cut-off.
  cutoff <- if (identical(variance, "reml")) {
    n_genes * log1p(stats::qf(level, 1, df) / df)
  } else {
    stats::qchisq(level, df = 1L)
  }
  intervals <- vapply(
    celltypes,
    function(celltype) {
      p_hat <- full$coefficients[[celltype]]
      deviance <- function(value) {
        restricted <- restricted_mle_single_cell(
          y,
          mean_signature_matrix,
          Sigma,
          fixed = stats::setNames(value, celltype),
          epsilon = epsilon,
          itmax = itmax
        )
        2 * (full$loglik - restricted$loglik) - cutoff
      }
      c(
        estimate = p_hat,
        lower = .profile_root(deviance, estimate = p_hat, face = 0),
        upper = .profile_root(deviance, estimate = p_hat, face = 1)
      )
    },
    numeric(3)
  )
  t(intervals)
}

#' Likelihood-ratio test for one single-cell composition
#'
#' @description
#' Compares the profiled maximum with the maximum on `fixed`
#' (passed as `null_value`). An interior null uses Wilks'
#' chi-square. A null on a face uses the chi-bar-square mixture of
#' [chi_bar_square_pvalue()]. The two technical scales are profiled
#' out under both hypotheses, so the degrees of freedom count
#' composition constraints only.
#'
#' @inheritParams restricted_mle_single_cell
#' @param null_value Named hypothesised ratios.
#' @param n_boundary,boundary_tol As in [lrt_decovart()].
#'
#' @return A one-row data frame, as for [lrt_decovart()].
#'
#' @seealso [confint_profile_single_cell()], [lrt_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(1, 0)) + c(0.2, -0.1, 0.05, 0.05)
#' lrt_single_cell(y, mu, sigma, null_value = c(ct2 = 0))
lrt_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  null_value,
  n_boundary = NULL,
  boundary_tol = 1e-6,
  epsilon = 1e-8,
  itmax = 100L
) {
  full <- restricted_mle_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  null <- restricted_mle_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    fixed = null_value,
    epsilon = epsilon,
    itmax = itmax
  )
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

# Firth penalty on the profile ---------------------------------------------

#' @noRd
.jeffreys_penalty_single_cell <- function(z, y, mean_signature_matrix, Sigma) {
  h <- tryCatch(
    hessian_profile_loglik_single_cell_ilr(
      z,
      y,
      mean_signature_matrix,
      Sigma
    ),
    error = function(e) NULL
  )
  if (is.null(h)) {
    return(-Inf)
  }
  ev <- eigen(-.symmetrise(h), symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(ev)) || min(ev) <= 1e-8) {
    return(-Inf)
  }
  0.5 * sum(log(ev))
}

#' Firth fit of one single-cell bulk sample
#'
#' @description
#' Maximises the profile log-likelihood plus half the log determinant
#' of the observed profile information. The penalty uses the negative
#' profile Hessian, not an expected Fisher matrix. Scales stay
#' concentrated out, so the penalty does not add a third free scale.
#'
#' @inheritParams deconvolute_single_cell
#'
#' @return A list with `coefficients`, `rho`, `loglik`, `penalty`,
#'   `penalised_loglik`, `mean_scale`, `variance_scale` and `convergence`.
#'
#' @seealso [confint_single_cell()], [fit_firth_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' fit_firth_single_cell(y, mu, sigma)$coefficients
fit_firth_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 1e-8,
  itmax = 80L
) {
  start <- deconvolute_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax,
    profile = TRUE
  )
  z0 <- isometric_log_ratio(start$composition)
  penalty <- function(z) {
    .jeffreys_penalty_single_cell(z, y, mean_signature_matrix, Sigma)
  }
  objective <- function(z) {
    ll <- profile_loglik_single_cell_ilr(
      z,
      y,
      mean_signature_matrix,
      Sigma
    )
    pen <- penalty(z)
    if (!is.finite(ll) || !is.finite(pen)) {
      return(-1e10)
    }
    ll + pen
  }
  gradient <- function(z) {
    g_ll <- gradient_profile_loglik_single_cell_ilr(
      z,
      y,
      mean_signature_matrix,
      Sigma
    )
    g_pen <- .central_gradient(penalty, z)
    if (any(!is.finite(g_pen))) {
      g_pen <- rep(0, length(z))
    }
    g_ll + g_pen
  }
  fit <- stats::optim(
    par = z0,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      fnscale = -1,
      reltol = epsilon,
      maxit = as.integer(itmax)
    )
  )
  nms <- colnames(mean_signature_matrix)
  p_hat <- repair_simplex(isometric_logistic(fit$par), open = TRUE, nms = nms)
  scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p_hat)
  list(
    coefficients = p_hat,
    rho = stats::setNames(fit$par, paste0("rho", seq_along(fit$par))),
    loglik = loglik_single_cell(
      p_hat,
      y,
      mean_signature_matrix,
      Sigma,
      scales$mean_scale,
      scales$variance_scale
    ),
    penalty = penalty(fit$par),
    penalised_loglik = fit$value,
    mean_scale = scales$mean_scale,
    variance_scale = scales$variance_scale,
    convergence = fit$convergence
  )
}

# Resampling ---------------------------------------------------------------

#' @noRd
.simulate_single_cell_bulk <- function(
  p,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  mean_bulk <- mean_scale * drop(mean_signature_matrix %*% p)
  covariance <- variance_scale *
    .mixture_factorisation(as.numeric(p), Sigma)$matrix
  stats::setNames(
    as.numeric(MASS::mvrnorm(1L, mu = mean_bulk, Sigma = covariance)),
    rownames(mean_signature_matrix)
  )
}

#' Parametric bootstrap for one single-cell bulk sample
#'
#' @description
#' Draws
#' \eqn{\boldsymbol{Y}^{(b)}\sim\mathcal{N}(\hat c\tilde{\boldsymbol{\mu}}
#' \hat{\boldsymbol{p}},\hat\phi\boldsymbol{S}_W(\hat{\boldsymbol{p}}))}
#' and refits the profile model. Percentile intervals are returned for
#' the composition and for the two scales. The reference moments stay
#' fixed.
#'
#' @inheritParams deconvolute_single_cell
#' @param n_boot Number of replicates.
#' @param level Confidence level.
#' @param optimiser Passed to [deconvolute_single_cell()] for the fit and
#'   every refit, so the bootstrap matches the estimator being assessed.
#'
#' @return A list with `estimates`, `interval`, `n_failed` (refits that
#'   errored, left as `NA` columns), `scales`, `scale_interval` and `fit`.
#'   `scales` and `scale_interval` have rows `mean_scale` (\eqn{c}),
#'   `variance_scale` (\eqn{\phi}), `platform_scale_a`
#'   (\eqn{a=\phi/c}) and `cell_count_n` (\eqn{n=c^{2}/\phi}).
#'
#' @seealso [confint_single_cell()], [bootstrap_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' set.seed(1)
#' bootstrap_single_cell(y, mu, sigma, n_boot = 5)$interval
bootstrap_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  n_boot = 200L,
  level = 0.95,
  epsilon = 1e-6,
  itmax = 40L,
  optimiser = c("newton", "lbfgsb")
) {
  optimiser <- match.arg(optimiser)
  level <- .check_level(level)
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  fit <- deconvolute_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax,
    profile = TRUE,
    optimiser = optimiser
  )
  nms <- names(fit$composition)
  estimates <- matrix(
    NA_real_,
    nrow = length(fit$composition),
    ncol = n_boot,
    dimnames = list(nms, NULL)
  )
  scales <- matrix(
    NA_real_,
    nrow = 2L,
    ncol = n_boot,
    dimnames = list(c("mean_scale", "variance_scale"), NULL)
  )
  for (b in seq_len(n_boot)) {
    y_b <- .simulate_single_cell_bulk(
      fit$composition,
      mean_signature_matrix,
      Sigma,
      fit$mean_scale,
      fit$variance_scale
    )
    fit_b <- tryCatch(
      deconvolute_single_cell(
        y_b,
        mean_signature_matrix,
        Sigma,
        epsilon = epsilon,
        itmax = itmax,
        profile = TRUE,
        optimiser = optimiser
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      next
    }
    estimates[, b] <- fit_b$composition
    scales[1L, b] <- fit_b$mean_scale
    scales[2L, b] <- fit_b$variance_scale
  }
  list(
    estimates = estimates,
    interval = .percentile_interval(estimates, level),
    n_failed = sum(is.na(estimates[1L, ])),
    scales = .derived_scales(scales),
    scale_interval = .percentile_interval(.derived_scales(scales), level),
    fit = fit
  )
}

#' Reference bootstrap for the single-cell mixture
#'
#' @description
#' Resamples the reference that produced
#' \eqn{(\tilde{\boldsymbol{\mu}},\boldsymbol{\Sigma}_j)} and refits
#' [deconvolute_single_cell()]. `resample = NULL` uses donors when every
#' cell type has at least two donor labels, and otherwise resamples
#' columns within each cell type. The within-type resample regenerates
#' the bulk from the fitted single-cell law, at the original
#' \eqn{(\hat{\boldsymbol{p}},\hat c,\hat\phi)} and the original
#' moments, and refits it with the resampled moments, so that both the
#' bulk noise and the reference error enter the interval.
#'
#' Shuffling cell-type names is not a reference bootstrap. It does not
#' redraw cells or donors, and the likelihood is equivariant under a
#' relabeling, so the draw would not be a source of uncertainty. That
#' check remains [equivariance_check_decovart()] for the population model.
#'
#' @inheritParams bootstrap_single_cell
#' @param reference_profiles Named list of gene-by-cell (or
#'   gene-by-replicate) matrices, one per cell type.
#' @param donor_ids Optional list of donor labels, aligned with
#'   `reference_profiles`.
#' @param resample `"donors"`, `"cells"`, or `NULL` for the automatic
#'   rule above.
#' @param regenerate_bulk `NULL` regenerates the bulk only for the
#'   within-type resample. `TRUE` or `FALSE` overrides that rule.
#' @param ridge Ridge passed to the sample covariance, relative to its
#'   mean diagonal.
#'
#' @return A list with `resample`, `regenerate_bulk`, `estimates`,
#'   `interval`, `n_failed`, `scales`, `scale_interval` and `fit`;
#'   `scales` and `scale_interval` also carry the derived rows
#'   `platform_scale_a` and `cell_count_n`.
#'
#' @seealso [reference_bootstrap_decovart()], [confint_single_cell()]
#' @family decovart_inference
#' @export
#' @examples
#' set.seed(1)
#' mu <- cbind(ct1 = c(2, 2, 0.2, 0.2), ct2 = c(0.2, 0.2, 3, 3))
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' refs <- lapply(seq_len(2), function(j) {
#'   t(MASS::mvrnorm(n = 6, mu = mu[, j], Sigma = sigma[, , j]))
#' })
#' names(refs) <- colnames(mu)
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + rnorm(4, sd = 0.2)
#' reference_bootstrap_single_cell(y, refs, n_boot = 4)$resample
reference_bootstrap_single_cell <- function(
  y,
  reference_profiles,
  donor_ids = NULL,
  resample = NULL,
  n_boot = 100L,
  level = 0.95,
  regenerate_bulk = NULL,
  ridge = 1e-6,
  epsilon = 1e-6,
  itmax = 40L
) {
  level <- .check_level(level)
  plan <- .reference_resample_plan(donor_ids)
  if (is.null(resample)) {
    resample <- plan$method
  } else {
    resample <- match.arg(resample, c("donors", "cells"))
  }
  if (is.null(regenerate_bulk)) {
    regenerate_bulk <- !identical(resample, "donors")
  }
  checked <- .check_reference_profiles(
    reference_profiles,
    donor_ids,
    resample
  )
  refs <- checked$refs
  donor_ids <- checked$donor_ids
  if (nrow(refs[[1L]]) <= length(refs)) {
    stop(
      "The single-cell model needs more genes than cell types.",
      call. = FALSE
    )
  }
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  moments <- .moments_from_profiles(refs, ridge = ridge)
  fit <- deconvolute_single_cell(
    y,
    moments$mean_signature_matrix,
    moments$Sigma,
    epsilon = epsilon,
    itmax = itmax,
    profile = TRUE
  )
  estimates <- matrix(
    NA_real_,
    nrow = length(refs),
    ncol = n_boot,
    dimnames = list(names(refs), NULL)
  )
  scales <- matrix(
    NA_real_,
    nrow = 2L,
    ncol = n_boot,
    dimnames = list(c("mean_scale", "variance_scale"), NULL)
  )
  for (b in seq_len(n_boot)) {
    refs_b <- .resample_reference_profiles(refs, resample, donor_ids)
    moments_b <- .moments_from_profiles(refs_b, ridge = ridge)
    # Drawn from the fitted model (original moments), fitted with the
    # resampled ones, so the reference error reaches the refit.
    y_b <- if (isTRUE(regenerate_bulk)) {
      .simulate_single_cell_bulk(
        fit$composition,
        moments$mean_signature_matrix,
        moments$Sigma,
        fit$mean_scale,
        fit$variance_scale
      )
    } else {
      y
    }
    fit_b <- tryCatch(
      deconvolute_single_cell(
        y_b,
        moments_b$mean_signature_matrix,
        moments_b$Sigma,
        epsilon = epsilon,
        itmax = itmax,
        profile = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      next
    }
    estimates[, b] <- fit_b$composition
    scales[1L, b] <- fit_b$mean_scale
    scales[2L, b] <- fit_b$variance_scale
  }
  list(
    resample = resample,
    regenerate_bulk = isTRUE(regenerate_bulk),
    estimates = estimates,
    interval = .percentile_interval(estimates, level),
    n_failed = sum(is.na(estimates[1L, ])),
    scales = .derived_scales(scales),
    scale_interval = .percentile_interval(.derived_scales(scales), level),
    fit = fit
  )
}

# Dispatcher ---------------------------------------------------------------

#' Confidence intervals for one single-cell bulk sample
#'
#' @description
#' The composition and the two technical scales are separated by the
#' Schur complement of the observed Hessian. `method` chooses the
#' interval for the composition. The default is the inverse ILR chart
#' applied to a Wald interval for \eqn{\boldsymbol{\rho}}, using the
#' profile information
#' \eqn{\{-\nabla^{2}_{\boldsymbol{\rho}}\ell^{\ast}\}^{-1}}.
#'
#' Intervals for \eqn{(c,\phi)} in `nuisance$interval` are Wald intervals
#' on the log scale, mapped back by the exponential, from the Schur
#' complement that removes \eqn{\boldsymbol{\rho}}. The rows
#' `platform_scale_a` (\eqn{a=\phi/c}) and `cell_count_n`
#' (\eqn{n=c^{2}/\phi}, a number of cells only when \eqn{a} is known)
#' follow by the delta method on \eqn{\log a=\log\phi-\log c} and
#' \eqn{\log n=2\log c-\log\phi}, so all four intervals are positive.
#' `nuisance$interval_wald` is the symmetric normal interval of
#' \eqn{(c,\phi)} on the original scale. Bootstrap methods replace both
#' with percentile intervals.
#'
#' \describe{
#'   \item{`ilr`}{Default. Inverse ILR chart.}
#'   \item{`wald`}{Delta method on the composition. Endpoints can leave
#'     \eqn{[0,1]}.}
#'   \item{`firth`}{Wald interval at the penalised profile estimate.
#'     `interval_ilr` maps that estimate through the inverse chart.}
#'   \item{`profile`}{Interior likelihood-ratio interval. A face null is
#'     [lrt_single_cell()].}
#'   \item{`bootstrap`}{Parametric percentile interval.}
#'   \item{`reference`}{Donor bootstrap, or a within-type cell resample
#'     that regenerates the bulk. See
#'     [reference_bootstrap_single_cell()].}
#' }
#'
#' These intervals treat the reference covariance as known. They do not
#' add a sandwich or higher-order correction for that plug-in.
#'
#' `variance = "reml"` corrects for the \eqn{\nu=G-J} residual degrees of
#' freedom of the variance scale: covariances are inflated by
#' \eqn{G/\nu}, Wald and ILR quantiles come from \eqn{t_{\nu}} (a
#' multivariate \eqn{t} for the ILR map when \eqn{J>2}), \eqn{\phi} gets
#' the chi-square interval \eqn{[G\hat\phi/\chi^2_{\nu;1-\alpha/2},
#' G\hat\phi/\chi^2_{\nu;\alpha/2}]}, and the profile interval uses the
#' cut-off of [confint_profile_single_cell()]. The maximum-likelihood
#' \eqn{\hat\phi} is biased by about \eqn{\nu/G}, so with few genes
#' (\eqn{G=10}, \eqn{J=2}) the `"ml"` intervals cover about 0.84 to 0.93
#' at the 0.95 level; the bootstrap methods are unaffected by `variance`.
#'
#' @inheritParams deconvolute_single_cell
#' @param method Interval for the composition. See Description.
#' @param level Confidence level.
#' @param variance `"ml"` (default) or `"reml"`; see Description.
#' @param n_draws Monte Carlo size for the ILR map when \eqn{J>2}.
#' @param n_boot Bootstrap size.
#' @param seed Optional seed, restored on exit, for the ILR Monte Carlo
#'   map and for the bootstrap methods.
#' @param reference_profiles,donor_ids,resample,regenerate_bulk Passed
#'   to [reference_bootstrap_single_cell()] when
#'   `method = "reference"`.
#' @param ... `epsilon` and `itmax` for the selected method.
#'
#' @return A list with `method`, `estimate` (the fitted composition),
#'   `interval` and `nuisance`. Bootstrap results also include the
#'   replicate matrices.
#'
#' @seealso [vcov_single_cell()], [confint_decovart()], [lrt_single_cell()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' colnames(mu) <- c("ct1", "ct2")
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' fit <- confint_single_cell(y, mu, sigma)
#' fit$interval
#' fit$nuisance$interval
confint_single_cell <- function(
  y,
  mean_signature_matrix = NULL,
  Sigma = NULL,
  method = c(
    "ilr",
    "wald",
    "firth",
    "profile",
    "bootstrap",
    "reference"
  ),
  level = 0.95,
  variance = c("ml", "reml"),
  n_draws = 2000L,
  n_boot = 200L,
  seed = NULL,
  reference_profiles = NULL,
  donor_ids = NULL,
  resample = NULL,
  regenerate_bulk = NULL,
  ...
) {
  method <- match.arg(method)
  variance <- match.arg(variance)
  level <- .check_level(level)
  if (identical(method, "reference")) {
    if (is.null(reference_profiles)) {
      stop(
        "method = \"reference\" needs `reference_profiles`.",
        call. = FALSE
      )
    }
    ref <- .with_seed(
      seed,
      reference_bootstrap_single_cell(
        y = y,
        reference_profiles = reference_profiles,
        donor_ids = donor_ids,
        resample = resample,
        n_boot = n_boot,
        level = level,
        regenerate_bulk = regenerate_bulk,
        ...
      )
    )
    return(list(
      method = "reference",
      resample = ref$resample,
      regenerate_bulk = ref$regenerate_bulk,
      estimate = ref$fit$composition,
      interval = ref$interval,
      nuisance = list(
        estimate = c(
          mean_scale = ref$fit$mean_scale,
          variance_scale = ref$fit$variance_scale
        ),
        interval = ref$scale_interval
      ),
      reference = ref
    ))
  }
  if (is.null(mean_signature_matrix) || is.null(Sigma)) {
    stop(
      "`mean_signature_matrix` and `Sigma` are required unless ",
      "method = \"reference\".",
      call. = FALSE
    )
  }
  if (identical(method, "bootstrap")) {
    boot <- .with_seed(
      seed,
      bootstrap_single_cell(
        y,
        mean_signature_matrix,
        Sigma,
        n_boot = n_boot,
        level = level,
        ...
      )
    )
    return(list(
      method = "bootstrap",
      estimate = boot$fit$composition,
      interval = boot$interval,
      nuisance = list(
        estimate = c(
          mean_scale = boot$fit$mean_scale,
          variance_scale = boot$fit$variance_scale
        ),
        interval = boot$scale_interval
      ),
      bootstrap = boot
    ))
  }
  calibration <- .single_cell_calibration(
    variance,
    n_genes = nrow(as.matrix(mean_signature_matrix)),
    n_types = ncol(as.matrix(mean_signature_matrix))
  )
  if (identical(method, "profile")) {
    profile <- confint_profile_single_cell(
      y,
      mean_signature_matrix,
      Sigma,
      level = level,
      variance = variance,
      ...
    )
    estimate <- stats::setNames(profile[, "estimate"], rownames(profile))
    interval <- profile[, c("lower", "upper"), drop = FALSE]
    colnames(interval) <- .level_bounds(level)
    pieces <- .composition_intervals(
      estimate,
      y,
      mean_signature_matrix,
      Sigma,
      level,
      how = "wald",
      n_draws = n_draws,
      seed = seed,
      calibration = calibration
    )
    return(list(
      method = "profile",
      estimate = estimate,
      interval = interval,
      variance = variance,
      nuisance = pieces$nuisance
    ))
  }
  if (identical(method, "firth")) {
    fit <- fit_firth_single_cell(
      y,
      mean_signature_matrix,
      Sigma,
      ...
    )
    pieces <- .composition_intervals(
      fit$coefficients,
      y,
      mean_signature_matrix,
      Sigma,
      level,
      how = "wald",
      n_draws = n_draws,
      seed = seed,
      calibration = calibration
    )
    mapped <- .composition_intervals(
      fit$coefficients,
      y,
      mean_signature_matrix,
      Sigma,
      level,
      how = "ilr",
      n_draws = n_draws,
      seed = seed,
      calibration = calibration
    )
    return(list(
      method = "firth",
      estimate = fit$coefficients,
      interval = pieces$interval,
      interval_ilr = mapped$interval,
      nuisance = pieces$nuisance,
      fit = fit
    ))
  }
  fit <- deconvolute_single_cell(
    y,
    mean_signature_matrix,
    Sigma,
    profile = TRUE,
    ...
  )
  pieces <- .composition_intervals(
    fit$composition,
    y,
    mean_signature_matrix,
    Sigma,
    level,
    how = method,
    n_draws = n_draws,
    seed = seed,
    calibration = calibration
  )
  c(list(method = method), pieces)
}
