# Resampling-based inference for DeCovarT: the parametric bootstrap (percentile
# intervals, and Monte Carlo calibration of boundary likelihood-ratio tests),
# the reference-sample bootstrap that propagates uncertainty in the plug-in
# moments, and the label-equivariance check of the estimator.
# See vignette("theory-population-inference").

# Parametric bootstrap -----------------------------------------------------

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
#' @param null_value Optional numeric vector of hypothesised ratios named
#'   by cell types; switches to the restricted bootstrap (see Details).
#' @param p Composition to simulate from. Defaults to the MLE computed
#'   from `y` (or the restricted MLE under `null_value`).
#' @param n_boot Number of bootstrap replicates.
#' @param level Confidence level for percentile intervals.
#' @param n_replicates Number of bulk columns simulated per replicate;
#'   defaults to the number of observed columns.
#' @param optimiser Solver for the unrestricted fits: `"newton"`
#'   (trust-region Newton in ILR coordinates) or `"lbfgsb"` (L-BFGS-B on
#'   the closed simplex, as in [deconvolute_ratios_L_BFGS_B()]), so that
#'   the bootstrap distribution matches the estimator being assessed.
#'   Restricted fits under `null_value` always use the Newton driver.
#'
#' @return A list with `estimates` (\eqn{J\times} `n_boot` matrix),
#'   `interval` (percentile bounds), `n_failed` (refits that errored,
#'   left as `NA` columns), `p_simulated`, and, when
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
  y,
  mean_signature_matrix,
  Sigma,
  p = NULL,
  null_value = NULL,
  n_boot = 500L,
  level = 0.95,
  n_replicates = NULL,
  epsilon = 10^-8,
  itmax = 500L,
  optimiser = c("newton", "lbfgsb")
) {
  optimiser <- match.arg(optimiser)
  y <- as.matrix(y)
  if (is.null(n_replicates)) {
    n_replicates <- ncol(y)
  }
  n_celltypes <- ncol(mean_signature_matrix)
  fit_bulk <- function(bulk, fixed = NULL) {
    if (is.null(fixed) && identical(optimiser, "lbfgsb")) {
      fit <- .lbfgsb_closure(
        rep(1 / n_celltypes, n_celltypes),
        loglik = function(p) {
          .pooled_population_loglik(p, bulk, mean_signature_matrix, Sigma)
        },
        score = function(p) {
          .pooled_population_score(p, bulk, mean_signature_matrix, Sigma)
        },
        epsilon = epsilon,
        itmax = itmax
      )
      return(list(
        coefficients = stats::setNames(
          fit$par,
          colnames(mean_signature_matrix)
        ),
        loglik = fit$loglik
      ))
    }
    restricted_mle_decovart(
      bulk,
      mean_signature_matrix,
      Sigma,
      fixed = fixed,
      epsilon = epsilon,
      itmax = itmax
    )
  }
  observed <- fit_bulk(y)
  statistic <- NULL
  p_simulated <- observed$coefficients
  if (!is.null(null_value)) {
    restricted <- fit_bulk(y, null_value)
    p_simulated <- restricted$coefficients
    statistic <- max(2 * (observed$loglik - restricted$loglik), 0)
  }
  if (!is.null(p)) {
    p_simulated <- p
  }
  estimates <- matrix(
    NA_real_,
    nrow = length(p_simulated),
    ncol = n_boot,
    dimnames = list(colnames(mean_signature_matrix), NULL)
  )
  null_statistics <- rep(NA_real_, n_boot)
  for (b in seq_len(n_boot)) {
    y_boot <- .simulate_bulk_replicates(
      p_simulated,
      mean_signature_matrix,
      Sigma,
      n_replicates
    )
    boot_full <- tryCatch(fit_bulk(y_boot), error = function(e) NULL)
    if (is.null(boot_full)) {
      next
    }
    estimates[, b] <- boot_full$coefficients
    if (!is.null(statistic)) {
      boot_null <- fit_bulk(y_boot, null_value)
      null_statistics[[b]] <- 2 *
        (max(boot_full$loglik, boot_null$loglik) - boot_null$loglik)
    }
  }
  out <- list(
    estimates = estimates,
    interval = .percentile_interval(estimates, level),
    n_failed = sum(is.na(estimates[1L, ])),
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

#' Check a named list of purified profiles and optional donor labels
#'
#' @noRd
.check_reference_profiles <- function(
  reference_profiles,
  donor_ids,
  method
) {
  method <- match.arg(method, c("donors", "cells"))
  refs <- lapply(reference_profiles, as.matrix)
  n_cols <- vapply(refs, ncol, integer(1))
  if (
    !is.list(reference_profiles) ||
      is.null(names(refs)) ||
      !all(nzchar(names(refs))) ||
      length(unique(vapply(refs, nrow, integer(1)))) != 1L ||
      any(n_cols < 2L)
  ) {
    stop(
      "`reference_profiles` must be a named list of gene x sample ",
      "matrices with the same genes and at least two columns per ",
      "cell type.",
      call. = FALSE
    )
  }
  if (identical(method, "donors")) {
    if (is.null(donor_ids)) {
      stop(
        "Supply `donor_ids` for method = \"donors\", ",
        "or use method = \"cells\".",
        call. = FALSE
      )
    }
    donor_ids <- lapply(donor_ids[names(refs)], as.character)
    if (!identical(lengths(donor_ids, use.names = FALSE), unname(n_cols))) {
      stop(
        "`donor_ids` must be named like `reference_profiles`, with one ",
        "label per purified column.",
        call. = FALSE
      )
    }
    if (
      any(vapply(donor_ids, function(d) length(unique(d)), integer(1)) < 2L)
    ) {
      stop(
        "method = \"donors\" needs at least two donors in every cell type.",
        call. = FALSE
      )
    }
  }
  list(refs = refs, donor_ids = donor_ids, method = method)
}

#' Draw donors or cells with replacement inside each labelled type
#'
#' @noRd
.resample_reference_profiles <- function(refs, method, donor_ids = NULL) {
  if (identical(method, "donors")) {
    return(Map(
      function(x, d) {
        donors <- unique(d)
        drawn <- sample(donors, size = length(donors), replace = TRUE)
        keep <- unlist(
          lapply(drawn, function(id) which(d == id)),
          use.names = FALSE
        )
        x[, keep, drop = FALSE]
      },
      refs,
      donor_ids
    ))
  }
  lapply(refs, function(x) {
    x[, sample.int(ncol(x), replace = TRUE), drop = FALSE]
  })
}

#' Donors when every type has two or more; otherwise cells, and regenerate
#'
#' Shuffling cell-type names is not a sampling model. One donor per type
#' falls back to resampling columns within the type.
#'
#' @noRd
.reference_resample_plan <- function(donor_ids) {
  if (is.null(donor_ids)) {
    return(list(method = "cells", regenerate_bulk = TRUE))
  }
  n_donors <- vapply(
    donor_ids,
    function(d) length(unique(as.character(d))),
    integer(1)
  )
  if (length(n_donors) == 0L || any(n_donors < 2L)) {
    return(list(method = "cells", regenerate_bulk = TRUE))
  }
  list(method = "donors", regenerate_bulk = FALSE)
}

# Reference-sample bootstrap -----------------------------------------------

#' Reference-based bootstrap for signature and composition uncertainty
#'
#' @description
#' The frequentist DeCovarT likelihood treats
#' \eqn{\boldsymbol{\mu}} and each \eqn{\boldsymbol{\Sigma}_j} as known.
#' In a labelled, reference-based model those moments come from purified
#' profiles whose *cell-type names stay attached*. This bootstrap
#' resamples the experimental units that generated the reference (donors,
#' by default), re-estimates the plug-in moments and refits
#' \insertCite{efronBootstrapMethodsAnother1979}{DeCovarT}.
#'
#' @details
#' Permuting gene order, or permuting cell-type labels of an already
#' averaged signature, is **not** a bootstrap: the first destroys the
#' gene-wise pairing of bulk and reference, and the second is algebraic
#' equivariance rather than a source of uncertainty (see
#' [equivariance_check_decovart()] and
#' `vignette("theory-population-inference")`). The two `method` options
#' instead resample, with replacement, units that actually vary.
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
#' }
#'
#' The observed bulk is kept fixed unless
#' `regenerate_bulk = TRUE`, in which case each replicate also draws
#' \eqn{\boldsymbol{Y}^{(b)}} from the fitted model, at `p` (the observed
#' MLE by default) and the moments of the full reference, and refits it
#' with the bootstrapped moments \eqn{S^{(b)}}. Drawing
#' \eqn{\boldsymbol{Y}^{(b)}} from \eqn{S^{(b)}} instead would make each
#' refit consistent with its own moments, so the reference error would not
#' reach \eqn{\hat{\boldsymbol{p}}^{(b)}}: with 60 cells per type, that
#' variant covered a known composition 76\% of the time at the 95\%
#' level, against 92\% for this one (60 replicates).
#'
#' When a cell type has fewer purified columns than genes, the sample
#' covariance is singular. A ridge `ridge * mean(diag(S))` is then added
#' on the diagonal so that \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} remains
#' factorisable. Set `ridge = 0` to disable that guard (and accept
#' possible Cholesky failures).
#'
#' @inheritParams bootstrap_decovart
#' @param reference_profiles Named list of \eqn{G\times n_j} matrices,
#'   one purified sample matrix per cell type (at least two columns each).
#'   The list order fixes the column order of the signature.
#' @param donor_ids Donor (or other clustering) labels: a list named like
#'   `reference_profiles`, with one label per purified column. Required for
#'   `method = "donors"`.
#' @param method One of `"donors"` (default) or `"cells"`.
#' @param regenerate_bulk If `TRUE`, simulate a new bulk from the fitted
#'   model and refit it with the bootstrapped moments, instead of keeping
#'   the observed \eqn{\boldsymbol{Y}} fixed.
#' @param ridge Non-negative ridge multiple for the sample covariance of
#'   each cell type. The default `1e-6` is relative to the mean diagonal.
#'
#' @return A list with `method`, `estimates` (\eqn{J\times} `n_boot`),
#'   `interval` (percentile bounds), `n_failed`, `mean_signature_matrix` /
#'   `Sigma`
#'   from the original reference, and `p_hat` from the original plug-in
#'   fit.
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
  y,
  reference_profiles,
  donor_ids = NULL,
  method = c("donors", "cells"),
  n_boot = 199L,
  level = 0.95,
  regenerate_bulk = FALSE,
  p = NULL,
  n_replicates = NULL,
  ridge = 1e-6,
  epsilon = 10^-8,
  itmax = 500L
) {
  checked <- .check_reference_profiles(
    reference_profiles,
    donor_ids,
    method
  )
  refs <- checked$refs
  donor_ids <- checked$donor_ids
  method <- checked$method
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 1L) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  y <- as.matrix(y)
  if (is.null(n_replicates)) {
    n_replicates <- ncol(y)
  }
  moments <- .moments_from_profiles(refs, ridge = ridge)
  observed <- restricted_mle_decovart(
    y,
    moments$mean_signature_matrix,
    moments$Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  p_simulated <- if (is.null(p)) observed$coefficients else p
  estimates <- matrix(
    NA_real_,
    nrow = length(refs),
    ncol = n_boot,
    dimnames = list(names(refs), NULL)
  )
  for (b in seq_len(n_boot)) {
    # Resample with replacement within each labelled cell type: whole donors
    # (cluster bootstrap) or individual purified columns.
    refs_b <- .resample_reference_profiles(refs, method, donor_ids)
    moments_b <- .moments_from_profiles(refs_b, ridge = ridge)
    # A regenerated bulk is drawn from the fitted model (original
    # moments) and fitted with the resampled moments; drawing it from the
    # resampled moments would make the refit consistent with them and
    # hide the reference error.
    y_b <- if (isTRUE(regenerate_bulk)) {
      .simulate_bulk_replicates(
        p_simulated,
        moments$mean_signature_matrix,
        moments$Sigma,
        n_replicates
      )
    } else {
      y
    }
    fit_b <- tryCatch(
      restricted_mle_decovart(
        y_b,
        moments_b$mean_signature_matrix,
        moments_b$Sigma,
        epsilon = epsilon,
        itmax = itmax
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) {
      estimates[, b] <- fit_b$coefficients
    }
  }
  list(
    method = method,
    estimates = estimates,
    interval = .percentile_interval(estimates, level),
    n_failed = sum(is.na(estimates[1L, ])),
    p_hat = observed$coefficients,
    mean_signature_matrix = moments$mean_signature_matrix,
    Sigma = moments$Sigma
  )
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

# Label-equivariance check -------------------------------------------------

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
#' `vignette("theory-population-inference")`. The labelled MLE on the
#' relabelled reference is then compared with
#' \eqn{\hat{\boldsymbol{p}}_{\mathrm{perm}}}.
#'
#' Use [reference_bootstrap_decovart()] for percentile intervals that
#' resample the experimental units of the reference (donors or cells).
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
  y,
  mean_signature_matrix,
  Sigma,
  perm = rev(seq_len(ncol(mean_signature_matrix))),
  epsilon = 10^-8,
  itmax = 500L
) {
  n_celltypes <- ncol(mean_signature_matrix)
  perm <- as.integer(perm)
  if (!identical(sort(perm), seq_len(n_celltypes))) {
    stop("`perm` must be a permutation of 1:J.", call. = FALSE)
  }
  # Reorder columns and covariance slices, keep the names in place: the
  # labelled composition of the relabelled reference is then p[perm].
  nms <- colnames(mean_signature_matrix)
  mu_perm <- mean_signature_matrix[, perm, drop = FALSE]
  colnames(mu_perm) <- nms
  sigma_perm <- Sigma[,, perm, drop = FALSE]
  if (!is.null(dimnames(Sigma)[[3L]])) {
    dimnames(sigma_perm)[[3L]] <- nms
  }
  observed <- restricted_mle_decovart(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  relabelled <- restricted_mle_decovart(
    y,
    mu_perm,
    sigma_perm,
    epsilon = epsilon,
    itmax = itmax
  )
  p_expected <- stats::setNames(observed$coefficients[perm], nms)
  list(
    perm = perm,
    p_hat = observed$coefficients,
    p_star = relabelled$coefficients,
    p_expected = p_expected,
    max_abs_diff = max(abs(relabelled$coefficients - p_expected)),
    loglik_diff = relabelled$loglik - observed$loglik
  )
}
