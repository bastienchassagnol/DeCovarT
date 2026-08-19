#' Validate generative-model parameters \eqn{\theta}
#'
#' Checks that `true_theta` stores the DeCovarT mixture parameters:
#' * `p`: mixing proportions of length \eqn{J}, or a matrix
#'   \eqn{J\times N} of sample-wise proportions (each column on the simplex);
#' * `mu`: mean signature \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}
#'   (also accepts \eqn{J\times G});
#' * `sigma` and/or `Theta`: array
#'   \eqn{G\times G\times J} of covariances and/or precision matrices.
#'
#' @param true_theta Named list with the fields above.
#' @param require_p Logical; if `TRUE` (default), `p` must be present and
#'   valid. If `FALSE`, a missing `p` is allowed (callers may default it).
#' @param J Optional expected number of cell types. When `NULL`, inferred
#'   from `mu` / the third dimension of `sigma` or `Theta`.
#' @param second_moment Which second-moment field to require:
#'   `"sigma"`, `"Theta"`, or `"either"` (default: at least one of
#'   `sigma` / `Theta`). Matching is case-insensitive.
#'
#' @srrstats {G2.3} Restricted character input (`second_moment`).
#' @srrstats {G2.3a} Validated via `.match_arg_ci()` (a `match.arg()`
#'   equivalent).
#' @srrstats {G2.3b} Matching is case-insensitive (`tolower()`).
#'
#' @return `TRUE` invisibly if the structure is valid; otherwise stops
#'   with an informative error.
#'
#' @seealso [compute_average_overlap()], [compute_average_jeffreys()]
#' @export
#' @examples
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' check_true_theta(theta)
check_true_theta <- function(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
) {
  .parse_true_theta(
    true_theta = true_theta,
    require_p = require_p,
    J = J,
    second_moment = second_moment
  )
  invisible(TRUE)
}

#' Parse and standardise \eqn{\theta} (internal)
#'
#' @return Named list with `p` (length-\eqn{J} vector or `NULL`),
#'   `mu` (\eqn{G\times J}), `sigma` and/or `Theta` (\eqn{G\times G\times J}),
#'   `G`, and `J`.
#' @keywords internal
#' @noRd
.parse_true_theta <- function(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
) {
  second_moment <- .match_arg_ci(second_moment, c("either", "sigma", "Theta"))
  if (!is.list(true_theta)) {
    stop("`true_theta` must be a list.", call. = FALSE)
  }
  if (is.null(true_theta$mu)) {
    stop("`true_theta` must contain `mu`.", call. = FALSE)
  }

  # --- Second moments: require Sigma and/or precision Theta -------------------
  # Each must be a G x G x J array (one PD matrix per cell type). When both are
  # supplied they must share the same (G, J).
  has_sigma <- !is.null(true_theta$sigma)
  has_theta <- !is.null(true_theta$Theta)
  if (identical(second_moment, "sigma") && !has_sigma) {
    stop("`true_theta` must contain `sigma`.", call. = FALSE)
  }
  if (identical(second_moment, "Theta") && !has_theta) {
    stop("`true_theta` must contain `Theta`.", call. = FALSE)
  }
  if (identical(second_moment, "either") && !has_sigma && !has_theta) {
    stop(
      "`true_theta` must contain `sigma` and/or `Theta`.",
      call. = FALSE
    )
  }

  # Helper: enforce cubic layout G x G x J and return (G, J).
  .check_ggj_array <- function(arr, nm) {
    if (is.null(dim(arr)) || length(dim(arr)) != 3L) {
      stop(
        "`true_theta$",
        nm,
        "` must be a G x G x J array.",
        call. = FALSE
      )
    }
    g1 <- dim(arr)[[1L]]
    g2 <- dim(arr)[[2L]]
    jj <- dim(arr)[[3L]]
    if (g1 != g2) {
      stop(
        "`true_theta$",
        nm,
        "` dims must be G x G x J.",
        call. = FALSE
      )
    }
    list(G = g1, J = jj)
  }

  dims_list <- list()
  if (has_sigma) {
    dims_list$sigma <- .check_ggj_array(true_theta$sigma, "sigma")
  }
  if (has_theta) {
    dims_list$Theta <- .check_ggj_array(true_theta$Theta, "Theta")
  }
  if (length(dims_list) == 2L) {
    if (
      !identical(dims_list$sigma$G, dims_list$Theta$G) ||
        !identical(dims_list$sigma$J, dims_list$Theta$J)
    ) {
      stop("`sigma` and `Theta` must share the same G x G x J dims.")
    }
  }
  G <- dims_list[[1L]]$G
  n_celltypes <- dims_list[[1L]]$J

  # --- Mean signature mu: genes x cell types ---------------------------------
  # Preferred storage is G x J; also accept MixSim-style J x G and transpose.
  mu <- as.matrix(true_theta$mu)
  if (!is.numeric(mu) || anyNA(mu)) {
    stop("`true_theta$mu` must be numeric without missing values.")
  }
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_gj <- mu
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_gj <- t(mu)
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  # --- Mixing proportions p: unit simplex ------------------------------------
  # Accept a single mixture weight vector (length J) or sample-wise matrices
  # (J x N or N x J). Matrix forms are reduced to length-J by averaging so that
  # MixSim / Jeffreys summaries still see one weight per cell type. Every
  # supplied simplex (vector, or each column/row of a matrix) must be
  # non-negative and sum to 1.
  p_vec <- NULL
  if (is.null(true_theta$p)) {
    if (isTRUE(require_p)) {
      stop("`true_theta` must contain `p`.", call. = FALSE)
    }
  } else {
    p_raw <- true_theta$p
    if (is.null(dim(p_raw))) {
      # Length-J vector on the simplex (shared mixture weights).
      p_vec <- as.numeric(p_raw)
      if (length(p_vec) != n_celltypes) {
        stop(
          "`true_theta$p` must have length J (= third dim of sigma/Theta).",
          call. = FALSE
        )
      }
    } else {
      p_mat <- as.matrix(p_raw)
      if (nrow(p_mat) == n_celltypes) {
        # J x N: each column is a sample-wise cellular ratio on the simplex;
        # average across samples for mixture-level summaries.
        if (
          any(abs(colSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as J x N must have non-negative columns ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- rowMeans(p_mat)
      } else if (ncol(p_mat) == n_celltypes) {
        # N x J: each row is a sample-wise ratio; average across samples.
        if (
          any(abs(rowSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as N x J must have non-negative rows ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- colMeans(p_mat)
      } else {
        stop(
          "`true_theta$p` must be length J, J x N, or N x J.",
          call. = FALSE
        )
      }
    }
    # Final simplex check on the (possibly averaged) length-J weight vector.
    if (
      any(p_vec < 0) ||
        anyNA(p_vec) ||
        abs(sum(p_vec) - 1) > 100 * .Machine$double.eps
    ) {
      stop(
        "`true_theta$p` must be non-negative and sum to 1 ",
        "(after averaging over samples if matrix).",
        call. = FALSE
      )
    }
  }

  # --- Cell-type count J -----------------------------------------------------
  # Infer from the third dim of sigma/Theta, or verify a caller-supplied J.
  if (is.null(J)) {
    J <- n_celltypes
  } else {
    J <- as.integer(J)
    if (length(J) != 1L || is.na(J) || J != n_celltypes) {
      stop(
        "`J` must match the third dimension of sigma/Theta.",
        call. = FALSE
      )
    }
  }
  if (J < 2L) {
    stop("`J` must be at least 2.", call. = FALSE)
  }

  # Standardised fields for downstream metrics (mu always G x J).
  list(
    p = p_vec,
    mu = mu_gj,
    sigma = if (has_sigma) true_theta$sigma else NULL,
    Theta = if (has_theta) true_theta$Theta else NULL,
    G = G,
    J = J
  )
}

#' Compute summary metrics for estimated proportions
#'
#' @description
#' When a ground truth \eqn{\boldsymbol{p}^{\star}} is supplied, scores compare
#' \eqn{\hat{\boldsymbol{p}}} to \eqn{\boldsymbol{p}^{\star}}. Otherwise scores
#' compare the reconstituted bulk
#' \eqn{\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}} to the
#' observed \eqn{\boldsymbol{y}}.
#'
#' @inheritParams deconvolute_ratios_Marquardt_Levenberg
#' @param estimated_p Estimated proportions
#'   \eqn{\hat{\boldsymbol{p}}\in\mathbb{R}^{J}}.
#' @param true_ratios Optional ground-truth proportions
#'   \eqn{\boldsymbol{p}^{\star}\in\mathbb{R}^{J}}. When supplied, metrics
#'   compare \eqn{\hat{\boldsymbol{p}}} to \eqn{\boldsymbol{p}^{\star}};
#'   otherwise they compare
#'   \eqn{\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}} to
#'   \eqn{\boldsymbol{y}}.
#'
#' @return A `tibble` with mse/rmse/mae, optionally
#'   \eqn{R^{2}} / adjusted \eqn{R^{2}}, and Pearson correlation.
#'
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2,
#'              dimnames = list(paste0("g", 1:2), paste0("ct", 1:2)))
#' y <- drop(mu %*% c(0.4, 0.6))
#' compute_benchmark_metrics(y, mu, estimated_p = c(0.45, 0.55),
#'                           true_ratios = c(0.4, 0.6))
#' @importFrom rlang .data
#' @export
compute_benchmark_metrics <- function(
  y,
  mean_signature_matrix,
  estimated_p,
  true_ratios = NULL
) {
  G <- nrow(mean_signature_matrix)
  J <- ncol(mean_signature_matrix)
  # free parameters: J - 1 proportions (simplex); residual df uses G genes
  df_res <- G - J + 1
  df_tot <- G - 1
  if (!is.null(true_ratios)) {
    # when the true parameters are known
    model_coef_determination <- max(
      0,
      1 - Metrics::rse(true_ratios, estimated_p)
    )
    scores <- tibble::tibble(
      model_mse = Metrics::mse(true_ratios, estimated_p),
      model_rmse = Metrics::rmse(true_ratios, estimated_p),
      model_mae = Metrics::mae(true_ratios, estimated_p),
      model_coef_determination = model_coef_determination,
      model_coef_determination_adjusted = max(
        0,
        1 - (1 - model_coef_determination) * df_tot / df_res
      ),
      model_cor = suppressWarnings(stats::cor(
        true_ratios,
        estimated_p,
        method = "pearson"
      ))
    )
  } else {
    # when they are unknown
    predicted_values <- as.vector(mean_signature_matrix %*% estimated_p)
    scores <- tibble::tibble(
      model_mse = Metrics::mse(y, predicted_values),
      model_rmse = Metrics::rmse(y, predicted_values),
      model_mae = Metrics::mae(y, predicted_values),
      model_cor = suppressWarnings(stats::cor(
        y,
        predicted_values,
        method = "pearson"
      ))
    )
  }
  return(scores)
}

#' Parallel deconvolution of a bulk expression matrix
#'
#' @description
#' For each column \eqn{\boldsymbol{y}_{\cdot i}} of the bulk matrix
#' \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}}, estimates
#' \eqn{\hat{\boldsymbol{p}}_{\cdot i}} with every supplied algorithm. When
#' covariance information is provided, DeCovarT methods maximise
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} under
#' \eqn{\boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(
#'   \boldsymbol{\mu}\boldsymbol{p},
#'   \boldsymbol{\Sigma}(\boldsymbol{p})
#' )}.
#'
#' @param signature_matrix Mean signature
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (rownames = genes,
#'   colnames = cell types). Numeric matrix of non-negative expression;
#'   \eqn{G \ge J} is required (the mixture is undetermined if
#'   \eqn{J > G}). Used as a frequentist **plug-in** for unobserved
#'   latent cell-type profiles \eqn{\boldsymbol{x}_{\cdot j}}.
#' @param bulk_expression Bulk matrix
#'   \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}}: numeric, non-negative,
#'   gene rownames matching `signature_matrix`.
#' @param true_ratios Optional ground-truth proportions: a length-\eqn{J}
#'   numeric vector (recycled over samples) or a numeric matrix
#'   \eqn{J\times N} (or \eqn{N\times J}).
#' @param Sigma Optional array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}
#'   of cell-type covariances (numeric; off-diagonals may be negative).
#' @param deconvolution_functions Named list; each element has `FUN` and
#'   optional `additional_parameters`.
#' @param scaled If `TRUE`, apply a log2 transform before estimation.
#' @param cores Number of parallel workers.
#'
#' @return A `tibble` of estimated \eqn{\hat{\boldsymbol{p}}} and metrics,
#'   after [repair_simplex()].
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
#' Y <- simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 2)$Y
#' deconvolute_ratios(
#'   signature_matrix = mu,
#'   bulk_expression = Y,
#'   true_ratios = c(0.5, 0.5),
#'   Sigma = Sigma,
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   cores = 1
#' )
#' @srrstats {G1.0} Closest published statistical method is DSection
#'   (Erkkila et al. 2010): univariate Gaussian convolution. DeMix, DeMixT
#'   and ISOpureR also treat cell-type profiles as latent; Ogundijo and Wang
#'   (2017) extend DSection by sequential Monte Carlo. DeCovarT generalises
#'   DSection to a multivariate convolution with sparse cell-type covariance.
#' @srrstats {G1.1} First implementation (in R or otherwise) of that
#'   multivariate Gaussian-convolution MLE for bulk deconvolution.
#' @srrstats {G2.8} Unique pre-processing gate:
#'   `.prepare_deconvolution_inputs()` then passes numeric matrices /
#'   arrays to every solver.
#' @srrstats {G2.13} Missing values in `y`, the signature, `Sigma` or
#'   `true_ratios` raise an error here, before any solver is called.
#' @srrstats {G2.14a} The missing-data policy is error-only (no imputation).
#' @srrstats {G2.15} After this check, `mean` / `cor` / `var` never receive
#'   incomplete expression data (default `na.rm = FALSE` is then safe).
#' @srrstats {G2.16} `NaN` / `Inf` / `-Inf` are rejected with the same error.
#' @srrstats {G5.8d} \eqn{J > G} raises an undetermined-mixture error.
#' @references
#' \insertRef{erkkilaProbabilisticAnalysisGene2010}{DeCovarT}
#' \insertRef{ahnDeMixDeconvolutionMixed2013}{DeCovarT}
#' \insertRef{wangTranscriptomeDeconvolutionHeterogeneous2018}{DeCovarT}
#' \insertRef{anghelISOpureRImplementationComputational2015}{DeCovarT}
#' \insertRef{ogundijoSequentialMonteCarlo2017}{DeCovarT}
#' \insertRef{hafemeisterNormalizationVarianceStabilization2019}{DeCovarT}
#' \insertRef{chionBayesianFrameworkMultivariate2023}{DeCovarT}
#' @importFrom rlang .data
#' @export
#' @seealso [deconvolute_ratios_Marquardt_Levenberg()]
deconvolute_ratios <- function(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  deconvolution_functions = NULL,
  scaled = FALSE,
  cores = ifelse(
    .Platform$OS.type == "unix",
    getOption("mc.cores", parallel::detectCores()),
    1
  )
) {
  aligned <- .prepare_deconvolution_inputs(
    signature_matrix = signature_matrix,
    bulk_expression = bulk_expression,
    true_ratios = true_ratios,
    Sigma = Sigma,
    scaled = scaled
  )
  mean_signature_matrix <- aligned$signature_matrix
  Y <- aligned$bulk_expression
  true_ratios <- aligned$true_ratios
  Sigma <- aligned$Sigma

  # estimation itself
  deconvolution_estimates <- purrr::imap_dfr(
    deconvolution_functions,
    function(deconvolution_function, algorithm) {
      additional_parameters <- deconvolution_function$additional_parameters
      metric_scores <- parallel::mclapply(
        seq_len(ncol(Y)),
        function(i) {
          y_i <- Y[, i, drop = TRUE]
          p_i <- if (is.null(true_ratios)) {
            NULL
          } else {
            true_ratios[, i, drop = TRUE]
          }
          list_arguments <- c(
            list(
              "y" = y_i,
              "mean_signature_matrix" = mean_signature_matrix,
              "Sigma" = Sigma
            ),
            additional_parameters
          )
          formal_args <- methods::formalArgs(deconvolution_function$FUN)
          success_estimation <- tryCatch(
            {
              estimated_p <- do.call(
                deconvolution_function$FUN,
                list_arguments[names(list_arguments) %in% formal_args]
              )
              compute_benchmark_metrics(
                y = y_i,
                mean_signature_matrix = mean_signature_matrix,
                estimated_p = estimated_p,
                true_ratios = p_i
              ) |>
                dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
            },
            error = function(e) {
              warning(conditionMessage(e), call. = FALSE)
              return(e)
            }
          )
          if (!inherits(success_estimation, "error")) {
            return(success_estimation)
          }
        },
        mc.cores = cores
      ) |>
        # One estimation with a given deconvolution algorithm terminated.
        dplyr::bind_rows()

      if (nrow(metric_scores) != 0) {
        metric_scores <- metric_scores |>
          dplyr::mutate(
            OMIC_ID = paste0("sample_", seq_len(nrow(metric_scores))),
            algorithm = algorithm
          )
      }
      return(metric_scores)
    }
  )
  return(deconvolution_estimates)
}
