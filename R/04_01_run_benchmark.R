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
#'   `sigma` / `Theta`).
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
  second_moment <- match.arg(second_moment)
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
        paste0("`true_theta$", nm, "` must be a G x G x J array."),
        call. = FALSE
      )
    }
    g1 <- dim(arr)[[1L]]
    g2 <- dim(arr)[[2L]]
    jj <- dim(arr)[[3L]]
    if (g1 != g2) {
      stop(
        paste0("`true_theta$", nm, "` dims must be G x G x J."),
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
        if (any(abs(colSums(p_mat) - 1) > 1e-8) || any(p_mat < 0)) {
          stop(
            paste(
              "`true_theta$p` as J x N must have non-negative columns",
              "summing to 1."
            ),
            call. = FALSE
          )
        }
        p_vec <- rowMeans(p_mat)
      } else if (ncol(p_mat) == n_celltypes) {
        # N x J: each row is a sample-wise ratio; average across samples.
        if (any(abs(rowSums(p_mat) - 1) > 1e-8) || any(p_mat < 0)) {
          stop(
            paste(
              "`true_theta$p` as N x J must have non-negative rows",
              "summing to 1."
            ),
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
    if (any(p_vec < 0) || abs(sum(p_vec) - 1) > 1e-8) {
      stop(
        paste(
          "`true_theta$p` must be non-negative and sum to 1",
          "(after averaging over samples if matrix)."
        ),
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
#'   colnames = cell types). Used as a frequentist **plug-in** for unobserved
#'   latent cell-type profiles \eqn{\boldsymbol{x}_{\cdot j}}.
#' @param bulk_expression Bulk matrix
#'   \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}}.
#' @param true_ratios Optional ground-truth proportions for scoring.
#' @param Sigma Optional array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}.
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
  # read in data
  if (!is.matrix(signature_matrix) || is.null(row.names(signature_matrix))) {
    stop(paste(
      "required format for signature is expression matrix,",
      "with rownames as genes"
    ))
  }
  mean_signature_matrix <- tibble::as_tibble(
    signature_matrix,
    rownames = "GENE_SYMBOL"
  )
  if (!is.matrix(bulk_expression) || is.null(row.names(bulk_expression))) {
    stop(
      "required format for mixture is expression matrix, with rownames as genes"
    )
  }
  Y <- tibble::as_tibble(bulk_expression, rownames = "GENE_SYMBOL")

  # remove potential missing data from Y database
  Y <- Y |> tidyr::drop_na()

  # intersect genes (we only keep genes that are common to both data bases)
  common_genes <- intersect(mean_signature_matrix$GENE_SYMBOL, Y$GENE_SYMBOL)

  if (length(common_genes) / dim(mean_signature_matrix)[1] < 0.5) {
    stop(paste(
      "Only",
      length(common_genes) / dim(mean_signature_matrix)[1],
      "fraction of genes are used in the signature matrix\n.
                  Half of common genes are required at least"
    ))
  }
  mean_signature_matrix <- mean_signature_matrix |>
    dplyr::filter(.data$GENE_SYMBOL %in% common_genes) |>
    dplyr::arrange(.data$GENE_SYMBOL) |>
    dplyr::select(dplyr::where(is.numeric)) |>
    as.matrix()
  Y <- Y |>
    dplyr::filter(.data$GENE_SYMBOL %in% common_genes) |>
    dplyr::arrange(.data$GENE_SYMBOL) |>
    dplyr::select(dplyr::where(is.numeric))

  if (scaled) {
    # log-2 normalise
    Y <- log2(Y)
    mean_signature_matrix <- log2(mean_signature_matrix)
  }

  # estimation itself
  deconvolution_estimates <- purrr::imap_dfr(
    deconvolution_functions,
    function(deconvolution_function, algorithm) {
      additional_parameters <- deconvolution_function$additional_parameters
      metric_scores <- parallel::mclapply(
        seq_len(ncol(Y)),
        function(i) {
          # metric_scores <- tibble::tibble(); for (i in 1:ncol(Y)) {
          list_arguments <- c(
            list(
              "y" = as.numeric(Y[[i]]),
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
                y = as.numeric(Y[[i]]),
                mean_signature_matrix = mean_signature_matrix,
                estimated_p = estimated_p,
                true_ratios = true_ratios
              ) |>
                dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
            },
            error = function(e) {
              warning(paste(e, "\n"))
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
