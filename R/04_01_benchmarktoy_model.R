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
#' @importFrom rlang .data
#' @export
compute_benchmark_metrics <- function(
  y,
  mean_signature_matrix,
  estimated_p,
  true_ratios = NULL
) {
  n <- nrow(mean_signature_matrix)
  k <- ncol(mean_signature_matrix)
  df_res <- n - k + 1 # number of free parameters: only k - 1 parameters must be learnt, with sum-to-one constraint
  df_tot <- n - 1 # no intercept for the moment, in the model (so n-1, or n?)
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
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))}.
#'
#' @param signature_matrix Mean signature
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (rownames = genes,
#'   colnames = cell types).
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
    stop(
      "required format for signature is expression matrix, with rownames as genes"
    )
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
        dplyr::bind_rows() # one estimation with a given deconvolution algorithm terminated

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
