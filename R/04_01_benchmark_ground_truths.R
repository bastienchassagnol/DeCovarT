

#' @title Compute summary metrics evaluating the quality of the estimate
#'
#' @description Compute metrics, either comparing th estimated ratios with a gold standard
#' or the divergence between the reconstituted virtual mixture,
#' using deterministic rule \eqn{\boldsymbol{\hat{y}}=\boldsymbol{mean_signature_matrix} \times \hat{\boldsymbol{p}}}
#' and the actual measured one
#'
#' @inheritParams deconvolute_ratios_Marquardt_Levenberg
#' @param estimated_p The ratios estimated by your favourite deconvolution algorithm
#'
#' @return A `tibble`, with the following scores:
#' * mse and rmse, for respectively \emph{mean} and \emph{root mean squared error}. See also the [Metrics::mse()] function.
#' * mae, for \emph{mean absolute error}. See also the [Metrics::mae()] function.
#' * \eqn{R^2} and adjusted \eqn{R^2}, corresponding to the percentage of variance
#' captured by the linear regression model. See also the [Metrics::rse()] function.
#' * cor, for the Pearson correlation between the estimated and true cellular ratios
#' giving the mean values of the variables within a given component. See also the [stats::cor()] function.
#' @export

compute_benchmark_metrics <- function(y, mean_signature_matrix, estimated_p, true_ratios = NULL) {
  n <- nrow(mean_signature_matrix)
  k <- ncol(mean_signature_matrix)
  df_res <- n - k + 1 # number of free parameters: only k - 1 parameters must be learnt, with sum-to-one constraint
  df_tot <- n - 1 # no intercept for the moment, in the model (so n-1, or n?)
  if (!is.null(true_ratios)) {
    # when the true parameters are known
    scores <- tibble::tibble(
      model_mse = Metrics::mse(true_ratios, estimated_p),
      model_rmse = Metrics::rmse(true_ratios, estimated_p),
      model_mae = Metrics::mae(true_ratios, estimated_p),
      model_coef_determination = max(
        0,
        1 - Metrics::rse(true_ratios, estimated_p)
      ),
      model_coef_determination_adjusted = max(
        0,
        1 - (1 - .data$model_coef_determination) * df_tot / df_res
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
    
    # invisible(utils::capture.output
    # model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_p)),
    # model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
    # # model_ccc= epiR::epi.ccc(y, predicted_values)
  }
  return(scores)
}

#' Main function of the package: deconvolute in parallel mixture samples
#'
#' @author Bastien CHASSAGNOL
#'
#' @param signature_matrix Parameter `mu`: \eqn{\boldsymbol{\mu}=(\mu_{g,j}) \in \mathbb{R}^{G \times J}},
#' storing in each each column the averaged expression of the `G` genes used to deconvolve all cell populations.
#' Name your `colnames` as the cell populations, and provide `rownames` argument for the name of the genes.
#' By convention, we use HGNC symbols.
#' @param bulk_expression Parameter `y`: \eqn{\boldsymbol{y}=(\mu_{g,i}) \in \mathbb{R}^{G \times I}},
#' storing in each each column the measured expression of the `G` genes in a heterogeneous sample, using any RNASeq or microarray technology.
#' Provide the sample ID for each of your samples in column, and the name of your genes in `rownames`.
#' @param true_ratios If available (for instance, in the context of a virtual benchmark, or if some standard cytometry techniques provide them),
#' vector of size \eqn{J}, storing the normalised proportions of the cell populations supposed present in the sample. Summary metrics
#' will then be computed against the ones returned by the deconvolution algorithms provided.
#' @param Sigma Only relevant for deconvolution algorithms which require a prior estimate
#' of the transcriptomic  covariance for each of the purified cell populations.
#' A 3-dimensional covariance matrix array is expected:
#'  \eqn{\mathrm{\Sigma}=(\Sigma_{l, k, j}) \in \mathbb{R}^{G \times G \times J}} to
#' parametrise the covariance transcriptomic structure of the \eqn{J} cell populations estimated.
#' @param deconvolution_functions The deconvolution functions themselves, a list
#' with for each item two attributes to be filled with:
#' * `FUN`: the function itself (not a string, but indeed any deconvolution function
#' integrating the default parameters listed in )
#' `additional_parameters` by default, set to NULL. If your deconvolution function
#' integrates any specific, additional parameter.
#' @param scaled Whether we should scale or not the dataset. By default, we consider that the provided dataset is in its original raw space,
#' and we do not scale the dataset, since our deconvolution algorithm assumes a multivariate Gaussian distribution on the raw counts themselves.
#' @param cores For a parallel estimation of ratios in a series of bulk samples,
#' assign a number of cores strictly inferior to the number of cores available
#' on your machine. By default, the maximum, minus in Unix systems, and for
#' OS compatibility, only one on Windows machines
#' @return A `tibble` storing for each row the measured cell proportions, as well as some summary metrics.
#' We ensure for each deconvolution algorithm that the returned estimates respect the unit simplex constraint,
#' with function [enforce_identifiability()].
#' @export
#'
#' @seealso [deconvolute_ratios_Marquardt_Levenberg()], to deconvolve a single, already normalised sample

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
  if (!is.matrix(signature_matrix) | is.null(row.names(signature_matrix))) {
    stop(
      "required format for signature is expression matrix, with rownames as genes"
    )
  }
  mean_signature_matrix <- tibble::as_tibble(signature_matrix, rownames = "GENE_SYMBOL")
  if (!is.matrix(bulk_expression) | is.null(row.names(bulk_expression))) {
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
  }
  mean_signature_matrix <- log2(mean_signature_matrix)
  
  # estimation itself
  deconvolution_estimates <- purrr::imap_dfr(
    deconvolution_functions,
    function(deconvolution_function, algorithm) {
      additional_parameters <- deconvolution_function$additional_parameters
      metric_scores <- parallel::mclapply(
        1:ncol(Y),
        function(i) {
          # metric_scores <- tibble::tibble(); for (i in 1:ncol(Y)) {
          success_estimation <- tryCatch(
            {
              list_arguments <- c(
                list(
                  "y" = Y[, i] |> as.matrix(),
                  "mean_signature_matrix" = mean_signature_matrix,
                  "Sigma" = Sigma,
                  "true_ratios" = true_ratios
                ),
                additional_parameters
              )
              # we only keep the arguments needed by the required function
              estimated_p <- do.call(
                deconvolution_function$FUN,
                list_arguments[methods::formalArgs(deconvolution_function$FUN)]
              )
            },
            error = function(e) {
              # dir.create("/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs", showWarnings = F)
              warning(paste(e, "\n"))
              dir.create(
                "./simulations/erreurs",
                showWarnings = F,
                recursive = TRUE
              )
              # "y"=Y[,i] |> as.matrix(), "mean_signature_matrix"=mean_signature_matrix, "Sigma"=Sigma,
              saveRDS(
                c(
                  list_arguments[methods::formalArgs(
                    deconvolution_function$FUN
                  )],
                  list(
                    "CALL" = call(
                      algorithm,
                      list_arguments[methods::formalArgs(
                        deconvolution_function$FUN
                      )]
                    ),
                    "error" = e,
                    "function_name" = algorithm
                  )
                ),
                file = paste0(
                  "/home/bncl_cb/rstudio/working/DeCovarT/simulations/erreurs/erreur_",
                  i,
                  "_function_",
                  algorithm,
                  ".rds"
                )
              )
              return(e)
            }
          )
          if (!inherits(success_estimation, "error")) {
            return(estimated_p)
          }
        },
        mc.cores = cores
      ) |>
        dplyr::bind_rows() # one estimation with a given deconvolution algorithm terminated
      
      if (nrow(metric_scores) != 0) {
        metric_scores <- metric_scores |>
          dplyr::mutate(
            OMIC_ID = paste0("sample_", 1:nrow(metric_scores)),
            algorithm = algorithm
          )
      }
      return(metric_scores)
    }
  )
  return(deconvolution_estimates)
}
