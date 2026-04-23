
############################################################################
############################################################################
###                                                                      ###
###                    OTHER DECONVOLUTION ALGORITHMS                    ###
###                                                                      ###
############################################################################
############################################################################

#' @describeIn deconvolute_ratios_Marquardt_Levenberg From this algorithm, providing
#' any explicit prior of the covariance matrix is pointless, since we are,
#' to our knowledge, the first ones to account for it explicitly. Here,
#' a custom implementation of the CIBERSORT algorithm.

deconvolute_ratios_CIBERSORT <- function(y, mean_signature_matrix, true_ratios = NULL) {
  # the set of nu values tested for best performance
  range.nu <- seq(0.2, 0.8, 0.3)
  model <- e1071::best.svm(
    mean_signature_matrix,
    y,
    type = "nu-regression",
    kernel = "linear",
    scale = F,
    nu = range.nu,
    tunecontrol = e1071::tune.control(sampling = "fix", fix = 0.75)
  ) # determine best fitted model
  
  e1071::tune.svm(
    mean_signature_matrix,
    y,
    type = "nu-regression",
    kernel = "linear",
    scale = F,
    nu = 0.25,
    tunecontrol = e1071::tune.control(sampling = "fix", fix = 0.75)
  )
  
  model <- e1071::svm(mean_signature_matrix, y = y)
  
  # get and normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- t(model$coefs) %*% model$SV |> enforce_identifiability()
  names(estimated_p) <- colnames(mean_signature_matrix)
  
  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @importFrom Rdpack reprompt
#' @describeIn deconvolute_ratios_Marquardt_Levenberg  Here,
#' a standard linear approach, as performed in \insertCite{abbas_etal09;textual}{DeCovarT} as it can be computed
#' with function [stats::lsfit()]. Nevertheless, similar to any other deconvolution methods,
#' inferred ratios are normalised back to the unit simplex space.

deconvolute_ratios_abbas <- function(y, mean_signature_matrix, true_ratios = NULL) {
  estimated_p <- stats::lsfit(mean_signature_matrix, y, intercept = F)$coefficients
  
  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- estimated_p |> enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg  Here,
#' a robust linear approach, as performed in \insertCite{monaco_etal19;textual}{DeCovarT}
#' as it can be computed with function [MASS::rlm()].

deconvolute_ratios_monaco <- function(y, mean_signature_matrix, true_ratios = NULL) {
  estimated_p <- MASS::rlm(y ~ mean_signature_matrix + 0, method = c("M"))$coefficients
  names(estimated_p) <- colnames(mean_signature_matrix)
  
  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_p <- estimated_p |> enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Here, function [nnls::nnls()]
#' is used as an interface to the Lawson-Hanson NNLS implementation, with the additional
#' constraints that the returned ratios can not be negative. 

deconvolute_ratios_nnls <- function(y, mean_signature_matrix, true_ratios = NULL) {
  estimated_p <- nnls::nnls(mean_signature_matrix, y)$mean_signature_matrix
  names(estimated_p) <- colnames(mean_signature_matrix)
  
  # normalize coefficients (sum to 1, as non-negativity is already enforced)
  estimated_p <- estimated_p |> enforce_identifiability()
  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}


#' @describeIn deconvolute_ratios_Marquardt_Levenberg  Here,
#' standard linear least squares optimisation, but accounting for
#' the two equality and inequality constraints of the unit simplex. Similar to the implementation
#' of the `deconRNASeq` algorithm, see also [limSolve::lsei()]
#' or the `lsqlin` function in Matlab for additional details.
#' @references
#' \insertAllCited{}

deconvolute_ratios_deconRNASeq <- function(y, mean_signature_matrix, true_ratios = NULL) {
  EE <- rep(1, ncol(mean_signature_matrix))
  FF <- 1 # encode the sum-to-one constraint
  GG <- diag(nrow = ncol(mean_signature_matrix))
  HH <- rep(0, ncol(mean_signature_matrix)) # encode the non-negativity constraint
  
  estimated_p <- limSolve::lsei(mean_signature_matrix, y, EE, FF, GG, HH)$mean_signature_matrix |>
    enforce_identifiability()
  names(estimated_p) <- colnames(mean_signature_matrix)
  metrics_scores <- compute_benchmark_metrics(
    y,
    mean_signature_matrix,
    estimated_p,
    true_ratios
  ) |>
    dplyr::bind_cols(tibble::as_tibble_row(estimated_p))
  return(metrics_scores)
}
