############################################################################
############################################################################
###                                                                      ###
###                    OTHER DECONVOLUTION ALGORITHMS                    ###
###                                                                      ###
############################################################################
############################################################################

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Linear baseline
#'   \eqn{\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}} via nu-SVR
#'   (CIBERSORT-style); no covariance prior is used.

deconvolute_ratios_CIBERSORT <- function(
  y,
  mean_signature_matrix,
  true_ratios = NULL
) {
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
#' @describeIn deconvolute_ratios_Marquardt_Levenberg Ordinary least squares
#'   for \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([stats::lsfit()]), following \insertCite{abbas_etal09;textual}{DeCovarT};
#'   estimates are projected back onto the simplex.

deconvolute_ratios_abbas <- function(
  y,
  mean_signature_matrix,
  true_ratios = NULL
) {
  estimated_p <- stats::lsfit(
    mean_signature_matrix,
    y,
    intercept = F
  )$coefficients

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

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Robust linear model
#'   \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([MASS::rlm()]), as in \insertCite{monaco_etal19;textual}{DeCovarT}.

deconvolute_ratios_monaco <- function(
  y,
  mean_signature_matrix,
  true_ratios = NULL
) {
  estimated_p <- MASS::rlm(
    y ~ mean_signature_matrix + 0,
    method = c("M")
  )$coefficients
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

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Non-negative least squares
#'   for \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([nnls::nnls()]), then simplex projection.

deconvolute_ratios_nnls <- function(
  y,
  mean_signature_matrix,
  true_ratios = NULL
) {
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


#' @describeIn deconvolute_ratios_Marquardt_Levenberg Equality- and
#'   inequality-constrained least squares on the simplex
#'   ([limSolve::lsei()]), in the spirit of `deconRNASeq`.
#' @references
#' \insertAllCited{}

deconvolute_ratios_deconRNASeq <- function(
  y,
  mean_signature_matrix,
  true_ratios = NULL
) {
  EE <- rep(1, ncol(mean_signature_matrix))
  FF <- 1 # encode the sum-to-one constraint
  GG <- diag(nrow = ncol(mean_signature_matrix))
  HH <- rep(0, ncol(mean_signature_matrix)) # encode the non-negativity constraint

  estimated_p <- limSolve::lsei(
    mean_signature_matrix,
    y,
    EE,
    FF,
    GG,
    HH
  )$mean_signature_matrix |>
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
