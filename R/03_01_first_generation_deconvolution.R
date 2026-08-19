############################################################################
############################################################################
###                                                                      ###
###                    First generation deconvolution algorithms         ###
###                                                                      ###
############################################################################
############################################################################

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Linear baseline
#'   \eqn{\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}} via nu-SVR
#'   (CIBERSORT-style); no covariance prior is used.
#' @export
deconvolute_ratios_cibersort <- function(y, mean_signature_matrix) {
  # nu grid as in CIBERSORT; skip tuning when too few genes for hold-out
  range_nu <- seq(0.2, 0.8, 0.3)
  n_obs <- nrow(mean_signature_matrix)
  model <- NULL
  if (n_obs >= 4L) {
    model <- tryCatch(
      e1071::best.svm(
        mean_signature_matrix,
        y,
        type = "nu-regression",
        kernel = "linear",
        scale = FALSE,
        nu = range_nu,
        tunecontrol = e1071::tune.control(sampling = "fix", fix = 0.75)
      ),
      error = function(e) NULL
    )
  }
  if (is.null(model)) {
    model <- e1071::svm(
      mean_signature_matrix,
      y = y,
      type = "nu-regression",
      kernel = "linear",
      scale = FALSE,
      nu = 0.5
    )
  }

  estimated_p <- as.numeric(t(model$coefs) %*% model$SV)
  estimated_p[estimated_p < 0] <- 0
  estimated_p <- repair_simplex(estimated_p)
  names(estimated_p) <- colnames(mean_signature_matrix)
  estimated_p
}

#' @importFrom Rdpack reprompt
#' @describeIn deconvolute_ratios_Marquardt_Levenberg Ordinary least squares
#'   for \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([stats::lsfit()]), following
#'   \insertCite{abbasDeconvolutionBloodMicroarray2009;textual}{DeCovarT};
#'   estimates are projected back onto the simplex.
#' @export
deconvolute_ratios_lsfit <- function(y, mean_signature_matrix) {
  estimated_p <- stats::lsfit(
    mean_signature_matrix,
    y,
    intercept = FALSE
  )$coefficients
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Robust linear model
#'   \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([MASS::rlm()]), as in
#'   \insertCite{monacoRNASeqSignaturesNormalized2019;textual}{DeCovarT}.
#' @export
deconvolute_ratios_rlm <- function(y, mean_signature_matrix) {
  estimated_p <- MASS::rlm(
    y ~ mean_signature_matrix + 0,
    method = "M"
  )$coefficients
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Non-negative least squares
#'   for \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}}
#'   ([nnls::nnls()]), then simplex projection.
#' @export
deconvolute_ratios_nnls <- function(y, mean_signature_matrix) {
  estimated_p <- nnls::nnls(mean_signature_matrix, as.numeric(y))$x
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Equality- and
#'   inequality-constrained least squares on the simplex
#'   ([limSolve::lsei()]), in the spirit of `deconRNASeq`.
#' @references
#' \insertAllCited{}
#' @export
deconvolute_ratios_deconrnaseq <- function(y, mean_signature_matrix) {
  n_celltypes <- ncol(mean_signature_matrix)
  estimated_p <- limSolve::lsei(
    A = mean_signature_matrix,
    B = as.numeric(y),
    E = rep(1, n_celltypes),
    F = 1,
    G = diag(nrow = n_celltypes),
    H = rep(0, n_celltypes),
    verbose = FALSE
  )$X
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}
