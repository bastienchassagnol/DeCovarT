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
  .check_suggested_package("e1071", "deconvolute_ratios_cibersort")
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
  .check_suggested_package("nnls", "deconvolute_ratios_nnls")
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
  .check_suggested_package("limSolve", "deconvolute_ratios_deconrnaseq")
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

#' Fixed residual covariance for a GLS competitor
#'
#' Evaluates the convolution covariance
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}
#' at a **known** composition (default the barycentre
#' \eqn{p_j=1/J}) and optionally keeps only the diagonal. That matrix
#' does **not** depend on the unknown \eqn{\boldsymbol{p}} at fit time.
#' Do not copy it into every slice of a DeCovarT tensor: that would
#' yield \eqn{\|p\|_2^2\boldsymbol{\Sigma}_{\mathrm{GLS}}}, which still
#' depends on \eqn{p}.
#'
#' @param Sigma Array \eqn{G\times G\times J} of cell-type covariances.
#' @param p Optional composition of length \eqn{J}; default
#'   \eqn{(1/J,\ldots,1/J)}.
#' @param diagonal If `TRUE` (default), return
#'   \eqn{\operatorname{diag}\{\boldsymbol{\Sigma}(\boldsymbol{p})\}}.
#'
#' @return A \eqn{G\times G} covariance matrix.
#'
#' @examples
#' Sigma <- array(c(diag(2), 2 * diag(2)), dim = c(2, 2, 2))
#' fixed_gls_covariance(Sigma)
#' @export
#' @seealso [deconvolute_ratios_gls()]
fixed_gls_covariance <- function(Sigma, p = NULL, diagonal = TRUE) {
  Sigma <- as.array(Sigma)
  n_celltypes <- dim(Sigma)[[3L]]
  if (is.null(p)) {
    p <- rep(1 / n_celltypes, n_celltypes)
  }
  p <- as.numeric(p)
  if (length(p) != n_celltypes) {
    stop("`p` must have length equal to dim(Sigma)[3].", call. = FALSE)
  }
  w <- .compute_global_variance(p, Sigma)
  if (isTRUE(diagonal)) {
    w <- diag(diag(w), nrow = nrow(w))
  }
  w
}

#' Generalised least squares with a fixed residual covariance
#'
#' Fits
#' \eqn{\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}} by
#' [MASS::lm.gls()] with a **known** \eqn{G\times G} residual
#' covariance `W` that does not depend on \eqn{\boldsymbol{p}}, then
#' [repair_simplex()]. Use [fixed_gls_covariance()] for the
#' global-diagonal competitor
#' \eqn{\operatorname{diag}\{\sum_j \bar p_j^2\boldsymbol{\Sigma}_j\}}.
#' This is the mean-only GLS gold standard for covariance-structure
#' benchmarks, not the DeCovarT convolution
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}.
#'
#' @param y Numeric bulk vector of length \eqn{G}.
#' @param mean_signature_matrix Numeric matrix \eqn{G\times J}.
#' @param W Residual covariance (when `inverse = TRUE`) or precision
#'   (when `inverse = FALSE`), of size \eqn{G\times G}.
#' @param inverse Passed to [MASS::lm.gls()]; `TRUE` means `W` is a
#'   covariance matrix.
#'
#' @return Named simplex vector of length \eqn{J}.
#'
#' @examples
#' mu <- cbind(ct1 = c(20, 40), ct2 = c(40, 20))
#' y <- drop(mu %*% c(0.4, 0.6))
#' Sigma <- array(c(diag(2), 2 * diag(2)), dim = c(2, 2, 2))
#' deconvolute_ratios_gls(y, mu, fixed_gls_covariance(Sigma))
#' @export
#' @seealso [fixed_gls_covariance()], [deconvolute_ratios_deconrnaseq()]
deconvolute_ratios_gls <- function(
  y,
  mean_signature_matrix,
  W,
  inverse = TRUE
) {
  y <- as.numeric(y)
  x <- as.matrix(mean_signature_matrix)
  n_genes <- nrow(x)
  if (!is.matrix(W) || nrow(W) != n_genes || ncol(W) != n_genes) {
    stop(
      "`W` must be a G x G matrix (fixed residual covariance), ",
      "not a cell-type tensor.",
      call. = FALSE
    )
  }
  ct_names <- colnames(x)
  if (is.null(ct_names)) {
    ct_names <- paste0("ct", seq_len(ncol(x)))
    colnames(x) <- ct_names
  }
  dat <- data.frame(y = y, as.data.frame(x), check.names = TRUE)
  form <- stats::reformulate(
    names(dat)[-1L],
    response = "y",
    intercept = FALSE
  )
  fit <- MASS::lm.gls(
    formula = form,
    data = dat,
    W = W,
    inverse = inverse
  )
  estimated_p <- as.numeric(stats::coef(fit))
  names(estimated_p) <- ct_names
  repair_simplex(estimated_p)
}
