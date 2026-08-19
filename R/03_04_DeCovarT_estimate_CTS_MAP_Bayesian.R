#' Maximum a posteriori purified profiles under a Gaussian convolution
#'
#' @description
#' For fixed plug-in
#' \eqn{\boldsymbol{\zeta}=(\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})} and a
#' bulk vector \eqn{\boldsymbol{y}}, returns cell-type-specific MAP estimates of
#' the **latent** purified profiles
#' \eqn{\boldsymbol{x}_{\cdot j}} in the additive Gaussian model of the article.
#' This is the Bayesian counterpart to the frequentist plug-in that replaces
#' \eqn{\boldsymbol{x}_{\cdot j}} by \eqn{\boldsymbol{\mu}_{\cdot j}} when
#' estimating proportions alone.
#'
#' @inheritParams loglik_multivariate
#'
#' @return List of length \eqn{J} with MAP vectors in \eqn{\mathbb{R}^{G}}.
#'
#' @keywords internal
#' @examples
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' y <- drop(mu %*% c(0.6, 0.4))
#' .map_gaussian_convolution(y, mu, Sigma)
#' @export
.map_gaussian_convolution <- function(y, mean_signature_matrix, Sigma) {
  J <- ncol(mean_signature_matrix)
  # compute intermediate calculations ---------------------------------------
  global_variance <- .compute_global_variance(
    p = c(1, 1),
    Sigma = Sigma
  )
  residual <- y - rowSums(mean_signature_matrix)
  global_precision_matrix <- solve(global_variance)

  # MAP, aka cell type specific expression, for each cell type --------------
  map_multivariate_gaussian <- lapply(seq_len(J), function(j) {
    mean_signature_matrix[, j, drop = FALSE] +
      Sigma[,, j] %*% global_precision_matrix %*% residual
  })

  return(map_multivariate_gaussian)
}
