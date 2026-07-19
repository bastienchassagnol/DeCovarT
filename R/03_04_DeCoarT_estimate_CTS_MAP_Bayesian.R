#' Maximum a posteriori purified profiles under a Gaussian convolution
#'
#' @description
#' For fixed plug-in
#' \eqn{\boldsymbol{\zeta}=(\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})} and a
#' bulk vector \eqn{\boldsymbol{y}}, returns cell-type-specific MAP estimates of
#' the latent purified profiles in the additive Gaussian model of the article.
#'
#' @param y Bulk vector \eqn{\boldsymbol{y}\in\mathbb{R}^{G}}.
#' @param mean_signature_matrix Mean matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}.
#' @param Sigma Array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}.
#'
#' @return List of length \eqn{J} with MAP vectors in \eqn{\mathbb{R}^{G}}.
#'
#' @keywords internal
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
