#' Compute the maximum a posteriori for a sum of Gaussian variables
#'
#' @param y Parameter `y`: \eqn{\boldsymbol{y}=(y_{g}) \in \mathbb{R}^{G}},
#' storing the measured expression of the `G` genes in the heterogeneous sample
#' @param mean_signature_matrix Parameter `mu`: \eqn{\boldsymbol{\mu}=(\mu_{g,j}) \in \mathbb{R}^{G \times J}},
#' storing in each each column the averaged expression of the `G` genes of the `J` cell populations.
#' @param Sigma A tensor storing for each cell type the GRN structure:
#'  \eqn{\mathrm{\Sigma}=(\Sigma_{l, k, j}) \in \mathbb{R}^{G \times G \times J}}, with each matrix
#' \eqn{\Sigma_{..j}, j \in \{ 1, \ldots, J\}} storing the covariance matrix
#' describing the covariance transcriptomic structure of a given cell population \eqn{j}

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
