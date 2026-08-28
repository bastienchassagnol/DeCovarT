#' Simulate bulk mixtures from a multivariate Gaussian convolution
#'
#' @description
#' For each bootstrap sample \eqn{i=1,\ldots,N}, draws **latent** purified
#' profiles
#' \eqn{\boldsymbol{x}_{\cdot j}^{(i)}\sim
#'   \mathcal{N}_{G}(\boldsymbol{\mu}_{\cdot j},
#' \boldsymbol{\Sigma}_j)} independently for each cell type
#' \eqn{j=1,\ldots,J}, then forms the bulk by the linear convolution
#' \deqn{
#'   \boldsymbol{y}_{\cdot i}
#'   =\boldsymbol{X}^{(i)}\boldsymbol{p}
#'   =\sum_{j=1}^{J}p_j\,\boldsymbol{x}_{\cdot j}^{(i)},
#' }
#' matching the article's conditional model
#' \eqn{\boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))} with
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=
#' \sum_j p_j^{2}\boldsymbol{\Sigma}_j}.
#'
#' The mean signature \eqn{\boldsymbol{\mu}} is shared across samples; only the
#' latent draws \eqn{\boldsymbol{x}_{\cdot j}^{(i)}} vary with \eqn{i}. Stacking
#' those draws into the three-way array
#' \eqn{\mathcal{X}=(x_{gji})\in\mathcal{M}_{G\times J\times N}}, the bulk
#' matrix is the mode-2 tensor–vector contraction
#' \deqn{
#'   \boldsymbol{Y}
#'   =\mathcal{X}\times_{2}\boldsymbol{p},
#'   \qquad
#'   y_{gi}=\sum_{j=1}^{J}x_{gji}\,p_{j}
#'   \quad(g=1,\ldots,G;\; i=1,\ldots,N),
#' }
#' which for each sample recovers
#' \eqn{\boldsymbol{y}_{\cdot i}=\boldsymbol{X}^{(i)}\,\boldsymbol{p}}
#' with \eqn{\boldsymbol{X}^{(i)}\in\mathcal{M}_{G\times J}}.
#'
#' @param signature_matrix Mean signature
#'   \eqn{\boldsymbol{\mu}=(\mu_{gj})\in\mathcal{M}_{G\times J}} (shared across
#'   samples; not the latent profiles).
#' @param Sigma Array of covariances
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}.
#' @param p Proportion vector \eqn{\boldsymbol{p}\in\Delta^{J-1}}
#'   (default: uniform).
#' @param n Number of bulk / bootstrap samples \eqn{N}.
#'
#' @return A list with:
#' * `latent_profiles`: array
#'   \eqn{\mathcal{X}=(x_{gji})\in\mathcal{M}_{G\times J\times N}} of
#'   **unobserved** cell-type-specific draws (one \eqn{G\times J} slice per
#'   sample \eqn{i});
#' * `Y`: matrix \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}} whose columns
#'   are bulk vectors \eqn{\boldsymbol{y}_{\cdot i}}, obtained as
#'   \eqn{\mathcal{X}\times_{2}\boldsymbol{p}}.
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
#' sim <- simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 5)
#' dim(sim$Y)
#' @export
#' @seealso [deconvolute_ratios()],
#'   [run_simulation_benchmark()]
simulate_bulk_mixture <- function(
  signature_matrix,
  Sigma,
  p = rep(1 / ncol(signature_matrix), ncol(signature_matrix)),
  n = 500
) {
  ##################################################################
  ##                        check validity                        ##
  ##################################################################
  gene_ok <- isTRUE(all.equal(
    rownames(signature_matrix),
    union(dimnames(Sigma)[[1]], dimnames(Sigma)[[2]])
  ))
  cell_ok <- isTRUE(all.equal(
    colnames(signature_matrix),
    dimnames(Sigma)[[3]]
  ))
  if (!gene_ok) {
    stop(
      "Some of the genes are distinct between expected",
      " and covariance expression"
    )
  }
  if (!cell_ok) {
    stop("Cell types differ between expected and covariance expression")
  }

  ##### simulation
  valid_celltypes <- colnames(signature_matrix)
  names_genes <- row.names(signature_matrix)
  Y <- matrix(
    0,
    nrow = nrow(signature_matrix),
    ncol = n,
    dimnames = list(names_genes, paste0("sample_", 1:n))
  )
  # Latent purified profiles X: G x J x N (not the mean signature mu)
  latent_profiles <- array(
    0,
    c(nrow(signature_matrix), ncol(signature_matrix), n),
    dimnames = list(names_genes, valid_celltypes, colnames(Y))
  )

  for (cell_name in valid_celltypes) {
    # Draw x_{·j}^{(i)} ~ N(mu_{·j}, Sigma_j) for i = 1, ..., N
    mean_parameter <- signature_matrix[, cell_name]
    covariance_parameter <- Sigma[,, cell_name]
    expression_per_celltype <- MASS::mvrnorm(
      n = n,
      mu = mean_parameter,
      Sigma = covariance_parameter,
      tol = 1e-12,
      empirical = FALSE
    )
    latent_profiles[, cell_name, ] <- t(expression_per_celltype)
  }

  ## Mode-2 contraction: Y = X ×_2 p, i.e. y_{gi} = Σ_j x_{gji} p_j
  Y <- tensor::tensor(p, B = latent_profiles, alongA = 1, alongB = 2)
  return(list(latent_profiles = latent_profiles, Y = Y))
}
