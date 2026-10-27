# Reference samples for benchmark methods that estimate the reference
# moments themselves (MEAD, RNA-Sieve) instead of taking them as plug-in
# values. Drawn independently of the bulk samples, so no bulk shares a
# purified profile with the reference.

#' Simulate an independent reference panel
#'
#' @description
#' Draws, for every cell type \eqn{j}, `n_individuals` \eqn{\times}
#' `n_cells` profiles from
#' \eqn{\mathcal{N}_G(\boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_j)}
#' and labels each with its individual. With `n_cells = 1` every column
#' is the purified profile of one individual, as in the population
#' convolution; with `n_cells > 1` every individual contributes
#' `n_cells` cells, as in a single-cell reference. The generative models
#' of this package have no between-individual layer, so individuals are
#' exchangeable labels.
#'
#' Benchmark methods that estimate the reference moments (MEAD averages
#' per-individual pseudo-bulks, RNA-Sieve pools the cells) receive this
#' panel; methods that take plug-in moments receive the true
#' \eqn{(\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})}. Call it with its own
#' `seed`, distinct from the seed of the bulk samples, so that no bulk is
#' built from a profile of the panel.
#'
#' @param mean_signature_matrix \eqn{G\times J} mean matrix with gene row
#'   names and cell-type column names.
#' @param Sigma \eqn{G\times G\times J} covariance array.
#' @param n_individuals Number of individuals \eqn{R\ge2}.
#' @param n_cells Profiles per individual and cell type.
#' @param seed Optional seed; the global RNG state is restored afterwards.
#'
#' @return A list with `reference`, a list named by cell type of
#'   \eqn{G\times(R\,n_{\mathrm{cells}})} matrices, and `individual`, a list
#'   named by cell type of integer vectors giving the individual of each
#'   column.
#'
#' @examples
#' mu <- matrix(
#'   c(5, 1, 1, 1, 5, 1),
#'   nrow = 3,
#'   dimnames = list(paste0("g", 1:3), c("ct1", "ct2"))
#' )
#' sigma <- array(
#'   c(diag(3), diag(3)),
#'   dim = c(3, 3, 2),
#'   dimnames = list(rownames(mu), rownames(mu), colnames(mu))
#' )
#' panel <- simulate_reference_replicates(mu, sigma, n_individuals = 4, seed = 1)
#' dim(panel$reference$ct1)
#' panel$individual$ct1
#' @export
#' @seealso [simulate_bulk_mixture()], [simulate_bulk_from_cells()]
simulate_reference_replicates <- function(
  mean_signature_matrix,
  Sigma,
  n_individuals,
  n_cells = 1L,
  seed = NULL
) {
  mean_signature_matrix <- as.matrix(mean_signature_matrix)
  n_genes <- nrow(mean_signature_matrix)
  nms <- colnames(mean_signature_matrix)
  if (is.null(nms)) {
    stop("`mean_signature_matrix` needs cell-type column names.", call. = FALSE)
  }
  if (!identical(dim(Sigma), c(n_genes, n_genes, length(nms)))) {
    stop(
      "`Sigma` must be a G x G x J array matching `mean_signature_matrix`.",
      call. = FALSE
    )
  }
  n_individuals <- as.integer(n_individuals)
  n_cells <- as.integer(n_cells)
  if (
    length(n_individuals) != 1L || is.na(n_individuals) || n_individuals < 2L
  ) {
    stop("`n_individuals` must be an integer of at least 2.", call. = FALSE)
  }
  if (length(n_cells) != 1L || is.na(n_cells) || n_cells < 1L) {
    stop("`n_cells` must be a positive integer.", call. = FALSE)
  }
  n_columns <- n_individuals * n_cells
  individual <- rep(seq_len(n_individuals), each = n_cells)
  draw <- function(j) {
    x <- MASS::mvrnorm(
      n = n_columns,
      mu = mean_signature_matrix[, j],
      Sigma = Sigma[,, j]
    )
    x <- matrix(t(x), nrow = n_genes)
    dimnames(x) <- list(
      rownames(mean_signature_matrix),
      paste0(
        "individual",
        individual,
        "_cell",
        sequence(rep(n_cells, n_individuals))
      )
    )
    x
  }
  reference <- .with_seed(seed, lapply(seq_along(nms), draw))
  list(
    reference = stats::setNames(reference, nms),
    individual = stats::setNames(rep(list(individual), length(nms)), nms)
  )
}
