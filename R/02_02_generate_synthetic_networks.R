## =========================================================================
## simulate_synthetic_GRNs.R
##
## Hierarchical generative model for Gaussian-based gene regulatory
## networks with structured mean profiles, negative-binomial-like
## marginal variances, and graph-constrained precision matrices.
## =========================================================================

#' Generate hierarchical mean profiles for parent and child populations
#'
#' Constructs two complementary parent mean vectors — one with high
#' expression in the first gene block and background in the second, the
#' other with the inverse pattern — then perturbs each parent mean to
#' produce two child subpopulations.
#'
#' @param n_expressed_genes Integer, number of expressed genes per block.
#' @param mean_lower_expressed,mean_upper_expressed Numeric bounds for
#'   the uniform distribution of expressed-gene means.
#' @param mean_lower_background,mean_upper_background Numeric bounds for
#'   the uniform distribution of background-gene means.
#' @param child_perturbation_sd Numeric, standard deviation of the
#'   centred Gaussian perturbation applied to parent means.
#' @param gene_names Character vector of gene identifiers.
#' @param parent_names Character vector of length 2.
#' @param child_names Character vector of length 4.
#'
#' @return A named list with elements \code{parent_means} (2 x p matrix)
#'   and \code{child_means} (4 x p matrix).
#'
#' @keywords internal
#' @noRd
generate_mean_profiles <- function(
  n_expressed_genes,
  mean_lower_expressed,
  mean_upper_expressed,
  mean_lower_background,
  mean_upper_background,
  child_perturbation_sd,
  gene_names,
  parent_names,
  child_names
) {
  n_total_genes <- 2L * n_expressed_genes
  block_1 <- seq_len(n_expressed_genes)
  block_2 <- n_expressed_genes + seq_len(n_expressed_genes)

  parent_1_mean <- numeric(n_total_genes)
  parent_2_mean <- numeric(n_total_genes)

  parent_1_mean[block_1] <- stats::runif(
    n_expressed_genes,
    mean_lower_expressed,
    mean_upper_expressed
  )
  parent_1_mean[block_2] <- stats::runif(
    n_expressed_genes,
    mean_lower_background,
    mean_upper_background
  )
  parent_2_mean[block_1] <- stats::runif(
    n_expressed_genes,
    mean_lower_background,
    mean_upper_background
  )
  parent_2_mean[block_2] <- stats::runif(
    n_expressed_genes,
    mean_lower_expressed,
    mean_upper_expressed
  )

  parent_means <- rbind(parent_1_mean, parent_2_mean)
  rownames(parent_means) <- parent_names
  colnames(parent_means) <- gene_names

  child_means <- matrix(
    NA_real_,
    nrow = length(child_names),
    ncol = n_total_genes,
    dimnames = list(child_names, gene_names)
  )

  parent_index <- c(1L, 1L, 2L, 2L)
  for (i in seq_along(child_names)) {
    perturbation <- stats::rnorm(
      n_total_genes,
      mean = 0,
      sd = child_perturbation_sd
    )
    child_means[i, ] <- parent_means[parent_index[i], ] + perturbation
  }

  child_means <- pmax(child_means, .Machine$double.eps)

  list(
    parent_means = parent_means,
    child_means = child_means
  )
}


#' Compute gene-wise marginal variances via a NB-like mean--variance law
#'
#' Applies \eqn{\sigma_g^2 = \mu_g + \mu_g^\alpha / L} where \eqn{L}
#' is the library size and \eqn{\alpha} the power parameter.  When
#' \eqn{\alpha = 2} this reduces to the classical negative-binomial
#' variance--mean relationship.
#'
#' @param mean_matrix Numeric matrix (populations x genes).
#' @param library_size Positive numeric scalar.
#' @param alpha Positive numeric scalar, power parameter.
#'
#' @return A numeric matrix of the same dimensions as \code{mean_matrix}.
#'
#' @keywords internal
#' @noRd
compute_marginal_variances <- function(mean_matrix, library_size, alpha) {
  variance_matrix <- mean_matrix + (mean_matrix^alpha) / library_size
  variance_matrix <- pmax(variance_matrix, .Machine$double.eps)
  dimnames(variance_matrix) <- dimnames(mean_matrix)
  variance_matrix
}


#' Sample a symmetric adjacency matrix from a random graph model
#'
#' Wraps \pkg{igraph} generators for either a Barabási--Albert
#' (preferential-attachment / power-law) graph or a stochastic block
#' model.
#'
#' @param n_genes Integer, number of nodes (genes).
#' @param graph_model Character, one of \code{"power_law"} or
#'   \code{"stochastic_block_model"}.
#' @param graph_params Named list of model-specific parameters (see
#'   \code{\link{simulate_hierarchical_grn_moments}} for details).
#'
#' @return A symmetric integer matrix of dimension \code{n_genes x
#'   n_genes} with zero diagonal.
#'
#' @keywords internal
#' @noRd
generate_random_network_skeleton <- function(
  n_genes,
  graph_model,
  graph_params
) {
  graph <- switch(
    graph_model,

    power_law = {
      power <- if (is.null(graph_params$power)) {
        1
      } else {
        graph_params$power
      }
      edges_per_node <- if (is.null(graph_params$edges_per_node)) {
        1L
      } else {
        as.integer(graph_params$edges_per_node)
      }
      igraph::sample_pa(
        n_genes,
        power = power,
        m = edges_per_node,
        directed = FALSE
      )
    },

    stochastic_block_model = {
      block_prob <- if (is.null(graph_params$block_prob)) {
        c(0.5, 0.25, 0.25)
      } else {
        graph_params$block_prob
      }
      p_within <- if (is.null(graph_params$p_within)) {
        0.25
      } else {
        graph_params$p_within
      }
      p_between <- if (is.null(graph_params$p_between)) {
        0.01
      } else {
        graph_params$p_between
      }

      n_blocks <- length(block_prob)
      pref_matrix <- matrix(p_between, n_blocks, n_blocks)
      diag(pref_matrix) <- p_within

      block_sizes <- as.integer(
        c(stats::rmultinom(1L, n_genes, block_prob))
      )
      block_sizes <- pmax(block_sizes, 1L)
      size_diff <- n_genes - sum(block_sizes)
      if (size_diff != 0L) {
        idx <- which.max(block_sizes)
        block_sizes[idx] <- block_sizes[idx] + size_diff
      }

      igraph::sample_sbm(
        n_genes,
        pref.matrix = pref_matrix,
        block.sizes = block_sizes
      )
    }
  )

  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(igraph::as_undirected(graph))
  )
  diag(adjacency) <- 0L
  adjacency
}


#' Build a normalised precision matrix from a binary adjacency matrix
#'
#' Applies an affine transformation to the adjacency matrix so that the
#' resulting precision matrix is symmetric positive-definite:
#' \deqn{
#'   \tilde\Omega = A \odot v, \quad
#'   \Omega = \tilde\Omega + (\lvert\lambda_{\min}(\tilde\Omega)\rvert + u)\,I
#' }
#'
#' @param adjacency_matrix Symmetric binary matrix.
#' @param precision_shift Positive numeric (\eqn{u}), additive diagonal
#'   shift above the absolute smallest eigenvalue.
#' @param precision_scale Positive numeric (\eqn{v}), multiplicative
#'   scale for off-diagonal entries.
#'
#' @return A symmetric positive-definite matrix of the same dimension.
#'
#' @keywords internal
#' @noRd
build_normalised_precision <- function(
  adjacency_matrix,
  precision_shift,
  precision_scale
) {
  omega <- adjacency_matrix * precision_scale
  eigenvalues <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
  diag(omega) <- abs(min(eigenvalues)) + precision_shift
  omega
}


#' Scale a shared correlation structure by population-specific variances
#'
#' Given a normalised precision matrix \eqn{\Omega} (shared across
#' populations), computes the correlation matrix
#' \eqn{R = \mathrm{cov2cor}(\Omega^{-1})} and then builds
#' population-specific covariance matrices
#' \eqn{\Sigma_k = D_k^{1/2}\,R\,D_k^{1/2}} where
#' \eqn{D_k = \mathrm{diag}(\sigma_{k,1}^2, \ldots, \sigma_{k,p}^2)}.
#'
#' @param normalised_precision Symmetric positive-definite matrix
#'   (p x p).
#' @param variance_matrix Numeric matrix (n_populations x p) of
#'   marginal variances.
#' @param population_names Character vector of population identifiers.
#' @param gene_names Character vector of gene identifiers.
#'
#' @return A three-dimensional numeric array of dimension
#'   (p x p x n_populations).
#'
#' @keywords internal
#' @noRd
build_precision_matrix <- function(
  normalised_precision,
  variance_matrix,
  population_names,
  gene_names
) {
  n_genes <- nrow(normalised_precision)
  n_populations <- nrow(variance_matrix)

  correlation_matrix <- stats::cov2cor(solve(normalised_precision))

  covariance_array <- array(
    NA_real_,
    dim = c(n_genes, n_genes, n_populations),
    dimnames = list(gene_names, gene_names, population_names)
  )

  for (k in seq_len(n_populations)) {
    sd_vector <- sqrt(variance_matrix[k, ])
    scale_matrix <- diag(sd_vector, nrow = n_genes)
    cov_k <- scale_matrix %*% correlation_matrix %*% scale_matrix
    cov_k <- (cov_k + t(cov_k)) / 2
    min_ev <- min(eigen(cov_k, symmetric = TRUE, only.values = TRUE)$values)
    if (min_ev < .Machine$double.eps) {
      diag(cov_k) <- diag(cov_k) + (abs(min_ev) + .Machine$double.eps)
    }
    covariance_array[,, k] <- cov_k
  }

  covariance_array
}


# ---- main exported function ----------------------------------------------

#' Simulate hierarchical gene regulatory network moments
#'
#' Generates a complete set of first- and second-order moments for a
#' hierarchical cell-population model.  Two parent cell lines are
#' defined with complementary expression profiles (block structure);
#' each parent gives rise to two child subpopulations via small Gaussian
#' perturbations.  Gene--gene dependencies are encoded in a
#' graph-constrained precision matrix, and marginal variances follow a
#' negative-binomial-like mean--variance law.
#'
#' The simulation proceeds in four stages:
#' \enumerate{
#'   \item \strong{Mean profiles.}
#'     Parent means are drawn from uniform distributions with distinct
#'     ranges for expressed and background genes.
#'     Child means are obtained by adding a centred Gaussian perturbation
#'     \eqn{\delta^{(k)} \sim \mathcal{N}(0, \sigma_\delta^2 I)} to
#'     the corresponding parent mean.
#'   \item \strong{Marginal variances.}
#'     Gene-wise variances are computed via
#'     \eqn{\sigma_g^2 = \mu_g + \mu_g^\alpha / L}, inspiring from the
#'      mean-variance relationship of a Negative Binomial distirbution
#'   \item \strong{Graph topology.}
#'     A random graph is sampled (power-law or stochastic block model)
#'     and converted to a binary adjacency matrix \eqn{A}.
#'   \item \strong{Covariance matrices.}
#'     A normalised precision matrix is built as
#'     \eqn{\Omega = A \odot v + (\lvert\lambda_{\min}\rvert + u)\,I},
#'     then inverted and scaled by population-specific variances to
#'     yield full covariance matrices.
#' }
#'
#' @param n_expressed_genes Integer, strictly positive.  Number of genes
#'   considered expressed per parent population.  The total gene count is
#'   \code{2 * n_expressed_genes} (one expressed block + one background
#'   block per parent).
#' @param mean_lower_expressed Numeric.  Lower bound of the uniform
#'   distribution for expressed-gene means.
#' @param mean_upper_expressed Numeric.  Upper bound of the uniform
#'   distribution for expressed-gene means.
#' @param mean_lower_background Numeric.  Lower bound of the uniform
#'   distribution for background-gene means.
#' @param mean_upper_background Numeric.  Upper bound of the uniform
#'   distribution for background-gene means.
#' @param library_size Numeric, positive.  Scaling factor in the
#'   mean--variance relationship
#'   \eqn{\sigma_g^2 = \mu_g + \mu_g^\alpha / L}.
#' @param alpha Numeric, positive.  Power parameter in the
#'   mean--variance relationship.  Setting \code{alpha = 2} recovers
#'   the classical negative-binomial variance.
#' @param precision_shift Numeric, positive.  Additive diagonal shift
#'   \eqn{u} ensuring positive-definiteness of the precision matrix.
#' @param precision_scale Numeric, positive.  Multiplicative scale
#'   \eqn{v} applied to off-diagonal precision entries.
#' @param child_perturbation_sd Numeric, positive (typically small).
#'   Standard deviation of the centred Gaussian perturbation added to
#'   parent means to create child subpopulation means.
#' @param graph_model Character.  Random graph generator to use; one of
#'   \code{"power_law"} (Barabási--Albert preferential attachment) or
#'   \code{"stochastic_block_model"}.
#' @param graph_params Named list of model-specific parameters:
#'   \describe{
#'     \item{For \code{"power_law"}:}{
#'       \code{power} (numeric, default 1) — attachment exponent;
#'       \code{edges_per_node} (integer, default 1) — edges added per
#'       new node.
#'     }
#'     \item{For \code{"stochastic_block_model"}:}{
#'       \code{block_prob} (numeric vector, default
#'       \code{c(0.5, 0.25, 0.25)}) — block-membership probabilities;
#'       \code{p_within} (numeric, default 0.25) — within-block edge
#'       probability;
#'       \code{p_between} (numeric, default 0.01) — between-block edge
#'       probability.
#'     }
#'   }
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{\code{parent_parameters}}{A named list containing:
#'       \describe{
#'         \item{\code{mean_profiles}}{Numeric matrix
#'           (2 x \code{2 * n_expressed_genes}).  Rows: \code{parent_1},
#'           \code{parent_2}; columns: \code{gene_1}, \ldots}
#'         \item{\code{covariance_matrices}}{Three-dimensional array
#'           (\code{2n x 2n x 2}), third axis indexed by parent name.}
#'       }
#'     }
#'     \item{\code{child_parameters}}{A named list containing:
#'       \describe{
#'         \item{\code{mean_profiles}}{Numeric matrix
#'           (4 x \code{2 * n_expressed_genes}).  Rows:
#'           \code{parent_1_child_a}, \code{parent_1_child_b},
#'           \code{parent_2_child_a}, \code{parent_2_child_b}.}
#'         \item{\code{covariance_matrices}}{Three-dimensional array
#'           (\code{2n x 2n x 4}), third axis indexed by child name.}
#'       }
#'     }
#'     \item{\code{graph_structure}}{A named list containing:
#'       \describe{
#'         \item{\code{adjacency_matrix}}{Binary symmetric matrix
#'           encoding the GRN topology.}
#'         \item{\code{normalised_precision}}{Positive-definite precision
#'           matrix derived from the adjacency matrix.}
#'       }
#'     }
#'   }
#'
#' @examples
#' set.seed(42)
#' moments <- simulate_hierarchical_grn_moments(
#'     n_expressed_genes     = 50,
#'     mean_lower_expressed  = 2,
#'     mean_upper_expressed  = 6,
#'     mean_lower_background = 0.1,
#'     mean_upper_background = 0.5,
#'     library_size          = 10000,
#'     alpha                 = 2,
#'     precision_shift       = 0.1,
#'     precision_scale       = 0.3,
#'     child_perturbation_sd = 0.1,
#'     graph_model           = "power_law",
#'     graph_params          = list(power = 1, edges_per_node = 2)
#' )
#' str(moments, max.level = 2)
#'
#' ## Verify positive-definiteness of a child covariance
#' eigen_vals <- eigen(
#'     moments$child_parameters$covariance_matrices[, , 1],
#'     only.values = TRUE
#' )$values
#' stopifnot(all(eigen_vals > 0))
#'
#' @importFrom igraph sample_pa sample_sbm as_adjacency_matrix
#'   as_undirected
#' @importFrom stats runif rnorm rmultinom cov2cor
#'
#' @export
simulate_hierarchical_grn_moments <- function(
  n_expressed_genes,
  mean_lower_expressed,
  mean_upper_expressed,
  mean_lower_background,
  mean_upper_background,
  library_size,
  alpha,
  precision_shift,
  precision_scale,
  child_perturbation_sd,
  graph_model = c("power_law", "stochastic_block_model"),
  graph_params = list()
) {
  ## --- input validation ------------------------------------------------

  graph_model <- match.arg(graph_model)

  if (
    !is.numeric(n_expressed_genes) ||
      length(n_expressed_genes) != 1L ||
      n_expressed_genes < 2L
  ) {
    stop("'n_expressed_genes' must be a single integer >= 2.", call. = FALSE)
  }
  n_expressed_genes <- as.integer(n_expressed_genes)

  if (mean_lower_expressed >= mean_upper_expressed) {
    stop(
      "'mean_lower_expressed' must be strictly less than ",
      "'mean_upper_expressed'.",
      call. = FALSE
    )
  }
  if (mean_lower_background >= mean_upper_background) {
    stop(
      "'mean_lower_background' must be strictly less than ",
      "'mean_upper_background'.",
      call. = FALSE
    )
  }
  if (library_size <= 0) {
    stop("'library_size' must be positive.", call. = FALSE)
  }
  if (alpha <= 0) {
    stop("'alpha' must be positive.", call. = FALSE)
  }
  if (precision_shift <= 0 || precision_scale <= 0) {
    stop(
      "'precision_shift' and 'precision_scale' must both be ",
      "positive.",
      call. = FALSE
    )
  }
  if (child_perturbation_sd <= 0) {
    stop("'child_perturbation_sd' must be positive.", call. = FALSE)
  }

  ## --- identifiers -----------------------------------------------------

  n_total_genes <- 2L * n_expressed_genes
  gene_names <- paste0("gene_", seq_len(n_total_genes))
  parent_names <- c("parent_1", "parent_2")
  child_names <- c(
    "parent_1_child_a",
    "parent_1_child_b",
    "parent_2_child_a",
    "parent_2_child_b"
  )

  ## --- step 1: hierarchical mean profiles ------------------------------

  mean_profiles <- generate_mean_profiles(
    n_expressed_genes = n_expressed_genes,
    mean_lower_expressed = mean_lower_expressed,
    mean_upper_expressed = mean_upper_expressed,
    mean_lower_background = mean_lower_background,
    mean_upper_background = mean_upper_background,
    child_perturbation_sd = child_perturbation_sd,
    gene_names = gene_names,
    parent_names = parent_names,
    child_names = child_names
  )

  ## --- step 2: NB-like marginal variances ------------------------------

  parent_variances <- compute_marginal_variances(
    mean_matrix = mean_profiles$parent_means,
    library_size = library_size,
    alpha = alpha
  )
  child_variances <- compute_marginal_variances(
    mean_matrix = mean_profiles$child_means,
    library_size = library_size,
    alpha = alpha
  )

  ## --- step 3: random graph adjacency matrix ---------------------------

  adjacency_matrix <- generate_random_network_skeleton(
    n_genes = n_total_genes,
    graph_model = graph_model,
    graph_params = graph_params
  )
  rownames(adjacency_matrix) <- gene_names
  colnames(adjacency_matrix) <- gene_names

  ## --- step 4: normalised precision → covariance arrays ----------------

  normalised_precision <- build_normalised_precision(
    adjacency_matrix = adjacency_matrix,
    precision_shift = precision_shift,
    precision_scale = precision_scale
  )

  parent_covariances <- build_precision_matrix(
    normalised_precision = normalised_precision,
    variance_matrix = parent_variances,
    population_names = parent_names,
    gene_names = gene_names
  )

  child_covariances <- build_precision_matrix(
    normalised_precision = normalised_precision,
    variance_matrix = child_variances,
    population_names = child_names,
    gene_names = gene_names
  )

  ## --- assemble output -------------------------------------------------

  list(
    parent_parameters = list(
      mean_profiles = mean_profiles$parent_means,
      covariance_matrices = parent_covariances
    ),
    child_parameters = list(
      mean_profiles = mean_profiles$child_means,
      covariance_matrices = child_covariances
    ),
    graph_structure = list(
      adjacency_matrix = adjacency_matrix,
      normalised_precision = normalised_precision
    )
  )
}
