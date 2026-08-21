#' Repair a numeric vector onto the unit simplex
#'
#' @description
#' Clips numerical under/overflow and renormalises so that
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1} and
#' \eqn{\boldsymbol{p}\ge\mathbf{0}}. This is a **repair /
#' renormalisation** step for estimated proportions, not a Euclidean
#' projection onto the simplex and not a statistical-identifiability
#' constraint.
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param tolerance Non-negative tolerance for treating entries as zero
#'   (default `100 * .Machine$double.eps`).
#'
#' @return Numeric vector on the simplex \eqn{\Delta^{J-1}}.
#'
#' @seealso [compositions::clo()] for compositional closure.
#'
#' @srrstats {G3.0} Default `tolerance` is `100 * .Machine$double.eps`;
#'   callers must not test floating-point equality with `==`.
#'
#' @examples
#' repair_simplex(c(0.2, 0.3, 0.5 + 1e-12))
#' @export
repair_simplex <- function(p, tolerance = 100 * .Machine$double.eps) {
  if (!is.numeric(p) || length(p) == 0L || anyNA(p)) {
    stop("`p` must be a non-empty numeric vector without missing values.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0) {
    stop("`tolerance` must be a single non-negative number.")
  }

  p[abs(p) <= tolerance] <- 0

  if (any(p < 0)) {
    stop("`p` contains negative values beyond numerical tolerance.")
  }

  total <- sum(p)
  if (!is.finite(total) || total <= 0) {
    stop("`p` must have a positive, finite sum.")
  }

  p <- p / total
  p[abs(p - 1) <= tolerance] <- 1
  p
}

#' Normalised Shannon entropy of a discrete distribution
#'
#' @description
#' For a probability vector \eqn{\boldsymbol{p}\in\Delta^{J-1}} over
#' \eqn{J} classes (here: cell types), the Shannon entropy is
#' \deqn{
#'   H(\boldsymbol{p})
#'   =
#'   -\sum_{j=1}^{J} p_j \log p_j,
#' }
#' with the convention \eqn{0\log 0 = 0}. Dividing by the maximum
#' entropy \eqn{\log J} (uniform over all \eqn{J} classes) yields
#' **Pielou's evenness**
#' \deqn{
#'   H^{\star}(\boldsymbol{p})
#'   =
#'   \frac{H(\boldsymbol{p})}{\log J}
#'   \in[0,1],
#' }
#' so \eqn{H^{\star}=0} for a Dirac mass on one type and
#' \eqn{H^{\star}=1} for the uniform distribution over the \eqn{J}
#' cell types. Zero masses are dropped only inside the sum; the
#' normaliser still uses the original class count \eqn{J}.
#'
#' This is preferable to changing the logarithm base to the number of
#' *positive* masses \eqn{J'}, which would renormalise only on the
#' support and hide sparsity relative to the full panel of cell types.
#'
#' @section Illustration:
#' Panel B of the figure contrasts a peaked (type-specific, low entropy)
#' share vector with a flat (multi-type, high entropy) one for \eqn{J=5}
#' classes.
#'
#' \if{html}{\figure{gini_vs_entropy_specificity.png}{options: width=700
#'   alt="Gini versus Shannon entropy"}}
#' \if{latex}{\figure{gini_vs_entropy_specificity.png}{options: width=5.5in}}
#'
#' @param ratios Numeric vector of cellular proportions
#'   \eqn{\boldsymbol{p}\in\Delta^{J-1}} (length \eqn{J}).
#' @return Scalar normalised entropy \eqn{H^{\star}\in[0,1]}.
#' @export
#' @examples
#' compute_shannon_entropy(c(1, 0, 0)) # Dirac -> 0
#' compute_shannon_entropy(rep(1 / 3, 3)) # uniform -> 1
compute_shannon_entropy <- function(ratios) {
  if (!is.numeric(ratios) || length(ratios) == 0L || anyNA(ratios)) {
    stop("`ratios` must be a non-empty numeric vector without missing values.")
  }
  if (min(ratios) < 0 || max(ratios) > 1) {
    stop("Probabilities must be strictly included between 0 and 1.")
  }

  # J = number of cell types in the full panel (before dropping zeros)
  J <- length(ratios)
  if (J == 1L) {
    return(0)
  }

  ratios <- ratios[ratios > 0]
  ratios <- ratios / sum(ratios)
  # Pielou evenness: H / log(J)
  -sum(ratios * log(ratios)) / log(J)
}

#' Average pairwise overlap of a Gaussian mixture
#'
#' @description
#' Returns MixSim's **BarOmega**: the unweighted mean of pairwise
#' overlaps
#' \deqn{
#'   \overline{\omega}
#'   =
#'   \frac{2}{J(J-1)}
#'   \sum_{1\le j<\ell\le J}
#'   \bigl(\Omega_{j\ell}+\Omega_{\ell j}\bigr)
#'   \in[0,1]
#' }
#' (up to the MixSim numerical convention), where
#' \eqn{\Omega_{j\ell}=\Pr_{X\sim f_j}(X\text{ classified as }\ell)}
#' already uses the mixture weights \eqn{\boldsymbol{p}} inside the
#' Bayes / MAP rule of [MixSim::overlap()]. Do **not** multiply the
#' directional masses by \eqn{p_j} again.
#'
#' @param true_theta List validated by [check_true_theta()]: `p` (length
#'   \eqn{J} or \eqn{J\times N}), `mu` (\eqn{G\times J}), `sigma`
#'   (\eqn{G\times G\times J}).
#' @param J Number of cell types (components). Defaults to the third
#'   dimension of `sigma`.
#' @return Scalar average pairwise overlap (MixSim `BarOmega`).
#' @export
#' @seealso [check_true_theta()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_overlap(theta)
compute_average_overlap <- function(true_theta, J = NULL) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = TRUE,
    J = J,
    second_moment = "sigma"
  )
  # MixSim::overlap expects Mu as J x G
  MixSim::overlap(
    Pi = theta$p,
    Mu = t(theta$mu),
    S = theta$sigma
  )$BarOmega
}

#' Gene scores from multinomial elastic-net cell-type classification
#'
#' @description
#' Fits a multinomial (or binomial) elastic net ([glmnet::glmnet()]) that
#' predicts cell type from expression features. Inputs are purified
#' expression profiles
#' \eqn{\boldsymbol{X}\in\mathcal{M}_{G\times J\times N}} (genes \eqn{\times}
#' cell types \eqn{\times} samples) and length-\eqn{J} cell-type labels.
#' Variability across samples replaces synthetic isotropic noise. Gene scores
#' are the sum over classes of absolute coefficients at a chosen
#' \eqn{\lambda} (intercept excluded). For nested / CV selection of
#' \eqn{\lambda}, see the experimental
#' \code{compute_glmnet_gene_scores_cv()} helper (not shipped in the package
#' build).
#'
#' @param expression_profiles Numeric array
#'   \eqn{G\times J\times N} of purified profiles.
#' @param celltype_labels Character or factor labels of length \eqn{J}
#'   (one per cell-type slice).
#' @param alpha Elastic-net mixing parameter in \eqn{[0,1]} (default `0.5`).
#' @param lambda Optional penalty value at which coefficients are extracted.
#'   When `NULL`, uses the smallest \eqn{\lambda} on the fitted path.
#' @param ... Additional arguments forwarded to [glmnet::glmnet()].
#'
#' @return Named numeric vector of length \eqn{G} (gene scores; larger means
#'   stronger multinomial signal).
#'
#' @export
#' @seealso [compute_average_jeffreys()], [compute_average_overlap()]
#' @examples
#' set.seed(1)
#' G <- 4L
#' J <- 3L
#' N <- 12L
#' profiles <- array(0, dim = c(G, J, N))
#' for (j in seq_len(J)) {
#'   profiles[j, j, ] <- 5 + stats::rnorm(N, sd = 0.2)
#' }
#' labels <- paste0("ct", seq_len(J))
#' scores <- compute_glmnet_gene_scores(profiles, labels)
#' names(scores)[which.max(scores)]
compute_glmnet_gene_scores <- function(
  expression_profiles,
  celltype_labels,
  alpha = 0.5,
  lambda = NULL,
  ...
) {
  if (
    is.null(dim(expression_profiles)) ||
      length(dim(expression_profiles)) != 3L
  ) {
    stop("`expression_profiles` must be a G x J x N array.")
  }
  if (!is.numeric(expression_profiles) || anyNA(expression_profiles)) {
    stop("`expression_profiles` must be numeric without missing values.")
  }
  G <- dim(expression_profiles)[[1L]]
  J <- dim(expression_profiles)[[2L]]
  N <- dim(expression_profiles)[[3L]]
  if (G < 1L || J < 2L || N < 1L) {
    stop("`expression_profiles` must be G x J x N with G >= 1, J >= 2, N >= 1.")
  }
  if (length(celltype_labels) != J) {
    stop("`celltype_labels` must have length J (second dim of profiles).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number in [0, 1].")
  }

  gene_names <- dimnames(expression_profiles)[[1L]]
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(G))
  }
  cell_names <- as.character(celltype_labels)

  # Design: one row per (cell type, sample) -> (J * N) x G
  x <- matrix(NA_real_, nrow = J * N, ncol = G)
  y <- factor(rep(cell_names, times = N), levels = unique(cell_names))
  row_id <- 1L
  for (n in seq_len(N)) {
    for (j in seq_len(J)) {
      x[row_id, ] <- expression_profiles[, j, n]
      row_id <- row_id + 1L
    }
  }
  colnames(x) <- gene_names

  family <- if (J == 2L) "binomial" else "multinomial"
  glm_args <- list(
    x = x,
    y = y,
    family = family,
    alpha = alpha
  )
  if (identical(family, "multinomial")) {
    glm_args$type.multinomial <- "grouped"
  }
  dots <- list(...)
  for (nm in names(dots)) {
    glm_args[[nm]] <- dots[[nm]]
  }

  fit <- do.call(glmnet::glmnet, glm_args)
  s_use <- if (is.null(lambda)) min(fit$lambda) else lambda
  beta <- stats::coef(fit, s = s_use)

  scores <- numeric(G)
  names(scores) <- gene_names
  if (identical(family, "binomial")) {
    cj <- as.matrix(beta)[-1L, 1L]
    scores <- abs(as.numeric(cj))
    names(scores) <- gene_names
  } else {
    for (j in seq_along(beta)) {
      cj <- as.matrix(beta[[j]])[-1L, 1L]
      scores <- scores + abs(as.numeric(cj))
    }
  }
  scores
}

#' Pairwise Jeffreys divergence between two multivariate Gaussians
#'
#' @description
#' Symmetrised KL divergence
#' \eqn{J(f_{j},f_{\ell})=D_{\mathrm{KL}}(f_{j}\parallel f_{\ell})+
#' D_{\mathrm{KL}}(f_{\ell}\parallel f_{j})} for
#' \eqn{f_{j}=\mathcal{N}_{G}(
#'   \boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_{j}
#' )}
#' and likewise for \eqn{\ell}, using the closed form in the feature-selection
#' vignette.
#'
#' @param mu_j,mu_l Numeric mean vectors of length \eqn{G}.
#' @param sigma_j,sigma_l Numeric \eqn{G\times G} covariances.
#'
#' @return Non-negative scalar Jeffreys divergence.
#'
#' @references
#' Kullback S, Leibler RA (1951).
#' "On Information and Sufficiency."
#' \emph{The Annals of Mathematical Statistics} 22(1), 79--86.
#' \doi{10.1214/aoms/1177729694}.
#'
#' Multivariate normal KL closed form:
#' \url{https://statproofbook.github.io/P/mvn-kl.html}.
#'
#' Symmetrised (Jeffreys) divergence:
#' \url{https://en.wikipedia.org/wiki/Kullback-Leibler_divergence}.
#'
#' @keywords internal
#' @examples
#' .jeffreys_gaussian(c(0, 0), c(1, 0), diag(2), diag(2))
#' @export
.jeffreys_gaussian <- function(mu_j, mu_l, sigma_j, sigma_l) {
  mu_j <- as.numeric(mu_j)
  mu_l <- as.numeric(mu_l)
  sigma_j <- as.matrix(sigma_j)
  sigma_l <- as.matrix(sigma_l)
  g <- length(mu_j)
  if (length(mu_l) != g) {
    stop("`mu_j` and `mu_l` must have the same length.", call. = FALSE)
  }
  if (!identical(dim(sigma_j), c(g, g)) || !identical(dim(sigma_l), c(g, g))) {
    stop("Covariances must be G x G.", call. = FALSE)
  }

  theta_j <- solve(sigma_j)
  theta_l <- solve(sigma_l)
  cov_term <- 0.5 *
    sum(diag(theta_l %*% sigma_j + theta_j %*% sigma_l - 2 * diag(g)))
  delta <- mu_j - mu_l
  mean_term <- 0.5 * as.numeric(crossprod(delta, (theta_j + theta_l) %*% delta))
  cov_term + mean_term
}

#' Average pairwise Jeffreys divergence of a Gaussian mixture
#'
#' @description
#' Returns the \eqn{p_j p_{\ell}}-weighted mean of pairwise Jeffreys
#' (symmetrised KL) divergences between purified Gaussians
#' \deqn{
#'   \overline{J}
#'   =
#'   \frac{
#'     \sum_{1\le j<\ell\le J} p_j p_{\ell}\, J(f_{j},f_{\ell})
#'   }{
#'     \sum_{1\le j<\ell\le J} p_j p_{\ell}
#'   }
#'   \in[0,\infty),
#' }
#' with
#' \eqn{f_{j}=\mathcal{N}_{G}(
#'   \boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_{j}
#' )}.
#' If `p` is omitted it defaults to the equi-balanced vector \eqn{1/J}, which
#' recovers the uniform pairwise average
#' \eqn{2/(J(J-1))\sum_{j<\ell}J(f_{j},f_{\ell})}.
#'
#' @param true_theta List validated by [check_true_theta()]: `mu`
#'   (\eqn{G\times J}), `sigma` (\eqn{G\times G\times J}), and optionally `p`
#'   (length \eqn{J} or \eqn{J\times N}). If `p` is missing it is set to
#'   \eqn{(1/J,\ldots,1/J)}.
#' @param J Number of cell types. Defaults to the third dimension of `sigma`.
#'
#' @return Scalar average pairwise Jeffreys divergence.
#' @export
#' @seealso [check_true_theta()], [compute_average_overlap()],
#'   [compute_glmnet_gene_scores()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_jeffreys(theta)
compute_average_jeffreys <- function(
  true_theta,
  J = NULL
) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = FALSE,
    J = J,
    second_moment = "sigma"
  )
  p <- theta$p
  if (is.null(p)) {
    p <- rep(1 / theta$J, theta$J)
  }
  mu_gj <- theta$mu
  sigma <- theta$sigma
  J <- theta$J

  pair_sum <- 0
  weight_sum <- 0
  for (j in seq_len(J - 1L)) {
    for (ell in seq.int(j + 1L, J)) {
      w <- p[[j]] * p[[ell]]
      pair_sum <- pair_sum +
        w *
          .jeffreys_gaussian(
            mu_gj[, j],
            mu_gj[, ell],
            sigma[,, j],
            sigma[,, ell]
          )
      weight_sum <- weight_sum + w
    }
  }
  if (weight_sum <= 0) {
    stop("Pair weights from `p` are degenerate.")
  }
  pair_sum / weight_sum
}
