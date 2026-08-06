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
#' \if{html}{\figure{gini_vs_entropy_specificity.png}{options: width="100\%"
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
#' @param true_theta List with `p` (length \eqn{J}), `mu` (\eqn{G\times J}
#'   mean matrix) and `sigma` (\eqn{G\times G\times J} covariance array),
#'   as in a GMM
#'   \eqn{(\boldsymbol{p},\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})}.
#' @param J Number of cell types (components). Defaults to
#'   `length(true_theta$p)`.
#' @return Scalar average pairwise overlap (MixSim `BarOmega`).
#' @export
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_overlap(theta)
compute_average_overlap <- function(true_theta, J = length(true_theta$p)) {
  stopifnot(is.list(true_theta))
  p <- true_theta$p
  mu <- as.matrix(true_theta$mu)
  sigma <- true_theta$sigma

  if (length(p) != J) {
    stop("`J` must equal length(true_theta$p).")
  }
  if (is.null(dim(sigma)) || length(dim(sigma)) != 3L) {
    stop("`true_theta$sigma` must be a G x G x J array.")
  }

  G <- dim(sigma)[[1L]]
  n_celltypes <- dim(sigma)[[3L]]
  if (dim(sigma)[[2L]] != G || n_celltypes != length(p)) {
    stop("`sigma` dims must be G x G x J with J = length(p).")
  }

  # DeCovarT stores mu as G x J; MixSim::overlap expects Mu as J x G
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_mixsim <- t(mu)
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_mixsim <- mu
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  MixSim::overlap(Pi = p, Mu = mu_mixsim, S = sigma)$BarOmega
}

#' Gene scores from multinomial elastic-net cell-type classification
#'
#' @description
#' Fits a multinomial elastic net ([glmnet::cv.glmnet()]) that predicts cell
#' type from expression features, using only the mean signature
#' \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (no covariances). Each
#' column \eqn{\boldsymbol{\mu}_{\cdot j}} is expanded into `n_rep` isotropic
#' Gaussian perturbations so that cross-validation is well-defined with one
#' prototype per type. Gene scores are the sum over classes of absolute
#' coefficients at `lambda.min` (intercept excluded).
#'
#' This is the `glmnet` screen in the four-score shortlist of the feature-selection
#' vignette (paired with Jeffreys, MixSim overlap and DEGs).
#'
#' @param mu Numeric mean signature \eqn{\boldsymbol{\mu}} with genes in rows
#'   and cell types in columns (\eqn{G\times J}).
#' @param alpha Elastic-net mixing parameter in \eqn{[0,1]} (default `0.5`).
#' @param n_rep Number of noisy replicates per cell-type mean (default `20L`).
#' @param noise_sd Isotropic Gaussian noise sd for replicates. Default is
#'   `1e-3` times the mean absolute column scale of `mu`.
#' @param nfolds Number of CV folds (default `min(10L, n_rep)`).
#' @param ... Additional arguments forwarded to [glmnet::cv.glmnet()].
#'
#' @return Named numeric vector of length \eqn{G} (gene scores; larger means
#'   stronger multinomial signal).
#'
#' @export
#' @seealso [compute_average_jeffreys()], [compute_average_overlap()]
#' @examples
#' set.seed(1)
#' mu <- cbind(c(0, 0, 5), c(5, 0, 0), c(0, 5, 0))
#' rownames(mu) <- paste0("g", seq_len(nrow(mu)))
#' scores <- compute_glmnet_gene_scores(mu)
#' names(scores)[which.max(scores)]
compute_glmnet_gene_scores <- function(
  mu,
  alpha = 0.5,
  n_rep = 20L,
  noise_sd = NULL,
  nfolds = NULL,
  ...
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for compute_glmnet_gene_scores().")
  }
  mu <- as.matrix(mu)
  if (!is.numeric(mu) || anyNA(mu)) {
    stop("`mu` must be a numeric matrix without missing values.")
  }
  G <- nrow(mu)
  J <- ncol(mu)
  if (is.null(G) || is.null(J) || G < 1L || J < 2L) {
    stop("`mu` must be G x J with G >= 1 and J >= 2.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number in [0, 1].")
  }
  n_rep <- as.integer(n_rep)
  if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L) {
    stop("`n_rep` must be a positive integer.")
  }

  gene_names <- rownames(mu)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(G))
  }
  cell_names <- colnames(mu)
  if (is.null(cell_names)) {
    cell_names <- paste0("celltype_", seq_len(J))
  }

  if (is.null(noise_sd)) {
    scale_mu <- mean(abs(mu))
    noise_sd <- if (scale_mu > 0) 1e-3 * scale_mu else 1e-3
  }
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L || noise_sd < 0) {
    stop("`noise_sd` must be a single non-negative number.")
  }

  # Design: n_rep replicates of each mean column -> (J * n_rep) x G
  x <- matrix(NA_real_, nrow = J * n_rep, ncol = G)
  y <- factor(rep(cell_names, each = n_rep), levels = cell_names)
  for (j in seq_len(J)) {
    idx <- ((j - 1L) * n_rep + 1L):(j * n_rep)
    x[idx, ] <- matrix(mu[, j], nrow = n_rep, ncol = G, byrow = TRUE) +
      matrix(stats::rnorm(n_rep * G, sd = noise_sd), nrow = n_rep, ncol = G)
  }
  colnames(x) <- gene_names

  if (is.null(nfolds)) {
    nfolds <- min(10L, n_rep)
  }
  nfolds <- as.integer(nfolds)
  if (nfolds < 2L || nfolds > nrow(x)) {
    stop("`nfolds` must satisfy 2 <= nfolds <= n_rep * J.")
  }

  family <- if (J == 2L) "binomial" else "multinomial"
  cv_args <- list(
    x = x,
    y = y,
    family = family,
    alpha = alpha,
    nfolds = nfolds
  )
  if (identical(family, "multinomial")) {
    cv_args$type.multinomial <- "grouped"
  }
  dots <- list(...)
  for (nm in names(dots)) {
    cv_args[[nm]] <- dots[[nm]]
  }

  fit <- do.call(glmnet::cv.glmnet, cv_args)
  beta <- stats::coef(fit, s = "lambda.min")

  scores <- numeric(G)
  names(scores) <- gene_names
  if (identical(family, "binomial")) {
    # single sparse matrix: intercept + G coefficients
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
#' \eqn{f_{j}=\mathcal{N}_{G}(\boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_{j})}
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
#' \url{https://en.wikipedia.org/wiki/Kullback\%E2\%80\%93Leibler_divergence#Symmetrised_divergence}.
#'
#' @keywords internal
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
#' \eqn{f_{j}=\mathcal{N}_{G}(\boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_{j})}.
#' If `p` is omitted it defaults to the equi-balanced vector \eqn{1/J}, which
#' recovers the uniform pairwise average
#' \eqn{2/(J(J-1))\sum_{j<\ell}J(f_{j},f_{\ell})}.
#'
#' @param true_theta List with `mu` (\eqn{G\times J}), `sigma`
#'   (\eqn{G\times G\times J}), and optionally `p` (length \eqn{J}). If `p` is
#'   missing it is set to \eqn{(1/J,\ldots,1/J)}.
#' @param J Number of cell types. Defaults to `length(p)` after the
#'   equi-balanced default is applied.
#'
#' @return Scalar average pairwise Jeffreys divergence.
#' @export
#' @seealso [compute_average_overlap()], [compute_glmnet_gene_scores()]
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
  stopifnot(is.list(true_theta))
  if (is.null(true_theta$mu) || is.null(true_theta$sigma)) {
    stop("`true_theta` must contain `mu` and `sigma`.")
  }
  mu <- as.matrix(true_theta$mu)
  sigma <- true_theta$sigma

  if (is.null(dim(sigma)) || length(dim(sigma)) != 3L) {
    stop("`true_theta$sigma` must be a G x G x J array.")
  }
  G <- dim(sigma)[[1L]]
  n_celltypes <- dim(sigma)[[3L]]
  if (dim(sigma)[[2L]] != G) {
    stop("`sigma` dims must be G x G x J.")
  }

  # Accept G x J (preferred) or J x G
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_gj <- mu
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_gj <- t(mu)
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  if (is.null(true_theta$p)) {
    p <- rep(1 / n_celltypes, n_celltypes)
  } else {
    p <- true_theta$p
  }
  if (is.null(J)) {
    J <- length(p)
  }
  if (length(p) != J || n_celltypes != J) {
    stop("`J` must match length(p) and the third dimension of sigma.")
  }
  if (J < 2L) {
    stop("`J` must be at least 2.")
  }
  if (any(p < 0) || abs(sum(p) - 1) > 1e-8) {
    stop("`true_theta$p` must be non-negative and sum to 1.")
  }

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
