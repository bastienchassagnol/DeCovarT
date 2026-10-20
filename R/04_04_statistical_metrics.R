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
#' When `open = TRUE`, entries on the closed-simplex boundary are nudged
#' into the relative interior so that [isometric_log_ratio()] is defined
#' (ILR / Marquardt / Newton starts). That open-simplex path is the
#' former role of [starting_simplex()].
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param tolerance Non-negative tolerance for treating entries as zero
#'   (default `100 * .Machine$double.eps`).
#' @param open Logical; if `TRUE`, push the repaired vector off the
#'   boundary into the open simplex.
#' @param nms Optional names attached to the returned vector (e.g.
#'   cell-type colnames).
#'
#' @return Numeric vector on the simplex \eqn{\Delta^{J-1}} (open when
#'   `open = TRUE`).
#'
#' @seealso [compositions::clo()] for compositional closure.
#'
#' @srrstats {G3.0} Default `tolerance` is `100 * .Machine$double.eps`;
#'   callers must not test floating-point equality with `==`.
#'
#' @examples
#' repair_simplex(c(0.2, 0.3, 0.5 + 1e-12))
#' repair_simplex(c(1, 0, 0), open = TRUE)
#' @export
repair_simplex <- function(
  p,
  tolerance = 100 * .Machine$double.eps,
  open = FALSE,
  nms = NULL
) {
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

  if (isTRUE(open)) {
    floor <- tolerance
    if (any(p < floor) || any(p > 1 - floor)) {
      p <- pmax(p, floor)
      p <- p / sum(p)
    }
  }

  if (!is.null(nms)) {
    names(p) <- nms
  }
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
#' normalised Shannon entropy
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
#' @seealso [composition_from_entropy()]
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
  # normalised Shannon entropy: H / log(J)
  -sum(ratios * log(ratios)) / log(J)
}

#' One-dominant composition with a target normalised Shannon entropy
#'
#' Returns a length-\eqn{J} simplex vector of the form
#' \eqn{(1-(J-1)q,\,q,\ldots,q)} whose normalised Shannon entropy
#' [compute_shannon_entropy()] equals `h_star`. The uniform composition
#' (`h_star = 1`) and a Dirac mass (`h_star = 0`) are returned exactly.
#'
#' @param h_star Target \eqn{H^{\star}\in[0,1]}.
#' @param n_celltypes Integer \eqn{J\ge 2}.
#' @param nms Optional names for the returned vector.
#'
#' @return Numeric simplex vector of length \eqn{J}.
#'
#' @examples
#' p <- composition_from_entropy(0.5, n_celltypes = 3L)
#' compute_shannon_entropy(p)
#' @export
#' @seealso [compute_shannon_entropy()]
composition_from_entropy <- function(
  h_star,
  n_celltypes = 3L,
  nms = NULL
) {
  j <- as.integer(n_celltypes)
  h_star <- as.numeric(h_star)
  if (length(j) != 1L || is.na(j) || j < 2L) {
    stop("`n_celltypes` must be an integer of at least 2.", call. = FALSE)
  }
  if (length(h_star) != 1L || is.na(h_star) || h_star < 0 || h_star > 1) {
    stop("`h_star` must lie in [0, 1].", call. = FALSE)
  }

  if (h_star >= 1 - 1e-12) {
    p <- rep(1 / j, j)
    return(repair_simplex(p, nms = nms))
  }
  if (h_star <= 1e-12) {
    p <- c(1, rep(0, j - 1L))
    return(repair_simplex(p, nms = nms))
  }

  h_of_q <- function(q) {
    p1 <- 1 - (j - 1) * q
    compute_shannon_entropy(c(p1, rep(q, j - 1L))) - h_star
  }
  q_max <- (1 / j) * (1 - 1e-8)
  q <- stats::uniroot(
    h_of_q,
    interval = c(1e-12, q_max),
    tol = 1e-12,
    maxiter = 200L
  )$root
  p <- c(1 - (j - 1) * q, rep(q, j - 1L))
  repair_simplex(p, nms = nms)
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
#' For \eqn{G\ge 4} MixSim's Davies quadrature becomes expensive.
#' The helper then switches to [overlap_gaussian_mc()]: stratified
#' Sobol draws, inverse-transform sampling through precomputed
#' Cholesky factors, and log-density MAP comparisons (`n_mc`
#' draws per component; default 10,000).
#'
#' @param true_theta List validated by [check_true_theta()]: `p` (length
#'   \eqn{J} or \eqn{J\times N}), `mu` (\eqn{G\times J}), `sigma`
#'   (\eqn{G\times G\times J}).
#' @param J Number of cell types (components). Defaults to the third
#'   dimension of `sigma`.
#' @param n_mc Monte Carlo draws per component when \eqn{G\ge 4}.
#' @param seed Optional seed forwarded to [overlap_gaussian_mc()].
#' @param verbose If `TRUE` (default), announce the MixSim-to-MC switch
#'   with `cli` when it is installed.
#' @return Scalar average pairwise overlap (MixSim `BarOmega`).
#' @export
#' @seealso [check_true_theta()], [overlap_gaussian_mc()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_overlap(theta)
compute_average_overlap <- function(
  true_theta,
  J = NULL,
  n_mc = 10000L,
  seed = NULL,
  verbose = TRUE
) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = TRUE,
    J = J,
    second_moment = "sigma"
  )
  if (theta$G >= 4L) {
    if (isTRUE(verbose)) {
      .ui_info(c(
        "G = {.val {theta$G}} >= 4: MixSim Davies overlap is replaced",
        "by stratified Sobol Monte Carlo ({.fn overlap_gaussian_mc},",
        "n_mc = {.val {n_mc}} draws per component)."
      ))
    }
    return(
      overlap_gaussian_mc(
        true_theta = theta,
        n_mc = n_mc,
        seed = seed,
        J = theta$J
      )$BarOmega
    )
  }
  .check_suggested_package("MixSim", "compute_average_overlap")
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
#' \dontrun{
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
#' }
compute_glmnet_gene_scores <- function(
  expression_profiles,
  celltype_labels,
  alpha = 0.5,
  lambda = NULL,
  ...
) {
  .check_suggested_package("glmnet", "compute_glmnet_gene_scores")
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
#' Sentinel for a score that does not apply in this metrics block
#'
#' Used instead of `NULL` / `NA` so tables stay joinable.
#'
#' @keywords internal
#' @noRd
.irrelevant_metric <- function() {
  "irrelevant"
}

#' Bind metric tibbles without collapsing to `NULL`
#'
#' @keywords internal
#' @noRd
.bind_metrics_rows <- function(xs) {
  out <- dplyr::bind_rows(xs)
  if (is.null(out)) {
    tibble::tibble(status = .irrelevant_metric())
  } else {
    out
  }
}

#' Mean squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar MSE.
#' @keywords internal
#' @noRd
.mse <- function(actual, predicted) {
  mean((actual - predicted)^2)
}

#' Root mean squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar RMSE.
#' @keywords internal
#' @noRd
.rmse <- function(actual, predicted) {
  sqrt(.mse(actual, predicted))
}

#' Mean absolute error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar MAE.
#' @keywords internal
#' @noRd
.mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}

#' Euclidean projection onto the probability simplex
#'
#' @param v Numeric vector.
#' @return Non-negative vector summing to 1.
#' @keywords internal
#' @noRd
.project_simplex <- function(v) {
  v <- as.numeric(v)
  n <- length(v)
  if (n == 0L) {
    stop("`v` must be non-empty.", call. = FALSE)
  }
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u > (cssv - 1) / seq_len(n)))
  theta <- (cssv[[rho]] - 1) / rho
  pmax(v - theta, 0)
}

#' Simplex KKT / projected-score residual
#'
#' \eqn{\|\Pi_\Delta(p + \alpha g) - p\|_2 / \alpha} for maximisation.
#'
#' @keywords internal
#' @noRd
.kkt_residual <- function(p, grad, alpha = 1) {
  p <- as.numeric(p)
  grad <- as.numeric(grad)
  if (length(p) != length(grad) || length(p) == 0L) {
    return(NA_real_)
  }
  if (!is.finite(alpha) || alpha <= 0) {
    alpha <- 1
  }
  step <- .project_simplex(p + alpha * grad) - p
  sqrt(sum(step^2)) / alpha
}

#' Total variation / normalised MAE on the simplex
#'
#' @keywords internal
#' @noRd
.tv <- function(p, p_hat) {
  0.5 * sum(abs(as.numeric(p) - as.numeric(p_hat)))
}

#' Aitchison (clr Euclidean) distance on the simplex
#'
#' @keywords internal
#' @noRd
.aitchison_distance <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  if (length(p) != length(p_hat) || length(p) < 2L) {
    return(NA_real_)
  }
  if (any(!is.finite(p)) || any(!is.finite(p_hat))) {
    return(NA_real_)
  }
  eps <- .Machine$double.eps
  p <- pmax(p, eps)
  p_hat <- pmax(p_hat, eps)
  p <- p / sum(p)
  p_hat <- p_hat / sum(p_hat)
  clr <- function(x) {
    log(x) - mean(log(x))
  }
  sqrt(sum((clr(p) - clr(p_hat))^2))
}

#' Vectorised Aitchison distance for two-part compositions
#'
#' @keywords internal
#' @noRd
.aitchison_pair <- function(p1, p2, h1, h2) {
  eps <- .Machine$double.eps
  p1 <- pmax(as.numeric(p1), eps)
  p2 <- pmax(as.numeric(p2), eps)
  h1 <- pmax(as.numeric(h1), eps)
  h2 <- pmax(as.numeric(h2), eps)
  sp <- p1 + p2
  sh <- h1 + h2
  p1 <- p1 / sp
  p2 <- p2 / sp
  h1 <- h1 / sh
  h2 <- h2 / sh
  mp <- 0.5 * (log(p1) + log(p2))
  mh <- 0.5 * (log(h1) + log(h2))
  sqrt((log(p1) - mp - (log(h1) - mh))^2 + (log(p2) - mp - (log(h2) - mh))^2)
}

#' \eqn{L_\infty} absolute error
#'
#' @keywords internal
#' @noRd
.maxae <- function(p, p_hat) {
  max(abs(as.numeric(p) - as.numeric(p_hat)))
}

#' Normalised angular distance on \eqn{[0, 1]}
#'
#' @keywords internal
#' @noRd
.angular_distance <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  denom <- sqrt(sum(p^2)) * sqrt(sum(p_hat^2))
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  cos_angle <- sum(p * p_hat) / denom
  cos_angle <- min(1, max(-1, cos_angle))
  (2 / pi) * acos(cos_angle)
}

#' HADACA-style SDID (\eqn{\sqrt{\sin\theta}})
#'
#' @keywords internal
#' @noRd
.sdid <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  denom <- sqrt(sum(p^2)) * sqrt(sum(p_hat^2))
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  cos_angle <- sum(p * p_hat) / denom
  sin_angle <- sqrt(1 - min(1, cos_angle^2))
  sqrt(sin_angle)
}

#' Safe Pearson correlation (NA if either vector is constant)
#'
#' @keywords internal
#' @noRd
.pearson_safe <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  sx <- stats::sd(x)
  sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) {
    return(NA_real_)
  }
  stats::cor(x, y, method = "pearson")
}

#' Presence / spillover diagnostics for one composition pair
#'
#' @keywords internal
#' @noRd
.presence_counts <- function(p, p_hat, threshold = 1e-4) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  # Failed fits store NA in p_hat; treat them as absent (mass 0).
  p[!is.finite(p)] <- 0
  p_hat[!is.finite(p_hat)] <- 0
  present_true <- p > threshold
  present_hat <- p_hat > threshold
  list(
    tp = sum(present_true & present_hat),
    fp = sum(!present_true & present_hat),
    fn = sum(present_true & !present_hat),
    tn = sum(!present_true & !present_hat),
    false_positive_mass = sum(p_hat[!present_true])
  )
}

#' F1 from presence counts
#'
#' @keywords internal
#' @noRd
.f1_from_counts <- function(tp, fp, fn) {
  tp <- as.numeric(tp)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  if (!is.finite(tp)) {
    tp <- 0
  }
  if (!is.finite(fp)) {
    fp <- 0
  }
  if (!is.finite(fn)) {
    fn <- 0
  }
  denom <- tp + fp + fn
  if (denom == 0) {
    return(NA_real_)
  }
  precision <- if ((tp + fp) == 0) {
    NA_real_
  } else {
    tp / (tp + fp)
  }
  recall <- if ((tp + fn) == 0) {
    NA_real_
  } else {
    tp / (tp + fn)
  }
  if (
    !is.finite(precision) || !is.finite(recall) || (precision + recall) == 0
  ) {
    return(0)
  }
  2 * precision * recall / (precision + recall)
}

#' Process memory in bytes (PSS on Linux when `ps` is available)
#'
#' Windows `ps` objects have no `pss` field; `x[["pss"]]` on a named
#' numeric vector then errors with subscript-out-of-bounds.
#'
#' @keywords internal
#' @noRd
.process_memory_bytes <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(NA_real_)
  }
  info <- tryCatch(
    ps::ps_memory_full_info(),
    error = function(e) {
      tryCatch(ps::ps_memory_info(), error = function(e2) NULL)
    }
  )
  if (is.null(info)) {
    return(NA_real_)
  }
  .one_named <- function(nm) {
    if (!nm %in% names(info)) {
      return(NULL)
    }
    val <- info[[nm]]
    if (!is.null(val) && length(val) == 1L && is.finite(val)) {
      return(as.numeric(val))
    }
    NULL
  }
  pss <- .one_named("pss")
  if (!is.null(pss)) {
    return(pss)
  }
  rss <- .one_named("rss")
  if (!is.null(rss)) {
    return(rss)
  }
  NA_real_
}

#' Monte Carlo SE of an empirical coverage probability
#'
#' @keywords internal
#' @noRd
.mcse_coverage <- function(covered) {
  coverage_mc_interval(covered, method = "wald")$mcse
}

#' Monte Carlo interval for an empirical coverage probability
#'
#' Interval for a binomial coverage rate
#' \eqn{\hat\pi=X/N} from independent coverage indicators.
#' Wilson is the default because the parameter is bounded in
#' \eqn{[0,1]} and Wald intervals behave poorly near 0 or 1.
#'
#' @param covered Logical or 0/1 coverage indicators (one per replicate).
#' @param conf_level Confidence level for the interval around
#'   \eqn{\hat\pi} (default `0.95`).
#' @param method `"wilson"` (default), `"wald"`, or `"agresti_coull"`.
#'
#' @return A list with `n`, `successes`, `coverage`, `mcse`, `lower`,
#'   `upper`, and `method`.
#'
#' @examples
#' coverage_mc_interval(c(TRUE, TRUE, TRUE, FALSE))
#' @export
coverage_mc_interval <- function(
  covered,
  conf_level = 0.95,
  method = c("wilson", "wald", "agresti_coull")
) {
  method <- .match_arg_case_insensitive(
    method,
    c(
      "wilson",
      "wald",
      "agresti_coull"
    )
  )
  covered <- as.numeric(covered)
  covered <- covered[is.finite(covered)]
  n <- length(covered)
  if (n == 0L) {
    return(list(
      n = 0L,
      successes = 0L,
      coverage = NA_real_,
      mcse = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      method = method
    ))
  }
  x <- sum(covered)
  p_hat <- x / n
  mcse <- sqrt(p_hat * (1 - p_hat) / n)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  bounds <- switch(
    method,
    wald = {
      c(
        max(0, p_hat - z * mcse),
        min(1, p_hat + z * mcse)
      )
    },
    wilson = {
      d <- 1 + z^2 / n
      centre <- (p_hat + z^2 / (2 * n)) / d
      half <- z /
        d *
        sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
      c(max(0, centre - half), min(1, centre + half))
    },
    agresti_coull = {
      n_tilde <- n + z^2
      p_tilde <- (x + z^2 / 2) / n_tilde
      half <- z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
      c(max(0, p_tilde - half), min(1, p_tilde + half))
    }
  )
  list(
    n = n,
    successes = as.integer(x),
    coverage = p_hat,
    mcse = mcse,
    lower = bounds[[1L]],
    upper = bounds[[2L]],
    method = method
  )
}

#' Align a proportion vector with signature colnames
#'
#' @keywords internal
#' @noRd
.as_named_p <- function(p, nms) {
  p <- as.numeric(p)
  if (!is.null(nms) && length(nms) == length(p)) {
    names(p) <- nms
  }
  p
}

#' Coerce a solver return value to a named proportion vector
#'
#' @keywords internal
#' @noRd
.coerce_estimated_p <- function(out, nms) {
  if (is.list(out) && !is.null(out$coefficients)) {
    p <- out$coefficients
  } else {
    p <- out
  }
  if (is.matrix(p) || is.data.frame(p)) {
    p <- as.numeric(p[, 1L, drop = TRUE])
  }
  .as_named_p(p, nms)
}

#' Optional standard errors from a solver return value
#'
#' @keywords internal
#' @noRd
.extract_estimated_se <- function(out, n_celltypes) {
  if (!is.list(out)) {
    return(rep(NA_real_, n_celltypes))
  }
  if (!is.null(out$se)) {
    return(as.numeric(out$se)[seq_len(n_celltypes)])
  }
  if (!is.null(out$vcov)) {
    return(sqrt(pmax(diag(as.matrix(out$vcov)), 0)))
  }
  rep(NA_real_, n_celltypes)
}

#' Map over bulk samples, optionally with `furrr`
#'
#' @keywords internal
#' @noRd
.map_samples <- function(
  n,
  fun,
  cores,
  verbose = FALSE,
  progress_every = 10L,
  progress_label = NULL
) {
  idx <- seq_len(n)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    .ui_abort("{.arg cores} must be a positive integer.")
  }
  if (cores == 1L) {
    if (!isTRUE(verbose) || n < 1L) {
      return(lapply(idx, fun))
    }
    every <- as.integer(progress_every)
    if (length(every) != 1L || is.na(every) || every < 1L) {
      every <- 10L
    }
    prefix <- if (is.null(progress_label) || !nzchar(progress_label)) {
      ""
    } else {
      paste0(progress_label, " \u00b7 ")
    }
    return(lapply(idx, function(i) {
      out <- fun(i)
      if (i == n || i == 1L || (i %% every == 0L)) {
        .ui_info("{prefix}inferred sample {.val {i}}/{.val {n}}.")
      }
      out
    }))
  }
  if (isTRUE(verbose)) {
    .ui_info(
      "Mapping {.val {n}} samples on {.val {cores}} workers (no per-sample ticks)."
    )
  }
  .check_suggested_package("furrr", "deconvolute_ratios")
  .check_suggested_package("future", "deconvolute_ratios")
  old_plan <- future::plan(future::multisession, workers = cores)
  on.exit(future::plan(old_plan), add = TRUE)
  furrr::future_map(
    idx,
    fun,
    .options = furrr::furrr_options(
      # L'Ecuyer-CMRG streams per worker (R >= 4.1 / future.seed).
      seed = TRUE,
      packages = "DeCovarT"
    )
  )
}

#' Closed-form (non-iterative) deconvolution keys
#'
#' @keywords internal
#' @noRd
.closed_form_algorithm <- function(algorithm) {
  .normalise_algorithm_key(algorithm) %in%
    c("nnls", "lsei", "cibersort", "deconrnaseq")
}

#' Empty homogenised optimiser record
#'
#' @keywords internal
#' @noRd
.empty_convergence <- function() {
  list(
    iterations = NA_integer_,
    criterion = NA_character_,
    code = NA_integer_,
    rdm = NA_real_
  )
}

#' First finite scalar of the requested type
#'
#' @keywords internal
#' @noRd
.first_scalar <- function(..., type = "integer") {
  vals <- list(...)
  for (v in vals) {
    if (is.null(v) || length(v) < 1L) {
      next
    }
    x <- v[[1L]]
    if (identical(type, "character")) {
      x <- as.character(x)
      if (is.na(x) || !nzchar(x)) {
        next
      }
      return(x)
    }
    if (identical(type, "double")) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 1L && is.finite(x)) {
        return(x)
      }
      next
    }
    x <- suppressWarnings(as.integer(x))
    if (length(x) == 1L && !is.na(x)) {
      return(x)
    }
  }
  if (identical(type, "character")) {
    return(NA_character_)
  }
  if (identical(type, "double")) {
    return(NA_real_)
  }
  NA_integer_
}

#' Homogenised optimiser diagnostics (iterative solvers only)
#'
#' Closed-form methods (NNLS, LSEI, CIBERSORT) return `NA` fields.
#'
#' @keywords internal
#' @noRd
.homogenise_convergence <- function(raw_out, algorithm) {
  if (.closed_form_algorithm(algorithm)) {
    return(.empty_convergence())
  }
  conv <- NULL
  if (is.list(raw_out) && !is.null(raw_out$convergence)) {
    conv <- raw_out$convergence
  }
  if (!is.list(conv)) {
    return(.empty_convergence())
  }
  list(
    iterations = .first_scalar(conv$iterations, conv$ni, type = "integer"),
    criterion = .first_scalar(
      conv$criterion,
      conv$message,
      type = "character"
    ),
    code = .first_scalar(
      conv$code,
      conv$istop,
      conv$convergence,
      type = "integer"
    ),
    rdm = .first_scalar(conv$rdm, type = "double")
  )
}
