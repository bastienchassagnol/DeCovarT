#' Hoyer sparsity of a nonnegative vector
#'
#' For a nonnegative weight vector \eqn{\boldsymbol{w}} of length \eqn{M}
#' (here: absolute off-diagonal partial correlations or
#' \eqn{\lvert r_{gh}\rvert}),
#' \deqn{
#'   S_{\mathrm{Hoyer}}(\boldsymbol{w})
#'   =
#'   \frac{
#'     \sqrt{M}-\lVert\boldsymbol{w}\rVert_1/\lVert\boldsymbol{w}\rVert_2
#'   }{
#'     \sqrt{M}-1
#'   }
#'   \in[0,1].
#' }
#' Zero is a flat spectrum of magnitudes; one is a single dominant edge.
#'
#' @keywords internal
#' @noRd
.hoyer_sparsity <- function(x) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) {
    return(NA_real_)
  }
  nrm2 <- sqrt(sum(x^2))
  if (nrm2 <= 0) {
    return(NA_real_)
  }
  (sqrt(n) - sum(x) / nrm2) / (sqrt(n) - 1)
}

#' Hellinger distance between two multivariate Gaussians
#'
#' Closed form through the Bhattacharyya coefficient
#' \deqn{
#'   \log\mathrm{BC}
#'   =
#'   \tfrac14\log\lvert\boldsymbol{\Sigma}_j\rvert
#'   +\tfrac14\log\lvert\boldsymbol{\Sigma}_\ell\rvert
#'   -\tfrac12\log\lvert\bar{\boldsymbol{\Sigma}}\rvert
#'   -\tfrac18\boldsymbol{\delta}^{\top}\bar{\boldsymbol{\Sigma}}^{-1}
#'   \boldsymbol{\delta},
#' }
#' with \eqn{\bar{\boldsymbol{\Sigma}}=(\boldsymbol{\Sigma}_j+
#' \boldsymbol{\Sigma}_\ell)/2} and \eqn{\boldsymbol{\delta}=
#' \boldsymbol{\mu}_j-\boldsymbol{\mu}_\ell}; then
#' \eqn{H=\sqrt{1-\mathrm{BC}}\in[0,1]}. It is symmetric, unlike the
#' MixSim misclassification overlap ([compute_average_overlap()]), and
#' bounded, unlike the Jeffreys divergence ([compute_average_jeffreys()]).
#'
#' @param mu_j,mu_l Mean vectors of length \eqn{G}.
#' @param sigma_j,sigma_l \eqn{G\times G} covariance matrices.
#'
#' @return The Hellinger distance, a number in \eqn{[0,1]}.
#'
#' @examples
#' hellinger_gaussian(c(0, 0), c(1, 1), diag(2), diag(2))
#' @export
hellinger_gaussian <- function(mu_j, mu_l, sigma_j, sigma_l) {
  mu_j <- as.numeric(mu_j)
  mu_l <- as.numeric(mu_l)
  sigma_j <- as.matrix(sigma_j)
  sigma_l <- as.matrix(sigma_l)
  sigma_bar <- (sigma_j + sigma_l) / 2
  log_det <- function(s) {
    as.numeric(determinant(s, logarithm = TRUE)$modulus)
  }
  log_bc <- 0.25 *
    log_det(sigma_j) +
    0.25 * log_det(sigma_l) -
    0.5 * log_det(sigma_bar)
  delta <- mu_j - mu_l
  mah <- tryCatch(
    as.numeric(crossprod(delta, solve(sigma_bar, delta))),
    error = function(e) Inf
  )
  log_bc <- log_bc - 0.125 * mah
  bc <- exp(min(log_bc, 0))
  sqrt(max(0, 1 - bc))
}

#' Off-diagonal density and mean degree of a symmetric matrix
#'
#' @keywords internal
#' @noRd
.undirected_density <- function(mat, eps = 1e-8) {
  mat <- as.matrix(mat)
  g <- nrow(mat)
  if (g < 2L) {
    return(list(density = NA_real_, mean_degree = NA_real_, n_edges = NA_real_))
  }
  upper <- abs(mat[upper.tri(mat)]) > eps
  n_edges <- sum(upper)
  n_possible <- g * (g - 1L) / 2
  list(
    density = n_edges / n_possible,
    mean_degree = 2 * n_edges / g,
    n_edges = n_edges
  )
}

#' Directional Fisher information along a simplex contrast
#'
#' For the Gaussian convolution
#' \eqn{\boldsymbol{y}\mid\boldsymbol{p}\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))}
#' with \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}
#' and \eqn{\boldsymbol{\Omega}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}},
#' the ambient Fisher matrix splits as
#' \deqn{
#'   I_{jk}(\boldsymbol{p})
#'   =
#'   \underbrace{\boldsymbol{\mu}_{\cdot j}^{\mathsf{T}}
#'     \boldsymbol{\Omega}\,
#'     \boldsymbol{\mu}_{\cdot k}}_{I^{\mathrm{mean}}_{jk}}
#'   +
#'   \underbrace{2 p_j p_k
#'     \mathrm{tr}\bigl(
#'       \boldsymbol{\Omega}\boldsymbol{\Sigma}_j
#'       \boldsymbol{\Omega}\boldsymbol{\Sigma}_k
#'     \bigr)}_{I^{\mathrm{cov}}_{jk}}.
#' }
#' Along a simplex contrast \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{d}=0},
#' \deqn{
#'   I_{\mathrm{mean}}(\boldsymbol{d})
#'   =
#'   (\boldsymbol{\mu}\boldsymbol{d})^{\mathsf{T}}
#'   \boldsymbol{\Omega}
#'   (\boldsymbol{\mu}\boldsymbol{d}),
#'   \qquad
#'   I_{\mathrm{cov}}(\boldsymbol{d})
#'   =
#'   \tfrac12
#'   \bigl\|
#'     \boldsymbol{\Omega}^{1/2}
#'     (D_{\boldsymbol{d}}\boldsymbol{\Sigma})
#'     \boldsymbol{\Omega}^{1/2}
#'   \bigr\|_{F}^{2}.
#' }
#' The covariance-information fraction is
#' \eqn{f_{\mathrm{cov}}=\mathrm{tr}(I^{\mathrm{cov}})/\mathrm{tr}(I)}.
#' The two matrix summands come from `.information_blocks()`.
#'
#' @keywords internal
#' @noRd
.fisher_directional <- function(
  d,
  p,
  mean_signature_matrix,
  Sigma,
  precision
) {
  d <- as.numeric(d)
  mu_d <- as.numeric(mean_signature_matrix %*% d)
  info_mean <- as.numeric(crossprod(mu_d, precision %*% mu_d))
  n_genes <- nrow(mean_signature_matrix)
  delta_sigma <- matrix(0, n_genes, n_genes)
  for (j in seq_along(p)) {
    delta_sigma <- delta_sigma + 2 * p[[j]] * d[[j]] * Sigma[,, j]
  }
  whitened <- precision %*% delta_sigma
  info_cov <- 0.5 * sum(diag(whitened %*% whitened))
  list(mean = info_mean, cov = info_cov, total = info_mean + info_cov)
}

#' Pairwise cosine summaries of mean-signature columns
#'
#' @keywords internal
#' @noRd
.mean_cosine_summaries <- function(mu) {
  mu <- as.matrix(mu)
  j <- ncol(mu)
  if (j < 2L) {
    return(list(
      mean_abs_cosine = NA_real_,
      min_cosine = NA_real_,
      max_cosine = NA_real_
    ))
  }
  norms <- sqrt(colSums(mu^2))
  gram <- crossprod(mu)
  cos_mat <- gram / tcrossprod(norms)
  off <- cos_mat[upper.tri(cos_mat)]
  list(
    mean_abs_cosine = mean(abs(off)),
    min_cosine = min(off),
    max_cosine = max(off)
  )
}

#' Mean pairwise Euclidean gap of signature columns (CLD)
#'
#' @keywords internal
#' @noRd
.mean_euclidean_distance <- function(mu) {
  obj <- compute_mean_profile_objectives(mu)
  n_pairs <- choose(ncol(as.matrix(mu)), 2L)
  if (!is.finite(n_pairs) || n_pairs < 1L) {
    return(NA_real_)
  }
  obj$sum_euclidean_distance / n_pairs
}

#' Describe a Gaussian-convolution simulation scenario
#'
#' Summarises the generating law
#' \eqn{\boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \sum_j p_j^2\boldsymbol{\Sigma}_j)} at five layers: composition,
#' mean geometry, covariance / SPD diagnostics, network sparsity, and
#' tangent Fisher information (mean versus covariance split). MixSim
#' BarOmega (`mixsim_baromega`), MaxOmega (`mixsim_maxomega`, the largest
#' pairwise overlap) and average pairwise Hellinger of the component
#' Gaussians are kept in `descriptors` (Hellinger is the unweighted pairwise mean;
#' `hellinger_weighted` is optional). Jeffreys / symmetrised KL is
#' returned in `supplementary`. Mean-geometry columns include the
#' highest pairwise cosine and the mean Euclidean (CLD) gap.
#'
#' @param true_theta List with `p`, `mu`, and `sigma` (and optionally
#'   `Theta` precision). Parsed by `.parse_true_theta()`.
#' @param adjacency Optional \eqn{G\times G\times J} adjacency array, or
#'   a length-\eqn{J} list of adjacencies, used for network density when
#'   supplied.
#' @param active_tol Threshold for counting an active simplex component.
#' @param include_mixsim Logical. If `TRUE` (default), compute MixSim
#'   BarOmega and MaxOmega in one overlap pass. Set `FALSE` to skip.
#'
#' @return A list with:
#' * `theta_true`: the convolution parameters `p`, `mu`, `sigma`;
#' * `descriptors`: one-row tibble of kept scenario statistics in six
#'   families (composition, mean geometry, SPD of
#'   \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}, tangent Fisher, network,
#'   component overlap). MixSim BarOmega uses [compute_average_overlap()]
#'   (Davies quadrature for \eqn{G<4}, Sobol Monte Carlo otherwise).
#'   Average pairwise Hellinger and affine-invariant Riemannian
#'   distance of the \eqn{\Sigma_j} are kept alongside.
#'   `hellinger_weighted` is optional. SPD columns include both
#'   \eqn{\kappa\{\boldsymbol{\Sigma}(\boldsymbol{p})\}}
#'   (`kappa_sigma_p`) and the reciprocal
#'   \eqn{\lambda_{\min}/\lambda_{\max}} (`kappa_sigma_reciprocal`);
#' * `supplementary`: one-row tibble of Jeffreys / symmetrised KL,
#'   recorded but not treated as a primary score;
#' * `call`: the matched call ([match.call()]).
#'
#' @examples
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(10, 0), c(0, 10)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' out <- describe_simulation_scenario(theta)
#' out$descriptors$h_star
#' out$call
#' @export
#' @seealso [run_simulation_benchmark()], [expected_fisher_unconstrained()],
#'   [compute_shannon_entropy()], [composition_from_entropy()],
#'   [helmert_basis()], [compute_average_overlap()],
#'   [compute_average_riemannian()]
describe_simulation_scenario <- function(
  true_theta,
  adjacency = NULL,
  active_tol = 1e-8,
  include_mixsim = TRUE
) {
  call <- match.call()
  theta <- .parse_true_theta(
    true_theta,
    require_p = TRUE,
    second_moment = "sigma"
  )
  p <- as.numeric(theta$p)
  mu <- theta$mu
  Sigma <- theta$sigma
  n_celltypes <- theta$J
  n_genes <- theta$G

  theta_true <- list(
    p = p,
    mu = mu,
    sigma = Sigma
  )

  h_star <- compute_shannon_entropy(p)
  h_nats <- if (n_celltypes <= 1L) {
    0
  } else {
    h_star * log(n_celltypes)
  }
  active <- p > active_tol
  n_active <- sum(active)
  min_active <- if (any(active)) {
    min(p[active])
  } else {
    NA_real_
  }

  sigma_p <- .compute_global_variance(p, Sigma)
  evals_sigma <- eigen(sigma_p, symmetric = TRUE, only.values = TRUE)$values
  lambda_min_sigma <- min(evals_sigma)
  lambda_max_sigma <- max(evals_sigma)
  kappa_sigma_reciprocal <- if (
    is.finite(lambda_max_sigma) && abs(lambda_max_sigma) > 0
  ) {
    lambda_min_sigma / lambda_max_sigma
  } else {
    NA_real_
  }

  sv <- svd(mu, nu = 0L, nv = 0L)$d
  sv <- sv[is.finite(sv) & sv > 0]
  gram_vol <- if (length(sv) == 0L) {
    0
  } else {
    exp(sum(log(sv)))
  }
  kappa_mu <- if (length(sv) >= 2L) {
    max(sv) / min(sv)
  } else if (length(sv) == 1L) {
    1
  } else {
    Inf
  }
  cos_summ <- .mean_cosine_summaries(mu)

  split <- .information_blocks(p, mu, Sigma)
  split$total <- split$mean + split$cov
  v_basis <- helmert_basis(n_celltypes)
  info_t <- crossprod(v_basis, split$total %*% v_basis)
  evals_t <- eigen(info_t, symmetric = TRUE, only.values = TRUE)$values
  lambda_min_it <- min(evals_t)
  lambda_max_it <- max(evals_t)
  kappa_it <- if (is.finite(lambda_min_it) && abs(lambda_min_it) > 0) {
    lambda_max_it / lambda_min_it
  } else {
    Inf
  }
  tr_total <- sum(diag(split$total))
  tr_cov <- sum(diag(split$cov))
  f_cov <- if (is.finite(tr_total) && tr_total > 0) {
    tr_cov / tr_total
  } else {
    NA_real_
  }

  f_cov_pairs <- numeric(0)
  if (n_celltypes >= 2L) {
    for (j in seq_len(n_celltypes - 1L)) {
      for (k in seq.int(j + 1L, n_celltypes)) {
        d <- numeric(n_celltypes)
        d[[j]] <- 1
        d[[k]] <- -1
        dir_info <- .fisher_directional(
          d = d,
          p = p,
          mean_signature_matrix = mu,
          Sigma = Sigma,
          precision = split$precision
        )
        if (is.finite(dir_info$total) && dir_info$total > 0) {
          f_cov_pairs <- c(f_cov_pairs, dir_info$cov / dir_info$total)
        }
      }
    }
  }
  f_cov_max <- if (length(f_cov_pairs) == 0L) {
    NA_real_
  } else {
    max(f_cov_pairs)
  }

  omega <- theta$Theta
  if (is.null(omega)) {
    omega <- array(NA_real_, dim = dim(Sigma))
    for (j in seq_len(n_celltypes)) {
      omega[,, j] <- tryCatch(
        solve(Sigma[,, j]),
        error = function(e) {
          matrix(NA_real_, n_genes, n_genes)
        }
      )
    }
  }

  dens_list <- vector("list", n_celltypes)
  if (!is.null(adjacency)) {
    if (is.list(adjacency)) {
      adj_array <- array(0, dim = c(n_genes, n_genes, n_celltypes))
      for (j in seq_len(n_celltypes)) {
        adj_array[,, j] <- as.matrix(adjacency[[j]])
      }
    } else {
      adj_array <- adjacency
    }
    for (j in seq_len(n_celltypes)) {
      dens_list[[j]] <- .undirected_density(adj_array[,, j])
    }
  } else {
    for (j in seq_len(n_celltypes)) {
      dens_list[[j]] <- .undirected_density(omega[,, j])
    }
  }
  network_density <- mean(vapply(dens_list, `[[`, numeric(1), "density"))
  network_mean_degree <- mean(vapply(
    dens_list,
    `[[`,
    numeric(1),
    "mean_degree"
  ))

  cor_off <- numeric(0)
  r_p <- tryCatch(
    stats::cov2cor(sigma_p),
    error = function(e) NULL
  )
  if (!is.null(r_p)) {
    cor_off <- abs(r_p[upper.tri(r_p)])
  }
  hoyer_r <- .hoyer_sparsity(cor_off)

  jeffreys <- tryCatch(
    compute_average_jeffreys(true_theta),
    error = function(e) NA_real_
  )
  hellinger_sum <- 0
  hellinger_w_sum <- 0
  weight_sum <- 0
  n_pairs_h <- 0L
  for (j in seq_len(n_celltypes - 1L)) {
    for (ell in seq.int(j + 1L, n_celltypes)) {
      h_jl <- hellinger_gaussian(
        mu[, j],
        mu[, ell],
        Sigma[,, j],
        Sigma[,, ell]
      )
      w_jl <- p[[j]] * p[[ell]]
      hellinger_sum <- hellinger_sum + h_jl
      hellinger_w_sum <- hellinger_w_sum + w_jl * h_jl
      weight_sum <- weight_sum + w_jl
      n_pairs_h <- n_pairs_h + 1L
    }
  }
  # Hellinger(f, g) = Hellinger(g, f); do not treat p_j p_k as a
  # directed divergence. The unweighted mean is the primary score.
  hellinger <- if (n_pairs_h > 0L) {
    hellinger_sum / n_pairs_h
  } else {
    NA_real_
  }
  hellinger_weighted <- if (weight_sum > 0) {
    hellinger_w_sum / weight_sum
  } else {
    NA_real_
  }

  mixsim_overlap <- NA_real_
  mixsim_max_overlap <- NA_real_
  if (isTRUE(include_mixsim)) {
    overlap <- tryCatch(
      .overlap_summary(
        true_theta,
        n_mc = 4000L,
        verbose = FALSE
      ),
      error = function(e) list(BarOmega = NA_real_, MaxOmega = NA_real_)
    )
    mixsim_overlap <- overlap$BarOmega
    mixsim_max_overlap <- overlap$MaxOmega
  }
  riemannian <- tryCatch(
    compute_average_riemannian(true_theta),
    error = function(e) NA_real_
  )

  descriptors <- tibble::tibble(
    n_genes = n_genes,
    n_celltypes = n_celltypes,
    h_star = h_star,
    n_eff = exp(h_nats),
    n_active = n_active,
    min_active = min_active,
    concentration = sum(p^2),
    mean_abs_cosine = cos_summ$mean_abs_cosine,
    min_cosine = cos_summ$min_cosine,
    max_cosine = cos_summ$max_cosine,
    mean_euclidean = .mean_euclidean_distance(mu),
    kappa_mu = kappa_mu,
    gram_volume = gram_vol,
    lambda_min_sigma_p = lambda_min_sigma,
    kappa_sigma_p = if (
      is.finite(lambda_min_sigma) && abs(lambda_min_sigma) > 0
    ) {
      lambda_max_sigma / lambda_min_sigma
    } else {
      Inf
    },
    kappa_sigma_reciprocal = kappa_sigma_reciprocal,
    lambda_min_it = lambda_min_it,
    kappa_it = kappa_it,
    f_cov = f_cov,
    f_cov_max = f_cov_max,
    network_density = network_density,
    network_mean_degree = network_mean_degree,
    hoyer_abs_correlation = hoyer_r,
    mixsim_baromega = mixsim_overlap,
    mixsim_maxomega = mixsim_max_overlap,
    hellinger = hellinger,
    hellinger_weighted = hellinger_weighted,
    riemannian_sigma = riemannian
  )

  supplementary <- tibble::tibble(
    jeffreys = jeffreys
  )

  list(
    theta_true = theta_true,
    descriptors = descriptors,
    supplementary = supplementary,
    call = call
  )
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
