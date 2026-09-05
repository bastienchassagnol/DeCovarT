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
#' @keywords internal
#' @noRd
.hellinger_gaussian <- function(mu_j, mu_l, sigma_j, sigma_l) {
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

#' Mean / covariance split of the ambient Fisher information
#'
#' For the Gaussian convolution
#' \eqn{\boldsymbol{y}\mid\boldsymbol{p}\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},\boldsymbol{\Sigma}(\boldsymbol{p}))}
#' with \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}
#' and \eqn{\boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}},
#' the ambient Fisher matrix splits as
#' \deqn{
#'   I_{jk}(\boldsymbol{p})
#'   =
#'   \underbrace{\boldsymbol{\mu}_{\cdot j}^{\mathsf{T}}
#'     \boldsymbol{\Theta}\,
#'     \boldsymbol{\mu}_{\cdot k}}_{I^{\mathrm{mean}}_{jk}}
#'   +
#'   \underbrace{2 p_j p_k
#'     \mathrm{tr}\bigl(
#'       \boldsymbol{\Theta}\boldsymbol{\Sigma}_j
#'       \boldsymbol{\Theta}\boldsymbol{\Sigma}_k
#'     \bigr)}_{I^{\mathrm{cov}}_{jk}}.
#' }
#' Along a simplex contrast \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{d}=0},
#' \deqn{
#'   I_{\mathrm{mean}}(\boldsymbol{d})
#'   =
#'   (\boldsymbol{\mu}\boldsymbol{d})^{\mathsf{T}}
#'   \boldsymbol{\Theta}
#'   (\boldsymbol{\mu}\boldsymbol{d}),
#'   \qquad
#'   I_{\mathrm{cov}}(\boldsymbol{d})
#'   =
#'   \tfrac12
#'   \bigl\|
#'     \boldsymbol{\Theta}^{1/2}
#'     (D_{\boldsymbol{d}}\boldsymbol{\Sigma})
#'     \boldsymbol{\Theta}^{1/2}
#'   \bigr\|_{F}^{2}.
#' }
#' The covariance-information fraction is
#' \eqn{f_{\mathrm{cov}}=\mathrm{tr}(I^{\mathrm{cov}})/\mathrm{tr}(I)}.
#'
#' @keywords internal
#' @noRd
.fisher_mean_cov_split <- function(p, mean_signature_matrix, Sigma) {
  n_celltypes <- length(p)
  theta <- .sigma_p_factorisation(p, Sigma)$inverse
  info_mean <- matrix(0, n_celltypes, n_celltypes)
  info_cov <- matrix(0, n_celltypes, n_celltypes)
  for (j in seq_len(n_celltypes)) {
    mu_j <- mean_signature_matrix[, j, drop = TRUE]
    sigma_j <- Sigma[,, j]
    for (k in seq_len(n_celltypes)) {
      mu_k <- mean_signature_matrix[, k, drop = TRUE]
      sigma_k <- Sigma[,, k]
      info_mean[j, k] <- .inner_product(mu_j, theta, mu_k)
      info_cov[j, k] <- 2 *
        p[[j]] *
        p[[k]] *
        sum(diag(theta %*% sigma_j %*% theta %*% sigma_k))
    }
  }
  list(
    mean = info_mean,
    cov = info_cov,
    total = info_mean + info_cov,
    precision = theta
  )
}

#' Directional Fisher information along a simplex contrast
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

#' Describe a Gaussian-convolution simulation scenario
#'
#' Summarises the generating law
#' \eqn{\boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \sum_j p_j^2\boldsymbol{\Sigma}_j)} at five layers: composition,
#' mean geometry, covariance / SPD diagnostics, network sparsity, and
#' tangent Fisher information (mean versus covariance split). MixSim
#' BarOmega and average pairwise Hellinger of the component Gaussians
#' are kept in `descriptors`. Jeffreys / symmetrised KL is returned in
#' `supplementary`.
#'
#' @param true_theta List with `p`, `mu`, and `sigma` (and optionally
#'   `Theta` precision). Parsed by [check_true_theta()] /
#'   `.parse_true_theta()`.
#' @param adjacency Optional \eqn{G\times G\times J} adjacency array, or
#'   a length-\eqn{J} list of adjacencies, used for network density when
#'   supplied.
#' @param active_tol Threshold for counting an active simplex component.
#'
#' @return A list with:
#' * `theta_true`: the convolution parameters `p`, `mu`, `sigma`;
#' * `descriptors`: one-row tibble of kept scenario statistics in six
#'   families (composition, mean geometry, SPD of
#'   \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}, tangent Fisher, network,
#'   component overlap). SPD columns include both
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
#'   [helmert_basis()]
describe_simulation_scenario <- function(
  true_theta,
  adjacency = NULL,
  active_tol = 1e-8
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

  split <- .fisher_mean_cov_split(p, mu, Sigma)
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
  hellinger <- 0
  hellinger_w <- 0
  for (j in seq_len(n_celltypes - 1L)) {
    for (ell in seq.int(j + 1L, n_celltypes)) {
      w <- p[[j]] * p[[ell]]
      hellinger <- hellinger +
        w *
          .hellinger_gaussian(
            mu[, j],
            mu[, ell],
            Sigma[,, j],
            Sigma[,, ell]
          )
      hellinger_w <- hellinger_w + w
    }
  }
  hellinger <- if (hellinger_w > 0) {
    hellinger / hellinger_w
  } else {
    NA_real_
  }

  mixsim_overlap <- NA_real_
  if (requireNamespace("MixSim", quietly = TRUE)) {
    mixsim_overlap <- tryCatch(
      MixSim::overlap(
        Pi = p,
        Mu = t(mu),
        S = Sigma
      )$BarOmega,
      error = function(e) NA_real_
    )
  }

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
    hellinger = hellinger
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
