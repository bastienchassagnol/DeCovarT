###############################################################################
# Figures for vignette 03-1, tier iii (bulk-versus-reference consistency,
# #sec-fs-range-tail and #sec-fs-subexp): five tail statements for the
# consistency statistic of one gene and of a module of dependent genes,
# against Monte Carlo.
#   vignettes/figures/fig_03_1_tail_bounds_gene.{svg,png}
#   vignettes/figures/fig_03_1_tail_bounds_module.{svg,png}
#
# Run from the repository root (system R, not renv):
#   Rscript --no-save --no-restore scripts/figures/fig_03_1_tail_bounds.R
# Runtime: about 90 s on one core (Monte Carlo dominates).
#
# Model (house notation of 03-1). One gene g, J = 3 types, bulk sample i with
# n_i = 1e4 cells and fixed composition m_ji = p_j n_i, p = (0.6, 0.3, 0.1);
# reference with C_j = (2000, 300, 30) cells. Per-cell counts of type j are
# NB(size psi = 2, mean mu_gj) in the failures parametrisation of
# eq-fs-nb-mgf; scale a_i = 1 known, no extra overdispersion.
#   y~_gi = T_gi / n_i,  mu^_gj = T_gj / C_j,
#   D_gi  = y~_gi - sum_j p_j mu^_gj,   E D = 0,
#   Var D = sum_j p_j sigma_gj^2 / n_i + sum_j p_j^2 sigma_gj^2 / C_j,
#   K_D(u) = sum_j m_ji log M_j(u / n_i) + sum_j C_j log M_j(-u p_j / C_j).
# log M_j(l) = -psi log1p(-mu_gj expm1(l) / psi), finite for l < -log pi_j.
#
# Statements (known parameters):
#   T1     Gaussian (Lindeberg-Feller) approximation of D.
#   T1-log log-normal approximation of R = y~ / sum_j p_j mu^_gj (delta method).
#   T3     Chernoff bound with Poisson CGFs (misspecified when counts are NB).
#   T4     Chernoff bound with the exact NB CGF K_D.
#   Module (N = 20 genes, sparse within-type cell-level correlation rho):
#   T1 independent CLT, T2 dependency-graph CLT, T3/T4 independent Chernoff,
#   T5 NB Chernoff after a proper colouring (Janson 2004), Jensen step.
# Ground truth: exact NB sums for one gene (J bulk and J reference NB draws per
# replicate); Gaussian copula on the NB totals for the module (approximation of
# the joint law of the sums). Checks are stopifnot() calls in section 6.
# Seed: 20261101.
###############################################################################

library(ggplot2)
library(DeCovarT) # theme_decovart_facets(), the repository theme

set.seed(20261101)

out_dir <- file.path("vignettes", "figures")
stopifnot(dir.exists(out_dir))

col_gauss <- "#1d4e89" # Gaussian family (T1, T2)
col_pois <- "#c45c26" # Poisson family (T3)
col_guide <- "#1c2b24" # Monte Carlo, guides
col_nb <- "#4f7f52" # muted green: NB family (T4, T5)
col_log <- "#7d9fc7" # muted light blue: log-normal (T1-log)

# ---- 1. setting -------------------------------------------------------------
alpha_tail <- 0.01 # alpha / 2 in each tail
n_i <- 1e4
p <- c(0.6, 0.3, 0.1)
m <- p * n_i
C <- c(2000, 300, 30)
psi <- 2
J <- length(p)
mu_gene <- list(low = c(0.2, 0.6, 0.1), high = c(20, 60, 10))
n_mc_gene <- 4e6
n_mc_module <- 5e5

var_D <- function(mu) {
  s2 <- mu + mu^2 / psi
  sum(p * s2) / n_i + sum(p^2 * s2 / C)
}
# third cumulant: per-cell NB kappa_3 = mu + 3 mu^2 / psi + 2 mu^3 / psi^2
kappa3_D <- function(mu) {
  k3 <- mu + 3 * mu^2 / psi + 2 * mu^3 / psi^2
  sum(p * k3) / n_i^2 - sum(p^3 * k3 / C^2)
}

# ---- 2. cumulant generating functions ---------------------------------------
# per-cell log-MGF and its derivative, NB (size psi) and Poisson
lmgf_nb <- function(l, mu) {
  x <- mu * expm1(l) / psi
  ifelse(x < 1, -psi * log1p(-x), Inf)
}
dlmgf_nb <- function(l, mu) mu * exp(l) / (1 - mu * expm1(l) / psi)
lmgf_pois <- function(l, mu) mu * expm1(l)
dlmgf_pois <- function(l, mu) mu * exp(l)

# CGF of D_g as a list (K, dK, domain lo < 0 < hi)
cgf_D <- function(mu, family = c("nb", "pois")) {
  family <- match.arg(family)
  if (family == "nb") {
    f <- lmgf_nb
    df <- dlmgf_nb
    b <- log1p(psi / mu) # -log pi_j
    hi <- n_i * min(b)
    lo <- max(-C * b / p)
  } else {
    f <- lmgf_pois
    df <- dlmgf_pois
    hi <- Inf
    lo <- -Inf
  }
  list(
    K = function(u) sum(m * f(u / n_i, mu)) + sum(C * f(-u * p / C, mu)),
    dK = function(u) sum(p * df(u / n_i, mu)) - sum(p * df(-u * p / C, mu)),
    lo = lo,
    hi = hi
  )
}

# Legendre transform sup_{side * u > 0} [side * u * t - K(u)], t >= 0, found
# from the stationarity condition K'(u) = side * t (K is convex).
chernoff_rate <- function(cgf, t, side = 1) {
  if (t <= 0) {
    return(0)
  }
  edge <- if (side > 0) cgf$hi else cgf$lo
  if (is.finite(edge)) {
    edge <- edge * (1 - 1e-12)
  } else {
    edge <- side
    while (side * cgf$dK(edge) < t) {
      edge <- 2 * edge
    }
  }
  g <- function(u) cgf$dK(u) - side * t
  stopifnot(side * g(edge) > 0)
  u_star <- stats::uniroot(
    g,
    sort(c(0, edge)),
    tol = 1e-13 * abs(edge),
    maxiter = 2000
  )$root
  side * u_star * t - cgf$K(u_star)
}

# smallest t with exp(-rate(t)) <= level
chernoff_threshold <- function(cgf, side, t0, level = alpha_tail) {
  target <- log(1 / level)
  hi <- t0
  while (chernoff_rate(cgf, hi, side) < target) {
    hi <- 2 * hi
  }
  stats::uniroot(
    function(t) chernoff_rate(cgf, t, side) - target,
    c(0, hi),
    tol = 1e-10 * hi
  )$root
}

# ---- 3. one gene: statements on a grid of t ---------------------------------
gene_statements <- function(mu, t_grid, side) {
  s <- sqrt(var_D(mu))
  mbar <- sum(p * mu)
  nb <- cgf_D(mu, "nb")
  po <- cgf_D(mu, "pois")
  t1_log <- if (side > 0) {
    stats::pnorm(mbar * log1p(t_grid / mbar) / s, lower.tail = FALSE)
  } else {
    ifelse(
      t_grid < mbar,
      stats::pnorm(mbar * log1p(-pmin(t_grid, mbar) / mbar) / s),
      0
    )
  }
  rbind(
    data.frame(
      t = t_grid,
      prob = stats::pnorm(t_grid / s, lower.tail = FALSE),
      statement = "T1"
    ),
    data.frame(t = t_grid, prob = t1_log, statement = "T1-log"),
    data.frame(
      t = t_grid,
      prob = exp(-vapply(t_grid, \(t) chernoff_rate(po, t, side), 0)),
      statement = "T3"
    ),
    data.frame(
      t = t_grid,
      prob = exp(-vapply(t_grid, \(t) chernoff_rate(nb, t, side), 0)),
      statement = "T4"
    )
  )
}

gene_thresholds <- function(mu, side, level = alpha_tail) {
  s <- sqrt(var_D(mu))
  mbar <- sum(p * mu)
  z <- stats::qnorm(1 - level)
  c(
    "T1" = z * s,
    "T1-log" = if (side > 0) {
      mbar * expm1(z * s / mbar)
    } else {
      -mbar * expm1(-z * s / mbar)
    },
    "T3" = chernoff_threshold(cgf_D(mu, "pois"), side, s, level),
    "T4" = chernoff_threshold(cgf_D(mu, "nb"), side, s, level)
  )
}

# ---- 4. one gene: Monte Carlo (exact NB sums) --------------------------------
# A sum of m i.i.d. NB(psi, mu) cells is NB(m psi, mean m mu) (common pi), so
# a replicate needs J bulk draws and J reference draws.
sim_D <- function(mu, n_mc) {
  bulk <- numeric(n_mc)
  ref <- numeric(n_mc)
  for (j in seq_len(J)) {
    bulk <- bulk + stats::rnbinom(n_mc, size = m[j] * psi, mu = m[j] * mu[j])
    ref <- ref +
      p[j] * stats::rnbinom(n_mc, size = C[j] * psi, mu = C[j] * mu[j]) / C[j]
  }
  list(D = bulk / n_i - ref, logR = log(bulk / n_i) - log(ref))
}

# empirical tail on a grid, from the sorted sample
mc_tail <- function(x_sorted, t_grid, side) {
  n <- length(x_sorted)
  if (side > 0) {
    (n - findInterval(t_grid, x_sorted, left.open = TRUE)) / n
  } else {
    findInterval(-t_grid, x_sorted) / n
  }
}

# empirical threshold (inverse ECDF) with a 95% order-statistic interval
mc_threshold <- function(x_sorted, side, level = alpha_tail) {
  n <- length(x_sorted)
  q <- if (side > 0) 1 - level else level
  k <- ceiling(n * q)
  h <- ceiling(1.96 * sqrt(n * q * (1 - q)))
  v <- x_sorted[c(k, max(1, k - h), min(n, k + h))]
  if (side > 0) v else -v[c(1, 3, 2)]
}

gene_curves <- list()
gene_thr <- list()
gene_mc <- list()
for (lev in names(mu_gene)) {
  mu <- mu_gene[[lev]]
  sim <- sim_D(mu, n_mc_gene)
  D_sorted <- sort(sim$D)
  gene_mc[[lev]] <- c(
    mean = mean(sim$D),
    var = stats::var(sim$D),
    var_formula = var_D(mu),
    sd_logR = stats::sd(sim$logR[is.finite(sim$logR)]),
    sd_over_mbar = sqrt(var_D(mu)) / sum(p * mu),
    skew = mean((sim$D - mean(sim$D))^3) / stats::var(sim$D)^1.5,
    skew_formula = kappa3_D(mu) / var_D(mu)^1.5
  )
  rm(sim)
  for (side in c(-1, 1)) {
    tail_name <- if (side > 0) "upper" else "lower"
    t_far <- max(gene_thresholds(mu, side, level = 1e-5))
    t_grid <- seq(0, 1.02 * t_far, length.out = 401)
    d <- gene_statements(mu, t_grid, side)
    d <- rbind(
      d,
      data.frame(
        t = t_grid,
        prob = mc_tail(D_sorted, t_grid, side),
        statement = "MC"
      )
    )
    d$expression <- lev
    d$tail <- tail_name
    gene_curves[[paste(lev, tail_name)]] <- d
    thr <- gene_thresholds(mu, side)
    mc_thr <- mc_threshold(D_sorted, side)
    gene_thr[[paste(lev, tail_name)]] <- data.frame(
      unit = paste(lev, "expression gene"),
      tail = tail_name,
      statement = c(names(thr), "MC"),
      t_001 = c(unname(thr), mc_thr[1]),
      mc_lo = mc_thr[2],
      mc_hi = mc_thr[3]
    )
  }
  rm(D_sorted)
}
gene_curves <- do.call(rbind, gene_curves)
gene_thr <- do.call(rbind, gene_thr)

# ---- 5. module of N dependent genes -----------------------------------------
N <- 20L
n_edges <- 30L # mean degree 2|E| / N = 3
max_degree <- 10L
rho <- 0.3 # within-type cell-level correlation on edges
sd_ln <- 0.5 # log-normal spread of the type means around the two genes
base_mu <- rbind(
  matrix(mu_gene$low, N / 2, J, byrow = TRUE),
  matrix(mu_gene$high, N / 2, J, byrow = TRUE)
)
mu_mod <- base_mu * exp(sd_ln * matrix(stats::rnorm(N * J), N, J))

all_pairs <- t(utils::combn(N, 2))
repeat {
  E <- all_pairs[sample.int(nrow(all_pairs), n_edges), , drop = FALSE]
  A <- matrix(0, N, N)
  A[E] <- 1
  A[E[, 2:1]] <- 1
  if (max(rowSums(A)) <= max_degree) break
}
deg <- rowSums(A)
R_cop <- diag(N) + rho * A # cell-level (and copula) correlation, every type
eig_R <- eigen(R_cop, symmetric = TRUE, only.values = TRUE)$values

# greedy (Welsh-Powell) proper colouring, vertices by decreasing degree
colour <- integer(N)
for (v in order(-deg, seq_len(N))) {
  colour[v] <- min(setdiff(seq_len(N), colour[A[v, ] == 1]))
}
k_col <- max(colour)
classes <- split(seq_len(N), colour)
bound_degree <- max(deg) + 1
bound_edges <- (1 + sqrt(8 * n_edges + 1)) / 2

# covariance of D across genes: Cov = sum_j (p_j/n + p_j^2/C_j) c_ghj,
# c_ghj = R_gh sigma_gj sigma_hj
sigma_mod <- sqrt(mu_mod + mu_mod^2 / psi)
w_j <- p / n_i + p^2 / C
cov_D <- (sigma_mod %*% diag(w_j) %*% t(sigma_mod)) * R_cop
sd_mod <- sqrt(diag(cov_D))
corr_D <- stats::cov2cor(cov_D)
var_S_t1 <- N
var_S_t2 <- sum(corr_D) # N + 2 sum_{edges} Corr(D_g, D_h)

# CGFs of S = sum_g D_g / sd_g
cgf_list_nb <- lapply(seq_len(N), \(g) cgf_D(mu_mod[g, ], "nb"))
cgf_list_po <- lapply(seq_len(N), \(g) cgf_D(mu_mod[g, ], "pois"))
cgf_S_indep <- function(cl) {
  list(
    K = function(u) sum(vapply(seq_len(N), \(g) cl[[g]]$K(u / sd_mod[g]), 0)),
    dK = function(u) {
      sum(vapply(seq_len(N), \(g) cl[[g]]$dK(u / sd_mod[g]) / sd_mod[g], 0))
    },
    lo = max(vapply(seq_len(N), \(g) cl[[g]]$lo * sd_mod[g], 0)),
    hi = min(vapply(seq_len(N), \(g) cl[[g]]$hi * sd_mod[g], 0))
  )
}
# T5: L(u) = log(k^-1 sum_c exp(sum_{g in I_c} K_{Z_g}(k u))) >= log E e^{uS}
cgf_S_colour <- function(cl) {
  a_c <- function(u) {
    vapply(
      classes,
      \(I) sum(vapply(I, \(g) cl[[g]]$K(k_col * u / sd_mod[g]), 0)),
      0
    )
  }
  da_c <- function(u) {
    vapply(
      classes,
      \(I) {
        sum(vapply(
          I,
          \(g) k_col * cl[[g]]$dK(k_col * u / sd_mod[g]) / sd_mod[g],
          0
        ))
      },
      0
    )
  }
  list(
    K = function(u) {
      a <- a_c(u)
      max(a) + log(mean(exp(a - max(a))))
    },
    dK = function(u) {
      a <- a_c(u)
      w <- exp(a - max(a))
      sum(w * da_c(u)) / sum(w)
    },
    lo = max(vapply(seq_len(N), \(g) cl[[g]]$lo * sd_mod[g], 0)) / k_col,
    hi = min(vapply(seq_len(N), \(g) cl[[g]]$hi * sd_mod[g], 0)) / k_col
  )
}
cgf_t3 <- cgf_S_indep(cgf_list_po)
cgf_t4 <- cgf_S_indep(cgf_list_nb)
cgf_t5 <- cgf_S_colour(cgf_list_nb)

module_thresholds <- function(side, level = alpha_tail) {
  z <- stats::qnorm(1 - level)
  c(
    "T1" = z * sqrt(var_S_t1),
    "T2" = z * sqrt(var_S_t2),
    "T3" = chernoff_threshold(cgf_t3, side, sqrt(N), level),
    "T4" = chernoff_threshold(cgf_t4, side, sqrt(N), level),
    "T5" = chernoff_threshold(cgf_t5, side, sqrt(N), level)
  )
}

# Monte Carlo: Gaussian copula (correlation R_cop) on the NB totals of each
# type, drawn independently for the bulk and for the reference.
qnb_from_z <- function(z, size, mu) {
  out <- numeric(length(z))
  neg <- z < 0
  out[neg] <- stats::qnbinom(stats::pnorm(z[neg]), size = size, mu = mu)
  out[!neg] <- stats::qnbinom(
    stats::pnorm(z[!neg], lower.tail = FALSE),
    size = size,
    mu = mu,
    lower.tail = FALSE
  )
  out
}
chol_R <- chol(R_cop)
draw_totals <- function(n_cells, mu_j) {
  z <- matrix(stats::rnorm(n_mc_module * N), n_mc_module, N) %*% chol_R
  for (g in seq_len(N)) {
    z[, g] <- qnb_from_z(z[, g], size = n_cells * psi, mu = n_cells * mu_j[g])
  }
  z
}
D_mod <- matrix(0, n_mc_module, N)
for (j in seq_len(J)) {
  D_mod <- D_mod + draw_totals(m[j], mu_mod[, j]) / n_i
  D_mod <- D_mod - p[j] * draw_totals(C[j], mu_mod[, j]) / C[j]
}
Z_mod <- sweep(D_mod, 2, sd_mod, "/")
rm(D_mod)
S_mc <- rowSums(Z_mod)
S_sorted <- sort(S_mc)
corr_mc <- stats::cor(Z_mod)

module_curves <- list()
module_thr <- list()
for (side in c(-1, 1)) {
  tail_name <- if (side > 0) "upper" else "lower"
  t_far <- max(module_thresholds(side, level = 1e-4))
  t_grid <- seq(0, 1.02 * t_far, length.out = 301)
  rate_curve <- function(cgf) {
    exp(-vapply(t_grid, \(t) chernoff_rate(cgf, t, side), 0))
  }
  d <- rbind(
    data.frame(
      t = t_grid,
      prob = stats::pnorm(t_grid / sqrt(var_S_t1), lower.tail = FALSE),
      statement = "T1"
    ),
    data.frame(
      t = t_grid,
      prob = stats::pnorm(t_grid / sqrt(var_S_t2), lower.tail = FALSE),
      statement = "T2"
    ),
    data.frame(t = t_grid, prob = rate_curve(cgf_t3), statement = "T3"),
    data.frame(t = t_grid, prob = rate_curve(cgf_t4), statement = "T4"),
    data.frame(t = t_grid, prob = rate_curve(cgf_t5), statement = "T5"),
    data.frame(
      t = t_grid,
      prob = mc_tail(S_sorted, t_grid, side),
      statement = "MC"
    )
  )
  d$tail <- tail_name
  module_curves[[tail_name]] <- d
  thr <- module_thresholds(side)
  mc_thr <- mc_threshold(S_sorted, side)
  module_thr[[tail_name]] <- data.frame(
    unit = "module",
    tail = tail_name,
    statement = c(names(thr), "MC"),
    t_001 = c(unname(thr), mc_thr[1]),
    mc_lo = mc_thr[2],
    mc_hi = mc_thr[3]
  )
}
module_curves <- do.call(rbind, module_curves)
module_thr <- do.call(rbind, module_thr)

# ---- 6. checks ---------------------------------------------------------------
rel <- function(a, b) abs(a - b) / abs(b)

# (a) numeric Legendre transform against closed forms, single sum of n cells
for (mu0 in c(0.2, 20)) {
  for (a in mu0 * c(0.3, 0.7, 1.4, 3)) {
    side <- sign(a - mu0)
    # Poisson: n (a log(a / mu) - a + mu)
    cg <- list(
      K = function(u) n_i * mu0 * expm1(u / n_i) - u * mu0,
      dK = function(u) mu0 * exp(u / n_i) - mu0,
      lo = -Inf,
      hi = Inf
    )
    stopifnot(
      rel(
        chernoff_rate(cg, abs(a - mu0), side),
        n_i * (a * log(a / mu0) - a + mu0)
      ) <
        1e-8
    )
    # NB: n I(a), eq-fs-nb-rate
    cg <- list(
      K = function(u) n_i * lmgf_nb(u / n_i, mu0) - u * mu0,
      dK = function(u) dlmgf_nb(u / n_i, mu0) - mu0,
      lo = -Inf,
      hi = n_i * log1p(psi / mu0)
    )
    I_a <- a *
      log(a * (mu0 + psi) / (mu0 * (a + psi))) +
      psi * log((mu0 + psi) / (a + psi))
    stopifnot(rel(chernoff_rate(cg, abs(a - mu0), side), n_i * I_a) < 1e-8)
  }
}

# (b) K_D(0) = 0, K_D'(0) = 0, K_D''(0) = Var D (central differences)
for (mu in mu_gene) {
  cg <- cgf_D(mu, "nb")
  h <- 1e-3 / sqrt(var_D(mu))
  stopifnot(
    abs(cg$K(0)) < 1e-12,
    abs(cg$dK(0)) < 1e-12,
    rel((cg$dK(h) - cg$dK(-h)) / (2 * h), var_D(mu)) < 1e-6,
    rel((cg$K(h) - 2 * cg$K(0) + cg$K(-h)) / h^2, var_D(mu)) < 1e-4,
    rel((cg$dK(h) - 2 * cg$dK(0) + cg$dK(-h)) / h^2, kappa3_D(mu)) < 1e-3
  )
}

# (c) Monte Carlo moments of D against the formulas
for (lev in names(gene_mc)) {
  g <- gene_mc[[lev]]
  stopifnot(
    abs(g["mean"]) < 5 * sqrt(g["var_formula"] / n_mc_gene),
    rel(g["var"], g["var_formula"]) < 0.01,
    abs(g["skew"] - g["skew_formula"]) < 0.01 # SE about sqrt(6 / n) = 0.001
  )
}

# (d) graph, positive definiteness, colouring
stopifnot(
  nrow(unique(t(apply(E, 1, sort)))) == n_edges,
  max(deg) <= max_degree,
  mean(deg) <= 3,
  min(eig_R) > 0,
  all(colour[E[, 1]] != colour[E[, 2]]),
  k_col <= min(bound_degree, bound_edges),
  # each colour class has an identity copula block: genes independent
  all(vapply(
    classes,
    \(I) isTRUE(all.equal(R_cop[I, I, drop = FALSE], diag(length(I)))),
    TRUE
  ))
)

# (e) Jensen step. Deterministic: e^{uS} <= k^-1 sum_c e^{k u S_c} on every
# draw. Exact Gaussian toy with the module correlation: exp(u^2 1'R1 / 2) <=
# k^-1 sum_c exp(|I_c| k^2 u^2 / 2). Monte Carlo: log E e^{uS} <= L(u).
S_class <- vapply(
  classes,
  \(I) rowSums(Z_mod[, I, drop = FALSE]),
  numeric(n_mc_module)
)
jensen_draws <- vapply(
  c(-0.2, -0.05, 0.05, 0.2),
  \(u) {
    a <- k_col * u * S_class
    rhs <- apply(a, 1, max) + log(rowMeans(exp(a - apply(a, 1, max))))
    all(u * S_mc <= rhs + 1e-12)
  },
  TRUE
)
u_toy <- seq(0.01, 1, by = 0.01)
jensen_toy <- all(
  u_toy^2 * sum(corr_D) / 2 <=
    log(vapply(
      u_toy,
      \(u) mean(exp(lengths(classes) * k_col^2 * u^2 / 2)),
      0
    ))
)
jensen_mc <- t(vapply(
  c(-0.1, -0.05, 0.05, 0.1),
  \(u) {
    c(
      u = u,
      log_mgf_mc = log(mean(exp(u * S_mc))),
      L_colour = cgf_t5$K(u),
      K_indep = cgf_t4$K(u)
    )
  },
  numeric(4)
))
stopifnot(
  all(jensen_draws),
  jensen_toy,
  all(jensen_mc[, "log_mgf_mc"] <= jensen_mc[, "L_colour"])
)

# (f) bounds T4 (gene) and T5 (module) are not below Monte Carlo at t_0.01
chk <- rbind(gene_thr, module_thr)
for (key in unique(paste(chk$unit, chk$tail))) {
  s <- chk[paste(chk$unit, chk$tail) == key, ]
  bnd <- if (grepl("module", key)) "T5" else "T4"
  stopifnot(s$t_001[s$statement == bnd] >= s$mc_lo[s$statement == "MC"])
}

# ---- 7. report ---------------------------------------------------------------
cat("\nGene Monte Carlo moments (n =", n_mc_gene, "):\n")
print(do.call(rbind, gene_mc), digits = 5)
cat(
  "\nModule: N =",
  N,
  " |E| =",
  n_edges,
  " max degree =",
  max(deg),
  " mean degree =",
  mean(deg),
  "\n"
)
cat("rho =", rho, " min eigenvalue of I + rho A =", min(eig_R), "\n")
cat(
  "greedy colouring: k =",
  k_col,
  " bounds Delta + 1 =",
  bound_degree,
  " (1 + sqrt(8|E| + 1))/2 =",
  round(bound_edges, 3),
  "\n"
)
cat("colour class sizes:", lengths(classes), "\n")
cat(
  "Var S: T1 =",
  var_S_t1,
  " T2 =",
  round(var_S_t2, 3),
  " Monte Carlo =",
  round(stats::var(S_mc), 3),
  "\n"
)
edge_idx <- E
cat(
  "edge correlations Corr(D_g, D_h): formula range",
  round(range(corr_D[edge_idx]), 3),
  "; Monte Carlo range",
  round(range(corr_mc[edge_idx]), 3),
  "\n"
)
cat(
  "max |off-edge Monte Carlo correlation| =",
  round(max(abs(corr_mc[upper.tri(corr_mc) & A == 0])), 4),
  "\n"
)
cat(
  "Jensen checks: per draw",
  all(jensen_draws),
  "; Gaussian toy",
  jensen_toy,
  "\n"
)
print(jensen_mc, digits = 5)
cat("\nModule type means (genes x types):\n")
print(round(mu_mod, 3))
cat("\nEdges:\n")
print(E)

thr_all <- rbind(gene_thr, module_thr)
thr_all$mc_value <- ave(
  ifelse(thr_all$statement == "MC", thr_all$t_001, NA),
  thr_all$unit,
  thr_all$tail,
  FUN = \(x) max(x, na.rm = TRUE)
)
thr_all$ratio_to_mc <- thr_all$t_001 / thr_all$mc_value
thr_all$under_covers <- thr_all$statement != "MC" &
  thr_all$t_001 < thr_all$mc_lo
cat("\nThresholds t_0.01 (MC with 95% order-statistic interval):\n")
print(thr_all, digits = 5, row.names = FALSE)

# ---- 8. figures ---------------------------------------------------------------
y_breaks <- 10^(0:-6)
y_labels <- c(
  "1",
  "0.1",
  "0.01",
  expression(10^-3),
  expression(10^-4),
  expression(10^-5),
  expression(10^-6)
)

theme_fig <- theme_decovart_facets(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = grid::unit(1.6, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e1d8", linewidth = 0.3),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "white", size = 9),
    plot.margin = margin(4, 8, 4, 4)
  )

save_figure <- function(plot, name, width, height) {
  svglite::svglite(
    file.path(out_dir, paste0(name, ".svg")),
    width = width,
    height = height,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    file.path(out_dir, paste0(name, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 220,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(NULL)
}

tail_plot <- function(
  curves,
  thr,
  keys,
  labels,
  colours,
  linetypes,
  y_min,
  facet,
  legend_order,
  legend_rows,
  legend_byrow
) {
  curves <- curves[curves$prob > 0, ]
  curves$statement <- factor(curves$statement, levels = keys, labels = labels)
  thr$statement <- factor(thr$statement, levels = keys, labels = labels)
  is_mc <- function(d) d$statement == labels[keys == "MC"]
  leg <- guide_legend(nrow = legend_rows, byrow = legend_byrow)
  ggplot(curves, aes(t, prob, colour = statement, linetype = statement)) +
    geom_hline(
      yintercept = alpha_tail,
      colour = col_guide,
      linewidth = 0.35,
      linetype = "22"
    ) +
    geom_line(data = curves[is_mc(curves), ], linewidth = 1.6, alpha = 0.35) +
    geom_line(data = curves[!is_mc(curves), ], linewidth = 0.75) +
    # thresholds t_0.01: Monte Carlo diamond underneath, statements on top
    geom_point(
      data = thr[is_mc(thr), ],
      aes(t_001, alpha_tail),
      shape = 23,
      size = 3.4,
      stroke = 0.4,
      fill = col_guide,
      colour = "white",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = thr[!is_mc(thr), ],
      aes(t_001, alpha_tail, fill = statement),
      shape = 21,
      size = 2,
      stroke = 0.45,
      colour = "white",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    facet +
    scale_colour_manual(values = colours, breaks = legend_order) +
    scale_fill_manual(values = colours, breaks = legend_order) +
    scale_linetype_manual(values = linetypes, breaks = legend_order) +
    scale_y_log10(breaks = y_breaks, labels = y_labels) +
    coord_cartesian(ylim = c(y_min, 1), expand = FALSE) +
    labs(x = "deviation t", y = "tail probability") +
    guides(colour = leg, linetype = leg) +
    theme_fig
}

le <- "≤"
ge <- "≥"
minus <- "−"
tail_label <- function(stat, tail) {
  ifelse(
    tail == "lower",
    paste0("lower tail P(", stat, " ", le, " ", minus, "t)"),
    paste0("upper tail P(", stat, " ", ge, " t)")
  )
}

# gene figure: rows low / high expression, columns lower / upper tail
gene_keys <- c("T1", "T1-log", "T3", "T4", "MC")
gene_labels <- c(
  "T1 Gaussian (CLT)",
  "T1-log log-normal (delta method)",
  "T3 Poisson Chernoff",
  "T4 NB Chernoff",
  "Monte Carlo (exact NB sums)"
)
panel_of <- function(expression, tail) {
  lab <- paste0(expression, " expression, ", tail_label("D", tail))
  lev <- paste0(
    rep(c("low", "high"), each = 2),
    " expression, ",
    tail_label("D", rep(c("lower", "upper"), 2))
  )
  factor(lab, levels = lev)
}
gene_curves$panel <- panel_of(gene_curves$expression, gene_curves$tail)
gene_thr$panel <- panel_of(
  sub(" expression gene", "", gene_thr$unit),
  gene_thr$tail
)
p_gene <- tail_plot(
  gene_curves,
  gene_thr,
  keys = gene_keys,
  labels = gene_labels,
  colours = stats::setNames(
    c(col_gauss, col_log, col_pois, col_nb, col_guide),
    gene_labels
  ),
  linetypes = stats::setNames(
    c("solid", "42", "solid", "solid", "solid"),
    gene_labels
  ),
  y_min = 1e-5,
  facet = facet_wrap(~panel, ncol = 2, scales = "free_x"),
  legend_order = gene_labels,
  legend_rows = 2,
  legend_byrow = TRUE
)
save_figure(p_gene, "fig_03_1_tail_bounds_gene", width = 7.4, height = 5.9)

# module figure: lower and upper tail of S
module_keys <- c("T1", "T2", "T3", "T4", "T5", "MC")
module_labels <- c(
  "T1 independent CLT, Var S = N",
  "T2 dependent CLT, Var S = N + 2Σ Corr",
  "T3 Poisson Chernoff, independent",
  "T4 NB Chernoff, independent",
  paste0("T5 NB Chernoff, colouring (k = ", k_col, ")"),
  "Monte Carlo (Gaussian copula)"
)
module_curves$panel <- factor(
  tail_label("S", module_curves$tail),
  levels = tail_label("S", c("lower", "upper"))
)
module_thr$panel <- factor(
  tail_label("S", module_thr$tail),
  levels = levels(module_curves$panel)
)
p_module <- tail_plot(
  module_curves,
  module_thr,
  keys = module_keys,
  labels = module_labels,
  colours = stats::setNames(
    c(col_gauss, col_gauss, col_pois, col_nb, col_nb, col_guide),
    module_labels
  ),
  linetypes = stats::setNames(
    c("22", "solid", "22", "22", "solid", "solid"),
    module_labels
  ),
  y_min = 1e-4,
  facet = facet_wrap(~panel, ncol = 2),
  # left column: independence assumed (dashed); right: dependence handled
  legend_order = module_labels[c(1, 3, 4, 2, 5, 6)],
  legend_rows = 3,
  legend_byrow = FALSE
) +
  labs(x = "deviation t of S (sum of the standardised D_g)")
save_figure(
  p_module,
  "fig_03_1_tail_bounds_module",
  width = 7.4,
  height = 4.1
)

utils::sessionInfo()
