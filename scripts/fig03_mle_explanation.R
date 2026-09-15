###############################################################################
###############################################################################
###                                                                         ###
###     FIGURE 03 – LLM pack: MLE explanation (balanced compositions)      ###
###                                                                         ###
###############################################################################
###############################################################################
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore scripts/fig03_mle_explanation.R \
#     > "logs/fig03_mle_explanation_$(date +%F).log" 2>&1 &
#
# Reads finished hybrid_*.rds (no ADEMP refit). Writes CSV + prompt under
# output/fig03/mle_explanation/ for an LLM agent (same layout as
# output/fig02/mle_explanation/).
###############################################################################

if (
  requireNamespace("devtools", quietly = TRUE) &&
    file.exists("DESCRIPTION")
) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(DeCovarT)
}
DeCovarT:::.ui_attach_script()

OUT_DIR <- file.path("output", "fig03", "mle_explanation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.ui_h1("Figure 03 · MLE explanation pack (equi-balanced p)")

artefacts <- read_simulation_artefacts(
  file.path("output", "fig03"),
  "hybrid",
  assemble = TRUE
)
artefacts$theta <- readRDS(file.path("output", "fig03", "hybrid_theta.rds"))
cfg <- artefacts$config
cfg <- cfg[as.character(cfg$proportions) == "balanced", , drop = FALSE]
ids <- as.character(cfg$ID)
.ui_info("Balanced scenarios: {.val {length(ids)}}.")

desc <- artefacts$descriptors
if (is.null(desc)) {
  desc <- readRDS(file.path("output", "fig03", "hybrid_descriptors.rds"))
}
desc <- desc[as.character(desc$ID) %in% ids, , drop = FALSE]

.ilr_info <- function(p, mu, Sigma) {
  info_p <- expected_fisher_unconstrained(p, mu, Sigma)
  z <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(z)
  info_z <- t(jac) %*% info_p %*% jac
  vcov_z <- tryCatch(
    solve(info_z),
    error = function(e) {
      matrix(NA_real_, length(z), length(z))
    }
  )
  vcov_p <- tryCatch(
    vcov_ilr_delta(p, mu, Sigma, warn = FALSE),
    error = function(e) {
      matrix(NA_real_, length(p), length(p))
    }
  )
  list(
    z = as.numeric(z),
    se_p = sqrt(pmax(diag(as.matrix(vcov_p)), 0)),
    se_z = sqrt(pmax(diag(as.matrix(vcov_z)), 0)),
    vcov_z = vcov_z,
    kappa_info_z = {
      ev <- eigen(info_z, symmetric = TRUE, only.values = TRUE)$values
      ev <- ev[is.finite(ev) & ev > 0]
      if (length(ev) < 1L) {
        NA_real_
      } else {
        max(ev) / min(ev)
      }
    }
  )
}

.theta_for_id <- function(id) {
  hit <- as.character(artefacts$theta$ID) == id
  DeCovarT:::.unwrap_true_theta(
    artefacts$theta$true_theta[which(hit)[[1L]]]
  )
}

.serialise_generators <- function(x) {
  if (is.null(x) || length(x) < 1L) {
    return(NA_character_)
  }
  if (is.list(x) && length(x) == 1L && is.list(x[[1L]])) {
    x <- x[[1L]]
  }
  jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")
}

# --------------------------------------------------------------------------
# Scenario descriptors + theoretical Wald SEs ----
# --------------------------------------------------------------------------

desc_rows <- vector("list", nrow(cfg))
wald_rows <- vector("list", nrow(cfg))
profile_chunks <- vector("list", nrow(cfg))
rho_lim <- 3
n_grid <- 25L
n_axis <- 81L
rho_grid <- seq(-rho_lim, rho_lim, length.out = n_grid)
rho_axis <- seq(-rho_lim, rho_lim, length.out = n_axis)

for (i in seq_len(nrow(cfg))) {
  id <- ids[[i]]
  row <- cfg[i, , drop = FALSE]
  th <- .theta_for_id(id)
  p <- as.numeric(th$p)
  nms <- names(th$p)
  if (is.null(nms)) {
    nms <- paste0("celltype_", seq_along(p))
  }
  names(p) <- nms
  mu <- th$mu
  Sigma <- th$sigma
  y <- drop(mu %*% p)
  inf <- .ilr_info(p, mu, Sigma)
  z_true <- inf$z
  d_hit <- desc[as.character(desc$ID) == id, , drop = FALSE]
  gen <- if ("graph_generators" %in% names(row)) {
    .serialise_generators(row$graph_generators)
  } else {
    NA_character_
  }
  desc_rows[[i]] <- data.frame(
    ID = id,
    n_genes = row$n_genes,
    n_celltypes = row$n_celltypes,
    proportions = as.character(row$proportions),
    entropy = row$entropy,
    graph_ct1 = as.character(row$graph_ct1),
    graph_ct2 = as.character(row$graph_ct2),
    graph_ct3 = as.character(row$graph_ct3),
    overlap_label = as.character(row$overlap_label),
    overlap_target = row$overlap_target,
    covariance_scale = row$covariance_scale,
    graph_generators = as.character(gen),
    p1 = p[[1L]],
    p2 = p[[2L]],
    p3 = p[[3L]],
    rho1_true = z_true[[1L]],
    rho2_true = z_true[[2L]],
    y_bar_l2 = sqrt(sum(y^2)),
    stringsAsFactors = FALSE
  )
  if (nrow(d_hit) > 0L) {
    keep <- setdiff(names(d_hit), names(desc_rows[[i]]))
    desc_rows[[i]] <- cbind(desc_rows[[i]], d_hit[, keep, drop = FALSE])
  }
  wald_rows[[i]] <- data.frame(
    ID = id,
    graph_ct1 = as.character(row$graph_ct1),
    graph_ct2 = as.character(row$graph_ct2),
    overlap_label = as.character(row$overlap_label),
    p1 = p[[1L]],
    p2 = p[[2L]],
    p3 = p[[3L]],
    rho1_true = z_true[[1L]],
    rho2_true = z_true[[2L]],
    theoretical_se_p1 = inf$se_p[[1L]],
    theoretical_se_p2 = inf$se_p[[2L]],
    theoretical_se_p3 = inf$se_p[[3L]],
    theoretical_se_rho1 = inf$se_z[[1L]],
    theoretical_se_rho2 = inf$se_z[[2L]],
    cov_rho1_rho2 = inf$vcov_z[1L, 2L],
    kappa_info_ilr = inf$kappa_info_z,
    stringsAsFactors = FALSE
  )
  ll_at <- function(z1, z2) {
    suppressWarnings(
      tryCatch(
        loglik_multivariate_constrained(c(z1, z2), y, mu, Sigma),
        error = function(e) NA_real_
      )
    )
  }
  g <- expand.grid(rho1 = rho_grid, rho2 = rho_grid, KEEP.OUT.ATTRS = FALSE)
  g$loglik <- mapply(ll_at, g$rho1, g$rho2)
  g$slice <- "grid"
  ax1 <- data.frame(
    rho1 = rho_axis,
    rho2 = z_true[[2L]],
    loglik = vapply(
      rho_axis,
      function(r) ll_at(r, z_true[[2L]]),
      numeric(1)
    ),
    slice = "rho1_axis"
  )
  ax2 <- data.frame(
    rho1 = z_true[[1L]],
    rho2 = rho_axis,
    loglik = vapply(
      rho_axis,
      function(r) ll_at(z_true[[1L]], r),
      numeric(1)
    ),
    slice = "rho2_axis"
  )
  pr <- rbind(g, ax1, ax2)
  pr$ID <- id
  pr$graph_ct1 <- as.character(row$graph_ct1)
  pr$graph_ct2 <- as.character(row$graph_ct2)
  pr$overlap_label <- as.character(row$overlap_label)
  pr$is_true_rho <- abs(pr$rho1 - z_true[[1L]]) < 1e-10 &
    abs(pr$rho2 - z_true[[2L]]) < 1e-10
  profile_chunks[[i]] <- pr[
    ,
    c(
      "ID",
      "graph_ct1",
      "graph_ct2",
      "overlap_label",
      "slice",
      "rho1",
      "rho2",
      "loglik",
      "is_true_rho"
    )
  ]
  .ui_info("Profiled {.val {id}} ({.val {i}}/{.val {nrow(cfg)}}).")
}

scenario_descriptors <- dplyr::bind_rows(desc_rows)
theoretical_wald_se <- dplyr::bind_rows(wald_rows)
ilr_loglik_profiles <- dplyr::bind_rows(profile_chunks)

# --------------------------------------------------------------------------
# Forest: Newton–Raphson and Marquardt–Levenberg ----
# --------------------------------------------------------------------------

mc <- artefacts$monte_carlo
mc <- mc[
  as.character(mc$proportions) == "balanced" &
    DeCovarT:::.algorithm_in(
      mc$algorithm,
      c("Newton-Raphson", "Marquardt-Levenberg")
    ),
  ,
  drop = FALSE
]
z <- stats::qnorm(0.975)
truth_map <- artefacts$theta
p_lookup <- lapply(ids, function(id) {
  th <- .theta_for_id(id)
  stats::setNames(as.numeric(th$p), names(th$p))
})
names(p_lookup) <- ids
mc$p_true <- NA_real_
for (r in seq_len(nrow(mc))) {
  id <- as.character(mc$ID[[r]])
  ct <- as.character(mc$cell_type[[r]])
  pv <- p_lookup[[id]]
  if (!is.null(pv) && ct %in% names(pv)) {
    mc$p_true[[r]] <- unname(pv[[ct]])
  } else if (!is.null(pv)) {
    j <- match(ct, paste0("celltype_", seq_along(pv)))
    if (is.finite(j)) {
      mc$p_true[[r]] <- unname(pv[[j]])
    }
  }
}
mc$mean_est <- mc$p_true + mc$bias
mc$wald_lo <- mc$mean_est - z * mc$theoretical_se
mc$wald_hi <- mc$mean_est + z * mc$theoretical_se
mc$emp_lo <- mc$mean_est - z * mc$empirical_sd
mc$emp_hi <- mc$mean_est + z * mc$empirical_sd
drop_list <- "graph_generators"
forest <- mc[, setdiff(names(mc), drop_list), drop = FALSE]
forest$algorithm <- as.character(forest$algorithm)
forest$graph_ct1 <- as.character(forest$graph_ct1)
forest$graph_ct2 <- as.character(forest$graph_ct2)
forest$graph_ct3 <- as.character(forest$graph_ct3)
forest$overlap_label <- as.character(forest$overlap_label)
forest$proportions <- as.character(forest$proportions)

utils::write.csv(
  scenario_descriptors,
  file.path(OUT_DIR, "scenario_descriptors.csv"),
  row.names = FALSE
)
utils::write.csv(
  ilr_loglik_profiles,
  file.path(OUT_DIR, "ilr_loglik_profiles.csv"),
  row.names = FALSE
)
utils::write.csv(
  theoretical_wald_se,
  file.path(OUT_DIR, "theoretical_wald_se.csv"),
  row.names = FALSE
)
utils::write.csv(
  forest,
  file.path(OUT_DIR, "forest_newton_marquardt.csv"),
  row.names = FALSE
)

.ui_success("Wrote CSVs under {.file {OUT_DIR}}.")
invisible(NULL)
