# --------------------------------------------------------------------------
# Tests for simulate_hierarchical_grn_moments() and mean-profile objectives
# --------------------------------------------------------------------------

default_args <- list(
  n_genes = 30L,
  n_celltypes = 3L,
  mean_scale = 10,
  target_cosine = 0.1,
  precision_shift = 0.1,
  precision_scale = 0.3,
  prop_inhibitory = 0.5,
  graph_model = "scale_free",
  graph_params = list()
)

test_that("simulate_hierarchical_grn_moments returns expected structure", {
  skip_if_not_installed("igraph")
  moments <- withr::with_seed(1L, {
    do.call(simulate_hierarchical_grn_moments, default_args)
  })

  expect_named(
    moments,
    c(
      "mean_profiles",
      "covariance_matrices",
      "precision_matrices",
      "graph_structure",
      "objectives"
    )
  )
  expect_identical(dim(moments$mean_profiles), c(30L, 3L))
  expect_identical(dim(moments$covariance_matrices), c(30L, 30L, 3L))
  expect_identical(dim(moments$precision_matrices), c(30L, 30L, 3L))
  expect_identical(
    dim(moments$graph_structure$adjacency_matrices),
    c(30L, 30L, 3L)
  )
  expect_identical(
    dim(moments$graph_structure$weighted_adjacencies),
    c(30L, 30L, 3L)
  )
  expect_identical(
    dim(moments$graph_structure$normalised_precision),
    c(30L, 30L, 3L)
  )
  expect_true(all(moments$mean_profiles >= 0))
  expect_named(
    moments$objectives,
    c("mean_abs_cosine", "sum_euclidean_distance")
  )
})

test_that("per-cell-type covariances are PD and match precision inverses", {
  moments <- withr::with_seed(2L, {
    do.call(simulate_hierarchical_grn_moments, default_args)
  })

  for (j in seq_len(3L)) {
    sigma_j <- moments$covariance_matrices[,, j]
    omega_j <- moments$precision_matrices[,, j]
    n_genes <- nrow(sigma_j)
    expect_equal(
      unname(sigma_j %*% omega_j),
      diag(n_genes),
      tolerance = 1e-8
    )
    eigen_vals <- eigen(sigma_j, only.values = TRUE)$values
    expect_true(all(eigen_vals > 0))
  }
  # Independent draws: cell-type precisions need not coincide.
  expect_false(isTRUE(all.equal(
    moments$precision_matrices[,, 1],
    moments$precision_matrices[,, 2]
  )))
})

test_that("higher target_cosine increases mean absolute cosine", {
  low_cos <- withr::with_seed(3L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(target_cosine = 0))
    )
  })
  high_cos <- withr::with_seed(3L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(target_cosine = 0.95))
    )
  })

  expect_lt(
    low_cos$objectives$mean_abs_cosine,
    high_cos$objectives$mean_abs_cosine
  )
  expect_lt(low_cos$objectives$mean_abs_cosine, 0.5)
  expect_gt(high_cos$objectives$mean_abs_cosine, 0.8)
})

test_that("higher mean_scale increases Euclidean separation", {
  near_origin <- withr::with_seed(4L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(mean_scale = 1, target_cosine = 0))
    )
  })
  far_apart <- withr::with_seed(4L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(mean_scale = 20, target_cosine = 0))
    )
  })

  expect_lt(
    near_origin$objectives$sum_euclidean_distance,
    far_apart$objectives$sum_euclidean_distance
  )
  expect_equal(
    far_apart$objectives$sum_euclidean_distance /
      near_origin$objectives$sum_euclidean_distance,
    20,
    tolerance = 1e-8
  )
})

test_that("precision_scale alters off-diagonal magnitude of Omega", {
  weak <- withr::with_seed(6L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(precision_scale = 0.05))
    )
  })
  strong <- withr::with_seed(6L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(precision_scale = 0.8))
    )
  })

  off_weak <- weak$precision_matrices[,, 1]
  off_strong <- strong$precision_matrices[,, 1]
  diag(off_weak) <- 0
  diag(off_strong) <- 0
  expect_lt(max(abs(off_weak)), max(abs(off_strong)))
})

test_that("precision_shift improves conditioning of Omega", {
  ill <- withr::with_seed(7L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(
        default_args,
        list(precision_shift = 0.01, precision_scale = 0.5)
      )
    )
  })
  well <- withr::with_seed(7L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(
        default_args,
        list(precision_shift = 1, precision_scale = 0.5)
      )
    )
  })

  expect_gt(
    kappa(ill$precision_matrices[,, 1], exact = TRUE),
    kappa(well$precision_matrices[,, 1], exact = TRUE)
  )
})

test_that("stochastic block, small-world, erdos and hub yield PD precision", {
  for (model in c(
    "erdos_renyi",
    "hub",
    "stochastic_block_model",
    "small_world"
  )) {
    moments <- withr::with_seed(8L, {
      do.call(
        simulate_hierarchical_grn_moments,
        modifyList(default_args, list(graph_model = model))
      )
    })
    for (j in seq_len(3L)) {
      eigen_vals <- eigen(
        moments$precision_matrices[,, j],
        only.values = TRUE
      )$values
      expect_true(all(eigen_vals > 0))
    }
  }
})

test_that("prop_inhibitory controls sign balance on weighted edges", {
  moments <- withr::with_seed(9L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(prop_inhibitory = 1))
    )
  })
  for (j in seq_len(3L)) {
    W <- moments$graph_structure$weighted_adjacencies[,, j]
    A <- moments$graph_structure$adjacency_matrices[,, j]
    upper <- W[upper.tri(W) & A == 1]
    expect_true(length(upper) == 0L || all(upper > 0))
  }
})

test_that("pre-built adjacency list yields cell-type-specific supports", {
  set.seed(10L)
  a1 <- generate_random_network_skeleton(12L, "hub", list(n_hubs = 1L))
  a2 <- generate_random_network_skeleton(12L, "erdos_renyi", list())
  a3 <- generate_random_network_skeleton(12L, "scale_free", list())
  moments <- simulate_hierarchical_grn_moments(
    n_genes = 12L,
    n_celltypes = 3L,
    mean_scale = 10,
    target_cosine = 0.2,
    precision_shift = 0.1,
    precision_scale = 0.3,
    adjacency = list(a1, a2, a3)
  )
  expect_equal(
    moments$graph_structure$adjacency_matrices[,, 1],
    a1,
    ignore_attr = TRUE
  )
  expect_equal(
    moments$graph_structure$adjacency_matrices[,, 2],
    a2,
    ignore_attr = TRUE
  )
  expect_equal(
    moments$graph_structure$adjacency_matrices[,, 3],
    a3,
    ignore_attr = TRUE
  )
})

test_that("generate_mean_signature_matrix is deterministic and scaled", {
  mu_a <- generate_mean_signature_matrix(
    n_genes = 2L,
    n_celltypes = 2L,
    mean_scale = 10,
    target_cosine = 0.3
  )
  mu_b <- generate_mean_signature_matrix(
    n_genes = 2L,
    n_celltypes = 2L,
    mean_scale = 10,
    target_cosine = 0.3
  )

  expect_identical(dim(mu_a), c(2L, 2L))
  expect_identical(mu_a, mu_b)
  expect_equal(
    unname(sqrt(colSums(mu_a^2))),
    rep(10, 2),
    tolerance = 1e-8
  )
  cos_hat <- sum(mu_a[, 1] * mu_a[, 2]) /
    sqrt(sum(mu_a[, 1]^2) * sum(mu_a[, 2]^2))
  expect_gt(cos_hat, 0)
  expect_lt(cos_hat, 1)
})

test_that("graph_model matching is case-insensitive", {
  a_lower <- generate_random_network_skeleton(8L, "erdos_renyi")
  a_upper <- generate_random_network_skeleton(8L, "ERDOS_RENYI")
  expect_identical(dim(a_lower), dim(a_upper))
})
