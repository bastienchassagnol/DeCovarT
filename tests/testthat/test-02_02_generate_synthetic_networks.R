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
  moments <- withr::with_seed(1L, {
    do.call(simulate_hierarchical_grn_moments, default_args)
  })

  expect_named(
    moments,
    c(
      "mean_profiles",
      "covariance_matrices",
      "graph_structure",
      "objectives"
    )
  )
  expect_equal(dim(moments$mean_profiles), c(30L, 3L))
  expect_equal(dim(moments$covariance_matrices), c(30L, 30L, 3L))
  expect_equal(
    dim(moments$graph_structure$adjacency_matrix),
    c(30L, 30L)
  )
  expect_equal(
    dim(moments$graph_structure$weighted_adjacency),
    c(30L, 30L)
  )
  expect_equal(
    dim(moments$graph_structure$normalised_precision),
    c(30L, 30L)
  )
  expect_true(all(moments$mean_profiles >= 0))
  expect_named(
    moments$objectives,
    c("mean_abs_cosine", "sum_euclidean_distance")
  )
})

test_that("shared covariances are PD and equal across cell types", {
  moments <- withr::with_seed(2L, {
    do.call(simulate_hierarchical_grn_moments, default_args)
  })

  sigma_1 <- moments$covariance_matrices[,, 1]
  sigma_2 <- moments$covariance_matrices[,, 2]
  expect_equal(sigma_1, sigma_2)
  expect_equal(
    sigma_1,
    solve(moments$graph_structure$normalised_precision),
    tolerance = 1e-8
  )
  eigen_vals <- eigen(sigma_1, only.values = TRUE)$values
  expect_true(all(eigen_vals > 0))
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

  off_weak <- weak$graph_structure$normalised_precision
  off_strong <- strong$graph_structure$normalised_precision
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
    kappa(ill$graph_structure$normalised_precision, exact = TRUE),
    kappa(well$graph_structure$normalised_precision, exact = TRUE)
  )
})

test_that("stochastic block and small-world models yield PD precision", {
  for (model in c("stochastic_block_model", "small_world")) {
    moments <- withr::with_seed(8L, {
      do.call(
        simulate_hierarchical_grn_moments,
        modifyList(default_args, list(graph_model = model))
      )
    })
    eigen_vals <- eigen(
      moments$graph_structure$normalised_precision,
      only.values = TRUE
    )$values
    expect_true(all(eigen_vals > 0))
  }
})

test_that("prop_inhibitory controls sign balance on weighted edges", {
  moments <- withr::with_seed(9L, {
    do.call(
      simulate_hierarchical_grn_moments,
      modifyList(default_args, list(prop_inhibitory = 1))
    )
  })
  W <- moments$graph_structure$weighted_adjacency
  upper <- W[upper.tri(W) & moments$graph_structure$adjacency_matrix == 1]
  expect_true(length(upper) == 0L || all(upper > 0))
})
