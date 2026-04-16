# --------------------------------------------------------------------------
# Tests for simulate_hierarchical_grn_moments()
# --------------------------------------------------------------------------

# Shared minimal valid arguments used across tests
default_args <- list(
  n_expressed_genes = 5L,
  mean_lower_expressed = 5,
  mean_upper_expressed = 20,
  mean_lower_background = 1,
  mean_upper_background = 3,
  library_size = 0,
  alpha = 2,
  precision_shift = 0.1,
  precision_scale = 0.3,
  child_perturbation_sd = 0.1,
  graph_model = "power_law",
  graph_params = list(power = 1, edges_per_node = 2)
)

run_default <- function(...) {
  args <- modifyList(default_args, list(...))
  do.call(simulate_hierarchical_grn_moments, args)
}

# ---- input validation ----------------------------------------------------

test_that("n_expressed_genes must be a single integer >= 2", {
  expect_error(
    run_default(n_expressed_genes = "a"),
    "must be a single integer >= 2"
  )
  expect_error(
    run_default(n_expressed_genes = 1L),
    "must be a single integer >= 2"
  )
  expect_error(
    run_default(n_expressed_genes = c(3L, 4L)),
    "must be a single integer >= 2"
  )
  expect_error(
    run_default(n_expressed_genes = -1),
    "must be a single integer >= 2"
  )
})

test_that("graph_model must be a valid choice", {
  expect_error(run_default(graph_model = "erdos_renyi"))
})

# ---- output structure (power-law graph) ----------------------------------

test_that("output has the expected top-level structure", {
  set.seed(1)
  res <- run_default()

  expect_type(res, "list")
  expect_named(
    res,
    c("parent_parameters", "child_parameters", "graph_structure")
  )
  expect_named(res$parent_parameters, c("mean_profiles", "covariance_matrices"))
  expect_named(res$child_parameters, c("mean_profiles", "covariance_matrices"))
  expect_named(
    res$graph_structure,
    c("adjacency_matrix", "normalised_precision")
  )
})

test_that("mean profiles have correct dimensions and names", {
  set.seed(2)
  n <- 5L
  res <- run_default(n_expressed_genes = n)
  p <- 2L * n
  gene_names <- paste0("gene_", seq_len(p))

  pm <- res$parent_parameters$mean_profiles
  expect_equal(dim(pm), c(2L, p))
  expect_equal(rownames(pm), c("parent_1", "parent_2"))
  expect_equal(colnames(pm), gene_names)

  cm <- res$child_parameters$mean_profiles
  expect_equal(dim(cm), c(4L, p))
  expect_equal(
    rownames(cm),
    c(
      "parent_1_child_a",
      "parent_1_child_b",
      "parent_2_child_a",
      "parent_2_child_b"
    )
  )
  expect_equal(colnames(cm), gene_names)
})

test_that("covariance arrays have correct dimensions and names", {
  set.seed(3)
  n <- 5L
  res <- run_default(n_expressed_genes = n)
  p <- 2L * n
  gene_names <- paste0("gene_", seq_len(p))

  pcov <- res$parent_parameters$covariance_matrices
  expect_equal(dim(pcov), c(p, p, 2L))
  expect_equal(dimnames(pcov)[[1]], gene_names)
  expect_equal(dimnames(pcov)[[2]], gene_names)
  expect_equal(dimnames(pcov)[[3]], c("parent_1", "parent_2"))

  ccov <- res$child_parameters$covariance_matrices
  expect_equal(dim(ccov), c(p, p, 4L))
  expect_equal(
    dimnames(ccov)[[3]],
    c(
      "parent_1_child_a",
      "parent_1_child_b",
      "parent_2_child_a",
      "parent_2_child_b"
    )
  )
})

test_that("adjacency and precision matrices have correct dimensions", {
  set.seed(4)
  n <- 5L
  res <- run_default(n_expressed_genes = n)
  p <- 2L * n
  gene_names <- paste0("gene_", seq_len(p))

  adj <- res$graph_structure$adjacency_matrix
  expect_equal(dim(adj), c(p, p))
  expect_equal(rownames(adj), gene_names)
  expect_equal(colnames(adj), gene_names)

  prec <- res$graph_structure$normalised_precision
  expect_equal(dim(prec), c(p, p))
})

# ---- mathematical properties --------------------------------------------

test_that("parent mean profiles are positive", {
  set.seed(5)
  res <- run_default()
  expect_true(all(res$parent_parameters$mean_profiles > 0))
})

test_that("child mean profiles are positive (clamped to machine eps)", {
  set.seed(6)
  res <- run_default()
  expect_true(all(res$child_parameters$mean_profiles > 0))
})

test_that("adjacency matrix is symmetric with zero diagonal", {
  set.seed(7)
  res <- run_default()
  adj <- res$graph_structure$adjacency_matrix
  expect_equal(adj, t(adj))
  expect_true(all(diag(adj) == 0))
})

test_that("adjacency matrix entries are binary", {
  set.seed(8)
  res <- run_default()
  adj <- res$graph_structure$adjacency_matrix
  expect_true(all(adj %in% c(0L, 1L)))
})

test_that("normalised precision is symmetric positive-definite", {
  set.seed(9)
  res <- run_default()
  prec <- res$graph_structure$normalised_precision
  expect_equal(prec, t(prec))
  eigenvalues <- eigen(prec, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("all parent covariance matrices are symmetric positive-definite", {
  set.seed(10)
  res <- run_default()
  pcov <- res$parent_parameters$covariance_matrices
  for (k in seq_len(dim(pcov)[3])) {
    mat <- pcov[,, k]
    expect_equal(mat, t(mat))
    ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
})

test_that("all child covariance matrices are symmetric positive-definite", {
  set.seed(11)
  res <- run_default()
  ccov <- res$child_parameters$covariance_matrices
  for (k in seq_len(dim(ccov)[3])) {
    mat <- ccov[,, k]
    expect_equal(mat, t(mat))
    ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
})

# ---- stochastic block model path ----------------------------------------

test_that("stochastic_block_model graph_model produces valid output", {
  set.seed(12)
  res <- run_default(
    graph_model = "stochastic_block_model",
    graph_params = list(
      block_prob = c(0.5, 0.25, 0.25),
      p_within = 0.3,
      p_between = 0.05
    )
  )

  expect_named(
    res,
    c("parent_parameters", "child_parameters", "graph_structure")
  )

  adj <- res$graph_structure$adjacency_matrix
  expect_equal(adj, t(adj))
  expect_true(all(diag(adj) == 0))
  expect_true(all(adj %in% c(0L, 1L)))

  prec <- res$graph_structure$normalised_precision
  ev <- eigen(prec, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

# ---- reproducibility -----------------------------------------------------

test_that("output is reproducible with set.seed", {
  set.seed(99)
  res1 <- run_default()
  set.seed(99)
  res2 <- run_default()

  expect_identical(
    res1$parent_parameters$mean_profiles,
    res2$parent_parameters$mean_profiles
  )
  expect_identical(
    res1$child_parameters$mean_profiles,
    res2$child_parameters$mean_profiles
  )
  expect_identical(
    res1$graph_structure$adjacency_matrix,
    res2$graph_structure$adjacency_matrix
  )
  expect_identical(
    res1$parent_parameters$covariance_matrices,
    res2$parent_parameters$covariance_matrices
  )
  expect_identical(
    res1$child_parameters$covariance_matrices,
    res2$child_parameters$covariance_matrices
  )
})

# ---- n_expressed_genes scaling -------------------------------------------

test_that("total gene count scales as 2 * n_expressed_genes", {
  set.seed(100)
  for (n in c(2L, 10L, 20L)) {
    res <- run_default(n_expressed_genes = n)
    p <- 2L * n
    expect_equal(ncol(res$parent_parameters$mean_profiles), p)
    expect_equal(nrow(res$graph_structure$adjacency_matrix), p)
    expect_equal(dim(res$parent_parameters$covariance_matrices)[1], p)
  }
})
