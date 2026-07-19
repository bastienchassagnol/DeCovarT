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

# ---- output structure (power-law graph) ----------------------------------

test_that("output has the expected top-level structure", {
  res <- withr::with_seed(1L, run_default())

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
  n <- 5L
  res <- withr::with_seed(2L, run_default(n_expressed_genes = n))
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
  n <- 5L
  res <- withr::with_seed(3L, run_default(n_expressed_genes = n))
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
  n <- 5L
  res <- withr::with_seed(4L, run_default(n_expressed_genes = n))
  p <- 2L * n
  gene_names <- paste0("gene_", seq_len(p))

  adj <- res$graph_structure$adjacency_matrix
  expect_equal(dim(adj), c(p, p))
  expect_equal(rownames(adj), gene_names)
  expect_equal(colnames(adj), gene_names)

  prec <- res$graph_structure$normalised_precision
  expect_equal(dim(prec), c(p, p))
})


test_that("normalised precision is symmetric positive-definite", {
  res <- withr::with_seed(9L, run_default())
  prec <- res$graph_structure$normalised_precision
  expect_equal(prec, t(prec))
  eigenvalues <- eigen(prec, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})

test_that("all child covariance matrices are symmetric positive-definite", {
  res <- withr::with_seed(11L, run_default())
  ccov <- res$child_parameters$covariance_matrices
  for (k in seq_len(dim(ccov)[3])) {
    mat <- ccov[,, k]
    expect_equal(mat, t(mat))
    ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
})
