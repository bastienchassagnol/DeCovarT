# Extended tests (G5.10). Run with DECOVART_EXTENDED_TESTS=true.
#
# Regular suite timings (desktop): all files < 1s; slowest is
# test-02_02_generate_synthetic_networks.R (hierarchical GRN, 30 genes).
# This file repeats that generator at a larger G so it stays out of CI.

skip_if_not_extended()

test_that("hierarchical GRN moments scale to more genes", {
  moments <- withr::with_seed(1L, {
    simulate_hierarchical_grn_moments(
      n_genes = 60L,
      n_celltypes = 3L,
      mean_scale = 10,
      target_cosine = 0.1,
      precision_shift = 0.1,
      precision_scale = 0.3,
      prop_inhibitory = 0.5,
      graph_model = "scale_free"
    )
  })
  expect_identical(dim(moments$mean_profiles), c(60L, 3L))
  expect_identical(dim(moments$covariance_matrices), c(60L, 60L, 3L))
})
