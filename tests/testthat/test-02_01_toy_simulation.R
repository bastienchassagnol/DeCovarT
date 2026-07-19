##################################################################
##                   bivariate simulations                   ##
##################################################################

test_that("small simulation testing", {
  simulation_two_genes <- withr::with_seed(
    seed = 3L,
    simulate_bulk_mixture(
      signature_matrix = matrix(
        c(20, 40, 40, 20),
        nrow = 2,
        dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))
      ),
      Sigma = array(
        c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
        dim = c(2, 2, 2),
        dimnames = list(
          paste0("genes_", 1:2),
          paste0("genes_", 1:2),
          paste0("cell_type_", 1:2)
        )
      ),
      n = 10
    )
  )
  expect_equal(
    simulation_two_genes$Y[, 1:2],
    matrix(
      c(30.16652, 29.72223, 30.31822, 29.41423),
      nrow = 2,
      dimnames = list(paste0("genes_", 1:2), paste0("sample_", 1:2))
    ),
    tolerance = 10^-3
  )
})
