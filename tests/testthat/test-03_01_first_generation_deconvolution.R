# Shared 2-gene / 2-celltype scenario for first-generation solvers.
.first_gen_setup <- function() {
  mean_signature_matrix <- matrix(
    c(20, 40, 40, 20),
    nrow = 2,
    dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2))
  )
  Sigma <- array(
    c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
    dim = c(2, 2, 2),
    dimnames = list(
      paste0("gene_", 1:2),
      paste0("gene_", 1:2),
      paste0("celltype_", 1:2)
    )
  )
  simulated_data <- withr::with_seed(
    3L,
    simulate_bulk_mixture(
      signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      p = c(0.5, 0.5),
      n = 1
    )
  )
  list(
    mean_signature_matrix = mean_signature_matrix,
    y = simulated_data$Y[, 1, drop = TRUE],
    p = c(0.5, 0.5)
  )
}

.expect_valid_simplex <- function(estimated_ratios, true_p, label) {
  testthat::expect_length(estimated_ratios, length(true_p))
  testthat::expect_false(anyNA(estimated_ratios), info = label)
  testthat::expect_true(
    all(estimated_ratios >= 0 & estimated_ratios <= 1),
    info = label
  )
  testthat::expect_equal(
    sum(estimated_ratios),
    1,
    tolerance = 1e-6,
    info = label
  )
}


test_that("First-generation deconvolution solvers return a valid simplex", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("limSolve")
  skip_if_not_installed("e1071")
  setup <- .first_gen_setup()

  # ------------------------------------------------------------------ #
  # Ordinary least squares (lsfit / Abbas-style)                       #
  # ------------------------------------------------------------------ #
  estimated_lsfit <- deconvolute_ratios_lsfit(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix
  )
  .expect_valid_simplex(estimated_lsfit, setup$p, "lsfit")

  # ------------------------------------------------------------------ #
  # Robust linear model (rlm / Monaco-style)                           #
  # ------------------------------------------------------------------ #
  estimated_rlm <- deconvolute_ratios_rlm(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix
  )
  .expect_valid_simplex(estimated_rlm, setup$p, "rlm")

  # ------------------------------------------------------------------ #
  # Non-negative least squares (nnls)                                  #
  # ------------------------------------------------------------------ #
  estimated_nnls <- deconvolute_ratios_nnls(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix
  )
  .expect_valid_simplex(estimated_nnls, setup$p, "nnls")

  # ------------------------------------------------------------------ #
  # Constrained least squares (lsei / deconRNASeq-style)               #
  # ------------------------------------------------------------------ #
  estimated_lsei <- deconvolute_ratios_deconrnaseq(
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix
  )
  .expect_valid_simplex(estimated_lsei, setup$p, "deconrnaseq")

  # ------------------------------------------------------------------ #
  # CIBERSORT-style nu-SVR (e1071); needs more genes than hold-out SVM #
  # ------------------------------------------------------------------ #
  set.seed(42L)
  n_genes <- 20L
  n_ct <- 2L
  cibersort_mu <- matrix(
    runif(n_genes * n_ct, min = 10, max = 50),
    nrow = n_genes,
    dimnames = list(
      paste0("gene_", seq_len(n_genes)),
      paste0("celltype_", seq_len(n_ct))
    )
  )
  true_p <- c(0.4, 0.6)
  cibersort_y <- as.numeric(cibersort_mu %*% true_p + rnorm(n_genes, sd = 1))
  estimated_cibersort <- deconvolute_ratios_cibersort(
    y = cibersort_y,
    mean_signature_matrix = cibersort_mu
  )
  .expect_valid_simplex(estimated_cibersort, true_p, "cibersort")
})
