test_that("Thomson projection returns two scores and a scree", {
  X <- withr::with_seed(1L, {
    MASS::mvrnorm(80, mu = c(0, 0, 0, 0, 1, 2), Sigma = diag(6))
  })
  proj <- DeCovarT:::.thomson_projection(X, n_factors = 2L)
  expect_equal(ncol(proj$scores), 2L)
  expect_equal(nrow(proj$scores), 80L)
  expect_true(length(proj$scree) <= 10L)
  expect_true(all(proj$scree >= 0))
})

test_that("labelled mixture draws match n and the simplex weights", {
  genes <- paste0("g", 1:4)
  cts <- paste0("ct", 1:3)
  mu <- matrix(
    c(0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0),
    nrow = 4L,
    dimnames = list(genes, cts)
  )
  sig <- diag(4)
  sigma <- array(c(sig, sig, sig), dim = c(4L, 4L, 3L))
  dimnames(sigma) <- list(genes, genes, cts)
  mix <- withr::with_seed(
    2L,
    DeCovarT:::.simulate_labelled_mixture(
      mu,
      sigma,
      p = c(0.5, 0.3, 0.2),
      n = 40L
    )
  )
  expect_equal(nrow(mix$Y), 40L)
  expect_length(mix$class, 40L)
  expect_true(all(mix$class %in% 1:3))
})

test_that("runtime book writes a PDF from optimisation rows", {
  skip_if_not_installed("ggdist")
  opt <- tibble::tibble(
    ID = "V1",
    sample_id = paste0("s", 1:30),
    algorithm = rep(c("lsei", "LBFGS", "cibersort"), each = 10),
    elapsed_sec = withr::with_seed(11L, runif(30, 0.01, 0.2)),
    graph_ct1 = "scale_free",
    graph_ct2 = "stochastic_block_model",
    overlap_label = rep(c("low", "moderate", "high"), each = 10),
    proportions = "balanced"
  )
  artefacts <- list(optimisation = opt)
  withr::with_tempfile("tf", fileext = ".pdf", {
    save_hybrid_runtime_book(artefacts, tf)
    expect_true(file.exists(tf))
    expect_gt(file.info(tf)$size, 0)
  })
})

test_that("bivariate runtime book writes a PDF from corner IDs", {
  skip_if_not_installed("ggdist")
  ids <- c("A", "B", "C", "D")
  opt <- tibble::tibble(
    ID = rep(ids, each = 8),
    sample_id = paste0("s", 1:32),
    algorithm = rep(c("lsei", "LBFGS"), 16),
    elapsed_sec = withr::with_seed(13L, runif(32, 0.01, 0.2))
  )
  cfg <- tibble::tibble(
    ID = ids,
    centroids = "small_CLD",
    variance = "homoscedastic",
    proportions = "balanced",
    correlation_celltype1 = c(0, -0.8, 0.8, -0.8),
    correlation_celltype2 = c(0, -0.8, 0.8, 0.8)
  )
  artefacts <- list(config = cfg, optimisation = opt)
  withr::with_tempfile("tf", fileext = ".pdf", {
    save_bivariate_runtime_book(artefacts, tf)
    expect_true(file.exists(tf))
    expect_gt(file.info(tf)$size, 0)
  })
})

test_that("memory book writes a PDF from optimisation rows", {
  skip_if_not_installed("ggdist")
  opt <- tibble::tibble(
    ID = "V1",
    sample_id = paste0("s", 1:30),
    algorithm = rep(c("lsei", "LBFGS", "cibersort"), each = 10),
    memory_bytes = withr::with_seed(12L, runif(30, 3.8e8, 4.2e8)),
    graph_ct1 = "scale_free",
    graph_ct2 = "scale_free",
    overlap_label = "low",
    proportions = "highly unbalanced"
  )
  artefacts <- list(optimisation = opt)
  withr::with_tempfile("tf", fileext = ".pdf", {
    save_hybrid_memory_book(artefacts, tf)
    expect_true(file.exists(tf))
    expect_gt(file.info(tf)$size, 0)
  })
})

test_that("latent projection book writes a one-page PDF", {
  skip_if_not_installed("EMMIXmfa")
  skip_if_not_installed("mclust")
  skip_if_not_installed("cowplot")
  genes <- paste0("g", seq_len(6L))
  cts <- paste0("celltype_", 1:3)
  mu <- matrix(
    c(8, 9, 9, 8, 2, 12, 2, 3, 11, 10, 4, 5, 5, 4, 12, 3, 6, 7),
    nrow = 6L,
    dimnames = list(genes, cts)
  )
  sig <- diag(6)
  sigma <- array(c(sig, sig, sig), dim = c(6L, 6L, 3L))
  dimnames(sigma) <- list(genes, genes, cts)
  th <- list(p = c(0.5, 0.3, 0.2), mu = mu, sigma = sigma)
  artefacts <- list(
    config = tibble::tibble(
      ID = "V1",
      proportions = "balanced",
      overlap_label = "low",
      graph_ct1 = "scale_free",
      graph_ct2 = "scale_free",
      graph_ct3 = "scale_free"
    ),
    theta = tibble::tibble(ID = "V1", true_theta = list(th))
  )
  withr::with_tempfile("tf", fileext = ".pdf", {
    save_hybrid_latent_projection_book(
      artefacts,
      tf,
      n = 40L,
      seed = 1L,
      itmax = 20L
    )
    expect_true(file.exists(tf))
    expect_gt(file.info(tf)$size, 0)
  })
})
