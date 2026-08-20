test_that("fit_decovart S3 accessors match the convolution (RE4, RE6)", {
  toy <- .toy_deconvolution()
  fit <- suppressWarnings(fit_decovart(
    signature_matrix = toy$signature_matrix,
    bulk_expression = toy$bulk_expression,
    true_ratios = toy$true_ratios,
    Sigma = toy$Sigma,
    method = "Marquardt-Levenberg",
    itmax = 60
  ))

  expect_s3_class(fit, "decovart_fit")
  expect_identical(rownames(coef(fit)), colnames(toy$signature_matrix))
  expect_identical(colnames(coef(fit)), colnames(toy$bulk_expression))
  expect_equal(
    fitted(fit),
    toy$signature_matrix %*% coef(fit),
    ignore_attr = TRUE,
    tolerance = .tol_srr
  )
  expect_equal(
    residuals(fit),
    toy$bulk_expression - fitted(fit),
    ignore_attr = TRUE,
    tolerance = .tol_srr
  )
  n <- nobs(fit)
  expect_identical(as.integer(n), ncol(toy$bulk_expression))
  expect_identical(attr(n, "n_genes"), nrow(toy$signature_matrix))
  expect_identical(attr(n, "n_celltypes"), ncol(toy$signature_matrix))
  expect_true(is.matrix(vcov(fit)) || is.list(vcov(fit)))
  expect_output(print(fit), "DeCovarT fit")
  expect_output(print(summary(fit)), "log-likelihood")
  grDevices::pdf(NULL)
  expect_silent(plot(fit))
  grDevices::dev.off()
  expect_false(is.function(utils::getS3method(
    "predict",
    "decovart_fit",
    optional = TRUE
  )))
  expect_false(is.function(utils::getS3method(
    "formula",
    "decovart_fit",
    optional = TRUE
  )))
})

test_that("gene-wise affine maps leave the MLE unchanged (RE2.3)", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- diag(c(1, 1.2, 0.8))
  Sigma[,, 2] <- diag(c(1.5, 0.9, 1.1))
  y <- matrix(
    drop(mu %*% c(0.65, 0.35)),
    ncol = 1,
    dimnames = list(genes, "s1")
  )

  p_raw <- suppressWarnings(deconvolute_ratios_Marquardt_Levenberg(
    drop(y),
    mu,
    Sigma,
    itmax = 80
  ))
  centre <- c(2, -1, 0.5)
  scale <- c(2, 0.5, 4)
  y_star <- (y - centre) / scale
  mu_star <- (mu - centre) / scale
  Sigma_star <- Sigma
  for (j in seq_len(2)) {
    Sigma_star[,, j] <- Sigma[,, j] * tcrossprod(1 / scale)
  }
  p_star <- suppressWarnings(deconvolute_ratios_Marquardt_Levenberg(
    drop(y_star),
    mu_star,
    Sigma_star,
    itmax = 80
  ))
  expect_equal(as.numeric(p_raw), as.numeric(p_star), tolerance = 1e-4)

  fit_z <- suppressWarnings(fit_decovart(
    mu,
    y,
    Sigma = Sigma,
    itmax = 80,
    standardise = TRUE
  ))
  fit_raw <- suppressWarnings(fit_decovart(
    mu,
    y,
    Sigma = Sigma,
    itmax = 80
  ))
  expect_equal(coef(fit_raw), coef(fit_z), tolerance = 1e-4)
})

test_that("log2 mixing is rejected; collinear signatures warn (RE2.2, RE2.4)", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- diag(3)
  Sigma[,, 2] <- diag(3)
  y <- matrix(drop(mu %*% c(0.5, 0.5)), ncol = 1, dimnames = list(genes, "s1"))

  expect_error(
    fit_decovart(mu, y, Sigma = Sigma, scaled = TRUE, itmax = 5),
    "log2"
  )

  mu_dup <- mu
  mu_dup[, 2] <- mu_dup[, 1]
  warn_dup <- character()
  withCallingHandlers(
    try(
      fit_decovart(mu_dup, y, Sigma = Sigma, itmax = 10),
      silent = TRUE
    ),
    warning = function(w) {
      warn_dup <<- c(warn_dup, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("collinear|identical", warn_dup)))

  cts3 <- c("child_1", "child_2", "parent")
  mu3 <- cbind(
    child_1 = c(20, 40, 15),
    child_2 = c(40, 20, 25),
    parent = 0.5 * c(20, 40, 15) + 0.5 * c(40, 20, 25)
  )
  rownames(mu3) <- genes
  Sigma3 <- array(0, dim = c(3, 3, 3), dimnames = list(genes, genes, cts3))
  for (j in 1:3) {
    Sigma3[,, j] <- diag(3)
  }
  y3 <- matrix(mu3 %*% c(0.2, 0.3, 0.5), ncol = 1)
  dimnames(y3) <- list(genes, "s1")
  expect_lt(qr(mu3)$rank, 3L)
  warn_parent <- character()
  withCallingHandlers(
    try(
      fit_decovart(mu3, y3, Sigma = Sigma3, itmax = 10),
      silent = TRUE
    ),
    warning = function(w) {
      warn_parent <<- c(warn_parent, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("collinear", warn_parent)))
})

test_that("Marquardt-Levenberg is stable under small bulk noise (RE5.9)", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- matrix(c(1, 0.2, 0, 0.2, 1.2, 0.1, 0, 0.1, 0.8), 3, 3)
  Sigma[,, 2] <- matrix(c(1.5, 0, 0.1, 0, 0.9, 0, 0.1, 0, 1.1), 3, 3)
  p_true <- c(0.6, 0.4)
  y0 <- matrix(drop(mu %*% p_true), ncol = 1, dimnames = list(genes, "s1"))

  fit_clean <- suppressWarnings(fit_decovart(
    mu,
    y0,
    Sigma = Sigma,
    method = "Marquardt-Levenberg",
    itmax = 80
  ))
  set.seed(2)
  y_pert <- y0 + rnorm(3, sd = 0.05)
  fit_pert <- suppressWarnings(fit_decovart(
    mu,
    y_pert,
    Sigma = Sigma,
    method = "Marquardt-Levenberg",
    itmax = 80
  ))
  expect_lt(
    sqrt(mean((drop(coef(fit_pert)) - drop(coef(fit_clean)))^2)),
    0.05
  )

  Sigma_diag <- Sigma
  for (j in 1:2) {
    Sigma_diag[,, j] <- diag(diag(Sigma[,, j]))
  }
  sigma_bar <- .compute_global_variance(p_true, Sigma)
  Sigma_global <- Sigma
  Sigma_global[,, 1] <- sigma_bar
  Sigma_global[,, 2] <- sigma_bar
  p_diag <- drop(coef(suppressWarnings(fit_decovart(
    mu,
    y_pert,
    Sigma = Sigma_diag,
    itmax = 80
  ))))
  p_global <- drop(coef(suppressWarnings(fit_decovart(
    mu,
    y_pert,
    Sigma = Sigma_global,
    itmax = 80
  ))))
  expect_equal(sum(p_diag), 1, tolerance = 1e-6)
  expect_equal(sum(p_global), 1, tolerance = 1e-6)
})

test_that("Marquardt and Newton recover p across seeds and starts (G5.6b, G5.9b)", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- matrix(c(1, 0.2, 0, 0.2, 1.2, 0.1, 0, 0.1, 0.8), 3, 3)
  Sigma[,, 2] <- matrix(c(1.5, 0, 0.1, 0, 0.9, 0, 0.1, 0, 1.1), 3, 3)
  p_true <- c(0.6, 0.4)
  solvers <- list(
    Marquardt = deconvolute_ratios_Marquardt_Levenberg,
    Newton = deconvolute_ratios_Newton_Raphson
  )
  seeds <- c(11L, 22L, 33L)

  dirichlet_start <- function(seed) {
    withr::with_seed(seed, {
      g <- stats::rgamma(length(p_true), shape = 1, rate = 1)
      g / sum(g)
    })
  }

  # G5.6b: three simulation seeds for the convolution draws.
  for (seed in seeds) {
    sim <- withr::with_seed(
      seed,
      simulate_bulk_mixture(mu, Sigma, p = p_true, n = 1L)
    )
    y <- sim$Y[, 1L, drop = TRUE]
    for (nm in names(solvers)) {
      p_hat <- suppressWarnings(solvers[[nm]](
        y,
        mu,
        Sigma,
        itmax = 80
      ))
      expect_equal(sum(p_hat), 1, tolerance = 1e-6, info = nm)
      expect_lt(
        sqrt(mean((as.numeric(p_hat) - p_true)^2)),
        0.2
      )
    }
  }

  # G5.9b: three random interior starts on a fixed perturbed bulk.
  y0 <- drop(mu %*% p_true)
  y_pert <- withr::with_seed(2L, y0 + stats::rnorm(3, sd = 0.05))
  for (nm in names(solvers)) {
    starts <- lapply(seeds, dirichlet_start)
    hats <- lapply(starts, function(p0) {
      suppressWarnings(solvers[[nm]](
        y_pert,
        mu,
        Sigma,
        itmax = 80,
        initial_p = p0
      ))
    })
    mat <- vapply(hats, as.numeric, numeric(2))
    expect_true(all(abs(colSums(mat) - 1) < 1e-6), info = nm)
    pairwise <- max(stats::dist(t(mat)))
    expect_lt(pairwise, 0.05)
  }
})

test_that("Diagonal and global weighted covariances stay on the simplex (RE7.1a)", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- matrix(c(1, 0.2, 0, 0.2, 1.2, 0.1, 0, 0.1, 0.8), 3, 3)
  Sigma[,, 2] <- matrix(c(1.5, 0, 0.1, 0, 0.9, 0, 0.1, 0, 1.1), 3, 3)
  p_true <- c(0.6, 0.4)
  y0 <- drop(mu %*% p_true)

  sigma_p <- p_true[[1L]]^2 * Sigma[,, 1] + p_true[[2L]]^2 * Sigma[,, 2]
  expect_equal(
    .compute_global_variance(p_true, Sigma),
    sigma_p,
    ignore_attr = TRUE,
    tolerance = 1e-10
  )

  Sigma_diag <- Sigma
  for (j in 1:2) {
    Sigma_diag[,, j] <- diag(diag(Sigma[,, j]))
  }
  Sigma_global <- Sigma
  Sigma_global[,, 1] <- sigma_p
  Sigma_global[,, 2] <- sigma_p

  specs <- list(
    full = Sigma,
    celltype_diagonal = Sigma_diag,
    global_weighted = Sigma_global
  )
  solvers <- list(
    Marquardt = deconvolute_ratios_Marquardt_Levenberg,
    Newton = deconvolute_ratios_Newton_Raphson
  )
  for (spec_nm in names(specs)) {
    for (sol_nm in names(solvers)) {
      p_hat <- suppressWarnings(solvers[[sol_nm]](
        y0,
        mu,
        specs[[spec_nm]],
        itmax = 80
      ))
      expect_equal(
        sum(p_hat),
        1,
        tolerance = 1e-6,
        info = paste(spec_nm, sol_nm)
      )
      expect_true(all(p_hat > 0), info = paste(spec_nm, sol_nm))
      expect_true(all(p_hat < 1), info = paste(spec_nm, sol_nm))
    }
  }
  p_full <- suppressWarnings(deconvolute_ratios_Marquardt_Levenberg(
    y0,
    mu,
    Sigma,
    itmax = 80
  ))
  expect_lt(sqrt(mean((as.numeric(p_full) - p_true)^2)), 5e-3)
})

test_that("L-BFGS-B and Newton-Raphson return models without repair_simplex", {
  genes <- paste0("g", 1:3)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 40, 25, 35, 30, 22),
    nrow = 3,
    dimnames = list(genes, cts)
  )
  Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
  Sigma[,, 1] <- diag(3)
  Sigma[,, 2] <- diag(3)
  y <- drop(mu %*% c(0.55, 0.45))

  for (method in c("L-BFGS-B", "Newton-Raphson")) {
    fit <- suppressWarnings(fit_decovart(
      mu,
      matrix(y, ncol = 1, dimnames = list(genes, "s1")),
      Sigma = Sigma,
      method = method,
      itmax = 40
    ))
    p <- drop(coef(fit))
    expect_equal(sum(p), 1, tolerance = 1e-6)
    expect_true(all(p >= 0 & p <= 1))
  }
})
