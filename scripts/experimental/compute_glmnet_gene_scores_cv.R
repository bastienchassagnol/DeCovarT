# Experimental: CV fold family for glmnet gene scores (not in package build).
# Source locally if needed: source("R/experimental/compute_glmnet_gene_scores_cv.R")

#' Gene scores via multinomial elastic net with CV lambda
#'
#' Cross-validated companion to [DeCovarT::compute_glmnet_gene_scores()].
#' Fits [glmnet::cv.glmnet()] on purified profiles
#' \eqn{G\times J\times N} and extracts absolute coefficient sums at
#' `lambda.min`.
#'
#' @inheritParams DeCovarT::compute_glmnet_gene_scores
#' @param nfolds Number of CV folds (default `min(10L, N)`).
#' @param ... Additional arguments forwarded to [glmnet::cv.glmnet()].
#'
#' @return Named numeric vector of length \eqn{G}.
#'
#' @keywords internal
compute_glmnet_gene_scores_cv <- function(
  expression_profiles,
  celltype_labels,
  alpha = 0.5,
  nfolds = NULL,
  ...
) {
  if (
    is.null(dim(expression_profiles)) ||
      length(dim(expression_profiles)) != 3L
  ) {
    stop("`expression_profiles` must be a G x J x N array.")
  }
  if (!is.numeric(expression_profiles) || anyNA(expression_profiles)) {
    stop("`expression_profiles` must be numeric without missing values.")
  }
  G <- dim(expression_profiles)[[1L]]
  J <- dim(expression_profiles)[[2L]]
  N <- dim(expression_profiles)[[3L]]
  if (G < 1L || J < 2L || N < 1L) {
    stop("`expression_profiles` must be G x J x N with G >= 1, J >= 2, N >= 1.")
  }
  if (length(celltype_labels) != J) {
    stop("`celltype_labels` must have length J (second dim of profiles).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number in [0, 1].")
  }

  gene_names <- dimnames(expression_profiles)[[1L]]
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(G))
  }
  cell_names <- as.character(celltype_labels)

  x <- matrix(NA_real_, nrow = J * N, ncol = G)
  y <- factor(rep(cell_names, times = N), levels = unique(cell_names))
  row_id <- 1L
  for (n in seq_len(N)) {
    for (j in seq_len(J)) {
      x[row_id, ] <- expression_profiles[, j, n]
      row_id <- row_id + 1L
    }
  }
  colnames(x) <- gene_names

  if (is.null(nfolds)) {
    nfolds <- min(10L, N)
  }
  nfolds <- as.integer(nfolds)
  if (nfolds < 2L || nfolds > nrow(x)) {
    stop("`nfolds` must satisfy 2 <= nfolds <= J * N.")
  }

  family <- if (J == 2L) "binomial" else "multinomial"
  cv_args <- list(
    x = x,
    y = y,
    family = family,
    alpha = alpha,
    nfolds = nfolds
  )
  if (identical(family, "multinomial")) {
    cv_args$type.multinomial <- "grouped"
  }
  dots <- list(...)
  for (nm in names(dots)) {
    cv_args[[nm]] <- dots[[nm]]
  }

  fit <- do.call(glmnet::cv.glmnet, cv_args)
  beta <- stats::coef(fit, s = "lambda.min")

  scores <- numeric(G)
  names(scores) <- gene_names
  if (identical(family, "binomial")) {
    cj <- as.matrix(beta)[-1L, 1L]
    scores <- abs(as.numeric(cj))
    names(scores) <- gene_names
  } else {
    for (j in seq_along(beta)) {
      cj <- as.matrix(beta[[j]])[-1L, 1L]
      scores <- scores + abs(as.numeric(cj))
    }
  }
  scores
}
