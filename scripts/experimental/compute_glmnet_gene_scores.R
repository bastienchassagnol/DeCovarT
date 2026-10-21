# EXPERIMENTAL -- not part of the DeCovarT package namespace.
#
# Gene scores from multinomial elastic-net cell-type classification
# (formerly R/04_04_statistical_metrics.R and a separate CV variant). They
# are not used by the current simulation or benchmark framework.
#
# Usage: source("scripts/experimental/compute_glmnet_gene_scores.R");
# requires the glmnet package.

# Fits a multinomial (binomial when J = 2) elastic net predicting the cell
# type from purified profiles X (G x J x N: genes x cell types x samples) and
# returns, per gene, the sum over classes of absolute coefficients (intercept
# excluded). `lambda = NULL` uses the smallest lambda of the path;
# `cross_validate = TRUE` selects lambda.min by cv.glmnet() with `nfolds`
# folds (default min(10, N)). `...` is forwarded to glmnet() / cv.glmnet().
compute_glmnet_gene_scores <- function(
  expression_profiles,
  celltype_labels,
  alpha = 0.5,
  lambda = NULL,
  cross_validate = FALSE,
  nfolds = NULL,
  ...
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("compute_glmnet_gene_scores() requires the glmnet package.")
  }
  dims <- dim(expression_profiles)
  if (
    length(dims) != 3L ||
      !is.numeric(expression_profiles) ||
      anyNA(expression_profiles)
  ) {
    stop("`expression_profiles` must be a numeric G x J x N array without NA.")
  }
  G <- dims[[1L]]
  J <- dims[[2L]]
  N <- dims[[3L]]
  if (J < 2L || length(celltype_labels) != J) {
    stop("Need J >= 2 cell types and one label per cell-type slice.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number in [0, 1].")
  }
  gene_names <- dimnames(expression_profiles)[[1L]]
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(G))
  }
  cell_names <- as.character(celltype_labels)

  # One row per (sample, cell type): (J * N) x G design, samples outermost.
  x <- t(matrix(expression_profiles, nrow = G, ncol = J * N))
  colnames(x) <- gene_names
  y <- factor(rep(cell_names, times = N), levels = unique(cell_names))

  family <- if (J == 2L) "binomial" else "multinomial"
  args <- list(x = x, y = y, family = family, alpha = alpha, ...)
  if (identical(family, "multinomial")) {
    args$type.multinomial <- "grouped"
  }
  if (isTRUE(cross_validate)) {
    args$nfolds <- if (is.null(nfolds)) min(10L, N) else as.integer(nfolds)
    if (args$nfolds < 2L || args$nfolds > nrow(x)) {
      stop("`nfolds` must satisfy 2 <= nfolds <= J * N.")
    }
    beta <- stats::coef(do.call(glmnet::cv.glmnet, args), s = "lambda.min")
  } else {
    fit <- do.call(glmnet::glmnet, args)
    beta <- stats::coef(
      fit,
      s = if (is.null(lambda)) min(fit$lambda) else lambda
    )
  }
  if (!is.list(beta)) {
    beta <- list(beta)
  }
  scores <- Reduce(
    `+`,
    lapply(beta, function(b) abs(as.numeric(as.matrix(b)[-1L, 1L])))
  )
  names(scores) <- gene_names
  scores
}
