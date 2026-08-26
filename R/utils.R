#' Case-insensitive match to a set of allowed strings
#'
#' Like [match.arg()], but `arg` is matched after `tolower()` so `"Sigma"`
#' and `"sigma"` are equivalent. If `arg` is the full `choices` vector (the
#' usual default when the caller omits the argument), the first choice is
#' returned.
#'
#' @param arg Character vector supplied by the user.
#' @param choices Character vector of allowed values (canonical case).
#'
#' @return A character vector of canonical `choices` entries.
#'
#' @srrstats {G2.3} Univariate character arguments are restricted to a
#'   documented set of values.
#' @srrstats {G2.3a} Validation uses this `match.arg()` equivalent.
#' @srrstats {G2.3b} Matching is case-insensitive via `tolower()`.
#' @srrstats {G3.0} Numeric "equality" elsewhere uses
#'   `100 * .Machine$double.eps`, not `==` on floats.
#'
#' @keywords internal
#' @noRd
.match_arg_case_insensitive <- function(arg, choices) {
  if (!is.character(arg) || length(arg) < 1L) {
    stop("Argument must be a non-empty character vector.", call. = FALSE)
  }
  if (anyNA(arg)) {
    stop("Argument must not contain missing values.", call. = FALSE)
  }
  # match.arg() default: omitted arg is the full choices vector.
  if (
    length(arg) == length(choices) &&
      setequal(tolower(arg), tolower(choices))
  ) {
    arg <- arg[[1L]]
  }
  vapply(
    arg,
    function(a) {
      hit <- choices[tolower(choices) == tolower(a)]
      if (length(hit) != 1L) {
        stop(
          "Argument must be one of ",
          toString(choices),
          ".",
          call. = FALSE
        )
      }
      hit[[1L]]
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

#' Error if an object contains NA, NaN, Inf or -Inf
#'
#' Transcriptomic counts after quantification are complete matrices:
#' zeros (dropouts) are valid observations, but missing or non-finite
#' values are not. Proteomic intensities, by contrast, often contain
#' structurally missing entries (limit of detection, DDA sampling);
#' that case is out of scope for the current API.
#'
#' @param x Atomic vector, matrix, or array.
#' @param name Argument name used in the error message.
#'
#' @srrstats {G2.13} Missing values are rejected before analytic routines.
#' @srrstats {G2.14} The only supported missing-data policy is to error
#'   (no ignore / impute options).
#' @srrstats {G2.14a} `NA` / `NaN` raise an error.
#' @srrstats {G2.15} Downstream `mean` / `var` / `cor` therefore never see
#'   incomplete data (they keep default `na.rm = FALSE`).
#' @srrstats {G2.16} `Inf` and `-Inf` are rejected as well.
#'
#' @keywords internal
#' @noRd
.assert_no_missing <- function(x, name) {
  if (anyNA(x)) {
    stop(
      "`",
      name,
      "` must not contain missing or NaN values.",
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop("`", name, "` must not contain Inf or -Inf.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Require a numeric matrix or array with a single storage mode
#'
#' @param x Object to check.
#' @param name Argument name used in the error message.
#' @param n_dim Expected `length(dim(x))` (`2` matrix, `3` array).
#' @param non_negative If `TRUE`, every entry must be `>= 0`.
#'
#' @srrstats {G2.1} Type assertions: numeric matrix / array only.
#' @srrstats {G2.6} Vectors are rejected here; see `.align_true_ratios()`
#'   for the one-dimensional `true_ratios` case.
#' @srrstats {G2.8} Callers receive objects of one defined class.
#' @srrstats {G2.9} Class mismatches error instead of converting (no
#'   information is added or dropped silently).
#'
#' @keywords internal
#' @noRd
.assert_numeric_array <- function(x, name, n_dim, non_negative = FALSE) {
  if (is.factor(x)) {
    stop("`", name, "` must not be a factor.", call. = FALSE)
  }
  if (is.data.frame(x)) {
    stop(
      "`",
      name,
      "` must be a numeric matrix or array, not a data.frame.",
      call. = FALSE
    )
  }
  dims <- dim(x)
  if (is.null(dims) || length(dims) != n_dim) {
    stop(
      "`",
      name,
      "` must be a numeric ",
      if (n_dim == 2L) "matrix" else "array",
      " with ",
      n_dim,
      " dimensions.",
      call. = FALSE
    )
  }
  if (!is.numeric(x)) {
    stop("`", name, "` must have numeric storage mode.", call. = FALSE)
  }
  .assert_no_missing(x, name)
  if (non_negative && any(x < 0)) {
    stop("`", name, "` must be non-negative.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Append a file suffix when the caller omitted it (G4.0)
#'
#' @param path File path, with or without an extension.
#' @param suffix Canonical suffix (`"pdf"`, `"csv"`, or `"rds"`), with or
#'   without a leading dot.
#'
#' @return `path` with the required suffix.
#'
#' @srrstats {G4.0} Missing suffixes are added; a mismatched suffix errors.
#'
#' @keywords internal
#' @noRd
.ensure_file_suffix <- function(path, suffix) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty character string.", call. = FALSE)
  }
  suffix <- sub("^\\.", "", suffix)
  suffix <- .match_arg_case_insensitive(suffix, c("pdf", "csv", "rds"))
  ext <- tools::file_ext(path)
  if (!nzchar(ext)) {
    return(paste0(path, ".", suffix))
  }
  if (!identical(tolower(ext), suffix)) {
    stop(
      "File '",
      path,
      "' must use suffix .",
      suffix,
      ".",
      call. = FALSE
    )
  }
  path
}

#' Write an artefact with a canonical suffix
#'
#' Intended for tests and scripts. Wrap calls in
#' [withr::with_tempfile()] so temporary files are cleaned up.
#' Plots use PDF; R objects use RDS; rectangular tables use CSV.
#'
#' @param x Object to write.
#' @param path Destination path.
#' @param kind One of `"rds"`, `"csv"`, `"pdf"`.
#'
#' @return The path actually written (invisible).
#'
#' @srrstats {G4.0} Delegates suffix handling to `.ensure_file_suffix()`.
#'
#' @keywords internal
#' @noRd
.write_artefact <- function(x, path, kind = c("rds", "csv", "pdf")) {
  kind <- .match_arg_case_insensitive(kind, c("rds", "csv", "pdf"))
  path <- .ensure_file_suffix(path, kind)
  switch(
    kind,
    rds = saveRDS(x, path),
    csv = utils::write.csv(x, path, row.names = TRUE),
    pdf = {
      grDevices::pdf(path)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (inherits(x, "Heatmap")) {
        .check_heatmap_dependencies()
        ComplexHeatmap::draw(x)
      } else if (
        is.list(x) &&
          all(vapply(x, inherits, logical(1L), "Heatmap"))
      ) {
        .check_heatmap_dependencies()
        lapply(x, ComplexHeatmap::draw)
      } else {
        print(x)
      }
    }
  )
  invisible(path)
}

#' Require a cubic \eqn{G\times G\times J} numeric array
#'
#' Shared by `.prepare_deconvolution_inputs()` (cell-type `Sigma`) and
#' `.parse_true_theta()` (generative `sigma` / `Theta`). Full merge of
#' those two gates is intentionally avoided: one validates deconvolution
#' inputs (\eqn{\boldsymbol{\mu}}, \eqn{\boldsymbol{Y}}, \eqn{\boldsymbol{\Sigma}}),
#' the other validates MixSim / Jeffreys generative parameters
#' \eqn{\theta=(p,\boldsymbol{\mu},\boldsymbol{\Sigma})}.
#'
#' @param arr Array to check.
#' @param name Argument name used in the error message.
#'
#' @return Named list with `G` and `J`.
#'
#' @keywords internal
#' @noRd
.assert_ggj_array <- function(arr, name) {
  if (is.null(dim(arr)) || length(dim(arr)) != 3L) {
    stop("`", name, "` must be a G x G x J array.", call. = FALSE)
  }
  g1 <- dim(arr)[[1L]]
  g2 <- dim(arr)[[2L]]
  jj <- dim(arr)[[3L]]
  if (g1 != g2) {
    stop("`", name, "` dims must be G x G x J.", call. = FALSE)
  }
  list(G = g1, J = jj)
}

#' Align optional ground-truth proportions to \eqn{J \times N}
#'
#' @param true_ratios `NULL`, length-\eqn{J} numeric vector, or numeric
#'   matrix with dimensions \eqn{J\times N} only (no \eqn{N\times J}
#'   transpose).
#' @param n_celltypes Number of cell types \eqn{J} (rows of the returned
#'   matrix; must match `ncol(signature_matrix)`).
#' @param n_samples Number of bulk samples \eqn{N} to deconvolve
#'   (columns of the returned matrix; must match
#'   `ncol(bulk_expression)`).
#'
#' @return `NULL` or a numeric matrix \eqn{J\times N}.
#'
#' @srrstats {G2.6} One-dimensional numeric vectors are accepted and
#'   recycled across samples; other 1-d classes error.
#'
#' @keywords internal
#' @noRd
.align_true_ratios <- function(true_ratios, n_celltypes, n_samples) {
  if (is.null(true_ratios)) {
    return(NULL)
  }
  if (is.factor(true_ratios) || is.data.frame(true_ratios)) {
    stop(
      "`true_ratios` must be a numeric vector or matrix.",
      call. = FALSE
    )
  }
  if (is.matrix(true_ratios)) {
    .assert_numeric_array(true_ratios, "true_ratios", 2L, non_negative = TRUE)
    if (nrow(true_ratios) == n_celltypes && ncol(true_ratios) == n_samples) {
      return(true_ratios)
    }
    stop(
      "`true_ratios` as a matrix must be J x N ",
      "(n_celltypes x n_samples).",
      call. = FALSE
    )
  }
  if (!is.numeric(true_ratios) || length(dim(true_ratios))) {
    stop(
      "`true_ratios` must be a numeric vector of length J or a ",
      "J x N matrix.",
      call. = FALSE
    )
  }
  .assert_no_missing(true_ratios, "true_ratios")
  if (any(true_ratios < 0)) {
    stop("`true_ratios` must be non-negative.", call. = FALSE)
  }
  if (length(true_ratios) != n_celltypes) {
    stop("`true_ratios` must have length J.", call. = FALSE)
  }
  matrix(
    true_ratios,
    nrow = n_celltypes,
    ncol = n_samples
  )
}

#' Validate and align inputs for [deconvolute_ratios()]
#'
#' Single pre-processing gate (G2.8): solvers downstream receive numeric
#' matrices / arrays only. Genes are matched by row names with
#' `drop = FALSE` (G2.10). When \eqn{J > G} the linear mixture is
#' undetermined and an error is raised.
#'
#' Complements [check_true_theta()] / `.parse_true_theta()`: those
#' validators check generative-model lists
#' \eqn{\theta=(p,\boldsymbol{\mu},\boldsymbol{\Sigma})} for MixSim /
#' Jeffreys metrics, while this gate checks the deconvolution call
#' signature (\eqn{\boldsymbol{\mu}}, \eqn{\boldsymbol{Y}}, optional
#' ground-truth ratios, \eqn{\boldsymbol{\Sigma}}). Shared cubic-array
#' checks use `.assert_ggj_array()`.
#'
#' @return Named list with aligned `signature_matrix`, `bulk_expression`,
#'   `true_ratios` (\eqn{J\times N} or `NULL`), `Sigma`, and optional
#'   gene-wise affine `centre` / `scale` when `standardise = TRUE`.
#'
#' @srrstats {G2.8} This is the unique conversion / dispatch point.
#' @srrstats {G2.10} Matrix subsetting uses `drop = FALSE`.
#' @srrstats {G5.8d} \eqn{J > G} (more cell types than genes) is rejected.
#' @srrstats {RE2.3} Gene-wise affine standardisation is the only supported
#'   scaling; log2 mixing (`scaled = TRUE`) is rejected.
#' @srrstats {RE2.4a} Rank-deficient or duplicate signature columns warn.
#'
#' @keywords internal
#' @noRd
.prepare_deconvolution_inputs <- function(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  standardise = FALSE,
  scaled = FALSE
) {
  .assert_numeric_array(
    signature_matrix,
    "signature_matrix",
    2L,
    non_negative = TRUE
  )
  .assert_numeric_array(
    bulk_expression,
    "bulk_expression",
    2L,
    non_negative = TRUE
  )
  if (
    is.null(rownames(signature_matrix)) ||
      is.null(rownames(bulk_expression))
  ) {
    stop(
      "`signature_matrix` and `bulk_expression` must have gene rownames.",
      call. = FALSE
    )
  }
  if (is.null(colnames(signature_matrix))) {
    stop("`signature_matrix` must have cell-type colnames.", call. = FALSE)
  }

  common_genes <- intersect(
    rownames(signature_matrix),
    rownames(bulk_expression)
  )
  if (length(common_genes) / nrow(signature_matrix) < 0.5) {
    stop(
      "Only ",
      length(common_genes) / nrow(signature_matrix),
      " fraction of genes are used in the signature matrix.\n",
      "Half of common genes are required at least",
      call. = FALSE
    )
  }
  common_genes <- sort(common_genes)
  signature_matrix <- signature_matrix[common_genes, , drop = FALSE]
  bulk_expression <- bulk_expression[common_genes, , drop = FALSE]

  n_genes <- nrow(signature_matrix)
  n_celltypes <- ncol(signature_matrix)
  n_samples <- ncol(bulk_expression)
  # J > G: undetermined linear mixture (arxiv:2310.14722, sec. 2.1.5).
  if (n_celltypes > n_genes) {
    stop(
      "Undetermined deconvolution: number of cell types J = ",
      n_celltypes,
      " exceeds number of genes G = ",
      n_genes,
      ".",
      call. = FALSE
    )
  }

  if (!is.null(Sigma)) {
    .assert_numeric_array(Sigma, "Sigma", 3L, non_negative = FALSE)
    dims_sigma <- .assert_ggj_array(Sigma, "Sigma")
    if (is.null(dimnames(Sigma)[[1L]]) || is.null(dimnames(Sigma)[[3L]])) {
      if (
        !identical(
          c(dims_sigma$G, dims_sigma$G, dims_sigma$J),
          c(n_genes, n_genes, n_celltypes)
        )
      ) {
        stop(
          "`Sigma` must be a G x G x J array matching the aligned ",
          "signature matrix.",
          call. = FALSE
        )
      }
    } else {
      missing_genes <- setdiff(common_genes, dimnames(Sigma)[[1L]])
      if (length(missing_genes)) {
        stop(
          "`Sigma` is missing genes present in the signature.",
          call. = FALSE
        )
      }
      Sigma <- Sigma[common_genes, common_genes, , drop = FALSE]
      if (dim(Sigma)[[3L]] != n_celltypes) {
        stop(
          "`Sigma` must have J slices (one covariance per cell type).",
          call. = FALSE
        )
      }
    }
  }

  true_ratios <- .align_true_ratios(
    true_ratios,
    n_celltypes = n_celltypes,
    n_samples = n_samples
  )

  .warn_collinear_signature(signature_matrix)

  if (isTRUE(scaled)) {
    stop(
      "`scaled = TRUE` (log2 mixing) is not supported: the logarithm ",
      "is a nonlinear map and, by Jensen's inequality, changes first ",
      "and second moments, breaking the linear convolution ",
      "y = mu p. CIBERSORT likewise requires non-negative expression, ",
      "no missing values, and a non-log linear scale ",
      "(Newman et al., 2015). Use `standardise = TRUE` for a gene-wise ",
      "affine z-score that leaves the theoretical MLE unchanged.",
      call. = FALSE
    )
  }

  centre <- NULL
  scale <- NULL
  if (isTRUE(standardise)) {
    std <- .standardise_gene_wise(
      signature_matrix,
      bulk_expression,
      Sigma
    )
    signature_matrix <- std$signature_matrix
    bulk_expression <- std$bulk_expression
    Sigma <- std$Sigma
    centre <- std$centre
    scale <- std$scale
  }

  list(
    signature_matrix = signature_matrix,
    bulk_expression = bulk_expression,
    true_ratios = true_ratios,
    Sigma = Sigma,
    centre = centre,
    scale = scale
  )
}

#' Warn when signature columns are collinear (RE2.4a, RE2.4b)
#'
#' @keywords internal
#' @noRd
.warn_collinear_signature <- function(signature_matrix) {
  n_celltypes <- ncol(signature_matrix)
  rank_mu <- qr(signature_matrix)$rank
  if (rank_mu < n_celltypes) {
    warning(
      "Signature columns are collinear (rank ",
      rank_mu,
      " < J = ",
      n_celltypes,
      "); mixture proportions are not identifiable.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  col_norm <- sqrt(colSums(signature_matrix^2))
  if (any(col_norm < .Machine$double.eps)) {
    warning(
      "A signature column is numerically zero.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  unit_cols <- sweep(signature_matrix, 2L, col_norm, "/")
  cosine <- abs(crossprod(unit_cols))
  diag(cosine) <- 0
  if (any(cosine > 1 - 1e-8)) {
    warning(
      "At least two signature columns are identical up to scaling.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Gene-wise affine z-score of bulk, means and covariances (RE2.3)
#'
#' Centre and scale are computed once from the reference means
#' \eqn{\boldsymbol{\mu}} (not per cell type or per bulk sample) and applied
#' as \(D^{-1}(\boldsymbol{x}-\boldsymbol{m})\) to \(\boldsymbol{Y}\) and
#' \(\boldsymbol{\mu}\), and as \(D^{-1}\boldsymbol{\Sigma}_j D^{-1}\) to each
#' covariance slice.
#'
#' @keywords internal
#' @noRd
.standardise_gene_wise <- function(
  signature_matrix,
  bulk_expression,
  Sigma = NULL
) {
  centre <- rowMeans(signature_matrix)
  scale <- apply(signature_matrix, 1L, stats::sd)
  if (any(!is.finite(scale) | scale < .Machine$double.eps)) {
    stop(
      "Gene-wise standardisation requires a positive finite standard ",
      "deviation for every gene in the signature.",
      call. = FALSE
    )
  }
  signature_star <- (signature_matrix - centre) / scale
  bulk_star <- (bulk_expression - centre) / scale
  sigma_star <- Sigma
  if (!is.null(Sigma)) {
    inv_scale <- 1 / scale
    scale_outer <- tcrossprod(inv_scale)
    for (j in seq_len(dim(Sigma)[[3L]])) {
      sigma_star[,, j] <- Sigma[,, j] * scale_outer
    }
  }
  list(
    signature_matrix = signature_star,
    bulk_expression = bulk_star,
    Sigma = sigma_star,
    centre = centre,
    scale = scale
  )
}
