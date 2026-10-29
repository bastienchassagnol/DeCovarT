# Second-generation competitors for the fig03 and fig05 benchmarks.
#
# Every method runs from its authors' code, with gene selection and
# preprocessing filters switched off (the scenarios are fully synthetic) and
# every other default kept.
#
# * MEAD (Xie, Gui and Wang): R package `MEAD` 1.0.2 from
#   https://github.com/DongyueXie/MEAD (commit d66fa10). Its NAMESPACE imports
#   Matrix, Rfast, SingleCellExperiment, cluster and vashr
#   (github::mengyin/vashr) without declaring them; install those first.
#   Called through `MEAD::MEAD()` with `preprocessing_control =
#   list(filter.gene = FALSE)`, so `MEAD_preprocessing()` only aligns genes,
#   then the official `MEAD_getX()` (per-cell unit-sum normalisation, mean
#   and between-individual covariance, vash weights) and `MEAD_est()` with
#   its defaults (hc.type = "hc3", use.QP = FALSE, truncation at zero).
#   R01 comes from `MEAD::get_R01()` on the scenario cohort of bulks
#   (alpha = 0.1, as in the package vignette); when it fails, R01 = NULL
#   (independent genes) and the route is recorded in `status`.
#   Intervals: `MEAD::get_ci()`. Because `MEAD_getX()` divides each reference
#   cell by its total, MEAD's p_hat = beta / sum(beta) weights each type by
#   its mean transcript total: it equals the cell fraction only when the
#   reference types have equal totals.
# * RNA-Sieve (Erdmann-Pham et al.): PyPI `rnasieve` 0.1.4
#   (https://github.com/songlab-cal/rna-sieve), run by `python/run_rnasieve.py`:
#   `model_from_raw_counts(trim_percent = 0, gene_thresh = 0,
#   normalization = False)`, joint `predict()` over the scenario cohort and
#   `compute_marginal_confidence_intervals()` (Godambe sandwich). When the
#   averaged Hessian has condition number >= 1e6 (stats.py:234), rnasieve
#   0.1.4 and the GitHub master subsample genes and return H J H' instead of
#   H^-1 J H^-T (stats.py:224), so its intervals collapse to [0, 1] (seen
#   at G = 20 on data drawn from RNA-Sieve's own model; when the branch does
#   not fire, the released and recomputed Godambe intervals coincide). The
#   runner therefore also returns the paper's Fisher (Wald) and Godambe
#   intervals, computed with rnasieve's own score and Hessian functions on
#   all genes (columns lower_fisher, upper_fisher, lower_godambe,
#   upper_godambe).
# * ReDeconv (Lu et al.): PyPI `redeconv` 1.2.1
#   (https://github.com/jyyulab/redeconv), run by `python/run_redeconv.py`:
#   the default `ReDeconv()` entry (MLMD, six reweighting loops) on a
#   signature written from the plug-in moments (true means and the square
#   roots of the diagonals of the true covariances), which bypasses its
#   signature-gene selection. No interval.
# * DECALS (Cai et al.): no public implementation, so it is not run. Its
#   point estimate, constrained OLS on the simplex solved as a quadratic
#   programme, is the `lsei` (DeconRNASeq) baseline already benchmarked; its
#   asymptotic interval is not evaluated.
#
# The two Python runners share the uv project in `python/` (Python 3.11,
# pinned by `uv.lock`). R calls them once per method for all scenarios, and
# exchanges CSV files through a temporary directory removed on exit.

.second_generation_methods <- c("MEAD", "RNA-Sieve", "ReDeconv")

#' Run MEAD, RNA-Sieve and ReDeconv on a batch of scenarios
#'
#' @param scenarios List of scenarios. Each is a list with `scenario_id`,
#'   `bulk` (G x N matrix, gene row names), `reference` (named list, one
#'   G x n_cells matrix per cell type), `individual` (named list, the
#'   individual of each reference column), `mean_signature_matrix` (G x J)
#'   and `Sigma` (G x G x J). The two last ones are the ReDeconv plug-in.
#' @param methods Subset of `c("MEAD", "RNA-Sieve", "ReDeconv")`.
#' @param level Confidence level of the intervals.
#' @param python_dir Directory of the uv project holding the runners.
#' @param workers Processes for RNA-Sieve and ReDeconv (one scenario per
#'   process at a time) and `future` workers for MEAD; the scenarios are
#'   independent fits.
#' @return A tibble with `scenario_id`, `algorithm`, `replicate`,
#'   `celltype`, `estimate`, `lower`, `upper`, `elapsed_sec` (scenario time
#'   divided by its number of bulk samples) and `status`; RNA-Sieve rows add
#'   `lower_fisher`, `upper_fisher`, `lower_godambe` and `upper_godambe`.
#' @noRd
run_second_generation_batch <- function(
  scenarios,
  methods = .second_generation_methods,
  level = 0.95,
  python_dir = file.path("scripts", "benchmark_metrics", "python"),
  workers = 1L
) {
  methods <- match.arg(methods, .second_generation_methods, several.ok = TRUE)
  stopifnot(is.numeric(level), length(level) == 1L, level > 0, level < 1)
  scenarios <- lapply(scenarios, .check_competitor_scenario)
  ids <- vapply(scenarios, `[[`, character(1), "scenario_id")
  if (anyDuplicated(ids)) {
    stop("`scenario_id` must be unique across scenarios.", call. = FALSE)
  }

  python_methods <- intersect(methods, c("RNA-Sieve", "ReDeconv"))
  results <- list()
  if (length(python_methods) > 0L) {
    work_dir <- withr::local_tempdir(pattern = "decovart_competitors_")
    in_dir <- file.path(work_dir, "input")
    .write_competitor_inputs(scenarios, in_dir)
    for (method in python_methods) {
      started <- proc.time()[["elapsed"]]
      results[[method]] <- .run_python_competitor(
        method,
        scenarios,
        in_dir,
        file.path(work_dir, method),
        level,
        python_dir,
        workers
      )
      .log_competitor(method, results[[method]], started)
    }
  }
  if ("MEAD" %in% methods) {
    started <- proc.time()[["elapsed"]]
    cli::cli_alert_info(
      "MEAD: {length(scenarios)} scenarios on {workers} worker(s)."
    )
    results[["MEAD"]] <- dplyr::bind_rows(furrr::future_map(
      scenarios,
      .run_mead_scenario,
      level = level,
      .options = furrr::furrr_options(seed = TRUE)
    ))
    .log_competitor("MEAD", results[["MEAD"]], started)
  }
  dplyr::bind_rows(results[intersect(.second_generation_methods, methods)])
}

.check_competitor_scenario <- function(scenario) {
  required <- c(
    "scenario_id",
    "bulk",
    "reference",
    "individual",
    "mean_signature_matrix",
    "Sigma"
  )
  missing <- setdiff(required, names(scenario))
  if (length(missing) > 0L) {
    stop(
      "Each scenario needs ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  scenario$scenario_id <- as.character(scenario$scenario_id)
  bulk <- as.matrix(scenario$bulk)
  if (is.null(rownames(bulk))) {
    rownames(bulk) <- paste0("gene_", seq_len(nrow(bulk)))
  }
  colnames(bulk) <- paste0("bulk_", seq_len(ncol(bulk)))
  mu <- as.matrix(scenario$mean_signature_matrix)
  celltypes <- colnames(mu)
  if (is.null(celltypes)) {
    stop("`mean_signature_matrix` needs cell-type column names.", call. = FALSE)
  }
  if (!setequal(names(scenario$reference), celltypes)) {
    stop(
      "`reference` must be named by the cell types of ",
      "`mean_signature_matrix`.",
      call. = FALSE
    )
  }
  if (nrow(mu) != nrow(bulk) || dim(scenario$Sigma)[1L] != nrow(bulk)) {
    stop("Bulk, means and covariances must share the genes.", call. = FALSE)
  }
  rownames(mu) <- rownames(bulk)
  scenario$reference <- lapply(scenario$reference[celltypes], function(x) {
    x <- as.matrix(x)
    if (nrow(x) != nrow(bulk)) {
      stop("Every reference matrix must have one row per gene.", call. = FALSE)
    }
    rownames(x) <- rownames(bulk)
    x
  })
  scenario$individual <- scenario$individual[celltypes]
  scenario$bulk <- bulk
  scenario$mean_signature_matrix <- mu
  scenario
}

.write_competitor_inputs <- function(scenarios, in_dir) {
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  ids <- vapply(scenarios, `[[`, character(1), "scenario_id")
  utils::write.csv(
    data.frame(scenario_id = ids),
    file.path(in_dir, "manifest.csv"),
    row.names = FALSE
  )
  for (scenario in scenarios) {
    scenario_dir <- file.path(in_dir, scenario$scenario_id)
    dir.create(scenario_dir, showWarnings = FALSE)
    celltypes <- colnames(scenario$mean_signature_matrix)
    .write_gene_table(scenario$bulk, file.path(scenario_dir, "bulk.csv"))
    utils::write.csv(
      data.frame(celltype = celltypes),
      file.path(scenario_dir, "celltypes.csv"),
      row.names = FALSE
    )
    for (k in seq_along(celltypes)) {
      .write_gene_table(
        scenario$reference[[k]],
        file.path(scenario_dir, sprintf("reference_%d.csv", k))
      )
    }
    sd <- sqrt(apply(scenario$Sigma, 3L, diag))
    dimnames(sd) <- dimnames(scenario$mean_signature_matrix)
    .write_gene_table(
      scenario$mean_signature_matrix,
      file.path(scenario_dir, "signature_mean.csv")
    )
    .write_gene_table(sd, file.path(scenario_dir, "signature_sd.csv"))
  }
  invisible(in_dir)
}

.write_gene_table <- function(x, path) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("column_", seq_len(ncol(x)))
  }
  out <- data.frame(gene = rownames(x), x, check.names = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
}

.run_python_competitor <- function(
  method,
  scenarios,
  in_dir,
  out_dir,
  level,
  python_dir,
  workers = 1L
) {
  uv <- Sys.which("uv")
  if (!nzchar(uv)) {
    stop(
      "`uv` is not on the PATH; install it to run ",
      method,
      " (https://docs.astral.sh/uv/).",
      call. = FALSE
    )
  }
  script <- switch(
    method,
    "RNA-Sieve" = "run_rnasieve.py",
    "ReDeconv" = "run_redeconv.py"
  )
  python_dir <- normalizePath(python_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "runner.log")
  args <- c(
    "run",
    "--frozen",
    "--quiet",
    "--project",
    shQuote(python_dir),
    "python",
    shQuote(file.path(python_dir, script)),
    shQuote(in_dir),
    shQuote(out_dir)
  )
  if (identical(method, "RNA-Sieve")) {
    args <- c(args, "--level", format(level, digits = 15))
  }
  args <- c(args, "--workers", as.integer(workers))
  # One line per finished scenario goes to this script's log (stdout);
  # tracebacks go to the runner log. One BLAS thread per process.
  exit <- system2(
    uv,
    args,
    stdout = "",
    stderr = log_file,
    env = c(
      "MPLBACKEND=Agg",
      "OPENBLAS_NUM_THREADS=1",
      "OMP_NUM_THREADS=1",
      "MKL_NUM_THREADS=1"
    )
  )
  status_file <- file.path(out_dir, "status.csv")
  if (!identical(exit, 0L) || !file.exists(status_file)) {
    log_tail <- utils::tail(readLines(log_file, warn = FALSE), 20L)
    stop(
      method,
      " runner failed (exit status ",
      exit,
      "):\n",
      paste(log_tail, collapse = "\n"),
      call. = FALSE
    )
  }
  status <- utils::read.csv(status_file, stringsAsFactors = FALSE)
  status$scenario_id <- as.character(status$scenario_id)
  status$message[is.na(status$message)] <- ""

  dplyr::bind_rows(lapply(scenarios, function(scenario) {
    row <- status[status$scenario_id == scenario$scenario_id, , drop = FALSE]
    n_bulk <- ncol(scenario$bulk)
    celltypes <- colnames(scenario$mean_signature_matrix)
    if (nrow(row) != 1L || row$status != "ok") {
      message <- if (nrow(row) == 1L) row$message else "no status reported"
      return(.competitor_rows(
        scenario$scenario_id,
        method,
        celltypes,
        n_bulk,
        estimate = NULL,
        elapsed = if (nrow(row) == 1L) row$elapsed_sec else NA_real_,
        status = paste("error:", message)
      ))
    }
    scenario_out <- file.path(out_dir, scenario$scenario_id)
    read_table <- function(name) {
      path <- file.path(scenario_out, name)
      if (!file.exists(path)) {
        return(NULL)
      }
      x <- utils::read.csv(path, row.names = 1L, check.names = FALSE)
      x <- as.matrix(x)
      x[celltypes, , drop = FALSE]
    }
    .competitor_rows(
      scenario$scenario_id,
      method,
      celltypes,
      n_bulk,
      estimate = read_table("estimates.csv"),
      lower = read_table("lower.csv"),
      upper = read_table("upper.csv"),
      extra = Filter(
        Negate(is.null),
        list(
          lower_fisher = read_table("lower_fisher.csv"),
          upper_fisher = read_table("upper_fisher.csv"),
          lower_godambe = read_table("lower_godambe.csv"),
          upper_godambe = read_table("upper_godambe.csv")
        )
      ),
      elapsed = row$elapsed_sec,
      status = "ok"
    )
  }))
}

# One log line per competitor: scenarios fitted, failures and wall time.
.log_competitor <- function(method, rows, started) {
  status <- unique(rows[c("scenario_id", "status")])
  # MEAD's status also records its route ("ok; get_R01(...)"); failures
  # start with "error".
  failed <- status$scenario_id[startsWith(status$status, "error")]
  elapsed <- round(proc.time()[["elapsed"]] - started)
  cli::cli_alert_success(
    "{method}: {nrow(status)} scenarios in {elapsed} s, {length(failed)} failed."
  )
  if (length(failed) > 0L) {
    cli::cli_alert_warning("{method} failed on {.val {failed}}.")
  }
}

.run_mead_scenario <- function(scenario, level) {
  .check_suggested_mead()
  celltypes <- colnames(scenario$mean_signature_matrix)
  n_bulk <- ncol(scenario$bulk)
  counts <- do.call(cbind, scenario$reference)
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  col_data <- S4Vectors::DataFrame(
    cell_type = factor(
      rep(celltypes, vapply(scenario$reference, ncol, integer(1))),
      levels = celltypes
    ),
    individual = factor(unlist(scenario$individual, use.names = FALSE))
  )
  ref <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  started <- proc.time()[["elapsed"]]
  r01 <- tryCatch(
    MEAD::get_R01(scenario$bulk, alpha = 0.1),
    error = function(e) e
  )
  r01_route <- "get_R01(alpha = 0.1)"
  if (inherits(r01, "error")) {
    r01_route <- paste0(
      "R01 = NULL (get_R01 failed: ",
      conditionMessage(r01),
      ")"
    )
    r01 <- NULL
  }
  # filter.gene = FALSE prints "Filtering ref genes recommended"; the
  # scenarios are synthetic, so the filter stays off and the note is muted.
  fit <- tryCatch(
    suppressMessages(MEAD::MEAD(
      scenario$bulk,
      ref,
      cell_types = celltypes,
      preprocessing_control = list(filter.gene = FALSE),
      R01 = r01
    )),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    return(.competitor_rows(
      scenario$scenario_id,
      "MEAD",
      celltypes,
      n_bulk,
      estimate = NULL,
      elapsed = elapsed,
      status = paste0("error: ", conditionMessage(fit), "; ", r01_route)
    ))
  }
  ci <- MEAD::get_ci(fit, alpha = 1 - level)
  .competitor_rows(
    scenario$scenario_id,
    "MEAD",
    celltypes,
    n_bulk,
    estimate = fit$p_hat[celltypes, , drop = FALSE],
    lower = t(ci$ci_l)[celltypes, , drop = FALSE],
    upper = t(ci$ci_r)[celltypes, , drop = FALSE],
    elapsed = elapsed,
    status = paste0("ok; ", r01_route)
  )
}

.check_suggested_mead <- function() {
  needed <- c("MEAD", "SingleCellExperiment", "S4Vectors", "vashr", "Rfast")
  available <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  missing <- needed[!available]
  if (length(missing) > 0L) {
    stop(
      "MEAD needs ",
      paste(missing, collapse = ", "),
      ". Install with pak::pak(c(\"Rfast\", \"bioc::SingleCellExperiment\", ",
      "\"mengyin/vashr\", \"DongyueXie/MEAD\")).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# One row per (bulk column, cell type); NULL matrices give NA columns.
.competitor_rows <- function(
  scenario_id,
  algorithm,
  celltypes,
  n_bulk,
  estimate,
  lower = NULL,
  upper = NULL,
  extra = list(),
  elapsed,
  status
) {
  as_values <- function(x) {
    if (is.null(x)) {
      return(rep(NA_real_, n_bulk * length(celltypes)))
    }
    as.numeric(x)
  }
  rows <- tibble::tibble(
    scenario_id = scenario_id,
    algorithm = algorithm,
    replicate = rep(seq_len(n_bulk), each = length(celltypes)),
    celltype = rep(celltypes, times = n_bulk),
    estimate = as_values(estimate),
    lower = as_values(lower),
    upper = as_values(upper),
    elapsed_sec = elapsed / n_bulk,
    status = status
  )
  for (name in names(extra)) {
    rows[[name]] <- as_values(extra[[name]])
  }
  rows
}
