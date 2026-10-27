# Benchmark artefacts shared by scripts/fig02_bivariate_toy.R,
# scripts/fig03_variance_covariance_driven.R and
# scripts/fig05_single_cell_hierarchical.R.
#
# Every figure keeps two files in output/figXX/:
# * figXX_benchmark.rds: list(meta, config, theta_true, descriptors,
#   supplementary, regression, monte_carlo, optimisation, intervals, data),
#   with every method (DeCovarT solvers, first-generation baselines and
#   second-generation competitors) scored on the same simulated bulks;
# * figXX_plot_data.rds: one entry per figure file stem, holding the data
#   handed to its plot builder (see scripts/visualisation/).
#
# The run size comes from environment variables, so a dry run and the
# production run share one script:
#   DECOVART_N_REPLICATES  Monte Carlo bulks per scenario (N)
#   DECOVART_N_BOOT        parametric-bootstrap refits per bulk (B)
#   DECOVART_N_REFERENCE   individuals in the reference panel (R; default N)
#   DECOVART_POSTPROCESS_ONLY=1 redraws the figures from figXX_benchmark.rds.

benchmark_settings <- function(default_replicates, default_boot) {
  read_count <- function(name, default) {
    value <- Sys.getenv(name, "")
    if (!nzchar(value)) {
      return(as.integer(default))
    }
    value <- suppressWarnings(as.integer(value))
    if (is.na(value) || value < 1L) {
      stop("`", name, "` must be a positive integer.", call. = FALSE)
    }
    value
  }
  n_replicates <- read_count("DECOVART_N_REPLICATES", default_replicates)
  list(
    n_replicates = n_replicates,
    n_boot = read_count("DECOVART_N_BOOT", default_boot),
    n_reference = read_count("DECOVART_N_REFERENCE", n_replicates),
    postprocess_only = identical(
      Sys.getenv("DECOVART_POSTPROCESS_ONLY", "0"),
      "1"
    )
  )
}

benchmark_meta <- function(seed, settings) {
  commit <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  list(
    seed = seed,
    n_replicates = settings$n_replicates,
    n_boot = settings$n_boot,
    n_reference = settings$n_reference,
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_commit = commit,
    session_info = utils::sessionInfo()
  )
}

# Competitor stage -----------------------------------------------------------

# Scenarios for run_second_generation_batch(): the simulated bulks of every
# scenario, its independent reference panel, and the plug-in moments.
competitor_scenarios <- function(benchmark, config) {
  lapply(seq_len(nrow(config)), function(i) {
    theta <- config$true_theta[[i]]
    data <- benchmark$data[[i]]
    list(
      scenario_id = config$ID[[i]],
      bulk = data$bulk,
      reference = data$reference$reference,
      individual = data$reference$individual,
      mean_signature_matrix = theta$mu,
      Sigma = theta$sigma
    )
  })
}

# Score the competitor estimates on the same bulks as the package solvers.
#
# Every competitor is scored against the target of the benchmark (the cell
# fractions). MEAD rescales every reference cell to a unit total, so its
# proportions are transcript fractions q_j = p_j S_j / sum(p S), with S_j the
# column total of the mean; its intervals are also scored against q, as
# `native (transcript fractions)`. RNA-Sieve has three intervals: the
# released code (`native (released code)`, [0, 1] whenever its averaged
# Hessian is ill-conditioned, rnasieve 0.1.4, stats.py:224 and 234) and the
# paper's Fisher and Godambe intervals recomputed by the runner on all genes
# (`fisher`, `godambe`).
score_competitors <- function(
  benchmark,
  competitors,
  config,
  model = c("population", "single_cell"),
  coverage_interval = "wilson"
) {
  model <- match.arg(model)
  meta_cols <- setdiff(names(config), c("true_theta", "n"))
  attach_meta <- function(tbl, i) {
    if (nrow(tbl) == 0L) {
      return(tbl)
    }
    meta <- config[rep(i, nrow(tbl)), meta_cols, drop = FALSE]
    dplyr::bind_cols(meta, tbl[setdiff(names(tbl), meta_cols)])
  }
  scored <- list()
  for (i in seq_len(nrow(config))) {
    theta <- config$true_theta[[i]]
    mu <- theta$mu
    data <- benchmark$data[[i]]
    functions <- DeCovarT:::.benchmark_model(model, mu, theta$sigma)
    rows_i <- competitors[competitors$scenario_id == config$ID[[i]], ]
    for (algorithm in unique(rows_i$algorithm)) {
      rows <- rows_i[rows_i$algorithm == algorithm, ]
      as_matrix <- function(column) {
        wide <- tidyr::pivot_wider(
          rows[c("celltype", "replicate", column)],
          names_from = "replicate",
          values_from = tidyr::all_of(column)
        )
        out <- as.matrix(wide[, -1L])
        rownames(out) <- wide$celltype
        out[colnames(mu), order(as.integer(colnames(out))), drop = FALSE]
      }
      p_hat <- as_matrix("estimate")
      bounds <- list(native = c("lower", "upper"))
      if (identical(algorithm, "RNA-Sieve")) {
        bounds <- list(
          `native (released code)` = c("lower", "upper"),
          fisher = c("lower_fisher", "upper_fisher"),
          godambe = c("lower_godambe", "upper_godambe")
        )
        bounds <- bounds[vapply(
          bounds,
          function(cols) all(cols %in% names(rows)),
          logical(1)
        )]
      }
      bounds <- lapply(bounds, function(cols) {
        list(lower = as_matrix(cols[[1L]]), upper = as_matrix(cols[[2L]]))
      })
      n_samples <- ncol(p_hat)
      per_sample <- function(fun) {
        vapply(
          seq_len(n_samples),
          function(s) {
            tryCatch(fun(s), error = function(e) NA_real_)
          },
          numeric(1)
        )
      }
      finite <- apply(p_hat, 2L, function(x) all(is.finite(x)))
      bounds <- Filter(function(b) any(is.finite(b$lower)), bounds)
      intervals <- if (length(bounds) > 0L) {
        lapply(seq_len(n_samples), function(s) {
          lapply(bounds, function(b) {
            cbind(lower = b$lower[, s], upper = b$upper[, s])
          })
        })
      }
      metrics <- compute_benchmark_metrics(
        y = data$bulk,
        mean_signature_matrix = mu,
        estimated_p = p_hat,
        true_ratios = data$truth,
        elapsed = tapply(rows$elapsed_sec, rows$replicate, mean)[
          as.character(seq_len(n_samples))
        ],
        numerical_converged = finite,
        kkt_residual = per_sample(function(s) {
          DeCovarT:::.kkt_residual(
            p_hat[, s],
            functions$score(p_hat[, s], data$bulk[, s])
          )
        }),
        loglik_hat = per_sample(function(s) {
          functions$loglik(p_hat[, s], data$bulk[, s])
        }),
        loglik_true = per_sample(function(s) {
          functions$loglik(data$truth[, s], data$bulk[, s])
        }),
        coverage_interval = coverage_interval,
        algorithm = algorithm,
        intervals = intervals
      )
      if (identical(algorithm, "MEAD") && !is.null(intervals)) {
        size <- colSums(mu)
        transcript <- data$truth * size
        transcript <- sweep(transcript, 2L, colSums(transcript), "/")
        native_q <- DeCovarT:::.interval_table(
          intervals,
          p_hat = p_hat,
          p_true = transcript,
          coverage_interval = coverage_interval,
          algorithm = algorithm
        )
        native_q$ci_method <- "native (transcript fractions)"
        metrics$intervals <- dplyr::bind_rows(metrics$intervals, native_q)
      }
      metrics$optimisation$status <- rows$status[[1L]]
      scored[[length(scored) + 1L]] <- list(
        regression_global = attach_meta(metrics$regression$global, i),
        regression_cell_type = attach_meta(metrics$regression$cell_type, i),
        monte_carlo = attach_meta(metrics$monte_carlo, i),
        optimisation = attach_meta(metrics$optimisation, i),
        intervals = attach_meta(metrics$intervals, i)
      )
    }
  }
  bind <- function(name) {
    dplyr::bind_rows(lapply(scored, `[[`, name))
  }
  benchmark$regression$global <- dplyr::bind_rows(
    benchmark$regression$global,
    bind("regression_global")
  )
  benchmark$regression$cell_type <- dplyr::bind_rows(
    benchmark$regression$cell_type,
    bind("regression_cell_type")
  )
  benchmark$monte_carlo <- dplyr::bind_rows(
    benchmark$monte_carlo,
    bind("monte_carlo")
  )
  benchmark$optimisation <- dplyr::bind_rows(
    benchmark$optimisation,
    bind("optimisation")
  )
  benchmark$intervals <- dplyr::bind_rows(
    benchmark$intervals,
    bind("intervals")
  )
  benchmark
}

# Files ------------------------------------------------------------------------

# The reference panels are large and only feed the competitor stage, so the
# saved data keep the bulks, their targets and the true scales.
save_benchmark <- function(benchmark, config, meta, file) {
  benchmark$data <- lapply(benchmark$data, function(x) {
    x[c("bulk", "truth", "nuisance")]
  })
  artefact <- c(
    list(
      meta = meta,
      config = config[setdiff(names(config), "true_theta")],
      theta_true = stats::setNames(config$true_theta, config$ID)
    ),
    benchmark[c(
      "descriptors",
      "supplementary",
      "regression",
      "monte_carlo",
      "optimisation",
      "intervals",
      "data"
    )]
  )
  names(artefact$data) <- config$ID
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(artefact, file)
  invisible(file)
}

read_benchmark <- function(file) {
  if (!file.exists(file)) {
    stop(
      "No benchmark at `",
      file,
      "`; run the figure script without DECOVART_POSTPROCESS_ONLY first.",
      call. = FALSE
    )
  }
  readRDS(file)
}
