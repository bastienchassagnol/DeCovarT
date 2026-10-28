# Likelihood objects of the benchmark, in the manner of stats::family():
# one constructor per model holds everything model-specific (simulator,
# log-likelihood and its derivatives, parameters, interval menu), so the
# benchmark code never branches on the name of the model.

#' Likelihood objects for the benchmark
#'
#' @description
#' `population_likelihood()` is the Gaussian convolution
#' \eqn{\mathcal{N}(\boldsymbol{\mu}\boldsymbol{p},
#' \sum_j p_j^2\boldsymbol{\Sigma}_j)};
#' `single_cell_likelihood()` is the mixture of cells, with its two scales
#' \eqn{(c,\phi)} profiled out, so every function is evaluated at the
#' closed-form scales of the composition. They play the part of
#' [stats::family()] in [stats::glm()]: [deconvolute_ratios()] and
#' [run_simulation_benchmark()] take one as `likelihood`, either the object
#' or the constructor itself.
#'
#' @details
#' Each object is a list of class `decovart_likelihood` with
#' * `name`: `"population"` or `"single_cell"`;
#' * `nuisance`: whether the model has the scales \eqn{(c,\phi)};
#' * `parameters(cell_names)`: the rows of every interval matrix, the
#'   cell types and, for the single-cell model, `mean_scale`,
#'   `variance_scale`, `platform_scale_a` and `cell_count_n`;
#' * `simulate(mean, Sigma, theta, n)`: `n` bulk samples drawn from the
#'   moments `mean` and `Sigma` with the composition and settings in
#'   `theta`; returns `Y`, the target `truth` and the true scales
#'   `nuisance` (`NULL` for the population model);
#' * `loglik(p, y, mu, Sigma)`, `score(p, y, mu, Sigma)`,
#'   `se(p, y, mu, Sigma)` (Wald standard errors of the composition) and
#'   `curvature(p, y, mu, Sigma)` (local-maximum flag and eigenvalue range
#'   of the Hessian in ILR coordinates, `NA` on the boundary);
#' * `interval(method, y, p_hat, p_true, mu, Sigma, ci)`: one method of
#'   the interval menu at the estimate `p_hat`, as a two-column matrix.
#'
#' The interval methods are
#' * `wald`: delta method on \eqn{\boldsymbol{p}}, symmetric;
#' * `ilr`: Wald interval in ILR coordinates mapped back to the simplex;
#' * `profile`: likelihood-ratio inversion; it refits from `y`, so it does
#'   not depend on the solver that produced `p_hat`;
#' * `bootstrap`: percentile interval of a parametric bootstrap refitted by
#'   `ci$optimiser`;
#' * `oracle`: `p_hat` plus or minus the Wald half-width evaluated at the
#'   true composition.
#'
#' For the single-cell model, `_reml` appended to `wald`, `ilr` or
#' `profile` applies the REML calibration (\eqn{\nu=G-J}), and Wald, ILR
#' and bootstrap intervals also cover the scales and the derived
#' \eqn{a=\phi/c} and \eqn{n=c^{2}/\phi}.
#'
#' Wald, ILR and oracle intervals are undefined when `p_hat` has a zero
#' coordinate: on a face (one zero), an edge (two zeros) or a vertex
#' (\eqn{J-1} zeros, one cell type making up the whole sample). L-BFGS-B
#' returns such estimates; Newton-Raphson in ILR coordinates does not.
#' Their rows are then `NA`, scales included, and the coverage of these
#' methods is conditional on an interior estimate.
#'
#' @return An object of class `decovart_likelihood`.
#'
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2, dimnames = list(NULL, c("a", "b")))
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' lik <- population_likelihood()
#' lik
#' y <- drop(mu %*% c(0.4, 0.6))
#' lik$loglik(c(0.4, 0.6), y, mu, Sigma)
#' kkt_residual(c(0.4, 0.6), lik$score(c(0.4, 0.6), y, mu, Sigma))
#' single_cell_likelihood()$parameters(c("a", "b"))
#' @seealso [deconvolute_ratios()], [run_simulation_benchmark()]
#' @export
population_likelihood <- function() {
  structure(
    list(
      name = "population",
      nuisance = FALSE,
      parameters = function(cell_names) cell_names,
      simulate = function(mean, Sigma, theta, n) {
        simulated <- simulate_bulk_mixture(
          signature_matrix = mean,
          Sigma = Sigma,
          p = theta$p,
          n = n,
          truncate_negative = FALSE
        )
        list(Y = simulated$Y, truth = theta$p, nuisance = NULL)
      },
      loglik = function(p, y, mu, Sigma) {
        loglik_multivariate(p, y, mu, Sigma)
      },
      score = function(p, y, mu, Sigma) {
        gradient_loglik_unconstrained(p, y, mu, Sigma)
      },
      se = function(p, y, mu, Sigma) {
        .ilr_wald_se(p, mu, Sigma, warn = FALSE)
      },
      curvature = function(p, y, mu, Sigma) {
        if (min(p) <= 0) {
          return(.no_curvature)
        }
        diagnostics <- boundary_diagnostics(p, y, mu, Sigma)
        as.list(diagnostics[c(
          "local_maximum",
          "min_eigenvalue",
          "max_eigenvalue"
        )])
      },
      interval = function(method, y, p_hat, p_true, mu, Sigma, ci) {
        switch(
          method,
          wald = .wald_bounds(p_hat, .ilr_wald_se(p_hat, mu, Sigma), ci$level),
          oracle = .wald_bounds(
            p_hat,
            .ilr_wald_se(p_true, mu, Sigma),
            ci$level
          ),
          ilr = suppressWarnings(confint_ilr_backtransform(
            y,
            mu,
            Sigma,
            p = p_hat,
            level = ci$level,
            n_draws = ci$n_draws
          ))$interval,
          profile = confint_profile_decovart(
            y,
            mu,
            Sigma,
            level = ci$level
          )[, c("lower", "upper"), drop = FALSE],
          bootstrap = bootstrap_decovart(
            y,
            mu,
            Sigma,
            p = p_hat,
            n_boot = ci$n_boot,
            level = ci$level,
            optimiser = ci$optimiser
          )$interval,
          cli::cli_abort(
            "Unknown interval method {.val {method}} for the population model."
          )
        )
      }
    ),
    class = "decovart_likelihood"
  )
}

#' @rdname population_likelihood
#' @export
single_cell_likelihood <- function() {
  scales_at <- function(p, y, mu, Sigma) {
    .single_cell_scales_at_p(y, mu, Sigma, p, strict = FALSE)
  }
  structure(
    list(
      name = "single_cell",
      nuisance = TRUE,
      parameters = function(cell_names) {
        c(
          cell_names,
          "mean_scale",
          "variance_scale",
          "platform_scale_a",
          "cell_count_n"
        )
      },
      simulate = function(mean, Sigma, theta, n) {
        simulated <- simulate_bulk_from_cells(
          mean,
          Sigma,
          p = theta$p,
          n_cells = theta$n_cells,
          platform_scale = if (is.null(theta$platform_scale)) {
            1
          } else {
            theta$platform_scale
          },
          n_samples = n
        )
        # The target of a single-cell bulk is its realised composition.
        list(
          Y = simulated$Y,
          truth = simulated$p_realised,
          nuisance = rbind(
            mean_scale = simulated$mean_scale,
            variance_scale = simulated$variance_scale
          )
        )
      },
      loglik = function(p, y, mu, Sigma) {
        scales <- scales_at(p, y, mu, Sigma)
        if (is.null(scales)) {
          return(NA_real_)
        }
        loglik_single_cell(
          p,
          y,
          mu,
          Sigma,
          scales$mean_scale,
          scales$variance_scale
        )
      },
      score = function(p, y, mu, Sigma) {
        scales <- scales_at(p, y, mu, Sigma)
        gradient_loglik_single_cell(
          p,
          y,
          mu,
          Sigma,
          scales$mean_scale,
          scales$variance_scale
        )
      },
      se = function(p, y, mu, Sigma) {
        vcov <- vcov_single_cell(p, y, mu, Sigma, warn = FALSE)$vcov
        stats::setNames(sqrt(pmax(diag(as.matrix(vcov)), 0)), names(p))
      },
      curvature = function(p, y, mu, Sigma) {
        if (min(p) <= 0) {
          return(.no_curvature)
        }
        # Same criterion as boundary_diagnostics(): small ILR score and a
        # negative-definite Hessian of the profile log-likelihood.
        z <- isometric_log_ratio(p)
        score <- gradient_profile_loglik_single_cell_ilr(z, y, mu, Sigma)
        ev <- eigen(
          hessian_profile_loglik_single_cell_ilr(z, y, mu, Sigma),
          symmetric = TRUE,
          only.values = TRUE
        )$values
        list(
          local_maximum = sqrt(sum(score^2)) < 1e-4 && max(ev) < 0,
          min_eigenvalue = min(ev),
          max_eigenvalue = max(ev)
        )
      },
      interval = function(method, y, p_hat, p_true, mu, Sigma, ci) {
        reml <- grepl("_reml$", method)
        base <- sub("_reml$", "", method)
        calibration <- .single_cell_calibration(
          if (reml) "reml" else "ml",
          n_genes = nrow(mu),
          n_types = length(p_hat)
        )
        switch(
          base,
          wald = ,
          ilr = {
            pieces <- .composition_intervals(
              p_hat,
              y,
              mu,
              Sigma,
              ci$level,
              how = base,
              n_draws = ci$n_draws,
              seed = NULL,
              calibration = calibration
            )
            rbind(pieces$interval, pieces$nuisance$interval)
          },
          oracle = {
            truth <- vcov_single_cell(p_true, y, mu, Sigma, warn = FALSE)
            .wald_bounds(p_hat, sqrt(pmax(diag(truth$vcov), 0)), ci$level)
          },
          profile = confint_profile_single_cell(
            y,
            mu,
            Sigma,
            level = ci$level,
            variance = calibration$variance
          )[, c("lower", "upper"), drop = FALSE],
          bootstrap = {
            boot <- bootstrap_single_cell(
              y,
              mu,
              Sigma,
              n_boot = ci$n_boot,
              level = ci$level,
              optimiser = ci$optimiser
            )
            rbind(boot$interval, boot$scale_interval)
          },
          cli::cli_abort(
            "Unknown interval method {.val {method}} for the single-cell model."
          )
        )
      }
    ),
    class = "decovart_likelihood"
  )
}

#' @export
print.decovart_likelihood <- function(x, ...) {
  cat(
    "<decovart_likelihood: ",
    x$name,
    if (isTRUE(x$nuisance)) ", scales (c, phi) profiled out" else "",
    ">\n",
    sep = ""
  )
  invisible(x)
}

#' @noRd
.no_curvature <- list(
  local_maximum = NA,
  min_eigenvalue = NA_real_,
  max_eigenvalue = NA_real_
)

#' Accept a likelihood object or its constructor, as glm() accepts a family
#' @noRd
.as_likelihood <- function(likelihood) {
  if (is.function(likelihood)) {
    likelihood <- likelihood()
  }
  if (!inherits(likelihood, "decovart_likelihood")) {
    cli::cli_abort(c(
      "{.arg likelihood} must be a likelihood object.",
      i = "Use {.fn population_likelihood} or {.fn single_cell_likelihood}."
    ))
  }
  likelihood
}

#' Rows of an interval matrix, padded and ordered to the model's parameters
#'
#' Every bulk then has the same rows, whatever the method returned; a
#' missing row is an `NA` interval.
#' @noRd
.pad_bounds <- function(x, parameters) {
  out <- matrix(
    NA_real_,
    length(parameters),
    2L,
    dimnames = list(parameters, c("lower", "upper"))
  )
  if (is.null(x)) {
    return(out)
  }
  x <- as.matrix(x)
  if (ncol(x) < 2L || is.null(rownames(x))) {
    return(out)
  }
  keep <- intersect(parameters, rownames(x))
  out[keep, ] <- x[keep, seq_len(2L), drop = FALSE]
  out
}

#' The interval menu of one fitted bulk sample
#'
#' Evaluates `ci$methods` at the estimate `p_hat`. Every result has the rows
#' `likelihood$parameters(names(p_hat))`. Wald, ILR and oracle intervals are
#' `NA` when `p_hat` has a zero coordinate (face, edge or vertex); a method
#' that fails gives `NA` bounds.
#'
#' @param ci List with `methods`, `level`, `n_boot`, `n_draws` and
#'   `optimiser`.
#' @return A named list of two-column matrices (lower, upper), one per
#'   method.
#' @noRd
.likelihood_intervals <- function(
  likelihood,
  y,
  p_hat,
  p_true,
  mean_signature_matrix,
  Sigma,
  ci
) {
  ci$level <- if (is.null(ci$level)) 0.95 else ci$level
  ci$n_boot <- if (is.null(ci$n_boot)) 200L else ci$n_boot
  ci$n_draws <- if (is.null(ci$n_draws)) 2000L else ci$n_draws
  ci$optimiser <- if (is.null(ci$optimiser)) "newton" else ci$optimiser
  parameters <- likelihood$parameters(names(p_hat))
  on_boundary <- any(p_hat <= 0)
  stats::setNames(
    lapply(ci$methods, function(method) {
      base <- sub("_reml$", "", method)
      if (on_boundary && base %in% c("wald", "ilr", "oracle")) {
        return(.pad_bounds(NULL, parameters))
      }
      out <- tryCatch(
        suppressWarnings(likelihood$interval(
          method,
          y,
          p_hat,
          p_true,
          mean_signature_matrix,
          Sigma,
          ci
        )),
        error = function(e) NULL
      )
      .pad_bounds(out, parameters)
    }),
    ci$methods
  )
}
