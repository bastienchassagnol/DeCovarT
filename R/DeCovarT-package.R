#' DeCovarT: covariance-aware bulk transcriptomic deconvolution
#'
#' @description
#' **DeCovarT** estimates cellular proportions from bulk RNA-seq by modelling
#' mixtures as Gaussian convolutions of purified cell-type means and
#' covariances. Notation follows the article: genes \eqn{g=1,\ldots,G}, cell
#' types \eqn{j=1,\ldots,J}, samples \eqn{i=1,\ldots,N}; bulk
#' \eqn{\boldsymbol{y}}, mean signature \eqn{\boldsymbol{\mu}}, proportions
#' \eqn{\boldsymbol{p}}, covariances / precisions
#' \eqn{\boldsymbol{\Sigma}_j}/\eqn{\boldsymbol{\Theta}_j}. Proportions live on
#' the open simplex and are optimised in unconstrained ALR coordinates
#' \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}
#' (`vignette("softmax-alr-derivatives", package = "DeCovarT")`).
#'
#' The frequentist API plugs in \eqn{\boldsymbol{\mu}} for unobserved latent
#' profiles \eqn{\boldsymbol{x}_{\cdot j}}; MAP recovery of those latents is
#' the Bayesian extension.
#'
#' To learn more, start with the package website vignettes, or in an R session:
#' `browseVignettes(package = "DeCovarT")`,
#' `?DeCovarT::fit_decovart` and
#' `?DeCovarT::deconvolute_ratios`.
#'
#' Reference manuals (PDF / HTML) can be regenerated with
#' `source("scripts/auxiliary/generate_package_manual.R")`.
#'
#' @srrstats {G1.4} Exported functions are documented with roxygen2.
#' @srrstats {G1.4a} Undocumented internals use `@noRd`; documented
#'   helpers are exported (with `@keywords internal` when low-level).
#' @srrstats {G1.3} Terminology: \eqn{\boldsymbol{p}} (simplex proportions),
#'   ALR coordinates \eqn{\boldsymbol{\rho}}, Gaussian convolution
#'   \eqn{\boldsymbol{y}\mid\boldsymbol{p}}, covariance
#'   \eqn{\boldsymbol{\Sigma}_j}. *Closed-reference*: every abundant type is a
#'   column of \eqn{\boldsymbol{\mu}}. *Gene-wise affine standardisation*:
#'   one centre/scale per gene, applied to bulk, means and covariances.
#'   A *sample-level covariate* \eqn{\boldsymbol{z}_{i}} is not a column of
#'   \eqn{\boldsymbol{\mu}}; it may shift \eqn{\boldsymbol{\mu}_{j}(\boldsymbol{z}_{i})}
#'   or the ALR of \eqn{\boldsymbol{p}_{\cdot i}}. *RNA fraction* vs *cell
#'   fraction*: transcriptome sizes \eqn{r_{j}} uncouple the two.
#'   Outlook (defined in
#'   `vignette("decovart-perspectives", package = "DeCovarT")`):
#'   *cell-type-specific (CTS)* sample-level latents
#'   \eqn{\boldsymbol{x}_{\cdot j}} (Bayesian / MAP);
#'   *spatial transcriptomics* as location-wise mixtures
#'   \eqn{\boldsymbol{y}(s)}; *log-normal* and *Poisson-log-normal*
#'   alternatives to the Gaussian bulk law.
#'
#' Paper-scale figures and large simulation assets are **not** bundled
#' here. They will live in a companion reproducibility repository,
#' following the CellRank 2 split between the maintained package
#' (`scverse/cellrank`) and `theislab/cellrank2_reproducibility`
#' (Weiler et al., *Nature Methods* 2024,
#' <https://doi.org/10.1038/s41592-024-02303-9>).
#'
#' Missing values: bulk and single-cell RNA-seq quantification yields
#' complete numeric matrices whose zeros are dropouts, not `NA`
#' \insertCite{hafemeisterNormalizationVarianceStabilization2019}{DeCovarT}.
#' Label-free proteomics, by contrast, routinely has 10-50% structurally
#' missing intensities (limit of detection, DDA ion sampling) that need
#' explicit missing-data models
#' \insertCite{chionBayesianFrameworkMultivariate2023}{DeCovarT}.
#' DeCovarT currently targets transcriptomes and therefore **errors** on
#' `NA` / `NaN` / `Inf` rather than imputing.
#'
#' @keywords internal
"_PACKAGE"
