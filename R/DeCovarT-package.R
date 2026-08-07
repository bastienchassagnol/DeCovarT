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
#' `browseVignettes(package = "DeCovarT")` and
#' `?DeCovarT::deconvolute_ratios`.
#'
#' Reference manuals (PDF / HTML) can be regenerated with
#' `source("scripts/auxiliary/generate_package_manual.R")`.
#'
#' @keywords internal
"_PACKAGE"
