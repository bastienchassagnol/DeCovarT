#' DeCovarT: covariance-aware bulk transcriptomic deconvolution
#'
#' @description
#' **DeCovarT** estimates cellular proportions from bulk RNA-seq by modelling
#' mixtures as Gaussian convolutions of purified cell-type means and
#' covariances. Proportions live on the open simplex and are optimised in an
#' unconstrained additive-logistic (ALR) coordinate system.
#'
#' To learn more, start with the package website vignettes, or in an R session:
#' `browseVignettes(package = "DeCovarT")` and
#' `?DeCovarT::deconvolute_ratios`.
#'
#' Reference manuals (PDF / HTML) can be regenerated with
#' `source("scripts/generate_package_manual.R")`.
#'
#' @keywords internal
"_PACKAGE"
