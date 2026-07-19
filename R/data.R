#' Bivariate simulation design grid
#'
#' Parameter configurations for the bivariate (\eqn{G=2}, \eqn{J=2}) Gaussian
#' convolution benchmark: each row stores entropy of \eqn{\boldsymbol{p}},
#' mixture overlap, and the GMM triple
#' \eqn{(\boldsymbol{p},\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\})}.
#'
#' @format A data frame with 392 rows and 7 columns:
#' \describe{
#'   \item{ID}{Scenario identifier}
#'   \item{overlap}{Average pairwise overlap}
#'   \item{entropy}{Normalised Shannon entropy of \eqn{\boldsymbol{p}}}
#'   \item{proportions, variance, centroids}{Scenario factors}
#'   \item{true_parameters}{List with `p`, `mu`, `sigma`}
#' }
"bivariate_configuration"


#' Bivariate deconvolution benchmark results
#'
#' Long-format scores from deconvolving simulated
#' \eqn{\boldsymbol{Y}} under the bivariate design in
#' [bivariate_configuration].
#'
#' @format A data frame with many rows and 8 columns:
#' \describe{
#'   \item{ID}{Links to [bivariate_configuration]}
#'   \item{correlation_celltype1, correlation_celltype2}{Off-diagonal
#'     correlations entering \eqn{\boldsymbol{\Sigma}_1},
#'     \eqn{\boldsymbol{\Sigma}_2}}
#'   \item{model_mse, \ldots}{Estimation metrics for
#'     \eqn{\hat{\boldsymbol{p}}} versus \eqn{\boldsymbol{p}^{\star}}}
#'   \item{p1, p2}{Coordinates of \eqn{\boldsymbol{p}}}
#'   \item{OMIC_ID}{Sample identifier}
#'   \item{algorithm}{Deconvolution method}
#' }
"bivariate_parameters"
