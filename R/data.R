#' Bivariate parameter configuration
#'
#' The complete parameter configurations used to model the bivariate framework
#'
#' @format ## `who`
#' A data frame with 392 rows and 7 columns:
#' \describe{
#'   \item{ID}{The unique ID identifying each general bivariate configuration tested}
#'   \item{overlap}{The average overlap}
#'   \item{entropy}{The Shannon entropy computed using the respective proportions of the two populations tested}
#'   \item{proportions, variance, centroids}{General descriptions about each parameter configuration}
#'   \item{true_parameters}{The actual list of parameters, p, mu and sigma, used to model the scenario}
#' }
"bivariate_configuration"



#' Bivariate parameter configuration
#'
#' The complete parameter configurations used to model the bivariate framework
#'
#' @format ## `bivariate_parameters`
#' A data frame with 7,838,260 rows and 8 columns:
#' \describe{
#'   \item{ID}{The unique ID identifying each general bivariate configuration tested}
#'   \item{correlation_celltype1,correlation_celltype2}{The pairwise correlation between the two genes, respectively 
#'   for population 1 and 2}
#'   \item{model_mse,...}{Some general metrics comparing the performance of the deconvolution algorithm 
#'   in inferring the cellular ratios, with respect to the ones provided in the simulation}
#'   \item{p1,p2}{Proportions of the two cell populations}
#'   \item{OMIC_ID}{ID characterinsg the sample}
#'   \item{algorithm}{The deconvolution algorithm used to infer the ratios}
#' }
"bivariate_parameters"