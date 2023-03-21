##################################################################
##                    boxplot representation                    ##
##################################################################

plot_boxplots_parameters <- function(distribution_parameters, score_variable=".data[[score_variable]]", with_outliers = FALSE) {

  if (!with_outliers) {
    # remove outlying points
    distribution_parameters <- distribution_parameters %>%
      dplyr::group_by(correlation_celltype1, correlation_celltype2) %>%
      dplyr::filter(.data[[score_variable]] > (quantile(.data[[score_variable]], probs = c(0.25)) - IQR(.data[[score_variable]])) &
                      .data[[score_variable]] < (quantile(.data[[score_variable]], probs = c(0.75)) +  IQR(.data[[score_variable]])))
  }

  # display parameters, removing outliers for visualisation purposes
  boxplot_parameters <- ggplot(distribution_parameters, aes(x = deconvolution_name, y = .data[[score_variable]], fill = variance)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 0.5, position = position_dodge(width = 0.9), width = 0.4) +
    stat_summary(
      fun = mean, geom = "point", shape = 3, size = 1, colour = "yellow",
      position = position_dodge(width = 0.9), show.legend = FALSE
    ) +
    facet_wrap(~correlation_celltype1 + correlation_celltype2) + # same correlation for cell type 1 in a given row
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 35, size = 12, vjust = 0.6),
      axis.ticks.length = unit(.1, "cm"),
      legend.text = element_text(size = 16),
      strip.text = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      title = element_blank(),
      panel.spacing = unit(.2, "pt")
    ) +
    scale_fill_viridis_d()

  # boxplot_parameters <- egg::tag_facet(boxplot_parameters,
  #                                      tag_pool = unique(distribution_parameters$name_parameter),
  #                                      open = "", close = "", hjust = -0.2, size = 5)
  return(boxplot_parameters)
}


##################################################################
##                    multivariate density     ##
##################################################################

plot_multivariate_density <- function(X, n=dim(X)[3]) {
  # concatenate gene expression from the two distinct cell populations
  X <- matrix(X, nrow = n *2, ncol = 2, byrow = T,
              dimnames = list(paste(dimnames(X)[[2]], rep(dimnames(X)[[3]], each=2), sep = "_"),dimnames(X)[[1]])) %>%
    tibble::as_tibble(rownames = "id") %>% dplyr::mutate(celltype=rep(dimnames(X)[[2]], n))

  isogradient <- ggplot(X, aes(x = gene_1, y = gene_2)) +
    stat_density_2d(geom = "tile",
                    aes(fill = ..density..),
                    contour = FALSE) +
    viridis::scale_fill_viridis() +
    labs(x ="Expression of gene 1",y = "Expression of gene 2") +
    theme_bw() +
    xlim(c(min(X[, c("gene_1")])-5, max(X[, c("gene_1")]) + 5)) +
    ylim(c(min(X[, c("gene_2")])-5, max(X[, c("gene_2")]) + 5)) +
    geom_point(mapping = aes(x=gene_1, y=gene_2, col=celltype, shape=celltype), data = mu_tibble,
               inherit.aes = F, show.legend = T, size=4) +
    scale_color_manual(values=c("red", "green"))


  return(isogradient)
}


# anti-log if max < 50 in mixture or signature file
# if(max (X)<50) {X=2**X; warning ("Log-scale spotted for reference signature")}
# if(max (Y)<50) {Y=2**Y; warning ("Log-scale spotted for mixture profile")}

# standardize both matrix, not really advised in benchmarking papers
# if (scaled) {
#   X <- scale(X);  Y <- scale(Y)
# }
# cibersort_ratios <- purrr::map_dfr(Y, ~ deconvolute_ratios_CIBERSORT(., X=X)) %>%
#     dplyr::mutate(OMIC_ID=colnames(Y))
