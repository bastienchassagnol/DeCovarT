
# Core algorithm (where svm regression is performed for each mixture)
deconvolute_ratios_CIBERSORT <- function(y, X, true_ratios=NULL, ...){
     #the set of nu values tested for best performance
    range.nu<-seq(0.2, 0.8, 0.3)
    model<-e1071::best.svm(X,y, type="nu-regression",kernel="linear",scale=F, 
                      nu = range.nu, tunecontrol=e1071::tune.control(sampling = "fix", fix=0.75)) #determine best fitted model
    
    e1071::tune.svm(X,y,type="nu-regression",kernel="linear",scale=F, 
                    nu = 0.25, tunecontrol=e1071::tune.control(sampling = "fix", fix=0.75))
    
    model <- e1071::svm(X, y=y)
    
    #get and normalize coefficients (sum to 1 and no negative coefficients)
    estimated_ratios <- t(model$coefs) %*% model$SV
    estimated_ratios[which(estimated_ratios<0)]<-0
    estimated_ratios <- estimated_ratios/sum(estimated_ratios)
    
    metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
      dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
    return(metrics_scores)
}

# linear regression
deconvolute_ratios_abbas <- function(y, X, true_ratios=NULL, ...) {
  estimated_ratios <- stats::lsfit(X, y, intercept = F)$coefficients
  
  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_ratios[which(estimated_ratios<0)]<-0; estimated_ratios <- estimated_ratios/sum(estimated_ratios) 
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
  return(metrics_scores)
}

# robust linear regression
deconvolute_ratios_monaco <- function(y, X, true_ratios=NULL, ...) {
  estimated_ratios <- MASS::rlm(y ~ X+ 0, method = c("M"))$coefficients; names(estimated_ratios) <- colnames(X)
  
  # normalize coefficients (sum to 1 and no negative coefficients)
  estimated_ratios[which(estimated_ratios<0)]<-0; estimated_ratios <- estimated_ratios/sum(estimated_ratios) 
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
  return(metrics_scores)
}


# quadratic programming (lsei function, other possibility is to use lsqlin function)
deconvolute_ratios_deconRNASeq <- function(y, X, true_ratios=NULL, ...) {

EE <- rep(1, ncol(X)); FF <- 1 # encode the sum-to-one constraint
GG <- diag(nrow=ncol(X)); HH <- rep(0, ncol(X)) # encode the non-negativity constraint

estimated_ratios <- limSolve::lsei(X, y, EE, FF, GG, HH)$X
metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
  dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
return(metrics_scores)
}

# nnls scoring
deconvolute_ratios_nnls <- function(y, X, true_ratios=NULL, ...) {
  estimated_ratios <- nnls::nnls(X, y)$x
  
  # normalize coefficients (sum to 1, as non-negativity is already enforced)
  estimated_ratios <- estimated_ratios/sum(estimated_ratios); names(estimated_ratios) <- colnames(X)
  metrics_scores <- compute_benchmark_metrics(y, X, estimated_ratios, true_ratios) %>%
    dplyr::bind_cols(tibble::as_tibble_row(estimated_ratios))
  return(metrics_scores)
}



# compute summary scores
compute_benchmark_metrics <- function(y, X, estimated_ratios, true_ratios=NULL) {
  n <- nrow(X); k <- ncol(X)
  df_res <- n - k + 1 # number of free parameters: only k - 1 parameters must be learnt, with sum-to-one constraint
  df_tot <- n - 1 # no intercept for the moment, in the model (so n-1, or n?)
  if(!is.null(true_ratios)) { # when the true parameters are known
    scores <- tibble::tibble(model_mse = Metrics::mse(true_ratios, estimated_ratios),
                             model_rmse = Metrics::rmse(true_ratios, estimated_ratios),
                             model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_ratios)),
                             model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
                             model_mae = Metrics::mae(true_ratios, estimated_ratios),
                             model_cor = suppressWarnings(cor(true_ratios, estimated_ratios, method = "pearson")))
   
  }
  else { # when they are unknown
    predicted_values <- as.vector(X%*%t(estimated_ratios))
    scores <- tibble::tibble(model_mse = Metrics::mse(y, predicted_values),
                             model_rmse = Metrics::rmse(y, predicted_values),
                             model_coef_determination = max(0, 1 - Metrics::rse(true_ratios, estimated_ratios)),
                             model_coef_determination_adjusted=max(0, 1 - (1-model_coef_determination) * df_tot/df_res),
                             model_mae = Metrics::mae(y, predicted_values),
                             model_cor = suppressWarnings(cor(y, predicted_values, method = "pearson")),
                             model_ccc= epiR::epi.ccc(y, predicted_values)) # only relevant comparing two sets of measured or simulated variables, not relevant in simulations)
    
  }
  return(scores)
}




# main function, launching all deconvolution algorithms
deconvolute_ratios <- function(signature_matrix, bulk_expression, scaled=F, true_ratios=NULL, Sigma=NULL,
                               deconvolution_functions=list("QP" =list(FUN=deconvolute_ratios_deconRNASeq,additionnal_parameters=NULL),
                                                           "lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL),
                                                           "rlm"=list(FUN=deconvolute_ratios_monaco,additionnal_parameters=NULL),
                                                           "nnls" = list(FUN=deconvolute_ratios_nnls,additionnal_parameters=NULL))){
    #read in data
    if (!is.matrix (signature_matrix) | is.null(row.names (signature_matrix))) 
      stop ("required format for signature is expression matrix, with rownames as genes")
    X <- tibble::as_tibble(signature_matrix, rownames="GENE_SYMBOL")
    if (!is.matrix (bulk_expression) | is.null(row.names (bulk_expression)))
      stop ("required format for mixture is expression matrix, with rownames as genes")
    Y <- tibble::as_tibble(bulk_expression, rownames="GENE_SYMBOL")
    
    # remove potential missing data from Y database
    Y <- Y %>% tidyr::drop_na()
    
    # intersect genes (we only keep genes that are common to both data bases)
    common_genes<-intersect(X$GENE_SYMBOL, Y$GENE_SYMBOL)
    
    if (length(common_genes)/dim(X)[1]<0.5) {
      stop (paste("Only", length(common_genes)/dim(X)[1],"fraction of genes are used in the signature matrix\n. 
                  Half of common genes are required at least"))}
    X <- X %>% dplyr::filter(GENE_SYMBOL %in% common_genes) %>% 
      dplyr::arrange(GENE_SYMBOL) %>% 
      dplyr::select(where(is.numeric)) %>% 
      as.matrix()
    Y <- Y %>% dplyr::filter(GENE_SYMBOL %in% common_genes) %>% 
      dplyr::arrange(GENE_SYMBOL) %>% 
      dplyr::select(where(is.numeric)) 


    # anti-log if max < 50 in mixture or signature file
    # if(max (X)<50) {X=2**X; warning ("Log-scale spotted for reference signature")}
    # if(max (Y)<50) {Y=2**Y; warning ("Log-scale spotted for mixture profile")} 
    
    # standardize both matrix, not really advised in benchmarking papers
    if (scaled) {
      X <- scale(X);  Y <- scale(Y)
    }
    # cibersort_ratios <- purrr::map_dfr(Y, ~ deconvolute_ratios_CIBERSORT(., X=X)) %>% 
    #     dplyr::mutate(OMIC_ID=colnames(Y))
    
    deconvolution_estimates <- purrr::imap_dfr(deconvolution_functions, function(deconvolution_function, deconvolution_name) {
      additionnal_parameters <- deconvolution_function$additionnal_parameters
      estimated_ratios <- purrr::map_dfr(Y, ~ do.call(deconvolution_function$FUN, c(list("y"=., "X"=X, "Sigma"=Sigma, "true_ratios"=true_ratios),
                                                                                        additionnal_parameters))) %>% 
        dplyr::mutate(OMIC_ID=colnames(Y), deconvolution_name=deconvolution_name)
      return(estimated_ratios)
    })
    return (deconvolution_estimates)
}




