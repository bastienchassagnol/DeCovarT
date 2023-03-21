############################################################################
############################################################################
###                                                                      ###
###                     TEST NUMERICALLY DERIVATIVES                     ###
###                                                                      ###
############################################################################
############################################################################

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

#################################################################
##                  first part of the Hessian                  ##
#################################################################
library(tensorA) # load package
J_log <- to.tensor(1:3,c("i"=3)) # Jacobian of the original log-likelihood
H_psi <- to.tensor(1:12,c(j=2,k=2,"i"=3)) # Hessian of the mapping function

# Einstein summation
hess_einstein <- J_log %e% H_psi # %e% refers to the Einstein operator, repeated indexes (here i)
# in both tensors are coupled and summed over, of note, any number of tensors can be passed
# %e% operator, provided you respect the shape configuration

# multivariate tensor product
hess_tensor_product <- mul.tensor(J_log,"i",H_psi,"i") # the standard tensor product,
# in which you have to indicate explicitly the indexes to aggregate

#  with the tensor package
library(tensor)
# in opposition to the mult product, you instead indicate the distinct indexes
hess_tensor <- tensor(A = J_log, B = H_psi, alongA = 1, alongB = 3)

# with base R operators
hess_base <- matrix(0, nrow = 2, ncol = 2, dimnames = )
for (i in 1:3) {
  hess_base <- hess_base +
    (J_log[i] %>% as.numeric()) * (H_psi[,,i] %>% to.matrix.tensor(j="I2"))
}

# check that all these elements are equal
testthat::expect_equal(hess_rieman, hess_tensor %>% to.tensor()) # pass
testthat::expect_equal(hess_rieman, hess_tensor_product) # pass
# pass but requires cumbersome operation
testthat::expect_equal(hess_rieman, hess_base %>%
                         to.tensor() %>% extract2(I1 =~"j", I2 =~"k"))

##################################################################
##                  second part of the Hessian                  ##
##################################################################

J_psi <- to.tensor(1:6,c("i"=3, "j"=2))
H_log <- to.tensor(1:9, c("i"=3, "l"=3))

# with base R
hess_base_R <- t(J_psi) %*% H_log %*% J_psi
# %*% designs here the dot product

# with the Einstein summation (needs an index change)
hess_einstein_second <- J_psi %e% H_log %e% J_psi[[i =~"l", j =~"k"]]
testthat::expect_equal(hess_base_R,
                       hess_einstein_second %>% matrix(nrow=2, ncol=2)) # pass



# testthat::expect_equal(J_psi %e% H_log %e% J_psi,
# mul.tensor(J_psi, i="i", H_log, j="i") %>% mul.tensor(i="k", J_psi, j="k"))
#
#
# J_psi <- to.tensor(1:6,c("i"=3, "k"=2))
# H_log <- to.tensor(1:9, c("i"=3, "j"=3))
# riemann.tensor(J_psi, H_log[[i =~"^i"]], J_psi[[k =~"^k"]])


##################################################################
##                  test numerical derivatives                  ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")
equibalanced_ratios <- c(1/3, 1/3, 1/3)
testthat::expect_equal(equibalanced_ratios,
                       equibalanced_ratios %>% inverse_mapping_function() %>% mapping_function())

# check Jacobian of the mapping function
jacobian_mapping_theoretical <- jacobian_mapping_function(c(0, 0, 0))
jacobian_mapping_numerical <- numDeriv::jacobian(mapping_function, c(0, 0, 0), method="simple", method.args=list(eps=1e-6))
testthat::expect_equal(jacobian_mapping_theoretical, jacobian_mapping_numerical, tolerance = 1e-6) # pass



jacobian_mapping_theoretical <- jacobian_mapping_function(inverse_mapping_function(c(0.90, 0.05, 0.05)))
# privilege the Ridchardon and its two-way estimation of derivative #
jacobian_mapping_numerical <- numDeriv::jacobian(mapping_function, inverse_mapping_function(c(0.90, 0.05, 0.05)),
                                                 method="Richardson", method.args=list(eps=1e-6))
testthat::expect_equal(jacobian_mapping_theoretical, jacobian_mapping_numerical, tolerance = 1e-6) # pass


##----------------------------------------------------------------
##                  test gradient log-likelihood                 -
##----------------------------------------------------------------

set.seed(3)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))


jacobian_mapping_numerical <- numDeriv::grad(loglik_multivariate, p,
                                             method="Richardson", method.args=list(eps=1e-4, r=6),
                                             y=y, X=X, Sigma=Sigma) # additional arguments

jacobian_mapping_theoretical <- gradient_loglik_unconstrained (p, y, X, Sigma)
testthat::expect_equal(jacobian_mapping_numerical, jacobian_mapping_theoretical)

#################################################################
##               test constrained log-likelihood               ##
#################################################################
set.seed(3); library(tensorA)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))
theta <- inverse_mapping_function(p)

Sigma_global <- .compute_global_variance(p, Sigma)

theta_vect <- seq(-30, 30, 0.5)
log_vect <- purrr::map_dbl(theta_vect, loglik_multivariate_constrained, y=y, X=X, Sigma=Sigma)
plot(log_vect ~theta_vect)

loglik_multivariate_constrained(theta_vect[1], )

# matrix product of the gradient unconstrained evaluated in p times Jacobian of p
grad_constrained_mapping_numerical <- numDeriv::grad(loglik_multivariate_constrained, theta,
                                                     method="Richardson", method.args=list(eps=1e-4, r=6),
                                                     y=y, X=X, Sigma=Sigma) # additional arguments
grad_constrained_mapping_theoretical <- gradient_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(grad_constrained_mapping_numerical, grad_constrained_mapping_theoretical %>% as.numeric())

###  more complex derivative
set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
# y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))
theta <- inverse_mapping_function(p)

# matrix product of the gradient unconstrained evaluated in p times Jacobian of p
grad_constrained_mapping_numerical <- numDeriv::grad(loglik_multivariate_constrained, theta,
                                                     method="Richardson", method.args=list(eps=1e-12, r=4),
                                                     y=y, X=X, Sigma=Sigma) # additional arguments
grad_constrained_mapping_theoretical <- gradient_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(grad_constrained_mapping_numerical, c(grad_constrained_mapping_theoretical))

##################################################################
##                 loglik_hessian_unconstrained                 ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")

set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X))# global gene expression, as linear combination
# y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))

hessian_loglik_numerical <- numDeriv::hessian(loglik_multivariate, p,
                                              method="Richardson", method.args=list(eps=1e-12, r=4),
                                              y=y, X=X, Sigma=Sigma) # additional arguments
# jacobian_grad_numerical <- numDeriv::jacobian(gradient_loglik_unconstrained, p,
#                                               method="Richardson", method.args=list(eps=1e-12, r=4),
#                                               y=y, X=X, Sigma=Sigma)
# testthat::expect_equal(hessian_loglik_numerical, jacobian_grad_numerical)


hessian_loglik_theoretical <- hessian_loglik_unconstrained (p, y, X, Sigma)


set.seed(3)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
# y <- X %*% p
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))

hessian_loglik_theoretical <- hessian_loglik_unconstrained (p, y, X, Sigma)
numDeriv::hessian(loglik_multivariate, p,
                  method="Richardson", method.args=list(eps=1e-12, r=4),
                  y=y, X=X, Sigma=Sigma)

##################################################################
##                 hessian mapping                 ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")
p <- c(0.2, 0.3, 0.4, 0.1); theta <- inverse_mapping_function(p)

# check Hessian of the mapping function
hessian_mapping_theoretical <- hessian_mapping_function(theta)
hessian_mapping_numerical <- array(0, dim = c(length(theta), length(theta), length(p)))
for (i in 1:length(p)) {
  hessian_mapping_numerical[,,i] <- numDeriv::hessian(mapping_function_univariate, theta,
                                                      method="Richardson", method.args=list(eps=1e-6), index=i)
}
testthat::expect_equal(hessian_mapping_theoretical, hessian_mapping_numerical, tolerance = 1e-6) # pass



set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
theta <- inverse_mapping_function(p)
num_genes <- nrow(X); num_celltypes <- ncol(X)
# y <- X %*% p + rnorm(nrow(X))# global gene expression, as linear combination
y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))



hessian_constrained_mapping_numerical <- numDeriv::hessian(loglik_multivariate_constrained, theta,
                                                           method="Richardson", method.args=list(eps=1e-12, r=4),
                                                           y=y, X=X, Sigma=Sigma) # additional arguments
hessian_constrained_mapping_theoretical <- hessian_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(hessian_constrained_mapping_numerical, hessian_constrained_mapping_theoretical)
