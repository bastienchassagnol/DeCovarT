# DeCovarT appendix: affine equivariance, collinearity and perturbations

This appendix records four numerical facts that sit behind the
`decovart_fit` API: gene-wise affine equivariance of the MLE, failure
under collinear signatures, small bulk perturbations when the covariance
is full, cell-type diagonal, or global, and the effect of random ALR
starts on second-order solvers. CPM / log2 normalisation, optimiser
tolerance, and scaling with G, J and overlap are left to a later
benchmark issue.

## Gene-wise affine equivariance

Let D=\mathrm{diag}(s_1,\ldots,s_G) with s_g\neq 0 and
\boldsymbol{m}\in\mathbb{R}^{G}. The same gene-wise map applied to bulk,
means and covariances, \boldsymbol{y}^\star\_{\cdot i} =
D^{-1}(\boldsymbol{y}\_{\cdot i}-\boldsymbol{m}), \quad
\boldsymbol{\mu}^\star =
D^{-1}(\boldsymbol{\mu}-\boldsymbol{m}\boldsymbol{1}\_J^\top), \quad
\boldsymbol{\Sigma}\_j^\star = D^{-1}\boldsymbol{\Sigma}\_j D^{-1},
preserves the linear convolution because
\boldsymbol{1}\_J^\top\boldsymbol{p}\_{\cdot i}=1:
\boldsymbol{y}^\star\_{\cdot i}
=\boldsymbol{\mu}^\star\boldsymbol{p}\_{\cdot i}
+\boldsymbol{\varepsilon}^\star\_{\cdot i}. The Mahalanobis residual is
invariant, so
\hat{\boldsymbol{p}}(\boldsymbol{Y},\boldsymbol{\mu},\boldsymbol{\Sigma})
=\hat{\boldsymbol{p}}(\boldsymbol{Y}^\star,\boldsymbol{\mu}^\star,
\boldsymbol{\Sigma}^\star) in exact arithmetic. A z-score computed
**once per gene from \boldsymbol{\mu}** (not per cell type, not per bulk
sample) is that map. Sample-level covariates that shift gene-wise means
are included provided the same D and \boldsymbol{m} are applied to the
adjusted \boldsymbol{\mu}\_j(\boldsymbol{z}\_i).

A logarithm is not affine. It is concave, so Jensen’s inequality changes
first and second moments, and
\log(\boldsymbol{\mu}\boldsymbol{p})\neq(\log\boldsymbol{\mu})\boldsymbol{p}.
`CIBERSORT` therefore requires non-negative expression, no missing
values, and a non-log linear scale ([Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)). The logarithm is monotone
on (0,\infty), so gene-wise orderings are preserved, but the
Gaussian-convolution MLE of \boldsymbol{p} is not.

``` r

library(DeCovarT)

set.seed(1)
genes <- paste0("g", 1:3)
cts <- paste0("ct", 1:2)
mu <- matrix(
  c(20, 40, 25, 35, 30, 22),
  nrow = 3,
  dimnames = list(genes, cts)
)
Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
Sigma[, , 1] <- diag(c(1.0, 1.2, 0.8))
Sigma[, , 2] <- diag(c(1.5, 0.9, 1.1))
p_true <- c(0.65, 0.35)
y <- matrix(
  drop(mu %*% p_true + rnorm(3, sd = 0.05)),
  ncol = 1,
  dimnames = list(genes, "s1")
)

fit_raw <- fit_decovart(
  mu,
  y,
  Sigma = Sigma,
  method = "Marquardt-Levenberg",
  itmax = 80
)

centre <- c(2, -1, 0.5)
scale <- c(2, 0.5, 4)
D_inv <- 1 / scale
y_star <- (y - centre) / scale
mu_star <- (mu - centre) / scale
Sigma_star <- Sigma
for (j in seq_len(2)) {
  Sigma_star[, , j] <- Sigma[, , j] * tcrossprod(D_inv)
}
# Solvers accept the starred arrays directly (z-scores need not be
# non-negative). The public gate still requires non-negative raw counts.
p_star <- deconvolute_ratios_Marquardt_Levenberg(
  drop(y_star),
  mu_star,
  Sigma_star,
  itmax = 80
)
fit_z <- fit_decovart(
  mu,
  y,
  Sigma = Sigma,
  method = "Marquardt-Levenberg",
  itmax = 80,
  standardise = TRUE
)

rbind(
  raw = drop(coef(fit_raw)),
  general_affine = p_star,
  gene_wise_z = drop(coef(fit_z))
)
#>                      ct1       ct2
#> raw            0.6489623 0.3510377
#> general_affine 0.6489623 0.3510377
#> gene_wise_z    0.6489623 0.3510377
```

The two rows should agree up to optimiser noise: `standardise = TRUE` is
a gene-wise affine map. Cell-type-wise
[`scale()`](https://rdrr.io/r/base/scale.html) of \boldsymbol{\mu} would
estimate a different D per column and is not supported.

## Collinearity

If two signature columns are identical, or a parent profile is an affine
combination of children, \boldsymbol{\mu} is rank-deficient and
\boldsymbol{p} is not identifiable (RE2.4b).

``` r

mu2 <- mu
mu2[, 2] <- mu2[, 1]
tryCatch(
  fit_decovart(
    mu2,
    y,
    Sigma = Sigma,
    method = "Marquardt-Levenberg",
    itmax = 20
  ),
  warning = function(w) conditionMessage(w),
  error = function(e) conditionMessage(e)
)
#> [1] "Signature columns are collinear (rank 1 < J = 2); mixture proportions are not identifiable."
```

``` r

genes3 <- paste0("g", 1:3)
cts3 <- c("child_1", "child_2", "parent")
mu3 <- cbind(
  child_1 = c(20, 40, 15),
  child_2 = c(40, 20, 25),
  parent = 0.5 * c(20, 40, 15) + 0.5 * c(40, 20, 25)
)
rownames(mu3) <- genes3
Sigma3 <- array(0, dim = c(3, 3, 3), dimnames = list(genes3, genes3, cts3))
for (j in 1:3) {
  Sigma3[, , j] <- diag(3)
}
y3 <- matrix(mu3 %*% c(0.2, 0.3, 0.5), ncol = 1)
dimnames(y3) <- list(genes3, "s1")
qr(mu3)$rank
#> [1] 2
tryCatch(
  fit_decovart(
    mu3,
    y3,
    Sigma = Sigma3,
    method = "Marquardt-Levenberg",
    itmax = 20
  ),
  warning = function(w) conditionMessage(w),
  error = function(e) conditionMessage(e)
)
#> [1] "Signature columns are collinear (rank 2 < J = 3); mixture proportions are not identifiable."
```

Missing values are not a separate predictor / response policy (RE2.2):
quantified bulk RNA-seq matrices are complete, and the API errors on
`NA` / `NaN` / `Inf`.

## Small perturbations of the bulk

With known \boldsymbol{p}, a small Gaussian perturbation of
\boldsymbol{y}\_{\cdot i} should move Marquardt–Levenberg estimates only
slightly. Three covariance specifications are compared:

1.  Full cell-type tensor \boldsymbol{\Sigma}\_j.
2.  Cell-type diagonal: off-diagonals of each \boldsymbol{\Sigma}\_j set
    to zero (gene-wise weighted convolution).
3.  Global competitor: every slice replaced by the convolution
    covariance \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
    p_j^{2}\boldsymbol{\Sigma}\_j evaluated at the known mixing weights.

``` r

set.seed(2)
p_true <- c(0.6, 0.4)
y0 <- matrix(drop(mu %*% p_true), ncol = 1, dimnames = list(genes, "s1"))
y_pert <- y0 + rnorm(3, sd = 0.2)

Sigma_diag <- Sigma
for (j in 1:2) {
  Sigma_diag[, , j] <- diag(diag(Sigma[, , j]))
}
sigma_bar <- .compute_global_variance(p_true, Sigma)
Sigma_global <- Sigma
Sigma_global[, , 1] <- sigma_bar
Sigma_global[, , 2] <- sigma_bar

fit_spec <- function(y, S) {
  drop(coef(fit_decovart(
    mu,
    y,
    Sigma = S,
    method = "Marquardt-Levenberg",
    itmax = 80
  )))
}

specs <- list(
  full = Sigma,
  celltype_diagonal = Sigma_diag,
  global_weighted = Sigma_global
)
est <- vapply(specs, function(S) fit_spec(y_pert, S), numeric(2))
rmse <- apply(est, 2, function(p_hat) sqrt(mean((p_hat - p_true)^2)))
list(estimates = est, rmse_to_truth = rmse)
#> $estimates
#>          full celltype_diagonal global_weighted
#> ct1 0.6112459         0.6112459       0.6115984
#> ct2 0.3887541         0.3887541       0.3884016
#> 
#> $rmse_to_truth
#>              full celltype_diagonal   global_weighted 
#>        0.01124588        0.01124588        0.01159839
```

### Random ALR starts

Marquardt–Levenberg and Newton–Raphson maximise the convolution
log-likelihood with analytic Hessians. That is cheaper per step than a
first-order quasi-Newton map, but a start near a simplex face can stall
on a plateau. Three Dirichlet(1,1) draws of \boldsymbol{p}^{(0)} (seeds
11, 22 and 33) are passed as `initial_p` on the same perturbed bulk. On
this well-conditioned toy the second-order estimates should agree; a
spread of several percentage points is a warning that the start has
trapped the solver.

``` r

dirichlet_start <- function(seed) {
  set.seed(seed)
  g <- stats::rgamma(2, shape = 1, rate = 1)
  g / sum(g)
}
seeds <- c(11L, 22L, 33L)
starts <- lapply(seeds, dirichlet_start)
names(starts) <- paste0("seed_", seeds)

fit_from_start <- function(fun, p0) {
  fun(drop(y_pert), mu, Sigma, itmax = 80, initial_p = p0)
}

ml <- vapply(
  starts,
  function(p0) fit_from_start(deconvolute_ratios_Marquardt_Levenberg, p0),
  numeric(2)
)
nr <- vapply(
  starts,
  function(p0) fit_from_start(deconvolute_ratios_Newton_Raphson, p0),
  numeric(2)
)
list(
  starts = starts,
  marquardt = ml,
  newton_raphson = nr
)
#> $starts
#> $starts$seed_11
#> [1] 0.09319046 0.90680954
#> 
#> $starts$seed_22
#> [1] 0.2746408 0.7253592
#> 
#> $starts$seed_33
#> [1] 0.1712805 0.8287195
#> 
#> 
#> $marquardt
#>      seed_11   seed_22   seed_33
#> ct1 0.611246 0.6112463 0.6112459
#> ct2 0.388754 0.3887537 0.3887541
#> 
#> $newton_raphson
#>      seed_11  seed_22  seed_33
#> ct1 0.611246 0.611246 0.611246
#> ct2 0.388754 0.388754 0.388754
```

A fuller comparison of normalisation (CPM, log2), tolerance, and
dimension (G,J) belongs in a later issue.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.
