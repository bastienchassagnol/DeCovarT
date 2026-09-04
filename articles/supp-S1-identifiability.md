# Appendix S1 — Regular-case MLE properties and identifiability

This appendix is the numerical counterpart of [MLE properties and
asymptotic
inference](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.md).
That article states what the Gaussian-convolution MLE is entitled to
claim. The chunks here check those claims on small, fully specified
scenarios.

Two regimes are kept distinct:

1.  **Regular (interior) case.** True \boldsymbol{p} sits in the open
    simplex. Affine equivariance, small bulk perturbations, and random
    ILR starts probe estimator invariance and solver stability when the
    Wilks \chi^2 calibration of [interior
    inference](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-interior)
    applies.
2.  **Identifiability edge cases.** Collinear or identical means,
    boundary zeros, and near-flat likelihood ridges, matching
    [identifiability](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-identifiability)
    and [boundary
    calibration](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-boundary).

Heavy Monte Carlo factorial runs live in
`scripts/supp_S1_identifiability.R` (`N_REPLICATES=2` for a smoke test;
outputs in `output/supp_S1/`). The live chunks below use G\le 3 and are
meant to be rebuilt with the vignette.

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
```

## Regular case: interior of the simplex

The generating composition below is strictly positive. That is the
setting in which
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
Wald intervals are defined and the score equations of [the convolution
log-likelihood](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#eq-loglik)
are regular.

### Gene-wise affine equivariance

Let D=\mathrm{diag}(s_1,\ldots,s_G) with s_g\neq 0 and
\boldsymbol{m}\in\mathbb{R}^{G}. The same gene-wise map applied to bulk,
means and covariances, \boldsymbol{y}^\star\_{\cdot i} =
D^{-1}(\boldsymbol{y}\_{\cdot i}-\boldsymbol{m}), \quad
\boldsymbol{\mu}^\star =
D^{-1}(\boldsymbol{\mu}-\boldsymbol{m}\boldsymbol{1}\_J^\top), \quad
\boldsymbol{\Sigma}\_j^\star = D^{-1}\boldsymbol{\Sigma}\_j D^{-1},
preserves the linear convolution because
\boldsymbol{1}\_J^\top\boldsymbol{p}\_{\cdot i}=1. The Mahalanobis
residual is invariant, so
\hat{\boldsymbol{p}}(\boldsymbol{Y},\boldsymbol{\mu},\boldsymbol{\Sigma})
=\hat{\boldsymbol{p}}(\boldsymbol{Y}^\star,\boldsymbol{\mu}^\star,
\boldsymbol{\Sigma}^\star) in exact arithmetic. A z-score computed
**once per gene from \boldsymbol{\mu}** (not per cell type, not per bulk
sample) is that map; see [equivariance in the MLE
article](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-mle-properties).

A logarithm is not affine. It is concave, so Jensen’s inequality changes
first and second moments, and
\log(\boldsymbol{\mu}\boldsymbol{p})\neq(\log\boldsymbol{\mu})\boldsymbol{p}.
`CIBERSORT` therefore requires non-negative expression, no missing
values, and a non-log linear scale ([Newman et al.
2015](#ref-newmanRobustEnumerationCell2015)).

``` r

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

### Small perturbations of the bulk

With known interior \boldsymbol{p}, a small Gaussian perturbation of
\boldsymbol{y}\_{\cdot i} should move Marquardt–Levenberg estimates only
slightly. Three covariance *models* are compared:

1.  Full cell-type tensor \boldsymbol{\Sigma}\_j (DeCovarT convolution).
2.  Cell-type diagonal: off-diagonals of each \boldsymbol{\Sigma}\_j set
    to zero (gene-wise weighted convolution).
3.  Fixed global diagonal GLS:
    W=\operatorname{diag}\\\Sigma(p^{\star})\\ held **constant**, fitted
    by
    [`deconvolute_ratios_gls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_gls.md)
    ([`MASS::lm.gls`](https://rdrr.io/pkg/MASS/man/lm.gls.html)). Do
    **not** copy that W into every tensor slice: \sum_j p_j^2
    W=\\p\\\_2^2 W still depends on p.

``` r

set.seed(2)
p_true <- c(0.6, 0.4)
y0 <- matrix(drop(mu %*% p_true), ncol = 1, dimnames = list(genes, "s1"))
y_pert <- y0 + rnorm(3, sd = 0.2)

Sigma_diag <- Sigma
for (j in 1:2) {
  Sigma_diag[, , j] <- diag(diag(Sigma[, , j]))
}
w_gls <- fixed_gls_covariance(Sigma, p = p_true, diagonal = TRUE)

fit_spec <- function(y, S) {
  drop(coef(fit_decovart(
    mu,
    y,
    Sigma = S,
    method = "Marquardt-Levenberg",
    itmax = 80
  )))
}

est <- cbind(
  full = fit_spec(y_pert, Sigma),
  celltype_diagonal = fit_spec(y_pert, Sigma_diag),
  gls_fixed = deconvolute_ratios_gls(
    drop(y_pert),
    mu,
    w_gls
  )
)
rmse <- apply(est, 2, function(p_hat) sqrt(mean((p_hat - p_true)^2)))
list(estimates = est, rmse_to_truth = rmse)
#> $estimates
#>          full celltype_diagonal gls_fixed
#> ct1 0.6112459         0.6112459 0.6121588
#> ct2 0.3887541         0.3887541 0.3878412
#> 
#> $rmse_to_truth
#>              full celltype_diagonal         gls_fixed 
#>        0.01124588        0.01124588        0.01215877
```

### Random ILR starts

Marquardt–Levenberg and Newton–Raphson maximise the convolution
log-likelihood with analytic Hessians. That is cheaper per step than a
first-order quasi-Newton map, but a start near a simplex face can stall
on a plateau. Three Dirichlet(1,1) draws
(`starting_simplex(..., "dirichlet")`, seeds 11, 22 and 33) are passed
as `initial_p` on the same perturbed bulk. Concentration \alpha=1 is
uniform on the simplex; \alpha\>1 is centre-biased and \alpha\<1 is
face-biased. Sequence several draws (or
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md))
rather than a single random start. A QP warm start is
`initial_p = "qp"`. On this well-conditioned interior toy the
second-order estimates should agree; a spread of several percentage
points is a warning that the start has trapped the solver.

``` r

seeds <- c(11L, 22L, 33L)
starts <- lapply(seeds, function(s) {
  set.seed(s)
  starting_simplex(2L, "dirichlet", nms = cts)
})
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
#>        ct1        ct2 
#> 0.09319046 0.90680954 
#> 
#> $starts$seed_22
#>       ct1       ct2 
#> 0.2746408 0.7253592 
#> 
#> $starts$seed_33
#>       ct1       ct2 
#> 0.1712805 0.8287195 
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
dimension (G,J) belongs in [Appendix
S3](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S3-scaling.md).

## Identifiability edge cases

Standard deconvolution theory assumes that the mean signature
\boldsymbol{\mu} has full column rank. That is **sufficient** for
identifiability of \boldsymbol{p} but **not necessary**: the convolution
law also depends on
\boldsymbol{\Sigma}(\boldsymbol{p})=\sum\_{j}p\_{j}^{2}\boldsymbol{\Sigma}\_{j}
([proposition](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#prp-identifiability)).
Identical means with distinct covariances remain identifiable; the
negative control is identical means *and* proportional covariances.

### Collinear mean columns

[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
flags collinear mean columns as a numerical hazard even when the
covariance tensor separates the types.

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

### Factorial Monte Carlo (script)

The remaining three sub-scenarios are too large for a vignette rebuild.
Run `Rscript scripts/supp_S1_identifiability.R` (or `N_REPLICATES=2` for
a smoke test).

| Sub-scenario | Factor | Levels |
|----|----|----|
| **S1a** Boundary | p_3 = 0 | Wald / profile / bootstrap CI on J=3, G=10 |
| **S1b** Same means | Covariance separation \rho\_{12} | 0, 0.5, 0.9 (J=2, G=5) |
| **S1c** Multimodal | Gene–gene correlation \rho | -0.8, 0, 0.8 (low CLD, J=2, G=2) |

N = 1000 Monte Carlo replicates per sub-scenario in a full run (smoke
test: n = 2). Solvers: Marquardt–Levenberg for S1b/S1c; S1a also records
interval coverage near a simplex face.

``` r

N_REPLICATES <- as.integer(Sys.getenv("N_REPLICATES", unset = "1000"))
# See scripts/supp_S1_identifiability.R for config_1a / config_1b /
# config_1c and the run_simulation_benchmark() calls.
```

#### Expected findings

- **S1a (boundary):** NNLS clips to zero by construction. Wald intervals
  from
  [`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  are undefined on a face; use
  [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
  /
  [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)
  with chi-bar-square calibration ([boundary
  nulls](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-boundary)).
- **S1b (same means):** RMSE decreases as covariance separation grows.
  Below \rho\_{12}=0 (identical diagonal \boldsymbol{\Sigma}\_j) the
  pair is a negative control. DeCovarT’s advantage is largest when means
  are uninformative.
- **S1c (multimodal):** The fraction of multistarts recovering the
  global mode decreases as the mean columns approach each other;
  Marquardt–Levenberg should explore more of the basin than a single
  L-BFGS-B start.

### See also

- [MLE
  properties](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.md)
- [Appendix S2 — Covariance
  inversion](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S2-covariance-inversion.md)
- [Appendix S3 — Solver
  scalability](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S3-scaling.md)
- [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.md)

### References

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.
