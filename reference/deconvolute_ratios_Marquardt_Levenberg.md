# DeCovarT MLE of cellular proportions for one bulk sample

Estimates \\\hat{\boldsymbol{p}}=\arg\max\_{\boldsymbol{p}}
\ell\_{\boldsymbol{y}\\\|\\\boldsymbol{\zeta}}(\boldsymbol{p})\\ under
the Gaussian convolution model \$\$
\boldsymbol{y}\\\|\\(\boldsymbol{\zeta},\boldsymbol{p})
\sim\mathcal{N}\_{G}\\\bigl(\boldsymbol{\mu}\boldsymbol{p},\\
\boldsymbol{\Sigma}(\boldsymbol{p})\bigr), \qquad
\boldsymbol{\Sigma}(\boldsymbol{p})=\sum\_{j=1}^{J}p_j^{2}\boldsymbol{\Sigma}\_j,
\$\$ subject to the simplex constraint
\\\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1\\,
\\\boldsymbol{p}\ge\mathbf{0}\\. Optimisation is performed in
unconstrained ILR coordinates \\\boldsymbol{z}\in\mathbb{R}^{J-1}\\ via
\\\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{z})\\
(Marquardt–Levenberg default; see other methods below and
[`vignette("theory-decovart-generative-model", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md)).

## Usage

``` r
deconvolute_ratios_cibersort(y, mean_signature_matrix)

deconvolute_ratios_rlm(y, mean_signature_matrix)

deconvolute_ratios_nnls(y, mean_signature_matrix)

deconvolute_ratios_deconrnaseq(y, mean_signature_matrix)

deconvolute_ratios_Marquardt_Levenberg(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
)

deconvolute_ratios_simulated_annealing(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  initial_p = NULL,
  dirichlet_alpha = 1
)

deconvolute_ratios_L_BFGS_B(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
)

deconvolute_ratios_Newton_Raphson(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL,
  dirichlet_alpha = 1
)

deconvolute_ratios_gradient_descent(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  initial_p = NULL,
  dirichlet_alpha = 1
)
```

## Arguments

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (plug-in means).

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

- epsilon, itmax:

  Absolute convergence tolerance and maximum number of iterations for
  the optimiser (same roles as `reltol` / `maxit` in
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html)).

- return_model:

  If `TRUE`, return a named list with coefficients, ILR coordinates,
  log-likelihood and optimiser diagnostics instead of the proportion
  vector.

- initial_p:

  Optional start: `NULL` / `"barycentre"` (default equi-balanced), a
  length-\\J\\ numeric vector, `"dirichlet"`, or `"qp"` (mean-only
  simplex QP). See
  [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md).
  Starts on a simplex face are nudged into the interior so the ILR map
  is defined (ILR methods) and so L-BFGS-B does not start on a
  degenerate \\\boldsymbol{\Sigma}(\boldsymbol{p})\\.

- dirichlet_alpha:

  Dirichlet concentration when `initial_p = "dirichlet"` (default `1`:
  uniform on the simplex; \\\alpha\>1\\ centre-biased; \\\alpha\<1\\
  face-biased). Independent draws are independent restarts; use
  [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  to sequence several.

## Value

Named numeric vector \\\hat{\boldsymbol{p}}\\ on the simplex (ILR
methods), or that list when `return_model = TRUE`. Benchmark metrics are
computed by
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

## Details

**Plug-in signature.** Argument `mean_signature_matrix` is the mean
\\\boldsymbol{\mu}\\, used as a proxy for the unobserved cell-type
profiles \\\boldsymbol{x}\_{\cdot j}\\. This is the frequentist plug-in;
recovering sample-specific latents is a Bayesian / MAP problem
([`.map_gaussian_convolution()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-map_gaussian_convolution.md)).

## Functions

- `deconvolute_ratios_cibersort()`: Linear baseline
  \\\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}\\ via
  nu-SVR (CIBERSORT-style); no covariance prior is used.

- `deconvolute_ratios_rlm()`: Robust linear model
  \\\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}\\
  ([`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html)), as in Monaco
  et al. (2019) .

- `deconvolute_ratios_nnls()`: Non-negative least squares for
  \\\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}\\
  ([`nnls::nnls()`](https://rdrr.io/pkg/nnls/man/nnls.html)), then
  simplex projection.

- `deconvolute_ratios_deconrnaseq()`: Equality- and
  inequality-constrained least squares on the simplex
  ([`limSolve::lsei()`](https://rdrr.io/pkg/limSolve/man/lsei.html)), in
  the spirit of `deconRNASeq`.

- `deconvolute_ratios_simulated_annealing()`: Simulated annealing on ILR
  coordinates \\\boldsymbol{z}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html) with
  `method = "SANN"`).

- `deconvolute_ratios_L_BFGS_B()`: Box-constrained L-BFGS-B directly in
  \\\boldsymbol{p}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  `method = "L-BFGS-B"`). The box keeps each coordinate in \\\[0,1\]\\;
  the returned vector is closed by \\p/\sum p\\ (no
  [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  clipping).

- `deconvolute_ratios_Newton_Raphson()`: Newton–Raphson / `nlminb` on
  ILR coordinates \\\boldsymbol{z}\\ using analytic gradient and Hessian
  ([`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)).

- `deconvolute_ratios_gradient_descent()`: BFGS quasi-Newton ascent on
  ILR coordinates \\\boldsymbol{z}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  `method = "BFGS"`).

## References

Monaco G, Lee B, Xu W, Mustafah S, Hwang YY, Carré C, Burdin N, Visan L,
Ceccarelli M, Poidinger M, Zippelius A, Pedro de Magalhães J, Larbi A
(2019). “RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute
Deconvolution of Human Immune Cell Types.” *Cell Reports*, **26**.
[doi:10.1016/j.celrep.2019.01.041](https://doi.org/10.1016/j.celrep.2019.01.041)
.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md),
[`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md)

## Examples

``` r
set.seed(1)
genes <- paste0("g", 1:2)
cts <- paste0("ct", 1:2)
mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
Sigma <- array(
  c(1, 0, 0, 1, 1, 0, 0, 1),
  dim = c(2, 2, 2),
  dimnames = list(genes, genes, cts)
)
y <- drop(mu %*% c(0.6, 0.4) + rnorm(2, sd = 0.1))
deconvolute_ratios_Marquardt_Levenberg(y, mu, Sigma, itmax = 50)
#>       ct1       ct2 
#> 0.5805233 0.4194767 
```
