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
unconstrained coordinates \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\ via
\\\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})\\
(Marquardt–Levenberg default; see other methods below and
[`vignette("softmax-alr-derivatives", package = "DeCovarT")`](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md)).

## Usage

``` r
deconvolute_ratios_cibersort(y, mean_signature_matrix)

deconvolute_ratios_lsfit(y, mean_signature_matrix)

deconvolute_ratios_rlm(y, mean_signature_matrix)

deconvolute_ratios_nnls(y, mean_signature_matrix)

deconvolute_ratios_deconrnaseq(y, mean_signature_matrix)

deconvolute_ratios_Marquardt_Levenberg(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200
)

deconvolute_ratios_simulated_annealing(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200
)

deconvolute_ratios_L_BFGS_B(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200
)

deconvolute_ratios_Newton_Raphson(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200
)

deconvolute_ratios_gradient_descent(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200
)
```

## Arguments

- y:

  Bulk expression vector \\\boldsymbol{y}\in\mathbb{R}^{G}\\ (one
  heterogeneous sample).

- mean_signature_matrix:

  Mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (columns = cell types; plug-in for latent profiles).

- Sigma:

  Array \\(\boldsymbol{\Sigma}\_j)\_{j=1}^{J}\in\mathcal{M}\_{G\times
  G\times J}\\ of cell-type covariances.

- epsilon, itmax:

  Absolute convergence tolerance and maximum number of iterations for
  the optimiser.

## Value

Named numeric vector \\\hat{\boldsymbol{p}}\\ on the simplex. Benchmark
metrics are computed by
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

- `deconvolute_ratios_lsfit()`: Ordinary least squares for
  \\\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}\\
  ([`stats::lsfit()`](https://rdrr.io/r/stats/lsfit.html)), following
  Abbas et al. (2009) ; estimates are projected back onto the simplex.

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

- `deconvolute_ratios_simulated_annealing()`: Simulated annealing on
  \\\boldsymbol{\rho}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html) with
  `method = "SANN"`).

- `deconvolute_ratios_L_BFGS_B()`: Box-constrained L-BFGS-B directly in
  \\\boldsymbol{p}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  `method = "L-BFGS-B"`).

- `deconvolute_ratios_Newton_Raphson()`: Newton–Raphson / `nlminb` on
  \\\boldsymbol{\rho}\\ using analytic gradient and Hessian
  ([`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)).

- `deconvolute_ratios_gradient_descent()`: BFGS quasi-Newton ascent on
  \\\boldsymbol{\rho}\\
  ([`stats::optim()`](https://rdrr.io/r/stats/optim.html)
  `method = "BFGS"`).

## References

Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z, Clark HF (2009).
“Deconvolution of Blood Microarray Data Identifies Cellular Activation
Patterns in Systemic Lupus Erythematosus.” *PloS One*, **4**.
[doi:10.1371/journal.pone.0006098](https://doi.org/10.1371/journal.pone.0006098)
.  
  
Monaco G, Lee B, Xu W, Mustafah S, Hwang YY, Carré C, Burdin N, Visan L,
Ceccarelli M, Poidinger M, Zippelius A, Pedro de Magalhães J, Larbi A
(2019). “RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute
Deconvolution of Human Immune Cell Types.” *Cell Reports*, **26**.
[doi:10.1016/j.celrep.2019.01.041](https://doi.org/10.1016/j.celrep.2019.01.041)
.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
