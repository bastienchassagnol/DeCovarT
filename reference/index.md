# Package index

## Main deconvolution function

Top-level function to deconvolve a matrix of bulk expression samples
using any combination of deconvolution algorithms in parallel.

- [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  : Parallel deconvolution of a bulk expression matrix

## DeCovarT algorithm

Covariance-aware optimisation of cellular ratios under a multivariate
Gaussian convolution model.

### Score equations and simplex maps

Analytic first- and second-order helpers for the unconstrained and
ALR-constrained log-likelihood (Jacobian / Hessian of the additive
logistic map; gradients and Hessians of ()). Documented for `?` / source
inspection; see also
[`vignette("generative-model-derivatives")`](https://bastienchassagnol.github.io/DeCovarT/articles/generative-model-derivatives.md).

- [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
  [`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
  : Additive logistic transform (unconstrained coordinates to the
  simplex)
- [`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md)
  : Jacobian \\\mathbf{J}\_{\boldsymbol{\psi}}\\ of the additive
  logistic map
- [`hessian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_additive_logistic.md)
  : Second derivatives (Hessian tensor) of the additive logistic map
- [`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md)
  : Gradient \\\nabla\_{\boldsymbol{p}}\ell\\ of the unconstrained
  log-likelihood
- [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md)
  : Hessian \\\mathbf{H}\\ of the unconstrained log-likelihood
- [`gradient_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_constrained.md)
  : Constrained gradient via the chain rule
- [`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md)
  : Constrained Hessian of \\\ell\circ\boldsymbol{\psi}\\

### Core optimiser

Marquardt–Levenberg deconvolution (exported entry point). Related first-
and second-order variants live in
`R/03_03_DeCovarT_estimate_ratios_frequentist.R`. The S3 model wrapper
is
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).

- [`deconvolute_ratios_cibersort()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_lsfit()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_rlm()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_nnls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_simulated_annealing()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_gradient_descent()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  : DeCovarT MLE of cellular proportions for one bulk sample
- [`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`print(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`summary(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`print(`*`<summary.decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`coef(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`fitted(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`residuals(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`vcov(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`nobs(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`confint(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  [`plot(`*`<decovart_fit>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
  : Fit the DeCovarT Gaussian-convolution model

### Fisher information and Wald intervals

Expected Fisher information of unconstrained proportions and ALR
delta-method covariance used by
[`vcov()`](https://rdrr.io/r/stats/vcov.html) /
[`confint()`](https://rdrr.io/r/stats/confint.html) on `decovart_fit`
objects (see `fit_decovart`).

- [`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
  : Expected Fisher information of unconstrained \\\boldsymbol{p}\\
- [`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md)
  : Cramer–Rao / ALR delta-method covariance of \\\hat{\boldsymbol{p}}\\

## Benchmark and evaluation

Quality metrics for estimated proportions, and the bivariate toy
convolution benchmark.

- [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  : Compute summary metrics for estimated proportions
- [`benchmark_bivariate_gaussian_convolutions()`](https://bastienchassagnol.github.io/DeCovarT/reference/benchmark_bivariate_gaussian_convolutions.md)
  : Benchmark bivariate Gaussian convolutions
- [`plot_correlation_Heatmap()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_correlation_Heatmap.md)
  : Plot deconvolution metric heatmaps

## Simulation

### Toy bulk mixtures

Simulate bulk mixtures as convolutions of multivariate Gaussians, and
hierarchical GRN moments with graph-constrained covariances.

- [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
  : Generate mean profiles with a target pairwise cosine
- [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
  : Simulate bulk mixtures from a multivariate Gaussian convolution
- [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
  : Simulate GRN first- and second-order moments

## Statistical metrics

Simplex repair, Shannon entropy, MixSim overlap, Jeffreys divergence,
and glmnet gene scores for signature design.

- [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  : Repair a numeric vector onto the unit simplex
- [`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)
  : Normalised Shannon entropy of a discrete distribution
- [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md)
  : Validate generative-model parameters \\\theta\\
- [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
  : Average pairwise overlap of a Gaussian mixture
- [`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)
  : Average pairwise Jeffreys divergence of a Gaussian mixture
- [`compute_glmnet_gene_scores()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_glmnet_gene_scores.md)
  : Gene scores from multinomial elastic-net cell-type classification
