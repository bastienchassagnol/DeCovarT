# Package index

## Main deconvolution function

Top-level function to deconvolve a matrix of bulk expression samples
using any combination of deconvolution algorithms in parallel.

- [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  : Parallel deconvolution of a bulk expression matrix

## DeCovarT algorithm

Covariance-aware optimisation of cellular ratios under a multivariate
Gaussian convolution model.

### Core optimiser

Marquardt–Levenberg deconvolution (exported entry point). Related first-
and second-order variants live in
`R/03_03_DeCovarT_estimate_ratios_frequentist.R`.

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

### Parametrisation

Additive logistic transform between the unconstrained parameter space
(^{J-1}) and the unit simplex (^{J-1}).

- [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
  [`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
  : Additive logistic transform (unconstrained coordinates to the
  simplex)

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

- [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
  : Simulate bulk mixtures from a Gaussian convolution
- [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
  : Simulate hierarchical GRN first- and second-order moments

## Utilities

Simplex projection, Shannon entropy, and average pairwise overlap used
when designing synthetic scenarios.

- [`enforce_identifiability()`](https://bastienchassagnol.github.io/DeCovarT/reference/enforce_identifiability.md)
  : Project estimated proportions onto the unit simplex
- [`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)
  : Normalised Shannon entropy of a discrete distribution
- [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
  : Average pairwise overlap of a Gaussian mixture
