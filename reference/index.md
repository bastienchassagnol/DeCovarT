# Package index

## Simplex chart

Shared by both generative models (`R/02_00_ilr_chart.R`): simplex
repair, the Helmert basis and the isometric log-ratio (ILR) chart with
its first and second derivatives.

- [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  : Repair a numeric vector onto the unit simplex
- [`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md)
  : Helmert contrast matrix for isometric log-ratio coordinates
- [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
  [`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
  : Isometric logistic transform (ILR coordinates to the simplex)
- [`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)
  : Jacobian of the isometric logistic map
- [`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md)
  : Hessian tensor of the isometric logistic map

## Population references: Gaussian convolution

Bulk covariance (\_j p_j^2_j) (article Section “Population-level
references”).

### Optimisation

Log-likelihood, score and Hessian, the solver family (Newton–Raphson by
default, with the PORT trust region of
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html)), starting values
and multistarts (`R/02_02_population_decovart_estimate_ratios.R`,
`R/02_03_population_decovart_fit.R`).

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
- [`deconvolute_ratios_cibersort()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_rlm()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_nnls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_Marquardt_Levenberg()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_simulated_annealing()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_Newton_Raphson()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  [`deconvolute_ratios_gradient_descent()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
  : DeCovarT MLE of cellular proportions for one bulk sample
- [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md)
  : Open-simplex start for ILR solvers
- [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  : Refit DeCovarT from several random starts to probe multimodality
- [`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
  : Unconstrained DeCovarT log-likelihood
- [`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md)
  : Gradient \\\nabla\_{\boldsymbol{p}}\ell\\ of the unconstrained
  log-likelihood
- [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md)
  : Hessian \\\mathbf{H}\\ of the unconstrained log-likelihood
- [`gradient_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_constrained.md)
  : Constrained gradient via the chain rule
- [`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md)
  : Constrained Hessian of \\\ell\circ\boldsymbol{\psi}\\

### Inference

Confidence intervals and tests (`R/03_01` to `R/03_03`): the
[`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md)
menu (ILR back-transformation by default, Wald, Firth, profile
likelihood, parametric and reference bootstraps), restricted fits on
simplex faces and chi-bar-square tests.

- [`confint_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_decovart.md)
  : Choose a confidence interval for a pooled DeCovarT composition
- [`confint_ilr_backtransform()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_ilr_backtransform.md)
  : Confidence intervals by mapping a Wald region through the ILR chart
- [`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
  : Expected Fisher information of unconstrained \\\boldsymbol{p}\\
- [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
  : Cramer–Rao / ILR delta-method covariance of \\\hat{\boldsymbol{p}}\\
- [`fit_firth_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_decovart.md)
  : Firth-penalised DeCovarT fit
- [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)
  : Profile-likelihood confidence intervals for cellular ratios
- [`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)
  : Maximum likelihood on the closed simplex, optionally on a face
- [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
  : Likelihood-ratio test for cellular ratios
- [`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md)
  : Chi-bar-square tail probability for boundary likelihood-ratio tests
- [`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md)
  : Boundary and stationarity diagnostics for one DeCovarT fit
- [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)
  : Parametric bootstrap for DeCovarT proportions and boundary tests
- [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  : Reference-based bootstrap for signature and composition uncertainty
- [`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
  : Permutation-equivariance check for a labelled reference

## Single-cell references: mixture of independent cells

Bulk covariance (*j p_j*{W,j}) with a mean scale (c) and a variance
scale () (article Section “Single-cell references”).

### Optimisation

Free-scale least-squares start, closed-form scales, and Newton–Raphson
on the profile log-likelihood
(`R/02_02_single_cell_decovart_estimate_ratios.R`).

- [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)
  : Single-cell deconvolution of one bulk sample
- [`single_cell_nuisance_scales()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_nuisance_scales.md)
  : Stage-i scales for one single-cell bulk sample
- [`loglik_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_single_cell.md)
  : Log-likelihood of one single-cell bulk sample
- [`profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md)
  [`gradient_profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md)
  [`hessian_profile_loglik_single_cell_ilr()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_single_cell_ilr.md)
  : Profile log-likelihood of one single-cell bulk sample in ILR
  coordinates

### Inference

The
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
menu, with intervals for the composition and for the two scales
(`R/03_04_single_cell_inference.R`).

- [`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
  : Confidence intervals for one single-cell bulk sample
- [`vcov_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_single_cell.md)
  : Observed information of one single-cell fit
- [`fit_firth_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_single_cell.md)
  : Firth fit of one single-cell bulk sample
- [`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md)
  : Profile-likelihood intervals for one single-cell bulk sample
- [`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md)
  : Maximum likelihood on a face of the single-cell model
- [`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)
  : Likelihood-ratio test for one single-cell composition
- [`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md)
  : Parametric bootstrap for one single-cell bulk sample
- [`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md)
  : Reference bootstrap for the single-cell mixture

## Generative models

### Population references

Mean signatures with prescribed pairwise cosines, graph-constrained
cell-type covariances, and bulk samples drawn from the Gaussian
convolution (`R/01_01_population_toy_simulation.R`,
`R/01_02_generate_synthetic_networks.R`).

- [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
  : Generate mean profiles from a target Gram matrix
- [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
  : Simulate GRN first- and second-order moments
- [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
  : Simulate bulk mixtures from a multivariate Gaussian convolution

### Single-cell references

Cell-size conventions, the four-scenario design, and bulk samples drawn
as sums of multinomial counts of cells
(`R/01_01_single_cell_toy_simulation.R`,
`R/01_02_single_cell_generate_synthetic_networks.R`).

- [`simulate_bulk_from_cells()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_from_cells.md)
  : Simulate a bulk sample as a sum of single cells
- [`rescale_single_cell_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/rescale_single_cell_moments.md)
  : Cell-size conventions for a single-cell reference
- [`single_cell_scenario_design()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_scenario_design.md)
  : Four single-cell simulation scenarios

## Benchmark

### Run

Fit every algorithm to every bulk sample, over a scenario grid
(`R/04_01_run_benchmark.R`, `R/04_02_run_simulation_benchmark.R`).

- [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  : Parallel deconvolution of a bulk expression matrix
- [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
  : Simulate bulk mixtures and benchmark deconvolution algorithms
- [`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md)
  : Describe a Gaussian-convolution simulation scenario

### Metrics

Global and per-cell-type performance scores, Monte Carlo coverage
(`R/04_03`, `R/04_04`), and scenario descriptors: entropy of the
composition, overlap and divergences between cell-type distributions
(`R/01_03` to `R/01_05`).

- [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  : Compute deconvolution benchmark metrics
- [`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md)
  : Monte Carlo interval for an empirical coverage probability
- [`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)
  : Normalised Shannon entropy of a discrete distribution
- [`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md)
  : One-dominant composition with a target normalised Shannon entropy
- [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
  : Average pairwise overlap of a Gaussian mixture
- [`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md)
  : MixSim-style pairwise overlap via stratified Sobol Monte Carlo
- [`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)
  : Mean pairwise AIRM distance of component covariances
- [`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)
  : Average pairwise Jeffreys divergence of a Gaussian mixture
- [`scale_covariance_array()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariance_array.md)
  : Scale component covariances, keeping precision zeros

### Visualisation

Monte Carlo plots and algorithm-similarity heatmaps of a benchmark
(`R/05_01_monte_carlo_plots.R`, `R/05_02_algorithm_similarity.R`).

- [`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md)
  : Pivot Monte Carlo proportion estimates to a long table
- [`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md)
  : Horizontal raincloud of Monte Carlo proportion estimates
- [`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md)
  : Forest plot of ADEMP Monte Carlo summaries
- [`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md)
  : Faceted dot plot of several ADEMP metrics
- [`algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/algorithm_similarity.md)
  : Algorithm-similarity correlation from a Monte Carlo benchmark
- [`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md)
  : Tile heatmap of algorithm-similarity correlations
- [`theme_decovart_facets()`](https://bastienchassagnol.github.io/DeCovarT/reference/theme_decovart_facets.md)
  : Faceted ggplot2 theme (black strips, panel border)
