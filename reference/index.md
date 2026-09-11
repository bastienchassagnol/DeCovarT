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
ILR-constrained log-likelihood (Helmert basis; Jacobian / Hessian of the
isometric logistic map; gradients and Hessians of ()). Additive
log-ratio maps are kept for the vignette appendix. Documented for `?` /
source inspection; see also
[`vignette("theory-decovart-generative-model")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md).

- [`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md)
  : Helmert contrast matrix for isometric log-ratio coordinates
- [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
  [`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
  : Isometric logistic transform (ILR coordinates to the simplex)
- [`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)
  : Jacobian of the isometric logistic map
- [`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md)
  : Hessian tensor of the isometric logistic map
- [`starting_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/starting_simplex.md)
  : Open-simplex start for ILR solvers
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

### Fixed-covariance GLS competitor

[`MASS::lm.gls()`](https://rdrr.io/pkg/MASS/man/lm.gls.html) wrapper
with a known residual covariance
([`fixed_gls_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/fixed_gls_covariance.md)
at a declared composition).

- [`deconvolute_ratios_gls()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_gls.md)
  : Generalised least squares with a fixed residual covariance
- [`fixed_gls_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/fixed_gls_covariance.md)
  : Fixed residual covariance for a GLS competitor

### Structure-aware covariance backends

Factorise (()=\_j p_j^2 \_j) by exploiting declared covariance structure
(block, band, sparse, diagonal-plus-low-rank) instead of a dense
Cholesky. Operators return the log-determinant and solves without a
materialised precision. See
[`vignette("theory-decovart-generative-model")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md).

- [`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  [`print(`*`<decovart_covariance>`*`)`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
  : Structure-aware covariance backend for
  \\\boldsymbol{\Sigma}(\boldsymbol{p})\\
- [`covariance_structure_from_graph_model()`](https://bastienchassagnol.github.io/DeCovarT/reference/covariance_structure_from_graph_model.md)
  : Recommend a covariance backend from a network topology
- [`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md)
  : Log-determinant of \\\boldsymbol{\Sigma}(\boldsymbol{p})\\
- [`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md)
  : Precision solve
  \\\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{B}\\
- [`sigma_quadform()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_quadform.md)
  : Mahalanobis quadratic form
  \\\mathbf{r}^{\mathsf{T}}\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{r}\\
- [`sigma_trace_precision_times()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_trace_precision_times.md)
  : Trace of the precision times a matrix,
  \\\operatorname{tr}(\boldsymbol{\Theta}(\boldsymbol{p})\mathbf{S})\\

### Fisher information and Wald intervals

Expected Fisher information of unconstrained proportions and ILR
delta-method covariance used by
[`vcov()`](https://rdrr.io/r/stats/vcov.html) /
[`confint()`](https://rdrr.io/r/stats/confint.html) on `decovart_fit`
objects (see `fit_decovart`).

- [`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
  : Expected Fisher information of unconstrained \\\boldsymbol{p}\\
- [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
  : Cramer–Rao / ILR delta-method covariance of \\\hat{\boldsymbol{p}}\\

### Likelihood-ratio and boundary inference

Profile likelihood, likelihood-ratio tests with chi-bar-square (Chernoff
/ Self–Liang) calibration on simplex faces, parametric bootstrap,
reference-sample bootstrap (donors, cells, or Dirichlet compositions),
and the boundary / multimodality diagnostics that a bare convergence
code cannot provide. See
[`vignette("theory-DeCovarT-MLE-properties")`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.md).

- [`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)
  : Restricted maximum likelihood with fixed cellular ratios
- [`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md)
  : Profile log-likelihood of one cellular ratio
- [`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)
  : Likelihood-ratio test for cellular ratios
- [`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md)
  : Profile-likelihood confidence intervals for cellular ratios
- [`chi_bar_square_pvalue()`](https://bastienchassagnol.github.io/DeCovarT/reference/chi_bar_square_pvalue.md)
  : Chi-bar-square tail probability for boundary likelihood-ratio tests
- [`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md)
  : Parametric bootstrap for DeCovarT proportions and boundary tests
- [`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md)
  : Reference-based bootstrap for signature and composition uncertainty
- [`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md)
  : Permutation-equivariance check for a labelled reference
- [`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md)
  : Boundary and stationarity diagnostics for one DeCovarT fit
- [`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md)
  : Refit DeCovarT from several random starts to probe multimodality

## Benchmark and evaluation

Quality metrics for estimated proportions, and the simulation benchmark
wrapper.

- [`compute_benchmark_metrics()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_benchmark_metrics.md)
  : Compute deconvolution benchmark metrics
- [`coverage_mc_interval()`](https://bastienchassagnol.github.io/DeCovarT/reference/coverage_mc_interval.md)
  : Monte Carlo interval for an empirical coverage probability
- [`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
  : Simulate bulk mixtures and benchmark deconvolution algorithms
- [`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
  : Parallel deconvolution of a bulk expression matrix
- [`pivot_mc_estimates()`](https://bastienchassagnol.github.io/DeCovarT/reference/pivot_mc_estimates.md)
  : Pivot Monte Carlo proportion estimates to a long table
- [`plot_mc_raincloud()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_raincloud.md)
  : Horizontal raincloud of Monte Carlo proportion estimates
- [`plot_mc_forest()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_forest.md)
  : Forest plot of ADEMP Monte Carlo summaries
- [`theme_decovart_facets()`](https://bastienchassagnol.github.io/DeCovarT/reference/theme_decovart_facets.md)
  : Faceted ggplot2 theme (black strips, panel border)
- [`algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/algorithm_similarity.md)
  : Algorithm-similarity correlation from a Monte Carlo benchmark
- [`plot_algorithm_similarity()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_algorithm_similarity.md)
  : Tile heatmap of algorithm-similarity correlations
- [`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md)
  : Faceted dot plot of several ADEMP metrics
- [`plot_correlation_Heatmap()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_correlation_Heatmap.md)
  : Plot deconvolution metric heatmaps
- [`write_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/write_simulation_artefacts.md)
  : Write split simulation artefacts (config, descriptors, theta,
  metrics)
- [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md)
  : Read split simulation artefacts and optionally reassemble a
  benchmark list
- [`slim_scenario_table()`](https://bastienchassagnol.github.io/DeCovarT/reference/slim_scenario_table.md)
  : Drop duplicated design columns from a scenario-tagged table
- [`plot_purified_density_2d()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_purified_density_2d.md)
  : 2-D density of purified Gaussians
- [`plot_bulk_convolution_density_2d()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bulk_convolution_density_2d.md)
  : 2-D density of the bulk convolution \\y\sim N(\mu p,\Sigma(p))\\
- [`gaussian_confidence_ellipse()`](https://bastienchassagnol.github.io/DeCovarT/reference/gaussian_confidence_ellipse.md)
  : Exact Gaussian probability ellipse (known mean and covariance)
- [`plot_bulk_loglik_surface_p()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bulk_loglik_surface_p.md)
  : Contour of the bulk log-likelihood on a \\(p_1,p_2)\\ lattice
- [`plot_bulk_loglik_ilr_profile()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bulk_loglik_ilr_profile.md)
  : Log-likelihood profile in the ILR coordinate
  \\\rho\in\mathbb{R}^{J-1}\\
- [`plot_bulk_loglik_rgl()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bulk_loglik_rgl.md)
  : rgl surface of the bulk log-likelihood on a \\(p_1,p_2)\\ lattice
- [`plot_bivariate_metric_tiles()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_bivariate_metric_tiles.md)
  : Tile heatmap of a bivariate metric on the (\\\rho_1\\,\\\rho_2\\)
  plane
- [`save_bivariate_metric_heatmaps()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_metric_heatmaps.md)
  : Write RMSE / MAE / Aitchison tile PDFs for the bivariate toy
- [`save_bivariate_similarity_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_similarity_book.md)
  : Multi-page PDF of clustered algorithm-similarity heatmaps
- [`save_bivariate_forest_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_forest_book.md)
  : Multi-page PDF of Wald forests at four correlation corners
- [`save_bivariate_raincloud_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_raincloud_book.md)
  : Multi-page PDF of rainclouds at four correlation corners
- [`save_bivariate_solver_dots_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_solver_dots_book.md)
  : 12-page PDF of RMSE/Aitchison solver dots on the correlation grid
- [`save_bivariate_bulk_density_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_bulk_density_book.md)
  : Multi-page PDF of bulk convolution densities
- [`save_bivariate_loglik_ilr_profile_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_loglik_ilr_profile_book.md)
  : Multi-page PDF of ILR log-likelihood profiles
  (\\\rho\in\mathbb{R}\\)
- [`save_bivariate_loglik_rgl_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_loglik_rgl_book.md)
  : Multi-page PDF (and optional HTML) of rgl bulk log-likelihood
  surfaces
- [`save_bivariate_loglik_surface_p_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_loglik_surface_p_book.md)
  : Multi-page PDF of bulk log-likelihood surfaces on \\(p_1,p_2)\\
- [`save_bivariate_purified_density_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_bivariate_purified_density_book.md)
  : Multi-page PDF of purified densities at four correlation corners

## Simulation

### Toy bulk mixtures

Simulate bulk mixtures as convolutions of multivariate Gaussians, and
hierarchical GRN moments with graph-constrained covariances.

- [`generate_mean_signature_matrix()`](https://bastienchassagnol.github.io/DeCovarT/reference/generate_mean_signature_matrix.md)
  : Generate mean profiles from a target Gram matrix
- [`equicorrelation_gram()`](https://bastienchassagnol.github.io/DeCovarT/reference/equicorrelation_gram.md)
  : Equicorrelation (constant-correlation) Gram matrix
- [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
  : Simulate bulk mixtures from a multivariate Gaussian convolution
- [`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
  : Simulate GRN first- and second-order moments
- [`describe_simulation_scenario()`](https://bastienchassagnol.github.io/DeCovarT/reference/describe_simulation_scenario.md)
  : Describe a Gaussian-convolution simulation scenario

## Statistical metrics

Simplex repair, Shannon entropy, MixSim overlap, Jeffreys divergence,
and glmnet gene scores for signature design.

- [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md)
  : Repair a numeric vector onto the unit simplex
- [`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)
  : Normalised Shannon entropy of a discrete distribution
- [`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md)
  : One-dominant composition with a target normalised Shannon entropy
- [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md)
  : Validate generative-model parameters \\\theta\\
- [`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
  : Average pairwise overlap of a Gaussian mixture
- [`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md)
  : MixSim-style pairwise overlap via stratified Sobol Monte Carlo
- [`spd_affine_invariant_distance()`](https://bastienchassagnol.github.io/DeCovarT/reference/spd_affine_invariant_distance.md)
  : Affine-invariant Riemannian (AIRM) distance between two SPD matrices
- [`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)
  : Mean pairwise AIRM distance of component covariances
- [`scale_covariance_array()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariance_array.md)
  : Scale component covariances, keeping precision zeros
- [`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md)
  : Calibrate a global covariance scale to a target MixSim BarOmega
- [`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)
  : Average pairwise Jeffreys divergence of a Gaussian mixture
- [`compute_glmnet_gene_scores()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_glmnet_gene_scores.md)
  : Gene scores from multinomial elastic-net cell-type classification
