# Simulating semi-synthetic pseudo-bulk mixtures for benchmarking

## Overview

This vignette documents the synthetic first- and second-order moment
generator
[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
and how its outputs feed bulk mixture simulation. Mean profiles are
constructed so that pairwise collinearity and centroid separation can be
dialled independently, following the AutoGeneS multi-objective rationale
of Aliee, Hananeh and Theis, Fabian J. ([2021](#ref-aliee_theis21)):
minimise correlation (here, cosine similarity) between cell-type columns
of \\\boldsymbol{\mu}\\ while maximising their pairwise Euclidean
distances. Covariances share a graph-constrained precision
\\\boldsymbol{\Omega}\\, mapped to
\\\boldsymbol{\Sigma}\_j=\boldsymbol{\Omega}^{-1}\\.

## Pipeline API

The exported entry point orchestrates internal helpers and returns
moments that can be passed to
[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)
and then to any deconvolution routine exposed by
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md).

``` mermaid
flowchart TD
  subgraph means ["Mean profiles (AutoGeneS-inspired)"]
    P1["target_cosine ρ / mean_scale s"] --> GM["generate_mean_signature_matrix()"]
    GM --> MU["μ"]
    MU --> OBJ["compute_mean_profile_objectives()"]
    OBJ --> SCORE["mean |cosine| and Euclidean distances"]
  end

  subgraph network ["Graph-constrained second order"]
    P2["graph_model / graph_params"] --> SK["generate_random_network_skeleton()"]
    SK --> A["A binary adjacency"]
    A --> BP["build_normalised_precision(u, v)"]
    BP --> OM["Omega positive definite"]
    OM --> CV["build_shared_covariance_array()"]
    CV --> SIG["Sigma_j = Omega inverse"]
  end

  MU --> MAIN["simulate_hierarchical_grn_moments()"]
  SIG --> MAIN
  SCORE --> MAIN
  MAIN --> MOM["moments: mu, Sigma, Omega, A"]
  MOM --> BULK["simulate_bulk_mixture()"]
  BULK --> Y["Y bulk matrix"]
  Y --> DEC["deconvolute_ratios()"]
```

``` r

library(DeCovarT)

set.seed(42)
moments <- simulate_hierarchical_grn_moments(
  n_genes = 40L,
  n_celltypes = 3L,
  mean_scale = 10,
  target_cosine = 0.1,
  precision_shift = 0.1,
  precision_scale = 0.3,
  graph_model = "power_law",
  graph_params = list(power = 1, edges_per_node = 2)
)

moments$objectives
#> $mean_abs_cosine
#> [1] 0.4144044
#> 
#> $sum_euclidean_distance
#> [1] 32.46647

bulk <- simulate_bulk_mixture(
  signature_matrix = moments$mean_profiles,
  Sigma = moments$covariance_matrices,
  p = c(0.5, 0.3, 0.2),
  n = 20
)
dim(bulk$Y)
#> [1] 40 20
```

Selecting genes (or, here, designing \\\boldsymbol{\mu}\\) by cosine
alone is insufficient for stable deconvolution ([Aliee, Hananeh and
Theis, Fabian J. 2021](#ref-aliee_theis21)): low Euclidean separation
makes ratio estimates noise-sensitive, and collinear centroids inflate
variance of \\\hat{\boldsymbol{p}}\\. The two knobs `target_cosine` and
`mean_scale` map onto that trade-off without running a genetic algorithm
at simulation time.

## Mathematical roles of the generators

### Mean signature: `generate_mean_signature_matrix()`

Each column \\\boldsymbol{\mu}\_{\cdot j}\\ is a normalised blend of a
shared direction \\\boldsymbol{u}\\ (uniform over genes) and a private
marker block \\\boldsymbol{v}\_{j}\\ (high on a gene partition for type
\\j\\, low elsewhere):

\\ \boldsymbol{\mu}\_{\cdot j} = s\\ \frac{ \sqrt{\rho}\\\boldsymbol{u}
+\sqrt{1-\rho}\\\boldsymbol{v}\_{j} }{ \bigl\\
\sqrt{\rho}\\\boldsymbol{u} +\sqrt{1-\rho}\\\boldsymbol{v}\_{j}
\bigr\\\_2 }, \qquad j=1,\ldots,J. \\

Here \\\rho\in\[0,1\]\\ is `target_cosine` and \\s\>0\\ is `mean_scale`.
Small \\\rho\\ yields complementary centroids; \\\rho\to 1\\ recovers a
shared direction (cosine one, vanishing Euclidean gap). Small \\s\\
places centroids near the origin and shrinks pairwise distances while
preserving angles.

### Objectives: `compute_mean_profile_objectives()`

For columns of \\\boldsymbol{\mu}\\,

\\ \overline{\lvert\cos\rvert} = \frac{2}{J(J-1)} \sum\_{1\le j\<k\le J}
\Biggl\| \frac{ \boldsymbol{\mu}\_{\cdot j}^{\mathsf{T}}
\boldsymbol{\mu}\_{\cdot k} }{ \\\boldsymbol{\mu}\_{\cdot j}\\\_2\\
\\\boldsymbol{\mu}\_{\cdot k}\\\_2 } \Biggr\|, \\

\\ D\_{\mathrm{Euc}} = \sum\_{1\le j\<k\le J} \\\boldsymbol{\mu}\_{\cdot
j}-\boldsymbol{\mu}\_{\cdot k}\\\_2. \\

AutoGeneS treats the first as a quantity to minimise and the second as a
quantity to maximise ([Aliee, Hananeh and Theis, Fabian J.
2021](#ref-aliee_theis21)). Cosine is preferred to Pearson correlation
so that collinear but unequally scaled profiles remain penalised.

### Network skeleton: `generate_random_network_skeleton()`

A binary adjacency \\\boldsymbol{A}\in\\0,1\\^{G\times G}\\ is drawn
from either a preferential-attachment (power-law / hub) model or a
stochastic block model (community structure). Both return an undirected
graph with zero diagonal; they encode presence of partial-correlation
edges, not their strengths.

### Positive-definite precision: `build_normalised_precision()`

Edge weights enter through an affine spectral shift of
\\\boldsymbol{A}\\, matching the construction used in the manuscript:

\\ \tilde{\boldsymbol{\Omega}} = \boldsymbol{A}\\v, \qquad
\boldsymbol{\Omega} = \tilde{\boldsymbol{\Omega}} + \bigl(
\lvert\lambda\_{\min}(\tilde{\boldsymbol{\Omega}})\rvert + u \bigr)
\mathbf{I}, \\

with off-diagonal scale \\v\\ (`precision_scale`) and diagonal cushion
\\u\\ (`precision_shift`). Larger \\v\\ strengthens interactions and
spreads the spectrum; larger \\u\\ improves conditioning. The map
\\Y=aX+bI\\ leaves eigenvectors unchanged and transforms eigenvalues as
\\\lambda_i(Y)=a\lambda_i(X)+b\\; choosing
\\b=\lvert\lambda\_{\min}(\tilde{\boldsymbol{\Omega}})\rvert+u\\ with
\\a=1\\ therefore guarantees \\\boldsymbol{\Omega}\succ 0\\ while
preserving the sparsity pattern of \\\boldsymbol{A}\\ off the diagonal.

### Shared covariances: `build_shared_covariance_array()`

Cell types share the same second-order structure,

\\ \boldsymbol{\Sigma}\_j = \boldsymbol{\Omega}^{-1}, \qquad
j=1,\ldots,J, \\

stacked as an array in \\\mathcal{M}\_{G\times G\times J}\\. Downstream
bulk simulation uses the Gaussian convolution
\\\boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p}, \sum_j
p_j^{2}\boldsymbol{\Sigma}\_j)\\.

## Perspectives: a single geometric score for \\\boldsymbol{\mu}\\

The present API exposes two AutoGeneS-style objectives. A natural
refinement is to replace them by one scalar that jointly rewards large
column norms (and hence Euclidean separation) and penalises alignment.

Two candidates are immediate. The Gram determinant

\\ V(\boldsymbol{\mu}) = \sqrt{ \det\bigl(
\boldsymbol{\mu}^{\mathsf{T}}\boldsymbol{\mu} \bigr) } \\

is the volume of the parallelepiped spanned by the columns of
\\\boldsymbol{\mu}\\. It vanishes under exact collinearity and grows
when columns lengthen or become more orthogonal—precisely the desired
MOO trade-off in one number.

Equivalently, the reciprocal condition number

\\ \kappa_2(\boldsymbol{\mu})^{-1} =
\frac{\sigma\_{\min}(\boldsymbol{\mu})}{\sigma\_{\max}(\boldsymbol{\mu})},
\\

available in R via `kappa(mean_profiles, exact = TRUE)` (or
`1 / kappa(...)`), measures how ill-posed linear recovery of
\\\boldsymbol{p}\\ from
\\\boldsymbol{y}\approx\boldsymbol{\mu}\boldsymbol{p}\\ is. Maximising
\\\kappa_2(\boldsymbol{\mu})^{-1}\\ (or minimising
\\\kappa_2(\boldsymbol{\mu})\\) collapses cosine and distance into a
deconvolution-centric criterion without a Pareto front.

A practical next step is to expose an optional
`mean_objective = c("autogenes", "volume", "condition")` in
[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md),
keep the constructive \\(\rho,s)\\ generator for controllable scenarios,
and report the chosen scalar alongside the existing pairwise
diagnostics.

## References

Aliee, Hananeh and Theis, Fabian J. 2021. ‘AutoGeneS: Automatic Gene
Selection Using Multi-Objective Optimization for RNA-seq Deconvolution’.
*Cell Systems*, ahead of print.
<https://doi.org/10.1016/j.cels.2021.05.006>.
