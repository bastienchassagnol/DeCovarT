# Appendix S2 — Covariance Inversion Strategies

``` r

library(DeCovarT)
```

> **Script:** `scripts/supp_S2_covariance_inversion.R` — benchmarks
> matrix inversion strategies. Outputs in `output/supp_S2/`.

------------------------------------------------------------------------

## 1️⃣ Generative Model

### Motivation

The DeCovarT log-likelihood requires repeated solves of the form
\boldsymbol{\Sigma}\_j^{-1} \boldsymbol{v}. Three inversion strategies
are compared for computational efficiency and numerical stability:

| Strategy | When optimal |
|----|----|
| **Dense Cholesky** | Small G or near-dense \boldsymbol{\Sigma}\_j |
| **Banded Cholesky** | \boldsymbol{\Sigma}\_j from a banded graph (local correlations only) |
| **Block-diagonal** | \boldsymbol{\Sigma}\_j with known block structure (connected components) |
| **Woodbury / Schur** | \boldsymbol{\Sigma}\_j = \boldsymbol{D} + \boldsymbol{U}\boldsymbol{C}\boldsymbol{U}^\top low-rank structure |

The benchmark measures **wall-clock time** and **memory** (peak RSS) as
functions of G and network structure. The same script first verifies
each structured backend (`block`, `band`, `sparse`, `diag_lowrank`)
against a dense Cholesky to machine precision on a small G.

### Factorial design

| Factor | Levels |
|----|----|
| G (genes) | 20, 50, 100, 250, 500 |
| Network structure | dense, banded (bandwidth 3), block-diagonal (4 blocks), low-rank (r = 5) |
| Backend | `dense_cholesky`, `banded_cholesky`, `block_diagonal`, `woodbury` |
| J (cell types) | 2, 3 |

N\_{\text{timing}} = 50 repeated solves per condition (median reported).

``` r

# -----------------------------------------------------------------
# SECTION 1: GENERATIVE MODEL
# Build precision matrices for each (G, structure) combination
# -----------------------------------------------------------------
SEED <- 20260807L

G_grid     <- c(20L, 50L, 100L, 250L, 500L)
structures <- c("dense", "banded", "block_diagonal", "low_rank")

build_precision <- function(G, structure, seed = SEED) {
  set.seed(seed)
  switch(structure,
    dense        = build_normalised_precision(
      generate_random_network_skeleton(G, model = "erdos_renyi", seed = seed)
    ),
    banded       = {
      adj <- bandSparse(G, G, k = -3:3, symmetric = TRUE,
                        diagonals = replicate(7, rep(1, G), simplify = FALSE))
      diag(adj) <- 0
      build_normalised_precision(as.matrix(adj))
    },
    block_diagonal = {
      block_size <- G %/% 4L
      adj_block  <- Matrix::bdiag(
        replicate(4L, {
          generate_random_network_skeleton(block_size, model = "erdos_renyi",
                                          seed = seed)
        }, simplify = FALSE)
      )
      build_normalised_precision(as.matrix(adj_block))
    },
    low_rank     = {
      # D + U C U^T structure (rank 5 perturbation)
      D_vals <- stats::runif(G, 0.5, 2)
      U_mat  <- matrix(rnorm(G * 5L), G, 5L)
      diag(D_vals) + tcrossprod(U_mat) * 0.1
    }
  )
}

precision_grids <- tidyr::expand_grid(G = G_grid, structure = structures) |>
  dplyr::mutate(
    Omega = purrr::map2(G, structure, build_precision)
  )
```

------------------------------------------------------------------------

## 2️⃣ Inference

### Timing benchmark

``` r

# -----------------------------------------------------------------
# SECTION 2: INFERENCE / TIMING
# Compare solve(Omega, v) using different backends
# -----------------------------------------------------------------
N_TIMING <- 50L
dir.create("output/supp_S2", recursive = TRUE, showWarnings = FALSE)

benchmark_results <- precision_grids |>
  dplyr::mutate(timing = purrr::map2(G, Omega, \(g, Om) {
    v <- rnorm(g)
    backends <- list(
      dense_cholesky  = \() solve(Om, v),
      block_diagonal  = \() {
        chol_Om <- chol(Om)
        backsolve(chol_Om, forwardsolve(t(chol_Om), v))
      }
    )
    purrr::map_dfr(backends, \(fn) {
      times <- replicate(N_TIMING, system.time(fn())[["elapsed"]])
      list(mean_ms = mean(times) * 1e3, sd_ms = sd(times) * 1e3)
    }, .id = "backend")
  })) |>
  tidyr::unnest(timing)

saveRDS(benchmark_results, "output/supp_S2/inversion_timing.rds")
```

------------------------------------------------------------------------

## 3️⃣ Visualisations

The four panels below split into two questions.
[Section 3.1](#sec-s2-inversion-panels) asks how the **linear algebra**
behind \boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\boldsymbol{v} scales
with G and graph sparsity. [Section 3.2](#sec-s2-solver-panels) asks how
the **iterative MLE solvers** spend that linear algebra: iterations to a
tolerance, wall-clock time, memory, simplex / KKT feasibility, and
agreement of analytic derivatives with Richardson extrapolation in
`numDeriv` ([Gilbert and Varadhan 2019](#ref-R-numDeriv)). The layout
follows the DICEPro supplementary solver and penalty panels (their
Figures S1, S2 and S5): log-scaled trajectories for cost, then a
separate check that first-order information and constraints are
numerically sound. The G–sparsity–time surface is the analogue of VIMixR
supplementary Figure S4.

| Output | Description |
|----|----|
| `output/supp_S2/s2_inversion_time.pdf` | Panel A: wall-clock solve time vs G |
| `output/supp_S2/s2_inversion_complexity.pdf` | Panel B: time and peak RSS vs G and sparsity |
| `output/supp_S2/s2_solver_cost.pdf` | Panel C: Marquardt, L-BFGS-B, Newton, gradient |
| `output/supp_S2/s2_kkt_gradient.pdf` | Panel D: KKT / simplex residual and score check |

### Covariance inversion: wall-clock time and complexity

Each likelihood, score, and Hessian evaluation factorises the mixture
covariance
\boldsymbol{\Sigma}(\boldsymbol{p})=\sum\_{j}p\_{j}^{2}\boldsymbol{\Sigma}\_{j}
once and then solves against one or more right-hand sides ([generative
model](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-numerical-speedups)).
Assembling \boldsymbol{\Sigma}(\boldsymbol{p}) costs O(J G^{2}). A dense
Cholesky factorisation of the G\times G result costs O(G^{3}) and
dominates when G is large. Structured backends (banded, block-diagonal,
Woodbury) replace that cubic with a cost that tracks bandwidth, block
size, or the low-rank residual r, which is why
[Figure 1](#fig-s2-inversion) must show **structure** and **backend**
together, not G alone.

**Panel A (running time).** Plot median wall-clock time per
[`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md)
against G on **log–log** axes, with one colour per backend and one
linetype per network structure. The log scale is the same device as
DICEPro Figure S1A: early-G cache effects stay visible while the
G=500–1000 regime does not collapse the faster solvers onto the axis. A
dense Cholesky should approach slope 3 on this plot. Banded and
block-diagonal backends should peel away below that line as soon as
fill-in stays cheaper than a full G^{3} factorisation. Ribbons are Monte
Carlo SD across the N\_{\text{timing}} repeated solves in Section 2, not
confidence intervals for a statistical estimand.

**Panel B (complexity and memory).** Following VIMixR Figure S4, treat
sparsity as a second experimental factor rather than a colour
annotation. Facet or group by edge density (or bandwidth / number of
blocks) and show both **seconds per solve** and **peak RSS**. The
intended reading is a three-way relationship: larger G inflates time
cubically for dense \boldsymbol{\Sigma}\_{j}; a sparser Markov network
reduces fill-in and peak memory; the structured backend only pays off
once that fill-in advantage exceeds the overhead of the specialised
kernel. At small G (G\le 50) dense Cholesky can still win because of
cache locality, as in the expected findings below.

``` r

# -----------------------------------------------------------------
# SECTION 3a: INVERSION PANELS (A, B)
# -----------------------------------------------------------------
benchmark_results <- readRDS("output/supp_S2/inversion_timing.rds")

p_time <- ggplot2::ggplot(
  benchmark_results,
  ggplot2::aes(G, mean_ms, colour = backend, linetype = structure)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = mean_ms - sd_ms, ymax = mean_ms + sd_ms, fill = backend),
    alpha = 0.15,
    colour = NA
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Number of genes G",
    y = "Solve time (ms, log10)"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  "output/supp_S2/s2_inversion_time.pdf",
  p_time,
  width = 8,
  height = 5
)
```

\(a\) Panel A: median wall-clock time per structured solve versus G on
log–log axes, coloured by backend and lined by network structure.

\(b\) Panel B: solve time and peak RSS versus G at several sparsity
levels (bandwidth, number of blocks, or edge density), the VIMixR-style
complexity surface.

Figure 1: Wall-clock cost of
\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\boldsymbol{v} as a function of
gene count G and Markov-network sparsity. Panel A compares inversion
backends; panel B relates the same times (and peak memory) to graph
sparsity. Assembly of \boldsymbol{\Sigma}(\boldsymbol{p}) is linear in J
and quadratic in G; dense Cholesky of the mixture covariance is cubic in
G.

### Constrained MLE solvers: cost, feasibility, and gradient checks

Covariance inversion is one callback. The MLE loop repeats that callback
until a convergence test fires.
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
exposes three second-order or quasi-Newton paths
([`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)):
Marquardt–Levenberg via `marqLevAlg` ([Philipps et al.
2023](#ref-R-marqLevAlg); [Commenges et al.
2006](#ref-commengesNewtonLikeAlgorithmLikelihood2006)), box-constrained
L-BFGS-B in \boldsymbol{p} via
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), and
Newton–Raphson / `nlminb` in ILR coordinates. The projected gradient in
ILR is the fourth comparator (the same quartet as the bivariate toy
vignette, excluding mean-only NNLS / `lsei` and simulated annealing).
Wall time for one sample therefore scales roughly as
(\text{iterations})\times\bigl(O(J G^{2})+O(G^{3})\bigr) for a dense
backend: **cubic in G**, and only **linear in J** through assembly of
\boldsymbol{\Sigma}(\boldsymbol{p}) and the J-vector score, not through
an extra matrix inverse per cell type.

**Panel C (solvers).** Three aligned metrics, in the spirit of DICEPro
Figure S1A–C:

- iterations (or function evaluations) until the package tolerance;
- total wall-clock time per sample;
- peak memory.

Do not map two of those scores to colour and size on one scatter. Facet
by metric, as in
[`plot_mc_metric_dots()`](https://bastienchassagnol.github.io/DeCovarT/reference/plot_mc_metric_dots.md).
Marquardt and Newton should reach the threshold in few iterations
because they use curvature; L-BFGS-B is cheaper per step but may need
more evaluations; projected gradient is the slow baseline (DICEPro’s
Adam / gradient-descent comparison). Report median and SD across Monte
Carlo samples, and keep J as a second x-axis or facet so the claimed
linear-J assembly cost is visible beside the cubic-G inversion cost.

**Panel D (simplex adherence and analytic gradients).** Two stacked
checks, combining DICEPro Figure S5 (unit-sum / reconstruction
trade-off) with their Figure S2A (analytical versus finite-difference
gradients).

The ILR / softmax chart is a C^{2} diffeomorphism onto the open simplex
([constrained
optimisation](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-constrained)),
so Marquardt, Newton, and the projected gradient keep \sum\_{j}\hat
p\_{j}=1 by construction. L-BFGS-B only boxes \[0,1\]^{J} and
renormalises after [`optim()`](https://rdrr.io/r/stats/optim.html);
during the line search \lvert 1-\sum\_{j}p\_{j}\rvert can leave zero,
which is why that path needs the guarded objective in [solver
safeguards](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-numerical-speedups).
The ADEMP field `optimisation$kkt_residual` ([synthetic
scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ademp))
is the projected-score residual
\lVert\Pi\_{\Delta}(\hat{\boldsymbol{p}}+\alpha\nabla\ell)-\hat{\boldsymbol{p}}\rVert\_{2}/\alpha.
It is the right analogue of DICEPro’s unit-sum violation for an interior
(or bound-constrained) stationary point: small residual means the score
has no feasible first-order improvement. Plot it on a log scale against
iteration or against reconstruction error, with a Pareto outline only if
a penalty-weight grid is actually run.

The companion scatter is analytic \nabla\ell (and, where stored, the
Hessian) against Richardson extrapolation in `numDeriv` ([Gilbert and
Varadhan 2019](#ref-R-numDeriv)), as already used in [the
generative-model
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.md)
and in `tests/testthat/test-03_03_DeCovarT.R`. DICEPro reported relative
errors below 6\times 10^{-6} with central differences
\varepsilon=10^{-5}. DeCovarT prefers Richardson’s tableau rather than a
single-step central difference, because the score mixes traces of
\boldsymbol{\Theta}(\boldsymbol{p}) with quadratic forms that have a
wide dynamic range. Points should fall on the identity; stratify by
parameter block (ILR coordinates versus unconstrained \boldsymbol{p} on
the L-BFGS-B path).

``` r

# -----------------------------------------------------------------
# SECTION 3b: SOLVER PANELS (C, D)
# -----------------------------------------------------------------
solver_results <- readRDS("output/supp_S2/solver_timing.rds")

p_solver <- ggplot2::ggplot(
  solver_results,
  ggplot2::aes(G, value, colour = solver)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_grid(metric ~ J, scales = "free_y", labeller = ggplot2::label_both) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Number of genes G",
    y = "Metric (log10)",
    colour = "Solver"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  "output/supp_S2/s2_solver_cost.pdf",
  p_solver,
  width = 9,
  height = 6
)
```

\(a\) Panel C: iterations to tolerance, wall-clock time, and peak memory
for Marquardt–Levenberg, L-BFGS-B, Newton–Raphson (`nlminb`), and
projected gradient, versus G and faceted by J.

\(b\) Panel D: log-scale KKT / simplex residual (feasibility of
\hat{\boldsymbol{p}}) and analytic score versus Richardson `numDeriv`
gradient, stratified by parameter block.

Figure 2: Computational profile of the four DeCovarT MLE solvers and
first-order diagnostics. Panel C separates speed of convergence from
wall-clock time and memory. Panel D checks unit-simplex / KKT adherence
and confirms that analytic derivatives match Richardson extrapolation.

### Expected findings

Dense Cholesky should dominate for G \le 50 where cache effects matter
more than sparsity. For G \ge 250 with banded or block-diagonal
structure, the corresponding sparse solvers should be 5–20\times faster.
Woodbury is competitive for low-rank perturbations (rank r = 5) at all
G. On log–log axes the dense inversion curve should approach slope 3 in
G; increasing J should shift times almost linearly through assembly of
\boldsymbol{\Sigma}(\boldsymbol{p}), not cubically.

Marquardt–Levenberg and Newton–Raphson should need far fewer iterations
than projected gradient, at a higher cost per step. L-BFGS-B should sit
between them in wall-clock time, with a larger KKT / unit-sum residual
before the final renormalisation. Analytic and Richardson gradients
should agree to relative error \ll 10^{-6} on interior points;
discrepancies flag a backend or ILR-chart bug, not optimiser tuning.

## Out-of-scope extensions

The following factors from the full simulation plan
(`temp_rag_simualtions_scenarios/detailled_siulation_code.md` §Full
factorial) are noted here but not yet implemented:

- **Low-rank residual r**: planned levels 0, 2, 8, 20 — would sharpen
  the Woodbury cross-over point.
- **Cholesky fill ratio** vs graph density: a finer grid would reveal
  when sparse-Cholesky reordering (AMD/COLAMD) pays off.
- **Condition number effect on solve accuracy**:
  \kappa(\boldsymbol{\Sigma}\_j) \in \\3, 10, 100, 1000\\ with iterative
  refinement.

## See also

- [Appendix S3 —
  Scaling](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S3-scaling.html)
  (RMSE and difficulty axes; this appendix is inversion and solver cost)
- [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-ademp)
- Structure-aware operators and solver safeguards: [generative
  model](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-numerical-speedups)

## References

Commenges, Daniel, Helene Jacqmin-Gadda, Cecile Proust, and Jeremie
Guedj. 2006. *A Newton-Like Algorithm for Likelihood Maximization: The
Robust-Variance Scoring Algorithm*. arXiv.
<https://doi.org/10.48550/arxiv.math/0610402>.

Gilbert, Paul, and Ravi Varadhan. 2019. *numDeriv: Accurate Numerical
Derivatives*. <http://optimizer.r-forge.r-project.org/>.

Philipps, Viviane, Cecile Proust-Lima, Melanie Prague, Boris Hejblum,
Daniel Commenges, and Amadou Diakite. 2023. *marqLevAlg: A Parallelized
General-Purpose Optimization Based on Marquardt-Levenberg Algorithm*.
