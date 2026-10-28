# How the figure benchmarks simulate samples

Three sampling regimes are involved. Notation: \(J\) cell types, \(G\)
genes, \(N\) bulks per scenario, \(B\) bootstrap refits per bulk, \(R\)
individuals in a reference panel.

## 1. No bootstrap (figures 2, 3 and 5; \(N = 200\), \(B = 1\))

Each scenario fixes one set of moments \((\boldsymbol\mu,
\boldsymbol\Sigma_1,\dots,\boldsymbol\Sigma_J)\) and draws \(N\)
independent bulks from it (`R/04_02_run_simulation_benchmark.R`). Only the
bulks vary, so \(B = 1\) means one draw of the moments, not a resampling.

- **Population model** (figures 2 and 3, `simulate_bulk_mixture()` in
  `R/01_01_population_toy_simulation.R`): for bulk \(i\) and type \(j\),
  \(\boldsymbol x_{ij}\sim\mathcal N_G(\boldsymbol\mu_j,
  \boldsymbol\Sigma_j)\), and \(\boldsymbol y_i=\sum_j p_j\boldsymbol
  x_{ij}\), so \(\boldsymbol y_i\sim\mathcal N(\boldsymbol\mu\boldsymbol p,
  \sum_j p_j^2\boldsymbol\Sigma_j)\). The target is \(\boldsymbol p\).
- **Single-cell model** (figure 5, `simulate_bulk_from_cells()` in
  `R/01_01_single_cell_toy_simulation.R`): cell counts \(\boldsymbol m_i\sim
  \mathrm{Multinomial}(n_i,\boldsymbol p)\); one Gaussian draw per type for
  the sum of its cells, \(\boldsymbol T_{ji}\sim\mathcal N(m_{ji}
  \boldsymbol\mu_j, m_{ji}\boldsymbol\Sigma_j)\), where \(\boldsymbol\mu_j\)
  carries the cell size \(S_j\) and \(\boldsymbol\Sigma_j\) carries
  \(S_j^2\); then \(\boldsymbol y_i=a\sum_j\boldsymbol T_{ji}\). \(n_i\) and
  the platform scale \(a\) are fixed per scenario. The target is the
  realised fraction \(\boldsymbol m_i/n_i\); the nuisance scales are
  \(c=a\,n_i\) and \(\phi=a^2n_i\).

DeCovarT's solvers receive the true \((\boldsymbol\mu,\boldsymbol\Sigma)\):
in the terms used by reference-based methods, a single reference whose
moments are known exactly, with no reference sampling error. Each bulk
gives one estimate and one interval per construction (Wald, ILR, profile,
Firth, oracle; REML variants in figure 5), and coverage is the share of
the \(N\) intervals that contain the target.

## 2. Parametric bootstrap (package API; no longer in the benchmarks)

For one observed bulk and its estimate \(\hat{\boldsymbol p}\), \(B\) new
bulks are drawn from the fitted model with the same \(\boldsymbol\mu\) and
\(\boldsymbol\Sigma\), each is refitted, and the interval is the percentile
interval of the \(B\) refits.

- Population (`bootstrap_decovart()`, `R/03_03_resampling_inference.R`):
  \(\boldsymbol y^{(b)}\sim\mathcal N(\boldsymbol\mu\hat{\boldsymbol p},
  \sum_j\hat p_j^2\boldsymbol\Sigma_j)\).
- Single cell (`bootstrap_single_cell()`, `R/03_04_single_cell_inference.R`):
  \(\boldsymbol y^{(b)}\sim\mathcal N(\hat c\,\boldsymbol\mu\hat{\boldsymbol
  p}, \hat\phi\sum_j\hat p_j\boldsymbol\Sigma_j)\), a Gaussian draw without
  the multinomial step.

A benchmark scenario then costs \(N\times B\) refits per solver (40,000 at
\(N=B=200\)), which dominated the 2026-10-27 runs. The reference bootstrap
(`reference_bootstrap_decovart()`, `reference_bootstrap_single_cell()`)
resamples donors of a reference panel instead; it was never part of the
benchmarks.

## 3. Reference panels for second-generation methods (\(R = 200\))

MEAD and RNA-Sieve estimate the reference moments themselves, so figures 3
and 5 also draw, for each scenario, an independent panel
(`simulate_reference_replicates()`, `R/01_06_reference_replicates.R`, seed
`reference$seed + i`). For each type \(j\) it holds \(R\times n_{\text{cells}}\)
independent columns from \(\mathcal N_G(\boldsymbol\mu_j,\boldsymbol
\Sigma_j)\), labelled by individual \(1,\dots,R\):

- figure 3: \(n_{\text{cells}} = 1\), i.e. \(R\) purified profiles per type;
- figure 5: \(n_{\text{cells}} = 20\) cells per individual, each a single
  cell at the per-cell moments.

All individuals share the same \(\boldsymbol\mu_j\) and
\(\boldsymbol\Sigma_j\): there is no between-individual effect, so the
labels are exchangeable. How the methods use the panel
(`scripts/benchmark_metrics/second_generation_methods.R`):

- **MEAD** receives every cell with its type and individual, estimates the
  means and the between-individual covariance after scaling each cell to
  a unit total, and so estimates transcript fractions (scored against
  them as `native (transcript fractions)`).
- **RNA-Sieve** receives one file per type and pools its cells; it has no
  individual labels.
- **ReDeconv** does not use the panel: it gets the true means and
  \(\sqrt{\operatorname{diag}\boldsymbol\Sigma_j}\).

DeCovarT never sees the panel; \(R\) matters only to MEAD and RNA-Sieve.
