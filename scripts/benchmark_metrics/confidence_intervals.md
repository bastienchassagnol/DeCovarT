# Confidence intervals

Two models, one menu. `confint_decovart()` is the population convolution
(covariance weights \(p_j^2\), no technical scales). `confint_single_cell()`
is the single-cell mixture (covariance weights \(p_j\), scales \(c\) and
\(\phi\) profiled out). Both default to `method = "ilr"`.

`confint.decovart_fit()` is unchanged. It is the per-column delta method
for a `decovart_fit`, and it can leave \([0,1]\). `confint_decovart()`
pools columns that share one composition and multiplies the Fisher
information by the number of columns.

The single-cell fit now returns `vcov_rho`, `vcov` and `vcov_nuisance`
from `deconvolute_single_cell()`. `vcov` is
\(\mathbf{J}_{\psi}(-\mathbf{H}_{\rho})^{-1}\mathbf{J}_{\psi}^{\top}\),
with \(\mathbf{H}_{\rho}\) the observed profile Hessian.

Higher-order corrections and a sandwich for plug-in covariances are not
implemented.

## Population convolution

| Request | Function | Default role |
|---|---|---|
| Wald, delta method on \(p\) | `confint_decovart(method = "wald")`, `vcov_ilr_delta()`, `confint.decovart_fit()` | Symmetric interval. Endpoints can leave \([0,1]\). Undefined on a face. |
| Firth-penalised Wald | `fit_firth_decovart()`, `confint_decovart(method = "firth")` | Jeffreys penalty \(\tfrac12\log\det I_{\rho}\) on the expected information. The reported interval is Wald at the penalised estimate. `interval_ilr` maps that estimate through the inverse chart. The penalty gradient is numerical; the likelihood score is analytic. Not an exact \(O(N^{-1})\) bias correction. |
| ILR inverse chart | `confint_ilr_backtransform()`, `confint_decovart()` | **Default.** For two cell types, the Wald endpoints of \(\rho\) are mapped by \(\psi\). For three or more, \(\rho^{(b)}\sim N(\hat\rho, I_{\rho}^{-1})\) and the interval is the marginal quantile of \(\psi(\rho^{(b)})\). `rho_interval` always holds the Wald bounds of \(\rho\). |
| Interior likelihood-ratio interval | `confint_profile_decovart()`, `confint_decovart(method = "profile")`, `lrt_decovart()` | Wilks \(\chi^2_1\). The test of an interior null is `lrt_decovart()`. |
| Boundary likelihood-ratio test | `lrt_decovart()` | Automatic chi-bar-square when a null value is on a face. This is a test, not a two-sided interval. A profile endpoint on a face should be read with this calibration. |
| Parametric bootstrap | `bootstrap_decovart()`, `confint_decovart(method = "bootstrap")` | Draws bulk profiles from the fitted convolution and refits. `null_value` adds a restricted Monte Carlo p-value. |
| Reference bootstrap | `reference_bootstrap_decovart()`, `confint_decovart(method = "reference")` | See the resampling rule below. |

## Single-cell mixture

The composition interval uses the profile information, the Schur
complement of the observed Hessian that removes \((c,\phi)\).

The scale interval uses the other Schur complement, which removes
\(\rho\). `nuisance$interval` is a Wald interval on \((\log c,\log\phi)\),
mapped back by the exponential, so it stays positive.
`nuisance$interval_wald` is the symmetric normal interval on \((c,\phi)\)
itself. Bootstrap methods replace both with percentiles.
\(I_{c\phi}=0\) in expectation does not make this block a pair of
one-parameter problems: the complement still accounts for estimation of
the composition.

| Request | Function |
|---|---|
| ILR inverse chart (default) | `confint_single_cell()` |
| Delta method | `confint_single_cell(method = "wald")`, `vcov_single_cell()` |
| Firth on the profile | `fit_firth_single_cell()`, `confint_single_cell(method = "firth")` |
| Interior profile interval | `confint_profile_single_cell()`, `confint_single_cell(method = "profile")` |
| Interior or boundary test | `lrt_single_cell()` |
| Parametric bootstrap | `bootstrap_single_cell()`, `confint_single_cell(method = "bootstrap")` |
| Reference bootstrap | `reference_bootstrap_single_cell()`, `confint_single_cell(method = "reference")` |
| Face-constrained fit | `restricted_mle_single_cell()` |

`estimate` is the composition that the likelihood optimises
(`composition` on the fit), before a cellular / transcript relabeling.

## Reference resampling

`confint_decovart(method = "reference")` and
`confint_single_cell(method = "reference")` choose as follows.

- Every cell type has at least two donor labels: resample donors within
  each type, rebuild \((\mu,\Sigma_j)\), and keep the observed bulk.
- Otherwise: resample columns within each cell type and regenerate the
  bulk from the fitted model (original moments), refitted with the
  resampled moments. (Drawing it from the resampled moments, as first
  implemented, hid the reference error: coverage 0.76 instead of 0.92.)

`resample` and `regenerate_bulk` override that rule.
`reference_bootstrap_decovart()` itself still defaults to
`method = "donors"` and `regenerate_bulk = FALSE`.

Shuffling cell-type names is not an interval. The fitted composition
moves with the columns, which is equivariance rather than sampling
variability. For the population model that check is
`equivariance_check_decovart()`.

## Checks (scripts/benchmark_metrics/check_confidence_interval_coverage.R)

Restricted fits on a face now use the trust-region Newton driver with the
pooled analytic Hessian (`.fit_on_face()`); the former BFGS start could
jump to a vertex with pooled columns. Run of 2026-10-24 (300 replicates,
level 0.95, seed 20261022), after the cell size entered the covariance
as S_j^2:

- Information matrices: the single-cell profile and nuisance
  covariances match the inverse numDeriv Hessian of l(rho, log c, log phi)
  to a relative error below 5e-7 (J = 2, 3, 4); the population expected
  information matches the Monte Carlo mean of -H to 1.5e-2 (J = 2) and
  5.7e-3 (J = 3), the Monte Carlo error of 4000 draws.
- Population intervals (J = 2, 3; N = 5): ILR 0.95 to 0.99, Wald 0.95
  to 0.97, profile 0.95 to 0.97, Firth 0.95 to 0.97. Parametric
  bootstrap 0.96, reference bootstrap 0.92 (J = 2, 50 replicates).
- Single-cell intervals (G = 10, J = 2, four scenarios of balanced or
  unbalanced composition and equal or unequal cell sizes): with
  `variance = "ml"` the composition intervals cover 0.84 to 0.93 and
  those of c and phi 0.82 to 0.90, because the variance scale has G - J
  degrees of freedom. `variance = "reml"` in `confint_single_cell()` and
  `confint_profile_single_cell()` inflates the covariances by
  G / (G - J), uses t quantiles on G - J degrees of freedom (a
  multivariate t for the ILR map when J > 2), the chi-square interval for
  phi, and the profile cut-off G log(1 + F(1, G - J) / (G - J)); it
  covers 0.93 to 0.97 for the composition and 0.94 to 0.96 for the
  scales.
- Boundary likelihood-ratio tests at the 0.05 level: rejection rates
  0.067, 0.060 and 0.043 for the population nulls (Monte Carlo SE about
  0.014), 0.093 for the single-cell null celltype_2 = 0 (SE 0.017),
  which is anti-conservative.

## Not in this pass

Score-test inversion, Bartlett correction, \(r^{\ast}\), the Godambe
sandwich, and any extra term for using \(\hat\Sigma_j\) in place of
\(\Sigma_j\). The reference bootstrap is the implemented route for
reference uncertainty. Coverage under the simulation benchmarks of
figures 2, 3 and 5 is in their `intervals` tables
(`output/figXX/figXX_benchmark.rds`).

The figure benchmarks run neither bootstrap: their menus are Wald, ILR,
profile, Firth and oracle (plus the REML variants in figure 5). The
parametric bootstrap cost \(N\times B\) refits per solver and scenario
(40,000 at \(N=B=200\)). The bootstrap functions remain in the package and
are covered by the checks above. See `sampling_design.md` for how the
bulks and reference panels are drawn.
