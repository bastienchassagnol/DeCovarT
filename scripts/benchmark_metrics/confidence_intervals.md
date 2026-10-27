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

## Checks (scripts/auxiliary/check_confidence_interval_coverage.R)

Restricted fits on a face now use the trust-region Newton driver with the
pooled analytic Hessian (`.fit_on_face()`); the former BFGS start could
jump to a vertex with pooled columns. Population intervals cover at
0.93 to 0.96 (J = 2, 3; N = 5; 300 replicates). Single-cell intervals
undercover at G = 10 (0.85 to 0.92), because the variance scale has
G - J degrees of freedom; at G = 40 the Wald interval covers at 0.94.
`variance = "reml"` in `confint_single_cell()` and
`confint_profile_single_cell()` inflates the covariances by G / (G - J),
uses t quantiles on G - J degrees of freedom (a multivariate t for the
ILR map when J > 2), the chi-square interval for phi, and the profile
cut-off G log(1 + F(1, G - J) / (G - J)).

## Not in this pass

Score-test inversion, Bartlett correction, \(r^{\ast}\), the Godambe
sandwich, and any extra term for using \(\hat\Sigma_j\) in place of
\(\Sigma_j\). The reference bootstrap is the implemented route for
reference uncertainty. Coverage of the new intervals has not been
checked by a Monte Carlo experiment.
