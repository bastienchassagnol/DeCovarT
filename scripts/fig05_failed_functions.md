# GitHub DeCovarT: functions exercised by the figure scripts

Install: `pak::pkg_install("bastienchassagnol/DeCovarT", dependencies = TRUE, upgrade = TRUE)`.

| Field | Value |
| --- | --- |
| Package | DeCovarT 2.3.1 |
| Remote | `github::bastienchassagnol/DeCovarT` |
| Commit | `164d3f81ecab10ed2b0a0860bb930a475f4340eb` |
| Library | `/home/bastienc/R/x86_64-pc-linux-gnu-library/4.6` |

The three figure scripts were run against this install (`library(DeCovarT)`), with local helpers from `scripts/benchmark_metrics/` and `scripts/visualisation/`. `N = 200` bulks and `B = 200` bootstrap refits.

## Missing names

None. Every function those scripts call is in the installed namespace. Six are not exported and were taken from that namespace with `get()`: `.hellinger_gaussian`, `.mixture_factorisation`, `.ui_h1`, `.ui_info`, `.ui_success`, `scale_covariances_to_overlap`.

## Runtime failure

`scripts/fig05_single_cell_hierarchical.R` stops while building the `G = 200` moments (`lineage_moments(200, seed = 20261024 + 200)`).

| Item | Detail |
| --- | --- |
| Exported caller | `build_covariance_array_from_precision()` |
| Internal callee | `.invert_precision()` in `R/01_02_generate_synthetic_networks.R` |
| Message | Could not invert a precision to a strictly positive-definite covariance by a uniform spectral shift. |
| `G = 20` | Same builder, seed `20261024 + 100`, succeeds (covariance array `20 x 20 x 4`). |
| `G = 200` | Fails. Scale-free sibling graph, `edges_per_node = 1`, `precision_shift = 0.2`, 5% sign flips. |

Figure 02 and figure 03 passed this step and entered the benchmark (`N = 200`, `B = 200`).
