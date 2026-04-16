# Links to the paper

- [Overleaf repository, shared with Anais Baudot](https://www.overleaf.com/project/697398ed149e81845cb5208e)
- [Hal Archive](https://arxiv.org/pdf/2309.09557)

# Project structure

## Simulations of synthetic networks

All six covariance matrices pass the positive-definiteness check for both graph models. Here is a summary of what was created.

---

**File:** `R/02_02_generate_synthetic_networks.R`

| Function | Role |
|---|---|
| `generate_mean_profiles` | Draws complementary parent means from U(lower, upper) for expressed/background blocks; perturbs via \(\delta \sim \mathcal{N}(0, \sigma_\delta^2 I)\) to produce 4 child means |
| `compute_marginal_variances` | Applies NB-like law \(\sigma_g^2 = \mu_g + \mu_g^\alpha / L\) |
| `generate_random_network_skeleton` | Wraps `igraph::sample_pa` (power-law) or `igraph::sample_sbm` (SBM) |
| `build_normalised_precision` | Affine transform \(\Omega = A \odot v + (\lvert\lambda_{\min}\rvert + u)I\) ensuring \(\Omega \succ 0\) |
| `build_precision_matrix` | Scales shared correlation \(R = \text{cov2cor}(\Omega^{-1})\) by population-specific variances into a 3D array |
| **`simulate_hierarchical_grn_moments`** | Orchestrates the full pipeline |

**Return value** — named list with three elements:

- `parent_parameters` — `mean_profiles` (2 x 2n matrix) + `covariance_matrices` (2n x 2n x 2 array)
- `child_parameters` — `mean_profiles` (4 x 2n matrix) + `covariance_matrices` (2n x 2n x 4 array)
- `graph_structure` — `adjacency_matrix` + `normalised_precision`

All row/column/slice names are explicit (`gene_1`, `parent_1`, `parent_1_child_a`, etc.). Full roxygen2 documentation with `@export`, `@importFrom`, `@examples`, and `@param` / `@return` blocks.