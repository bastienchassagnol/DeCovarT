# Appendix S4 — Covariance Modelling Assumptions

``` r

library(DeCovarT)
```

> **Script:** `scripts/supp_S4_covariance_modeling.R` — compares
> covariance modelling assumptions. Outputs in `output/supp_S4/`.

------------------------------------------------------------------------

## :one: Generative Model

### Motivation

The DeCovarT likelihood requires specifying \boldsymbol{\Sigma}\_j — the
within-cell-type covariance. In practice, these are estimated from
external (reference) data and may be mis-specified. This appendix
compares four **inference models** that differ in how
\boldsymbol{\Sigma}\_j is supplied to the solver:

| Model | Description | Oracle? |
|----|----|----|
| **True full** | Exact \boldsymbol{\Sigma}\_j from the data-generating process | Yes |
| **Cell-type diagonal** | \operatorname{diag}(\boldsymbol{\Sigma}\_j) — ignores off-diagonal | Partial |
| **Global weighted diagonal** | \boldsymbol{D}\_\text{global} = \operatorname{diag}\bigl(\sum_j \frac{1}{J^2}\boldsymbol{\Sigma}\_j\bigr), replicated across types | Partial |
| **Identity** | \boldsymbol{I}\_G — mean-only equivalent | No |

The global weighted diagonal uses balanced weights p_j = 1/J and is
non-oracular, matching the standard practice when only pooled
variability is known.

### Factorial design

| Factor | Levels |
|----|----|
| Inference model | True full, cell-type diagonal, global weighted diagonal, identity |
| Mean cosine \rho\_\mu | 0, 0.5, 0.9 |
| **Total** | 4 \times 3 = \mathbf{12} scenarios |

G = 50 genes, J = 3 cell types; N = 200 Monte Carlo replicates (smoke
test: n = 2).

``` r

# -----------------------------------------------------------------
# SECTION 1: GENERATIVE MODEL
# Build true Sigma matrices and compute mis-specified alternatives
# -----------------------------------------------------------------
SEED <- 20260807L
G    <- 50L;  J <- 3L

cov_models_grid <- tidyr::expand_grid(
  cosine_mu   = c(0.0, 0.5, 0.9),
  cov_model   = c("true_full", "ct_diagonal", "global_diagonal", "identity")
) |>
  dplyr::mutate(scenario_id = dplyr::row_number())

build_true_cov <- function(G, J, seed = SEED) {
  purrr::map(seq_len(J), \(j) {
    net <- generate_random_network_skeleton(G, model = "erdos_renyi",
                                           seed = seed + j)
    build_normalised_precision(net) |> solve()
  })
}

true_covs <- build_true_cov(G, J, SEED)

# Mis-specified alternatives
ct_diag_covs <- purrr::map(true_covs, \(Sigma) diag(diag(Sigma)))

# Global weighted diagonal D_global = diag( sum_j (1/J^2) Sigma_j )
D_global_vec <- Reduce(`+`, purrr::map(true_covs, \(S) S / J^2)) |>
  diag()
global_diag_covs <- purrr::map(seq_len(J), \(.) diag(D_global_vec))

identity_covs <- purrr::map(seq_len(J), \(.) diag(G))

cov_specs <- list(
  true_full       = true_covs,
  ct_diagonal     = ct_diag_covs,
  global_diagonal = global_diag_covs,
  identity        = identity_covs
)
```

------------------------------------------------------------------------

## :two: Inference

``` r

# -----------------------------------------------------------------
# SECTION 2: INFERENCE
# -----------------------------------------------------------------
N_REPLICATES <- as.integer(Sys.getenv("N_REPLICATES", unset = "200"))
ALGORITHMS   <- c("Marquardt-Levenberg")   # fix solver, vary Sigma input
dir.create("output/supp_S4", recursive = TRUE, showWarnings = FALSE)

s4_results <- purrr::pmap(
  list(cov_models_grid$cosine_mu, cov_models_grid$cov_model),
  \(rho, model) {
    mu <- generate_mean_signature_matrix(J, G, cosine_similarity = rho,
                                         mean_scale = 100, seed = SEED)
    run_simulation_benchmark(
      n_replicates          = N_REPLICATES,
      true_proportions      = rep(1 / J, J),
      mean_signature_matrix = mu,
      covariance_matrices   = cov_specs[[model]],
      algorithm_names       = ALGORITHMS
    )
  },
  .progress = TRUE
)

saveRDS(s4_results, "output/supp_S4/cov_model_benchmark.rds")
saveRDS(cov_models_grid, "output/supp_S4/cov_models_grid.rds")
```

------------------------------------------------------------------------

## :three: Visualisations

| Output | Description |
|----|----|
| `output/supp_S4/s4_rmse_by_model.pdf` | RMSE by covariance model and mean cosine |
| `output/supp_S4/s4_bias_by_model.pdf` | Bias comparison across models |
| `output/supp_S4/s4_delta_vs_truefull.pdf` | RMSE gap vs. true-full oracle |

``` r

# -----------------------------------------------------------------
# SECTION 3: VISUALISATIONS
# -----------------------------------------------------------------
s4_results      <- readRDS("output/supp_S4/cov_model_benchmark.rds")
cov_models_grid <- readRDS("output/supp_S4/cov_models_grid.rds")

if (N_REPLICATES >= 10) {
  summary_df <- cov_models_grid |>
    dplyr::mutate(
      rmse = purrr::map_dbl(s4_results, \(res) {
        mean(res[["Marquardt-Levenberg"]][["RMSE"]], na.rm = TRUE)
      })
    )

  p_rmse <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(factor(cosine_mu), rmse, fill = cov_model)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(x = "Mean cosine ρ_μ", y = "Mean RMSE",
                  fill = "Covariance model") +
    ggplot2::theme_minimal()

  ggplot2::ggsave("output/supp_S4/s4_rmse_by_model.pdf", p_rmse, width = 8, height = 5)
}
```

### Expected findings

- **True full** should have the lowest RMSE at all cosine levels.
- **Cell-type diagonal** degrades gracefully: off-diagonal elements add
  variance but the diagonal captures heteroscedasticity.
- **Global weighted diagonal** loses cell-type specificity; its RMSE
  tracks the identity model when cell types are heteroscedastic.
- **Identity** (mean-only) should perform comparably to NNLS, confirming
  that ignoring second-order structure removes the DeCovarT advantage.
- The **RMSE gap** (true full minus alternatives) should increase with
  cosine \rho\_\mu, because near-collinear means force greater reliance
  on the covariance for discrimination.

## Out-of-scope: reference noise

An important practical extension (noted in the full factorial plan) is
testing **reference covariance estimation noise** — supplying
\hat{\boldsymbol{\Sigma}}\_j estimated from a finite sample of
n\_\text{ref} \in \\20, 50, 200\\ single-cell donors instead of the
oracle. This adds a shrinkage/regularisation dimension to the modelling
comparison and is planned for a future update.

## See also

- [Appendix S1 —
  Identifiability](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.html)
- [Appendix S3 —
  Scaling](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S3-scaling.html)
- [Appendix S5 — Model
  Mis-specification](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S5-misspecification.html)

## References
