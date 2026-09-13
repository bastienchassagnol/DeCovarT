# 16-page 2D latent-space book for the covariance-driven scenario

Extreme MixSim overlap (low / high) crossed with extreme Shannon
composition (balanced / highly unbalanced) and the four CT1/CT2 graph
assignments (\\2 \times 2 \times 4 = 16\\ pages). Each page is a 3-by-2
of projections: independent Thomson factor analyses (Thomson 1938) on
the purified types (top), then a shared plane from mixtures of common
factor analysers (McLachlan and Peel 2000) on convolution draws, the
same MCFA on an unsupervised mixture with weights \\\boldsymbol{p}\\,
and supervised
[`mclust::MclustDR()`](https://mclust-org.github.io/mclust/reference/MclustDR.html)
(Scrucca 2010, 2015) (bottom).

## Usage

``` r
save_hybrid_latent_projection_book(
  artefacts,
  file,
  n = 300L,
  seed = NULL,
  itmax = 150L,
  data_rds = NULL
)
```

## Arguments

- artefacts:

  Split or assembled fig03 artefacts from
  [`read_simulation_artefacts()`](https://bastienchassagnol.github.io/DeCovarT/reference/read_simulation_artefacts.md).

- file:

  Output PDF path.

- n:

  Draws per cell type / bulk sample (default 300).

- seed:

  Optional seed; `NULL` leaves the RNG state unchanged.

- itmax:

  EM iterations for
  [`EMMIXmfa::mcfa()`](https://rdrr.io/pkg/EMMIXmfa/man/mcfa.html).

- data_rds:

  Optional directory for ggplot `data` RDS files.

## Value

`file`, invisibly.

## References

McLachlan GJ, Peel D (2000). “Mixtures of Factor Analyzers.” In *Finite
Mixture Models*, 238–256. John Wiley & Sons, Ltd.
[doi:10.1002/0471721182.ch8](https://doi.org/10.1002/0471721182.ch8) .  
  
Scrucca L (2010). “Dimension Reduction for Model-Based Clustering.”
*Statistics and Computing*, **20**(4), 471–484.
[doi:10.1007/s11222-009-9138-7](https://doi.org/10.1007/s11222-009-9138-7)
.  
  
Scrucca L (2015). “Graphical Tools for Model-Based Mixture Discriminant
Analysis.” https://arxiv.org/abs/1508.01695v1.
[doi:10.1007/s11634-013-0147-1](https://doi.org/10.1007/s11634-013-0147-1)
.  
  
Thomson GH (1938). “Methods of Estimating Mental Factors.” *Nature*,
**141**(3562), 246–246.
[doi:10.1038/141246a0](https://doi.org/10.1038/141246a0) .

## See also

[`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md),
[`save_hybrid_runtime_book()`](https://bastienchassagnol.github.io/DeCovarT/reference/save_hybrid_runtime_book.md)

## Examples

``` r
skip <- !requireNamespace("EMMIXmfa", quietly = TRUE) ||
  !requireNamespace("mclust", quietly = TRUE) ||
  !requireNamespace("cowplot", quietly = TRUE)
if (!skip) {
  genes <- paste0("g", seq_len(6L))
  cts <- paste0("celltype_", 1:3)
  mu <- matrix(
    c(8, 9, 9, 8, 2, 12, 2, 3, 11, 10, 4, 5, 5, 4, 12, 3, 6, 7),
    nrow = 6L,
    dimnames = list(genes, cts)
  )
  sig <- diag(6)
  sigma <- array(c(sig, sig, sig), dim = c(6L, 6L, 3L))
  dimnames(sigma) <- list(genes, genes, cts)
  th <- list(p = c(0.5, 0.3, 0.2), mu = mu, sigma = sigma)
  artefacts <- list(
    config = tibble::tibble(
      ID = "V1",
      proportions = "balanced",
      overlap_label = "low",
      graph_ct1 = "scale_free",
      graph_ct2 = "scale_free",
      graph_ct3 = "scale_free"
    ),
    theta = tibble::tibble(ID = "V1", true_theta = list(th))
  )
  tf <- withr::local_tempfile(fileext = ".pdf")
  save_hybrid_latent_projection_book(
    artefacts,
    tf,
    n = 40L,
    seed = 1L,
    itmax = 20L
  )
}
```
