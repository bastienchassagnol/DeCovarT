# Appendix S5 — Model Mis-specification (Description Only)

> **Note**
>
> This vignette describes **planned but not yet implemented**
> mis-specification experiments. Implementation is deferred pending
> realistic network generators and validated noise models for bulk
> RNA-seq; see §Blocking issues below.

------------------------------------------------------------------------

## Overview

The DeCovarT model assumes that bulk expression arises from a Gaussian
mixture:

\boldsymbol{y} \sim \sum\_{j=1}^{J} p_j \\
\mathcal{N}(\boldsymbol{\mu}\_j,\\ \boldsymbol{\Sigma}\_j).

Real bulk RNA-seq violates this assumption in at least five ways:

1.  **Observation noise model**: counts follow a negative binomial (NB)
    or Poisson-log-normal (PLN) distribution rather than Gaussian.
2.  **Heavy-tailed or contaminated bulk**: outlier samples (technical
    artefacts, cell doublets) inflate tails.
3.  **Missing reference component**: an unmodelled cell type contributes
    to the bulk; all estimated proportions are redistributed.
4.  **Inter-sample (copula) correlations**: genes share latent batch
    structure, inducing correlation across observations.
5.  **Wrong covariance reference**: the within-cell-type covariance
    \boldsymbol{\Sigma}\_j is estimated from a different cohort or
    technology, introducing systematic bias.

------------------------------------------------------------------------

## Planned experiments

### S5a — Observation noise model

| Factor | Levels |
|----|----|
| Data-generating law | Gaussian (baseline); t\_\nu (\nu = 3); log-normal/PLN; contaminated Gaussian (\varepsilon = 0.05) |
| Algorithm | NNLS, DeconRNASeq, Marquardt–Levenberg (full \boldsymbol{\Sigma}) |
| G, J | 50 genes, 3 cell types |

**Hypothesis:** DeCovarT should retain lower bias under log-normal noise
than NNLS because the covariance captures more of the signal geometry,
but may be more sensitive than NNLS to heavy-tail contamination.

Connection to HADACA3 ([Barbot and Richard
2026](#ref-barbotPromisesLimitsMultimodal2026)): the `SDC5` dataset uses
a copula + NB noise model; `SDE5` uses an EM-based noise model. These
represent the practical domain of mis-specification most relevant to
pancreatic tumour deconvolution.

### S5b — Missing reference component

| Factor                | Levels                                          |
|-----------------------|-------------------------------------------------|
| Missing-type fraction | None (0%); rare (2%); moderate (10%)            |
| Detection strategy    | None (ignore); re-normalise; add null cell type |
| Metric                | RMSE of remaining types; false inclusion rate   |

**Hypothesis:** closure forces the missing mass onto the J modelled cell
types; the distribution of the error across types depends on the
mean-collinearity structure. DeCovarT may exhibit a different
redistribution pattern than NNLS because the second-order term couples
the proportions estimation.

### S5c — Inter-sample copula correlations

| Factor               | Levels                                       |
|----------------------|----------------------------------------------|
| Copula family        | Gaussian; Clayton (positive tail dependence) |
| Marginal             | NB (r = 5, \mu = \bar{\mu}\_g)               |
| Correlation strength | \tau_K \in \\0, 0.2, 0.5\\                   |

This directly mirrors the `SDC5` HADACA3 scenario ([Barbot and Richard
2026](#ref-barbotPromisesLimitsMultimodal2026)).

------------------------------------------------------------------------

## Blocking issues (why not implemented yet)

1.  **Gaussian-only
    [`simulate_bulk_mixture()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_bulk_mixture.md)**:
    the current package implementation generates bulk expression as a
    Gaussian mixture. NB/PLN data generation would require a dedicated
    simulator (e.g. integration with `splatter` ([Zappia et al.
    2017](#ref-zappiaSplatterSimulationSinglecell2017a)) or `SymSim`).

2.  **Copula sampling**: copula-correlated NB samples require either a
    dedicated `VineCopula`-based generator or the HADACA3 framework code
    (currently not a DeCovarT dependency).

3.  **No count-data likelihood**: DeCovarT fits a Gaussian
    log-likelihood. Comparing it under NB data requires either a
    properly mis-specified analysis (fit Gaussian to count data) or a
    new NB likelihood branch.

These blockers are tracked as GitHub issues; this vignette will be
updated once a count-data simulation pathway is merged.

------------------------------------------------------------------------

## Connection to the full factorial plan

The full factorial plan (`detailled_siulation_code.md`) includes these
additional factors that touch on mis-specification and are also out of
scope:

| Factor | Planned levels | Notes |
|----|----|----|
| Data-generating law | Gaussian; t\_\nu; PLN; contaminated Gaussian | S5a above |
| Missing-reference component | None; 2%; 10% | S5b above |
| Reference covariance sample size | Oracle; 20; 50; 200 donors | See Appendix S4 |
| Initial composition (multistart) | Equi-balanced; Dirichlet(1); 10 multistarts | See Appendix S1 |

------------------------------------------------------------------------

## See also

- [Appendix S1 —
  Identifiability](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S1-identifiability.html)
- [Appendix S4 — Covariance
  Modelling](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S4-covariance-modeling.html)
- [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/how-to-build-synthetic-scenarios-mean-covariance.md)

## References

Barbot, Hugo, and Magali Richard. 2026. ‘On the Promises and Limits of
Multimodal Integration for Deconvolution: The HADACA3 Benchmark’.
*NeurIPS*.

Zappia, Luke, Belinda Phipson, and Alicia Oshlack. 2017. ‘Splatter:
Simulation of Single-Cell RNA Sequencing Data’. *Genome Biology* 18:
174. <https://doi.org/10.1186/s13059-017-1305-0>.
