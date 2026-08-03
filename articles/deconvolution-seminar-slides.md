# Doctorants seminar slides (embedded)

## Overview

This page embeds the Quarto reveal.js deck built from the doctorants
deconvolution seminar (sources under `slides/` in the repository).
Pattern inspired by
[GTBioSS_slides.qmd](https://github.com/bastienchassagnol/COPASI_Team216_Lo2016/blob/master/vignettes/GTBioSS_slides.qmd)
and the [quarto-webr slide-embed
example](https://quarto-webr.thecoatlessprofessor.com/examples/book/slide-embed.html#presentation).

[Open the deck full
screen](https://bastienchassagnol.github.io/DeCovarT/slides/index.html).

# Cellular deconvolution with DeCovarT

A network-centric view of bulk mixtures

Bastien Chassagnol

MMG / LPSM / INSERM

Gregory Nuel

Etienne Becht

Anaïs Baudot

2026-08-03

## 

Cellular deconvolution

DeCovarT: a network-centric approach to bulk mixtures

------------------------------------------------------------------------

Bastien Chassagnol  
*MMG*

Gregory Nuel  
*LPSM*

Etienne Becht  
*INSERM*

Anaïs Baudot  
*MMG*

Slides:
[GitHub](https://github.com/bastienchassagnol/DeCovarT/tree/main/slides)
· Full screen after render: `slides/index.html`

## Outline

1.  Why deconvolution?
2.  A review of deconvolution techniques
3.  DeCovarT: a network-centric approach
4.  Biological applications and perspectives
5.  Quiz

# Why deconvolution?

## Quantify components of biological systems

Bulk tissues mix cell populations. Autoimmunity and cancer change both
**composition** and **state** of those populations.

Schnell et al., *Nature*, 2020 —
[doi:10.1038/s41422-020-0277-x](https://doi.org/10.1038/s41422-020-0277-x)

## Holistic complexity

Genes cooperate in pathways; metabolites and cell populations interact
in bio-molecular networks. Deconvolution alone is not enough — we need
**structure** on the gene side as well.

Schnell et al., *Nature*, 2020 —
[doi:10.1038/s41422-020-0277-x](https://doi.org/10.1038/s41422-020-0277-x)

## Why bulk is limited for biomarkers

Shoemaker et al. (2012): the same bulk fold-change can come from

- **Scenario A** — activation within an existing cell type
- **Scenario B** — arrival of a new cell population

Shoemaker et al., *BMC Genomics*, 2012 —
[doi:10.1186/1471-2164-13-460](https://doi.org/10.1186/1471-2164-13-460)

## Technical alternatives

**Bulk RNA-seq**

- Cheaper, larger cohorts
- Longitudinal studies
- Aggregates heterogeneous signals

**Single-cell / spatial**

- Resolve rare types & lineages
- Costly, noisy, heavy compute
- Complements bulk rather than replaces it

Deconvolution recovers cell-type proportions from bulk when single-cell
is impractical at scale.

# Review of techniques

## The ecosystem of algorithms

Shen-Orr & Gannett (2013) taxonomy:

- **Partial deconvolution** — estimate proportions \\p\\ given signature
  \\X\\
- **Complete deconvolution** — jointly infer \\p\\ and \\X\\ (ill-posed
  without priors)

## Granularity

Deconvolution targets vary:

- Mixtures of **cell populations** (e.g. immune fractions)
- Mixtures of **tissues** (e.g. tumour purity)
- Mixtures of **cell-cycle phases**

## Reference-based principle

1.  Purified cellular profiles (signature matrix \\X\\)
2.  Bulk mixture \\y\\
3.  Estimate proportions \\p\\ (often with simplex constraints)

Classic tools: DeconRNASeq, CIBERSORT, EPIC, quanTIseq, FARDEEP, …

## Signature-based families

**Regression**

- Quadratic programming (DeconRNASeq)
- Robust / SVR (CIBERSORT, FARDEEP)
- Weighted QP (EPIC, quanTIseq)

**Probabilistic**

- Continuous: DSection, DeMix
- Discrete: PERT, TEMT

# DeCovarT

## Global pipeline

1.  Collect / curate cell-type datasets
2.  Learn features & GRN structure
3.  **DeCovarT** estimation of \\p\\
4.  Biological exploitation

Pipeline framing inspired by Avila Cobos et al., *Genome Biology*, 2018
—
[doi:10.1186/s13059-018-1479-0](https://doi.org/10.1186/s13059-018-1479-0)

## Multivariate Gaussian convolution

Cell type \\j\\: \\x_j \sim \mathcal{N}\_G(\mu_j, \Sigma_j)\\.

Bulk observation:

\\ y \mid p \\\sim\\ \mathcal{N}\_G\\\Bigl( \mu p,\\ \sum_j
p_j^{2}\\\Sigma_j \Bigr). \\

Covariance structure enters the likelihood — gene–gene dependence is
first-class, not noise.

## Unit-simplex constraint

Proportions live on \\\Delta^{J-1}\\. DeCovarT uses an unconstrained map
\\\rho \mapsto p\\ (additive logistic / soft-max) that is a
\\C^{2}\\-diffeomorphism, then optimises with Marquardt–Levenberg.

## Toy model (2 cells × 2 genes)

Bivariate convolutions make the geometry visible: entropy of \\p\\,
overlap of \\\Sigma_j\\, and solver behaviour (DeconRNASeq vs DeCovarT).

See package vignettes [synthetic
scenarios](../vignettes/synthetic-scenarios.qmd) and
[`benchmark_bivariate_gaussian_convolutions()`](https://bastienchassagnol.github.io/DeCovarT/reference/benchmark_bivariate_gaussian_convolutions.html).

# Applications & perspectives

## Sjögren’s disease

Transcriptomic patient clusters and IFN modules motivate
composition-aware analysis of immune infiltration.

Soret et al., *Nature Communications*, 2021 —
[doi:10.1038/s41467-021-23949-4](https://doi.org/10.1038/s41467-021-23949-4)

## DTOO / gastruloid perspectives

Digital twins and multimodal gastruloid maps need deconvolution across
bulk, single-cell, and ATAC layers — composition + covariance again.

## Spatial deconvolution (outlook)

Emerging toolkits mix:

- Graph / optimal transport pairing of cells to spots
- Bayesian count models (NB, Poisson–gamma, topic models)
- Regression with simplex / robust / ensemble constraints

Gaspard-Boulinc & Cavalli, *Nature Reviews Genetics*, 2025 —
[doi:10.1038/s41576-025-00845-y](https://doi.org/10.1038/s41576-025-00845-y)

## Perspectives for DeCovarT

- Rare / unknown cell content and spillover
- Multi-omics and spatio-temporal resolution
- Hierarchical ontologies (coarse vs fine labels)
- Compositional statistics on \\p\\

# Quiz

## What does reference-based deconvolution estimate?

- Gene regulatory network topology only
- Cell-type proportions \\p\\ given a signature \\X\\
- Single-cell UMAP embeddings
- Batch-correction factors alone

## Which statement about DeCovarT is true?

- It ignores gene–gene covariance
- It models bulk as a Gaussian convolution with covariance \\\sum_j
  p_j^2\Sigma_j\\
- It only works on DNA methylation
- It requires spatial coordinates for every cell

## Scenario A vs B (Shoemaker 2012)

A bulk gene increases. Scenario B means:

- The same cell type up-regulates the gene
- A new cell population arrives in the mixture
- Sequencing depth doubled
- The gene was a spike-in control

# Acknowledgements

## Thanks

**Supervision & co-authors**

- Anaïs Baudot
- Gregory Nuel
- Etienne Becht
- Marielle Péré

**Collaborators**

- SysBioMed / MMG colleagues
- Doctorants seminar organisers

## Take-away messages

- Bulk mixes cell types; deconvolution recovers proportions \\p\\.
- Reference-based methods need a signature \\X\\; unsupervised methods
  need strong priors.
- DeCovarT models gene–gene covariance via multivariate Gaussian
  convolutions.
- Always validate with aliquots (bulk + FACS / scRNA-seq) when possible.

DeCovarT — doctorants seminar
