# External resources: manuscript, presentations, and preprints

This page collects **material outside the package vignette tutorials**:
the working methods manuscript, seminar slides, and canonical URLs. The
in-package articles under *Vignettes* document implementation,
simulation, and mathematical detail; start here when you need the paper
PDF or talk deck.

| Resource | Link |
|:---|:---|
| Source repository | [bastienchassagnol/DeCovarT](https://github.com/bastienchassagnol/DeCovarT) |
| Methods preprint | [arXiv:2309.09557](https://arxiv.org/abs/2309.09557) ([Chassagnol et al. 2023](#ref-chassagnolDeCovarTMultidimensionalProbalistic2023)) |
| Package site | [bastienchassagnol.github.io/DeCovarT](https://bastienchassagnol.github.io/DeCovarT/) |
| Function call network | [decovart_function_network.html](https://bastienchassagnol.github.io/DeCovarT/package_network/decovart_function_network.html) |

## Methods manuscript

The manuscript describes covariance-aware deconvolution of bulk
transcriptomes as a multivariate Gaussian convolution. Cell-type
proportions lie on the simplex; they are recovered by constrained
maximum likelihood with analytic gradients and Hessians. The PDF below
is the working build of `article/main.pdf` (Springer Nature `sn-jnl`
template).

Open [full
screen](https://bastienchassagnol.github.io/DeCovarT/article/main.pdf)
or [download a
copy](https://bastienchassagnol.github.io/DeCovarT/articles/decovart-manuscript.pdf).

## Doctorants seminar slides

This embeds the Quarto reveal.js deck from the doctorants deconvolution
seminar (sources under `slides/` in the repository). Pattern inspired by
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

2026-08-31

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
**composition** and **state** of those populations ([Schnell et al.
2020](#/references)).

Figure 1: Immune balance as a seesaw between regulation and
inflammation.

Schnell et al. ([2020](#/references))

## Immune balance (schematic)

Complementary schematic of [Figure 1](#/fig-immune-balance):

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TB
  H["Balanced<br/>Regulation ≈ Inflammation<br/>Healthy"]
  U["Under-reaction<br/>Regulation ≫ Inflammation<br/>Cancer risk"]
  O["Over-reaction<br/>Inflammation ≫ Regulation<br/>Autoimmunity"]
  H -->|"imbalance?"| U
  H -->|"imbalance?"| O
```

Figure 2: Healthy balance vs under- and over-reaction.

Schnell et al. ([2020](#/references))

## Holistic complexity

Genes cooperate in pathways; metabolites and cell populations interact
in bio-molecular networks. Deconvolution alone is not enough — we need
**structure** on the gene side as well ([Schnell et al.
2020](#/references)).

Figure 3: Holistic view: transcriptome and cell populations feed immune
balance.

Schnell et al. ([2020](#/references))

## Why bulk is limited for biomarkers

Shoemaker et al. ([2012](#/references)): the same bulk fold-change can
come from

- **Scenario A** — activation within an existing cell type
  1.  
- **Scenario B** — arrival of a new cell population
  ([Figure 5](#/fig-bulk-b))

Figure 4: Scenario A: within-type activation.

Figure 5: Scenario B: new population.

Shoemaker et al. ([2012](#/references))

## Technical alternatives

Figure 6: Bulk versus single-cell trade-offs.

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

Taxonomy after Shen-Orr & Gannett (2013) and reviews such as Avila Cobos
et al. ([2018](#/references)), Gaspard-Boulinc et al.
([2025](#/references)):

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

Figure 7: Signature \\X\\, bulk \\y\\, proportions \\p\\.

1.  Purified cellular profiles (signature matrix \\X\\)
2.  Bulk mixture \\y\\
3.  Estimate proportions \\p\\ (often with simplex constraints)

Classic tools include DeconRNASeq ([Gong and Szustakowski
2013](#/references)), CIBERSORT ([Newman et al. 2015](#/references)),
and EPIC ([Racle et al. 2017](#/references)).

Avila Cobos et al. ([2018](#/references))

## Signature-based families

Compact survey of common reference-based tools
([Table 1](#/tbl-methods)):

Table 1: Selected reference-based deconvolution methods.

**Regression**

- Quadratic programming (DeconRNASeq)
- Robust / SVR (CIBERSORT, FARDEEP)
- Weighted QP (EPIC, quanTIseq)

**Probabilistic**

- Continuous: DSection, DeMix
- Discrete: PERT, TEMT

# DeCovarT

## Global pipeline

Figure 8: From references to biological exploitation.

1.  Collect / curate cell-type datasets
2.  Learn features & GRN structure
3.  **DeCovarT** estimation of \\p\\
4.  Biological exploitation

Avila Cobos et al. ([2018](#/references))

## Multivariate Gaussian convolution

Cell type \\j\\: \\x_j \sim \mathcal{N}\_G(\mu_j, \Sigma_j)\\.

Bulk observation ([Equation 1](#/eq-bulk-gauss)):

\\ y \mid p \\\sim\\ \mathcal{N}\_G\\\Bigl( \mu p,\\ \sum_j
p_j^{2}\\\Sigma_j \Bigr). \tag{1}\\

- \\G\\: number of genes (dimension of \\y\\)
- \\J\\: number of cell types; \\p \in \Delta^{J-1}\\
- \\\mu\\: \\G \times J\\ mean signature; \\\Sigma_j\\: gene–gene
  covariance of type \\j\\

Covariance structure enters the likelihood — gene–gene dependence is
first-class, not noise.

Figure 9: Another view of the convolution principle.

## Unit-simplex constraint

Proportions live on \\\Delta^{J-1}\\. DeCovarT uses an unconstrained map
\\\rho \mapsto p\\ (additive logistic / soft-max) that is a
\\C^{2}\\-diffeomorphism, then optimises with Marquardt–Levenberg.

## Toy model (2 cells × 2 genes)

Bivariate convolutions make the geometry visible: entropy of \\p\\,
overlap of \\\Sigma_j\\, and solver behaviour (DeconRNASeq vs DeCovarT).
See [Figure 10](#/fig-toy-model) and [Figure 11](#/fig-toy-results).

Figure 10: Toy bivariate setup.

Figure 11: Toy estimation results.

Package vignettes: [synthetic
scenarios](../vignettes/synthetic-scenarios.qmd) and
[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.html)
with scenario grids from `scripts/configure_bivariate_toy_scenarios.R`.

# Applications & perspectives

## Sjögren’s disease

Transcriptomic patient clusters and IFN modules motivate
composition-aware analysis of immune infiltration ([Soret et al.
2021](#/references)).

Figure 12: Sjögren patient stratification motivating deconvolution.

Soret et al. ([2021](#/references))

## DTOO / gastruloid perspectives

Digital twins and multimodal gastruloid maps need deconvolution across
bulk, single-cell, and ATAC layers — composition + covariance again.

## Spatial deconvolution (outlook)

Emerging toolkits mix ([Gaspard-Boulinc et al. 2025](#/references)):

- Graph / optimal transport pairing of cells to spots
- Bayesian count models (NB, Poisson–gamma, topic models)
- Regression with simplex / robust / ensemble constraints

Gaspard-Boulinc et al. ([2025](#/references))

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
  need strong priors ([Table 1](#/tbl-methods)).
- DeCovarT models gene–gene covariance via multivariate Gaussian
  convolutions ([Equation 1](#/eq-bulk-gauss)).
- Always validate with aliquots (bulk + FACS / scRNA-seq) when possible.

## References

Avila Cobos, Francisco, Jo Vandesompele, Pieter Mestdagh, and Katleen De
Preter. 2018. ‘Computational Deconvolution of Transcriptomics Data from
Mixed Cell Populations’. *Bioinformatics (Oxford, England)* 34.
<https://doi.org/10.1093/bioinformatics/bty019>.

Gaspard-Boulinc, Lucie C., Luca Gortana, Thomas Walter, Emmanuel
Barillot, and Florence M. G. Cavalli. 2025. ‘Cell-Type Deconvolution
Methods for Spatial Transcriptomics’. *Nature Reviews Genetics* 26.
<https://doi.org/10.1038/s41576-025-00845-y>.

Gong, Ting, and Joseph D. Szustakowski. 2013. ‘DeconRNASeq: A
Statistical Framework for Deconvolution of Heterogeneous Tissue Samples
Based on mRNA-Seq Data’. *Bioinformatics (Oxford, England)* 29.
<https://doi.org/10.1093/bioinformatics/btt090>.

Newman, Aaron, Chih Liu, Michael Green, et al. 2015. ‘Robust Enumeration
of Cell Subsets from Tissue Expression Profiles’. *Nature Methods* 12.
<https://doi.org/10.1038/nmeth.3337>.

Racle, Julien, Kaat de Jonge, Petra Baumgaertner, Daniel E Speiser, and
David Gfeller. 2017. ‘Simultaneous Enumeration of Cancer and Immune Cell
Types from Bulk Tumor Gene Expression Data’. *eLife* 6.
<https://doi.org/10.7554/elife.26476>.

Schnell, Alexandra, Lloyd Bod, Asaf Madi, and Vijay K. Kuchroo. 2020.
‘The Yin and Yang of Co-Inhibitory Receptors: Toward Anti-Tumor Immunity
Without Autoimmunity’. *Cell Research* 30.
<https://doi.org/10.1038/s41422-020-0277-x>.

Shoemaker, Jason E., Tiago JS Lopes, Samik Ghosh, Yukiko Matsuoka,
Yoshihiro Kawaoka, and Hiroaki Kitano. 2012. ‘CTen: A Web-Based Platform
for Identifying Enriched Cell Types from Heterogeneous Microarray Data’.
*BMC Genomics* 13. <https://doi.org/10.1186/1471-2164-13-460>.

Soret, Perrine, Christelle Le Dantec, Emiko Desvaux, et al. 2021. ‘A New
Molecular Classification to Drive Precision Treatment Strategies in
Primary Sjögren’s Syndrome’. *Nature Communications* 12.
<https://doi.org/10.1038/s41467-021-23949-4>.

DeCovarT — doctorants seminar

## References

Chassagnol, Bastien, Grégory Nuel, and Etienne Becht. 2023. *DeCovarT, a
Multidimensional Probabilistic Model for the Deconvolution of
Heterogeneous Transcriptomic Samples*. arXiv.
<https://doi.org/10.48550/arxiv.2309.09557>.
