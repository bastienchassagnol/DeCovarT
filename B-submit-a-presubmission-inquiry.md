# NA

Submitting Author Name: Bastien Chassagnol Submitting Author Github
Handle: @bastienchassagnol Other Package Authors Github handles: (comma
separated, delete if none) Repository:
<https://github.com/bastienchassagnol/DeCovarT> Submission type:
Pre-submission Language: en

> **ORCID:**
> [0000-0002-8955-2391](https://orcid.org/0000-0002-8955-2391)
> **Documentation:** <https://bastienchassagnol.github.io/DeCovarT/>
> **Method preprint:** <https://arxiv.org/abs/2309.09557>

------------------------------------------------------------------------

- Paste the full DESCRIPTION file inside a code block below:

&nbsp;

    Package: DeCovarT
    Title: Covariance-Aware Deconvolution of Bulk Transcriptomic Samples
    Version: 2.0.1
    Authors@R:
        person("Bastien", "Chassagnol",
               email = "bastien_chassagnol@laposte.net",
               role = c("aut", "cre", "cph"),
               comment = c(
                 ORCID = "0000-0002-8955-2391",
                 LinkedIn = "https://www.linkedin.com/in/bastien-chassagnol-677b67140/?locale=fr"
               ))
    Description: Estimates cell-type proportions in bulk transcriptomic samples
        with a probabilistic convolution model that integrates the gene-gene
        covariance structure of purified reference profiles. Cellular ratios
        are recovered by maximum likelihood under a multivariate Gaussian
        convolution, using analytic gradients and Hessians, an additive
        log-ratio reparametrisation that enforces the simplex constraint, and
        Marquardt-Levenberg or Newton-type optimisers. The methodology is
        described in Chassagnol, Nuel and Becht (2023)
        <doi:10.48550/arXiv.2309.09557>.
    License: MIT + file LICENSE
    Additional_repositories: https://bioconductor.org/packages/release/bioc
    Encoding: UTF-8
    Roxygen: list(markdown = TRUE)
    Suggests:
        compositions,
        flextable,
        ggplot2,
        knitr,
        litedown,
        numDeriv,
        pkgdown,
        quarto,
        reactable (>= 0.4.4),
        readr,
        rmarkdown,
        stringr,
        testthat (>= 3.0.0),
        tinytable,
        withr
    VignetteBuilder: knitr, quarto
    Config/testthat/edition: 3
    Imports:
        circlize,
        ComplexHeatmap,
        dplyr,
        e1071,
        igraph,
        limSolve,
        marqLevAlg,
        MASS,
        methods,
        Metrics,
        MixSim,
        nnls,
        glmnet,
        purrr,
        Rdpack (>= 0.7),
        rlang,
        tensor,
        tibble,
        tidyr,
        viridis
    Depends:
        R (>= 4.1.0)
    RdMacros: Rdpack
    Language: en-GB
    URL: https://github.com/bastienchassagnol/DeCovarT, https://bastienchassagnol.github.io/DeCovarT/
    BugReports: https://github.com/bastienchassagnol/DeCovarT/issues
    Config/roxygen2/version: 8.0.0
    RoxygenNote: 7.3.3

> *DeCovarT* estimates cell-type proportions in bulk RNA-seq samples
> using a probabilistic convolution model that incorporates the
> gene–gene covariance structure of cellular profiles. Ratios are
> recovered by maximum likelihood under a multivariate Gaussian
> convolution, with analytic gradients and Hessians, additive log-ratio
> reparametrisation on the simplex, and Marquardt–Levenberg or
> Newton-type optimisers. The package also ships classical deconvolution
> baselines (OLS, NNLS, DeconRNASeq-style QP, robust linear regression,
> CIBERSORT-style ν-SVR) for benchmarking and evaluation utilities.

## Scope

Please indicate which category or categories from our [package fit
policies](https://ropensci.github.io/dev_guide/policies.html#package-categories)
or [statistical package
categories](https://stats-devguide.ropensci.org/overview.html#overview-categories)
this package falls under. (Please check one or more appropriate boxes
below):

**Data Lifecycle Packages**

data retrieval

data extraction

data munging

data deposition

data validation and testing

workflow automation

version control

citation management and bibliometrics

scientific software wrappers

field and lab reproducibility tools

database software bindings

geospatial data

translation

**Statistical Packages**

Bayesian and Monte Carlo Routines

Dimensionality Reduction, Clustering, and Unsupervised Learning

Machine Learning

Regression and Supervised Learning

Exploratory Data Analysis (EDA) and Summary Statistics

Spatial Analyses

Time Series Analyses

Probability Distributions

Explain how and why the package falls under these categories (briefly,
1-2 sentences). Please note any areas you are unsure of:

> *I believe DeCovarT falls under **statistical software** rather than
> the general data-lifecycle categories above.* The current release fits
> **Regression and Supervised Learning** (constrained MLE of cell-type
> proportions on the simplex, mapping bulk transcriptomic mixtures to
> reference signatures) and **Probability Distributions** (multivariate
> Gaussian convolution with cell-type-specific covariance structure,
> including the derivation of asymptotic confidence intervals from the
> expected Fisher information matrix).
>
> *I am aware that **Probability Distributions** is currently marked
> “not yet in scope” in the statistical software guide. I have ticked it
> as the most accurate description of a core component of the method,
> but defer to the editors on whether it should be treated as a
> secondary descriptor or deferred.*
>
> *A planned extension will model latent, unobserved realisations of
> stochastic cellular profiles \boldsymbol{X} in a Bayesian framework
> (posterior inference over cell-type compositions), which would
> additionally fall under **Bayesian and Monte Carlo Routines**. I have
> not checked that box for the current release, but would welcome
> guidance on whether to anticipate it here.*

- If submitting a statistical package, have you already [incorporated
  documentation of standards into your code via the **srr**
  package](https://stats-devguide.ropensci.org/pkgdev.html#pkgdev-srr)?

> *Partially.* The `srr::srr_stats_roclet` is enabled in `DESCRIPTION`,
> and `R/srr-stats-standards.R` documents all 116 **Regression and
> Supervised Learning** standards (v0.2.0): **66** tagged `@srrstats`
> (addressed), **1** `@srrstatsNA` (G4.0, no file-writing wrappers), and
> **49** `@srrstatsTODO` (to be completed during review). Running
> `srr::srr_stats_pre_submit()` reports **57% compliance**, above the
> ≥50% threshold for full submission. Representative standards already
> met include:
>
> - *Input/output checks and documentation:* `checkmate` validation and
>   roxygen2 parameter documentation (G2.0–G2.3b).
> - *Model specification and convergence:* Marquardt–Levenberg / Newton
>   optimisers with convergence diagnostics and configurable thresholds
>   (RE3.0–RE3.3).
> - *Uncertainty quantification:* asymptotic intervals from the Fisher
>   information matrix (RE4.3, RE4.6).
> - *Testing:* correctness, edge-case, and parameter-recovery tests with
>   fixed seeds (G5.0–G5.9a, RE7.0–RE7.3).
>
> Remaining gaps include missing-data handling options (G2.14), extended
> test suites (G5.10–G5.12), goodness-of-fit reporting (RE4.10–RE4.11),
> and distributing `@srrstats` tags from the central standards file into
> individual `R/` source files — all tractable during review.

- Who is the target audience and what are scientific applications of
  this package?

> **Target audience:** computational biologists and immunologists
> working with bulk RNA-seq who need reproducible cell-type
> deconvolution that accounts for gene–gene covariance in reference
> profiles.
>
> **Scientific applications:** estimating cell composition in bulk
> transcriptomes; comparative benchmarking of deconvolution methods; and
> realistic synthetic validation scenarios that model gene interactions
> via random network generators (see the
> [synthetic-scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/synthetic-scenarios.html)
> vignette).

- Are there other R packages that accomplish the same thing? If so, how
  does yours differ or meet [our criteria for
  best-in-category](https://ropensci.github.io/dev_guide/policies.html#overlap)?

> Several R/Bioconductor packages address bulk deconvolution, including
> reference-based methods such as **MuSiC**, **BayesPrism**,
> **BisqueRNA**, **SCDC**, **InstaPrism**, and signature-based tools
> (e.g. **EPIC**, **CIBERSORT** implementations). A recent benchmarking
> study ([Dietrich et al., *Genome Biology*
> 2026](https://doi.org/10.1186/s13059-026-03955-w)) provides a unifying
> framework (`omnideconv`) for second-generation methods that build
> signatures from single-cell profiles rather than cell-population
> averages; the associated
> [`deconvBench`](https://github.com/omnideconv/deconvbench) pipeline
> standardises comparisons across tumour microenvironments,
> developmental cell lines, and organoid settings.
>
> *DeCovarT differs from all of the above by:*
>
> 1.  **Modelling:** multivariate Gaussian convolution MLE with explicit
>     covariance structure, whereas most standard approaches assume
>     independence between genes (including all tools currently
>     benchmarked in `omnideconv`).
> 2.  **Optimisation:** constrained MLE with analytic derivatives and
>     simplex reparametrisation, enabling asymptotic intervals from an
>     explicit expected Fisher information matrix — a form of
>     uncertainty quantification absent from most existing tools.
>
> *I plan to run a comprehensive benchmark comparing DeCovarT against
> second-generation single-cell-informed methods across diverse
> biological settings (tumour microenvironment, organoids) using the
> `deconvBench` framework. A key challenge for genuine comparison is
> that most existing algorithms rely solely on mean expression profiles,
> whereas DeCovarT additionally requires a covariance structure — making
> performance strongly dependent on the quality of inferred gene
> regulatory networks and covariance matrices. Reconstructing those
> networks accurately remains an open problem: expression data alone
> largely yields correlative rather than causal edges, and adding
> chromatin accessibility or other modalities is typically required for
> mechanistic inference ([Maizels & Briscoe, *Nature Reviews Genetics*
> 2026](https://doi.org/10.1038/s41576-026-00939-1)); furthermore, even
> with richer inputs, GRN inference can fail systematically under common
> regimes ([Arnold et al. 2025,
> arXiv:2605.04930](https://arxiv.org/abs/2605.04930)). I will address
> this dependency explicitly in the benchmark design and manuscript.*
>
> *I am happy to expand the README “similar packages” section if editors
> consider that helpful before a full submission.*

- (If applicable) Does your package comply with our [guidance around
  *Ethics, Data Privacy and Human Subjects
  Research*](https://devguide.ropensci.org/policies.html#ethics-data-privacy-and-human-subjects-research)?

> *The package operates on user-supplied expression matrices and
> reference signatures; it does not collect, store, or distribute
> human-subject data. I am not aware of ethics or privacy issues
> specific to the software itself, but I will follow any additional
> guidance editors recommend for applied vignettes that use public
> cohorts.*

- Any other questions or issues we should be aware of?:

> **Publication and review timing.** I plan to submit a **methods
> article to *Nature Methods*** describing the methodology and
> validation. That manuscript concerns the **scientific method**, not a
> separate software paper about the R package itself.
>
> - The package is **not** currently under review at another
>   software-review venue.
> - The *Nature Methods* manuscript is **not yet submitted**; I intend
>   to pursue rOpenSci software review **first**, in line with the
>   [author
>   guide](https://devguide.ropensci.org/softwarereview_author.html),
>   and will not open the journal submission while an rOpenSci review is
>   active.
>
> *Please advise whether this sequencing is acceptable.*

> **Package readiness (`pkgcheck`).** All **mandatory** `pkgcheck` items
> pass, including CI ([latest
> run](https://github.com/bastienchassagnol/DeCovarT/actions/runs/31252824642)),
> **86.1%** test coverage, and clean `R CMD check`. Remaining optional
> `:eyes:` items include a relatively large number of Imports (20),
> mostly from the tidyverse ecosystem rather than base R; I can address
> these during review if the package is in scope. **Lifecycle:**
> [Active](https://www.repostatus.org/#active) / stable — **License:**
> MIT — **Maintainer commitment:** I expect to maintain the package for
> ≥2 years.

### Miscellaneous questions

> *The following points address the editor’s follow-up questions from
> [issue
> \#798](https://github.com/ropensci/software-review/issues/798#issuecomment-5226685188).*

**`srr` compliance estimate.** As detailed above under the `srr`
question, I have read the **Regression and Supervised Learning**
standards and estimate compliance with well over half of them in the
current codebase. The principal remaining gaps are systematic
missing-data handling and goodness-of-fit reporting for the convolution
model; both are tractable and I commit to addressing them during review.

**Scope of planned benchmarking.** I plan a comprehensive benchmark of
DeCovarT against first- and second-generation deconvolution methods,
including those unified by
[`omnideconv`](https://pmc.ncbi.nlm.nih.gov/articles/PMC12837286/)
([Dietrich et al., *Genome Biology*
2026](https://doi.org/10.1186/s13059-026-03955-w)), using the
[`deconvBench`](https://github.com/omnideconv/deconvbench) Nextflow
pipeline as a standardised evaluation platform. Planned biological
contexts include tumour microenvironment mixtures and developmental
settings such as organoids.

*A critical caveat:* virtually all competing methods require only mean
expression profiles per cell type, whereas DeCovarT additionally
requires a covariance (or precision) matrix encoding gene–gene
interactions. Method performance therefore depends jointly on
deconvolution accuracy *and* on how well the input covariance is
reconstructed from scRNA-seq data or inferred via gene regulatory
network (GRN) methods. Recovering accurate GRNs is an active and
challenging research problem: expression data alone largely yields
correlative rather than causal regulatory edges, and adding chromatin
accessibility (ATAC-seq) or other modalities is typically required for
mechanistic inference ([Maizels & Briscoe, *Nature Reviews Genetics*
2026](https://doi.org/10.1038/s41576-026-00939-1)); furthermore, GRN
inference can fail systematically even with richer inputs ([Arnold et
al. 2025, arXiv:2605.04930](https://arxiv.org/abs/2605.04930)). An open
question — which the benchmark will address explicitly — is whether
recovering the true causal regulatory architecture is *necessary* for
improved compositional recovery, or whether statistically consistent
covariance estimation alone is sufficient for the gains in precision
that DeCovarT targets. The benchmark design will evaluate performance
across GRN reconstruction strategies of varying quality, and I intend to
be transparent about this dependency in both the benchmark vignette and
the manuscript.

Beyond transcriptomics, the same Gaussian-convolution framework extends
naturally to **mixture experimental design** for compositional data —
Scheffé-type models with simplex-constrained proportions — with
potential applications in formulation chemistry, food science, and other
domains where component blending is the design variable (outlined in the
[DeCovarT
perspectives](https://bastienchassagnol.github.io/DeCovarT/articles/decovart-perspectives.html)
vignette, using `AlgDesign` / `skpr` for designed simplex points).

**Questions for the editors.**

1.  Should **Probability Distributions** be treated as a secondary
    descriptor for the Gaussian convolution component, or should the
    submission focus on **Regression and Supervised Learning** only?
2.  Should I check **Bayesian and Monte Carlo Routines** now (given the
    planned latent-profile extension), or only once that API ships?

## Use of Generative AI

Generative AI tools were used to produce some of the material in this
submission.

If so, please describe usage, and include links to any relevant aspects
of your repository. See [our blog
post](https://ropensci.org/blog/2026/02/26/ropensci-ai-policy/) for
background. (Explicit advice is not yet included in our *Dev Guide*; we
are hoping to update very soon, and ask your cooperation and
transparency in the meantime.)

> *As a non-native English speaker, I have used generative AI tools
> primarily to refine grammar and scientific prose in documentation,
> vignettes, and this inquiry — not to generate statistical results or
> undocumented code paths.* In the repository, project rules in
> [`.cursorrules`](https://github.com/bastienchassagnol/DeCovarT/blob/main/.cursorrules)
> encode good-practice expectations (comprehensive roxygen2
> documentation, systematic testing of core functionality, and
> linting/formatting conventions). *AI-assisted editing is distinct from
> AI-generated methodology; all estimators, tests, and numerical claims
> are implemented and verified in R.*
