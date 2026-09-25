Submitting Author Name: Bastien Chassagnol
Submitting Author Github Handle: <!--author1-->@bastienchassagnol<!--end-author1-->
Other Package Authors Github handles: (comma separated, delete if none) <!--author-others-->[@mariellePr](https://github.com/mariellePr)<!--end-author-others-->
Repository:  <!--repourl-->https://github.com/bastienchassagnol/DeCovarT<!--end-repourl-->
Version submitted: 2.3.1
Submission type: <!--submission-type-->Stats<!--end-submission-type-->
Badge grade: <!--statsgrade-->silver<!--end-statsgrade-->
Editor: <!--editor--> TBD <!--end-editor-->
Reviewers: <!--reviewers-list--> TBD <!--end-reviewers-list-->
<!--due-dates-list--><!--end-due-dates-list-->
Archive: TBD
Version accepted: TBD
Language: <!--language-->en<!--end-language-->

This follows the approved pre-submission inquiry [#798](https://github.com/ropensci/software-review/issues/798).

-   Paste the full DESCRIPTION file inside a code block below:

```
Package: DeCovarT
Title: Covariance-Aware Deconvolution of Bulk Transcriptomic Samples
Version: 2.3.1
Authors@R:
    person("Bastien", "Chassagnol", , "bastien_chassagnol@laposte.net", role = c("aut", "cre", "cph"),
           comment = c(ORCID = "0000-0002-8955-2391"))
Description: Estimates cell-type proportions in bulk transcriptomic
    samples with a probabilistic convolution model that integrates the
    gene-gene covariance structure of purified reference profiles.
    Cellular ratios are recovered by maximum likelihood under a
    multivariate Gaussian convolution, using analytic gradients and
    Hessians, an isometric log-ratio reparametrisation (Helmert basis)
    that enforces the simplex constraint, and Marquardt-Levenberg or
    Newton-type optimisers. The methodology is described in Chassagnol,
    Nuel and Becht (2023) <doi:10.48550/arXiv.2309.09557>.
License: MIT + file LICENSE
URL: https://github.com/bastienchassagnol/DeCovarT,
    https://bastienchassagnol.github.io/DeCovarT/
BugReports: https://github.com/bastienchassagnol/DeCovarT/issues
Depends:
    R (>= 4.1.0)
Imports:
    dplyr,
    ggplot2 (>= 3.5.0),
    marqLevAlg,
    MASS,
    Matrix,
    purrr,
    Rdpack (>= 0.7),
    rlang,
    tensor,
    tibble,
    tidyr
Suggests:
    circlize,
    cli,
    ComplexHeatmap,
    compositions,
    cowplot,
    e1071,
    EMMIXmfa,
    flextable,
    forcats,
    funkyheatmap,
    furrr,
    future,
    ggdendro,
    ggdist,
    ggtext,
    glmnet,
    gridExtra,
    htmlwidgets,
    igraph,
    knitr,
    limSolve,
    litedown,
    mclust,
    Metrics,
    MixSim,
    mvtnorm,
    nnls,
    numDeriv,
    patchwork,
    pkgdown,
    plotly,
    png,
    ps,
    qrng,
    quarto,
    qqplotr,
    reactable (>= 0.4.4),
    readr,
    rgl,
    rmarkdown,
    spelling,
    stringr,
    testthat (>= 3.0.0),
    tinytable (>= 0.15.0),
    viridis,
    withr
VignetteBuilder: knitr, quarto
RdMacros: Rdpack
Config/roxygen2/version: 8.1.0
Config/testthat/edition: 3
Encoding: UTF-8
Language: en-GB
Roxygen: list(markdown = TRUE, roclets = c("namespace", "rd",
    "srr::srr_stats_roclet"))
RoxygenNote: 8.1.0
```

## Scope

- Please indicate which of our [statistical package categories](https://stats-devguide.ropensci.org/overview.html#overview-categories) this package falls under.

     **Statistical Packages**

	- [ ] Bayesian and Monte Carlo Routines
	- [ ] Dimensionality Reduction, Clustering, and Unsupervised Learning
	- [ ] Machine Learning
	- [x] Regression and Supervised Learning
	- [ ] Exploratory Data Analysis (EDA) and Summary Statistics
	- [ ] Spatial Analyses
	- [ ] Time Series Analyses
	- [x] Probability Distributions

## Pre-submission Inquiry

- [x] A pre-submission inquiry has been approved in [issue#798](https://github.com/ropensci/software-review/issues/798)

## General Information

-   Who is the target audience and what are scientific applications of this package?

Computational biologists and immunologists who estimate cell-type proportions from bulk RNA-seq when a purified or single-cell reference is available. Typical uses are covariance-aware deconvolution of mixed tissues, simulation studies that need a known Gaussian convolution, and method comparison against mean-only baselines that ship with the package (NNLS, simplex QP, GLS, `rlm`, CIBERSORT-style ν-SVR).

Manuscript-scale biological pipelines (feature selection, sparse precision estimation, population alignment, gastruloid case study) live in the companion repository [DeCovarT_reproducibility](https://github.com/bastienchassagnol/DeCovarT_reproducibility), not in the CRAN tarball.

-   Paste your responses to our [*General Standard* **G1.1** here](https://stats-devguide.ropensci.org/standards.html#general-standards):

**The first implementation of a novel algorithm**, with related software documented rather than treated as the same estimator.

DeCovarT maximises a *multivariate* Gaussian convolution likelihood for bulk mixtures, with cell-type-specific covariances entering as \(\sum_j p_j^2 \Sigma_j\), an isometric log-ratio chart on the simplex, and analytic score and Hessian. The closest published statistical method is the *univariate* Gaussian convolution of DSection ([Erkkilä *et al.* 2010](https://doi.org/10.1093/bioinformatics/btq470)); DeMixT ([Wang *et al.* 2018](https://www.bioconductor.org/packages/DeMixT); [GitHub](https://github.com/wwylab/DeMixT)) and [ISOpureR](https://cran.r-project.org/package=ISOpureR) also treat profiles as latent, but none of these is a sparse multivariate convolution MLE on the simplex.

Mean-only or weighted-least-squares engines that share a linear mean and a simplex constraint, but not this likelihood, include [MuSiC](https://github.com/xuranw/MuSiC) / [MuSiC2](https://github.com/Jiaxin-Fan/MuSiC2), [DWLS](https://github.com/dtsoucas/DWLS), [BisqueRNA](https://github.com/cozygene/bisque), [CIBERSORTx](https://cibersortx.stanford.edu), [EPIC](https://github.com/GfellerLab/EPIC), [SCDC](https://github.com/meichendong/SCDC), [BayesPrism](https://github.com/Danko-Lab/BayesPrism), and the spatial analogue [RCTD / spacexr](https://github.com/dmcable/spacexr). Second-generation probabilistic tools that *do* put a law on the bulk, but a different one, are [RNA-Sieve](https://github.com/songlab-cal/rna-sieve) (gene-wise CLT, linear variances), MEAD (errors-in-variables, gene–gene correlation), DECALS (subject-specific \(\Sigma_i\), constrained least squares), and ReDeconv (univariate variances linear in \(p_{ji}\), CLTS scale). The statistical contrast is written out in the pkgdown article [Tracks for a three-layer DeCovarT](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-tracks).

Since [#798](https://github.com/ropensci/software-review/issues/798), the package has also gained:

- Dirichlet / QP / multi-start initialisation and explicit local-vs-global diagnostics (`multistart_decovart()`, `boundary_diagnostics()`), because the finite-sample log-likelihood is not globally concave;
- profile likelihood-ratio intervals, chi-bar-square calibration for active zeros, and parametric / reference bootstrap, so Wald is no longer the only interior interval;
- the Godambe sandwich covariance remains an **outlook** item in the statistical-perspectives vignette, not a shipped estimator;

- the companion reproducibility repository for high-dimensional bottlenecks that do not belong in the package API (feature selection, GGM comparison, population alignment, gastruloid pipeline).

-   (If applicable) Does your package comply with our [guidance around *Ethics, Data Privacy and Human Subjects Research*](https://devguide.ropensci.org/policies.html#ethics-data-privacy-and-human-subjects-research)?

Yes. The package ships only simulated and toy matrices. The companion repository uses publicly released GEO single-cell data (Suppinger gastruloids); no additional human-subjects collection is required to run the package.

## Badging

-    What grade of badge are you aiming for? **silver**

-    If aiming for silver or gold, describe which of the [four aspects](https://stats-devguide.ropensci.org/pkgdev.html#pkgdev-silver) the package fulfils:

1. **Standards beyond a minimal subset.** All applicable G / RE / PD standards are tagged `@srrstats` or justified `@srrstatsNA` (162 / 37 / 0 TODO). G5.7 and RE5.0 are N/A by editorial agreement in #798.
2. **Documentation and testing.** Help pages, theory vignettes (generative model, MLE properties, identifiability, statistical perspectives), S3 `decovart_fit` methods, and testthat coverage of the estimator, intervals, and edge cases (collinear means, active zeros, multi-start).
3. **Generality.** The API is any bulk matrix plus cell-type means and covariances, not one tissue atlas. First-generation mean-only solvers are included for comparison; the gastruloid pipeline is one use case in the companion repository, not the only supported workflow.

(Internal design of the convolution MLE — analytic derivatives, ILR chart, cached factorisation of \(V(p)\) — is a possible fourth aspect; I do not rely on it alone.)

## Technical checks

- [x] I have read the [rOpenSci packaging guide](https://devguide.ropensci.org/building.html).
- [x] I have read the [author guide](https://devguide.ropensci.org/softwarereview_author.html) and I expect to maintain this package for at least 2 years or have another maintainer identified.
- [x] I/we have read the [*Statistical Software Peer Review* Guide for Authors](https://stats-devguide.ropensci.org/pkgdev.html).
- [x] I/we have run [`autotest`](https://github.com/ropensci-review-tools/autotest) checks on the package *(local run; log under `logs/autotest_*.log` — I will paste a summary of failures, if any, when posting)*.
- [x] The [`srr_stats_pre_submit()` function](https://ropensci-review-tools.github.io/srr/reference/srr_stats_pre_submit.html) confirms this package may be submitted.
- [ ] The [`pkgcheck()` function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) confirms this package may be submitted *(coverage was ~54% while manuscript figure-book helpers still lived in `R/`; those helpers have been moved to `scripts/visualisation/`. I will re-run `pkgcheck` before posting.)*.

This package:

- [x] does not violate the Terms of Service of any service it interacts with.
- [x] has a CRAN and OSI accepted license.
- [x] contains a [README with instructions for installing the development version](https://devguide.ropensci.org/building.html#readme).

## Use of Generative AI

- [x] Generative AI tools were used to produce some of the material in this submission.

Cursor / LLMs were used for documentation, SRR tagging, packaging hygiene, and this issue text. The estimator, proofs, and simulation design are the author’s. Project-local programming conventions were enforced via the repository [`.cursorrules`](https://github.com/bastienchassagnol/DeCovarT/blob/main/.cursorrules). Background: [rOpenSci AI policy](https://ropensci.org/blog/2026/02/26/ropensci-ai-policy/).

## Publication options

- [x] Do you intend for this package to go on CRAN?
- [ ] Do you intend for this package to go on Bioconductor?

I intend to submit to CRAN. I am also balancing a later relocation to Bioconductor: the companion chapter [From single cells to type-level Gaussians](https://github.com/bastienchassagnol/DeCovarT_reproducibility/blob/main/docs/06-population-alignment.qmd) discusses parameterisations that are quite specific to single-cell and bulk RNA-seq technical modalities (multiplicative rather than additive noise, library-depth normalisation, and related assay structure). Those pipelines live in the companion repository today; a Bioconductor home would make more sense if that assay-specific layer is later folded into the package.

## Code of conduct

- [x] I agree to abide by [rOpenSci's Code of Conduct](https://ropensci.org/code-of-conduct/) during the review process and in maintaining my package should it be accepted.
