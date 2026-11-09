# 3.2 Cell-type covariance for the DeCovarT plug-in

> **Notation used in this vignette**
>
> A symbol absent from these groups is defined in this vignette, with
> one meaning. The full table is the [notation
> article](https://bastienchassagnol.github.io/DeCovarT/articles/notation.md).
>
> | Symbol | Meaning |
> |----|----|
> | **Indices** |  |
> | g = 1, \ldots, G | Gene. |
> | j = 1, \ldots, J; \ell | Cell type, and a second cell type. |
> | i = 1, \ldots, N | Bulk sample. The index is dropped inside a derivation and restored in the model statements. |
> | r = 1, \ldots, R | Reference donor (replicate). |
> | c; C_j | Cell, and the number of reference cells of type j. |
>
> Indices. The scalar R is the number of reference donors. {.table
> .caption-top}
>
> | Symbol | Meaning |
> |----|----|
> | **Reference means** |  |
> | \boldsymbol{x}\_{\cdot j,i} | Purified profile of type j in sample i. |
> | x\_{gc}, \boldsymbol{x}\_{\cdot c} | Raw UMI count of gene g in reference cell c, and that cell’s count vector. |
> | \tilde{\boldsymbol{\mu}}\_{\cdot j} = C_j^{-1} \sum\_{c \in j} \boldsymbol{x}\_{\cdot c}^(a) | Per-cell mean of the raw UMI counts of the C_j reference cells of type j. |
> | S_j = \mathbf{1}^{\top} \tilde{\boldsymbol{\mu}}\_{\cdot j}^(b) | Mean transcriptome size of a type-j cell. |
> | ^(a) A tilde marks a quantity rescaled to a reference unit: per cell, or relative to the reference protocol. This mean is a mean of raw counts. It is not a library-size proportion. |  |
> | ^(b)S_j summarises the raw-count mean. It is not a divisor applied to the counts. \boldsymbol{\theta} is not used for this split. |  |
>
> Purified profiles and per-cell means of raw counts. {.table
> .caption-top}
>
> | Symbol | Meaning |
> |----|----|
> | **Within one reference** |  |
> | \boldsymbol{\mu}\_{\cdot j}, \boldsymbol{\Sigma}\_j^(a) | Mean and G \times G gene–gene covariance of a type-j profile. |
> | \boldsymbol{\Omega}\_j = \boldsymbol{\Sigma}\_j^{-1} | Precision matrix of type j. |
> | \boldsymbol{\Sigma}\_{W,j} | Within-donor cell-to-cell covariance of type j (one single-cell sample). |
> | \boldsymbol{V}\_j | Empirical covariance of type j, from its reference profiles or cells. |
> | \boldsymbol{R}\_j | Empirical correlation matrix of type j. Bold: the scalar R is the number of donors. |
> | \boldsymbol{\Delta}\_j = \operatorname{diag}(\boldsymbol{V}\_j) | Diagonal matrix of the empirical variances of type j. |
> | \lambda^(b) | Penalty of a regularised covariance or precision estimator. |
> | **Across donors or types** |  |
> | \boldsymbol{\Sigma}\_{B,j} | Between-donor covariance of type j (needs R \ge 2 donors). |
> | \boldsymbol{\Gamma}\_g^(c) | J \times J covariance across cell types for gene g. |
> | \boldsymbol{S}\_W(\boldsymbol{v}) = \sum_j v_j \boldsymbol{\Sigma}\_{W,j} | Linear covariance mixture with weights \boldsymbol{v}. |
> | \boldsymbol{S}\_B(\boldsymbol{v}) = \sum_j v_j^2 \boldsymbol{\Sigma}\_{B,j} | Quadratic covariance mixture with weights \boldsymbol{v}. |
> | ^(a)\boldsymbol{\Sigma}\_j couples genes within one cell type. |  |
> | ^(b) A penalty. The size-free expression of a cell is \boldsymbol{\lambda}\_{\cdot c}, in the technical-scale table. |  |
> | ^(c)\boldsymbol{\Gamma}\_g couples cell types within one gene. The opposite orientation from \boldsymbol{\Sigma}\_j. |  |
>
> Gene–gene covariance within a type, and variation across donors or
> types. {.table .caption-top}
>
> | Symbol | Meaning |
> |----|----|
> | **Technical scales** |  |
> | s_c | Size factor of reference cell c. |
> | \boldsymbol{\lambda}\_{\cdot c} | Size-free expression of cell c. |
> | \boldsymbol{z}\_{\cdot c} | Latent log-expression of cell c in the Poisson log-normal model. |
> | a_i | Sample-wide scale between the reference and the bulk library. |
> | c_i = a_i n_i, \phi_i = a_i^2 n_i | Mean scale and variance scale of the single-cell model. |
> | b_g | Capture efficiency of gene g in one protocol. |
> | \tilde{b}\_g = b_g^{\mathrm{bulk}} / b_g^{\mathrm{ref}}^(a) | Gene-wise platform ratio. \tilde{\boldsymbol{B}} = \operatorname{diag}(\tilde{\boldsymbol{b}}). |
> | \boldsymbol{D}\_i^(b) | Diagonal, heteroscedastic bulk measurement noise. |
> | ^(a) A tilde: relative to the reference protocol, not per cell. |  |
> | ^(b) Measurement noise. Not the deviance \mathcal{D}\_g of a count model. |  |
>
> Size factors, capture efficiencies, and measurement noise. {.table
> .caption-top}

> **Scope**
>
> This vignette states which covariance DeCovarT needs from a
> single-cell reference, and how to estimate it. The target is the
> within-type covariance \boldsymbol{\Sigma}\_{W,j} that enters the
> likelihood, not a biologically interpretable gene regulatory network.
> The default estimator is linear shrinkage. Three sparse-precision
> alternatives are derived in [3.3 Precision
> networks](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.md)
> and scored by the same cross-validated criterion. The application to
> the Suppinger gastruloid reference (solver settings, outputs, plots)
> is in the [covariance chapter of the reproducibility
> book](https://github.com/bastienchassagnol/DeCovarT_reproducibility/blob/main/docs/03-covariance-estimation.qmd).
>
> Local symbols, defined once here: \omega\_{gg'} is an entry of
> \boldsymbol{\Omega}\_j; d is the maximum degree of a graph; K is the
> number of cross-validation folds; \alpha is an FDR level; \gamma is
> the EBIC parameter; \circ is the elementwise product and
> e^{\boldsymbol{A}} an elementwise exponential.

## The target: a covariance, not a graph

DeCovarT evaluates a Gaussian likelihood whose covariance is a sum of
cell-type covariances. With one single-cell specimen, the identified
object is the within-type cell-to-cell covariance
\boldsymbol{\Sigma}\_{W,j}, which enters linearly in p\_{ji}
([one-replicate
regime](https://bastienchassagnol.github.io/DeCovarT/articles/02-1-single-cell-generative-optimisation.html#sec-one-replicate)).
With purified population profiles, \boldsymbol{\Sigma}\_j enters through
\sum_jp\_{ji}^{2}\boldsymbol{\Sigma}\_j
([convolution](https://bastienchassagnol.github.io/DeCovarT/articles/05-1-review-distribution-aware.html#eq-gaussian-convolution)).
Either way the estimate must live on the linear count scale of the mean
signature
([normalisation](https://bastienchassagnol.github.io/DeCovarT/articles/02-1-single-cell-generative-optimisation.html#sec-norm-linear)),
not on \log(1+x) or on nonparanormal scores.

Two properties make the covariance, not the precision, the natural
object.

- **Marginalisation.** If the cells of type j have covariance
  \boldsymbol{\Sigma}\_j over all genes, the covariance of a panel
  \mathcal{P} is the sub-matrix
  (\boldsymbol{\Sigma}\_j)\_{\mathcal{P}\mathcal{P}}, whatever other
  genes are measured. The precision of the panel is
  \\(\boldsymbol{\Sigma}\_j)\_{\mathcal{P}\mathcal{P}}\\^{-1}, a Schur
  complement of \boldsymbol{\Omega}\_j, not the sub-matrix
  (\boldsymbol{\Omega}\_j)\_{\mathcal{P}\mathcal{P}}. Adding or removing
  a gene changes every partial correlation and every conditional
  variance 1/\omega\_{gg}; it leaves each covariance entry unchanged. A
  covariance entry is a property of a gene pair; a graph is a property
  of the panel.
- **The mixture has no sparse precision.** Sums such as
  \sum_jp\_{ji}\boldsymbol{\Sigma}\_{W,j} or
  \sum_jp\_{ji}^{2}\boldsymbol{\Sigma}\_j have dense inverses even when
  every \boldsymbol{\Omega}\_j is sparse. The likelihood never uses the
  zeros of \boldsymbol{\Omega}\_j.

Sparsity of \boldsymbol{\Omega}\_j therefore buys DeCovarT only a lower
variance, through \|E_j\|+G free parameters instead of G(G+1)/2.
Existence for C_j\le G does not require it: any regulariser provides
that. Interpretation of the zeros and the recovery theory of
[Sec. 2](#sec-cov-regimes) matter for a network, not for the
deconvolution. The default plug-in is therefore a **dense shrinkage
covariance** ([Sec. 3.1](#sec-cov-shrinkage)). Sparse precisions are
alternatives, judged on the same held-out score ([Sec. 4](#sec-cov-cv)).
When the proportions are known, `CSNet` estimates type-level
co-expression from the same convolution
([CSNet](https://bastienchassagnol.github.io/DeCovarT/articles/05-1-review-distribution-aware.html#sec-csnet))
([Su et al. 2024](#ref-suEstimatingCellTypeSpecificGene2024)).

| Estimator | Data it sees | Tuning parameter | Object exported to DeCovarT | Graph |
|----|----|----|----|----|
| Shrinkage, target D ([Sec. 3.1](#sec-cov-shrinkage)) | linear counts | \lambda\in\[0,1\]; closed-form \hat\lambda^{\ast} as start | dense \widehat{\boldsymbol{\Sigma}}\_j, variances kept | none |
| Graphical lasso, unpenalised diagonal ([graphical lasso](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-cov-glasso)) | correlation \boldsymbol{R}\_j of linear counts | \lambda | \boldsymbol{\Delta}\_j^{1/2}\widehat{\boldsymbol{\Sigma}}\_{R}\boldsymbol{\Delta}\_j^{1/2} | support of \widehat{\boldsymbol{\Omega}}\_j |
| `SILGGM` after the nonparanormal map ([nonparanormal](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-cov-npn)) | nonparanormal scores | FDR level \alpha | linear-scale refit under the zeros of the support | FDR-controlled support |
| `PLNnetwork` ([`PLNnetwork`](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-cov-pln)) | raw UMIs, offset \log s_c | \lambda | linear moments of the latent (\boldsymbol{\mu}\_z,\boldsymbol{\Sigma}\_z) ([linear-moment map](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#eq-pln-linear)) | support of the latent precision |

Table 1: Four covariance estimators. Each has one tuning parameter,
chosen by the criterion of [Sec. 4](#sec-cov-cv).

## Sample-size regimes

Report C_j against G for every type (and time point) before reading any
estimate. In the Suppinger reference several types have C_j\<G, so the
sample covariance \boldsymbol{V}\_j has rank at most C_j-1 and no
unpenalised estimator exists for them. Three thresholds are often
conflated.

1.  **The penalised estimate exists.** With an off-diagonal penalty, the
    graphical lasso has a unique positive definite solution for any
    \lambda\>0 once every V\_{gg}\>0, even when C_j\<G ([Carter
    2025](#ref-carterExistencePenalisedLikelihood2025)). Shrinkage
    towards a positive definite target exists for any C_j\ge2. No
    sparsity of the truth is needed.

2.  **The maximum-likelihood estimate on a fixed graph exists.**
    Dempster’s refit under the zeros of a graph E,

    \widehat{\boldsymbol{\Omega}}\_{E}
    =\arg\max\_{\boldsymbol{\Omega}\succ0}
    \log\det\boldsymbol{\Omega}-\operatorname{tr}(\boldsymbol{V}\_j\boldsymbol{\Omega})
    \quad\text{s.t.}\quad \omega\_{gg'}=0\\ \\ \forall\\(g,g')\notin E,
    \tag{1}

    exists for almost every sample once C_j reaches the
    maximum-likelihood threshold \operatorname{mlt}(E). That threshold
    lies between the clique number and the treewidth plus one, and
    equals the clique number for chordal graphs (Uhler 2012; Bernstein &
    Outlaw 2023). TODO-cite: Uhler 2012, doi:10.1214/11-AOS957
    TODO-cite: Bernstein & Outlaw 2023, arXiv:2312.03145 Bernstein and
    Outlaw simulate graphs selected by the graphical lasso from
    independent Gaussian data (G\le9): the probability that
    C_j\ge\operatorname{mlt} of the selected graph rises with \lambda
    and with C_j; with a single observation the graphical lasso never
    returned the empty graph, the only graph whose estimate then exists
    (pp. 4–5).

3.  **The graph is recovered.** Support recovery by the graphical lasso
    needs C_j\gtrsim d^{2}\log G, an irrepresentability condition,
    bounded eigenvalues and a minimum \|\omega\_{gg'}\| (Ravikumar et
    al. 2011). TODO-cite: Ravikumar, Wainwright, Raskutti & Yu 2011,
    doi:10.1214/11-EJS631 Neighbourhood selection needs d\log G
    ([Chiquet 2015, 41](#ref-chiquetContributionsSparseMethods2015)),
    under its own neighbourhood stability condition (Meinshausen &
    Bühlmann 2006, assumption 6). TODO-cite: Meinshausen & Bühlmann
    2006, doi:10.1214/009053606000000281 With G=500 and d=10, d^{2}\log
    G\approx620 before constants, and hubs raise d. On models derived
    from gene expression the irrepresentability conditions never held
    ([Heinävaara et al.
    2016](#ref-hein%C3%A4vaaraInconsistency%E2%84%931penalisedSparse2016)),
    and they involve the unknown \boldsymbol{\Sigma}\_j: “we can never
    verify them on real data” ([Chiquet 2015,
    77](#ref-chiquetContributionsSparseMethods2015)).

The consequence for DeCovarT: threshold 1 is all the plug-in needs.
Threshold 2 governs the refits used by `SILGGM`
([nonparanormal](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-cov-npn)).
Threshold 3 governs any biological reading of the zeros, which is out of
scope here; it can be probed by simulating precisions at the observed
C_j with
[`simulate_hierarchical_grn_moments()`](https://bastienchassagnol.github.io/DeCovarT/reference/simulate_hierarchical_grn_moments.md)
and measuring recovery.

## Estimators

### Linear shrinkage (default)

Schäfer and Strimmer shrink the sample covariance towards a target
\boldsymbol{T},
\widehat{\boldsymbol{\Sigma}}=\lambda\boldsymbol{T}+(1-\lambda)\boldsymbol{V}\_j
([Schäfer and Strimmer
2005](#ref-sch%C3%A4ferShrinkageApproachLargescale2005), eq. 3). Their
target D (“diagonal, unequal variance”) keeps the variances and shrinks
the correlations towards zero (Table 1, eq. 10):

\hat\Sigma\_{gg}=V\_{gg}, \qquad
\hat\Sigma\_{gg'}=(1-\lambda)\\R\_{j,gg'}\sqrt{V\_{gg}V\_{g'g'}}, \qquad
\hat\lambda^{\ast} =\frac{\sum\_{g\neq
g'}\widehat{\operatorname{Var}}(R\_{j,gg'})} {\sum\_{g\neq
g'}R\_{j,gg'}^{2}}, \tag{2}

with \hat\lambda^{\ast} truncated to \[0,1\]. The closed form minimises
the Frobenius risk (Ledoit–Wolf lemma; Ledoit & Wolf 2004). TODO-cite:
Ledoit & Wolf 2004, doi:10.1016/S0047-259X(03)00096-4 The estimate is
positive definite for any C_j and G, invariant to the scale of each
gene, and costs about two sample covariances. Its precision is dense;
partial correlations and a local-FDR graph can be read from it if wanted
(eqs. 11–15).

The Frobenius risk is not the Gaussian likelihood that DeCovarT
evaluates. \hat\lambda^{\ast} is therefore only the starting point of
the cross-validated search of [Sec. 4](#sec-cov-cv). In
`corpcor::cov.shrink()`, `lambda.var = 0` reproduces target D; the
default also shrinks the variances towards their median.

Sparse precisions are not the default plug-in. The graphical lasso, the
nonparanormal `SILGGM` estimator and `PLNnetwork` are derived in [3.3
Precision
networks](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.md).
Each still exports a covariance on the linear scale, and each is scored
by the same held-out criterion as shrinkage ([Sec. 4](#sec-cov-cv)).

## One criterion for all four: cross-validated Gaussian log-score

Each estimator of [Table 1](#tbl-cov-estimators) has one tuning
parameter: \lambda for shrinkage (searched around \hat\lambda^{\ast}),
\lambda for the graphical lasso, \alpha for `SILGGM`, \lambda for
`PLNnetwork`. All four are tuned by the same score.

**Score.** Split the C_j cells of type j into K folds. Fit
(\widehat{\boldsymbol{\mu}}^{(-k)},\widehat{\boldsymbol{\Sigma}}^{(-k)})
on the other folds, exported on the linear scale, and score the C\_{j,k}
held-out cells by their mean Gaussian log-density,

\ell^{(k)} =-\frac12\Bigl\[\log\det\widehat{\boldsymbol{\Sigma}}^{(-k)}
+\operatorname{tr}\bigl\\(\widehat{\boldsymbol{\Sigma}}^{(-k)})^{-1}\boldsymbol{V}\_{j,k}\bigr\\
+G\log2\pi\Bigr\], \tag{3}

with \boldsymbol{V}\_{j,k}=C\_{j,k}^{-1}\sum\_{c\in
k}(\boldsymbol{x}\_{\cdot
c}-\widehat{\boldsymbol{\mu}}^{(-k)})(\boldsymbol{x}\_{\cdot
c}-\widehat{\boldsymbol{\mu}}^{(-k)})^{\top}. The tuning value maximises
K^{-1}\sum_k\ell^{(k)} (or the one standard-error rule). Minus the
cross-validated score estimates the Kullback–Leibler risk of the fitted
law up to a constant ([Commenges and Jacqmin-Gadda 2015, sect.
2.4.2](#ref-commengesDynamicalBiostatisticalModels2015)). It applies to
maximum-likelihood, penalised and shrinkage fits alike, whereas AIC and
BIC are justified for the maximum-likelihood estimate. The folds split
cells of one specimen. The score therefore measures prediction of new
cells of that specimen, which is the population of
\boldsymbol{\Sigma}\_{W,j}, and says nothing about new specimens
([Zimmerman et al.
2021](#ref-zimmermanPracticalSolutionPseudoreplication2021)).

**Why a Gaussian score is legitimate for non-Gaussian cells.** For a
test cell with mean \boldsymbol{\mu}^{\ast} and covariance
\boldsymbol{\Sigma}^{\ast}, whatever its law,

\mathbb{E}\bigl\\-\log\phi_G(\boldsymbol{x};\boldsymbol{\mu},\boldsymbol{\Sigma})\bigr\\
=\tfrac12\bigl\\\log\det\boldsymbol{\Sigma}
+\operatorname{tr}(\boldsymbol{\Sigma}^{-1}\boldsymbol{\Sigma}^{\ast})
+(\boldsymbol{\mu}-\boldsymbol{\mu}^{\ast})^{\top}\boldsymbol{\Sigma}^{-1}(\boldsymbol{\mu}-\boldsymbol{\mu}^{\ast})
+G\log2\pi\bigr\\. \tag{4}

Only the first two moments enter. The excess over its value at
(\boldsymbol{\mu}^{\ast},\boldsymbol{\Sigma}^{\ast}) is, at
\boldsymbol{\mu}=\boldsymbol{\mu}^{\ast}, half of Stein’s loss
\operatorname{tr}(\boldsymbol{\Sigma}^{-1}\boldsymbol{\Sigma}^{\ast})-\log\det(\boldsymbol{\Sigma}^{-1}\boldsymbol{\Sigma}^{\ast})-G\ge0,
which is zero only at \boldsymbol{\Sigma}=\boldsymbol{\Sigma}^{\ast}.
The Gaussian log-score is thus a proper scoring rule for the mean and
covariance, which is all that DeCovarT’s bulk likelihood uses: the bulk
law is Gaussian by the central limit theorem over many cells, not
because single cells are Gaussian.

**Comparing estimators: nested cross-validation.** Tuning needs one
loop. Comparing the four tuned estimators needs two: an inner loop that
chooses each tuning parameter on the training part of an outer fold, and
an outer loop that scores the tuned estimators on cells the tuning never
saw. Reporting the inner-loop maximum as the estimator’s performance is
optimistic, because it is the maximum over a grid. Compare estimators by
their paired outer-fold differences, with the diagonal covariance
\boldsymbol{\Delta}\_j as the baseline.

**LOOCV is not nested cross-validation.** Leave-one-out is the choice
K=C_j in a single loop. It can serve as the inner or the outer loop of a
nested scheme but is not one, and it costs C_j fits per grid value.

**Fast approximation.** The universal approximate cross-validation
criterion approximates leave-one-out from one fit,

\mathrm{UACV}=\Psi(\hat{\boldsymbol{\theta}})+\operatorname{tr}\bigl(\boldsymbol{H}\_{\Phi}^{-1}\boldsymbol{U}\bigr),
\tag{5}

where \Psi is the mean assessment loss (here minus the log-score),
\boldsymbol{H}\_{\Phi} the Hessian of the estimating function and
\boldsymbol{U} (their \boldsymbol{K}) the mean cross-product of the
per-cell assessment and estimating scores; the error is O_p(C_j^{-2})
(Commenges et al. 2012, eqs. 6–7). TODO-cite: Commenges, Proust-Lima,
Samieri & Liquet 2012, arXiv:1206.1753 For a maximum-likelihood fit
assessed by its own log-likelihood it is the Takeuchi information
criterion divided by 2C_j (p. 7). A smooth penalty enters
\boldsymbol{H}\_{\Phi} and lowers the correction (p. 8). An \ell_1
penalty is not differentiable at zero, so for the graphical lasso and
`PLNnetwork` the formula is applied on the active set, with the support
held fixed; that restriction is our reading, not the paper’s.

**Diagnostics, not selectors.** EBIC ([Foygel and Drton
2010](#ref-foygelExtendedBayesianInformation2010)) charges \log
C_j+4\gamma\log G per edge and relies on the maximum-likelihood degrees
of freedom; at small C_j it returns empty graphs ([screening
threshold](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-grn-screening)).
StARS ([Liu et al. 2010](#ref-liuStabilityApproachRegularization2010))
chooses the least regularised graph whose edge instability across
subsamples stays below a threshold. It costs B\times L fits for B
subsamples and L penalties, targets the support rather than the
covariance, and tends to select too many edges in Gaussian graphical
models. Thresholded adaptive validation targets the support as well
([Laszkiewicz et al.
2021](#ref-laszkiewiczThresholdedAdaptiveValidation2021)). All three are
reported beside the cross-validated choice, not instead of it.

## The diagonal: conditional variance, not self-regulation

In a Gaussian graphical model,
\omega\_{gg}=1/\operatorname{Var}(x_g\mid\boldsymbol{x}\_{-g}): the part
of gene g’s variance that the other panel genes do not explain. On
counts it mixes transcriptional bursting and auto-regulation, regulators
outside the panel, and Poisson sampling, and \boldsymbol{V}\_j alone
cannot separate them; at droplet depths sampling is the largest share
for most genes ([Sarkar and Stephens
2021](#ref-sarkarSeparatingMeasurementExpression2021)). By
marginalisation ([Sec. 1](#sec-cov-target)) \omega\_{gg} also changes
with the panel.

**Penalising the diagonal biases the variances.** With a penalised
diagonal the KKT conditions give \hat\Sigma\_{gg}=V\_{gg}+\lambda
([Friedman et al. 2008](#ref-friedmanSparseInverseCovariance2008);
[Mazumder and Hastie 2011](#ref-mazumderGraphicalLassoNew2011)). On the
correlation scale every implied correlation shrinks by about
1/(1+\lambda), in the same direction as Poisson attenuation ([count
attenuation](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-grn-attenuation)).
The penalty does not move the empty-graph threshold of [screening
threshold](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-grn-screening).
DeCovarT’s likelihood uses the variances directly, so this bias enters
the deconvolution.

**Self-regulation belongs to directed dynamic models.** CardamomOT
estimates a directed self-interaction \theta\_{ii} for each gene and
gives it no specific penalty beyond its common \ell_1 term towards a
prior network; a positive auto-activation with a negative basal term is
“the most parsimonious way to produce bistability for an isolated gene”
(Maugé & Ventre 2026, p. 28, citing Ferrell 2002). TODO-cite: Maugé &
Ventre 2026, doi:10.1371/journal.pcbi.1014838 TODO-cite: Ferrell 2002,
doi:10.1016/S0955-0674(02)00314-9 (journal line on p. 1 disagrees with
the DOI; verify on Crossref) Such a parameter is linked to a precision
only through the stationary covariance of a linearised model. For
\mathrm{d}\boldsymbol{x}=\boldsymbol{A}(\boldsymbol{x}-\boldsymbol{\mu})\\\mathrm{d}t+\mathrm{d}\boldsymbol{w}
with noise covariance \boldsymbol{Q}, the stationary covariance solves
the Lyapunov equation
\boldsymbol{A}\boldsymbol{\Sigma}+\boldsymbol{\Sigma}\boldsymbol{A}^{\top}+\boldsymbol{Q}=\boldsymbol{0}.
For two genes with self-loops -a_1,-a_2, an edge 1\to2 of weight k and
\boldsymbol{Q}=\operatorname{diag}(q_1,q_2), writing
\kappa=(a_1+a_2)^{2}q_2+k^{2}q_1,

\omega\_{11}=\frac{2(a_1+a_2)\\a_1(a_1+a_2)q_2+k^{2}q_1\\}{q_1\\\kappa},
\qquad \omega\_{22}=\frac{2a_2(a_1+a_2)^{2}}{\kappa}, \qquad
\omega\_{12}=-\frac{2a_2(a_1+a_2)k}{\kappa}. \tag{6}

Each diagonal entry depends on both self-loops, the edge and both noise
levels; only for k=0 does it reduce to 2a_g/q_g, which still mixes the
self-loop with the noise. For a three-gene chain 1\to2\to3 the
stationary \omega\_{13} is non-zero for generic parameters although
genes 1 and 3 share no edge. A GGM diagonal is therefore neither a
self-loop nor its sign. Strong positive self-loops also produce
bistable, hence bimodal, within-type laws, which a Gaussian
\boldsymbol{\Sigma}\_{W,j} cannot represent; the DeCovarT counterpart is
to split the type into subtypes.

> **Decision**
>
> Never penalise the diagonal of \boldsymbol{\Omega}\_j: `glasso` with
> `penalize.diagonal = FALSE`, `PLNnetwork` with an unpenalised
> diagonal, shrinkage target D with `lambda.var = 0`. Schemes that
> re-estimate the diagonal after the off-diagonal step are compatible
> with this rule: SPACE estimates conditional variances from residuals
> ([Peng et al. 2009](#ref-pengPartialCorrelationEstimation2009)); the
> partial-correlation graphical lasso penalises only the partial
> correlations and is scale invariant, computed by alternating the two
> blocks ([Carter and Molinari
> 2025](#ref-carterExistenceOptimisationPartial2025)); the simplest
> version is a refit without penalty on a selected support
> ([Eq. 1](#eq-hard-cov-selection)).

A near-diagonal precision is not, by itself, evidence that genes are
independent. The five competing explanations, the screening threshold
and the count attenuation are in [3.3 Precision
networks](https://bastienchassagnol.github.io/DeCovarT/articles/03-3-precision-networks.html#sec-grn-diagonal).
Directed Gaussian networks and joint estimation across types are in the
same vignette. Graph scores (Hamming distance, structural Hamming
distance, BEELINE, stability selection used as a graph rule) do not
choose \widehat{\boldsymbol{\Sigma}}\_j. The plug-in is the covariance
with the best held-out Gaussian log-score ([Sec. 4](#sec-cov-cv)). Those
score definitions stay in the reproducibility file
`docs/05-grn-metrics.qmd`, which is not a listed chapter.

## References

Carter, Jack Storror. 2025. ‘Existence of Penalised Likelihood Estimates
and Posterior Propriety of Separable Prior Distributions for Gaussian
Precision Matrices’. In *arXiv.org*. Https://arxiv.org/abs/2505.20005v2.

Carter, Jack Storror, and Cesare Molinari. 2025. *Existence and
Optimisation of the Partial Correlation Graphical Lasso*.
arXiv:2510.25712. arXiv. <https://doi.org/10.48550/arxiv.2510.25712>.

Chiquet, Julien. 2015. ‘Contributions to Sparse Methods for Complex Data
Analysis’. PhD thesis, Université Paris-Saclay / AgroParisTech.

Commenges, Daniel, and Helene Jacqmin-Gadda. 2015. *Dynamical
Biostatistical Models*. Chapman; Hall/CRC.
<https://doi.org/10.1201/b19109>.

Foygel, Rina, and Mathias Drton. 2010. *Extended Bayesian Information
Criteria for Gaussian Graphical Models*. arXiv:1011.6640. arXiv.
<https://doi.org/10.48550/arxiv.1011.6640>.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. ‘Sparse
Inverse Covariance Estimation with the Graphical Lasso’. *Biostatistics
(Oxford, England)* 9. <https://doi.org/10.1093/biostatistics/kxm045>.

Heinävaara, Otte, Janne Leppä-aho, Jukka Corander, and Antti Honkela.
2016. ‘On the Inconsistency of \ell1-Penalised Sparse Precision Matrix
Estimation’. *BMC Bioinformatics* 17 (Suppl 16): 448.
<https://doi.org/10.1186/s12859-016-1309-x>.

Laszkiewicz, Mike, Asja Fischer, and Johannes Lederer. 2021.
‘Thresholded Adaptive Validation: Tuning the Graphical Lasso for Graph
Recovery’. *Proceedings of The 24th International Conference on
Artificial Intelligence and Statistics*, 1864–72.

Liu, Han, Kathryn Roeder, and Larry Wasserman. 2010. ‘Stability Approach
to Regularization Selection (StARS) for High Dimensional Graphical
Models’. *Advances in Neural Information Processing Systems* 24 (2):
1432–40.

Mazumder, Rahul, and Trevor Hastie. 2011. ‘The Graphical Lasso: New
Insights and Alternatives’. *Electronic Journal of Statistics* 6.
<https://doi.org/10.1214/12-ejs740>.

Peng, Jie, Pei Wang, Nengfeng Zhou, and Ji Zhu. 2009. ‘Partial
Correlation Estimation by Joint Sparse Regression Models’. *Journal of
the American Statistical Association* 104 (486): 735–46.
<https://doi.org/10.1198/jasa.2009.0126>.

Sarkar, Abhishek, and Matthew Stephens. 2021. ‘Separating Measurement
and Expression Models Clarifies Confusion in Single-Cell RNA Sequencing
Analysis’. *Nature Genetics* 53 (6): 770–77.
<https://doi.org/10.1038/s41588-021-00873-4>.

Schäfer, Juliane, and Korbinian Strimmer. 2005. ‘A Shrinkage Approach to
Large-Scale Covariance Matrix Estimation and Implications for Functional
Genomics’. *Statistical Applications in Genetics and Molecular Biology*
4: Article32. <https://doi.org/10.2202/1544-6115.1175>.

Su, Chang, Jingfei Zhang, and Hongyu Zhao. 2024. ‘Estimating
Cell-Type-Specific Gene Co-Expression Networks from Bulk Gene Expression
Data with an Application to Alzheimer’s Disease’. *Journal of the
American Statistical Association* 119 (546): 811–24.
<https://doi.org/10.1080/01621459.2023.2297467>.

Zimmerman, Kip D., Mark A. Espeland, and Carl D. Langefeld. 2021. ‘A
Practical Solution to Pseudoreplication Bias in Single-Cell Studies’.
*Nature Communications* 12 (1): 738.
<https://doi.org/10.1038/s41467-021-21038-1>.
