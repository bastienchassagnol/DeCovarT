# Notation

A symbol absent from these tables is defined in the vignette that uses
it, with one meaning there. An article embeds only the groups it uses.
\boldsymbol{\theta} is not a relative-abundance symbol: it is the
parameter of interest in a profile likelihood. Other methods may still
write a per-cell mean as S_j\boldsymbol{\theta}\_{\cdot j}; that is
their split, quoted when a method is reviewed.

## Indices

| Symbol | Meaning |
|----|----|
| Indices | Indices |
| \$g = 1, \ldots, G\$ | Gene. |
| \$j = 1, \ldots, J\$; \$\ell\$ | Cell type, and a second cell type. |
| \$i = 1, \ldots, N\$ | Bulk sample. The index is dropped inside a derivation and restored in the model statements. |
| \$r = 1, \ldots, R\$ | Reference donor (replicate). |
| \$c\$; \$C_j\$ | Cell, and the number of reference cells of type \$j\$. |

Indices. The scalar \$R\$ is the number of reference donors.
{#tinytable_0hg7hskt86ksms3p27iw .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Cellular ratios

| Symbol | Meaning |
|----|----|
| Cell fractions | Cell fractions |
| \$m\_{ji}\$, \$n_i = \sum_j m\_{ji}\$ | Cells of type \$j\$ in bulk sample \$i\$, and their total. |
| \$\boldsymbol{p}\_{\cdot i} \in \Delta^{J-1}\$ | Cell fractions, \$p\_{ji} = m\_{ji} / n_i\$. |
| \$\boldsymbol{\rho}\_i \in \mathbb{R}^{J-1}\$^(a) | Isometric log-ratio coordinates of \$\boldsymbol{p}\_{\cdot i}\$, \$\boldsymbol{p}\_{\cdot i} = \operatorname{softmax}(\mathbf{V}\boldsymbol{\rho}\_i)\$ with a Helmert basis \$\mathbf{V}\$. |
| Transcript content | Transcript content |
| \$\boldsymbol{q}\_{\cdot i}\$^(b) | Transcript fractions, \$q\_{ji} = p\_{ji} S_j / \sum\_\ell p\_{\ell i} S\_\ell\$. |
| ^(a) An isometric log-ratio coordinate of the cell fractions. Not a correlation. |  |
| ^(b) Content weight used when a method normalises by transcriptome size. Not the cell-fraction estimand, which remains \$\boldsymbol{p}\$. Dividing counts by library size is not an observation model: genes that share a divisor become correlated. |  |

Cell fractions and transcript fractions.
{#tinytable_rx1v0xhz1qlr5gz1sh45 .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Reference means

| Symbol | Meaning |
|----|----|
| Reference means | Reference means |
| \$\boldsymbol{x}\_{\cdot j,i}\$ | Purified profile of type \$j\$ in sample \$i\$. |
| \$x\_{gc}\$, \$\boldsymbol{x}\_{\cdot c}\$ | Raw UMI count of gene \$g\$ in reference cell \$c\$, and that cell's count vector. |
| \$\tilde{\boldsymbol{\mu}}\_{\cdot j} = C_j^{-1} \sum\_{c \in j} \boldsymbol{x}\_{\cdot c}\$^(a) | Per-cell mean of the raw UMI counts of the \$C_j\$ reference cells of type \$j\$. |
| \$S_j = \mathbf{1}^{\top} \tilde{\boldsymbol{\mu}}\_{\cdot j}\$^(b) | Mean transcriptome size of a type-\$j\$ cell. |
| ^(a) A tilde marks a quantity rescaled to a reference unit: per cell, or relative to the reference protocol. This mean is a mean of raw counts. It is not a library-size proportion. |  |
| ^(b) \$S_j\$ summarises the raw-count mean. It is not a divisor applied to the counts. \$\boldsymbol{\theta}\$ is not used for this split. |  |

Purified profiles and per-cell means of raw counts.
{#tinytable_9yw0rpka2hk70mh4fb39 .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Biological variability

| Symbol | Meaning |
|----|----|
| Within one reference | Within one reference |
| \$\boldsymbol{\mu}\_{\cdot j}\$, \$\boldsymbol{\Sigma}\_j\$^(a) | Mean and \$G \times G\$ gene–gene covariance of a type-\$j\$ profile. |
| \$\boldsymbol{\Omega}\_j = \boldsymbol{\Sigma}\_j^{-1}\$ | Precision matrix of type \$j\$. |
| \$\boldsymbol{\Sigma}\_{W,j}\$ | Within-donor cell-to-cell covariance of type \$j\$ (one single-cell sample). |
| \$\boldsymbol{V}\_j\$ | Empirical covariance of type \$j\$, from its reference profiles or cells. |
| \$\boldsymbol{R}\_j\$ | Empirical correlation matrix of type \$j\$. Bold: the scalar \$R\$ is the number of donors. |
| \$\boldsymbol{\Delta}\_j = \operatorname{diag}(\boldsymbol{V}\_j)\$ | Diagonal matrix of the empirical variances of type \$j\$. |
| \$\lambda\$^(b) | Penalty of a regularised covariance or precision estimator. |
| Across donors or types | Across donors or types |
| \$\boldsymbol{\Sigma}\_{B,j}\$ | Between-donor covariance of type \$j\$ (needs \$R \ge 2\$ donors). |
| \$\boldsymbol{\Gamma}\_g\$^(c) | \$J \times J\$ covariance across cell types for gene \$g\$. |
| \$\boldsymbol{S}\_W(\boldsymbol{v}) = \sum_j v_j \boldsymbol{\Sigma}\_{W,j}\$ | Linear covariance mixture with weights \$\boldsymbol{v}\$. |
| \$\boldsymbol{S}\_B(\boldsymbol{v}) = \sum_j v_j^2 \boldsymbol{\Sigma}\_{B,j}\$ | Quadratic covariance mixture with weights \$\boldsymbol{v}\$. |
| ^(a) \$\boldsymbol{\Sigma}\_j\$ couples genes within one cell type. |  |
| ^(b) A penalty. The size-free expression of a cell is \$\boldsymbol{\lambda}\_{\cdot c}\$, in the technical-scale table. |  |
| ^(c) \$\boldsymbol{\Gamma}\_g\$ couples cell types within one gene. The opposite orientation from \$\boldsymbol{\Sigma}\_j\$. |  |

Gene–gene covariance within a type, and variation across donors or
types. {#tinytable_rdln4unpuub1yylw523o .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Technical scales

| Symbol | Meaning |
|----|----|
| Technical scales | Technical scales |
| \$s_c\$ | Size factor of reference cell \$c\$. |
| \$\boldsymbol{\lambda}\_{\cdot c}\$ | Size-free expression of cell \$c\$. |
| \$\boldsymbol{z}\_{\cdot c}\$ | Latent log-expression of cell \$c\$ in the Poisson log-normal model. |
| \$a_i\$ | Sample-wide scale between the reference and the bulk library. |
| \$c_i = a_i n_i\$, \$\phi_i = a_i^2 n_i\$ | Mean scale and variance scale of the single-cell model. |
| \$b_g\$ | Capture efficiency of gene \$g\$ in one protocol. |
| \$\tilde{b}\_g = b_g^{\mathrm{bulk}} / b_g^{\mathrm{ref}}\$^(a) | Gene-wise platform ratio. \$\tilde{\boldsymbol{B}} = \operatorname{diag}(\tilde{\boldsymbol{b}})\$. |
| \$\boldsymbol{D}\_i\$^(b) | Diagonal, heteroscedastic bulk measurement noise. |
| ^(a) A tilde: relative to the reference protocol, not per cell. |  |
| ^(b) Measurement noise. Not the deviance \$\mathcal{D}\_g\$ of a count model. |  |

Size factors, capture efficiencies, and measurement noise.
{#tinytable_ej182c18nqy4b0ou9mf5 .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Bulk observation

| Symbol | Meaning |
|----|----|
| Observation | Observation |
| \$y\_{gi}\$, \$\boldsymbol{y}\_{\cdot i}\$ | Bulk expression on a linear scale. |
| \$\tilde{\boldsymbol{y}}\_{\cdot i} = \boldsymbol{y}\_{\cdot i} / n_i\$^(a) | Bulk expression per cell. |
| \$\boldsymbol{T}\_{ji} = \sum\_{c=1}^{m\_{ji}} \boldsymbol{X}\_{jic}\$ | Total expression contributed by type \$j\$ to sample \$i\$. |
| Single-cell signatures | Single-cell signatures |
| \$\sum_j m_j \boldsymbol{\Sigma}\_j\$^(b) | Covariance when counts are fixed and cells of a type are independent. |
| Purified profiles | Purified profiles |
| \$\sum_j p_j^2 \boldsymbol{\Sigma}\_j\$^(c) | Covariance when each type contributes one population profile, as in a FACS or cytometry reference. |
| ^(a) Per cell, by the number of cells. Not a division of each gene by the library total: that division induces a spurious correlation between genes. |  |
| ^(b) Single-cell signature paradigm: the reference is the cells. |  |
| ^(c) Purified-profile paradigm: one profile per type. The likelihood in the package fits this covariance. It is not a preferred reference over the single-cell signatures. |  |

The bulk profile and the two reference paradigms. Neither is a
gold-standard reference. {#tinytable_frblsvm8de5n2472ag4y .table
.tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Feature selection

| Symbol | Meaning |
|----|----|
| Gene screens | Gene screens |
| \$\psi_g\$ | Negative binomial dispersion of gene \$g\$. |
| \$\mathcal{D}\_g\$ | Deviance of gene \$g\$ in a count model. |
| \$\delta\_{gj}\$ | Detection rate of gene \$g\$ in type \$j\$: the fraction of cells with a non-zero count. |

Dispersion, deviance, and detection rate.
{#tinytable_r9elketpv93r2xc1jub4 .table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}

## Overlap

| Symbol | Meaning |
|----|----|
| Overlap | Overlap |
| \$\varpi\_{\ell \mid j}\$, \$\varpi\_{j\ell}\$, \$\bar{\varpi}\$ | Probability of misclassifying a draw of type \$j\$ into type \$\ell\$; pairwise overlap \$\varpi\_{j\ell} = \varpi\_{\ell \mid j} + \varpi\_{j \mid \ell}\$; and its average over pairs. |

MixSim overlap between cell types. {#tinytable_oa6iteudl5l2mousuf5v
.table .tinytable
style="table-layout: fixed; width: 100% !important; margin-left: auto; margin-right: auto;"
quarto-disable-processing="true"}
