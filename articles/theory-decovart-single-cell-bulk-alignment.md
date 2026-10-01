# Aligning single-cell references with bulk: three generative formulations

> **Scope**
>
> Released DeCovarT fits the convolution \boldsymbol{y}\_{\cdot
> i}\mid\boldsymbol{p}\_{\cdot i}
> \sim\mathcal{N}\_G(\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i}, \sum_j
> p\_{ji}^{2}\boldsymbol{\Sigma}\_j) of the [derivatives
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-unconstrained).
> That law is exact when each reference column is one random
> **population profile**. A single-cell reference supplies a different
> object: many cells per type, a within-type cell-to-cell covariance,
> and a bulk sample that is a physical sum of cells measured on another
> platform. This vignette writes down three generative formulations that
> start from cells rather than from profiles, derives the gradient and
> Hessian of each in unconstrained coordinates, and places `MuSiC`,
> `RNA-Sieve`, `DWLS`, `Bisque` and `DECALS` on the same map. The
> simplex charts of the derivatives vignette (ILR, or ALR in its
> appendix) apply unchanged to the proportion block of every
> formulation, so they are stated once and not repeated.
>
> Notation follows the [perspectives
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-notation).
> Section cross-references to that vignette and to the derivatives
> vignette are given as links. Nothing here is implemented in
> [`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md);
> the equations are the specification for a later release.

> **Note 1: Notation added to the perspectives table**
>
> | Symbol | Meaning |
> |----|----|
> | m\_{ji}\in\mathbb{N} | Number of cells of type j in bulk sample i (`MuSiC` writes m\_{kj}) |
> | n_i=\sum_j m\_{ji} | Total number of cells in bulk sample i |
> | \boldsymbol{m}\_{\cdot i}=(m\_{1i},\ldots,m\_{Ji})^{\top} | Cell-count vector; \boldsymbol{p}\_{\cdot i}=\boldsymbol{m}\_{\cdot i}/n_i |
> | \kappa_i=1/n_i | Per-cell averaging factor (the effective cells-per-bulk scale of the perspectives vignette) |
> | \boldsymbol{w}\_{\cdot i}=a_i\boldsymbol{m}\_{\cdot i} | Effective counts on the bulk molecule scale |
> | \boldsymbol{D}\_i=\boldsymbol{\Sigma}\_{\mathrm{obs},i}=\operatorname{diag}(\tau\_{1i}^{2},\ldots,\tau\_{Gi}^{2}) | Diagonal, heteroscedastic bulk measurement noise |
> | \boldsymbol{T}\_{ji}=\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic} | Total expression contributed by type j to sample i |
> | \boldsymbol{S}\_W(\boldsymbol{v})=\sum_j v_j\boldsymbol{\Sigma}\_{W,j} | Linear covariance mixture with weights \boldsymbol{v} |
> | \boldsymbol{S}\_B(\boldsymbol{v})=\sum_j v_j^{2}\boldsymbol{\Sigma}\_{B,j} | Quadratic covariance mixture |
> | \varrho_j\in\[-1/(m\_{ji}-1),1\] | Cell-to-cell correlation inside type j ([Sec. 7.1](#sec-persp-equicorrelation)); written \varrho because \boldsymbol{\rho} is the ILR coordinate of the derivatives vignette |
> | \boldsymbol{\vartheta} | Generic parameter vector; \boldsymbol{\theta}\_{\cdot j} keeps its meaning of relative abundance |
> | \boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k,\boldsymbol{\mu}\_{kl},\boldsymbol{\Sigma}\_{kl} | First and second partial derivatives of the mean and covariance in \vartheta_k,\vartheta_l |
>
> Bulk sample i is fixed throughout, so the index i is dropped inside
> derivations and restored in the model statements.

``` mermaid
%%{init: {"theme": "sandstone", "flowchart": {"curve": "basis"}}}%%
flowchart TD
  ref["Single-cell reference: mu_j and Sigma_W,j from R = 1, Sigma_B,j needs R at least 2"] --> q{"What is the bulk observation?"}
  q -->|"absolute molecules"| S1["Strategy 1: cell-count convolution. Target m_i in N^J, no simplex"]
  q -->|"per-cell average"| S2["Strategy 2: averaged mixture. Target p_i on the simplex plus n_i"]
  q -->|"library-normalised"| S3["Strategy 3: compositional bulk. Target RNA fractions q_i"]
  S1 --> L11["1. linear: sum_j m_ji Sigma_W,j"]
  L11 --> L12["2. m_ji non-negative integers"]
  L12 --> L13["3. sample scale a_i, enters squared"]
  L13 --> L14["4. diagonal bulk noise D_i"]
  L14 --> L15["5. R replicates: add m_ji squared Sigma_B,j"]
  L15 --> out1["p_ji = m_ji / n_i by division"]
  S2 --> L21["1. kappa_i sum_j p_ji Sigma_W,j, with n_i from RNA-Sieve-type variance or Bisque-type pairing"]
  L21 --> L22["2. scale a_i, per sample or shared"]
  L22 --> L23["3. diagonal bulk noise D_i"]
  L23 --> L24["4. R replicates: add p_ji squared Sigma_B,j"]
  L24 --> out2["DeCovarT today = quadratic layer only"]
  S3 --> L31["closure Z_i = T_i / total, delta-method Gaussian on the tangent space"]
  L31 --> out3["q_i differs from p_i unless S_j constant"]
```

Figure 1: Three generative formulations of a single-cell-referenced
bulk, each with its layers of complexity. Strategy 1 keeps absolute cell
counts and needs no simplex. Strategy 2 averages per cell and recovers
DeCovarT’s quadratic layer when biological replicates exist. Strategy 3
closes the bulk to a composition and is kept as a note.

[Figure 1](#fig-three-strategies) is the map of this vignette.
[Sec. 1](#sec-one-formula) shows that the linear and quadratic
covariance weights are two endpoints of one identity.
[Sec. 2](#sec-general) gives a single score and Hessian theorem that
every later model specialises. [Sec. 3](#sec-cell-counts),
[Sec. 4](#sec-averaged) and [Sec. 5](#sec-compositional) are the three
strategies. [Sec. 6](#sec-summary) compares them, and
[Sec. 7](#sec-perspectives) lists extensions.

## Convolution and weighted sum are one formula

Let
\boldsymbol{X}\sim\mathcal{N}\_G(\boldsymbol{\mu},\boldsymbol{\Sigma}).
Two statements look contradictory,

\operatorname{Cov}(a\boldsymbol{X})=a^{2}\boldsymbol{\Sigma}
\qquad\text{versus}\qquad
\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)=a\boldsymbol{\Sigma},
\tag{1}

and both are correct. The difference lies entirely in the dependence
between the a summands, not in two rules for Gaussian vectors.

### Characteristic-function proof

The characteristic function of \boldsymbol{X} is
\phi\_{\boldsymbol{X}}(\boldsymbol{t})
=\exp(i\boldsymbol{t}^{\top}\boldsymbol{\mu}
-\tfrac12\boldsymbol{t}^{\top}\boldsymbol{\Sigma}\boldsymbol{t}).

**Scaling** substitutes a\boldsymbol{t} for \boldsymbol{t}:

\phi\_{a\boldsymbol{X}}(\boldsymbol{t})
=\phi\_{\boldsymbol{X}}(a\boldsymbol{t})
=\exp\bigl(i\boldsymbol{t}^{\top}(a\boldsymbol{\mu})
-\tfrac12\boldsymbol{t}^{\top}(a^{2}\boldsymbol{\Sigma})\boldsymbol{t}\bigr),
\qquad
a\boldsymbol{X}\sim\mathcal{N}\_G(a\boldsymbol{\mu},a^{2}\boldsymbol{\Sigma}).
\tag{2}

**Independent convolution** multiplies a identical factors:

\phi\_{\boldsymbol{X}\_1+\cdots+\boldsymbol{X}\_a}(\boldsymbol{t})
=\phi\_{\boldsymbol{X}}(\boldsymbol{t})^{a}
=\exp\bigl(i\boldsymbol{t}^{\top}(a\boldsymbol{\mu})
-\tfrac12\boldsymbol{t}^{\top}(a\boldsymbol{\Sigma})\boldsymbol{t}\bigr),
\qquad
\sum\_{r}\boldsymbol{X}\_r\sim\mathcal{N}\_G(a\boldsymbol{\mu},a\boldsymbol{\Sigma}).
\tag{3}

The first route squares the quadratic exponent; the second adds a copies
of it. The density convolution (f\*f)(\boldsymbol{s})=\int
f(\boldsymbol{x})f(\boldsymbol{s}-\boldsymbol{x})\\d\boldsymbol{x} is
the law of \boldsymbol{X}\_1+\boldsymbol{X}\_2 only when the two draws
are independent. Saying “a convolution of Gaussians” therefore already
asserts independent summands.

### The master identity and its two endpoints

For arbitrary random vectors with
\operatorname{Cov}(\boldsymbol{X}\_r)=\boldsymbol{\Sigma} and
cross-covariances
\boldsymbol{C}\_{rs}=\operatorname{Cov}(\boldsymbol{X}\_r,\boldsymbol{X}\_s),

\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)
=\sum\_{r,s}\operatorname{Cov}(\boldsymbol{X}\_r,\boldsymbol{X}\_s)
=a\boldsymbol{\Sigma}+\sum\_{r\neq s}\boldsymbol{C}\_{rs}. \tag{4}

Under **equicorrelation**,
\boldsymbol{C}\_{rs}=\varrho\boldsymbol{\Sigma} for r\neq s, the a(a-1)
off-diagonal blocks collapse to

\operatorname{Cov}\Bigl(\sum\_{r=1}^{a}\boldsymbol{X}\_r\Bigr)
=a\bigl\[1+(a-1)\varrho\bigr\]\boldsymbol{\Sigma}, \qquad
-\frac{1}{a-1}\le\varrho\le 1 . \tag{5}

Setting \varrho=1 (the same realised vector reused a times) gives
a^{2}\boldsymbol{\Sigma}; setting \varrho=0 (fresh independent draws)
gives a\boldsymbol{\Sigma}. [Table 1](#tbl-dependence-regimes) lists the
cases. A sum of merely *marginally* Gaussian vectors need not be
Gaussian; [Eq. 5](#eq-equicorrelation) is a Gaussian law only when the
stacked vector
(\boldsymbol{X}\_1^{\top},\ldots,\boldsymbol{X}\_a^{\top})^{\top} is
**jointly** Gaussian, which is the assumption made whenever an
intermediate \varrho is used below.

| Construction | Dependence | Mean | Covariance | Law |
|----|----|----|----|----|
| Same \boldsymbol{X} repeated | \varrho=1 | a\boldsymbol{\mu} | a^{2}\boldsymbol{\Sigma} | \mathcal{N}\_G |
| a i.i.d. copies | \varrho=0 | a\boldsymbol{\mu} | a\boldsymbol{\Sigma} | \mathcal{N}\_G |
| Correlated copies | general \boldsymbol{C}\_{rs} | a\boldsymbol{\mu} | a\boldsymbol{\Sigma}+\sum\_{r\neq s}\boldsymbol{C}\_{rs} | Gaussian if jointly Gaussian |
| Equicorrelated copies | \boldsymbol{C}\_{rs}=\varrho\boldsymbol{\Sigma} | a\boldsymbol{\mu} | a\[1+(a-1)\varrho\]\boldsymbol{\Sigma} | Gaussian if jointly Gaussian |

Table 1: Dependence regimes for a sum of a Gaussian vectors with common
marginal covariance. The first two rows are the endpoints \varrho=1 and
\varrho=0 of the last.

The biological reading is the one that matters for a single-cell
reference. If \boldsymbol{\Sigma}\_j is the covariance **between cells**
of type j, summing m_j independent cells yields
m_j\boldsymbol{\Sigma}\_j. If \boldsymbol{\Sigma}\_j is the covariance
of **one latent type-level profile** that is multiplied by a weight, the
contribution is p_j^{2}\boldsymbol{\Sigma}\_j. Released DeCovarT is the
second case. The random-effects decomposition of the [perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-three-regimes),

\boldsymbol{X}\_{jc}=\boldsymbol{M}\_j+\boldsymbol{\varepsilon}\_{jc},
\qquad
\operatorname{Cov}\Bigl(\sum\_{c=1}^{m_j}\boldsymbol{X}\_{jc}\Bigr)
=m_j\boldsymbol{\Sigma}\_{W,j}+m_j^{2}\boldsymbol{\Sigma}\_{B,j},
\tag{6}

contains both endpoints at once: independent cell fluctuations
accumulate linearly, a shared type-level fluctuation is amplified
quadratically. Equicorrelated cells are the special case
\boldsymbol{\Sigma}\_{B,j}=\varrho_j\boldsymbol{\Sigma}\_j and
\boldsymbol{\Sigma}\_{W,j}=(1-\varrho_j)\boldsymbol{\Sigma}\_j
([Sec. 7.1](#sec-persp-equicorrelation)).

## One score and Hessian for every formulation

All models below are Gaussian with a parametrised mean and covariance.
Write \boldsymbol{y}\mid\boldsymbol{\vartheta}\sim
\mathcal{N}\_G\bigl(\boldsymbol{\mu}(\boldsymbol{\vartheta}),
\boldsymbol{\Sigma}(\boldsymbol{\vartheta})\bigr),
\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}(\boldsymbol{\vartheta}),
\boldsymbol{\Omega}=\boldsymbol{\Sigma}^{-1}, and up to a constant

\ell(\boldsymbol{\vartheta})
=\tfrac12\log\det\boldsymbol{\Omega}(\boldsymbol{\vartheta})
-\tfrac12\boldsymbol{r}^{\top}\boldsymbol{\Omega}(\boldsymbol{\vartheta})\boldsymbol{r}.
\tag{7}

The identities of the derivatives vignette ([matrix
calculus](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-matrix-calculus))
give the following.

**Theorem 1 (Score, Hessian and expected information of a parametrised
Gaussian)** With
\boldsymbol{\mu}\_k=\partial\boldsymbol{\mu}/\partial\vartheta_k,
\boldsymbol{\Sigma}\_k=\partial\boldsymbol{\Sigma}/\partial\vartheta_k,
\boldsymbol{\mu}\_{kl}=\partial^{2}\boldsymbol{\mu}/\partial\vartheta_k\partial\vartheta_l
and
\boldsymbol{\Sigma}\_{kl}=\partial^{2}\boldsymbol{\Sigma}/\partial\vartheta_k\partial\vartheta_l,

\frac{\partial\ell}{\partial\vartheta_k}
=\underbrace{-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k)}\_{\text{determinant}}
+\underbrace{\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{r}}\_{\text{mean
residual}}
+\underbrace{\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}}\_{\text{covariance
quadratic}}, \tag{8}

\begin{aligned}
\frac{\partial^{2}\ell}{\partial\vartheta_k\partial\vartheta_l}
&=\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l)
-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl})
-\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_l
+\boldsymbol{\mu}\_{kl}^{\top}\boldsymbol{\Omega}\boldsymbol{r} \\
&\quad
-\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{\mu}\_l^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl}\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r},
\end{aligned} \tag{9}

I(\boldsymbol{\vartheta})\_{kl}
=\mathbb{E}\Bigl\[-\frac{\partial^{2}\ell}{\partial\vartheta_k\partial\vartheta_l}\Bigr\]
=\boldsymbol{\mu}\_k^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_l
+\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l).
\tag{10}

*Proof*. The determinant term follows from
\partial\log\det\boldsymbol{\Omega}/\partial\vartheta_k
=-\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k) and
\partial\boldsymbol{\Omega}/\partial\vartheta_k
=-\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}.
Differentiating
-\tfrac12\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{r} uses
\partial\boldsymbol{r}/\partial\vartheta_k=-\boldsymbol{\mu}\_k and the
symmetry of \boldsymbol{\Omega}. Differentiating
[Eq. 8](#eq-score-general) once more and using
\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{r}
=\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}
(transpose of a scalar) gives [Eq. 9](#eq-hess-general). For
[Eq. 10](#eq-fisher-general) substitute
\mathbb{E}\[\boldsymbol{r}\]=\boldsymbol{0} and
\mathbb{E}\[\boldsymbol{r}\boldsymbol{r}^{\top}\]=\boldsymbol{\Sigma}:
the \boldsymbol{\Sigma}\_{kl} traces cancel and
\mathbb{E}\[\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l\boldsymbol{\Omega}\boldsymbol{r}\]
=\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_k\boldsymbol{\Omega}\boldsymbol{\Sigma}\_l).

> **Tip 2: Check against released DeCovarT**
>
> Take \boldsymbol{\vartheta}=\boldsymbol{p},
> \boldsymbol{\mu}(\boldsymbol{p})=\boldsymbol{\mu}\boldsymbol{p} and
> \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
> p_j^{2}\boldsymbol{\Sigma}\_j. Then
> \boldsymbol{\mu}\_{p_k}=\boldsymbol{\mu}\_{\cdot k},
> \boldsymbol{\mu}\_{p_kp_l}=\boldsymbol{0}, the covariance derivative
> is \boldsymbol{\Sigma}\_{p_k}=2p_k\boldsymbol{\Sigma}\_k with
> \boldsymbol{\Sigma}\_k the type covariance,
> \boldsymbol{\Sigma}\_{p_kp_k}=2\boldsymbol{\Sigma}\_k and
> \boldsymbol{\Sigma}\_{p_kp_l}=\boldsymbol{0} for k\neq l. Substituting
> into [Eq. 8](#eq-score-general) to [Eq. 10](#eq-fisher-general)
> recovers the [unconstrained
> gradient](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-score-unconstrained),
> [Hessian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-hess-unconstrained)
> and [expected Fisher
> information](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-fisher-p)
> of the derivatives vignette term by term, including the factor 2p_kp_l
> in front of the covariance trace.

Every model below is therefore specified by its tuple
(\boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k,\boldsymbol{\mu}\_{kl},\boldsymbol{\Sigma}\_{kl}).
When a parameter block lives on the simplex, the ILR chain rule
([constrained score and
Hessian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#thm-chain-rule))
is applied to that block only; the Jacobian of the full map is block
diagonal with identity blocks for the scalar parameters a_i, \kappa_i
and \tau\_{gi}^{2}. Positivity of scalars is handled by a log chart
(a_i=e^{\eta_i}), which multiplies the corresponding row of the score by
a_i and adds the score itself to the diagonal of the Hessian.

## Strategy 1: cell-count convolution

The first formulation keeps the bulk on an absolute molecule scale and
asks how many cells of each type it contains. It is the cell-level
reading of `MuSiC`’s equation (2), X\_{jg}=\sum_k
m\_{kj}S\_{kj}\theta\_{kjg}, where m\_{kj} is the number of cells of
type k in subject j ([Wang et al. 2019](#ref-wangBulkTissueCell2019)).
With \boldsymbol{\mu}\_{\cdot j}=S_j\boldsymbol{\theta}\_{\cdot j} the
mean below is exactly that identity before `MuSiC` divides by the
library total. The estimand is \boldsymbol{m}\_{\cdot
i}\in\mathbb{N}^{J}, proportions follow by division, and no simplex
constraint is imposed.

### Layer 1: independent cells, linear covariance

Assume one single-cell replicate (R=1) and, conditional on the counts,
independent cells

\boldsymbol{X}\_{jic}\overset{\text{i.i.d.}}{\sim}
\mathcal{N}\_G(\boldsymbol{\mu}\_{\cdot j},\boldsymbol{\Sigma}\_{W,j}),
\qquad c=1,\ldots,m\_{ji}, \tag{11}

independent across types. The bulk is the physical sum

\boldsymbol{y}\_{\cdot i} =\sum\_{j=1}^{J}\boldsymbol{T}\_{ji}
=\sum\_{j=1}^{J}\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}. \tag{12}

Applying [Eq. 4](#eq-master-identity) with
\boldsymbol{C}\_{rs}=\boldsymbol{0} inside each type and independence
across types,

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i}
\sim\mathcal{N}\_G\Bigl( \boldsymbol{\mu}\boldsymbol{m}\_{\cdot i},\\
\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i}) \Bigr), \qquad
\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i})=\sum_j
m\_{ji}\boldsymbol{\Sigma}\_{W,j}. \tag{13}

Each term has one meaning. \boldsymbol{\mu}\boldsymbol{m}\_{\cdot i} is
the expected molecule count of a tissue made of m\_{ji} average cells of
each type. m\_{ji}\boldsymbol{\Sigma}\_{W,j} is the accumulated
cell-to-cell scatter of type j: independent deviations add, so the
weight is linear. There is no \boldsymbol{\Sigma}\_{B,j} because one
replicate cannot identify it. The law is exact under Gaussian cells.
Under non-Gaussian cells it is the central-limit approximation of
`RNA-Sieve`, applied here to the full vector rather than gene by gene
([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)).

Treat \boldsymbol{m}\_{\cdot i} as a point of (0,\infty)^{J} for the
moment. The tuple of [Theorem 1](#thm-score-general) is

\boldsymbol{\mu}\_k=\boldsymbol{\mu}\_{\cdot k}, \qquad
\boldsymbol{\Sigma}\_k=\boldsymbol{\Sigma}\_{W,k}, \qquad
\boldsymbol{\mu}\_{kl}=\boldsymbol{0}, \qquad
\boldsymbol{\Sigma}\_{kl}=\boldsymbol{0}. \tag{14}

**Theorem 2 (Score, Hessian and information for cell counts)** With
\boldsymbol{\Omega}=\boldsymbol{S}\_W(\boldsymbol{m})^{-1} and
\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{m},

\frac{\partial\ell}{\partial m_k}
=-\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k})
+\boldsymbol{\mu}\_{\cdot k}^{\top}\boldsymbol{\Omega}\boldsymbol{r}
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{r},
\tag{15}

\frac{\partial^{2}\ell}{\partial m_k\partial m_l}
=\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,l})
-\boldsymbol{\mu}\_{\cdot
k}^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_{\cdot l}
-\boldsymbol{\mu}\_{\cdot
k}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,l}\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{\mu}\_{\cdot
l}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{r}
-\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,l}\boldsymbol{\Omega}\boldsymbol{r},
\tag{16}

I(\boldsymbol{m})\_{kl} =\boldsymbol{\mu}\_{\cdot
k}^{\top}\boldsymbol{\Omega}\boldsymbol{\mu}\_{\cdot l}
+\tfrac12\\\mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,l}).
\tag{17}

Compared with released DeCovarT the Hessian loses two terms, because
\boldsymbol{\Sigma}(\boldsymbol{m}) is **linear** in the parameter:
there is no \mathrm{tr}(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{kl})
and no second residual quadratic. The covariance derivative is a
constant matrix, so one Cholesky factorisation of
\boldsymbol{S}\_W(\boldsymbol{m}) per iteration serves objective, score
and Hessian, as in
[`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-numerical-speedups).

Proportions and their uncertainty follow by the delta method. With \hat
n=\mathbf{1}^{\top}\hat{\boldsymbol{m}} and
\hat{\boldsymbol{p}}=\hat{\boldsymbol{m}}/\hat n,

\frac{\partial\boldsymbol{p}}{\partial\boldsymbol{m}^{\top}}
=\frac{1}{n}\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr),
\qquad \operatorname{Var}(\hat{\boldsymbol{p}}) \approx\frac{1}{n^{2}}
\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr)
I(\boldsymbol{m})^{-1}
\bigl(\mathbf{I}\_J-\boldsymbol{p}\mathbf{1}^{\top}\bigr)^{\top}.
\tag{18}

That covariance is singular along \mathbf{1}, as it must be for a
composition, without any log-ratio chart. A type with \hat m\_{ji}=0 is
an ordinary point of the parameter space, not a boundary of a logarithm;
the structural-zero machinery of the [derivatives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-structural-zeros)
is not needed here.

> **Important 3: Where the GRN precision goes**
>
> The reproducibility book estimates one sparse precision per type from
> single cells and exports \\\hat{\boldsymbol{\mu}}\_{\cdot
> j},\hat{\boldsymbol{\Sigma}}\_j,\hat{\boldsymbol{\Omega}}\_j\\ on the
> linear count scale ([GGM
> chapter](https://github.com/bastienchassagnol/DeCovarT_reproducibility/blob/main/docs/02-grn-graphical-lasso.qmd)).
> The rows of that Gaussian table are **cells** from one experiment, so
> the estimate is (\boldsymbol{\mu}\_{\cdot
> j},\boldsymbol{\Sigma}\_{W,j},\boldsymbol{\Omega}\_{W,j}), not a
> population-profile covariance. Its slot is the linear term
> m\_{ji}\boldsymbol{\Sigma}\_{W,j} of [Eq. 13](#eq-cc-law), and the
> per-cell average of C_j such cells has covariance
> \boldsymbol{\Sigma}\_{W,j}/C_j and precision
> C_j\boldsymbol{\Omega}\_{W,j}: the zero pattern survives averaging.
> The precision of the **sum** \boldsymbol{S}\_W(\boldsymbol{m}) is
> dense in general even when every \boldsymbol{\Omega}\_{W,j} is sparse.
> The structure-aware backends of the derivatives vignette act on
> \boldsymbol{\Sigma}, not on \boldsymbol{\Omega}, so a dense Cholesky
> is the default here unless the \boldsymbol{\Sigma}\_{W,j} share a
> block or low-rank structure.

### Layer 2: counts are non-negative integers

m\_{ji}\in\\0,1,2,\ldots\\, so the likelihood is a function on the
lattice \mathbb{N}^{J} and [Eq. 15](#eq-cc-score) is the score of its
continuous relaxation. Three facts make the relaxation usable.

1.  The Gaussian law [Eq. 13](#eq-cc-law) is well defined on the whole
    lattice: \boldsymbol{S}\_W(\boldsymbol{m}) is positive definite as
    soon as one m\_{ji}\>0 and every \boldsymbol{\Sigma}\_{W,j} is
    positive definite. Zero counts cost nothing.
2.  For abundant types the lattice spacing is one cell out of thousands;
    rounding the relaxed maximiser changes \log-likelihood by O(1/n_i)
    relative terms.
3.  For rare types (m\_{ji}\le 5) the covariance changes in visible
    steps, and the relaxed optimum can sit between integers. A local
    lattice search over the 3^{J'} neighbours \\\lfloor\hat
    m\rfloor-1,\lfloor\hat m\rfloor,\lceil\hat m\rceil\\ of the J' small
    coordinates is cheap because each candidate needs one Cholesky.
    Branch-and-bound on the concave relaxation is the exact alternative.

A soft version keeps the Gaussian observation law and adds a count
prior, m\_{ji}\sim\mathrm{Poisson}(n_i\pi\_{ji}) or a multinomial over
types, and maximises the posterior over the lattice. The observation
model stays in the Gaussian regime; only the parameter space becomes
discrete.

### Layer 3: a sample-wide capture scale

Bulk and reference are measured on different platforms. Let a_i\>0 be
the sample-wide ratio of molecule capture between the bulk library and
the reference cells, shared by all genes. The whole realised signal is
amplified,

\boldsymbol{y}\_{\cdot i}
=a_i\sum_j\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}, \qquad
\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl( a_i\boldsymbol{\mu}\boldsymbol{m}\_{\cdot i},\\
a_i^{2}\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot i}) \Bigr). \tag{19}

The two scalings of [Sec. 1](#sec-one-formula) appear side by side:
m\_{ji} enters linearly because m\_{ji} independent cells are convolved;
a_i enters squared because one realised fluctuation is amplified a_i
times. Reparametrise to effective counts \boldsymbol{w}\_{\cdot
i}=a_i\boldsymbol{m}\_{\cdot i}:

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{w}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl( \boldsymbol{\mu}\boldsymbol{w}\_{\cdot i},\\
a_i\boldsymbol{S}\_W(\boldsymbol{w}\_{\cdot i}) \Bigr). \tag{20}

The mean identifies \boldsymbol{w}\_{\cdot i} alone; the covariance
magnitude then identifies a_i, and \boldsymbol{m}\_{\cdot
i}=\boldsymbol{w}\_{\cdot i}/a_i. That is the route by which `RNA-Sieve`
infers its total cell number n: the variance of a sum of n cells grows
like n while the mean grows like n too, so their ratio pins the count
([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)). Its
supplement states the caveat that applies here as well: n is physically
meaningful within one protocol, or when relative amplification factors
are known, and loses that meaning when the cross-protocol scale is
unknown. In [Eq. 20](#eq-cc-w-law) the same caveat reads: a_i is
identified only if \boldsymbol{\Sigma}\_{W,j} is on the molecule scale
of the reference cells, which is where the GRN estimates of
[Note 3](#nte-sparse-precision-slot) live.

In the (\boldsymbol{w},a) chart the tuple is

\boldsymbol{\mu}\_{w_k}=\boldsymbol{\mu}\_{\cdot k},\quad
\boldsymbol{\mu}\_{a}=\boldsymbol{0},\quad
\boldsymbol{\Sigma}\_{w_k}=a\boldsymbol{\Sigma}\_{W,k},\quad
\boldsymbol{\Sigma}\_{a}=\boldsymbol{S}\_W(\boldsymbol{w}),\quad
\boldsymbol{\Sigma}\_{w_ka}=\boldsymbol{\Sigma}\_{W,k},\quad
\boldsymbol{\Sigma}\_{w_kw_l}=\boldsymbol{\Sigma}\_{aa}=\boldsymbol{0},
\tag{21}

and all \boldsymbol{\mu}\_{kl}=\boldsymbol{0}. The scale score is

\frac{\partial\ell}{\partial a}
=-\tfrac12\\\mathrm{tr}\bigl(\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})\bigr)
+\tfrac12\\\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})\boldsymbol{\Omega}\boldsymbol{r}
=\frac{1}{2a}\bigl(\boldsymbol{r}^{\top}\boldsymbol{\Omega}\boldsymbol{r}-G\bigr),
\tag{22}

because
\boldsymbol{\Omega}\boldsymbol{S}\_W(\boldsymbol{w})=a^{-1}\mathbf{I}\_G
when \boldsymbol{D}\_i=\boldsymbol{0}. Setting it to zero gives the
closed form \hat
a=\boldsymbol{r}^{\top}\boldsymbol{S}\_W(\boldsymbol{w})^{-1}\boldsymbol{r}/G
at fixed \boldsymbol{w}: the scale is the mean squared Mahalanobis
residual, and it can be profiled out. The Fisher block for a is
I\_{aa}=G/(2a^{2}), the usual variance-parameter information, so the
relative standard error of \hat a is \sqrt{2/G}. With G in the hundreds
that is a few percent, which is the “fairly weak dependence” `RNA-Sieve`
reports for n.

### Layer 4: diagonal bulk measurement noise

Add technical bulk noise that is independent across genes,

\boldsymbol{y}\_{\cdot i}
=a_i\sum_j\boldsymbol{T}\_{ji}+\boldsymbol{\epsilon}\_i, \qquad
\boldsymbol{\epsilon}\_i\sim\mathcal{N}\_G(\boldsymbol{0},\boldsymbol{D}\_i),
\qquad \boldsymbol{\Sigma}(\boldsymbol{w},a,\boldsymbol{\tau}^{2})
=a\boldsymbol{S}\_W(\boldsymbol{w})+\boldsymbol{D}\_i . \tag{23}

\boldsymbol{D}\_i is diagonal by construction: gene to gene correlation
is already carried by the \boldsymbol{\Sigma}\_{W,j}, and an
off-diagonal residual would count the network twice ([perspectives,
three-layer
law](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-enough-replicates)).
The new derivatives are
\boldsymbol{\Sigma}\_{\tau_g^{2}}=\boldsymbol{e}\_g\boldsymbol{e}\_g^{\top}
and all second derivatives involving \boldsymbol{\tau}^{2} vanish, so

\frac{\partial\ell}{\partial\tau_g^{2}}
=-\tfrac12\\\Omega\_{gg}+\tfrac12\\(\boldsymbol{\Omega}\boldsymbol{r})\_g^{2}.
\tag{24}

The closed form of [Eq. 22](#eq-cc-score-a) no longer holds, because
\boldsymbol{\Omega}\boldsymbol{S}\_W is not a multiple of the identity.
Identifiability now rests on structure rather than on magnitude: a_i
inflates every entry of a_i\boldsymbol{S}\_W(\boldsymbol{w}),
\boldsymbol{D}\_i inflates the diagonal only. The **off-diagonal**
entries of the bulk covariance, which are the gene network, identify
a_i; the diagonal then identifies \boldsymbol{D}\_i. A gene-wise
deconvolution has no off-diagonal entries and cannot separate the two;
`MuSiC` and `RNA-Sieve` accordingly keep a single gene-wise variance.
With one bulk column and G free \tau\_{gi}^{2} the model is
over-parametrised; constrain \boldsymbol{D}\_i=\tau_i^{2}\mathbf{I}\_G,
or \tau\_{gi}^{2}=\tau_i^{2}\\\mu_g(\boldsymbol{w})^{\gamma} with a
mean–variance exponent, or estimate \boldsymbol{D} from technical bulk
replicates.

### Layer 5: several single-cell replicates

With R\ge 2 biological replicates in the reference, write
\boldsymbol{X}\_{jrc}=\boldsymbol{\mu}\_{\cdot
j}+\boldsymbol{B}\_{jr}+\boldsymbol{\varepsilon}\_{jrc} as in the
perspectives vignette, with
\boldsymbol{B}\_{jr}\sim\mathcal{N}(\boldsymbol{0},\boldsymbol{\Sigma}\_{B,j})
shared by all cells of type j in replicate r. A new bulk donor i draws
its own \boldsymbol{B}\_{ji}, so

\boldsymbol{T}\_{ji} =m\_{ji}(\boldsymbol{\mu}\_{\cdot
j}+\boldsymbol{B}\_{ji})
+\sum\_{c=1}^{m\_{ji}}\boldsymbol{\varepsilon}\_{jic}, \qquad
\operatorname{Cov}(\boldsymbol{T}\_{ji})
=m\_{ji}\boldsymbol{\Sigma}\_{W,j}+m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j},
\tag{25}

and the full five-layer law is

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{m}\_{\cdot i},a_i
\sim\mathcal{N}\_G\Bigl( a_i\boldsymbol{\mu}\boldsymbol{m}\_{\cdot i},\\
a_i^{2}\bigl\[\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
i})+\boldsymbol{S}\_B(\boldsymbol{m}\_{\cdot i})\bigr\]
+\boldsymbol{D}\_i \Bigr). \tag{26}

The quadratic term is DeCovarT’s original covariance written in cell
counts. In the \boldsymbol{m} chart the new derivatives are

\boldsymbol{\Sigma}\_{m_k}
=a^{2}\bigl(\boldsymbol{\Sigma}\_{W,k}+2m_k\boldsymbol{\Sigma}\_{B,k}\bigr),
\qquad \boldsymbol{\Sigma}\_{m_km_k}=2a^{2}\boldsymbol{\Sigma}\_{B,k},
\qquad \boldsymbol{\Sigma}\_{m_km_l}=\boldsymbol{0}\\ (k\neq l),
\tag{27}

so the two Hessian terms that vanished in [Sec. 3.1](#sec-cc-linear)
return, carried by \boldsymbol{\Sigma}\_{B,k} alone. Estimation of
\boldsymbol{\Sigma}\_{B,j} follows the method-of-moments contrast of the
[perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-enough-replicates):
the replicate means have covariance
\boldsymbol{\Sigma}\_{B,j}+\boldsymbol{\Sigma}\_{W,j}/C\_{jr}. A dense
G\times G between-replicate covariance needs many more replicates than a
single-cell study provides; with R around ten, independent donors
justify a **diagonal** \boldsymbol{\Sigma}\_{B,j} whose G marginal
variances are estimable, and the gene network stays inside
\boldsymbol{\Sigma}\_{W,j}.

### Comparison with MuSiC, RNA-Sieve and DWLS

**`MuSiC`.** Equation (2) of Wang et al.
([2019](#ref-wangBulkTissueCell2019)) is the mean of
[Eq. 13](#eq-cc-law) with \boldsymbol{\mu}\_{\cdot
j}=S_j\boldsymbol{\theta}\_{\cdot j}. `MuSiC` then divides by the
library total, which removes m_j and leaves relative abundances; it
never estimates counts. Its variance, equations (9) and (10),
\operatorname{Var}(Y\_{jg}\mid p_j)=C_j^{2}\delta\_{jg}^{2}
+C_j^{2}\sum_k p\_{jk}^{2}S_k^{2}\sigma\_{gk}^{2}, has a normalising
constant squared and **quadratic** weights p\_{jk}^{2}. The quadratic
weights arise because \theta\_{kjg} is a subject-level random profile
with cross-subject variance \sigma\_{gk}^{2}: that is the
\boldsymbol{\Sigma}\_{B,j} layer of [Eq. 26](#eq-cc-five-layer),
restricted to its diagonal, and \delta\_{jg}^{2} is \boldsymbol{D}\_i.
The within-type cell scatter \boldsymbol{\Sigma}\_{W,j} has no
counterpart in the `MuSiC` variance; it is used only to select
consistent genes.

**`RNA-Sieve`.** The bulk is a sum of n cells and the likelihood,
equation (8) of Erdmann-Pham et al.
([2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)), is
\prod_g\mathcal{N}\bigl(\tilde
b_g;\\n(M\alpha)\_g,\\n\sigma_g^{2}(M,\alpha,S)\bigr) with
\sigma_g^{2}=\sum_k\alpha_k\[s\_{g,k}+m\_{g,k}^{2}\]-b_g^{2}. Two
differences with [Eq. 13](#eq-cc-law) matter. First, the variance is
linear in n, as here, but it is the variance of a **randomly labelled**
cell: the term \sum_k\alpha_k m\_{g,k}^{2}-b_g^{2} is the spread of type
means and comes from drawing each cell’s type at random.
[Eq. 13](#eq-cc-law) conditions on the counts and therefore omits that
term; the multinomial-label version adds n\sum_j
p_j(\boldsymbol{\mu}\_{\cdot
j}-\bar{\boldsymbol{\mu}})(\boldsymbol{\mu}\_{\cdot
j}-\bar{\boldsymbol{\mu}})^{\top}, a rank-(J-1) matrix that the
cell-count model deliberately does not carry. Second, `RNA-Sieve` treats
genes as independent and uses a Godambe sandwich for intervals;
[Eq. 13](#eq-cc-law) keeps the full \boldsymbol{\Sigma}\_{W,j}.
`RNA-Sieve` also models the reference means as noisy, \tilde
M=M+\epsilon_M with variance s\_{g,k}/c_k. The analogue here is an
errors-in-variables term \sum_j
m\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/C_j on the bulk covariance,
quadratic in the counts and shrinking with the number of reference
cells, which [Sec. 6](#sec-summary) lists as an alternative.

**`DWLS`.** Tsoucas et al.
([2019](#ref-tsoucasAccurateEstimationCelltype2019)) average cells to a
signature, build an artificial bulk by summing single-cell profiles, and
fit weighted least squares with weights
1/(\boldsymbol{\mu}\hat{\boldsymbol{p}})\_g^{2} plus a dampening
constant. In the present notation that is [Eq. 23](#eq-cc-noise-law)
with \boldsymbol{\Sigma}\_{W,j}=\boldsymbol{0} and
\boldsymbol{D}\_i\propto\operatorname{diag}(\boldsymbol{\mu}\boldsymbol{p})^{2},
a multiplicative noise model without a cell-count scale. The paper
states the equal-total-RNA-per-cell assumption under which its weights
are cell fractions rather than RNA fractions; that is the S_j question
of [Sec. 5](#sec-compositional).

## Strategy 2: averaged mixture on the simplex

The second formulation divides both sides by the number of cells and
targets the composition \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1}. The
reference is already a per-cell average ([mean
signature](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-mean-signature)).
The bulk is not: a library is a sum, and the number of cells behind it
is unknown and usually far from the number of reference cells. Each
layer below adds one way of closing that gap.

### Layer 1: the per-cell average and the missing n_i

Divide [Eq. 12](#eq-cc-sum) by n_i:

\bar{\boldsymbol{y}}\_{\cdot i}
=\frac{1}{n_i}\sum_j\sum\_{c=1}^{m\_{ji}}\boldsymbol{X}\_{jic}, \qquad
\bar{\boldsymbol{y}}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},n_i
\sim\mathcal{N}\_G\Bigl( \boldsymbol{\mu}\boldsymbol{p}\_{\cdot i},\\
\kappa_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i}) \Bigr),
\qquad\kappa_i=\frac{1}{n_i}. \tag{28}

The mean is the DeCovarT mean. The covariance is the linear mixture of
the [perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-one-replicate)
with the ideal value \kappa_i=1/n_i: averaging n_i independent cells
shrinks the cell-level covariance by 1/n_i. The quadratic weights also
appear algebraically,
\operatorname{Cov}(p\_{ji}\bar{\boldsymbol{X}}\_{ji})
=p\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/m\_{ji}=p\_{ji}\boldsymbol{\Sigma}\_{W,j}/n_i,
but the covariance of the type-specific average depends on its own cell
count, which turns the square into a linear weight.

n_i is unobserved. Three ways to supply it:

- **Variance route (`RNA-Sieve`).** Estimate \kappa_i jointly with
  \boldsymbol{p}\_{\cdot i} from the covariance magnitude, as in
  1.  This needs \boldsymbol{\Sigma}\_{W,j} on the molecule scale of the
      bulk.
- **Pairing route (`Bisque`).** Jew et al.
  ([2020](#ref-jewAccurateEstimationCell2020)) form a pseudo-bulk
  \boldsymbol{Y}=\boldsymbol{Z}\boldsymbol{p} from the reference profile
  and the proportions counted in the single-cell data of paired
  individuals, then fit a per-gene linear map
  Y_j=\beta_jX'\_j+\epsilon_j (their equation 1) from observed bulk to
  pseudo-bulk, or match first and second moments gene by gene (their
  equation 3). A multivariate Gaussian version of that step is natural
  here: for a paired sample i' with counted \boldsymbol{m}\_{\cdot i'},
  the pseudo-bulk law is [Eq. 13](#eq-cc-law) and the observed bulk is
  an affine image of it, \boldsymbol{y}\_{\cdot
  i'}\sim\mathcal{N}\_G\bigl(\boldsymbol{c}+\boldsymbol{B}\boldsymbol{\mu}\boldsymbol{m}\_{\cdot
  i'}, \boldsymbol{B}\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
  i'})\boldsymbol{B}^{\top}+\boldsymbol{D}\bigr) with diagonal
  \boldsymbol{B}; (\boldsymbol{c},\boldsymbol{B},\boldsymbol{D}) are
  estimated once on the paired cohort and \kappa_i for an unpaired
  sample follows from the fitted mean scale
  ([Sec. 7.2](#sec-persp-affine)).
- **External counts.** Nuclei counts, DNA content or tissue mass give
  n_i directly; then \kappa_i is a known offset.

With \boldsymbol{\vartheta}=(\boldsymbol{p},\kappa) the tuple is

\boldsymbol{\mu}\_{p_k}=\boldsymbol{\mu}\_{\cdot k},\quad
\boldsymbol{\Sigma}\_{p_k}=\kappa\boldsymbol{\Sigma}\_{W,k},\quad
\boldsymbol{\Sigma}\_{\kappa}=\boldsymbol{S}\_W(\boldsymbol{p}),\quad
\boldsymbol{\Sigma}\_{p_k\kappa}=\boldsymbol{\Sigma}\_{W,k},\quad
\boldsymbol{\Sigma}\_{p_kp_l}=\boldsymbol{\Sigma}\_{\kappa\kappa}=\boldsymbol{0},
\tag{29}

with all \boldsymbol{\mu}\_{kl}=\boldsymbol{0}. The simplex is imposed
on the \boldsymbol{p} block by the ILR chart,
\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{\rho}), with
Jacobian \mathbf{J}\_\psi=\mathbf{S}(\boldsymbol{p})\mathbf{V} and the
score contraction of the map Hessian exactly as in the [chain
rule](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-generative-model.html#sec-chain-rule);
the \kappa row is untouched. In particular the stationarity condition
\mathbf{V}^{\top}\mathbf{S}(\boldsymbol{p})\nabla\_{\boldsymbol{p}}\ell=\boldsymbol{0}
is the same equation with \nabla\_{\boldsymbol{p}}\ell taken from
[Eq. 8](#eq-score-general) and the tuple [Eq. 29](#eq-av-tuple).

### Layer 2: a global scale, per sample or shared

On an absolute scale the bulk is \boldsymbol{y}\_{\cdot
i}=a_i\bar{\boldsymbol{y}}\_{\cdot i} with the same a_i as in 1, now
absorbing n_i as well as capture:

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},a_i,\kappa_i
\sim\mathcal{N}\_G\Bigl( a_i\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i},\\
a_i^{2}\kappa_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i}) \Bigr).
\tag{30}

Here the simplex constraint buys identifiability. Because
\mathbf{1}^{\top}\boldsymbol{p}\_{\cdot i}=1, the mean
a_i\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i} pins a_i from first moments
alone, and \kappa_i is then read from the covariance magnitude. In
Strategy 1 the mean fixes only the product a_i\boldsymbol{m}\_{\cdot i}.
The tuple adds

\boldsymbol{\mu}\_{a}=\boldsymbol{\mu}\boldsymbol{p},\quad
\boldsymbol{\mu}\_{p_k}=a\boldsymbol{\mu}\_{\cdot k},\quad
\boldsymbol{\mu}\_{p_ka}=\boldsymbol{\mu}\_{\cdot k},\quad
\boldsymbol{\Sigma}\_{a}=2a\kappa\boldsymbol{S}\_W(\boldsymbol{p}),\quad
\boldsymbol{\Sigma}\_{aa}=2\kappa\boldsymbol{S}\_W(\boldsymbol{p}),\quad
\boldsymbol{\Sigma}\_{p_ka}=2a\kappa\boldsymbol{\Sigma}\_{W,k},\quad
\boldsymbol{\Sigma}\_{a\kappa}=2a\boldsymbol{S}\_W(\boldsymbol{p}),
\tag{31}

and the p and \kappa entries of [Eq. 29](#eq-av-tuple) are multiplied by
a^{2}. When N bulk samples share one protocol, a **shared** scale a
replaces a_i and the pooled log-likelihood
\sum_i\ell_i(\boldsymbol{p}\_{\cdot i},a,\kappa_i) is maximised; its a
score is the sum of the per-sample scores, and the Hessian is block
arrow-shaped with the a row coupling all samples. This pooling is what
`MuSiC` achieves with its sample-level normalising constant C_j and what
`RCTD` achieves with a known depth offset ([probabilistic
engines](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-probabilistic-engines)).

### Layer 3: a residual bulk term

Adding
\boldsymbol{\epsilon}\_i\sim\mathcal{N}\_G(\boldsymbol{0},\boldsymbol{D}\_i)
gives

\boldsymbol{\Sigma}(\boldsymbol{p},a,\kappa,\boldsymbol{\tau}^{2})
=a^{2}\kappa\boldsymbol{S}\_W(\boldsymbol{p})+\boldsymbol{D}\_i,
\tag{32}

with
\boldsymbol{\Sigma}\_{\tau_g^{2}}=\boldsymbol{e}\_g\boldsymbol{e}\_g^{\top}
as in [Eq. 24](#eq-cc-score-tau). \kappa_i and \boldsymbol{D}\_i both
inflate the variance and are separated by structure, not by magnitude:
the network entries of \boldsymbol{S}\_W(\boldsymbol{p}) identify
\kappa_i, the diagonal excess identifies \boldsymbol{D}\_i. If the
reference network is weak, a free \boldsymbol{D}\_i will absorb
\kappa_i, and n_i becomes unidentifiable; constrain \boldsymbol{D}\_i as
in [Sec. 3.4](#sec-cc-noise) or fix \kappa_i externally.

### Layer 4: within and between covariance from R replicates

With R\ge 2 replicates, divide [Eq. 25](#eq-cc-replicate-type) by n_i:
the between term
m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}/n_i^{2}=p\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}
keeps its quadratic weight, the within term becomes
\kappa_ip\_{ji}\boldsymbol{\Sigma}\_{W,j}. The four-layer law is the
three-layer law of the perspectives vignette with the scale written out,

\boldsymbol{y}\_{\cdot i}\mid\boldsymbol{p}\_{\cdot i},a_i,\kappa_i
\sim\mathcal{N}\_G\Bigl( a_i\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i},\\
a_i^{2}\bigl\[\boldsymbol{S}\_B(\boldsymbol{p}\_{\cdot
i})+\kappa_i\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot i})\bigr\]
+\boldsymbol{D}\_i \Bigr), \tag{33}

and the covariance derivatives in \boldsymbol{p} are

\boldsymbol{\Sigma}\_{p_k}
=a^{2}\bigl(2p_k\boldsymbol{\Sigma}\_{B,k}+\kappa\boldsymbol{\Sigma}\_{W,k}\bigr),
\qquad \boldsymbol{\Sigma}\_{p_kp_k}=2a^{2}\boldsymbol{\Sigma}\_{B,k},
\qquad \boldsymbol{\Sigma}\_{p_kp_l}=\boldsymbol{0}\\ (k\neq l).
\tag{34}

Setting \kappa_i\to 0 (a bulk of infinitely many cells),
\boldsymbol{D}\_i=\boldsymbol{0} and a_i=1 returns the released DeCovarT
derivatives with \boldsymbol{\Sigma}\_j=\boldsymbol{\Sigma}\_{B,j}.
Released DeCovarT is therefore the large-n_i, noise-free limit of this
formulation, which is why its \boldsymbol{\Sigma}\_j must be a
between-replicate covariance to carry a sampling interpretation. The
diagonal-\boldsymbol{\Sigma}\_{B,j} remark of
[Sec. 3.5](#sec-cc-replicates) applies: a diagonal quadratic layer plus
a sparse-precision linear layer has a dense precision, so the backends
of the derivatives vignette see only the structure of
\boldsymbol{\Sigma}.

### Comparison with DECALS

`DECALS` is the method closest to released DeCovarT ([Cai et al.
2024](#ref-caiStatisticalInferenceCelltype2024)). Its model (their
equation 2) is y\_{ij}=\sum_k\pi\_{ik}w\_{kj}+\epsilon\_{ij} on FPKM,
fitted by constrained ordinary least squares (their equation 3), and its
subject-specific covariance under independent cell-type profiles is
\boldsymbol{\Sigma}\_i=\operatorname{Cov}(\sum_k\pi\_{ik}\boldsymbol{x}\_i^{(k)})
=\sum_k\pi\_{ik}^{2}\boldsymbol{\Sigma}^{(k)} (their Section 2.3): the
quadratic layer \boldsymbol{S}\_B(\boldsymbol{p}) and nothing else.
Three contrasts with [Eq. 33](#eq-av-four-layer):

- **Where \boldsymbol{\Sigma}^{(k)} comes from.** `DECALS` regresses the
  residual outer products \hat z\_{ij}\hat z\_{ij'} on
  (\hat\pi\_{i1}^{2},\ldots,\hat\pi\_{iK}^{2}) across the n bulk samples
  (their equation 6). The design matrix
  \boldsymbol{H}=\[\hat\pi\_{ik}^{2}\] must have full column rank, so it
  needs N\ge J bulk samples with varying composition. DeCovarT plugs in
  reference estimates and needs no bulk replication for the covariance.
  The two sources are complementary: a residual regression on N bulk
  samples is a route to \boldsymbol{\Sigma}\_{B,j} when R=1 in the
  reference.
- **Estimator.** `DECALS` keeps constrained OLS for the point estimate
  and uses \boldsymbol{\Sigma}\_i only for the asymptotic covariance
  V_i=UDU^{\top} with
  D=(\tfrac1pW^{\top}W)^{-1}\tfrac1pW^{\top}\boldsymbol{\Sigma}\_iW(\tfrac1pW^{\top}W)^{-1}.
  It considers constrained GLS with \boldsymbol{\Sigma}\_i^{-1/2} and
  rejects it because the estimated inverse square root is too noisy per
  sample. DeCovarT maximises the full Gaussian likelihood, whose
  determinant term and covariance-quadratic term are absent from any GLS
  criterion ([GLS
  competitor](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-robust-gls)).
- **Finite-sample correction.** The `DECALS` bias correction (their
  equation 7, Proposition 1) targets the moment estimator of
  \boldsymbol{\Sigma}^{(k)}, not the proportion estimator.
  [Sec. 7.3](#sec-persp-firth) sets it beside Firth’s adjustment.

> **Important 4: Strategies 1 and 2 are one model with one redundant
> coordinate**
>
> Start from the five-layer law [Eq. 26](#eq-cc-five-layer) and set
> \boldsymbol{m}\_{\cdot i}=n_i\boldsymbol{p}\_{\cdot i}:
>
> \begin{aligned} a_i\boldsymbol{\mu}\boldsymbol{m}\_{\cdot i}
> &=(a_in_i)\\\boldsymbol{\mu}\boldsymbol{p}\_{\cdot i}, \\
> a_i^{2}\bigl\[\boldsymbol{S}\_W(\boldsymbol{m}\_{\cdot
> i})+\boldsymbol{S}\_B(\boldsymbol{m}\_{\cdot i})\bigr\]
> &=(a_in_i)^{2}\Bigl\[\frac{1}{n_i}\boldsymbol{S}\_W(\boldsymbol{p}\_{\cdot
> i})+\boldsymbol{S}\_B(\boldsymbol{p}\_{\cdot i})\Bigr\]. \end{aligned}
> \tag{35}
>
> Renaming a_in_i\mapsto a_i and 1/n_i\mapsto\kappa_i gives
> [Eq. 33](#eq-av-four-layer) exactly. The map (\boldsymbol{m}\_{\cdot
> i},a_i)\mapsto(\boldsymbol{p}\_{\cdot i},\kappa_i,a_in_i) is a
> bijection from (0,\infty)^{J}\times(0,\infty) onto
> \Delta^{J-1}\times(0,\infty)^{2}, so the two strategies are
> reparametrisations of the same Gaussian family. The cellular ratio
> p\_{ji} is the count m\_{ji} averaged by n_i, and the quadratic
> DeCovarT weight is the between-population layer written per cell:
> m\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}/n_i^{2}=p\_{ji}^{2}\boldsymbol{\Sigma}\_{B,j}.
>
> What differs is how the redundant degree of freedom is spent. Strategy
> 1 has J+1 free parameters (\boldsymbol{m},a) and the mean identifies
> only a\boldsymbol{m}; the count scale comes from the covariance.
> Strategy 2 has (J-1)+2 free parameters (\boldsymbol{p},\kappa,a); the
> simplex removes one coordinate and lets the mean identify a outright,
> leaving \kappa to the covariance. The integer constraint of
> [Sec. 3.2](#sec-cc-integer) exists only in Strategy 1, because
> n_i\boldsymbol{p}\_{\cdot i} is not constrained to the lattice once
> n_i is a continuous nuisance. Choose by the question: counts and a
> scale calibrated to the reference platform (Strategy 1), or a
> composition with n_i as a nuisance (Strategy 2).

## Strategy 3: library-size normalisation gives a composition

The third formulation divides the bulk by its own total rather than by a
cell count, as CPM, TPM, RPKM and FPKM do. Let
\boldsymbol{T}\_i=\sum_j\sum_c\boldsymbol{X}\_{jic} be the unnormalised
vector and

\boldsymbol{Z}\_i=\frac{\boldsymbol{T}\_i}{\mathbf{1}^{\top}\boldsymbol{T}\_i},
\qquad \mathbf{1}^{\top}\boldsymbol{Z}\_i=1 . \tag{36}

\boldsymbol{Z}\_i is a random ratio. Even when \boldsymbol{T}\_i is
exactly Gaussian, \boldsymbol{Z}\_i is not; its support is the
composition hyperplane and a Gaussian description lives on the
(G-1)-dimensional tangent space. With
\boldsymbol{t}\_i=\mathbb{E}\[\boldsymbol{T}\_i\],
\boldsymbol{V}\_i=\operatorname{Cov}(\boldsymbol{T}\_i),
L_i=\mathbf{1}^{\top}\boldsymbol{t}\_i and
\boldsymbol{z}\_i=\boldsymbol{t}\_i/L_i, the delta method gives

\boldsymbol{J}\_i=\frac{1}{L_i}\bigl(\mathbf{I}\_G-\boldsymbol{z}\_i\mathbf{1}^{\top}\bigr),
\qquad
\boldsymbol{Z}\_i\approx\mathcal{N}\bigl(\boldsymbol{z}\_i,\boldsymbol{J}\_i\boldsymbol{V}\_i\boldsymbol{J}\_i^{\top}\bigr),
\qquad \mathbf{1}^{\top}\boldsymbol{J}\_i=\boldsymbol{0}. \tag{37}

The covariance is singular in \mathbb{R}^{G}. Closure also **creates**
correlation. For a diagonal
\boldsymbol{V}\_i=\operatorname{diag}(\boldsymbol{v}),

\bigl(\boldsymbol{J}\_i\boldsymbol{V}\_i\boldsymbol{J}\_i^{\top}\bigr)\_{gg'}
=\frac{1}{L_i^{2}}\Bigl(
v_g\delta\_{gg'}-z_gv\_{g'}-z\_{g'}v_g+z_gz\_{g'}\sum_hv_h \Bigr),
\tag{38}

so independent genes become negatively correlated after normalisation,
and the gene network estimated on raw counts is no longer the network of
the normalised bulk. The target also changes: with
\boldsymbol{\mu}\_{\cdot j}=S_j\boldsymbol{\theta}\_{\cdot j} the
normalised mean is \sum_jq\_{ji}\boldsymbol{\theta}\_{\cdot j} with RNA
fractions q\_{ji}=m\_{ji}S_j/\sum_km\_{ki}S_k, which equal cell
fractions only when S_j is constant ([RNA fraction versus cell
fraction](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-uncoupling)).
This vignette does not develop Strategy 3 further: it removes the count
scale that Strategies 1 and 2 estimate, breaks the additive mixture that
the convolution needs ([normalisation that preserves the linear
mixture](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-norm-linear)),
and is the regime in which `MuSiC` and `DECALS` operate by choice.

## Summary and critical comparison

|  | Strategy 1: cell counts | Strategy 2: averaged mixture | Strategy 3: compositional | Released DeCovarT |
|----|----|----|----|----|
| Target | \boldsymbol{m}\_{\cdot i}\in\mathbb{N}^{J}, then \boldsymbol{p}=\boldsymbol{m}/n | \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1}, n_i nuisance | RNA fractions \boldsymbol{q}\_{\cdot i} | \boldsymbol{p}\_{\cdot i}\in\Delta^{J-1} |
| Covariance weights | linear m\_{ji} in \boldsymbol{\Sigma}\_{W,j}; quadratic m\_{ji}^{2} in \boldsymbol{\Sigma}\_{B,j} | linear \kappa_ip\_{ji}; quadratic p\_{ji}^{2} | delta-method image of either | quadratic p\_{ji}^{2} only |
| Constraint | none; lattice | simplex via ILR | closure, singular | simplex via ILR |
| Scale identification | mean fixes a\boldsymbol{m}; covariance fixes a | mean fixes a; covariance fixes \kappa | scale removed | a absorbed in \boldsymbol{y} |
| Advantages | counts are the physical quantity; zeros are ordinary points; Hessian simpler (linear covariance); one Cholesky per iteration | DeCovarT is its noise-free limit; simplex buys identifiability of a; direct comparison with MuSiC, DECALS | matches CPM/TPM pipelines; no scale to estimate | implemented; exact Gaussian; sparse \boldsymbol{\Omega}\_j as declared input |
| Drawbacks | a needs \boldsymbol{\Sigma}\_{W,j} on the reference molecule scale; cross-protocol a loses physical meaning; integer search for rare types | \kappa and \boldsymbol{D} separated only through off-diagonal structure; n_i weakly identified when the network is weak | induced correlations; \boldsymbol{q}\neq\boldsymbol{p} unless S_j constant; Gaussian only approximately | \boldsymbol{\Sigma}\_j must be a between-replicate covariance; silences \boldsymbol{\Sigma}\_{W,j} for rare types |

Table 2: The three formulations and the released model.

The honest reading of [Table 2](#tbl-strategies) is that Strategy 2 with
the four layers of [Sec. 4.4](#sec-av-replicates) is the model DeCovarT
should fit when a single-cell reference is used, and Strategy 1 is its
count-scale chart. Four further alternatives stay inside the
multivariate Gaussian regime and deserve a place in a simulation study
before any is implemented.

1.  **Errors in variables for the signature.** The reference mean is
    itself an average of C_j cells, so \hat{\boldsymbol{\mu}}\_{\cdot
    j}\sim\mathcal{N}(\boldsymbol{\mu}\_{\cdot
    j},\boldsymbol{\Sigma}\_{W,j}/C_j). Marginalising it adds
    \sum_jp\_{ji}^{2}\boldsymbol{\Sigma}\_{W,j}/C_j to
    [Eq. 33](#eq-av-four-layer): a quadratic within-type term that is
    Monte Carlo error of the signature, not biology, and that vanishes
    as the reference grows. `RNA-Sieve` and `MEAD` carry it
    ([Erdmann-Pham et al.
    2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021); [Xie
    and Wang 2023](#ref-xieRobustStatisticalInference2023)). It is the
    one quadratic term that a one-replicate reference does identify.
2.  **Residual regression for \boldsymbol{\Sigma}\_{B,j}.** When R=1 but
    N\ge J bulk samples share the reference, the `DECALS` moment
    regression of residual outer products on \hat p\_{ji}^{2} estimates
    a between-population covariance without single-cell replicates. A
    diagonal restriction keeps it estimable.
3.  **Multi-sample pooling.** Shared a across samples
    ([Sec. 4.2](#sec-av-scale)), shared \boldsymbol{D} across technical
    replicates, or shared \kappa across samples of the same tissue mass
    turn weakly identified per-sample scalars into pooled ones.
4.  **Compositional Gaussian on the bulk side.** If a pipeline insists
    on normalised bulk, model \operatorname{clr}(\boldsymbol{Z}\_i) or
    \operatorname{ilr}(\boldsymbol{Z}\_i) as Gaussian (a logistic-normal
    bulk) with mean the log-ratio image of the mixture ([Aitchison
    1982](#ref-aitchisonStatisticalAnalysisCompositional1982);
    [Pawlowsky-Glahn and Buccianti
    2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)). The mean
    map is no longer linear in \boldsymbol{p}, so the analytic score of
    [Theorem 1](#thm-score-general) still applies but
    \boldsymbol{\mu}\_{kl}\neq\boldsymbol{0}.

Outside the Gaussian regime, a Poisson–log-normal read-level mixture
(`RCTD`) is the natural competitor and is already described in the
[perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-probabilistic-engines).

## Perspectives

### Correlated cells inside a type

[Sec. 3](#sec-cell-counts) assumed independent cells within a type.
Cells of one type in one tissue share micro-environment, cell cycle
phase and paracrine signals, so a positive within-type correlation is
plausible. The jointly Gaussian equicorrelated model of
[Eq. 5](#eq-equicorrelation) with a type-specific \varrho_j gives

\operatorname{Cov}(\boldsymbol{T}\_{ji})
=m\_{ji}\bigl\[1+(m\_{ji}-1)\varrho_j\bigr\]\boldsymbol{\Sigma}\_{W,j}
=m\_{ji}(1-\varrho_j)\boldsymbol{\Sigma}\_{W,j}+m\_{ji}^{2}\varrho_j\boldsymbol{\Sigma}\_{W,j}.
\tag{39}

The second form shows that equicorrelation is the random-effects model
[Eq. 6](#eq-random-effects) with proportional layers,
\boldsymbol{\Sigma}\_{B,j}=\varrho_j\boldsymbol{\Sigma}\_{W,j} and a
within covariance deflated to (1-\varrho_j)\boldsymbol{\Sigma}\_{W,j}.
The covariance derivative in \varrho_j is
\boldsymbol{\Sigma}\_{\varrho_j}=m\_{ji}(m\_{ji}-1)\boldsymbol{\Sigma}\_{W,j},
so [Theorem 1](#thm-score-general) applies. Identifiability is the
obstacle: with one bulk sample, \varrho_j and m\_{ji} enter only through
the scalar m\_{ji}\[1+(m\_{ji}-1)\varrho_j\] in front of
\boldsymbol{\Sigma}\_{W,j} and cannot be separated. Several samples with
different compositions, or an external estimate of \varrho_j from
spatial neighbourhoods of the same type, are needed. The lower bound
\varrho_j\ge-1/(m\_{ji}-1) is the positive-semidefiniteness constraint
and shrinks to zero for abundant types, so negative within-type
correlation is a rare-type phenomenon only. The claim that a sum of
equicorrelated cells is Gaussian rests on **joint** Gaussianity of the
stacked cells, which is an assumption and not a consequence of Gaussian
marginals.

### Gene-wise and type-wise platform scaling

1 used one scalar a_i. Capture efficiency in single-cell protocols is
gene dependent and technology dependent, so the alignment between
reference and bulk is in general an affine map. Let
\boldsymbol{c}\in\mathbb{R}^{G} be a background intercept and
\boldsymbol{B}\_j=\operatorname{diag}(\boldsymbol{d}\_j) a gene-wise,
type-wise scale. Affine invariance of the Gaussian family gives exactly

\boldsymbol{y}\_{\cdot i}
=\boldsymbol{c}+\sum_j\boldsymbol{B}\_j\boldsymbol{T}\_{ji}+\boldsymbol{\epsilon}\_i
\sim\mathcal{N}\_G\Bigl(
\boldsymbol{c}+\sum_jm\_{ji}\boldsymbol{B}\_j\boldsymbol{\mu}\_{\cdot
j},\\
\sum_jm\_{ji}\boldsymbol{B}\_j\boldsymbol{\Sigma}\_{W,j}\boldsymbol{B}\_j^{\top}+\boldsymbol{D}\_i
\Bigr), \tag{40}

and the tuple of [Theorem 1](#thm-score-general) follows by replacing
\boldsymbol{\mu}\_{\cdot j} with
\boldsymbol{B}\_j\boldsymbol{\mu}\_{\cdot j} and
\boldsymbol{\Sigma}\_{W,j} with
\boldsymbol{B}\_j\boldsymbol{\Sigma}\_{W,j}\boldsymbol{B}\_j. Three
special cases are already in the literature.
\boldsymbol{B}\_j=a_i\mathbf{I}\_G, \boldsymbol{c}=\boldsymbol{0} is 1.
\boldsymbol{c}=(\alpha\_{jg})\_g with \boldsymbol{B}\_j=\mathbf{I}\_G is
the gene- and subject-specific intercept `MuSiC` adds to its equation
(8) to absorb protocol bias ([Wang et al.
2019](#ref-wangBulkTissueCell2019)).
\boldsymbol{B}\_j=\operatorname{diag}(\boldsymbol{d}) shared across
types, with a per-gene intercept, is the `Bisque` transformation learned
on paired individuals ([Jew et al.
2020](#ref-jewAccurateEstimationCell2020)) and the \boldsymbol{d} of
`MEAD` and `RCTD`
([alignment](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-alignment)).

Two warnings carry over. A free \boldsymbol{d} fitted on the sample
being deconvolved is not identified from \boldsymbol{p}\_{\cdot i}; it
must come from a paired cohort, a pseudo-bulk of all samples, or an
external measurement model. A type-wise \boldsymbol{d}\_j is weaker
still and is a candidate only when sorted bulk profiles of type j exist
on the bulk platform. The quantity that an external model could supply
is the gene-specific capture probability of the single-cell protocol.
Sarkar and Stephens separate that protocol from the biological state: an
observation model is a measurement model given the latent abundance,
plus an expression model for the abundance ([Sarkar and Stephens
2021](#ref-sarkarSeparatingMeasurementExpression2021)). `bayNorm` is a
binomial measurement model for the unobserved original count, so the sum
of posterior means over genes is a cell’s inferred transcript total and
the within-type mean of those totals is S_j ([Tang et al.
2020](#ref-tangBayNormBayesianGene2020)). `Sanity` is the
expression-model step: a posterior on transcriptional activity after
Poisson sampling noise, not a molecule count ([Breda et al.
2021](#ref-bredaBayesianInferenceGene2021)). The ratios of those S_j are
what the RNA-to-cell correction needs; one shared mean capture
efficiency cancels. Any such plug-in is protocol specific, so
\boldsymbol{B}\_j and S_j would be stored with the reference, never
re-estimated per bulk sample. The longer statement is in the
[perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-uncoupling).

### Firth penalisation and DECALS-type finite-sample correction

The Firth adjustment of the [perspectives
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-decovart-statistical-perspectives.html#sec-firth)
maximises
\ell\_{\mathrm{F}}(\boldsymbol{\vartheta})=\ell(\boldsymbol{\vartheta})+\tfrac12\log\det
I(\boldsymbol{\vartheta}) ([Firth
1993](#ref-firthBiasReductionMaximum1993)). With
[Eq. 10](#eq-fisher-general) the adjusted score is

\frac{\partial\ell\_{\mathrm{F}}}{\partial\vartheta_q}
=\frac{\partial\ell}{\partial\vartheta_q}
+\tfrac12\\\mathrm{tr}\Bigl(I(\boldsymbol{\vartheta})^{-1}\frac{\partial
I(\boldsymbol{\vartheta})}{\partial\vartheta_q}\Bigr). \tag{41}

For the cell-count model of [Sec. 3.1](#sec-cc-linear) the information
depends on \boldsymbol{m} only through \boldsymbol{\Omega}, and
\partial\boldsymbol{\Omega}/\partial
m_q=-\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,q}\boldsymbol{\Omega},
so

\frac{\partial I(\boldsymbol{m})\_{kl}}{\partial m_q}
=-\boldsymbol{\mu}\_{\cdot
k}^{\top}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,q}\boldsymbol{\Omega}\boldsymbol{\mu}\_{\cdot
l}
-\mathrm{tr}\bigl(\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,k}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,q}\boldsymbol{\Omega}\boldsymbol{\Sigma}\_{W,l}\bigr),
\tag{42}

where the two trace terms of the product rule coincide by transposition
and cyclicity. Both summands are negative semidefinite contributions:
information about counts **decreases** as counts grow, because more
cells inflate \boldsymbol{S}\_W(\boldsymbol{m}) and deflate
\boldsymbol{\Omega}. Jeffreys’ prior \sqrt{\det I} therefore favours
small totals, and the Firth penalty on the cell-count chart pulls \hat
n_i downward. That is a prior statement, not a proven bias reduction:
Firth’s exact removal of the O(N^{-1}) term holds for the canonical
parameter of a full exponential family, whereas a Gaussian whose mean
and covariance both depend on \boldsymbol{\vartheta} is a curved family,
where the Jeffreys penalty is one of several bias-reducing modified
scores and does not cancel the first-order bias exactly. On the simplex
chart of Strategy 2 the same derivative applies to the \boldsymbol{p}
block with \kappa fixed, and the collapse of \det I near faces and under
collinear signatures is the regularising effect already argued in the
perspectives vignette. A simulation under the ADEMP template of the
[synthetic-scenarios
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.md)
is the only way to decide whether the pull is a correction or a bias.

`DECALS` corrects a different estimator. Its moment regression estimates
the type-level covariances from residuals computed with
\hat{\boldsymbol{\pi}}\_i instead of \boldsymbol{\pi}\_i, and the
finite-sample gap between \hat{\boldsymbol{H}}=\[\hat\pi\_{ik}^{2}\] and
\boldsymbol{H}=\[\pi\_{ik}^{2}\] deflates the estimated covariance and
the coverage of its intervals. Proposition 1 of Cai et al.
([2024](#ref-caiStatisticalInferenceCelltype2024)) quantifies
\mathbb{E}\[\hat{\boldsymbol{H}}^{\top}\hat{\boldsymbol{H}}-\boldsymbol{H}^{\top}\boldsymbol{H}\]
and \mathbb{E}\[\hat{\boldsymbol{H}}-\boldsymbol{H}\] in terms of the
asymptotic covariance V_i of \hat{\boldsymbol{\pi}}\_i and subtracts
them (their equation 7). The DeCovarT analogue is the plug-in of
estimated covariances. If \hat{\boldsymbol{\Sigma}}\_{W,j} is the sample
covariance of C_j Gaussian cells, \nu_j=C_j-1, then
\mathbb{E}\[\hat{\boldsymbol{\Sigma}}\_{W,j}^{-1}\]=\tfrac{\nu_j}{\nu_j-G-1}\boldsymbol{\Sigma}\_{W,j}^{-1}
for \nu_j\>G+1: the plug-in precision is inflated by a factor that is
far from one for rare types with few cells and many genes. The
first-order correction multiplies \hat{\boldsymbol{\Omega}}\_{W,j} by
(\nu_j-G-1)/\nu_j; the graphical-lasso shrinkage of the reproducibility
book is the regularised alternative when \nu_j\le G+1. The two
corrections are complementary: Firth acts on the likelihood geometry of
the composition, the Wishart or `DECALS` correction acts on the
covariance that the likelihood plugs in. A `DECALS`-style iteration,
alternating between \hat{\boldsymbol{p}}\_{\cdot i} and a bias-corrected
\hat{\boldsymbol{\Sigma}}\_{B,j} regressed on \hat p\_{ji}^{2} across N
bulk samples, is the route to a between-population covariance when the
reference has a single replicate.

## References

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Breda, Jérémie, Mihaela Zavolan, and Erik van Nimwegen. 2021. ‘Bayesian
Inference of Gene Expression States from Single-Cell RNA-Seq Data’.
*Nature Biotechnology* 39 (8): 1008–16.
<https://doi.org/10.1038/s41587-021-00875-x>.

Cai, Biao, Emma Jingfei Zhang, Hongyu Li, Chang Su, and Hongyu Zhao.
2024. ‘Statistical Inference of Cell-Type Proportions Estimated from
Bulk Expression Data’. *Journal of the American Statistical Association*
119 (548): 2521–32. <https://doi.org/10.1080/01621459.2024.2382435>.

Erdmann-Pham, Dan D., Jonathan Fischer, Justin Hong, and Yun S. Song.
2021. ‘Likelihood-Based Deconvolution of Bulk Gene Expression Data Using
Single-Cell References’. *Genome Research* 31 (10): 1794–806.
<https://doi.org/10.1101/gr.272344.120>.

Firth, David. 1993. ‘Bias Reduction of Maximum Likelihood Estimates’.
*Biometrika* 80 (1): 27–38. <https://doi.org/10.1093/biomet/80.1.27>.

Jew, Brandon, Marcus Alvarez, Elior Rahmani, et al. 2020. ‘Accurate
Estimation of Cell Composition in Bulk Expression Through Robust
Integration of Single-Cell Information’. *Nature Communications* 11.
<https://doi.org/10.1038/s41467-020-15816-6>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Sarkar, Abhishek, and Matthew Stephens. 2021. ‘Separating Measurement
and Expression Models Clarifies Confusion in Single-Cell RNA Sequencing
Analysis’. *Nature Genetics* 53 (6): 770–77.
<https://doi.org/10.1038/s41588-021-00873-4>.

Tang, Wenhao, François Bertaux, Philipp Thomas, et al. 2020. ‘bayNorm:
Bayesian Gene Expression Recovery, Imputation and Normalization for
Single-Cell RNA-Sequencing Data’. *Bioinformatics* 36 (4): 1174–81.
<https://doi.org/10.1093/bioinformatics/btz726>.

Tsoucas, Daphne, Rui Dong, Haide Chen, Qian Zhu, Guoji Guo, and
Guo-Cheng Yuan. 2019. ‘Accurate Estimation of Cell-Type Composition from
Gene Expression Data’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-019-10802-z>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.

Xie, Dongyue, and Jingshu Wang. 2023. *Robust Statistical Inference for
Cell Type Deconvolution*. arXiv.
<https://doi.org/10.48550/arxiv.2202.06420>.
