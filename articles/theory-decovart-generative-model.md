# Derivatives of the DeCovarT generative model under simplex transforms

DeCovarT maximises a Gaussian-convolution log-likelihood for cellular
ratios \boldsymbol{p}\in\Delta^{J-1}. This vignette derives the analytic
first- and second-order derivatives used by the frequentist solvers,
first with respect to unconstrained ambient proportions \boldsymbol{p},
then after the isometric log-ratio (ILR) map that enforces the simplex.
It closes with the expected Fisher information and the ILR delta-method
covariance that feed
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
and
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).
Package entry points are
[`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md),
[`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md),
[`gradient_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_constrained.md),
[`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md),
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
and
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md).
The additive log-ratio chart is retained as an appendix
([Sec. 8](#sec-alr)) and as
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md)
for reference-invariance checks.

> **Note 1: Notation**
>
> - \boldsymbol{\mu}\in\mathbb{R}^{G\times J}: mean signature.
> - \boldsymbol{\Sigma}\_j\in\mathrm{SPD}(G): cell-type covariance;
>   \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
>   p_j^{2}\boldsymbol{\Sigma}\_j and
>   \boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}.
> - \boldsymbol{p}: proportions on the open simplex (ambient score
>   equations treat them as free coordinates in (0,1)^{J}).
> - \boldsymbol{z}\in\mathbb{R}^{J-1}: ILR coordinates on a Helmert
>   basis \mathbf{V}; \boldsymbol{p}=\operatorname{softmax}(\mathbf{V}
>   \boldsymbol{z}).

## Matrix calculus reminders

The score and Hessian of the convolution log-likelihood follow from a
handful of matrix-derivative identities. The statements below mirror the
matrix-calculus remarks in the manuscript Methods appendix.

*Remark 1* (First-order matrix calculus). Let
\boldsymbol{A}=\boldsymbol{A}(p) and \boldsymbol{B}=\boldsymbol{B}(p) be
invertible matrix-valued functions of a scalar p, and let
\boldsymbol{U}, \boldsymbol{V} be constant matrices of compatible size.
Then

\begin{aligned} \frac{\partial\det(\boldsymbol{A})}{\partial p} &=
\det(\boldsymbol{A})\\ \mathrm{tr}\\\left( \boldsymbol{A}^{-1}
\frac{\partial\boldsymbol{A}}{\partial p} \right), \\
\frac{\partial(\boldsymbol{U}\boldsymbol{A}\boldsymbol{V})}{\partial p}
&= \boldsymbol{U} \frac{\partial\boldsymbol{A}}{\partial p}
\boldsymbol{V}, \\ \frac{\partial\boldsymbol{A}^{-1}}{\partial p} &= -
\boldsymbol{A}^{-1} \frac{\partial\boldsymbol{A}}{\partial p}
\boldsymbol{A}^{-1}. \end{aligned} \tag{1}

The chain rule on \log\det yields

\frac{\partial\log\det(\boldsymbol{A})}{\partial p} =
\mathrm{tr}\\\left( \boldsymbol{A}^{-1}
\frac{\partial\boldsymbol{A}}{\partial p} \right), \qquad
\frac{\partial\log\det(\boldsymbol{A}^{-1})}{\partial p} = -
\mathrm{tr}\\\left( \boldsymbol{A}^{-1}
\frac{\partial\boldsymbol{A}}{\partial p} \right). \tag{2}

For a residual quadratic form with symmetric precision
\boldsymbol{\Theta} and mean map linear in p,

\frac{\partial}{\partial p} \bigl(
(\boldsymbol{y}-\boldsymbol{x}p)^{\top} \boldsymbol{\Theta}
(\boldsymbol{y}-\boldsymbol{x}p) \bigr) = -2\\ \boldsymbol{x}^{\top}
\boldsymbol{\Theta} (\boldsymbol{y}-\boldsymbol{x}p). \tag{3}

*Remark 2* (Second-order matrix calculus). Differentiating the inverse
once more gives

\frac{\partial^{2}\boldsymbol{A}^{-1}}{\partial p\_{i}\partial p\_{j}} =
\boldsymbol{A}^{-1} \Biggl( \frac{\partial\boldsymbol{A}}{\partial
p\_{i}} \boldsymbol{A}^{-1} \frac{\partial\boldsymbol{A}}{\partial
p\_{j}} - \frac{\partial^{2}\boldsymbol{A}}{\partial p\_{i}\partial
p\_{j}} + \frac{\partial\boldsymbol{A}}{\partial p\_{j}}
\boldsymbol{A}^{-1} \frac{\partial\boldsymbol{A}}{\partial p\_{i}}
\Biggr) \boldsymbol{A}^{-1}, \tag{4}

and the trace is linear, \partial\mathrm{tr}(\boldsymbol{A})/\partial
p\_{i} =\mathrm{tr}(\partial\boldsymbol{A}/\partial p\_{i}). Combining
[Eq. 1](#eq-matrix-first-order) with that linearity yields

\frac{\partial^{2}\log\det(\boldsymbol{A}^{-1})}{\partial p\_{i}\partial
p\_{j}} = - \mathrm{tr}\\\left\[ \boldsymbol{A}^{-1}
\frac{\partial^{2}\boldsymbol{A}}{\partial p\_{i}\partial p\_{j}}
\right\] + \mathrm{tr}\\\left\[ \left( \boldsymbol{A}^{-1}
\frac{\partial\boldsymbol{A}}{\partial p\_{i}} \right) \left(
\boldsymbol{A}^{-1} \frac{\partial\boldsymbol{A}}{\partial p\_{j}}
\right) \right\]. \tag{5}

> **Tip 2: Inner products in the score**
>
> Several terms are \boldsymbol{\Theta}-inner products of signature
> columns, \boldsymbol{\mu}\_{\cdot
> j}^{\top}\boldsymbol{\Theta}\boldsymbol{\mu}\_{\cdot k}. Package
> helper
> [`.inner_product()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-inner_product.md)
> evaluates that bilinear form; see also the squared Mahalanobis
> distance helper in the same family.

## Unconstrained score and Hessian

Conditional on the reference moments
\boldsymbol{\zeta}=(\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_{j}\\), the
bulk follows \boldsymbol{y}\mid\boldsymbol{\zeta},\boldsymbol{p}
\sim\mathcal{N}\_{G}\bigl(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p})\bigr). Up to an additive constant
independent of \boldsymbol{p},

\ell\_{\boldsymbol{y}\mid\boldsymbol{\zeta}}(\boldsymbol{p}) =
\tfrac{1}{2}\log\det\boldsymbol{\Theta}(\boldsymbol{p}) - \tfrac{1}{2}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta}(\boldsymbol{p})
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}). \tag{6}

> **Important 3: The determinant term carries a factor \tfrac{1}{2}**
>
> Both terms of [Eq. 6](#eq-loglik) are halved, because
> \log\det\boldsymbol{\Theta}=-\log\det\boldsymbol{\Sigma} enters the
> Gaussian log-density as
> -\tfrac{1}{2}\log\det\boldsymbol{\Sigma}(\boldsymbol{p}). DeCovarT
> releases before 2.3.0 implemented
> -\log\det\boldsymbol{\Sigma}(\boldsymbol{p}), which doubled the
> determinant contribution. The objective was then *not* the stated
> Gaussian model, and it disagreed with the expected Fisher information
> of [Sec. 4](#sec-fisher-wald), which had always been the one implied
> by the correct density. With the factor restored,
> [`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
> matches `mvtnorm::dmvnorm(..., log = TRUE)` up to the constant
> \tfrac{G}{2}\log(2\pi), and
> \mathbb{E}\[-\mathbf{H}\]=I(\boldsymbol{p}) holds exactly (see
> [Note 7](#nte-obs-vs-exp)). Point estimates move only slightly,
> because the two terms trade off differently, but likelihood
> **values**, AIC, likelihood-ratio statistics, and any curvature-based
> standard error do change.

> **Warning 4: What “unconstrained” means here**
>
> [Sec. 2](#sec-unconstrained) treats each p\_{j} as an ambient
> coordinate in (0,1) and ignores \sum\_{j}p\_{j}=1. The resulting
> gradient and Hessian are the building blocks of the ILR chain rule in
> [Sec. 3](#sec-constrained). They are **not** the coordinates used by
> Marquardt–Levenberg or Newton–Raphson in the package (those optimise
> in \boldsymbol{z}-space). The box-constrained L-BFGS-B path is the
> exception: it works directly in \boldsymbol{p} with \[0,1\]^{J} boxes
> and does not enforce the simplex exactly.

Applying [Eq. 2](#eq-logdet) and [Eq. 3](#eq-quad-form-deriv) with
\partial\boldsymbol{\Sigma}/\partial
p\_{j}=2p\_{j}\boldsymbol{\Sigma}\_{j} gives the ambient score. Colour
coding matches the manuscript: purple determinant / precision terms,
blue mean residual, orange covariance quadratic.

**Theorem 1 (Unconstrained gradient)** \begin{aligned}
\frac{\partial\ell}{\partial p\_{j}} &=
\underbrace{-p\_{j}\\\mathrm{tr}\bigl(
\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{j} \bigr)}\_{\text{purple}} +
\underbrace{ (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\mu}\_{\cdot j} }\_{\text{blue}} \\
&\quad+ \underbrace{ p\_{j}\\
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{j} \boldsymbol{\Theta}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}) }\_{\text{orange}}.
\end{aligned} \tag{7}

Differentiating [Eq. 7](#eq-grad-unconstrained) once more yields the
ambient Hessian \mathbf{H}\in\mathcal{M}\_{J\times J}. Diagonal and
off-diagonal blocks keep the same colour pairing.

**Theorem 2 (Unconstrained Hessian)** For i\in\\1,\ldots,J\\,

\begin{aligned} \mathbf{H}\_{ii} &= \frac{\partial^{2}\ell}{\partial
p\_{i}^{2}} = \underbrace{
-\\\mathrm{tr}(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{i}) +
2p\_{i}^{2}\\\mathrm{tr}\bigl(
(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{i})^{2} \bigr)
}\_{\text{purple}} \\ &\quad \underbrace{ -2p\_{i}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{i} \boldsymbol{\Theta}
\boldsymbol{\mu}\_{\cdot i} - \boldsymbol{\mu}\_{\cdot i}^{\top}
\boldsymbol{\Theta} \boldsymbol{\mu}\_{\cdot i} }\_{\text{blue}} \\
&\quad \underbrace{ -2p\_{i}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{i} \boldsymbol{\Theta}
\boldsymbol{\mu}\_{\cdot i} -
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \bigl(
4p\_{i}^{2}\boldsymbol{\Sigma}\_{i}\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{i} -
\boldsymbol{\Sigma}\_{i} \bigr) \boldsymbol{\Theta}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}) }\_{\text{orange}},
\end{aligned} \tag{8}

and for i\neq j,

\begin{aligned} \mathbf{H}\_{ij} &= \frac{\partial^{2}\ell}{\partial
p\_{i}\partial p\_{j}} = \underbrace{ 2p\_{i}p\_{j}\\ \mathrm{tr}\bigl(
\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{j}
\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{i} \bigr) }\_{\text{purple}} \\
&\quad \underbrace{ -2p\_{i}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{i} \boldsymbol{\Theta}
\boldsymbol{\mu}\_{\cdot j} - \boldsymbol{\mu}\_{\cdot i}^{\top}
\boldsymbol{\Theta} \boldsymbol{\mu}\_{\cdot j} }\_{\text{blue}} \\
&\quad \underbrace{ -2p\_{j}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{j} \boldsymbol{\Theta}
\boldsymbol{\mu}\_{\cdot i} - 4p\_{i}p\_{j}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta} \boldsymbol{\Sigma}\_{i} \boldsymbol{\Theta}
\boldsymbol{\Sigma}\_{j} \boldsymbol{\Theta}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}) }\_{\text{orange}}.
\end{aligned} \tag{9}

The stationary equation \nabla\ell(\boldsymbol{p})=\boldsymbol{0} has no
closed root, so DeCovarT uses iterative solvers. Analytic gradients and
Hessians are checked against Richardson extrapolation in `numDeriv`
([Gilbert and Varadhan 2019](#ref-R-numDeriv)).

## Constrained optimisation on the simplex

Constrained optimisation may eliminate constraints by reparametrisation,
enforce them softly with penalties or barriers, or enforce them exactly
with projections / KKT methods. DeCovarT uses a C^{2} diffeomorphism
from free Euclidean coordinates \boldsymbol{z}\in\mathbb{R}^{J-1} onto
the open simplex, so every Newton step stays feasible by construction.

### Log-sum-exp and softmax

The softmax of scores \boldsymbol{z}\in\mathbb{R}^{K} is evaluated
through the log-sum-exp (LSE) to avoid overflow,

s\_{i} = \frac{e^{z\_{i}}}{\sum\_{k}e^{z\_{k}}} =
e^{z\_{i}-\operatorname{LSE}(\boldsymbol{z})}, \qquad
\operatorname{LSE}(\boldsymbol{z}) = m+\log\sum\_{k}e^{z\_{k}-m}, \quad
m=\max\_{k}z\_{k}. \tag{10}

Its Jacobian is
\mathbf{J}\_{\mathrm{softmax}}=\operatorname{diag}(\boldsymbol{s})
-\boldsymbol{s}\boldsymbol{s}^{\top}, which is also
\nabla^{2}\operatorname{LSE}(\boldsymbol{z}). The second derivative is
the third-order tensor

H\_{ijk} = s\_{i}\bigl\[ (\delta\_{ij}-s\_{j})(\delta\_{ik}-s\_{k}) -
s\_{j}(\delta\_{jk}-s\_{k}) \bigr\]. \tag{11}

### Isometric log-ratio map

Choose a Helmert sub-matrix \mathbf{V}\in\mathbb{R}^{J\times(J-1)} with

\mathbf{V}^{\mathsf{T}}\mathbf{V}=\mathbf{I}\_{J-1}, \qquad
\mathbf{V}^{\mathsf{T}}\mathbf{1}=\mathbf{0}. \tag{12}

ILR coordinates and their inverse are ([Pawlowsky-Glahn and Buccianti
2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011))

\boldsymbol{z} = \mathbf{V}^{\mathsf{T}}\log\boldsymbol{p} =
\mathbf{V}^{\mathsf{T}}\operatorname{clr}(\boldsymbol{p}), \qquad
\boldsymbol{p} = \operatorname{softmax}(\mathbf{V}\boldsymbol{z}) =
\mathcal{C}\bigl\\\exp(\mathbf{V}\boldsymbol{z})\bigr\\, \tag{13}

where \mathcal{C} denotes closure to unit sum. No cell type is pinned as
a reference (contrast [Sec. 8](#sec-alr)). Package helpers
[`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md),
[`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md)
and
[`isometric_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_log_ratio.md)
implement [Eq. 12](#eq-helmert)–[Eq. 13](#eq-ilr-maps) with a
log-sum-exp softmax, matching the
[`compositions::ilr()`](https://rdrr.io/pkg/compositions/man/ilr.html) /
`ilrInv()` convention on a Helmert basis ([van den Boogaart et al.
2025](#ref-R-compositions); [Egozcue et al.
2003](#ref-egozcueIsometricLogratioTransformations2003)).

**Proposition 1 (ILR is a C^{2} diffeomorphism)**
\boldsymbol{\psi}:\mathbb{R}^{J-1}\to(0,1)^{J} with \sum\_{j}p\_{j}=1 is
bijective and twice continuously differentiable, so every critical point
in \boldsymbol{z}-space maps to a unique open simplex composition and
conversely. Any other valid ILR basis is \mathbf{V}^{\star}=\mathbf{V}Q
with Q^{\mathsf{T}}Q=\mathbf{I}, which rotates coordinates as
\boldsymbol{z}^{\star}=Q^{\mathsf{T}} \boldsymbol{z} and leaves
eigenvalues, condition numbers and traces of properly transformed
quadratic forms unchanged.

> **Note 5: Why ILR is the default chart**
>
> Changing an ALR denominator changes the coordinates, not the
> composition. If derivatives, Fisher information and the delta method
> are transformed correctly, \hat{\boldsymbol{p}},
> \ell(\hat{\boldsymbol{p}}) and
> \operatorname{Var}(\hat{\boldsymbol{p}}) must agree after mapping back
> to the simplex. A material change in \hat{\boldsymbol{p}} when the
> reference cell is permuted is therefore numerical or optimisation
> sensitivity, not a statistical reference dependence.
>
> ILR still improves the chart because distinct orthonormal bases differ
> by an orthogonal rotation, whereas changing an ALR denominator is a
> generally non-orthogonal reparametrisation. Euclidean gradient norms,
> Hessian condition numbers and Newton step geometry are then natural up
> to rotation. A fixed Helmert basis introduces an ordering convention,
> not a denominator. The ALR formulae are kept in [Sec. 8](#sec-alr) for
> comparison and for the structural-zero discussion
> ([Sec. 8.1](#sec-structural-zeros)).

### Efficient derivatives of \boldsymbol{\psi}

Write \mathbf{S}(\boldsymbol{p})=\operatorname{diag}(\boldsymbol{p})
-\boldsymbol{p}\boldsymbol{p}^{\mathsf{T}}. The Jacobian is

\mathbf{J}\_{\boldsymbol{\psi}}(\boldsymbol{z}) =
\frac{\partial\boldsymbol{p}}{\partial\boldsymbol{z}^{\mathsf{T}}} =
\mathbf{S}(\boldsymbol{p})\mathbf{V}. \tag{14}

[`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md)
returns this J\times(J-1) matrix. Let
\bar{\boldsymbol{v}}=\mathbf{V}^{\mathsf{T}}\boldsymbol{p},
\boldsymbol{q}\_{i}=\mathbf{V}\_{i\cdot}-\bar{\boldsymbol{v}} and
\mathbf{C}\_{V}=\mathbf{V}^{\mathsf{T}}\mathbf{S}(\boldsymbol{p})
\mathbf{V}. Slice i of the map Hessian is

\mathbf{H}\_{\psi\_{i}} = p\_{i}\bigl(
\boldsymbol{q}\_{i}\boldsymbol{q}\_{i}^{\mathsf{T}} -\mathbf{C}\_{V}
\bigr). \tag{15}

[`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md)
stores the tensor of shape (J-1)\times(J-1)\times J.

### Change of variables: ambient \boldsymbol{p} to ILR \boldsymbol{z}

Let \ell\_{\boldsymbol{z}}(\boldsymbol{z})
=\ell\_{\boldsymbol{y}\mid\boldsymbol{\zeta}}
\bigl(\boldsymbol{\psi}(\boldsymbol{z})\bigr). The chain rule maps the
unconstrained derivatives of [Sec. 2](#sec-unconstrained) into free
Euclidean coordinates.

**Theorem 3 (Constrained score and Hessian)**
\nabla\_{\boldsymbol{z}}\ell =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}\nabla\_{\boldsymbol{p}}\ell
=
\mathbf{V}^{\mathsf{T}}\mathbf{S}(\boldsymbol{p})\nabla\_{\boldsymbol{p}}\ell,
\tag{16}

\mathbf{H}\_{\boldsymbol{z}} =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\mathbf{H}\_{\boldsymbol{p}} \mathbf{J}\_{\boldsymbol{\psi}} +
\sum\_{i=1}^{J} \frac{\partial\ell}{\partial p\_{i}}
\mathbf{H}\_{\psi\_{i}}. \tag{17}

At an interior KKT point every active score equals the same multiplier,
\partial\ell/\partial p\_{i}=\lambda. Because
\sum\_{i}p\_{i}(\boldsymbol{z})=1 for every \boldsymbol{z},
\sum\_{i}\mathbf{H}\_{\psi\_{i}}=\mathbf{0}, so the second summand of
[Eq. 17](#eq-chain-second) vanishes and \mathbf{H}\_{\boldsymbol{z}}
=\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\mathbf{H}\_{\boldsymbol{p}} \mathbf{J}\_{\boldsymbol{\psi}}. Away from
stationarity that contraction must be kept.

> **Important 6: Implementation map**
>
> | Object | Package function |
> |:---|:---|
> | Ambient score / Hessian | [`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md), [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md) |
> | Map Jacobian / Hessian | [`jacobian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_isometric_logistic.md), [`hessian_isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_isometric_logistic.md) |
> | ILR score / Hessian | [`gradient_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_constrained.md), [`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md) |
>
> [Eq. 17](#eq-chain-second) is the sum of a covariant pull-back of
> \mathbf{H} and a contraction of the ambient score against the map
> Hessian. Omitting the second summand is a common bug when porting
> first-order code to Newton methods.

CLR coordinates and Bayesian formulations sit on the same footing
([Pawlowsky-Glahn and Buccianti
2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)): a Dirichlet
or logistic-normal prior on \boldsymbol{p} pairs naturally with the
Gaussian bulk likelihood (the experimental MAP `CTS` path in
`R/03_04_DeCovarT_estimate_CTS_MAP_Bayesian.R` is one starting point).
Regardless of coordinate system, the workflow is the same: optimise in
an unconstrained parameterisation, then map back to the simplex for
interpretation.

## Expected Fisher information and Wald inference

Point estimation alone does not report uncertainty. Under a regular
large-sample regime the MLE is asymptotically normal with covariance
given by the inverse Fisher information. DeCovarT builds Wald standard
errors from the *expected* Fisher information of the multivariate-normal
convolution, then maps that information through the ILR chart.

### Unconstrained expected Fisher information

For
\boldsymbol{y}\sim\mathcal{N}\_{G}\bigl(\boldsymbol{\mu}\boldsymbol{p},
\boldsymbol{\Sigma}(\boldsymbol{p})\bigr) the expected Fisher
information has entries

I(\boldsymbol{p})\_{jk} = \boldsymbol{\mu}\_{\cdot j}^{\top}
\boldsymbol{\Theta}(\boldsymbol{p}) \boldsymbol{\mu}\_{\cdot k} +
2p\_{j}p\_{k}\\ \mathrm{tr}\bigl(
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_{j}
\boldsymbol{\Theta}(\boldsymbol{p})\boldsymbol{\Sigma}\_{k} \bigr).
\tag{18}

The first summand is the mean contribution (a \boldsymbol{\Theta}-inner
product of signature columns); the second is the covariance contribution
of the quadratic map
\boldsymbol{p}\mapsto\boldsymbol{\Sigma}(\boldsymbol{p}). See the
multivariate-normal formula on [Wikipedia: Fisher
information](https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution).
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
implements [Eq. 18](#eq-fisher-p).

> **Note 7: Observed versus expected information**
>
> The analytic Hessian of [Sec. 2](#sec-unconstrained) is the *observed*
> information (up to sign). Wald intervals in `decovart_fit` use the
> *expected* information [Eq. 18](#eq-fisher-p), which does not depend
> on the residual \boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p} and is
> positive definite under the usual full-rank assumptions on
> (\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_{j}\\).
>
> The two are consistent by construction. Substituting
> \mathbb{E}\[\boldsymbol{r}\]=\boldsymbol{0} and
> \mathbb{E}\[\boldsymbol{r}\boldsymbol{r}^{\top}\]
> =\boldsymbol{\Sigma}(\boldsymbol{p}) into [Eq. 8](#eq-hess-diag) and
> [Eq. 9](#eq-hess-off) cancels the determinant traces against the
> residual traces and returns [Eq. 18](#eq-fisher-p) exactly, i.e.
> \mathbb{E}\[-\mathbf{H}\]=I(\boldsymbol{p}). That identity only holds
> with the factor \tfrac{1}{2} of [Note 3](#nte-halfdet) in place.

### Constrained case: ILR delta method

Fisher information transforms as a covariant quadratic form under the
ILR chart. With \mathbf{J}\_{\boldsymbol{\psi}}
=\mathbf{S}(\boldsymbol{p})\mathbf{V},

I\_{\boldsymbol{z}} = \mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
I(\boldsymbol{p}) \mathbf{J}\_{\boldsymbol{\psi}} =
\mathbf{V}^{\mathsf{T}} \mathbf{S}(\boldsymbol{p}) I(\boldsymbol{p})
\mathbf{S}(\boldsymbol{p}) \mathbf{V}. \tag{19}

The mean/covariance split of [Eq. 18](#eq-fisher-p) pulls back
separately, I\_{\boldsymbol{z}}=I\_{\boldsymbol{z},\mathrm{mean}}
+I\_{\boldsymbol{z},\mathrm{cov}}. Under regularity,
\hat{\boldsymbol{z}} \overset{a}{\sim}
\mathcal{N}(\boldsymbol{z}\_{0},I\_{\boldsymbol{z}}^{-1}) ([Vaart
2000](#ref-vaartAsymptoticStatistics2000)). The first-order delta method
([Oehlert 1992](#ref-oehlertNoteDeltaMethod1992)) then returns the
simplex covariance used by
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),

\mathrm{Var}(\hat{\boldsymbol{p}}) \approx
\mathbf{J}\_{\boldsymbol{\psi}} I\_{\boldsymbol{z}}^{-1}
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}. \tag{20}

Diagonal square roots are asymptotic standard errors; Wald intervals at
level 1-\alpha are \hat{p}\_{j}\pm z\_{1-\alpha/2}\\\mathrm{SE}\_{j}
with z\_{q}=\Phi^{-1}(q), as in
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
implements [Eq. 19](#eq-fisher-z)–[Eq. 20](#eq-delta-p). Orthogonal
rotations of \mathbf{V} leave [Eq. 20](#eq-delta-p) unchanged. The
construction is undefined on the simplex boundary (the log-ratio chart
blows up); the helper then returns `NA` with a warning.
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md)
is the ALR analogue of [Sec. 8](#sec-alr): the reconstructed simplex
covariance must agree with [Eq. 20](#eq-delta-p) up to numerical error.

> **Warning 8: Wald is the cheapest, not the safest, interval**
>
> [Eq. 20](#eq-delta-p) linearises the ILR chart at
> \hat{\boldsymbol{p}}, so it is neither invariant to a non-orthogonal
> reparametrisation nor confined to \[0,1\], and it degenerates entirely
> when a proportion reaches a simplex face. The companion vignette [MLE
> properties and asymptotic
> inference](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.md)
> derives the alternatives: profile likelihood-ratio intervals, the
> chi-bar-square calibration of Chernoff and of Self and Liang for nulls
> on the boundary, a parametric bootstrap, and a reference-sample
> bootstrap of donors or cells (or a Dirichlet composition sweep). It
> also shows why these limits need replicate samples that share one
> composition. Shuffling genes or cell-type names of an averaged
> signature is not a bootstrap: the maximised likelihood is equivariant
> under relabelling.

> **Numerical speed-ups and solver safeguards**
>
> The analytic maps above are only half of a usable optimiser. Practical
> bottlenecks that showed up on the hybrid scenario in [Variance-driven
> hybrid
> scenario](https://bastienchassagnol.github.io/DeCovarT/articles/fig03-variance-driven.html#sec-hybrid-deconvolution)
> live in `R/03_03_DeCovarT_estimate_ratios_frequentist.R`.
>
> ------------------------------------------------------------------------
>
> #### Cache a Cholesky factorisation of \boldsymbol{\Sigma}(\boldsymbol{p})
>
> Within one iteration, objective / gradient / Hessian callbacks hit the
> same trial \boldsymbol{p}. The helper
> [`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
> caches a single Cholesky factor \boldsymbol{R} with
> \boldsymbol{\Sigma}(\boldsymbol{p})=\boldsymbol{R}^{\mathsf{T}}\boldsymbol{R}
> and returns \log\det\boldsymbol{\Sigma}(\boldsymbol{p})=2\sum_g\log
> R\_{gg} together with the precision
> \boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{R}^{-1}\boldsymbol{R}^{-\mathsf{T}}
> (via `chol2inv`) without repeating an O(G^{3}) factorisation.
>
> [`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
> then evaluates the Mahalanobis term exactly as
> [`mvtnorm::dmvnorm`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
> ([Genz et al. 2026](#ref-R-mvtnorm)): it solves
> \boldsymbol{R}^{\mathsf{T}}\boldsymbol{z}
> =\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p} by
> `backsolve(..., transpose = TRUE)` and takes
> \lVert\boldsymbol{z}\rVert^{2}, rather than forming
> \boldsymbol{\Theta}(\boldsymbol{p}) and a dense quadratic form. The
> two routes agree to machine precision; the backsolve path is the one
> used for the objective, while the cached inverse is retained for the
> analytic score and Hessian, which need
> \boldsymbol{\Theta}(\boldsymbol{p}) in several trace and inner-product
> terms. A QR factorisation of the covariance itself would recover the
> same log-determinant and quadratic form at a larger O(G^{3}) constant
> and is not used: Cholesky is the natural factorisation of a symmetric
> positive-definite matrix.
>
> The implemented log-density omits the additive -\tfrac{G}{2}\log(2\pi)
> of the full Gaussian, which does not depend on \boldsymbol{p} and
> therefore cannot change the MLE, the score, or the Hessian. Tests
> check that subtracting that constant recovers
> `dmvnorm(..., log = TRUE)` exactly.
>
> ------------------------------------------------------------------------
>
> #### Guard the box-constrained L-BFGS-B path
>
> [`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
> optimises in \boldsymbol{p} with boxes \[0,1\]^{J} that do not enforce
> \mathbf{1}^{\top}\boldsymbol{p}=1. Near \sum\_{j}p\_{j}\approx 0 (or
> on a failed [`chol()`](https://rdrr.io/r/base/chol.html)), safeguarded
> wrappers return a finite penalty and a zero gradient rather than
> aborting [`optim()`](https://rdrr.io/r/stats/optim.html).
>
> ------------------------------------------------------------------------
>
> #### `marqLevAlg` Hessian sign under `minimize = FALSE`
>
> [`marqLevAlg`](https://cran.r-project.org/package=marqLevAlg)
> ([Philipps et al. 2023](#ref-R-marqLevAlg)) implements
> Marquardt–Levenberg with the relative-distance-to-minimum (RDM) rule
> of Commenges et al. ([Commenges et al.
> 2006](#ref-commengesNewtonLikeAlgorithmLikelihood2006)). With
> `minimize = FALSE`, the package negates only `fn` and `gr`; the
> analytic `hess` is passed through unchanged while internal Cholesky
> routines assume a **positive-definite** matrix at the optimum. When
> maximising a log-likelihood, pass
> `hess = function(...) -hessian_loglik_constrained(...)`. Reported
> upstream as
> [VivianePhilipps/marqLevAlgParallel#3](https://github.com/VivianePhilipps/marqLevAlgParallel/issues/3).

## Model objects and S3 accessors

[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
returns a `decovart_fit` object with standard S3 methods for
coefficients, fitted values, residuals, and (under regularity) a
delta-method covariance. The chunk below mirrors the appendix
equivariance setup:

``` r

library(DeCovarT)
genes <- paste0("g", 1:3)
cts <- paste0("ct", 1:2)
mu <- matrix(
  c(20, 40, 25, 35, 30, 22),
  nrow = 3,
  dimnames = list(genes, cts)
)
Sigma <- array(0, dim = c(3, 3, 2), dimnames = list(genes, genes, cts))
Sigma[, , 1] <- diag(c(1.0, 1.2, 0.8))
Sigma[, , 2] <- diag(c(1.5, 0.9, 1.1))
p_true <- c(0.65, 0.35)
y <- matrix(
  drop(mu %*% p_true + rnorm(3, sd = 0.05)),
  ncol = 1L,
  dimnames = list(genes, "s1")
)

fit <- fit_decovart(
  mu,
  y,
  Sigma = Sigma,
  method = "Marquardt-Levenberg",
  itmax = 80
)
coef(fit)
nobs(fit)
c(
  n_genes = attr(nobs(fit), "n_genes"),
  n_celltypes = attr(nobs(fit), "n_celltypes")
)
fitted(fit)
residuals(fit)
vcov(fit)
```

[`formula()`](https://rdrr.io/r/stats/formula.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) are not implemented:
DeCovarT does not forecast bulk expression.

## Numerical consistency of the likelihood and its derivatives

Analytic formulae can be differentiated correctly *and* still implement
the wrong objective. The checks below separate those two questions. They
are the same identities that `tests/testthat/test-03_03_DeCovarT.R` and
`tests/testthat/test-03_06_decovart_inference.R` enforce on every
package build; the chunks here are a readable, small-(G) illustration.

### Density against `mvtnorm`

The first check is against an independent implementation of the
multivariate Gaussian log-density, not against our own inverse. With the
\tfrac{1}{2}\log\det factor restored, DeCovarT and
[`mvtnorm::dmvnorm`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html) ([Genz
et al. 2026](#ref-R-mvtnorm)) differ by exactly \tfrac{G}{2}\log(2\pi)
(DeCovarT omits the negative constant):

``` r

n_genes <- length(y_chk)
ours <- loglik_multivariate(p_chk, y_chk, mu_chk, Sigma_chk)
reference <- if (requireNamespace("mvtnorm", quietly = TRUE)) {
  as.numeric(mvtnorm::dmvnorm(
    y_chk,
    mean = drop(mu_chk %*% p_chk),
    sigma = cov_chk,
    log = TRUE
  ))
} else {
  residual <- y_chk - drop(mu_chk %*% p_chk)
  as.numeric(
    -0.5 * determinant(cov_chk, logarithm = TRUE)$modulus -
      0.5 * drop(crossprod(residual, solve(cov_chk, residual))) -
      0.5 * n_genes * log(2 * pi)
  )
}
c(
  loglik_multivariate = ours,
  dmvnorm = reference,
  difference = ours - 0.5 * n_genes * log(2 * pi) - reference
)
```

A finite-difference check of the *same* objective cannot detect a shared
factor-of-two error in the log-likelihood and its derivatives. That is
why this density comparison is independent of `numDeriv`.

### Richardson extrapolation of the score and Hessian

Once the objective is the intended Gaussian, the analytic score and
Hessian are checked against Richardson extrapolation in `numDeriv`
([Gilbert and Varadhan 2019](#ref-R-numDeriv)). Richardson’s tableau
extrapolates central differences to cancel successive even powers of the
step size, which is more accurate than a single two-point stencil at the
same O(\varepsilon^{-1}) cost in function evaluations.

``` r

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  numderiv_report <- c(gradient_rel = NA_real_, hessian_rel = NA_real_)
} else {
  grad_n <- numDeriv::grad(
    loglik_multivariate,
    p_chk,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = y_chk,
    mean_signature_matrix = mu_chk,
    Sigma = Sigma_chk
  )
  hess_n <- numDeriv::hessian(
    loglik_multivariate,
    p_chk,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 4),
    y = y_chk,
    mean_signature_matrix = mu_chk,
    Sigma = Sigma_chk
  )
  grad_a <- gradient_loglik_unconstrained(
    p_chk, y_chk, mu_chk, Sigma_chk
  )
  hess_a <- hessian_loglik_unconstrained(
    p_chk, y_chk, mu_chk, Sigma_chk
  )
  numderiv_report <- c(
    gradient_rel = max(abs(grad_a - grad_n)) / max(abs(grad_n)),
    hessian_rel = max(abs(hess_a - hess_n)) / max(abs(hess_n))
  )
}
numderiv_report
```

Relative discrepancies of order (10^{-8}) are typical on this three-gene
toy; that is the same tolerance used in the test suite.

### Expected information versus Monte Carlo (\[-\])

Matching `numDeriv` shows that the derivatives belong to the implemented
objective. Matching the *expected* Fisher information
[Eq. 18](#eq-fisher-p) to a Monte Carlo average of (-(; )) shows that
the objective is the Gaussian whose information matrix we use for Wald
intervals. Draw (^{(b)}\_G(, ())) at the *true* () (the parametric
bootstrap of
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
used here as a diagnostic rather than as an inferential procedure):

``` r

set.seed(2)
n_mc <- 40L
info <- expected_fisher_unconstrained(p_chk, mu_chk, Sigma_chk)
neg_hess <- array(NA_real_, c(3L, 3L, n_mc))
for (b in seq_len(n_mc)) {
  y_b <- MASS::mvrnorm(
    n = 1L,
    mu = drop(mu_chk %*% p_chk),
    Sigma = cov_chk
  )
  neg_hess[, , b] <- -hessian_loglik_unconstrained(
    p_chk, y_b, mu_chk, Sigma_chk
  )
}
mc_info <- apply(neg_hess, c(1L, 2L), mean)
c(
  rel_frobenius = sqrt(sum((mc_info - info)^2)) / sqrt(sum(info^2))
)
```

The same identity is what closed the factor-of-two inconsistency: with
the old (-) objective, (\[-\]) could not equal (I()). A few dozen draws
already bring the relative Frobenius gap below (10^{-2}) on this toy;
the test suite uses a tighter Monte Carlo.

A related, weaker check compares the empirical covariance of
parametric-bootstrap MLEs to the Cramér–Rao map
[`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md).
That comparison is noisy at small replication (see [MLE
properties](https://bastienchassagnol.github.io/DeCovarT/articles/theory-DeCovarT-MLE-properties.html#sec-replication)
for why one bulk column is not an asymptotic sample) and is therefore
left to the inference vignette and to
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md).

## Structure-aware covariance backends

Every solver call factorises \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^{2}\boldsymbol{\Sigma}\_j. The dense Cholesky in
[`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
is the universal default and the only method that needs no assumptions.
When the covariance inherits structure from the underlying gene network,
the same log-determinant and solves can be obtained far more cheaply.
Exploiting band, sparse and diagonal-plus-low-rank structure in a
Newton-type solver is standard practice ([Boyd et al. 2004, sect.
9.7](#ref-boydConvexOptimization2004)); the backends below make that
structure explicit for DeCovarT.

[`new_decovart_covariance()`](https://bastienchassagnol.github.io/DeCovarT/reference/new_decovart_covariance.md)
wraps the covariance array with a declared `structure`, and the
operators
[`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md),
[`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md),
[`sigma_quadform()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_quadform.md)
and
[`sigma_trace_precision_times()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_trace_precision_times.md)
return exactly the quantities the objective, score and Hessian need.
They expose *operations*, never a materialised dense precision: a sparse
or Woodbury solve must not be silently undone by a dense `chol2inv`.

### Decision tree

Structure is supplied explicitly, never inferred from floating-point
near-zeros, so a numerical near-zero is never mistaken for an exact
conditional independence.
[`covariance_structure_from_graph_model()`](https://bastienchassagnol.github.io/DeCovarT/reference/covariance_structure_from_graph_model.md)
maps a network topology to the recommended backend. A sparse *precision*
\boldsymbol{\Omega}\_j does not imply a sparse *covariance*
\boldsymbol{\Sigma}\_j; only structure that survives in
\boldsymbol{\Sigma}(\boldsymbol{p}) itself accelerates the likelihood.

``` mermaid
%%{init: {"theme": "sandstone"}}%%
flowchart TD
  A["Sigma(p) = sum_j p_j^2 Sigma_j"] --> B{"Declared covariance structure?"}
  B -->|"none (default)"| D["Dense Cholesky\nO(G^3)"]
  B -->|"disconnected modules\n(stochastic block model)"| E["Block Cholesky\nO(sum_b G_b^3)"]
  B -->|"ordered local dependence\n(band / AR)"| F["Banded Cholesky\nO(G b^2)"]
  B -->|"fixed sparse pattern\n(Erdos-Renyi, scale-free, small-world)"| G["Sparse Cholesky\nsymbolic ordering computed once"]
  B -->|"shared programs / bridge hubs\n(pathways, housekeeping genes)"| H["Woodbury / Schur\nO(G r^2)"]
```

Figure 1: Choosing a covariance backend; dense Cholesky is the default
and structured backends are exact refinements.

Each refinement is exact, not approximate:

- **Block** (disconnected gene modules, e.g. a stochastic block model
  with no between-block edges ([Holland et al.
  1983](#ref-hollandStochasticBlockmodelsFirst1983))).
  \boldsymbol{\Sigma}(\boldsymbol{p}) is block-diagonal, so the
  log-likelihood separates over blocks and one Cholesky per block
  replaces one dense Cholesky.
- **Band / AR.** Only genes within bandwidth b covary; a banded Cholesky
  costs O(G b^{2}) ([Golub and Van Loan
  2013](#ref-golubMatrixComputations2013)).
- **Sparse** (Erdős–Rényi, scale-free ([Barabási and Albert
  1999](#ref-barabasiEmergenceScalingRandom1999)), small-world ([Watts
  and Strogatz 1998](#ref-wattsCollectiveDynamicsSmallworld1998))). The
  fill-reducing ordering is a *symbolic factorisation* computed once and
  reused for every numeric refactor as \boldsymbol{p} changes; only the
  numerical Cholesky repeats ([Chen et al.
  2008](#ref-chenAlgorithm887CHOLMOD2008)).
- **Diagonal + low rank.** Genes interact through r\ll G shared
  regulatory programs,
  \boldsymbol{\Sigma}\_j=\mathbf{D}\_j+\mathbf{U}\mathbf{C}\_j
  \mathbf{U}^{\mathsf{T}} — the factor-analysis form
  \boldsymbol{\Sigma}=\mathbf{W}\mathbf{W}^{\mathsf{T}}+\boldsymbol{\Psi}.
  The Woodbury identity and the matrix-determinant lemma ([Hager
  1989](#ref-hagerUpdatingInverseMatrix1989)) move the O(G^{3}) cost to
  the low rank r. A few bridge / housekeeping hub genes give the
  arrow-matrix special case handled by the same Schur-complement
  algebra.

``` r

# Declare structure explicitly, then use the operator interface.
cov <- new_decovart_covariance(Sigma, structure = "diag_lowrank",
  diagonal = D, loadings = U, core = C)
ld <- sigma_logdet(cov, p)            # log det Sigma(p)
z  <- sigma_solve(cov, p, residual)   # Sigma(p)^{-1} r, no dense inverse
q  <- sigma_quadform(cov, p, residual)
```

### Computational complexity

For a dense symmetric positive-definite G\times G matrix the Cholesky
factorisation costs about G^{3}/3 flops before the cheaper triangular
solves ([Golub and Van Loan 2013](#ref-golubMatrixComputations2013)).
The structured backends reduce that leading term:

| Structure | When it applies | Factorisation cost | Backend argument |
|:---|:---|:---|:---|
| Dense | universal default | O(G^{3}) | (default) |
| Block | disconnected modules (SBM) | O(\sum_b G_b^{3}); \approx B^{2} saving for B equal blocks | `structure = "block"` |
| Band / AR | bandwidth b | O(G b^{2}) | `structure = "band"` |
| Sparse | fixed sparse pattern | pattern-dependent; symbolic ordering reused across \boldsymbol{p} | `structure = "sparse"` |
| Diagonal + low rank | r shared programs | O(G r^{2}) solve, O(r^{3}) core | `structure = "diag_lowrank"` |

Table 1: Leading factorisation cost per backend. G genes, J cell types,
B blocks of size G_b, bandwidth b, rank r\ll G.

These are worst-case flop counts and only rough predictors of wall-clock
time on cache-based hardware with tuned BLAS; the block and low-rank
savings (\approx B^{2} and O(G^{2})\\\to\\O(Gr^{2}) respectively) are
the most robust. The `"band"` and `"sparse"` backends use the imported
`Matrix` package (supernodal sparse Cholesky in CHOLMOD, ([Chen et al.
2008](#ref-chenAlgorithm887CHOLMOD2008))).
`scripts/supp_S2_covariance_inversion.R` verifies each backend against
the dense path to machine precision. Approximate sparse Cholesky of
otherwise dense kernel covariances is a separate research line ([Schäfer
et al. 2020](#ref-sch%C3%A4ferCompressionInversionApproximate2020)); the
`"sparse"` backend here is exact for a declared sparse pattern, not that
approximation.

## Appendix: additive log-ratio chart

Pin the last category as reference (\rho\_{J}\equiv 0). The *additive
logistic* transform \boldsymbol{\psi} (inverse ALR) and its inverse are
([Aitchison 1982](#ref-aitchisonStatisticalAnalysisCompositional1982))

p\_{j} = \frac{e^{\rho\_{j}}}{\sum\_{k\<J}e^{\rho\_{k}}+1}\\ (j\<J),
\qquad p\_{J} = \frac{1}{\sum\_{k\<J}e^{\rho\_{k}}+1}, \qquad \rho\_{j}
= \log\\\left(\frac{p\_{j}}{p\_{J}}\right). \tag{21}

Equivalently,
\boldsymbol{p}=\operatorname{softmax}(\rho\_{1},\ldots,\rho\_{J-1},0).
Package helpers
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
and
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_log_ratio.md)
implement [Eq. 21](#eq-alr-maps); they match
[`compositions::alrInv()`](https://rdrr.io/pkg/compositions/man/alr.html)
/ `alr()` with the last part as reference ([van den Boogaart et al.
2025](#ref-R-compositions)). The Jacobian is the softmax Jacobian with
the reference column dropped, \partial
p\_{i}/\partial\rho\_{a}=p\_{i}(\delta\_{ia}-p\_{a}) for a=1,\ldots,J-1
([`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md);
[`hessian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_additive_logistic.md)).

> **Note 9: Why ALR is not the default**
>
> ALR is a valid C^{2} diffeomorphism onto the open simplex, and the
> reconstructed \operatorname{Var}(\hat{\boldsymbol{p}}) of
> [`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md)
> must match
> [`vcov_ilr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_ilr_delta.md)
> when both charts are implemented correctly. The last cell type is
> nevertheless a privileged denominator: permuting which type occupies
> position J changes the Euclidean metric on coordinate space. ILR /
> Helmert ([Sec. 3.2](#sec-ilr)) removes that privilege up to orthogonal
> rotation, which is why the solvers and
> [`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
> use ILR. Keep ALR when comparing against a multinomial-logit
> reference-category parameterisation, or when checking that a change of
> denominator does not move \hat{\boldsymbol{p}}.

### Structural zeros and folded simplex models

Neither ALR nor ILR can represent an exact zero at a finite coordinate,
because both contain a logarithm. For a true simplex face the package
therefore does **not** send \log 0 through the chart. It fixes the
inactive set, optimises ILR coordinates on the remaining active face
([`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)),
and checks KKT inequalities for the zero coordinates
([`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md);
see also the MLE-properties vignette).

ALR remains in this appendix as the classical reference-category chart
([Aitchison 1982](#ref-aitchisonStatisticalAnalysisCompositional1982)).
Folded \alpha-models on the simplex ([Tsagris and Stewart
2020](#ref-tsagrisFoldedModelCompositional2020)) and nonparametric
regressions that admit exact zeros without a log-ratio chart ([Tsagris
et al. 2023](#ref-tsagrisFlexibleNonparametricRegression2023)) are
related *data-side* constructions. They are not used in the DeCovarT
convolution, whose mean map is still linear in \boldsymbol{p}. Rounded
zeros (detection limits) remain an imputation problem; true structural
zeros belong on a face and should be fitted as an active-set / KKT
problem, not by substituting a small positive constant into
[Eq. 13](#eq-ilr-maps).

## References

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Barabási, Albert-László, and Réka Albert. 1999. ‘Emergence of Scaling in
Random Networks’. *Science* 286.
<https://doi.org/10.1126/science.286.5439.509>.

Boyd, Stephen, Stephen P. Boyd, and Lieven Vandenberghe. 2004. *Convex
Optimization*. 1st edn. Cambridge University Press; Cambridge University
Press. <https://doi.org/10.1017/cbo9780511804441>.

Chen, Yanqing, Timothy A. Davis, William W. Hager, and Sivasankaran
Rajamanickam. 2008. ‘Algorithm 887: CHOLMOD, Supernodal Sparse Cholesky
Factorization and Update/Downdate’. *ACM Transactions on Mathematical
Software* 35 (3): 22:1–14. <https://doi.org/10.1145/1391989.1391995>.

Commenges, Daniel, Helene Jacqmin-Gadda, Cecile Proust, and Jeremie
Guedj. 2006. *A Newton-Like Algorithm for Likelihood Maximization: The
Robust-Variance Scoring Algorithm*. arXiv.
<https://doi.org/10.48550/arxiv.math/0610402>.

Egozcue, J. J., V. Pawlowsky-Glahn, G. Mateu-Figueras, and C.
Barceló-Vidal. 2003. ‘Isometric Logratio Transformations for
Compositional Data Analysis’. *Mathematical Geology* 35 (3): 279–300.
<https://doi.org/10.1023/a:1023818214614>.

Genz, Alan, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, and Torsten Hothorn.
2026. *Mvtnorm: Multivariate Normal and t Distributions*.
<http://mvtnorm.R-forge.R-project.org>.

Gilbert, Paul, and Ravi Varadhan. 2019. *numDeriv: Accurate Numerical
Derivatives*. <http://optimizer.r-forge.r-project.org/>.

Golub, Gene H., and Charles F. Van Loan. 2013. *Matrix Computations*.
JHU Press.

Hager, William W. 1989. ‘Updating the Inverse of a Matrix’. *SIAM
Review* 31 (2): 221–39. <https://doi.org/10.1137/1031049>.

Holland, Paul W., Kathryn Blackmond Laskey, and Samuel Leinhardt. 1983.
‘Stochastic Blockmodels: First Steps’. *Social Networks* 5.
<https://doi.org/10.1016/0378-8733(83)90021-7>.

Oehlert, Gary W. 1992. ‘A Note on the Delta Method’. *The American
Statistician* 46 (1): 27–29.
<https://doi.org/10.1080/00031305.1992.10475842>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Philipps, Viviane, Cecile Proust-Lima, Melanie Prague, Boris Hejblum,
Daniel Commenges, and Amadou Diakite. 2023. *marqLevAlg: A Parallelized
General-Purpose Optimization Based on Marquardt-Levenberg Algorithm*.

Schäfer, Florian, T. J. Sullivan, and Houman Owhadi. 2020. *Compression,
Inversion, and Approximate PCA of Dense Kernel Matrices at Near-Linear
Computational Complexity*. arXiv.
<https://doi.org/10.48550/arxiv.1706.02205>.

Tsagris, Michail, Abdulaziz Alenazi, and Connie Stewart. 2023. ‘Flexible
Non-Parametric Regression Models for Compositional Response Data with
Zeros’. *Statistics and Computing* 33 (5): 106.
<https://doi.org/10.1007/s11222-023-10277-5>.

Tsagris, Michail, and Connie Stewart. 2020. ‘A Folded Model for
Compositional Data Analysis’. *Australian & New Zealand Journal of
Statistics* 62 (2): 249–77. <https://doi.org/10.1111/anzs.12289>.

Vaart, A. W. van der. 2000. *Asymptotic Statistics*. Cambridge
University Press.

van den Boogaart, K. Gerald, Raimon Tolosana-Delgado, and Matevz Bren.
2025. *Compositions: Compositional Data Analysis*.
<http://www.stat.boogaart.de/compositions/>.

Watts, Duncan J., and Steven H. Strogatz. 1998. ‘Collective Dynamics of
“Small-World” Networks’. *Nature* 393 (6684): 440–42.
<https://doi.org/10.1038/30918>.
