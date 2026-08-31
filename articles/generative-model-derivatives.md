# Derivatives of the DeCovarT generative model under simplex transforms

DeCovarT maximises a Gaussian-convolution log-likelihood for cellular
ratios \boldsymbol{p}\in\Delta^{J-1}. This vignette derives the analytic
first- and second-order derivatives used by the frequentist solvers,
first with respect to unconstrained ambient proportions \boldsymbol{p},
then after the additive log-ratio (ALR) map that enforces the simplex.
It closes with the expected Fisher information and the ALR delta-method
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
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md).

> **Note 1: Notation**
>
> - \boldsymbol{\mu}\in\mathbb{R}^{G\times J}: mean signature.
> - \boldsymbol{\Sigma}\_j\in\mathrm{SPD}(G): cell-type covariance;
>   \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
>   p_j^{2}\boldsymbol{\Sigma}\_j and
>   \boldsymbol{\Theta}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}.
> - \boldsymbol{p}: proportions on the open simplex (ambient score
>   equations treat them as free coordinates in (0,1)^{J}).
> - \boldsymbol{\rho}\in\mathbb{R}^{J-1}: ALR coordinates with reference
>   part J; \boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho}).

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
> [Note 6](#nte-obs-vs-exp)). Point estimates move only slightly,
> because the two terms trade off differently, but likelihood
> **values**, AIC, likelihood-ratio statistics, and any curvature-based
> standard error do change.

> **Warning 4: What “unconstrained” means here**
>
> [Sec. 2](#sec-unconstrained) treats each p\_{j} as an ambient
> coordinate in (0,1) and ignores \sum\_{j}p\_{j}=1. The resulting
> gradient and Hessian are the building blocks of the ALR chain rule in
> [Sec. 3](#sec-constrained). They are **not** the coordinates used by
> Marquardt–Levenberg or Newton–Raphson in the package (those optimise
> in \boldsymbol{\rho}-space). The box-constrained L-BFGS-B path is the
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
from free Euclidean coordinates \boldsymbol{\rho}\in\mathbb{R}^{J-1}
onto the open simplex, so every Newton step stays feasible by
construction.

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

### Additive log-ratio map

Pin the last category as reference (\rho\_{J}\equiv 0). The *additive
logistic* transform \boldsymbol{\psi} (inverse ALR) and its inverse are
([Aitchison 1982](#ref-aitchisonStatisticalAnalysisCompositional1982))

p\_{j} = \frac{e^{\rho\_{j}}}{\sum\_{k\<J}e^{\rho\_{k}}+1}\\ (j\<J),
\qquad p\_{J} = \frac{1}{\sum\_{k\<J}e^{\rho\_{k}}+1}, \qquad \rho\_{j}
= \log\\\left(\frac{p\_{j}}{p\_{J}}\right). \tag{12}

Equivalently,
\boldsymbol{p}=\operatorname{softmax}(\rho\_{1},\ldots,\rho\_{J-1},0).
Package helpers
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
and
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_log_ratio.md)
implement [Eq. 12](#eq-alr-maps); they match
[`compositions::alrInv()`](https://rdrr.io/pkg/compositions/man/alr.html)
/ `alr()` with the last part as reference ([van den Boogaart et al.
2025](#ref-R-compositions)).

**Proposition 1 (ALR is a C^{2} diffeomorphism)**
\boldsymbol{\psi}:\mathbb{R}^{J-1}\to(0,1)^{J} with \sum\_{j}p\_{j}=1 is
bijective and twice continuously differentiable, so every critical point
in \boldsymbol{\rho}-space maps to a unique open simplex composition and
conversely.

### Efficient derivatives of \boldsymbol{\psi}

Write A=\sum\_{j'\<J}e^{\rho\_{j'}}+1 and B\_{i}=A-e^{\rho\_{i}}. The
Jacobian \mathbf{J}\_{\boldsymbol{\psi}}\in\mathcal{M}\_{J\times(J-1)}
is the softmax Jacobian with the reference column dropped,

\frac{\partial p\_{i}}{\partial\rho\_{j}} = \begin{cases}
e^{\rho\_{i}}B\_{i}/A^{2}, & i=j,\\ i\<J,\\
-e^{\rho\_{j}}e^{\rho\_{i}}/A^{2}, & i\neq j,\\ i\<J,\\
-e^{\rho\_{j}}/A^{2}, & i=J, \end{cases} \tag{13}

or compactly \partial
p\_{i}/\partial\rho\_{a}=p\_{i}(\delta\_{ia}-p\_{a}) for a=1,\ldots,J-1.
The Hessian is a third-order tensor of shape (J-1)\times(J-1)\times J;
explicit case distinctions appear in the manuscript appendix and in
[`hessian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_additive_logistic.md).
In practice DeCovarT evaluates the Jacobian / Hessian from the softmax
restriction rather than expanding every Boolean case by hand
([`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md)).

### Change of variables: ambient \boldsymbol{p} to ALR \boldsymbol{\rho}

Let \ell\_{\boldsymbol{\rho}}(\boldsymbol{\rho})
=\ell\_{\boldsymbol{y}\mid\boldsymbol{\zeta}}
\bigl(\boldsymbol{\psi}(\boldsymbol{\rho})\bigr). The chain rule maps
the unconstrained derivatives of [Sec. 2](#sec-unconstrained) into free
Euclidean coordinates.

**Theorem 3 (Constrained score and Hessian)**
\frac{\partial\ell\_{\boldsymbol{\rho}}}{\partial\rho\_{j}} =
\sum\_{i=1}^{J} \frac{\partial\ell}{\partial p\_{i}} \frac{\partial
p\_{i}}{\partial\rho\_{j}}, \tag{14}

\begin{aligned}
\frac{\partial^{2}\ell\_{\boldsymbol{\rho}}}{\partial\rho\_{k}\partial\rho\_{j}}
&= \sum\_{i=1}^{J} \sum\_{l=1}^{J} \frac{\partial
p\_{i}}{\partial\rho\_{j}} \frac{\partial^{2}\ell}{\partial
p\_{i}\partial p\_{l}} \frac{\partial p\_{l}}{\partial\rho\_{k}} \\
&\quad+ \sum\_{i=1}^{J} \frac{\partial\ell}{\partial p\_{i}}
\frac{\partial^{2}p\_{i}}{\partial\rho\_{k}\partial\rho\_{j}}.
\end{aligned} \tag{15}

> **Important 5: Implementation map**
>
> | Object | Package function |
> |:---|:---|
> | Ambient score / Hessian | [`gradient_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_unconstrained.md), [`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md) |
> | Map Jacobian / Hessian | [`jacobian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/jacobian_additive_logistic.md), [`hessian_additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_additive_logistic.md) |
> | ALR score / Hessian | [`gradient_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/gradient_loglik_constrained.md), [`hessian_loglik_constrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_constrained.md) |
>
> [Eq. 15](#eq-chain-second) is the sum of a covariant pull-back of
> \mathbf{H} and a contraction of the ambient score against the map
> Hessian. Omitting the second summand is a common bug when porting
> first-order code to Newton methods.

Extensions on the same footing include **ILR** or **CLR** coordinates
([Pawlowsky-Glahn and Buccianti
2011](#ref-pawlowsky-glahnCompositionalDataAnalysis2011)) and Bayesian
formulations: a Dirichlet or logistic-normal prior on \boldsymbol{p}
pairs naturally with the Gaussian bulk likelihood (the experimental MAP
`CTS` path in `R/03_04_DeCovarT_estimate_CTS_MAP_Bayesian.R` is one
starting point). Regardless of coordinate system, the workflow is the
same: optimise in an unconstrained parameterisation, then map back to
the simplex for interpretation.

## Expected Fisher information and Wald inference

Point estimation alone does not report uncertainty. Under a regular
large-sample regime the MLE is asymptotically normal with covariance
given by the inverse Fisher information. DeCovarT builds Wald standard
errors from the *expected* Fisher information of the multivariate-normal
convolution, then maps that information through the ALR chart.

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
\tag{16}

The first summand is the mean contribution (a \boldsymbol{\Theta}-inner
product of signature columns); the second is the covariance contribution
of the quadratic map
\boldsymbol{p}\mapsto\boldsymbol{\Sigma}(\boldsymbol{p}). See the
multivariate-normal formula on [Wikipedia: Fisher
information](https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution).
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md)
implements [Eq. 16](#eq-fisher-p).

> **Note 6: Observed versus expected information**
>
> The analytic Hessian of [Sec. 2](#sec-unconstrained) is the *observed*
> information (up to sign). Wald intervals in `decovart_fit` use the
> *expected* information [Eq. 16](#eq-fisher-p), which does not depend
> on the residual \boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p} and is
> positive definite under the usual full-rank assumptions on
> (\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_{j}\\).
>
> The two are consistent by construction. Substituting
> \mathbb{E}\[\boldsymbol{r}\]=\boldsymbol{0} and
> \mathbb{E}\[\boldsymbol{r}\boldsymbol{r}^{\top}\]
> =\boldsymbol{\Sigma}(\boldsymbol{p}) into [Eq. 8](#eq-hess-diag) and
> [Eq. 9](#eq-hess-off) cancels the determinant traces against the
> residual traces and returns [Eq. 16](#eq-fisher-p) exactly, i.e.
> \mathbb{E}\[-\mathbf{H}\]=I(\boldsymbol{p}). That identity only holds
> with the factor \tfrac{1}{2} of [Note 3](#nte-halfdet) in place.

### Constrained case: ALR delta method

Fisher information transforms as a covariant quadratic form under the
ALR chart. With \mathbf{J}\_{\boldsymbol{\psi}}
=\partial\boldsymbol{\psi}/\partial\boldsymbol{\rho}^{\top},

I\_{\boldsymbol{\rho}} = \mathbf{J}\_{\boldsymbol{\psi}}^{\top}
I(\boldsymbol{p}) \mathbf{J}\_{\boldsymbol{\psi}}. \tag{17}

Under regularity, \hat{\boldsymbol{\rho}} \overset{a}{\sim}
\mathcal{N}(\boldsymbol{\rho}\_{0},I\_{\boldsymbol{\rho}}^{-1}) ([Vaart
2000](#ref-vaartAsymptoticStatistics2000)). The first-order delta method
([Oehlert 1992](#ref-oehlertNoteDeltaMethod1992)) then returns the
simplex covariance used by
[`vcov.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),

\mathrm{Var}(\hat{\boldsymbol{p}}) \approx
\mathbf{J}\_{\boldsymbol{\psi}} I\_{\boldsymbol{\rho}}^{-1}
\mathbf{J}\_{\boldsymbol{\psi}}^{\top}. \tag{18}

Diagonal square roots are asymptotic standard errors; Wald intervals at
level 1-\alpha are \hat{p}\_{j}\pm z\_{1-\alpha/2}\\\mathrm{SE}\_{j}
with z\_{q}=\Phi^{-1}(q), as in
[`confint.decovart_fit()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md).
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md)
implements [Eq. 17](#eq-fisher-rho)–[Eq. 18](#eq-delta-p). The
construction is undefined on the simplex boundary (the ALR chart blows
up); the helper then returns `NA` with a warning. The ALR chart itself
follows Aitchison’s compositional geometry ([Aitchison
1982](#ref-aitchisonStatisticalAnalysisCompositional1982)).

> **Warning 7: Wald is the cheapest, not the safest, interval**
>
> [Eq. 18](#eq-delta-p) linearises the ALR chart at
> \hat{\boldsymbol{p}}, so it is neither invariant to reparametrisation
> nor confined to \[0,1\], and it degenerates entirely when a proportion
> reaches a simplex face. The companion vignette [MLE properties and
> asymptotic
> inference](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-MLE-properties.html)
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
> bottlenecks that showed up on the hybrid scenario in [Manuscript
> synthetic simulation
> scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-manuscript-scenarios.html#sec-hybrid-deconvolution)
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
#>           s1
#> ct1 0.648736
#> ct2 0.351264
nobs(fit)
#> [1] 1
#> attr(,"n_genes")
#> [1] 3
#> attr(,"n_celltypes")
#> [1] 2
#> attr(,"n_samples")
#> [1] 1
c(
  n_genes = attr(nobs(fit), "n_genes"),
  n_celltypes = attr(nobs(fit), "n_celltypes")
)
#>     n_genes n_celltypes 
#>           3           2
fitted(fit)
#>          s1
#> g1 25.26896
#> g2 36.48736
#> g3 23.94621
residuals(fit)
#>             s1
#> g1 -0.03549669
#> g2 -0.01537423
#> g3  0.09402526
vcov(fit)
#>              ct1          ct2
#> ct1  0.001804864 -0.001804864
#> ct2 -0.001804864  0.001804864
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
#> loglik_multivariate             dmvnorm          difference 
#>           0.9560062          -1.8008094           0.0000000
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
#> gradient_rel  hessian_rel 
#> 3.893819e-11 4.275761e-12
```

Relative discrepancies of order (10^{-8}) are typical on this three-gene
toy; that is the same tolerance used in the test suite.

### Expected information versus Monte Carlo (\[-\])

Matching `numDeriv` shows that the derivatives belong to the implemented
objective. Matching the *expected* Fisher information
[Eq. 16](#eq-fisher-p) to a Monte Carlo average of (-(; )) shows that
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
#> rel_frobenius 
#>   0.004978603
```

The same identity is what closed the factor-of-two inconsistency: with
the old (-) objective, (\[-\]) could not equal (I()). A few dozen draws
already bring the relative Frobenius gap below (10^{-2}) on this toy;
the test suite uses a tighter Monte Carlo.

A related, weaker check compares the empirical covariance of
parametric-bootstrap MLEs to the Cramér–Rao map
[`vcov_alr_delta()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_alr_delta.md).
That comparison is noisy at small replication (see [MLE
properties](https://bastienchassagnol.github.io/DeCovarT/articles/DeCovarT-MLE-properties.html#sec-replication)
for why one bulk column is not an asymptotic sample) and is therefore
left to the inference vignette and to
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md).

## References

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Commenges, Daniel, Helene Jacqmin-Gadda, Cecile Proust, and Jeremie
Guedj. 2006. *A Newton-Like Algorithm for Likelihood Maximization: The
Robust-Variance Scoring Algorithm*. arXiv.
<https://doi.org/10.48550/arxiv.math/0610402>.

Genz, Alan, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, and Torsten Hothorn.
2026. *Mvtnorm: Multivariate Normal and t Distributions*.
<http://mvtnorm.R-forge.R-project.org>.

Gilbert, Paul, and Ravi Varadhan. 2019. *numDeriv: Accurate Numerical
Derivatives*. <http://optimizer.r-forge.r-project.org/>.

Oehlert, Gary W. 1992. ‘A Note on the Delta Method’. *The American
Statistician* 46 (1): 27–29.
<https://doi.org/10.1080/00031305.1992.10475842>.

Pawlowsky-Glahn, Vera, and Antonella Buccianti, eds. 2011.
*Compositional Data Analysis: Theory and Applications*. Wiley.

Philipps, Viviane, Cecile Proust-Lima, Melanie Prague, Boris Hejblum,
Daniel Commenges, and Amadou Diakite. 2023. *marqLevAlg: A Parallelized
General-Purpose Optimization Based on Marquardt-Levenberg Algorithm*.

Vaart, A. W. van der. 2000. *Asymptotic Statistics*. Cambridge
University Press.

van den Boogaart, K. Gerald, Raimon Tolosana-Delgado, and Matevz Bren.
2025. *Compositions: Compositional Data Analysis*.
<http://www.stat.boogaart.de/compositions/>.
