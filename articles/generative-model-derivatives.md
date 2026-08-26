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
\log\det\boldsymbol{\Theta}(\boldsymbol{p}) - \tfrac{1}{2}
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\top}
\boldsymbol{\Theta}(\boldsymbol{p})
(\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}). \tag{6}

> **Warning 3: What “unconstrained” means here**
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
\underbrace{-2p\_{j}\\\mathrm{tr}\bigl(
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
-2\\\mathrm{tr}(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_{i}) +
4p\_{i}^{2}\\\mathrm{tr}\bigl(
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
p\_{i}\partial p\_{j}} = \underbrace{ 4p\_{i}p\_{j}\\ \mathrm{tr}\bigl(
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

> **Important 4: Implementation map**
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

> **Note 5: Observed versus expected information**
>
> The analytic Hessian of [Sec. 2](#sec-unconstrained) is the *observed*
> information (up to sign). Wald intervals in `decovart_fit` use the
> *expected* information [Eq. 16](#eq-fisher-p), which does not depend
> on the residual \boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p} and is
> positive definite under the usual full-rank assumptions on
> (\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_{j}\\).

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

## Numerical speed-ups and solver safeguards

The analytic maps above are only half of a usable optimiser. Practical
bottlenecks that showed up on the hybrid scenario in [Deconvolution use
cases](https://bastienchassagnol.github.io/DeCovarT/articles/DeCoVart-use-cases.html#sec-hybrid-deconvolution)
live in `R/03_03_DeCovarT_estimate_ratios_frequentist.R`.

### Cache a Cholesky factorisation of \boldsymbol{\Sigma}(\boldsymbol{p})

Within one iteration, objective / gradient / Hessian callbacks hit the
same trial \boldsymbol{p}. The helper
[`.sigma_p_factorisation()`](https://bastienchassagnol.github.io/DeCovarT/reference/dot-sigma_p_factorisation.md)
caches a single Cholesky factor and returns
\log\det\boldsymbol{\Sigma}(\boldsymbol{p}) and
\boldsymbol{\Sigma}(\boldsymbol{p})^{-1} without repeating an O(G^{3})
factorisation.

### Guard the box-constrained L-BFGS-B path

[`deconvolute_ratios_L_BFGS_B()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
optimises in \boldsymbol{p} with boxes \[0,1\]^{J} that do not enforce
\mathbf{1}^{\top}\boldsymbol{p}=1. Near \sum\_{j}p\_{j}\approx 0 (or on
a failed [`chol()`](https://rdrr.io/r/base/chol.html)), safeguarded
wrappers return a finite penalty and a zero gradient rather than
aborting [`optim()`](https://rdrr.io/r/stats/optim.html).

### `marqLevAlg` Hessian sign under `minimize = FALSE`

> **Caution 6: `marqLevAlg(minimize = FALSE)` does not flip `hess`**
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

## References

Aitchison, J. 1982. ‘The Statistical Analysis of Compositional Data’.
*Journal of the Royal Statistical Society: Series B (Methodological)* 44
(2): 139–60. <https://doi.org/10.1111/j.2517-6161.1982.tb01195.x>.

Commenges, Daniel, Helene Jacqmin-Gadda, Cecile Proust, and Jeremie
Guedj. 2006. *A Newton-Like Algorithm for Likelihood Maximization: The
Robust-Variance Scoring Algorithm*. arXiv.
<https://doi.org/10.48550/arxiv.math/0610402>.

Gilbert, Paul, and Ravi Varadhan. 2019. *numDeriv: Accurate Numerical
Derivatives*. <http://optimizer.r-forge.r-project.org/>.

Oehlert, Gary W. 1992. ‘A Note on the Delta Method’. *The American
Statistician* 46 (1): 27–29.
<https://doi.org/10.1080/00031305.1992.10475842>.

Philipps, Viviane, Cecile Proust-Lima, Melanie Prague, Boris Hejblum,
Daniel Commenges, and Amadou Diakite. 2023. *marqLevAlg: A Parallelized
General-Purpose Optimization Based on Marquardt-Levenberg Algorithm*.

Vaart, A. W. van der. 2000. *Asymptotic Statistics*. Cambridge
University Press.

van den Boogaart, K. Gerald, Raimon Tolosana-Delgado, and Matevz Bren.
2025. *Compositions: Compositional Data Analysis*.
<http://www.stat.boogaart.de/compositions/>.
