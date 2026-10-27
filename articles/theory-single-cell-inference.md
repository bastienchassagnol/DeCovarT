# Single-cell DeCovarT inference: profile information, intervals and the REML calibration

> **Scope**
>
> This vignette derives the intervals and tests for the single-cell
> mixture fitted by
> [`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md):
> one bulk sample, one single-cell reference, a composition on the
> simplex and two technical scales estimated with it. The model, the
> closed-form scales and the profile Newton fit are in the [single-cell
> generative and optimisation
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-single-cell-generative-optimisation.html#sec-profile-newton).
> The population-level counterpart, with the regularity conditions, the
> three first-order tests and the boundary theory that this vignette
> reuses, is the [population inference
> vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.md).
>
> 1.  [Sec. 1](#sec-sc-setting) states the model and its scale
>     equivariance.
> 2.  [Sec. 2](#sec-sc-information) gives the information matrix and the
>     Schur complement that carries the cost of estimating the scales.
> 3.  [Sec. 3](#sec-sc-small-dispersion) states the limit that justifies
>     a first-order interval for one bulk sample.
> 4.  [Sec. 4](#sec-sc-intervals) derives the six interval methods of
>     [`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md),
>     and [Sec. 4.1](#sec-sc-reml) the calibration for few genes.
> 5.  [Sec. 5](#sec-sc-nuisance) gives the intervals for the scales,
>     [Sec. 6](#sec-sc-tests) the likelihood-ratio tests and the
>     boundary, and [Sec. 7](#sec-sc-limits) what none of these
>     intervals covers.

## Model and scale equivariance

For bulk sample i with n_i cells and capture ratio a_i, the averaged
mixture of independent cells is

\boldsymbol{y}\mid\boldsymbol{p}
\sim\mathcal{N}\_G\bigl(c\\\tilde{\boldsymbol{\mu}}\boldsymbol{p},\\\phi\\\boldsymbol{S}\_W(\boldsymbol{p})\bigr),
\qquad
\boldsymbol{S}\_W(\boldsymbol{p})=\sum\_{j=1}^{J}p_j\boldsymbol{\Sigma}\_{W,j},
\qquad c=a_in_i,\quad\phi=a_i^{2}n_i, \tag{1}

with \tilde{\boldsymbol{\mu}}\_{\cdot j} the mean of a type-j cell on
raw counts and \boldsymbol{\Sigma}\_{W,j} its within-type covariance.
The composition is charted by
\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})=\operatorname{softmax}(\mathbf{V}\boldsymbol{\rho})
with a Helmert basis \mathbf{V} ([ILR
chart](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-generative-optimisation.html#sec-ilr)).
Concentrating the scales out gives the profile log-likelihood

\ell^{\ast}(\boldsymbol{\rho}) =-\frac{G}{2}\log Q(\boldsymbol{p})
-\frac12\log\det\boldsymbol{S}\_W(\boldsymbol{p})+\text{const}, \qquad
Q(\boldsymbol{p})=\hat{\boldsymbol{r}}^{\top}\boldsymbol{S}\_W(\boldsymbol{p})^{-1}\hat{\boldsymbol{r}},
\tag{2}

with \hat{\boldsymbol{r}}=\boldsymbol{y}-\hat
c(\boldsymbol{p})\tilde{\boldsymbol{\mu}}\boldsymbol{p} the generalised
least-squares residual ([closed-form
scales](https://bastienchassagnol.github.io/DeCovarT/articles/theory-single-cell-generative-optimisation.html#sec-profile-scales)).

**Scale equivariance.** Multiply the bulk by k\>0. The closed forms give
\hat c(k\boldsymbol{y})=k\hat c(\boldsymbol{y}), \hat{\boldsymbol{r}}
scales by k, Q by k^{2} and \hat\phi by k^{2}, so

\ell^{\ast}(\boldsymbol{\rho};k\boldsymbol{y})=\ell^{\ast}(\boldsymbol{\rho};\boldsymbol{y})-G\log
k . \tag{3}

The shift does not depend on \boldsymbol{\rho}. The estimate
\hat{\boldsymbol{p}}, the profile Hessian, every interval built from
them and every likelihood ratio are unchanged, while (\hat c,\hat\phi)
move to (k\hat c,k^{2}\hat\phi). The capture ratio and the number of
cells enter only through c and \phi. A fixed capture ratio is therefore
a relabelling of the scales and cannot change the composition, which is
why a simulation study need not cross it with the other factors.

## Information and the Schur complement

Write \boldsymbol{\eta}=(c,\phi),
\boldsymbol{W}=\boldsymbol{S}\_W(\boldsymbol{p})^{-1} and
\boldsymbol{m}=\tilde{\boldsymbol{\mu}}\boldsymbol{p}. The expected
information of [Eq. 1](#eq-sc-model) has the blocks

I\_{c\phi}=0,\qquad I\_{\phi\phi}=\frac{G}{2\phi^{2}},\qquad
I\_{cc}=\frac{\boldsymbol{m}^{\top}\boldsymbol{W}\boldsymbol{m}}{\phi},\qquad
I\_{p_jc}=\frac{c\\\tilde{\boldsymbol{\mu}}\_{\cdot
j}^{\top}\boldsymbol{W}\boldsymbol{m}}{\phi},\qquad
I\_{p_j\phi}=\frac{\operatorname{tr}(\boldsymbol{W}\boldsymbol{\Sigma}\_{W,j})}{2\phi}.
\tag{4}

The two scales are orthogonal to each other, but neither is orthogonal
to the composition. The asymptotic covariance of \hat{\boldsymbol{\rho}}
is therefore the inverse of the Schur complement

I\_{\boldsymbol{\rho}\mid\boldsymbol{\eta}}
=I\_{\boldsymbol{\rho}\boldsymbol{\rho}}
-I\_{\boldsymbol{\rho}\boldsymbol{\eta}}I\_{\boldsymbol{\eta}\boldsymbol{\eta}}^{-1}I\_{\boldsymbol{\eta}\boldsymbol{\rho}},
\tag{5}

which is never larger than I\_{\boldsymbol{\rho}\boldsymbol{\rho}}:
estimating the scales widens the composition intervals, and the Schur
complement accounts for it without a separate correction. When an
interest parameter is orthogonal to the nuisance, its standard error is
the same whether the nuisance is known or estimated ([Cox and Reid
1987](#ref-coxParameterOrthogonalityApproximate1987));
[Eq. 4](#eq-sc-information) shows that this does not happen here. A
scalar interest parameter can always be made orthogonal by
reparametrisation; a vector one in general cannot.

The fit uses the observed counterpart. By the envelope theorem the
negative Hessian of the profile log-likelihood
[Eq. 2](#eq-sc-profile-recap) is the Schur complement of the observed
information at the closed-form scales ([profile
Hessian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-single-cell-generative-optimisation.html#sec-profile-derivatives)),
so profile-based and joint Wald intervals coincide.
[`vcov_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/vcov_single_cell.md)
and the fit return

\widehat{\operatorname{Var}}(\hat{\boldsymbol{\rho}})=\bigl(-\nabla^{2}\ell^{\ast}(\hat{\boldsymbol{\rho}})\bigr)^{-1},
\qquad \widehat{\operatorname{Var}}(\hat{\boldsymbol{p}})
=\mathbf{J}\_{\boldsymbol{\psi}}\\\widehat{\operatorname{Var}}(\hat{\boldsymbol{\rho}})\\\mathbf{J}\_{\boldsymbol{\psi}}^{\top},
\tag{6}

as `vcov_rho` and `vcov`, with \mathbf{J}\_{\boldsymbol{\psi}} the [ILR
Jacobian](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-generative-optimisation.html#sec-ilr-deriv).
The complementary Schur complement, which removes \boldsymbol{\rho} from
the scale block, gives `vcov_nuisance` for (c,\phi).

## Small-dispersion asymptotics in the number of cells

One bulk sample has no replicate, so the usual limit in the number of
independent observations does not apply, and the G genes are not
independent replicates either. Divide the bulk by c: the per-cell bulk
\tilde{\boldsymbol{y}}=\boldsymbol{y}/c satisfies

\tilde{\boldsymbol{y}}\sim\mathcal{N}\_G\bigl(\tilde{\boldsymbol{\mu}}\boldsymbol{p},\\n^{-1}\boldsymbol{S}\_W(\boldsymbol{p})\bigr),
\qquad n=c^{2}/\phi, \tag{7}

a family whose dispersion 1/n shrinks as the number of cells grows. That
is the small-dispersion limit of Jørgensen
([1987](#ref-j%C3%B8rgensenSmallDispersionAsymptotics1987)). Under it
the maximum-likelihood estimate is asymptotically normal around the
realised composition, with covariance the inverse of
[Eq. 5](#eq-sc-schur), and the likelihood-ratio statistic is
asymptotically \chi^{2}. The approximation improves with the number of
cells of every type, m_j=np_j, not with G. The number of genes enters
through the residual degrees of freedom of the covariance scale
([Sec. 4.1](#sec-sc-reml)). With Gaussian cells the sum is exactly
Gaussian at any n; with counts, [Eq. 7](#eq-sc-per-cell) is a
central-limit approximation that needs every m_j to be large.

## Intervals for the composition

[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
takes the bulk and the reference moments, refits, and returns one of six
intervals. Its `estimate` is the composition that the likelihood
optimises, before any conversion between cell and transcript fractions.

| `method` | Interval | Inside \[0,1\] | On a face |
|:---|:---|:---|:---|
| `"wald"` | \hat p_j\pm z\_{1-\alpha/2}\\\mathrm{SE}\_j from [Eq. 6](#eq-sc-vcov) | no | undefined |
| `"ilr"` (default) | Wald in \boldsymbol{\rho}, mapped by \boldsymbol{\psi} | yes | undefined |
| `"profile"` | likelihood-ratio inversion, scales concentrated out | yes | needs the boundary test |
| `"firth"` | Wald at the penalised estimate, also mapped by \boldsymbol{\psi} | ILR version only | estimate stays interior |
| `"bootstrap"` | parametric percentile | yes | automatic |
| `"reference"` | percentile over resampled references | yes | automatic |

Table 1: Composition intervals of the single-cell mixture.

**Delta method.** The `"wald"` interval uses the diagonal of
[Eq. 6](#eq-sc-vcov). It is symmetric and can leave \[0,1\] for a rare
type.

**Inverse ILR chart.** The `"ilr"` interval forms the Wald interval in
\boldsymbol{\rho} and maps it back. For J=2 the endpoints of the single
coordinate are mapped by \boldsymbol{\psi}. For J\>2, draws
\boldsymbol{\rho}^{(b)}\sim\mathcal{N}\bigl(\hat{\boldsymbol{\rho}},\widehat{\operatorname{Var}}(\hat{\boldsymbol{\rho}})\bigr),
b=1,\ldots,`n_draws`, are mapped to
\boldsymbol{\psi}(\boldsymbol{\rho}^{(b)}) and the interval is the
marginal quantile of each proportion. Every draw is a composition, so
the interval stays inside the simplex. Pass `seed` for a reproducible
Monte Carlo map.

**Profile likelihood.** The `"profile"` interval
([`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md))
fixes p_j at a trial value, refits the other proportions and both scales
with
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md),
and keeps the values whose deviance stays below the \chi^{2}\_{1}
quantile ([Wilks
1938](#ref-wilksLargeSampleDistributionLikelihood1938)). It does not
depend on the chart and stays in \[0,1\]. An endpoint that lands on a
face is not a Wilks interval and must be read with the boundary test of
[Sec. 6](#sec-sc-tests).

**Firth penalty.** The `"firth"` interval maximises
\ell^{\ast}(\boldsymbol{\rho})+\tfrac12\log\det\\-\nabla^{2}\ell^{\ast}(\boldsymbol{\rho})\\
([`fit_firth_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_firth_single_cell.md))
([Firth 1993](#ref-firthBiasReductionMaximum1993)), then reports the
Wald interval at the penalised estimate. The penalty uses the observed
profile information, with the scales concentrated out, so it adds no
free scale. The information in \boldsymbol{\rho} inherits the factor
\mathbf{S}(\boldsymbol{p}) of the ILR Jacobian, whose determinant
vanishes at a face, so the penalty tends to -\infty there and the
penalised estimate is always interior. Firth’s exact removal of the
first-order bias holds for the canonical parameter of an exponential
family; this curved Gaussian family needs a simulation check ([Firth
penalisation](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-firth)).

**Parametric bootstrap.** The `"bootstrap"` interval
([`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md))
draws \boldsymbol{Y}^{(b)}\sim\mathcal{N}\bigl(\hat
c\\\tilde{\boldsymbol{\mu}}\hat{\boldsymbol{p}},\hat\phi\\\boldsymbol{S}\_W(\hat{\boldsymbol{p}})\bigr),
refits the profile model, and returns percentile intervals for the
composition and for both scales ([Efron
1979](#ref-efronBootstrapMethodsAnother1979)). The reference moments
stay fixed.

**Reference bootstrap.** The `"reference"` interval
([`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md))
propagates the uncertainty of the reference itself. When every type has
at least two donors it resamples donors within each type, rebuilds
(\tilde{\boldsymbol{\mu}},\boldsymbol{\Sigma}\_{W,j}) and refits the
observed bulk. With one donor it resamples cells within each type and
regenerates the bulk from the fitted model with the original moments,
then refits with the resampled moments. Drawing that bulk from the
resampled moments instead would hide the reference error, since the bulk
and the reference would then agree by construction. Shuffling cell-type
labels is not an interval: the fitted composition moves with the labels,
which is equivariance, not sampling variability.

### Few genes: the REML calibration

The covariance scale is estimated from the whitened residual,
\hat\phi=Q(\hat{\boldsymbol{p}})/G, after fitting J mean parameters (J-1
coordinates of \boldsymbol{\rho} and c). It has \nu=G-J residual degrees
of freedom, and the maximum-likelihood divisor G biases it downwards by
a factor of about \nu/G. Every interval that plugs in \hat\phi is then
too narrow when G is small: with G=20 and J=4, \nu/G=0.8.

`variance = "reml"` in
[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
and
[`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md)
corrects for those degrees of freedom, as restricted maximum likelihood
does for a variance component:

- the covariances of [Eq. 6](#eq-sc-vcov) are inflated by G/\nu;
- Wald and ILR quantiles come from t\_{\nu} instead of the normal, with
  a multivariate t for the Monte Carlo ILR map when J\>2;
- \phi gets the chi-square interval
  \bigl\[G\hat\phi/\chi^{2}\_{\nu;1-\alpha/2},\\G\hat\phi/\chi^{2}\_{\nu;\alpha/2}\bigr\];
- the profile interval uses the cut-off

D\_{\alpha}=G\log\Bigl(1+\frac{F\_{1,\nu;1-\alpha}}{\nu}\Bigr), \tag{8}

in place of \chi^{2}\_{1;1-\alpha}.

[Eq. 8](#eq-sc-reml-cutoff) is exact in a linear model with an unknown
variance: with the variance profiled out, twice the log likelihood ratio
for one linear constraint equals G\log(1+F/\nu) for an F\_{1,\nu}
statistic. As \nu\to\infty it tends to \chi^{2}\_{1}. The bootstrap
methods ignore `variance`, because they draw \hat\phi afresh in every
replicate. The Monte Carlo check of these calibrations is
`scripts/benchmark_metrics/check_confidence_interval_coverage.R`: the
maximum-likelihood intervals under-cover with G=10 genes, and the REML
calibration brings them back to the nominal level.

## Intervals for the two scales

[`confint_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_single_cell.md)
also returns `nuisance`, intervals for c and \phi. The default
`nuisance$interval` is a Wald interval on (\log c,\log\phi) from the
complementary Schur complement, which removes \boldsymbol{\rho}, mapped
back by the exponential, so both limits are positive.
`nuisance$interval_wald` is the symmetric normal interval on (c,\phi)
itself. Under `variance = "reml"`, \phi gets the chi-square interval of
[Sec. 4.1](#sec-sc-reml). The bootstrap methods replace both with
percentile intervals.

I\_{c\phi}=0 does not reduce this block to two one-parameter problems:
the complement still accounts for the estimation of the composition. The
derived quantities \hat a=\hat\phi/\hat c and \hat n=\hat c^{2}/\hat\phi
inherit the downward bias of \hat\phi, so \hat n is biased upwards. They
describe the sample; \hat n is a number of cells only when the capture
ratio is known, from the same platform or from paired samples
([Erdmann-Pham et al.
2021](#ref-erdmann-phamLikelihoodbasedDeconvolutionBulk2021)).

## Likelihood-ratio tests and the boundary

[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)
tests a value of the composition with the scales profiled out under both
hypotheses. An interior null, p_j=c_0\in(0,1), is calibrated by \chi^{2}
with as many degrees of freedom as constraints. A null on a face, p_j=0,
puts the true parameter on the boundary of the simplex: the limit is the
projection of a Gaussian score on the tangent cone, and for one absent
type with orthogonal information the statistic follows
\tfrac12\chi^{2}\_{0}+\tfrac12\chi^{2}\_{1} ([Chernoff
1954](#ref-hermanDistributionLikelihoodRatio1954); [Self and Liang
1987](#ref-selfAsymptoticPropertiesMaximum1987)).
[`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md)
counts the null proportions on the boundary and switches to the
chi-bar-square calibration automatically. The restricted fit on a face,
[`restricted_mle_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_single_cell.md),
uses the same trust-region Newton driver on the remaining coordinates,
with the scales concentrated out. The theory and the chi-bar-square
weights are in the [population inference
vignette](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-boundary).

The L-BFGS-B branch of
[`deconvolute_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_single_cell.md)
can return an exact zero. There the ILR coordinates, the Wald
covariances and the Wald and ILR intervals are undefined (`NA`); the
profile interval and the boundary test are the procedures that remain
valid.

## What these intervals do not cover

- **The reference is a plug-in.** Every method except the reference
  bootstrap treats (\tilde{\boldsymbol{\mu}},\boldsymbol{\Sigma}\_{W,j})
  as known. No sandwich or higher-order term corrects for the estimation
  of \boldsymbol{\Sigma}\_{W,j}; the reference bootstrap is the
  implemented route.
- **One donor has no between-donor layer.** \phi inflates the covariance
  globally but cannot correct its shape, so with a reference from one
  donor the intervals describe within-donor variation only.
- **The variance function must be right.** The determinant term of the
  likelihood adds information only if
  \phi\\\boldsymbol{S}\_W(\boldsymbol{p}) is the true covariance.
  Otherwise the quasi-score of the [generalised linear model
  view](https://bastienchassagnol.github.io/DeCovarT/articles/theory-population-inference.html#sec-glm)
  is the safer estimator.
- **Gaussian cells.** The sum of Gaussian cells is exactly Gaussian. For
  count data, [Eq. 7](#eq-sc-per-cell) is a central-limit approximation
  that needs many cells of every type.

| Question | Procedure | Function |
|:---|:---|:---|
| interval for a proportion, interior estimate | ILR inverse chart, or profile | `confint_single_cell(method = "ilr")`, [`confint_profile_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_single_cell.md) |
| few genes, G-J small | the same, with the REML calibration | `variance = "reml"` |
| proportion close to zero | Firth penalty, or profile | `confint_single_cell(method = "firth")` |
| is type j absent? | boundary likelihood-ratio test | [`lrt_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_single_cell.md) |
| no trust in the asymptotics | parametric bootstrap | [`bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_single_cell.md) |
| uncertainty of the reference | reference bootstrap | [`reference_bootstrap_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_single_cell.md) |
| capture ratio and number of cells | log-scale Wald or REML chi-square | `confint_single_cell()$nuisance` |

Table 2: Which single-cell procedure answers which question.

## References

Chernoff, Herman. 1954. ‘On the Distribution of the Likelihood Ratio’.
*The Annals of Mathematical Statistics* 25 (3): 573–78.
<https://doi.org/10.1214/aoms/1177728725>.

Cox, D. R., and N. Reid. 1987. ‘Parameter Orthogonality and Approximate
Conditional Inference’. *Journal of the Royal Statistical Society:
Series B (Methodological)* 49 (1): 1–18.
<https://doi.org/10.1111/j.2517-6161.1987.tb01422.x>.

Efron, B. 1979. ‘Bootstrap Methods: Another Look at the Jackknife’. *The
Annals of Statistics* 7 (1): 1–26.
<https://doi.org/10.1214/aos/1176344552>.

Erdmann-Pham, Dan D., Jonathan Fischer, Justin Hong, and Yun S. Song.
2021. ‘Likelihood-Based Deconvolution of Bulk Gene Expression Data Using
Single-Cell References’. *Genome Research* 31 (10): 1794–806.
<https://doi.org/10.1101/gr.272344.120>.

Firth, David. 1993. ‘Bias Reduction of Maximum Likelihood Estimates’.
*Biometrika* 80 (1): 27–38. <https://doi.org/10.1093/biomet/80.1.27>.

Jørgensen, Bent. 1987. ‘Small Dispersion Asymptotics’. *Brazilian
Journal of Probability and Statistics* 1 (1): 59–90.

Self, Steven G., and Kung-Yee Liang. 1987. ‘Asymptotic Properties of
Maximum Likelihood Estimators and Likelihood Ratio Tests Under
Nonstandard Conditions’. *Journal of the American Statistical
Association* 82 (398): 605–10.
<https://doi.org/10.1080/01621459.1987.10478472>.

Wilks, S. S. 1938. ‘The Large-Sample Distribution of the Likelihood
Ratio for Testing Composite Hypotheses’. *The Annals of Mathematical
Statistics* 9 (1): 60–62. <https://doi.org/10.1214/aoms/1177732360>.
