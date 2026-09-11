# Distances, overlap, and covariance information

``` r

library(DeCovarT)
```

This note collects the **component-separation** scores used in
[synthetic scenario
descriptors](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-scenario-descriptors)
and the **covariance-driven** factorial
([§2.2](https://bastienchassagnol.github.io/DeCovarT/articles/fig03-covariance-driven.md)).
The scores describe the J purified Gaussians
\mathcal{N}(\boldsymbol{\mu}\_{\cdot j},\boldsymbol{\Sigma}\_j) that
enter the convolution, not the deconvolution estimator.

## Overlap, total variation, and MixSim

Pairwise MixSim overlap is the sum of two Bayes misclassification
probabilities ([Maitra and Melnykov
2010](#ref-maitraSimulatingDataStudy2010); [Melnykov et al.
2012](#ref-melnykovMixSimPackageSimulating2012)). For components
j\neq\ell with weights \pi_j,\pi\_\ell and densities f_j,f\_\ell,

\Omega\_{j\ell} = \Pr\_{X\sim f_j}\bigl(\pi\_\ell f\_\ell(X)\>\pi_j
f_j(X)\bigr), \qquad \omega\_{j\ell} = \Omega\_{j\ell}+\Omega\_{\ell j}.
\tag{1}

`BarOmega` is the unweighted mean of \omega\_{j\ell} over j\<\ell. The
FSDA implementation of the same construction ([MixSim
help](https://rosa.unipr.it/FSDA/MixSim.html); [source
pointer](https://github.com/UniprJRC/FSDA/blob/ccc82f8ee453b4cd9e3818512e8ea220b94564b3/toolbox/helpfiles/pointersHTML/MixSim.html);
([Riani et al. 2015](#ref-rianiSimulatingMixturesMultivariate2015)))
targets **average** overlap, not the maximum. DeCovarT follows that
choice.

The **histogram similarity** (density overlap) of two laws P,Q is

\operatorname{OVL}(P,Q) = \int \min\\p,q\\\\\mathrm{d}\mu =
1-\operatorname{TV}(P,Q), \tag{2}

with total variation \operatorname{TV}=\tfrac12\int\lvert p-q\rvert
([Nielsen and Sun 2018](#ref-nielsenGuaranteedDeterministicBounds2018)).
For **two equal-weight** Gaussians the MixSim pairwise overlap \omega
equals \operatorname{OVL} and therefore 1-\operatorname{TV}. Unequal
\boldsymbol{\pi} replace f_j by \pi_j f_j in the MAP rule, so the
identity is then approximate.

Bayes error in the equal-prior two-class problem is
P_e=\tfrac12(1-\operatorname{TV})=\tfrac12\operatorname{OVL}.

## Hellinger, Bhattacharyya, and Fisher–Rao

Hellinger distance H_2(f,g)=H_2(g,f) is a **metric** on densities. For
Gaussians it has a closed form through the Bhattacharyya coefficient
\mathrm{BC}=\int\sqrt{fg} ([Hellinger
distance](https://en.wikipedia.org/wiki/Hellinger_distance);
[Bhattacharyya
properties](https://en.wikipedia.org/wiki/Bhattacharyya_distance#Properties)).
The package stores the unweighted pairwise mean as `hellinger`;
`hellinger_weighted` reweights by p_j p_k (a prevalence weight, not a
directed divergence).

On the unit interval H\in\[0,1\] with H^2=1-\mathrm{BC}, Hellinger and
total variation satisfy ([connection with
TV](https://en.wikipedia.org/wiki/Hellinger_distance#Connection_with_total_variation_distance))

H^2 \le \operatorname{TV} \le H\sqrt{2-H^2}. \tag{3}

Because \operatorname{OVL}=1-\operatorname{TV}, those inequalities
reverse for overlap. Locally, Hellinger is equivalent to the Fisher–Rao
metric on parametric families. For univariate normals the Fisher
information metric is
\mathrm{d}s^2=\mathrm{d}\mu^2/\sigma^2+2(\mathrm{d}\sigma/\sigma)^2
([normal
distribution](https://en.wikipedia.org/wiki/Fisher_information_metric#Normal_distribution)).
For covariances at equal means that geodesic coincides (up to scale)
with the affine-invariant Riemannian metric of [Sec. 3](#sec-airm).

Jeffreys (symmetrised KL) remains a supplementary descriptor
([`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)):
it is not bounded and is not a metric.

## Affine-invariant Riemannian distance

The space of symmetric positive-definite (SPD) matrices is a curved
cone. The **affine-invariant Riemannian metric** (AIRM) distance

d_R(A,B) = \bigl\lVert \log\bigl(A^{-1/2}BA^{-1/2}\bigr) \bigr\rVert_F =
\Bigl(\sum_g \log^2\lambda_g(A^{-1}B)\Bigr)^{1/2} \tag{4}

is congruence-invariant
(d_R(XAX^{\mathsf{T}},XBX^{\mathsf{T}})=d_R(A,B)) and
inversion-invariant (d_R(A^{-1},B^{-1})=d_R(A,B)). The geodesic stays
inside the SPD cone; the distance to a singular matrix is infinite.

We **do not** use the Frobenius chord \lVert A-B\rVert_F. That Euclidean
metric treats SPD matrices as a flat vector space, is not
inversion-invariant, can interpolate off the cone (negative
eigenvalues), and produces the swelling effect (determinants of
midpoints larger than both endpoints).
[`spd_affine_invariant_distance()`](https://bastienchassagnol.github.io/DeCovarT/reference/spd_affine_invariant_distance.md)
and
[`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)
implement only AIRM. The descriptor column `riemannian_sigma` is the
unweighted mean of d_R(\Sigma_j,\Sigma\_\ell) over j\<\ell. The distance
is *scale-invariant*: d_R(sA,sB)=d_R(A,B). In fig03 it therefore tracks
the **shape** of the completed graphs, not the overlap knob s.

## Computing overlap in moderate dimension

[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html)
evaluates the exact pairwise misclassification integrals by Davies’
algorithm for linear combinations of chi-squares ([Maitra and Melnykov
2010](#ref-maitraSimulatingDataStudy2010)). That quadrature scales
poorly with the gene dimension G.
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
therefore uses Davies for G\<4 and stratified Sobol Monte Carlo
otherwise.

### MixSim / FSDA (exact, G\<4)

For G\in\\1,2,3\\ DeCovarT calls
[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html) and
returns `BarOmega` (the unweighted mean of the J(J-1)/2 pairwise
overlaps \omega\_{j\ell}). The MixSim **generator** `MixSim()` can
simulate mixtures with a prescribed *average* overlap ([CRAN
MixSim](https://cran.r-project.org/web/packages/MixSim/refman/MixSim.html#MixSim))
but it draws unconstrained means and covariances. It cannot keep the
structural zeros of a graph-constrained precision. Fig03 therefore does
**not** call `MixSim()` as a generator; it only matches `BarOmega` after
the graph is fixed ([Sec. 5](#sec-fig03-scale)).

### Stratified Sobol Monte Carlo (G\ge 4)

[`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md)
estimates \Omega\_{j\ell} by drawing n points **from each component**
(default n=10\\000), then comparing log-densities under the two-class
MAP rule. Two devices are stacked, and they do different jobs:

- **Sobol** fills the unit cube (0,1)^G more evenly than i.i.d. uniforms
  (low-discrepancy / quasi-Monte Carlo). It is still a sample in the
  *uniform* cube, not a Gaussian sample.
- **Inverse-transform sampling** then maps each coordinate
  u_g\mapsto\Phi^{-1}(u_g) and left-multiplies by the Cholesky factor of
  \Sigma_j, so the image is exactly \mathcal{N}(\boldsymbol{\mu}\_{\cdot
  j},\boldsymbol{\Sigma}\_j). Without that map, Sobol points would not
  have the right Gaussian law.

The remaining implementation details are precomputed Cholesky factors,
log-density evaluation, and stratification by component (not draws from
the mixture). When G\ge 4,
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
emits a `cli` note and uses this estimator. Hellinger of the purified
Gaussians remains a cheap closed-form *descriptor*, not a MixSim
`BarOmega` substitute.

``` default
flowchart TD
  A["Component Gaussians"] --> B{"G < 4?"}
  B -->|yes| C["MixSim::overlap Davies"]
  B -->|no| D["overlap_gaussian_mc Sobol"]
```

``` mermaid
flowchart TD
  A["Component Gaussians"] --> B{"G < 4?"}
  B -->|yes| C["MixSim::overlap Davies"]
  B -->|no| D["overlap_gaussian_mc Sobol"]
```

Figure 1: Routes to MixSim-style average overlap. Exact Davies
quadrature is used only for G \< 4; otherwise stratified Sobol Monte
Carlo.

## Fig03: topology fixed, overlap as the knob

The covariance-driven script `scripts/fig03_variance_driven.R` still
draws one sparse graph per cell type (scale-free or cluster SBM on the
mean-collinear types; type 3 held at scale-free). After SPD completion,
each \Sigma_j is replaced by s\Sigma_j. Then \Omega_j(s)=\Omega_j/s:
**every structural zero is preserved**. Larger s inflates the Gaussians,
raises MixSim `BarOmega`, and raises the covariance-information fraction
f\_{\mathrm{cov}} (mean Fisher scales as 1/s; the whitened Frobenius
covariance block is scale-invariant). See the tangent-Fisher callout in
[scenario
descriptors](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#nte-desc-fisher).

[`scale_covariances_to_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/scale_covariances_to_overlap.md)
binary-searches s so that `BarOmega` at the balanced composition matches
a target in \\\text{low},\text{moderate},\text{high}\\. That is the
graph-constrained analogue of MixSim / FSDA average-overlap simulation:
topology is an experimental factor; overlap is a second, stronger
predictor of deconvolution error when means are collinear.

### See also

- Scenario descriptors: [How to build synthetic
  scenarios](https://bastienchassagnol.github.io/DeCovarT/articles/theory-synthetic-scenarios-mean-covariance.html#sec-scenario-descriptors)
- Covariance-driven factorial:
  [§2.2](https://bastienchassagnol.github.io/DeCovarT/articles/fig03-covariance-driven.md)
- Feature-selection overlap monitors: [Appendix
  S6](https://bastienchassagnol.github.io/DeCovarT/articles/supp-S6-feature-selection.md)

### References

Maitra, Ranjan, and Volodymyr Melnykov. 2010. ‘Simulating Data to Study
Performance of Finite Mixture Modeling and Clustering Algorithms’.
*Journal of Computational and Graphical Statistics* 19 (2): 354–76.
<https://doi.org/10.1198/jcgs.2009.08054>.

Melnykov, Volodymyr, Wei-Chen Chen, and Ranjan Maitra. 2012. ‘MixSim: An
R Package for Simulating Data to Study Performance of Clustering
Algorithms’. *Journal of Statistical Software* 51.
<https://doi.org/10.18637/jss.v051.i12>.

Nielsen, Frank, and Ke Sun. 2018. *Guaranteed Deterministic Bounds on
the Total Variation Distance Between Univariate Mixtures*.
arXiv:1806.11311. <https://doi.org/10.48550/arxiv.1806.11311>.

Riani, Marco, Andrea Cerioli, Domenico Perrotta, and Francesca Torti.
2015. ‘Simulating Mixtures of Multivariate Data with Fixed Cluster
Overlap in FSDA Library’. *Advances in Data Analysis and Classification*
9 (4): 461–81. <https://doi.org/10.1007/s11634-015-0223-9>.
