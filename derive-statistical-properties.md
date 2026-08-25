# DeCovarT likelihood geometry, constrained inference, and rOpenSci regression standards

## Principal conclusion

The central theoretical claim you proposed — that the DeCovarT sample
log-likelihood is globally concave and therefore has a unique global
maximum because its Hessian is everywhere strictly negative definite —
is **false, even under the favourable assumptions you specified**:
strictly positive expression values, strictly positive-definite
cell-type covariance matrices, (G\>J), linearly independent signature
columns, and true proportions strictly inside the simplex.

There are, however, two strong results that *can* be proved and are, in
my view, the right theoretical statements for the article:

1.  **Global identifiability and a unique population optimum.** If the
    signature matrix is injective on the simplex — full column rank of
    () is a sufficient, though stronger-than-necessary condition — then
    \[ p (p, (p)), (p)=\_jp_j^2_j , \] is identifiable. Consequently the
    expected log-likelihood under the true (p_0) has a **unique global
    maximum at (p_0)**, because its difference from the optimum is a
    Kullback–Leibler divergence.

2.  **Strict local curvature in expectation.** For the correctly
    normalised Gaussian likelihood, the expected Fisher information is
    strictly positive definite under your full-rank assumption. Hence
    the **expected Hessian is strictly negative definite** in
    identifiable directions: \[ E_p\[-H_p(Y)\]=I_p. \] After restriction
    to the simplex or transformation to ALR coordinates, the
    corresponding Fisher information remains strictly positive definite.

These are considerably stronger and cleaner results than merely saying
that an optimiser usually converges, but they must not be confused with
global concavity of the realised, sample-specific likelihood.

More strongly, there is an explicit (J=2,G=3) counterexample satisfying
all your assumptions in which the correctly specified DeCovarT
likelihood has **two distinct, equal, interior global maxima**. The ALR
transformation preserves those two maxima. Thus neither positivity,
covariance positive-definiteness, nor absence of collinearity can
establish uniqueness of the finite-sample MLE.

Before developing that theory, there is also an important
implementation/manuscript correction: the Gaussian determinant term in
both the current article and the current R implementation is off by a
factor of two. Your model is stated as \[ YpN_G!(p,(p)), \] but the
current code evaluates \[ -(p) -12r^((p)){-1}r, \] whereas the
multivariate-normal log-density requires \[ -(p) -12r^((p)){-1}r. \] The
current repository explicitly implements the former expression in
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md).
fileciteturn10file0 The manuscript derives the same doubled
determinant contribution. fileciteturn12file0 This is especially
important because your expected-Fisher formula is actually the one
associated with the **correct** Gaussian likelihood, so at present the
objective and the Wald/Fisher inference are not internally consistent.
The current fitting code gives \[ I\_{jk} =\_j^\_k+ 2p_jp_k,(\_j_k), \]
which is indeed the standard Gaussian mean–covariance Fisher
information. fileciteturn1file0

The Matrix Cookbook identities concerning derivatives of determinants
and inverses are appropriate tools for correcting this calculation.
Petersen and Pedersen’s current Matrix Cookbook is specifically a
reference for matrix identities, inverse derivatives and determinant
derivatives. citeturn14search5

## Geometry of the DeCovarT likelihood

### The exact Gaussian convolution

For one bulk sample, suppress the sample index and write

\[ X_jN_G(\_j,\_j),j=1,,J, \]

independently, and

\[ Y=\_{j=1}^Jp_jX_j. \]

Since every (\_j), choose a Cholesky or symmetric square root
(\_j=L_jL_j^), and write

\[ X_j=\_j+L_jZ_j,Z_jN_G(0,I_G). \]

Then

\[ Y =\_jp_j_j+\_jp_jL_jZ_j =p+B(p)Z, \]

where

\[ B(p)= \begin{bmatrix} p_1L_1 & \cdots & p_JL_J \end{bmatrix} , Z=
\begin{bmatrix} Z_1^\top&\cdots&Z_J^\top \end{bmatrix}

^{!}. \]

Therefore

\[ B(p)B(p)^ =\_jp_j^(2L_jL_j) =\_jp_j^2_j V(p), \]

and

\[ \]

This follows directly from affine closure of Gaussian random variables;
the same Gaussian-affine construction is used explicitly in Chiquet,
Mariadassou and Robin’s probabilistic-PCA formulation, where an affine
Gaussian latent vector induces a Gaussian marginal with covariance
(BB⁺2I). citeturn16view0 Your own repository uses exactly
(V(p)=\_jp_j^2_j). fileciteturn10file0

This observation also answers the Gauss–Hermite part of the question:
**Gauss–Hermite quadrature is not needed for the Gaussian DeCovarT
likelihood**. The latent integral is analytic. Gaussian quadrature,
including Gauss–Hermite rules, is valuable when the latent expectation
is no longer closed form; Golub and Welsch give the classical
construction of Gaussian quadrature rules from orthogonal-polynomial
recurrences. citeturn15search1 In the present Gaussian convolution,
introducing quadrature would replace an exact expression by an
approximation and would not create a concavity property absent from the
exact likelihood.

The Chiquet–Mariadassou–Robin appendix cannot be transferred directly
either. Their variational criterion is proved *biconcave*: concave in
one parameter block with the other held fixed, but explicitly not
jointly concave in general. Their Gaussian-affine convexity lemma works
because the relevant integrand is a convex scalar function evaluated at
an expression affine in each parameter block. citeturn16view0 In
DeCovarT, (p) simultaneously controls the mean and the
**scale/covariance**, and marginalisation produces (V(p)) and
(V(p)^{-1}). That is a different geometry.

### Correct likelihood, score and Hessian

Let

\[ r(p)=y-p, V(p)=\_jp_j^2_j, T(p)=V(p)^{-1}. \]

Up to the (p)-independent constant (-G(2)/2), the true Gaussian
log-likelihood is

\[ \]

By contrast,
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
currently uses (-V). fileciteturn10file0 The consequence propagates
into the analytic gradient: the current manuscript has the
determinant-score contribution (-2p_j(T_j)), whereas the Gaussian
likelihood requires only half of it. fileciteturn12file0

Since

\[ V_j =2p_j_j, \]

the corrected score is

\[ \]

A particularly useful form of the Hessian is the **directional
Hessian**, because its sign can be studied directly. For any direction
(aR^J), define

\[ m_a=a, \]

\[ W_a =D V(p)\[a\] = 2_jp_ja_j_j, \]

and

\[ U_a =D^2V(p)\[a,a\] = 2_ja_j^2_j. \]

Matrix inverse differentiation gives (DT\[a\]=-TW_aT), one of the
standard identities collected in the Matrix Cookbook.
citeturn14search5 Direct differentiation gives

\[ \]

There is **no fixed sign** here. The terms involving (r=y-p) can
dominate with either sign. Positive definiteness of (T), (\_j), or (V)
cannot remove those terms.

This is precisely why standard matrix-convexity facts do not give the
proposed theorem. (X) is concave on the positive-definite cone and the
matrix-fractional quadratic (x^X^{-1}x) has well-understood convexity
properties, but those preservation rules require the arguments to enter
in suitable affine ways. In DeCovarT, \[ r(p)=y-p \] is affine but \[
V(p)=\_jp_j^2_j \] is quadratic. Boyd and Vandenberghe’s
convex-optimisation framework makes exactly this distinction between
known curvature of a matrix function and the validity of the composition
rule used to substitute parameter-dependent arguments.
citeturn18search0

### The two-cell unconstrained case

Now take (J=2), initially without imposing (p_1+p_2=1):

\[ p=(p_1,p_2)^, V=p_1^(2_1+p_2)2_2. \]

Even in an exceptionally simple isotropic setting, \[ \_1=\_2=I_G, \]
one has

\[ V=(p_1^(2+p_2)2)I_G, \]

and therefore

\[ (p) = -G2(p_1^(2+p_2)2) - {2(p_1^(2+p_2)2)} +C. \]

This is plainly not an ordinary concave Gaussian regression likelihood:
its noise variance itself changes with the coefficients.

A concrete example satisfying your stated positivity and rank
assumptions is

\[ G=3, = \begin{pmatrix} 2&1\\ 1&2\\ 1&6/5 \end{pmatrix}

, \_1=\_2=I_3, \]

with

\[ y= \begin{pmatrix}1/2\\1/2\end{pmatrix} = \begin{pmatrix}
3/2\\3/2\\11/10 \end{pmatrix}

. \]

All entries of (y) and () are strictly positive, () has rank two,
(G=3\>J=2), and both covariance matrices are SPD.

At the strictly interior unconstrained point

\[ p=(1/5,3/10), \]

direct differentiation of the **correct** Gaussian likelihood gives

\[ H(p) = \begin{pmatrix} -458908&-751880\\ -751880&-860592
\end{pmatrix}

. \]

Its determinant is

\[ H = -\<0, \]

and its eigenvalues are approximately

\[ -654.53,+53.93. \]

Hence (H) is indefinite. This one calculation is sufficient to disprove
global concavity of the unconstrained likelihood.

### The two-cell simplex and ALR case

Impose the unit simplex:

\[ p_1=t,p_2=1-t,\<t\<1. \]

It is useful to give an even stronger counterexample, because this one
proves not merely non-concavity but **non-uniqueness of the global
maximum**.

Take

\[ G=3, = \begin{pmatrix} 2&1\\ 1&2\\ 1&1 \end{pmatrix}

, \_1=\_2=I_3, \]

and

\[ y= \begin{pmatrix} 2\\2\\5/2 \end{pmatrix}

. \]

Again:

\[ y_g\>0, \_{gj}\>0, \_j, G\>J, ()=2. \]

Put

\[ s(t)=t^(2+(1-t))2. \]

Then

\[ \|y-p(t)\|^2 = 2t^2-2t+ =s(t)+, \]

so

\[ (t) = -s(t) - +C. \]

Its first derivative simplifies exactly to

\[ \]

There are three stationary points:

\[ t_0=, \]

and

\[ \]

Numerically,

\[ t\_-, t\_+. \]

At the central point,

\[ ’’(1/2)=6\>0, \]

so (t=1/2) is a local **minimum**. At the other two stationary points,

\[ ’’(t\_)=-\<0. \]

Both have the same likelihood,

\[ (t\_-)=(t\_+), \]

and both exceed the limiting likelihood at (t) or (t). They are
consequently two distinct **global maxima**, and both are strictly
inside (\]0,1\[).

Thus:

\[ \]

The ALR transformation cannot repair this. For (J=2),

\[ =, t=()=. \]

Because (:R(0,1)) is a bijection, every maximiser in (t)-space
corresponds one-to-one to a maximiser in ()-space. Therefore the two
global maxima remain two global maxima.

More explicitly,

\[ (()) = \[t(1-t)\]^2’‘(t) + t(1-t)(1-2t)’(t). \]

At (), corresponding to (t=1/2), the second term vanishes and

\[ = = \>0. \]

The ALR likelihood therefore has a local **minimum** at (). It is
manifestly not concave.

This scalar identity is the (J=2) version of the chain rule already
correctly documented in your manuscript and package:

\[ H\_ = J\_^H_pJ\_+ *i H*{\_i}. \]

Your softmax/ALR vignette correctly implements both the Jacobian and the
second-derivative tensor. fileciteturn9file0 The second summand is
precisely why a nonlinear reparameterisation does not generally preserve
Hessian definiteness.

### Extension to arbitrary numbers of cell types

The preceding failure is not a peculiarity of (J=2). A two-cell
counterexample embeds into a higher-dimensional model, and small
perturbations of the remaining signatures and covariance matrices can
retain the local non-concavity. Therefore there can be no theorem saying
that (G\>J), SPD covariances and linearly independent columns imply
global concavity for arbitrary (J).

What *does* generalise beautifully is the expected information.

Because (rN(0,V)) under the model,

\[ E(r)=0, E(rr^)=V. \]

Taking expectations in the directional Hessian causes the (U_a)-terms to
cancel:

\[ E(r^TU_aTr) = (TU_a), \]

while

\[ E(r^TW_aTW_aTr) = (TW_aTW_a). \]

Hence

\[ \]

Equivalently,

\[ \]

That is exactly the Fisher formula already encoded in your current
fitting implementation. fileciteturn1file0

The strict-definiteness proof is now simple. For any (a),

\[ a^I(p)a = *{} +* {}. \]

If () has full column rank,

\[ a a. \]

Since (T),

\[ \|T^{1/2}a\|^2\>0. \]

Therefore

\[ \]

and

\[ \]

This is, I think, the appropriate strict-negative-definiteness theorem
for the article: **the expected Hessian, not the observed Hessian, is
strictly negative definite under full-rank signatures**.

On the simplex, choose for example the linear contrast basis

\[ C= \begin{pmatrix} I\_{J-1}\\ -\mathbf1^\top \end{pmatrix}

, \]

so every tangent perturbation is (a=Cu). Then

\[ I\_{} = C^I(p)C. \]

If (I(p)) and (C) has full column rank,

\[ u^C^I C u = (Cu)^I(Cu)\>0, \]

for every (u). Thus

\[ \]

For ALR coordinates the same argument yields

\[ \]

because (J\_) has rank (J-1) at every point of the open simplex. Your
implementation already performs precisely this information
transformation. fileciteturn1file0 Aitchison’s work provides the
classical statistical rationale for log-ratio coordinates on the
simplex. citeturn11search13

Two refinements are worth making to the manuscript.

First, **(G\>J) is sufficient but not necessary**. For unconstrained
full-column-rank (), (GJ) suffices; equality (G=J) is perfectly
compatible with invertibility. On the simplex, only (J-1) independent
directions exist. The weakest mean-identifiability condition is

\[ (C)=J-1, \]

that is, affine independence of the cell-type signatures. Full column
rank of () is stronger than needed. In fact, (J=G+1) signatures can in
principle be affinely independent in (R^G). Your current code
deliberately rejects (J\>G), so it enforces a conservative sufficient
condition rather than the weakest mathematical condition. The current
SRR annotations explicitly describe this rejection.
fileciteturn14file0

Second, the strict positivity of (y) and () **does not help with
concavity**. It is biologically sensible, but the Hessian is governed by
(r), (V), and their derivatives. The counterexample above has strictly
positive (y) and () already. Moreover, if “strictly positive (X_j)” is
intended as an almost-sure support assumption rather than a property of
the observed reference profiles, that assumption is mathematically
incompatible with an untruncated multivariate Gaussian, whose support is
all of (R^G). A truncated Gaussian or log-normal model would be a
different generative likelihood.

## What can be proved about uniqueness and identifiability

Although sample-wise unique MLEs cannot be guaranteed, you can obtain a
particularly clean **population uniqueness theorem**.

Let

\[ P_p=N_G(p,V(p)). \]

For a true (p_0), define the population criterion

\[ Q(p) = E\_{p_0}{(p;Y)}. \]

Then

\[ Q(p)-Q(p_0) = -(P\_{p_0}P_p). \]

Since KL divergence is non-negative,

\[ Q(p)Q(p_0), \]

with equality exactly when (P_p=P\_{p_0}).

Therefore, if the mapping

\[ p(p,V(p)) \]

is injective on the parameter space,

\[ \]

Under your full-column-rank assumption this injectivity follows
immediately from the mean:

\[ p=q (p-q)=0 p=q. \]

So you do **not even need the covariance map** to establish
identifiability under full rank.

This is substantially different from sample concavity:

\[ , \]

and

\[ y. \]

That distinction should be explicit in the paper and in issue \#8. Your
existing issue correctly identifies consistency, identifiability, Fisher
information and possible concavity as separate theoretical questions;
the analysis above shows why they should remain separate.
fileciteturn3file0

There is also a useful existence result on the **closed simplex**. Since
every (\_j),

\[ \_{}V(p) *jp_j^2*{}(\_j). \]

If

\[ = *j*{}(\_j)\>0, \]

then on (\_jp_j=1),

\[ \_jp_j^21J, \]

so

\[ \_{}V(p) \>0. \]

Thus (V(p)) is uniformly positive definite on the entire closed simplex.
The Gaussian likelihood is continuous there, and because the closed
simplex is compact, **at least one global MLE exists**.

That statement no longer follows merely from compactness if you insist
that the estimator live only in the *open* simplex. A true parameter may
be interior while a particular sample has its maximum on a face. ALR
then represents the limiting solution only through (\_j-). This matters
directly for both optimisation and the null-proportion tests discussed
below.

### Affine invariance

There is nevertheless a genuine and useful invariance result that can be
retained.

For any nonsingular (A), set

\[ y^=Ay, ^=A, \_j^(=A_jA). \]

Then

\[ r^=A r, V^(=AVA). \]

Consequently

\[ (r⁾(V⁾{-1}r^ = r^V^{-1}r, \]

while

\[ V^ = V+2\|A\|. \]

The transformed log-likelihood differs only by a constant independent of
(p), so

\[ \]

For an affine transformation (xA(x-b)) on the **simplex**, the common
translation can also be absorbed into each signature column because
(\_jp_j=1):

\[ A!(y-b-\_jp_j(\_j-b)) = A(y-p). \]

This provides a rigorous foundation for applying the same gene-wise
centring/scaling transformation to bulk (y), signatures (), and
covariances. Your current fitting documentation already describes this
affine equivariance. fileciteturn1file0

It does **not** imply invariance under nonlinear transformations such as
(y). Your manuscript already notes that \[ (p)()p \] in general. That
distinction is correct.

## Likelihood-ratio inference and confidence intervals

### The constrained likelihood-ratio statistic

A likelihood-ratio test does not actually require a delta-method
transformation. A smooth one-to-one reparameterisation leaves likelihood
ratios invariant.

For testing whether cell type (j) is absent,

\[ H_0:p_j=0 \]

against

\[ H_1:p_j, \]

define

\[ p = \_{p^{J-1}} (p), \]

and the reduced estimate

\[ p^{(0)} = \_{} (p). \]

Equivalently, remove cell type (j) and fit a ((J-1))-part composition on
its own simplex. Then

\[ \]

The reduced model has (J-2) free compositional coordinates versus (J-1)
in the full model.

For an **ordinary interior, regular equality restriction**, Wilks’
theorem would give

\[ D_j_1^2. \]

Wilks’ original theorem is precisely the classical large-sample
likelihood-ratio result for regular composite hypotheses.
citeturn14search1

But (p_j=0) is **not an interior point of the simplex**. It is on its
boundary. Moreover, ordinary ALR coordinates cannot represent it at a
finite value:

\[ p_j \_j- \]

when (j) is not the ALR reference part.

Therefore, using a standard (\_1^2) Wilks calibration for the null
(p_j=0) is generally incorrect.

For a single active non-negativity constraint, with all other
proportions strictly positive and all nuisance directions identifiable,
the standard tangent-cone limit is the Chernoff/Self–Liang mixture

\[ \]

Chernoff introduced the likelihood-ratio geometry at boundaries, and
Self and Liang develop the general MLE/LRT theory for parameters on the
boundary, showing that mixtures of chi-squared laws arise in many such
cases. citeturn14search8turn11search0

For (D_j\>0), the corresponding simple one-boundary p-value is therefore

\[ \]

At the conventional 5% level the one-boundary critical value is

\[ ^2\_{1,0.90}, \]

not the regular-Wilks value (3.841).

That (50{:}50) result should **not** be hard-coded as a universal rule.
If two or more components are simultaneously zero, the limiting law is
generally a chi-bar-square distribution whose weights depend on the
local information geometry. Self and Liang also show that some
nuisance-boundary configurations do not even yield the simplest
chi-square-mixture form. citeturn11search0

### Why ALR should not be used to encode the null

The full interior model can be optimised in ALR coordinates, but the
null model (p_j=0) should be fitted separately on the lower-dimensional
simplex.

That gives a clean implementation:

\[ q^{J-2} p^{(0)}=(q_1,,0,,q\_{J-1}), \]

with the remaining (J-1) cell types given their own ((J-2))-dimensional
ALR map.

This is preferable to setting

\[ p_j= \]

for a tiny (), because the latter tests (p_j=), not (p_j=0), and the
result becomes tolerance-dependent.

It is also preferable to the present direct L-BFGS-B route if the latter
is intended to represent an exact simplex likelihood: the current
implementation optimises over box-constrained (p_j) and only afterwards
closes the result by dividing by (\_jp_j). The source itself notes that
box constraints do not impose (\_jp_j=1). fileciteturn6file0 That
post-hoc normalisation is not generally equivalent to maximising the
likelihood on the simplex, because your likelihood is not scale
invariant in (p).

For formal reduced/full LRT fitting I would therefore use either:

\[ , \]

or an optimisation method that imposes

\[ p_j,^p=1 \]

simultaneously.

### Delta method and Fisher information

For an interior true composition, let

\[ p=(), J\_=. \]

The expected Fisher information transforms exactly as a covariant
quadratic form:

\[ \]

Under a genuine regular large-sample regime,

\[ N!( *0,, I*^{-1} ), \]

and the first-order delta method yields

\[ \]

This is the formula currently used in
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md),
and the current code deliberately warns that the Wald covariance becomes
undefined at the simplex boundary. fileciteturn1file0

There is an important distinction from the **observed Hessian**. For a
scalar objective composed with (),

\[ H\_^{} = J\_^H_p^{}J\_+ *i g*{p,i}H\_{\_i}. \]

For Fisher information, by contrast, regular score identities remove the
second term in expectation, giving simply (J^I J). Confusing these two
transformation rules is one of the main delta-method pitfalls in this
problem.

At a constrained stationary point there is an additional simplification.
The first-order KKT condition in (p)-space is

\[ g_p= \]

for an interior simplex optimum. Because

\[ \_ip_i()=1 \]

identically,

\[ *iH*{\_i}=0, \]

so

\[ *ig*{p,i}H\_{\_i} = *iH*{\_i}=0. \]

Hence at such an optimum,

\[ H\_ = J\_^H_pJ\_. \]

That establishes **local** curvature equivalence at a constrained
stationary point, but not global ALR concavity.

### The most important asymptotic caveat

There is a deeper problem with interpreting the current Fisher matrix as
an ordinary asymptotic `vcov`.

For each bulk sample (i), DeCovarT observes **one (G)-dimensional random
vector**

\[ Y_iN_G(p_i,V(p_i)) \]

and estimates a separate (p_i).

Increasing the number of bulk samples (N) does not give more independent
replicates for any one (p_i), because the number of unknown compositions
also increases with (N).

Likewise, the (G) genes cannot automatically be treated as (G)
independent observations: the entire point of DeCovarT is that they are
correlated through (V(p)).

Therefore the usual argument

\[ n(p-p_0)N(0,I^{-1}) \]

does **not** follow merely by taking (N) in the current per-sample
deconvolution design.

This changes how I would describe the current result. For one bulk
vector,

\[ I(p)^{-1} \]

is a model-based Fisher-information/Cramér–Rao quantity and a local
quadratic-likelihood approximation; it is **not automatically the actual
sampling covariance of (p)**.

Your manuscript currently correctly calls the Cramér–Rao result a lower
bound, but then changes to an equality when mapping it back to (p). The
latter should instead be presented as an asymptotic/local approximation
unless an explicit asymptotic regime is supplied.

A classical MLE/Wilks/Wald limit becomes rigorous, for example, if one
observes

\[ Y_1,,Y_n N_G(p,V(p)) \]

with the **same** composition (p). Your usual application does not have
those repeated observations. Alternatively, an increasing-(G) theorem
would require a separate triangular-array/dependence argument; it cannot
be borrowed from ordinary iid MLE theory.

This is probably the most important theoretical qualification to add to
issue \#8.

### Better confidence intervals

For DeCovarT, I would provide at least three inferential modes.

**Profile-likelihood intervals** are the most natural analytic
improvement. For candidate (c),

\[ *{p_j}(c) =* {}(p), \]

and form

\[ R_j(c) = 2{(p)-\_{p_j}(c)}. \]

The confidence set is obtained by retaining (c) for which (R_j(c)) is
below the chosen critical value. It automatically respects the simplex
and can be strongly asymmetric near zero. Interior calibration may use
the regular (\_1^2) approximation; boundary endpoints need boundary or
bootstrap calibration. Wilks’ theorem supplies the regular interior
reference distribution, while Chernoff/Self–Liang explain why the
boundary is different. citeturn14search1turn11search0

**Parametric-bootstrap intervals/tests** are especially attractive here
because your entire conditional distribution is available analytically.
Given a fitted candidate (p), simulate

\[ Y^{(b)} N_G(p,V(p)), \]

refit DeCovarT, and use the empirical distribution of (p^{(b)}) or of
the LR statistic. Bootstrap methodology originates with Efron’s
resampling framework. citeturn14search6 For a boundary LRT, a
restricted parametric bootstrap under (H_0) avoids relying on the simple
(50{:}50) asymptotic approximation and automatically incorporates the
actual (G,J,,) geometry.

There is even a useful finite-sample possibility: conditional on (,) and
a **fully specified** (p_0), Monte Carlo simulation from (N(p_0,V(p_0)))
can calibrate a test statistic to arbitrary simulation precision. For a
composite null such as (p_j=c), nuisance proportions remain and must be
profiled or maximised over to guarantee frequentist coverage; a plug-in
restricted-MLE bootstrap is therefore approximate rather than literally
exact.

**Reference-sample bootstrap** addresses a different missing source of
uncertainty. Your frequentist likelihood treats () and every (\_j) as
known plug-in quantities. In reality they are estimated from purified
reference populations. Resampling the purified observations,
re-estimating (,\_j), and then refitting the bulk composition propagates
this uncertainty. This is particularly important for
covariance/precision estimates, where finite-reference-sample
variability and regularisation can be considerable.

Other defensible options include inversion of a one-sided score test,
higher-order likelihood corrections, or a Bayesian
logistic-normal/Dirichlet compositional model when prior regularisation
is scientifically acceptable. A Bayesian credible interval is of course
not an “exact frequentist CI”, but it is often more stable at the
boundary than an ALR Wald interval.

### Delta/Fisher pitfalls specific to DeCovarT

The main pitfalls are therefore:

| Issue | Consequence |
|----|----|
| Current likelihood has (-V) rather than (-V) | Current score/Hessian objective is not the stated Gaussian log-likelihood, while the Fisher formula corresponds to the proper Gaussian model. fileciteturn10file0turn1file0 |
| Only one multivariate bulk observation per (p_i) | (I^{-1}) is not automatically an asymptotic covariance; ordinary iid Wald/Wilks theory needs an explicit replication/asymptotic regime. |
| (p_j=0) is a boundary hypothesis | Ordinary Wald normality and ordinary (\_1^2) Wilks calibration fail; tangent-cone/boundary or bootstrap inference is required. citeturn11search0 |
| ALR only parameterises the open simplex | A zero proportion occurs at infinite ALR coordinate; exact reduced models must be fitted on simplex faces. |
| (J\_) degenerates near a face | ALR standard errors can explode; delta linearisation deteriorates as (p_j). |
| (J)-dimensional covariance of (p) is singular | (^p=1) forces rank at most (J-1); inversion belongs in tangent/ALR coordinates. |
| ALR reference part is arbitrary | ()-covariance depends on the chosen reference; inference mapped back to (p) should be reference-invariant up to numerical error. |
| (,\_j) are plug-in estimates | Conditional Fisher intervals omit reference-estimation uncertainty. |
| Sparse/regularised precision estimation | Penalisation bias and tuning uncertainty are not represented by the nominal Gaussian Fisher matrix. |
| Model misspecification | A Gaussian information matrix need not give robust frequentist coverage under non-Gaussian transcriptomic distributions. |

## Boundary and identifiability edge cases

### A cell type is a linear combination of the others

Suppose there is (a) with

\[ a=0. \]

Then the first term in the information quadratic form disappears in that
direction:

\[ a^I_pa = \|T^({1/2}W_aT){1/2}\|\_F^2. \]

So collinearity of mean signatures does **not automatically imply total
non-identifiability** in DeCovarT. The covariance structure may rescue
it.

If

\[ W_a = 2_jp_ja_j_j, \]

the Fisher information may still be locally positive in that direction.

The genuinely dangerous local condition is

\[ \]

for some admissible tangent direction (a), because then

\[ a^I_pa=0. \]

The Fisher matrix is singular and first-order local identification
fails.

For global identifiability, the relevant criterion is stronger:

\[ \]

If both

\[ p=q \]

and

\[ \_jp_j^2_j = \_jq_j^2_j \]

hold for two distinct compositions, then (P_p=P_q) exactly and no
estimator can distinguish them from the bulk distribution alone.

That distinction is worth encoding in your collinearity diagnostics. The
current package already warns about rank-deficient signature columns and
tests collinear configurations. fileciteturn14file0turn16file0
Rather than describing every rank deficiency as strict
non-identifiability, the warning could say **“mean-based identifiability
is lost; covariance structure may still identify the composition, but
Fisher conditioning should be checked.”**

### A true cell fraction is zero

Let

\[ p_j=0, \]

while at least two other components remain positive.

The covariance remains strictly positive definite:

\[ V(p) = \_{kj}p_k^2_k. \]

So the Gaussian model itself is perfectly well-defined.

What fails is interior regularity:

\[ p^{J-1}. \]

ALR has no finite representation, the ordinary delta method becomes
singular, and likelihood-ratio limits become boundary distributions.
Your current code already detects near-boundary compositions and warns
that Wald `vcov` is undefined. fileciteturn1file0 This is
statistically appropriate.

One subtlety is that

\[ . \|\_{p_j=0} =0. \]

Thus, to first order, an absent cell type enters through the mean, not
its covariance contribution. On the simplex, introducing (p_j) while
reducing a reference component (p_r) creates mean tangent

\[ \_j-\_r. \]

If this contrast is non-zero, the absent type can still be first-order
identifiable despite its variance contribution being second-order in
(p_j).

### Several proportions are zero

If (m\>1) components are absent, the null is at the intersection of
several simplex faces. Multiple inequality constraints become active,
and the simple

\[ \_0^(2+\_1)2 \]

law no longer suffices. The limiting LR distribution is determined by
projection of a Gaussian score onto the local tangent cone, giving a
chi-bar-square law in common regular cases. Self and Liang’s theory is
the appropriate foundation. citeturn11search0

At a simplex vertex, (J-1) inequalities are active. That is the most
non-regular configuration and should be tested separately.

### Nearly collinear signatures

Even if () has exact full rank, a tiny smallest singular value makes

\[ ^T \]

poorly conditioned. Then theoretical identifiability coexists with high
practical uncertainty.

A useful diagnostic for each fitted sample is therefore not merely

\[ (), \]

but also:

\[ (^T), *{}(I*), (I\_). \]

The last two directly describe the inferential geometry in the simplex
tangent space.

### Positivity and the open simplex

One should distinguish:

\[ p_0^{J-1} \]

as an assumption about the **true biological parameter**, from

\[ p^{J-1} \]

as a claim about a **sample estimator**.

The first does not imply the second.

Consequently, restricting all optimisation to an ALR-open simplex can
prevent attainment of a genuine boundary MLE. For routine estimation
this may be a deliberate biological regularisation; for hypothesis
testing (p_j=0), however, the closure of the simplex must explicitly be
available.

## Changes I would make in the DeCovarT code and manuscript

The first changes are theoretical correctness issues rather than
optional enhancements.

**Correct the Gaussian likelihood before adding further inference.** In
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md),

``` r

log_lik <- -0.5 * sigma_p$log_det -
  0.5 * .bilinear_form(residual, sigma_p$inverse)
```

must replace the current doubled determinant term if equation \[
YpN_G(p,V(p)) \] is indeed the intended generative model. The current
implementation and documentation explicitly use `-sigma_p$log_det`.
fileciteturn10file0

Then update the analytic gradient and Hessian. In particular,

\[ -2p_j(T_j) \]

must become

\[ -p_j(T_j). \]

Your existing numerical derivative tests are very useful — they compare
the analytic unconstrained and constrained gradients/Hessians against
`numDeriv`. fileciteturn16file0 But those tests currently only prove
that the derivatives are consistent with the **implemented objective**.
They cannot detect a common error in both the mathematical objective and
its derivatives. Add a separate direct density test against, for
example, the closed-form `mvtnorm`/manual Gaussian log-density for fixed
(p).

**Replace any proposed “global strict concavity” theorem by three
separate propositions:**

\[ \]

from compactness and uniform (V(p));

\[ \]

from injectivity/KL divergence; and

\[ \]

under full-rank/affine-independent signatures.

Then include the (J=2) counterexample above as a remark explicitly
showing that realised likelihoods may be multimodal.

That is scientifically stronger than making an incorrect convexity
claim: it tells readers exactly which properties are global, which are
local/asymptotic, and which depend on the observed bulk vector.

**Add a multimodality diagnostic.** Since global concavity is false,
optimiser convergence alone cannot certify global optimality. Run
several dispersed ALR starts and retain all converged
solutions/log-likelihoods; warn when distinct solutions have likelihood
differences below tolerance. Your SRR annotation already states that
tests vary ALR starts and acknowledges that second-order solvers can
stall from poor starts. fileciteturn14file0 That should become part
of the returned fit diagnostics rather than just a test concern.

A simple output could contain

``` r

fit$optimisation$starts
fit$optimisation$solutions
fit$optimisation$logLik
fit$optimisation$n_distinct_modes
fit$optimisation$globality_certified <- FALSE
```

unless the problem is in a special case for which globality can actually
be proved.

**Do not describe Levenberg–Marquardt or Newton convergence as proof of
the maximum.** Because an observed Hessian can be indefinite and
multiple modes can exist, a solver may converge to a local maximum,
minimum or saddle unless its stopping rule checks curvature/objective
value. Your tests already compare several optimisation methods, which is
useful. fileciteturn16file0

**Implement an exact simplex-face fitting helper**, for example
conceptually:

``` r
fit_decovart_reduced <- function(
  y,
  mean_signature_matrix,
  Sigma,
  absent,
  ...
)
```

which removes the absent columns/slices, fits the remaining composition,
and reconstructs the (J)-vector with zeros. This becomes the foundation
for profile likelihoods and LRTs.

**Add a `lrtest.decovart_fit()` or dedicated `test_proportion()`
interface** with clearly differentiated calibrations:

``` r

test_proportion(
  object,
  cell_type,
  null = 0,
  method = c("boundary_lrt", "parametric_bootstrap", "wald"),
  B = 1999L
)
```

For `null = 0`, default to the boundary LRT or bootstrap, not Wald.

**Add profile-likelihood confidence intervals** to
[`confint()`](https://rdrr.io/r/stats/confint.html):

``` r

confint(
  object,
  method = c("profile", "wald_alr", "parametric_bootstrap"),
  level = 0.95
)
```

I would make `profile` the statistically preferred frequentist method
and retain the current Fisher/delta intervals as the fast approximation.

**Rename or qualify the Fisher covariance.** The current documentation
says [`vcov()`](https://rdrr.io/r/stats/vcov.html) is the “Cramer–Rao /
expected-Fisher bound mapped to the simplex”, and
[`confint()`](https://rdrr.io/r/stats/confint.html) uses the ALR
delta-Wald approximation. fileciteturn1file0 That is acceptable if it
is explicitly documented as a *model-based local/asymptotic
approximation*, but not as an exact covariance from one (G)-variate
observation.

**Return information diagnostics**, for each sample:

\[ *{}(I*), *{}(I*), (I\_), \]

plus a boundary-distance statistic such as

\[ \_jp_j. \]

This provides direct warnings when the Wald approximation is least
trustworthy.

**Add the two-maxima example as a permanent regression test.** This is
particularly valuable because it prevents a future documentation claim
that the optimiser is globally concave. The test should verify that
starts near the two sides recover

\[ t\_= \]

for the corrected Gaussian likelihood and that both attain the same
objective.

Finally, after correcting the likelihood, regenerate all stored fixture
values in `tests/testthat/fixtures/`; the current tests deliberately
compare against stored bivariate benchmark results.
fileciteturn16file0

## rOpenSci regression standards and implementation priorities

Your attached SRR report dated 21 August 2026 shows **85/116 standards
tagged as satisfied**, with 28 tagged not applicable and only three
explicit TODOs: **G1.5, G5.7 and RE5.0**. The current repository’s
`srr-stats-standards.R` confirms those outstanding
performance/reproducibility items and has already documented many
regression-specific NAs, including the absence of a formula interface,
forecast methods and separate predictor/response missing-value
semantics. fileciteturn14file0

The official rOpenSci standards say all statistical packages must
satisfy the general standards plus at least one applicable category, and
explicitly permit standards to be marked non-applicable where justified.
citeturn12view0 Your decision not to create an artificial `formula`
interface is therefore defensible: RE1.0 itself permits omission when
reasons are explicitly documented. citeturn13view0 Your current SRR
file now does exactly that, explaining that DeCovarT is a
convolution/deconvolution model rather than a conventional
predictor-response formula model. fileciteturn14file0

The highest-value remaining work is not to add `lm`-like surface
features. It is to strengthen statistical correctness and benchmarking.

### Performance scaling should become an actual package artefact

RE5.0 asks regression software to document scaling between input size
and computational speed. citeturn13view2 G5.7 similarly asks for
algorithm-performance tests as properties of the data change. rOpenSci’s
general testing standards also ask for parameter-recovery,
edge-condition and noise-susceptibility tests. citeturn12view0

Your current SRR file still labels G5.7 and RE5.0 TODO.
fileciteturn14file0 I would not defer these entirely to a later
benchmark paper, because a modest in-package/extended benchmark would
satisfy the standard much more cleanly.

Benchmark separately against:

\[ G{20,50,100,200,}, \]

\[ J{2,3,5,10,}, \]

and covariance structure/condition number.

For one likelihood evaluation, covariance assembly is approximately

\[ O(JG^2), \]

and dense Cholesky factorisation is

\[ O(G^3). \]

Your implementation already recognises the (O(G^3)) cost and caches the
covariance factorisation between objective/gradient/Hessian calls at the
same trial point. fileciteturn10file0 That is exactly the kind of
implementation feature RE5.0 should quantify.

Report not merely elapsed time, but:

\[ , , , , \]

versus (G,J,(V)), and perhaps signature overlap.

### Theoretical tests should be upgraded in light of non-concavity

RE7.0 and RE7.1 ask for noiseless exact predictor relationships and
exact predictor-response relationships, including expected behaviour in
degenerate cases. citeturn13view3 Your current source now tags these
and contains noiseless recovery tests. fileciteturn14file0

For DeCovarT, I would extend that test family to distinguish:

- identifiable/full-rank configurations;
- nearly collinear but technically identifiable configurations;
- mean-collinear but covariance-identifiable configurations;
- jointly non-identifiable configurations;
- interior unique numerical modes;
- intentionally multimodal configurations;
- one active zero proportion;
- multiple active zero proportions;
- covariance eigenvalues approaching zero without actually violating
  SPD.

This makes the tests genuinely specific to the mathematics of your
package instead of merely translating linear-regression standards
mechanically.

### Convergence handling is almost compliant, but make suppression explicit

RE3.0–RE3.3 require warnings on failure, suppressibility of those
messages, sensible thresholds and user-controllable thresholds.
citeturn13view1 Your current SRR annotations say convergence
diagnostics are stored and that warnings may be suppressed using
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html).
fileciteturn14file0

That technically addresses the standard, but a package-level argument is
cleaner:

``` r

fit_decovart(..., warn_convergence = TRUE)
```

while always retaining

``` r

fit$convergence
```

irrespective of warning suppression. This makes the behaviour
discoverable rather than requiring users to know about the generic
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html) mechanism.

More importantly, given the non-concavity result, `$convergence` should
distinguish:

\[ \]

from

\[ \]

and from

\[ . \]

An optimiser stopping successfully is not a proof of a unique MLE.

### The covariance standard is substantively addressed

G3.1 says software relying on covariance calculations should permit
alternatives to one fixed covariance estimator and document arbitrary
covariance methods. citeturn12view0 Your current SRR source says
DeCovarT accepts user-provided covariance matrices generated by
arbitrary estimation methods and documents a diagonal approximation as
one option. fileciteturn14file0

That is a reasonable interpretation of G3.1. I would reinforce it with
one vignette example showing the exact same
[`fit_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/fit_decovart.md)
call using:

\[ \_j^{}, \_j^{}, \_j^{}, \]

without forcing those estimators to become package dependencies.

This also gives you a natural place to warn that inferential
[`vcov()`](https://rdrr.io/r/stats/vcov.html) conditions on the supplied
(\_j); it does not include uncertainty in their estimation.

### `confint()` deserves higher priority now that RE4.3 is claimed

rOpenSci RE4.3 explicitly asks for confidence intervals on coefficients,
and RE4.6 asks for parameter variance-covariance output.
citeturn13view0turn13view1 Your current code claims both via
expected-Fisher/delta Wald inference. fileciteturn1file0

Given the theoretical findings here, satisfying the *letter* of RE4.3 is
less important than ensuring that its statistical interpretation
satisfies the general rOpenSci requirement that statistical terminology
and assumptions underlying confidence intervals be unambiguously
documented. The general standards explicitly state that the assumptions
and implications of reported confidence intervals should be clarified.
citeturn12view0

I therefore regard the following as a pre-submission issue:

> **Document `confint.decovart_fit(..., method="wald")` as a
> conditional, model-based Fisher/delta approximation valid away from
> simplex boundaries; do not imply exact frequentist coverage from one
> bulk vector. Add profile/bootstrap alternatives before presenting
> `p_j=0` inference as standard Wald inference.**

### Reproducibility and performance claims remain real TODOs

G1.5 asks the software to include code required to reproduce performance
claims made in associated publications. citeturn12view0 Your SRR
source currently says that article-figure reproduction will live in a
companion repository. fileciteturn14file0 That can be reasonable
organisationally, but the package documentation should then provide an
explicit, versioned link and identify which package release/commit
reproduces which manuscript version.

For a statistical peer review, I would make the companion repository
part of the release process rather than leave G1.5 indefinitely marked
TODO.

### A proposed implementation priority order

The highest-priority work is therefore:

| Priority | Change | Reason |
|----|----|----|
| **Critical** | Correct (-V(p)) in likelihood, gradient, Hessian and manuscript | Generative model/objective/Fisher consistency. fileciteturn10file0turn1file0 |
| **Critical** | Remove any claim of global concavity/unique finite-sample MLE | Explicit admissible counterexample disproves it. |
| **Critical** | Recast theorem as identifiability + unique population optimum + (I\_) | This is the rigorous result your assumptions actually support. |
| **High** | Implement profile/reduced-face fitting and boundary LRT | Correct test for (p_j=0); ordinary Wilks/Wald is non-regular. citeturn11search0 |
| **High** | Add parametric-bootstrap inference | Finite-sample calibration is particularly natural for your known Gaussian convolution. citeturn14search6 |
| **High** | Clarify Fisher `vcov` as conditional/local/asymptotic approximation | There is one correlated (G)-vector per sample-specific (p_i), not an iid replication sequence. |
| **High** | Add multimodal/multi-start diagnostics | Non-concavity means optimiser convergence does not certify globality. |
| **High** | Complete G5.7 and RE5.0 scaling benchmarks | These are the main remaining SRR TODOs; RE5.0 explicitly requires input-size/runtime scaling. citeturn13view2 |
| **Medium** | Bootstrap purified references to propagate (,) uncertainty | Current Fisher treats plug-in moments as known. |
| **Medium** | Return (*{}(I*)), ((I\_)), boundary distance | Gives users direct identifiability/Wald-quality diagnostics. |
| **Medium** | Reconsider the hard (J\>G) rejection | (G\>J) is sufficient, not the weakest simplex-identifiability criterion. |
| **Medium** | Add explicit `warn_convergence` and richer mode diagnostics | Better RE3.1 interface and clearer interpretation of convergence. citeturn13view1 |
| **Medium** | Versioned reproducibility companion for G1.5 | Needed to substantiate manuscript performance claims. citeturn12view0 |

The most defensible theoretical statement for the revised article would
therefore be something close to:

> **Proposition.** Assume (\_j) for all (j), and let (p) belong to the
> interior of the unit simplex. If the cell-type signature means are
> affinely independent on the simplex (in particular, if () has full
> column rank), then the DeCovarT Gaussian family (N_G(p,\_jp_j^2_j)) is
> identifiable. Its population log-likelihood has a unique global
> maximum at the true composition. Moreover, the expected Fisher
> information restricted to the simplex tangent space, equivalently the
> Fisher information in ALR coordinates, is strictly positive definite.
> These properties imply local regular identifiability at interior
> points, but **do not imply global concavity or uniqueness of the
> realised finite-sample likelihood**, which can be multimodal even for
> (J=2), (G\>J), positive data and positive-definite covariances.

That theorem is fully compatible with the Gaussian affine decomposition,
your existing ALR differential machinery, the matrix-calculus identities
in the Matrix Cookbook, and the Fisher expression already present in
DeCovarT. It also cleanly explains why the convexity lemma from the
Poisson-PCA appendix is inspirational but cannot furnish the stronger
global-concavity result here. citeturn14search5turn16view0
