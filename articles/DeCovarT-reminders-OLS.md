# OLS reminders for linear deconvolution

> **Scope**
>
> Material moved out of the Nature Methods Article appendix to keep the
> main manuscript within editorial length. Notation matches the article
> Methods section and MuSiC-style mean references ([Wang et al.
> 2019](#ref-wangBulkTissueCell2019)).

## Notation

| Symbol | Role |
|:---|:---|
| g=1,\ldots,G | genes |
| j=1,\ldots,J | cell types |
| i=1,\ldots,N | samples (subjects); independence across i |
| \boldsymbol{y}=(y\_{gi})\in\mathbb{R}\_{+}^{G\times N} | bulk; column \boldsymbol{y}\_{\cdot i} |
| \boldsymbol{\mu}=(\mu\_{gj})\in\mathbb{R}^{G\times J} | mean signature; column \boldsymbol{\mu}\_{\cdot j} |
| \boldsymbol{p}=(p\_{ji})\in\\\]0,1\[^{J\times N} | proportions; column \boldsymbol{p}\_{\cdot i} |

Dot notation: \boldsymbol{\mu}\_{\cdot j} is the j-th column of
\boldsymbol{\mu}; \boldsymbol{y}\_{g\cdot} would be the g-th row.
Classical texts often write the design matrix as \boldsymbol{X}; here
the fixed design is the **mean** signature \boldsymbol{\mu}, as in MuSiC
cross-subject means.

Ideal noiseless model for sample i:

\boldsymbol{y}\_{\cdot i}=\boldsymbol{\mu}\\\boldsymbol{p}\_{\cdot i},
\qquad y\_{gi}=\sum\_{j=1}^{J}\mu\_{gj}p\_{ji}. \tag{1}

With G\>J and \operatorname{rank}(\boldsymbol{\mu})=J, the system is
overdetermined ([Abbas et al.
2009](#ref-abbasDeconvolutionBloodMicroarray2009)); uniqueness in the
square full-rank case follows from the Rouché–Capelli theorem
([Shafarevich and Remizov 2013](#ref-shafarevichLinearEquations2013)).

## Ordinary least squares

With additive residual \epsilon\_{gi}, OLS minimises squared error for
each sample (independently):

\hat{\boldsymbol{p}}\_{\cdot i}^{\mathrm{OLS}} \equiv
\arg\min\_{\boldsymbol{p}\_{\cdot i}} \bigl\\\boldsymbol{y}\_{\cdot
i}-\boldsymbol{\mu}\\\boldsymbol{p}\_{\cdot i}\bigr\\\_{2}^{2} =
\arg\min\_{\boldsymbol{p}\_{\cdot i}} \sum\_{g=1}^{G}
\Bigl(y\_{gi}-\sum\_{j=1}^{J}\mu\_{gj}p\_{ji}\Bigr)^{2}. \tag{2}

When \boldsymbol{\mu}^{\top}\boldsymbol{\mu} is invertible, the normal
equations give

\hat{\boldsymbol{p}}\_{\cdot i}^{\mathrm{OLS}} =
\bigl(\boldsymbol{\mu}^{\top}\boldsymbol{\mu}\bigr)^{-1}
\boldsymbol{\mu}^{\top}\boldsymbol{y}\_{\cdot i}. \tag{3}

Existence of this inverse requires full column rank J (no cell-type
profile is an exact linear combination of the others; parent/child
lineages that are collinear cannot both be estimated).

## Homoscedastic Gaussian noise and MLE

Under the classical linear model (Article Fig.~DAG panel **a**),

y\_{gi}=\sum\_{j=1}^{J}\mu\_{gj}p\_{ji}+\epsilon\_{gi}, \qquad
\epsilon\_{gi}\sim\mathcal{N}(0,\sigma\_{i}^{2}), \tag{4}

i.e. y\_{gi}\sim\mathcal{N}\bigl(\sum\_{j}\mu\_{gj}p\_{ji},\\\sigma\_{i}^{2}\bigr)
with fixed (exogenous) \boldsymbol{\mu}. Then the MLE for
\boldsymbol{p}\_{\cdot i} coincides with [Equation 3](#eq-ols-estimate).

## Gauss–Markov assumptions (sketch)

Under weak exogeneity of \boldsymbol{\mu}, homoscedasticity
\operatorname{Var}(\epsilon\_{gi})=\sigma\_{i}^{2} for all g, zero mean
\mathbb{E}(\epsilon\_{gi})=0, and uncorrelated residuals across genes,
the OLS estimator is the BLUE (best linear unbiased estimator). With
i.i.d. Gaussian errors the sample log-likelihood for sample i is (up to
constants)

\ell(\boldsymbol{p}\_{\cdot i},\sigma\_{i}\\\|\\\boldsymbol{y}\_{\cdot
i},\boldsymbol{\mu}) = -G\log\sigma\_{i} -\frac{1}{2\sigma\_{i}^{2}}
\sum\_{g=1}^{G} \Bigl(y\_{gi}-\sum\_{j=1}^{J}\mu\_{gj}p\_{ji}\Bigr)^{2},

so maximising \ell in \boldsymbol{p}\_{\cdot i} recovers
[Equation 2](#eq-ols-task). The residual variance MLE is
\hat\sigma\_{i}^{2}=G^{-1}\sum\_{g}\hat\epsilon\_{gi}^{2}. Independence
across samples i=1,\ldots,N is the modelling assumption used throughout
DeCovarT; microarray literature sometimes questions gene-wise
independence ([Efron 2009](#ref-efronAreSetMicroarrays2009)).

> **What DeCovarT changes**
>
> DeCovarT replaces gene-wise scalar noise with a multivariate
> convolution: latent \boldsymbol{x}\_{\cdot
> j}\sim\mathcal{N}\_{G}(\boldsymbol{\mu}\_{\cdot j},
> \boldsymbol{\Sigma}\_{j}), and \boldsymbol{y}\_{\cdot
> i}\\\|\\\boldsymbol{p}\_{\cdot i}
> \sim\mathcal{N}\_{G}\bigl(\boldsymbol{\mu}\\\boldsymbol{p}\_{\cdot i},
> \sum\_{j}p\_{ji}^{2}\boldsymbol{\Sigma}\_{j}\bigr), with precision
> \boldsymbol{\Theta}\_{j}=\boldsymbol{\Sigma}\_{j}^{-1} (typically from
> `gLasso`). See the article Methods and [softmax / ALR
> derivatives](https://bastienchassagnol.github.io/DeCovarT/articles/softmax-alr-derivatives.md).

Abbas, Alexander R., Kristen Wolslegel, Dhaya Seshasayee, Zora Modrusan,
and Hilary F. Clark. 2009. ‘Deconvolution of Blood Microarray Data
Identifies Cellular Activation Patterns in Systemic Lupus
Erythematosus’. *PloS One* 4.
<https://doi.org/10.1371/journal.pone.0006098>.

Efron, Bradley. 2009. ‘Are a Set of Microarrays Independent of Each
Other?’ *The Annals of Applied Statistics* 3.
<https://doi.org/10.1214/09-aoas236>.

Shafarevich, Igor R., and Alexey O. Remizov. 2013. ‘Linear Equations’.
In *Linear Algebra and Geometry*. Springer Berlin Heidelberg.
<https://doi.org/10.1007/978-3-642-30994-6_1>.

Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao
Li. 2019. ‘Bulk Tissue Cell Type Deconvolution with Multi-Subject
Single-Cell Expression Reference’. *Nature Communications* 10.
<https://doi.org/10.1038/s41467-018-08023-x>.
