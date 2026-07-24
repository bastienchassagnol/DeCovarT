# Simplex coordinate maps in DeCovarT: additive logistic and additive log-ratio

DeCovarT estimates the cellular ratios \\\boldsymbol{p}\in\Delta^{J-1}\\
of a bulk sample by maximising a log-likelihood over an *unconstrained*
parameter \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\. Two reciprocal maps
move between the two representations, and both are standard transforms
from the analysis of compositional data.

## The two maps

The forward map
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
sends unconstrained coordinates \\\boldsymbol{\rho}\\ to the open
simplex,

\\ p_j=\frac{e^{\rho_j}}{\sum\_{k\<J}e^{\rho_k}+1}\\ (j\<J),\qquad
p_J=\frac{1}{\sum\_{k\<J}e^{\rho_k}+1}. \tag{1}\\

[Eq. 1](#eq-additive-logistic) is a softmax in which the last category
\\J\\ is pinned as a reference (\\\rho_J\equiv 0\\); in Aitchison’s
framework it is the *additive logistic transform*, i.e. the inverse
additive log-ratio map (\\\mathrm{alr}^{-1}\\).

The inverse map
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
sends a composition back to \\\mathbb{R}^{J-1}\\,

\\ \rho_j=\log\\\left(\frac{p_j}{p_J}\right),\qquad j=1,\ldots,J-1.
\tag{2}\\

[Eq. 2](#eq-additive-log-ratio) is the *additive log-ratio*
(\\\mathrm{alr}\\) transform with reference part \\p_J\\, equivalently
the multinomial-logit link with reference category \\J\\.

## Correspondence with compositional transforms

[Table 1](#tbl-coordinate-maps) summarises the two DeCovarT helpers, the
direction of each map, and their equivalent names across the
statistical, machine-learning, and compositional-data-analysis
literatures.

| DeCovarT function | Direction | Transform | Equivalent formulations |
|:---|:---|:---|:---|
| [`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md) | \\\boldsymbol{\rho}\mapsto\boldsymbol{p}\\ (to \\\Delta^{J-1}\\) | inverse additive log-ratio (\\\mathrm{alr}^{-1}\\) | softmax with reference category \\J\\; [`compositions::alrInv()`](https://rdrr.io/pkg/compositions/man/alr.html); `skbio.stats.composition.alr_inv()` |
| [`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md) | \\\boldsymbol{p}\mapsto\boldsymbol{\rho}\\ (to \\\mathbb{R}^{J-1}\\) | additive log-ratio (\\\mathrm{alr}\\) | multinomial-logit link (reference \\J\\); [`compositions::alr()`](https://rdrr.io/pkg/compositions/man/alr.html); `skbio.stats.composition.alr()` |

Table 1: DeCovarT simplex coordinate maps and their equivalents.

## The `compositions` R package

The [`compositions`](https://cran.r-project.org/package=compositions)
package provides a general toolbox for Aitchison geometry on the
simplex. Beyond the additive log-ratio pair `alr()` / `alrInv()` used
here, it implements the centred log-ratio `clr()` / `clrInv()` and the
isometric log-ratio `ilr()` / `ilrInv()`, together with composition
classes (`acomp()`, `rcomp()`) and the perturbation and powering
operations of the Aitchison simplex. DeCovarT relies only on the
additive log-ratio pair because it yields the sparsest
\\(J-1)\\-dimensional parametrisation with an interpretable reference
category, but any of the alternative bases would define an equally valid
unconstrained coordinate system.

The closed-form Jacobians and Hessians of these maps, together with
efficient `R` and `Python` implementations, are derived in the companion
vignette on softmax and additive log-ratio derivatives.
