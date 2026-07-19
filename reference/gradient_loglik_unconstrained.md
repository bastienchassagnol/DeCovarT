# Gradient \\\nabla\_{\boldsymbol{p}}\ell\\ of the unconstrained log-likelihood

Analytic gradient of
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md)
with respect to \\\boldsymbol{p}\\. Writing
\\\boldsymbol{\Theta}=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\\ and
\\\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}\\, the
\\j\\-th coordinate is \$\$ \frac{\partial\ell}{\partial p_j} =
-2p_j\\\mathrm{Tr}\\\bigl(\boldsymbol{\Theta}\boldsymbol{\Sigma}\_j\bigr)
+\boldsymbol{r}^{\mathsf{T}}\boldsymbol{\Theta}\boldsymbol{\mu}\_{\cdot
j} +p_j\\\boldsymbol{r}^{\mathsf{T}}
\boldsymbol{\Theta}\boldsymbol{\Sigma}\_j\boldsymbol{\Theta}\boldsymbol{r}.
\$\$

## Usage

``` r
gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Numeric vector \\\nabla\_{\boldsymbol{p}}\ell\in\mathbb{R}^{J}\\.

## Details

Unit tests compare this analytic gradient to a numerical reference from
[`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html) applied
to
[`loglik_multivariate()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_multivariate.md).
For that check the Richardson method is preferred; main `method.args`
knobs:

- `eps`:

  Initial finite-difference step (default `1e-4`).

- `r`:

  Number of Richardson extrapolations (default `4`; tests use `6`).
  Raising `r` usually improves accuracy more safely than shrinking `eps`
  alone.

- `d`, `v`:

  Relative step factor and geometric reduction between extrapolations
  (default `v = 2`).

- `zero.tol`, `show.details`:

  See [`?numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html).

Alternative `method` values: `"simple"` and `"complex"`.

## See also

[`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html),
[`hessian_loglik_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/hessian_loglik_unconstrained.md)
