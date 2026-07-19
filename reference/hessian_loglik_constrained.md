# Constrained Hessian of \\\ell\circ\boldsymbol{\psi}\\

Second-order chain rule \$\$ \mathbf{H}\_{\boldsymbol{\rho}} =
\mathbf{J}\_{\boldsymbol{\psi}}^{\mathsf{T}}
\mathbf{H}\_{\boldsymbol{p}} \mathbf{J}\_{\boldsymbol{\psi}}
+\sum\_{i=1}^{J} \frac{\partial\ell}{\partial p_i}\\
\frac{\partial^{2}p_i}{\partial\boldsymbol{\rho}\partial\boldsymbol{\rho}^{\mathsf{T}}}.
\$\$

## Usage

``` r
hessian_loglik_constrained(rho, y, mean_signature_matrix, Sigma)
```

## Arguments

- rho:

  Numeric vector \\\boldsymbol{\rho}\in\mathbb{R}^{J-1}\\.

- y:

  Numeric vector (or one-column matrix)
  \\\boldsymbol{y}\in\mathbb{R}^{G}\\.

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\.

- Sigma:

  Array of cell-type covariances in \\\mathcal{M}\_{G\times G\times
  J}\\.

## Value

Symmetric matrix in \\\mathcal{M}\_{(J-1)\times(J-1)}\\.
