# Profile log-likelihood of one single-cell bulk sample in ILR coordinates

\\\ell^{\ast}(\boldsymbol{\rho})
=\ell(\boldsymbol{\psi}(\boldsymbol{\rho}),\hat c,\hat\phi)\\, with the
scales recomputed in closed form at every composition (see
[`single_cell_nuisance_scales()`](https://bastienchassagnol.github.io/DeCovarT/reference/single_cell_nuisance_scales.md)).
By the envelope theorem its score is the score of
[`loglik_single_cell()`](https://bastienchassagnol.github.io/DeCovarT/reference/loglik_single_cell.md)
at those scales, mapped by \\\mathbf{J}\_{\psi}^{\top}\\. Its Hessian is
the Schur complement of the observed Hessian in
\\(\boldsymbol{p},c,\phi)\\, mapped by the second-order chain rule.

The profile objective returns `-Inf` where a closed-form scale is not
positive, so the trust-region driver rejects that trial point.

## Usage

``` r
profile_loglik_single_cell_ilr(z, y, mean_signature_matrix, Sigma, V = NULL)

gradient_profile_loglik_single_cell_ilr(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
)

hessian_profile_loglik_single_cell_ilr(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
)
```

## Arguments

- z:

  ILR coordinate of length \\J-1\\.

- y:

  Bulk vector of length \\G\\.

- mean_signature_matrix:

  Reference mean, \\G\times J\\.

- Sigma:

  Cell-level covariances, \\G\times G\times J\\.

- V:

  Optional ILR basis. See
  [`isometric_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/isometric_logistic.md).

## Value

`profile_loglik_single_cell_ilr()`: scalar.
`gradient_profile_loglik_single_cell_ilr()`: vector of length \\J-1\\.
`hessian_profile_loglik_single_cell_ilr()`: symmetric
\\(J-1)\times(J-1)\\ matrix.
