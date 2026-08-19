# Validate generative-model parameters \\\theta\\

Checks that `true_theta` stores the DeCovarT mixture parameters:

- `p`: mixing proportions of length \\J\\, or a matrix \\J\times N\\ of
  sample-wise proportions (each column on the simplex);

- `mu`: mean signature \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (also accepts \\J\times G\\);

- `sigma` and/or `Theta`: array \\G\times G\times J\\ of covariances
  and/or precision matrices.

## Usage

``` r
check_true_theta(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
)
```

## Arguments

- true_theta:

  Named list with the fields above.

- require_p:

  Logical; if `TRUE` (default), `p` must be present and valid. If
  `FALSE`, a missing `p` is allowed (callers may default it).

- J:

  Optional expected number of cell types. When `NULL`, inferred from
  `mu` / the third dimension of `sigma` or `Theta`.

- second_moment:

  Which second-moment field to require: `"sigma"`, `"Theta"`, or
  `"either"` (default: at least one of `sigma` / `Theta`). Matching is
  case-insensitive.

## Value

`TRUE` invisibly if the structure is valid; otherwise stops with an
informative error.

## See also

[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
[`compute_average_jeffreys()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_jeffreys.md)

## Examples

``` r
theta <- list(
  p = c(0.5, 0.5),
  mu = cbind(c(0, 0), c(3, 0)),
  sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
)
check_true_theta(theta)
```
