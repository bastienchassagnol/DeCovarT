# Repair a numeric vector onto the unit simplex

Clips numerical under/overflow and renormalises so that
\\\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1\\ and
\\\boldsymbol{p}\ge\mathbf{0}\\. This is a **repair / renormalisation**
step for estimated proportions, not a Euclidean projection onto the
simplex and not a statistical-identifiability constraint.

When `open = TRUE`, entries on the closed-simplex boundary are nudged
into the relative interior so that
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
is defined (ALR / Marquardt / Newton starts). That open-simplex path is
the former role of `.starting_simplex()`.

## Usage

``` r
repair_simplex(
  p,
  tolerance = 100 * .Machine$double.eps,
  open = FALSE,
  nms = NULL
)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- tolerance:

  Non-negative tolerance for treating entries as zero (default
  `100 * .Machine$double.eps`).

- open:

  Logical; if `TRUE`, push the repaired vector off the boundary into the
  open simplex.

- nms:

  Optional names attached to the returned vector (e.g. cell-type
  colnames).

## Value

Numeric vector on the simplex \\\Delta^{J-1}\\ (open when
`open = TRUE`).

## See also

[`compositions::clo()`](https://rdrr.io/pkg/compositions/man/clo.html)
for compositional closure.

## Examples

``` r
repair_simplex(c(0.2, 0.3, 0.5 + 1e-12))
#> [1] 0.2 0.3 0.5
repair_simplex(c(1, 0, 0), open = TRUE)
#> [1] 1.000000e+00 2.220446e-14 2.220446e-14
```
