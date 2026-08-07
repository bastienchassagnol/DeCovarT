# Repair a numeric vector onto the unit simplex

Clips numerical under/overflow and renormalises so that
\\\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1\\ and
\\\boldsymbol{p}\ge\mathbf{0}\\. This is a **repair / renormalisation**
step for estimated proportions, not a Euclidean projection onto the
simplex and not a statistical-identifiability constraint.

## Usage

``` r
repair_simplex(p, tolerance = 100 * .Machine$double.eps)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

- tolerance:

  Non-negative tolerance for treating entries as zero (default
  `100 * .Machine$double.eps`).

## Value

Numeric vector on the simplex \\\Delta^{J-1}\\.

## See also

[`compositions::clo()`](https://rdrr.io/pkg/compositions/man/clo.html)
for compositional closure.

## Examples

``` r
repair_simplex(c(0.2, 0.3, 0.5 + 1e-12))
#> [1] 0.2 0.3 0.5
```
