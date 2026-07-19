# Project estimated proportions onto the unit simplex

Clips numerical under/overflow and renormalises so that
\\\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1\\ and
\\\boldsymbol{p}\ge\mathbf{0}\\.

## Usage

``` r
enforce_identifiability(p)
```

## Arguments

- p:

  Numeric vector \\\boldsymbol{p}\in\mathbb{R}^{J}\\.

## Value

Numeric vector on the simplex \\\Delta^{J-1}\\.
