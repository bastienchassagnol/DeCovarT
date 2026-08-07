# Score mean profiles with AutoGeneS-style objectives

Computes the two objectives used by AutoGeneS for signature quality. For
columns \\\boldsymbol{\mu}\_{\cdot j},\boldsymbol{\mu}\_{\cdot k}\\ of
\\\boldsymbol{\mu}\\, the **cosine** (angle similarity, akin to Pearson
correlation of centred, unit-norm vectors but here without
mean-centring—only \\\ell_2\\ normalisation) is \$\$
\cos(\boldsymbol{\mu}\_{\cdot j},\boldsymbol{\mu}\_{\cdot k}) = \frac{
\boldsymbol{\mu}\_{\cdot j}^{\mathsf{T}} \boldsymbol{\mu}\_{\cdot k} }{
\\\boldsymbol{\mu}\_{\cdot j}\\\_2\\ \\\boldsymbol{\mu}\_{\cdot k}\\\_2
}, \$\$ and the **Euclidean** separation is \$\$
\\\boldsymbol{\mu}\_{\cdot j}-\boldsymbol{\mu}\_{\cdot k}\\\_2 = \sqrt{
\\\boldsymbol{\mu}\_{\cdot j}\\\_2^{2} +\\\boldsymbol{\mu}\_{\cdot
k}\\\_2^{2} -2\\ \boldsymbol{\mu}\_{\cdot j}^{\mathsf{T}}
\boldsymbol{\mu}\_{\cdot k} }. \$\$ The returned scores are the mean
absolute pairwise cosine (to minimise) and the sum of pairwise Euclidean
distances (to maximise).

## Usage

``` r
compute_mean_profile_objectives(mean_signature_matrix)
```

## Arguments

- mean_signature_matrix:

  Numeric matrix \\\boldsymbol{\mu}\in\mathcal{M}\_{G\times J}\\
  (columns = cell types).

## Value

A named list with `mean_abs_cosine` and `sum_euclidean_distance`.
