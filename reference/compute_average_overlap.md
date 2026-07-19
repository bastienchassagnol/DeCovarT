# Average pairwise overlap of a Gaussian mixture

Approximates mixture overlap from pairwise misclassification
probabilities returned by
[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html),
weighted by \\p_i\Omega\_{ij}+p_j\Omega\_{ji}\\ and averaged over pairs.

## Usage

``` r
compute_average_overlap(true_theta, k = length(true_theta$p))
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` as in a GMM
  \\(\boldsymbol{p},\boldsymbol{\mu},\\\boldsymbol{\Sigma}\_j\\)\\.

- k:

  Number of components \\J\\ (default `length(true_theta$p)`).

## Value

Scalar average overlap.
