# Complete a weighted adjacency to an SPD precision

Applies a uniform spectral shift that preserves off-diagonal signs and
support: \$\$ \boldsymbol{\Omega} = \boldsymbol{W} +
\bigl(\lvert\lambda\_{\min}(\boldsymbol{W})\rvert + u\bigr) \mathbf{I}.
\$\$

## Usage

``` r
build_normalised_precision(weighted_adjacency, precision_shift)
```

## Arguments

- weighted_adjacency:

  Symmetric numeric matrix (zero diagonal).

- precision_shift:

  Positive diagonal cushion \\u\\.

## Value

Symmetric positive-definite matrix.
