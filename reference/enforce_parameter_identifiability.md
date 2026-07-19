# Order GMM parameters for unique labelling

Lexicographically orders mixture components by columns of
\\\boldsymbol{\mu}\\ and renormalises \\\boldsymbol{p}\\.

## Usage

``` r
enforce_parameter_identifiability(theta)
```

## Arguments

- theta:

  List with elements `p` (\\\boldsymbol{p}\\), `mu`
  (\\\boldsymbol{\mu}\\), and `sigma` (array of
  \\\boldsymbol{\Sigma}\_j\\).

## Value

Relabelled list with the same structure.
