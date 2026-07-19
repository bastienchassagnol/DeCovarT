# Normalised Shannon entropy of a discrete distribution

For \\\boldsymbol{p}\\ on the simplex (after dropping zeros), \$\$
H(\boldsymbol{p}) = -\sum\_{j}p_j\log\_{J'}p_j \in\[0,1\], \$\$ with
\\J'\\ the number of positive masses (1 = uniform).

## Usage

``` r
compute_shannon_entropy(ratios)
```

## Arguments

- ratios:

  Numeric vector \\\boldsymbol{p}\\.

## Value

Scalar entropy in \\\[0,1\]\\.
