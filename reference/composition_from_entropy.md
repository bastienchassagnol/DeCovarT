# One-dominant composition with a target normalised Shannon entropy

Returns a length-\\J\\ simplex vector of the form
\\(1-(J-1)q,\\q,\ldots,q)\\ whose Pielou evenness
[`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)
equals `h_star`. The uniform composition (`h_star = 1`) and a Dirac mass
(`h_star = 0`) are returned exactly.

## Usage

``` r
composition_from_entropy(h_star, n_celltypes = 3L, nms = NULL)
```

## Arguments

- h_star:

  Target \\H^{\star}\in\[0,1\]\\.

- n_celltypes:

  Integer \\J\ge 2\\.

- nms:

  Optional names for the returned vector.

## Value

Numeric simplex vector of length \\J\\.

## See also

[`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md)

## Examples

``` r
p <- composition_from_entropy(0.5, n_celltypes = 3L)
compute_shannon_entropy(p)
#> [1] 0.5
```
