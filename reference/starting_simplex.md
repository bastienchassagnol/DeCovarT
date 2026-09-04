# Open-simplex start for ILR solvers

The convolution log-likelihood is not globally concave, so the start can
change which basin the solver enters. Draw several independent Dirichlet
starts (different RNG streams or
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md))
rather than a single draw when probing multimodality.

## Usage

``` r
starting_simplex(
  n_celltypes,
  initial_p = NULL,
  nms = NULL,
  y = NULL,
  mean_signature_matrix = NULL,
  dirichlet_alpha = 1
)
```

## Arguments

- n_celltypes:

  Integer \\J\\.

- initial_p:

  One of:

  - `NULL` or `"barycentre"`: equi-balanced \\(1/J,\ldots,1/J)\\;

  - a numeric vector of length \\J\\: used as-is after
    [`repair_simplex()`](https://bastienchassagnol.github.io/DeCovarT/reference/repair_simplex.md);

  - `"dirichlet"`: one Dirichlet\\(\alpha,\ldots,\alpha)\\ draw (see
    `dirichlet_alpha`);

  - `"qp"` (aliases `"deconrnaseq"`, `"lsei"`): mean-only simplex QP
    from
    [`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)
    (`y` and `mean_signature_matrix` required).

- nms:

  Optional names (cell-type colnames).

- y, mean_signature_matrix:

  Bulk and signature, required for a QP start.

- dirichlet_alpha:

  Positive concentration, recycled to length \\J\\. The default `1` is
  uniform on the simplex. \\\alpha\>1\\ concentrates mass near the
  barycentre; \\\alpha\<1\\ puts extra mass near faces (boundary-biased
  restarts).

## Value

A length-\\J\\ open-simplex vector.

## See also

[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`deconvolute_ratios_deconrnaseq()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios_Marquardt_Levenberg.md)

## Examples

``` r
starting_simplex(3L)
#> [1] 0.3333333 0.3333333 0.3333333
set.seed(1)
starting_simplex(3L, "dirichlet")
#>        ct1        ct2        ct3 
#> 0.04037978 0.48994649 0.46967373 
```
