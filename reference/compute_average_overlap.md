# Average pairwise overlap of a Gaussian mixture

Returns MixSim's **BarOmega**: the unweighted mean of pairwise overlaps
\$\$ \overline{\omega} = \frac{2}{J(J-1)} \sum\_{1\le j\<\ell\le J}
\bigl(\Omega\_{j\ell}+\Omega\_{\ell j}\bigr) \in\[0,1\] \$\$ (up to the
MixSim numerical convention), where \\\Omega\_{j\ell}=\Pr\_{X\sim
f_j}(X\text{ classified as }\ell)\\ already uses the mixture weights
\\\boldsymbol{p}\\ inside the Bayes / MAP rule of
[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html). Do
**not** multiply the directional masses by \\p_j\\ again.

## Usage

``` r
compute_average_overlap(true_theta, J = NULL)
```

## Arguments

- true_theta:

  List validated by
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md):
  `p` (length \\J\\ or \\J\times N\\), `mu` (\\G\times J\\), `sigma`
  (\\G\times G\times J\\).

- J:

  Number of cell types (components). Defaults to the third dimension of
  `sigma`.

## Value

Scalar average pairwise overlap (MixSim `BarOmega`).

## See also

[`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md)

## Examples

``` r
set.seed(1)
theta <- list(
  p = c(0.5, 0.5),
  mu = cbind(c(0, 0), c(3, 0)),
  sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
)
compute_average_overlap(theta)
#> [1] 0.1336144
```
