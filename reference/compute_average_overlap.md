# Average pairwise overlap of a Gaussian mixture

Returns MixSim's **BarOmega**: the unweighted mean of pairwise overlaps
\$\$ \overline{\omega} = \frac{2}{J(J-1)} \sum\_{1\le j\<\ell\le J}
\bigl(\Omega\_{j\ell}+\Omega\_{\ell j}\bigr) \in\[0,1\] \$\$ (up to the
MixSim numerical convention), where \\\Omega\_{j\ell}=\Pr\_{X\sim
f_j}(X\text{ classified as }\ell)\\ already uses the mixture weights
\\\boldsymbol{p}\\ inside the Bayes / MAP rule of
[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html). Do
**not** multiply the directional masses by \\p_j\\ again.

For \\G\ge 4\\ MixSim's Davies quadrature becomes expensive. The helper
then switches to
[`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md):
stratified Sobol draws, inverse-transform sampling through precomputed
Cholesky factors, and log-density MAP comparisons (`n_mc` draws per
component; default 10,000).

## Usage

``` r
compute_average_overlap(
  true_theta,
  J = NULL,
  n_mc = 10000L,
  seed = NULL,
  verbose = TRUE
)
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

- n_mc:

  Monte Carlo draws per component when \\G\ge 4\\.

- seed:

  Optional seed forwarded to
  [`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md).

- verbose:

  If `TRUE` (default), announce the MixSim-to-MC switch with `cli` when
  it is installed.

## Value

Scalar average pairwise overlap (MixSim `BarOmega`).

## See also

[`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md),
[`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md)

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
