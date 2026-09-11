# MixSim-style pairwise overlap via stratified Sobol Monte Carlo

Estimates the MixSim Omega map \\\Omega\_{j\ell}=\Pr\_{X\sim
f_j}(\pi\_\ell f\_\ell(X)\>\pi_j f_j(X))\\ by drawing `n_mc` quasi-Monte
Carlo points **from each component** (inverse-transform sampling of a
Sobol sequence through the Cholesky factor) and comparing
**log-densities**. The average pairwise overlap `BarOmega` is the
unweighted mean of \\\Omega\_{j\ell}+\Omega\_{\ell j}\\ over
\\j\<\ell\\, matching
[`MixSim::overlap()`](https://rdrr.io/pkg/MixSim/man/overlap.html) /
[FSDA MixSim](https://rosa.unipr.it/FSDA/MixSim.html).

## Usage

``` r
overlap_gaussian_mc(true_theta, n_mc = 10000L, seed = NULL, J = NULL)
```

## Arguments

- true_theta:

  List with `p`, `mu`, `sigma` as in
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md).

- n_mc:

  Draws **per component** (default `10000`).

- seed:

  Optional RNG seed (passed to Sobol randomisation and to the uniform
  fallback).

- J:

  Optional number of components.

## Value

A list with `BarOmega`, `MaxOmega`, and `OmegaMap`.

## Details

For two equal-weight components this pairwise overlap equals the
histogram similarity \\\int\min(f_j,f\_\ell)\\ and therefore
\\1-\mathrm{TV}(f_j,f\_\ell)\\ (Nielsen and Sun 2018) .

## References

Nielsen F, Sun K (2018). “Guaranteed Deterministic Bounds on the Total
Variation Distance between Univariate Mixtures.” arXiv:1806.11311.
[doi:10.48550/arxiv.1806.11311](https://doi.org/10.48550/arxiv.1806.11311)
. <https://arxiv.org/abs/1806.11311>.

## See also

[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md),
[`compute_average_riemannian()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_riemannian.md)

## Examples

``` r
set.seed(1)
theta <- list(
  p = c(0.5, 0.5),
  mu = cbind(c(0, 0), c(3, 0)),
  sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
)
overlap_gaussian_mc(theta, n_mc = 2000L)$BarOmega
#> [1] 0.1335
```
