# Calibrate a global covariance scale to a target MixSim BarOmega

Binary-searches \\s\>0\\ so that
[`compute_average_overlap()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_average_overlap.md)
of \\\\s\Sigma_j\\\\ matches `target`. Larger \\s\\ inflates the
Gaussians and **increases** average overlap while leaving graph zeros
intact. This is the graph-constrained analogue of MixSim / FSDA's
average-overlap generator, which cannot take a precision support as
input.

## Usage

``` r
scale_covariances_to_overlap(
  mu,
  sigma,
  p,
  target,
  n_mc = 4000L,
  s_range = c(0.001, 1000),
  tol = 0.005,
  seed = NULL,
  verbose = FALSE,
  pairs = NULL
)
```

## Arguments

- mu:

  \\G\times J\\ mean signature.

- sigma:

  Base \\G\times G\times J\\ covariances.

- p:

  Simplex weights (length \\J\\).

- target:

  Target average pairwise overlap in \\(0,1)\\.

- n_mc:

  Monte Carlo size forwarded to
  [`overlap_gaussian_mc()`](https://bastienchassagnol.github.io/DeCovarT/reference/overlap_gaussian_mc.md)
  when \\G\ge 4\\.

- s_range:

  Length-2 search interval for \\s\\.

- tol:

  Absolute overlap tolerance.

- seed:

  Optional seed for overlap evaluations.

- verbose:

  If `TRUE`, print the MixSim-to-MC switch.

- pairs:

  Optional two-column matrix of cell-type indices. When given, the
  calibrated quantity is the mean overlap
  \\\omega\_{j\ell}+\omega\_{\ell j}\\ over these pairs only (for
  example the sibling types of each lineage), instead of the average
  over all pairs.

## Value

A list with `sigma`, `Theta`, `scale`, `baromega` (the calibrated
overlap at the returned scale), and `target`.
