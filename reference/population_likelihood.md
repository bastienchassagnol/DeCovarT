# Likelihood objects for the benchmark

`population_likelihood()` is the Gaussian convolution
\\\mathcal{N}(\boldsymbol{\mu}\boldsymbol{p}, \sum_j
p_j^2\boldsymbol{\Sigma}\_j)\\; `single_cell_likelihood()` is the
mixture of cells, with its two scales \\(c,\phi)\\ profiled out, so
every function is evaluated at the closed-form scales of the
composition. They play the part of
[`stats::family()`](https://rdrr.io/r/stats/family.html) in
[`stats::glm()`](https://rdrr.io/r/stats/glm.html):
[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md)
and
[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)
take one as `likelihood`, either the object or the constructor itself.

## Usage

``` r
population_likelihood()

single_cell_likelihood()
```

## Value

An object of class `decovart_likelihood`.

## Details

Each object is a list of class `decovart_likelihood` with

- `name`: `"population"` or `"single_cell"`;

- `nuisance`: whether the model has the scales \\(c,\phi)\\;

- `parameters(cell_names)`: the rows of every interval matrix, the cell
  types and, for the single-cell model, `mean_scale`, `variance_scale`,
  `platform_scale_a` and `cell_count_n`;

- `simulate(mean, Sigma, theta, n)`: `n` bulk samples drawn from the
  moments `mean` and `Sigma` with the composition and settings in
  `theta`; returns `Y`, the target `truth` and the true scales
  `nuisance` (`NULL` for the population model);

- `loglik(p, y, mu, Sigma)`, `score(p, y, mu, Sigma)`,
  `se(p, y, mu, Sigma)` (Wald standard errors of the composition) and
  `curvature(p, y, mu, Sigma)` (local-maximum flag and eigenvalue range
  of the Hessian in ILR coordinates, `NA` on the boundary);

- `interval(method, y, p_hat, p_true, mu, Sigma, ci)`: one method of the
  interval menu at the estimate `p_hat`, as a two-column matrix.

The interval methods are

- `wald`: delta method on \\\boldsymbol{p}\\, symmetric;

- `ilr`: Wald interval in ILR coordinates mapped back to the simplex;

- `profile`: likelihood-ratio inversion; it refits from `y`, so it does
  not depend on the solver that produced `p_hat`;

- `bootstrap`: percentile interval of a parametric bootstrap refitted by
  `ci$optimiser`;

- `oracle`: `p_hat` plus or minus the Wald half-width evaluated at the
  true composition.

For the single-cell model, `_reml` appended to `wald`, `ilr` or
`profile` applies the REML calibration (\\\nu=G-J\\), and Wald, ILR and
bootstrap intervals also cover the scales and the derived \\a=\phi/c\\
and \\n=c^{2}/\phi\\.

Wald, ILR and oracle intervals are undefined when `p_hat` has a zero
coordinate: on a face (one zero), an edge (two zeros) or a vertex
(\\J-1\\ zeros, one cell type making up the whole sample). L-BFGS-B
returns such estimates; Newton-Raphson in ILR coordinates does not.
Their rows are then `NA`, scales included, and the coverage of these
methods is conditional on an interior estimate.

## See also

[`deconvolute_ratios()`](https://bastienchassagnol.github.io/DeCovarT/reference/deconvolute_ratios.md),
[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md)

## Examples

``` r
mu <- matrix(c(20, 22, 22, 20), 2, dimnames = list(NULL, c("a", "b")))
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
lik <- population_likelihood()
lik
#> <decovart_likelihood: population>
y <- drop(mu %*% c(0.4, 0.6))
lik$loglik(c(0.4, 0.6), y, mu, Sigma)
#> [1] 0.6539265
kkt_residual(c(0.4, 0.6), lik$score(c(0.4, 0.6), y, mu, Sigma))
#> [1] 0.5439283
single_cell_likelihood()$parameters(c("a", "b"))
#> [1] "a"                "b"                "mean_scale"       "variance_scale"  
#> [5] "platform_scale_a" "cell_count_n"    
```
