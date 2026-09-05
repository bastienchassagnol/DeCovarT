# Describe a Gaussian-convolution simulation scenario

Summarises the generating law
\\\boldsymbol{y}\mid(\boldsymbol{\zeta},\boldsymbol{p})\sim
\mathcal{N}\_{G}(\boldsymbol{\mu}\boldsymbol{p}, \sum_j
p_j^2\boldsymbol{\Sigma}\_j)\\ at five layers: composition, mean
geometry, covariance / SPD diagnostics, network sparsity, and tangent
Fisher information (mean versus covariance split). MixSim BarOmega and
average pairwise Hellinger of the component Gaussians are kept in
`descriptors`. Jeffreys / symmetrised KL is returned in `supplementary`.

## Usage

``` r
describe_simulation_scenario(true_theta, adjacency = NULL, active_tol = 1e-08)
```

## Arguments

- true_theta:

  List with `p`, `mu`, and `sigma` (and optionally `Theta` precision).
  Parsed by
  [`check_true_theta()`](https://bastienchassagnol.github.io/DeCovarT/reference/check_true_theta.md)
  / `.parse_true_theta()`.

- adjacency:

  Optional \\G\times G\times J\\ adjacency array, or a length-\\J\\ list
  of adjacencies, used for network density when supplied.

- active_tol:

  Threshold for counting an active simplex component.

## Value

A list with:

- `theta_true`: the convolution parameters `p`, `mu`, `sigma`;

- `descriptors`: one-row tibble of kept scenario statistics in six
  families (composition, mean geometry, SPD of
  \\\boldsymbol{\Sigma}(\boldsymbol{p})\\, tangent Fisher, network,
  component overlap). SPD columns include both
  \\\kappa\\\boldsymbol{\Sigma}(\boldsymbol{p})\\\\ (`kappa_sigma_p`)
  and the reciprocal \\\lambda\_{\min}/\lambda\_{\max}\\
  (`kappa_sigma_reciprocal`);

- `supplementary`: one-row tibble of Jeffreys / symmetrised KL, recorded
  but not treated as a primary score;

- `call`: the matched call
  ([`match.call()`](https://rdrr.io/r/base/match.call.html)).

## See also

[`run_simulation_benchmark()`](https://bastienchassagnol.github.io/DeCovarT/reference/run_simulation_benchmark.md),
[`expected_fisher_unconstrained()`](https://bastienchassagnol.github.io/DeCovarT/reference/expected_fisher_unconstrained.md),
[`compute_shannon_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/compute_shannon_entropy.md),
[`composition_from_entropy()`](https://bastienchassagnol.github.io/DeCovarT/reference/composition_from_entropy.md),
[`helmert_basis()`](https://bastienchassagnol.github.io/DeCovarT/reference/helmert_basis.md)

## Examples

``` r
theta <- list(
  p = c(0.5, 0.5),
  mu = cbind(c(10, 0), c(0, 10)),
  sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
)
out <- describe_simulation_scenario(theta)
out$descriptors$h_star
#> [1] 1
out$call
#> describe_simulation_scenario(true_theta = theta)
```
