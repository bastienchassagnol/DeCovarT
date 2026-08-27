# Chi-bar-square tail probability for boundary likelihood-ratio tests

Upper-tail probability of the mixture \$\$ \bar{\chi}^{2} =
\sum\_{i=0}^{q} w_i\\\chi^{2}\_{s+i}, \qquad w_i=\binom{q}{i}2^{-q},
\$\$ with \\q\\ constraints active on the boundary of the parameter
space and \\s\\ interior constraints. This is Case 9 of (Self and Liang
1987) , which generalises the
\\\tfrac{1}{2}\chi^{2}\_{0}+\tfrac{1}{2}\chi^{2}\_{1}\\ result of
(Chernoff 1954) . For \\q=0\\ the mixture collapses to the ordinary
\\\chi^{2}\_{s}\\ of (Wilks 1938) .

## Usage

``` r
chi_bar_square_pvalue(
  statistic,
  n_boundary = 1L,
  df_interior = 0L,
  weights = NULL
)
```

## Arguments

- statistic:

  Numeric vector of observed likelihood-ratio statistics
  \\D=2(\hat{\ell}\_1-\hat{\ell}\_0)\\.

- n_boundary:

  Integer \\q\\: number of constraints active on the boundary (e.g. null
  ratios set to zero).

- df_interior:

  Integer \\s\\: number of interior (two-sided) constraints tested
  simultaneously.

- weights:

  Optional numeric vector of length `n_boundary + 1` giving the mixing
  weights \\w_0,\ldots,w_q\\. Defaults to the binomial weights.

## Value

Numeric vector of upper-tail probabilities.

## Details

The binomial weights hold when the block of the Fisher information
corresponding to the active constraints is diagonal after
orthogonalisation against the nuisance block. Under correlated
constraints the mixing weights depend on the geometry of the
approximating tangent cone and the binomial choice is only an
approximation; Self and Liang also exhibit a nuisance-parameter
configuration whose limit is **not** a chi-square mixture. Supply
`weights` explicitly, or calibrate by simulation with
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
whenever \\q\>1\\ and the constrained cell types share signal.

## Atom at zero

\\\chi^{2}\_{0}\\ is a point mass at the origin, so its strict upper
tail \\\Pr(\chi^{2}\_{0}\>D)\\ is zero for every \\D\ge 0\\. The
component therefore only removes weight: with one active constraint and
no interior constraint the p-value is
\\\tfrac{1}{2}\Pr(\chi^{2}\_{1}\>D)\\, capped at \\0.5\\ when \\D=0\\.
Note [`stats::pchisq()`](https://rdrr.io/r/stats/Chisquare.html)
switches to the *closed* tail at exactly `q = 0` for `df = 0`; the
strict convention is used here so the p-value stays continuous in
`statistic`.

## References

Self SG, Liang K (1987). “Asymptotic Properties of Maximum Likelihood
Estimators and Likelihood Ratio Tests under Nonstandard Conditions.”
*Journal of the American Statistical Association*, **82**(398), 605–610.
ISSN 0162-1459.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)
.

Chernoff H (1954). “On the Distribution of the Likelihood Ratio.” *The
Annals of Mathematical Statistics*, **25**(3), 573–578. ISSN 0003-4851.
[doi:10.1214/aoms/1177728725](https://doi.org/10.1214/aoms/1177728725) .

Wilks SS (1938). “The Large-Sample Distribution of the Likelihood Ratio
for Testing Composite Hypotheses.” *The Annals of Mathematical
Statistics*, **9**(1), 60–62. ISSN 0003-4851.
[doi:10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360) .

## See also

[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md)

Other decovart_inference:
[`bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/bootstrap_decovart.md),
[`boundary_diagnostics()`](https://bastienchassagnol.github.io/DeCovarT/reference/boundary_diagnostics.md),
[`confint_profile_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/confint_profile_decovart.md),
[`equivariance_check_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/equivariance_check_decovart.md),
[`lrt_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/lrt_decovart.md),
[`multistart_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/multistart_decovart.md),
[`profile_loglik_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/profile_loglik_decovart.md),
[`reference_bootstrap_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/reference_bootstrap_decovart.md),
[`restricted_mle_decovart()`](https://bastienchassagnol.github.io/DeCovarT/reference/restricted_mle_decovart.md)

## Examples

``` r
# One active boundary constraint: 50:50 mixture of chi2_0 and chi2_1.
chi_bar_square_pvalue(2.71, n_boundary = 1L)
#> [1] 0.0498605
# Interior (Wilks) calibration.
chi_bar_square_pvalue(3.84, n_boundary = 0L, df_interior = 1L)
#> [1] 0.05004352
```
