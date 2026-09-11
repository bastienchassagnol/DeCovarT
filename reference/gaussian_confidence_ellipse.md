# Exact Gaussian probability ellipse (known mean and covariance)

Boundary of the set
\\(x-\mu)^{\mathsf{T}}\Sigma^{-1}(x-\mu)\le\chi^2\_{2,\alpha}\\ for a
bivariate Gaussian with **known** \\\mu\\ and \\\Sigma\\. The squared
Mahalanobis distance is exactly \\\chi^2_2\\.

## Usage

``` r
gaussian_confidence_ellipse(mu, Sigma, level = 0.95, n = 200L)
```

## Arguments

- mu:

  Length-2 mean.

- Sigma:

  \\2 \times 2\\ covariance.

- level:

  Probability content (default 0.95).

- n:

  Number of boundary points.

## Value

A data frame with `x` and `y`.

## Examples

``` r
mu <- c(20, 22)
Sigma <- matrix(c(1, 0.4, 0.4, 1), 2)
head(gaussian_confidence_ellipse(mu, Sigma))
#>          x        y
#> 1 22.39611 22.50010
#> 2 22.41071 22.57549
#> 3 22.42290 22.65031
#> 4 22.43267 22.72448
#> 5 22.44002 22.79793
#> 6 22.44494 22.87058
```
