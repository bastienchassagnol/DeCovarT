# Monte Carlo interval for an empirical coverage probability

Interval for a binomial coverage rate \\\hat\pi=X/N\\ from independent
coverage indicators. Wilson is the default because the parameter is
bounded in \\\[0,1\]\\ and Wald intervals behave poorly near 0 or 1.

## Usage

``` r
coverage_mc_interval(
  covered,
  conf_level = 0.95,
  method = c("wilson", "wald", "agresti_coull")
)
```

## Arguments

- covered:

  Logical or 0/1 coverage indicators (one per replicate).

- conf_level:

  Confidence level for the interval around \\\hat\pi\\ (default `0.95`).

- method:

  `"wilson"` (default), `"wald"`, or `"agresti_coull"`.

## Value

A list with `n`, `successes`, `coverage`, `mcse`, `lower`, `upper`, and
`method`.

## Examples

``` r
coverage_mc_interval(c(TRUE, TRUE, TRUE, FALSE))
#> $n
#> [1] 4
#> 
#> $successes
#> [1] 3
#> 
#> $coverage
#> [1] 0.75
#> 
#> $mcse
#> [1] 0.2165064
#> 
#> $lower
#> [1] 0.3006418
#> 
#> $upper
#> [1] 0.9544127
#> 
#> $method
#> [1] "wilson"
#> 
```
