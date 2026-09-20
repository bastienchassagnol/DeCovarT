# Q-Q plot of whitened ILR coordinates versus N(0, 1)

Tail-sensitive simultaneous envelope (`qqplotr` `bandType = "ts"`)
calibrated to the theoretical \\N(0,1)\\ law (`identity = TRUE`,
`mu = 0`, `sigma = 1`). Colour and fill distinguish solvers. The plot is
**not** detrended.

## Usage

``` r
plot_mc_qq_normal(df, title = NULL)
```

## Arguments

- df:

  Long table with `whitened`, `algorithm`, and `panel`.

- title:

  Page title.

## Value

A `ggplot`.
