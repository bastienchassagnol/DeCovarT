# Equicorrelation (constant-correlation) Gram matrix

\\R=(1-\rho)I+\rho\mathbf{1}\mathbf{1}^{\mathsf{T}}\\. Positive
semidefinite for \\\rho\in\[-1/(J-1),1\]\\.

## Usage

``` r
equicorrelation_gram(n_celltypes, target_cosine = 0)
```

## Arguments

- n_celltypes:

  Integer \\J\ge 2\\.

- target_cosine:

  Pairwise cosine \\\rho\\.

## Value

Symmetric \\J\times J\\ matrix with unit diagonal.

## Examples

``` r
equicorrelation_gram(3L, 0.4)
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.4  0.4
#> [2,]  0.4  1.0  0.4
#> [3,]  0.4  0.4  1.0
```
