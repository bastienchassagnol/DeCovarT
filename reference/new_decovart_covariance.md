# Structure-aware covariance backend for \\\boldsymbol{\Sigma}(\boldsymbol{p})\\

Wraps the cell-type covariance array together with a declared
`structure`, so that the bulk covariance
\\\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j
p_j^2\boldsymbol{\Sigma}\_j\\ can be factorised with the cheapest exact
method. The object exposes operators via
[`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md),
[`sigma_solve()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_solve.md),
[`sigma_quadform()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_quadform.md)
and
[`sigma_trace_precision_times()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_trace_precision_times.md);
the dense array is always accepted and wrapped as `"dense"`.

## Usage

``` r
new_decovart_covariance(
  Sigma = NULL,
  structure = c("dense", "block", "band", "sparse", "diag_lowrank"),
  blocks = NULL,
  bandwidth = NULL,
  diagonal = NULL,
  loadings = NULL,
  core = NULL,
  tol = 1e-08
)

# S3 method for class 'decovart_covariance'
print(x, ...)
```

## Arguments

- Sigma:

  Numeric array of cell-type covariances in \\\mathcal{M}\_{G\times
  G\times J}\\. For `"diag_lowrank"` it may be omitted and is then
  reconstructed from `diagonal`, `loadings`, `core`.

- structure:

  One of `"dense"`, `"block"`, `"band"`, `"sparse"`, `"diag_lowrank"`.
  Matching is case-insensitive.

- blocks:

  For `"block"`: an integer / factor vector of length \\G\\ giving block
  membership, or a list of index vectors.

- bandwidth:

  For `"band"`: non-negative integer bandwidth \\b\\.

- diagonal:

  For `"diag_lowrank"`: numeric \\G\times J\\ matrix whose column \\j\\
  is the diagonal of \\\mathbf{D}\_j\\.

- loadings:

  For `"diag_lowrank"`: shared loadings
  \\\mathbf{U}\in\mathcal{M}\_{G\times r}\\.

- core:

  For `"diag_lowrank"`: array \\\mathcal{M}\_{r\times r\times J}\\ of
  cores \\\mathbf{C}\_j\\.

- tol:

  Relative tolerance for validating that off-structure entries vanish
  (block / band).

- x:

  A `decovart_covariance` object (print method).

- ...:

  Ignored.

## Value

An object of class `"decovart_covariance"`.

## Details

Supported structures:

- `"dense"`: universal fallback; one Cholesky of the assembled matrix.

- `"block"`: exactly block-diagonal covariance (disconnected gene
  modules); one Cholesky per block, so the factorisation cost drops by
  roughly the squared number of equal blocks.

- `"band"`: covariance with bandwidth \\b\\ (\\\Sigma\_{jk}=0\\ for
  \\\|j-k\|\>b\\); a banded / sparse Cholesky.

- `"sparse"`: sparse covariance with a fixed nonzero pattern across
  \\\boldsymbol{p}\\; the fill-reducing ordering (symbolic
  factorisation) is computed once and reused for every numeric refactor.

- `"diag_lowrank"`: \\\boldsymbol{\Sigma}\_j=\mathbf{D}\_j+
  \mathbf{U}\mathbf{C}\_j\mathbf{U}^{\mathsf{T}}\\ with shared loadings
  \\\mathbf{U}\\; solves use the Woodbury identity and the
  log-determinant uses the matrix-determinant lemma, transferring the
  \\O(G^3)\\ cost to the low rank \\r\ll G\\.

Band and sparse factorisations use Matrix (CHOLMOD; (Chen et al. 2008)
). Dense Cholesky complexity follows (Golub and Van Loan 2013) ; the
low-rank path uses the Woodbury identity (Hager 1989) .

## References

Chen Y, Davis TA, Hager WW, Rajamanickam S (2008). “Algorithm 887:
CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate.”
*ACM Transactions on Mathematical Software*, **35**(3), 22:1–22:14.
[doi:10.1145/1391989.1391995](https://doi.org/10.1145/1391989.1391995) .

Golub GH, Van Loan CF (2013). *Matrix Computations*. JHU Press. ISBN
978-1-4214-0859-0.

Hager WW (1989). “Updating the Inverse of a Matrix.” *SIAM Review*,
**31**(2), 221–239.
[doi:10.1137/1031049](https://doi.org/10.1137/1031049) .

## See also

[`sigma_logdet()`](https://bastienchassagnol.github.io/DeCovarT/reference/sigma_logdet.md),
[`covariance_structure_from_graph_model()`](https://bastienchassagnol.github.io/DeCovarT/reference/covariance_structure_from_graph_model.md)

## Examples

``` r
Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
cov_dense <- new_decovart_covariance(Sigma, "dense")
sigma_logdet(cov_dense, c(0.6, 0.4))
#> [1] -1.307853
```
