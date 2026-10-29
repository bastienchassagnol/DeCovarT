# OpenBLAS 0.3.20 selects a faulty kernel on some AVX-512 / AMX CPUs
# (Intel Xeon w9-3575X): solve() is wrong from n = 12, chol() from n = 40.
# Every likelihood in the package relies on both. Setting
# OPENBLAS_CORETYPE=Haswell before R starts avoids it. SkylakeX fixes
# chol() and solve() but not every dsyevr call.

test_that("BLAS / LAPACK return exact Cholesky factors and inverses", {
  for (n in c(12L, 20L, 64L)) {
    a <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
    r <- chol(a)
    expect_lt(max(abs(crossprod(r) - a)), 1e-10)
    expect_lt(max(abs(solve(a) %*% a - diag(n))), 1e-8)
    expect_lt(max(abs(chol2inv(r) %*% a - diag(n))), 1e-8)
  }
})
