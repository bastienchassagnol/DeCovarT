# Numerically stable softmax and additive log-ratio derivatives

This vignette collects the closed-form derivatives that DeCovarT needs
to move between the unconstrained coordinates
$`\boldsymbol{\rho}\in\mathbb{R}^{J-1}`$ and the cellular ratios
$`\boldsymbol{p}\in\Delta^{J-1}`$. For every map we give the explicit
tensor formula, an efficient base-`R` implementation (evaluated and
checked against
[`numDeriv`](https://cran.r-project.org/package=numDeriv)), and a
batched `PyTorch` counterpart. The naming follows the DeCovarT API:
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
for $`\boldsymbol{\rho}\mapsto\boldsymbol{p}`$ and
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
for $`\boldsymbol{p}\mapsto\boldsymbol{\rho}`$.

## Log-sum-exp and softmax

The softmax of a score vector $`\boldsymbol{z}\in\mathbb{R}^{K}`$ is
evaluated through the log-sum-exp (LSE) to avoid overflow,

``` math
s_i=\operatorname{softmax}(\boldsymbol{z})_i=\frac{e^{z_i}}{\sum_{k}e^{z_k}}
=e^{z_i-\operatorname{LSE}(\boldsymbol{z})},
\qquad
\operatorname{LSE}(\boldsymbol{z})=m+\log\!\sum_{k}e^{z_k-m},
\quad m=\max_k z_k.
 \qquad(1)
```

``` r

log_sum_exp <- function(z) {
  m <- max(z)
  m + log(sum(exp(z - m)))
}

softmax <- function(z) {
  exp(z - log_sum_exp(z))
}
```

Listing 1: Stable log-sum-exp and softmax in R

``` python
import torch
import torch.nn.functional as F


def softmax(z: torch.Tensor, dim: int = -1) -> torch.Tensor:
    return F.softmax(z, dim=dim)


def log_sum_exp(z: torch.Tensor, dim: int = -1) -> torch.Tensor:
    return torch.logsumexp(z, dim=dim)
```

Listing 2: Softmax and log-sum-exp in PyTorch

### Softmax Jacobian and Hessian

The softmax Jacobian is the difference between a diagonal and a rank-one
term,

``` math
\frac{\partial s_i}{\partial z_j}=s_i(\delta_{ij}-s_j),
\qquad
\mathbf{J}_{\mathrm{softmax}}=\operatorname{diag}(\boldsymbol{s})
-\boldsymbol{s}\boldsymbol{s}^{\top}.
 \qquad(2)
```

This matrix also equals the Hessian of the log-sum-exp,
$`\nabla^2\operatorname{LSE}(\boldsymbol{z})=\operatorname{diag}(\boldsymbol{s})
-\boldsymbol{s}\boldsymbol{s}^{\top}`$. Because softmax is
vector-valued, its second derivative is a third-order tensor
$`H_{ijk}=\partial^2 s_i/(\partial z_j\partial z_k)`$,

``` math
H_{ijk}=s_i\bigl[(\delta_{ij}-s_j)(\delta_{ik}-s_k)-s_j(\delta_{jk}-s_k)\bigr],
 \qquad(3)
```

which, slice by slice, is the rank-one correction

``` math
\mathbf{H}^{(i)}=s_i\Bigl[(\boldsymbol{e}_i-\boldsymbol{s})
(\boldsymbol{e}_i-\boldsymbol{s})^{\top}
-\bigl(\operatorname{diag}(\boldsymbol{s})-\boldsymbol{s}\boldsymbol{s}^{\top}\bigr)
\Bigr],
 \qquad(4)
```

with $`\boldsymbol{e}_i`$ the $`i`$-th canonical basis vector.
[Listing 3](#lst-softmax-deriv-r) materialises [Eq. 2](#eq-softmax-jac)
and [Eq. 4](#eq-softmax-hess-slice) in $`O(K^2)`$ per slice.

``` r

softmax_jacobian <- function(z) {
  s <- softmax(z)
  diag(s) - tcrossprod(s)
}

softmax_hessian <- function(z) {
  s <- softmax(z)
  k <- length(s)
  jac <- diag(s) - tcrossprod(s)
  hessian <- array(0, dim = c(k, k, k))
  for (i in seq_len(k)) {
    di <- (seq_len(k) == i) - s # e_i - s
    hessian[i, , ] <- s[i] * (tcrossprod(di) - jac)
  }
  hessian
}
```

Listing 3: Softmax Jacobian and Hessian tensor in R

``` python
def softmax_jacobian(z: torch.Tensor) -> torch.Tensor:
    """z: [..., K] -> [..., K, K]."""
    s = softmax(z, dim=-1)
    return torch.diag_embed(s) - s.unsqueeze(-1) * s.unsqueeze(-2)


def softmax_hessian(z: torch.Tensor) -> torch.Tensor:
    """z: [K] -> H[i, j, k] of shape [K, K, K]."""
    s = softmax(z, dim=0)
    eye = torch.eye(s.numel(), dtype=s.dtype, device=s.device)
    diff = eye - s  # (e_i - s) stacked over i
    outer = diff.unsqueeze(-1) * diff.unsqueeze(-2)
    jac = torch.diag(s) - torch.outer(s, s)
    return s[:, None, None] * (outer - jac.unsqueeze(0))
```

Listing 4: Batched softmax Jacobian and Hessian in PyTorch

## Additive log-ratio coordinates

Let $`\boldsymbol{p}=(p_1,\ldots,p_J)`$ be a composition with $`p_j>0`$
and $`\sum_j p_j=1`$, and take the last part $`J`$ as reference. The
*additive log-ratio* (alr) transform sends the composition to
$`\mathbb{R}^{J-1}`$,

``` math
\rho_i=\log\!\Bigl(\frac{p_i}{p_J}\Bigr),\qquad i=1,\ldots,J-1,
 \qquad(5)
```

and its inverse, the *additive logistic* transform
($`\mathrm{alr}^{-1}`$), pins a reference logit at zero,
$`\tilde{\boldsymbol{\rho}}=(\rho_1,\ldots,\rho_{J-1},0)`$, and applies
a softmax,

``` math
\boldsymbol{p}=\operatorname{softmax}(\tilde{\boldsymbol{\rho}}),
\qquad
p_i=\frac{e^{\rho_i}}{1+\sum_{j<J}e^{\rho_j}}\ (i<J),
\qquad
p_J=\frac{1}{1+\sum_{j<J}e^{\rho_j}}.
 \qquad(6)
```

``` r

# rho -> p  (additive logistic, alr^{-1}); reference logit appended as 0
additive_logistic <- function(rho) {
  softmax(c(rho, 0))
}

# p -> rho  (additive log-ratio, alr); last part as reference
additive_log_ratio <- function(p) {
  j <- length(p)
  log(p[-j] / p[j])
}
```

Listing 5: Additive logistic and additive log-ratio maps in R

``` python
def additive_logistic(rho: torch.Tensor) -> torch.Tensor:
    """rho: [..., J-1] -> p: [..., J] (reference category = last)."""
    reference = torch.zeros_like(rho[..., :1])
    return torch.softmax(torch.cat([rho, reference], dim=-1), dim=-1)


def additive_log_ratio(p: torch.Tensor) -> torch.Tensor:
    """p: [..., J] -> rho: [..., J-1]."""
    log_p = torch.log(p)
    return log_p[..., :-1] - log_p[..., -1:]
```

Listing 6: Additive logistic and additive log-ratio maps in PyTorch

### Derivatives of the additive logistic map ($`\boldsymbol{\rho}\mapsto\boldsymbol{p}`$)

The additive logistic map is the softmax restricted to the first $`J-1`$
inputs, so its Jacobian
$`\mathbf{J}_{\boldsymbol{\psi}}\in\mathcal{M}_{J\times(J-1)}`$ is the
softmax Jacobian with the reference column dropped,

``` math
\frac{\partial p_i}{\partial\rho_a}=p_i(\delta_{ia}-p_a),
\qquad
\mathbf{J}_{\boldsymbol{\psi}}
=\bigl[\operatorname{diag}(\boldsymbol{p})
-\boldsymbol{p}\boldsymbol{p}^{\top}\bigr]_{:,\,1:(J-1)},
\qquad a=1,\ldots,J-1.
 \qquad(7)
```

Its Hessian is a tensor
$`\mathbf{H}_{\boldsymbol{\psi}}\in\mathcal{M}_{J\times(J-1)\times(J-1)}`$,
the softmax Hessian [Eq. 4](#eq-softmax-hess-slice) restricted to the
first $`J-1`$ input indices,

``` math
\frac{\partial^2 p_i}{\partial\rho_a\partial\rho_b}
=p_i\bigl[(\delta_{ia}-p_a)(\delta_{ib}-p_b)-p_a(\delta_{ab}-p_b)\bigr],
\qquad a,b=1,\ldots,J-1.
 \qquad(8)
```

``` r

additive_logistic_jacobian <- function(rho) {
  p <- additive_logistic(rho)
  j <- length(p)
  (diag(p) - tcrossprod(p))[, -j, drop = FALSE]
}

additive_logistic_hessian <- function(rho) {
  p <- additive_logistic(rho)
  j <- length(p)
  jac <- diag(p) - tcrossprod(p)
  hessian <- array(0, dim = c(j, j - 1, j - 1))
  for (i in seq_len(j)) {
    di <- (seq_len(j) == i) - p
    slice <- p[i] * (tcrossprod(di) - jac)
    hessian[i, , ] <- slice[-j, -j]
  }
  hessian
}
```

Listing 7: Jacobian and Hessian of the additive logistic map in R

``` python
def additive_logistic_jacobian(rho: torch.Tensor) -> torch.Tensor:
    """rho: [J-1] -> Jacobian [J, J-1]."""
    p = additive_logistic(rho)
    return (torch.diag(p) - torch.outer(p, p))[:, :-1]


def additive_logistic_hessian(rho: torch.Tensor) -> torch.Tensor:
    """rho: [J-1] -> Hessian [J, J-1, J-1]."""
    reference = torch.zeros_like(rho[:1])
    augmented = torch.cat([rho, reference])
    return softmax_hessian(augmented)[:, :-1, :-1]
```

Listing 8: Jacobian and Hessian of the additive logistic map in PyTorch

### Derivatives of the additive log-ratio map ($`\boldsymbol{p}\mapsto\boldsymbol{\rho}`$)

Treating $`\boldsymbol{p}`$ as ambient coordinates in $`\mathbb{R}^{J}`$
(admissible tangent directions satisfy $`\sum_i \mathrm{d}p_i=0`$), the
alr Jacobian $`\mathbf{J}_{\mathrm{alr}}\in\mathcal{M}_{(J-1)\times J}`$
is sparse,

``` math
\frac{\partial\rho_i}{\partial p_j}=\frac{\delta_{ij}}{p_i}-\frac{\delta_{jJ}}{p_J},
\qquad
\mathbf{J}_{\mathrm{alr}}
=\Bigl[\operatorname{diag}\!\bigl(p_1^{-1},\ldots,p_{J-1}^{-1}\bigr)
\ \big|\ -p_J^{-1}\boldsymbol{1}_{J-1}\Bigr],
 \qquad(9)
```

and its Hessian
$`\mathbf{H}_{\mathrm{alr}}\in\mathcal{M}_{(J-1)\times J\times J}`$ is
diagonal in the last two modes,

``` math
\frac{\partial^2\rho_i}{\partial p_j\partial p_k}
=-\frac{\delta_{ij}\delta_{ik}}{p_i^{2}}
+\frac{\delta_{jJ}\delta_{kJ}}{p_J^{2}}.
 \qquad(10)
```

``` r

additive_log_ratio_jacobian <- function(p) {
  j <- length(p)
  idx <- seq_len(j - 1)
  jac <- matrix(0, j - 1, j)
  jac[cbind(idx, idx)] <- 1 / p[idx]
  jac[, j] <- -1 / p[j]
  jac
}

additive_log_ratio_hessian <- function(p) {
  j <- length(p)
  hessian <- array(0, dim = c(j - 1, j, j))
  for (i in seq_len(j - 1)) {
    hessian[i, i, i] <- -1 / p[i]^2
    hessian[i, j, j] <- 1 / p[j]^2
  }
  hessian
}
```

Listing 9: Jacobian and Hessian of the additive log-ratio map in R

``` python
def additive_log_ratio_jacobian(p: torch.Tensor) -> torch.Tensor:
    """p: [J] -> Jacobian [J-1, J]."""
    j = p.numel()
    jac = torch.zeros(j - 1, j, dtype=p.dtype, device=p.device)
    idx = torch.arange(j - 1, device=p.device)
    jac[idx, idx] = 1.0 / p[:-1]
    jac[:, -1] = -1.0 / p[-1]
    return jac


def additive_log_ratio_hessian(p: torch.Tensor) -> torch.Tensor:
    """p: [J] -> Hessian [J-1, J, J]."""
    j = p.numel()
    hessian = torch.zeros(j - 1, j, j, dtype=p.dtype, device=p.device)
    idx = torch.arange(j - 1, device=p.device)
    hessian[idx, idx, idx] = -1.0 / p[:-1].square()
    hessian[:, -1, -1] = 1.0 / p[-1].square()
    return hessian
```

Listing 10: Jacobian and Hessian of the additive log-ratio map in
PyTorch

## Reference implementations

Both transforms are standard in compositional data analysis. In `R`, the
[`compositions`](https://cran.r-project.org/web/packages/compositions/refman/compositions.html#alr)
package exposes `alr()` / `alrInv()`. In `Python`,
[`scikit-bio`](https://scikit.bio/docs/dev/generated/skbio.stats.composition.alr_inv.html)
provides `skbio.stats.composition.alr()` and `alr_inv()`. DeCovarT’s
[`additive_logistic()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
and
[`additive_log_ratio()`](https://bastienchassagnol.github.io/DeCovarT/reference/additive_logistic.md)
are the same maps with the last part fixed as reference.

``` r

p <- c(0.2, 0.3, 0.5)
rho <- additive_log_ratio(p)

if (requireNamespace("compositions", quietly = TRUE)) {
  rho_pkg <- as.numeric(compositions::alr(compositions::acomp(p)))
  cat("alr agreement:", isTRUE(all.equal(as.numeric(rho), rho_pkg)), "\n")
  p_pkg <- as.numeric(compositions::alrInv(rho_pkg))
  cat("alrInv round-trip:", isTRUE(all.equal(p, p_pkg)), "\n")
}
```

Listing 11: Cross-check against the compositions package

## Checking the closed forms

The analytic derivatives are validated against Richardson extrapolation
in `R` ([Listing 12](#lst-check-r)).

``` r

rho <- c(0.4, -1.1)
p <- c(0.2, 0.3, 0.5)

# additive logistic: Jacobian and per-output Hessian
jac_ok <- all.equal(
  numDeriv::jacobian(additive_logistic, rho),
  additive_logistic_jacobian(rho)
)
hess_ok <- all.equal(
  numDeriv::hessian(function(r) additive_logistic(r)[1], rho),
  additive_logistic_hessian(rho)[1, , ]
)

# additive log-ratio: ambient Jacobian
alr_jac_ok <- all.equal(
  numDeriv::jacobian(additive_log_ratio, p),
  additive_log_ratio_jacobian(p)
)

c(logistic_jacobian = isTRUE(jac_ok),
  logistic_hessian = isTRUE(hess_ok),
  alr_jacobian = isTRUE(alr_jac_ok))
```

Listing 12: Validation with numDeriv
