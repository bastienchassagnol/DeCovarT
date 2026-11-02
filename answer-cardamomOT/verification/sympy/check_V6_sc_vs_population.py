# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V6: bulk covariance under single-cell references (linear in p) against
population references (quadratic in p), at J = 2 and G = 2.

Derivation route: covariances are read off the Hessian of the cumulant
generating function (log of the characteristic function, with t -> -i t),
built from the independence structure of each model:
- population: y = sum_j p_j x_j, x_j independent N(mu_j, Sigma_j), so
  K_y(t) = sum_j K_j(p_j t);
- single-cell, composition given: y = a * (sum of m_j cells of type j), so
  K_y(t) = sum_j m_j K_j(a t), with m_j = n p_j;
- single-cell, cells drawn i.i.d. from the mixture: K_y(t) = n log sum_j p_j
  exp K_j(a t) (the phi_mix(t)^n form of the note).
Linearity and quadratic dependence in p are then tested exactly.

V6e reaches the i.i.d.-draw model by a second, independent route: the type
counts m ~ Multinomial(n, p) are averaged out of the conditional characteristic
function prod_j phi_j^{m_j} (exact for n = 1..8 at J = 2 and n = 1..5 at
J = 3), and the covariance follows from the law of total covariance with the
multinomial moments built from one-hot type indicators (symbolic n).
"""

from __future__ import annotations

import json

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line, detail_full=None):
    """detail: short text for the ledger table, with $...$ LaTeX (sympy.latex);
    detail_full: the full engine output, kept in ledger.json."""
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def _hess0(K, ts):
    return sp.hessian(K, ts).subs({t: 0 for t in ts}).applyfunc(sp.simplify)


def _zero(M):
    return all(sp.simplify(e) == 0 for e in M)


def run() -> list[dict]:
    out: list[dict] = []
    p = sp.Symbol("p", positive=True)
    a, n = sp.symbols("a n", positive=True)
    ts = sp.symbols("t1 t2", real=True)
    tv = sp.Matrix(ts)
    mus, Sigs = [], []
    for j in (1, 2):
        mus.append(sp.Matrix(sp.symbols(f"m{j}_1 m{j}_2", real=True)))
        s11, s22 = sp.symbols(f"S{j}_11 S{j}_22", positive=True)
        s12 = sp.Symbol(f"S{j}_12", real=True)
        Sigs.append(sp.Matrix([[s11, s12], [s12, s22]]))

    def K(j, tt):  # Gaussian CGF of one profile / one cell of type j
        return (tt.T * mus[j])[0] + (tt.T * Sigs[j] * tt)[0] / 2

    w = [p, 1 - p]

    # population references
    K_pop = sum(K(j, w[j] * tv) for j in range(2))
    C_pop = _hess0(K_pop, ts)
    ok_pop = _zero(C_pop - (p**2 * Sigs[0] + (1 - p) ** 2 * Sigs[1]))
    curv = C_pop.diff(p, 2).applyfunc(sp.simplify)
    out.append(
        _rec(
            "V6a",
            r"Population references: $\operatorname{Cov}\bm y=\sum_j p_j^2\bm\Sigma_j$, quadratic in $\bm p$",
            "exact",
            "pass" if ok_pop and _zero(curv - 2 * (Sigs[0] + Sigs[1])) else "fail",
            r"Hessian at $\bm 0$ of $K_{\bm y}(\bm t)=\sum_j K_j(p_j\bm t)$ is $p^2\bm\Sigma_1+(1-p)^2\bm\Sigma_2$; "
            r"$\partial_p^2=2(\bm\Sigma_1+\bm\Sigma_2)\succ0$",
            75,
            detail_full=(
                "Hessian at 0 of K_y(t) = sum_j K_j(p_j t) (independent profiles) equals "
                "p^2 Sigma_1 + (1-p)^2 Sigma_2; d2/dp2 = 2 (Sigma_1 + Sigma_2), positive definite"
            ),
        )
    )

    # single-cell references, composition given (m_j = n p_j)
    K_sc = sum(n * w[j] * K(j, a * tv) for j in range(2))
    C_sc = _hess0(K_sc, ts)
    phi = a**2 * n
    ok_sc = _zero(C_sc - phi * (p * Sigs[0] + (1 - p) * Sigs[1]))
    lin = _zero(C_sc.diff(p, 2))
    out.append(
        _rec(
            "V6b",
            r"Single-cell references, given $m_j=n_ip_j$ cells per type: $\operatorname{Cov}\bm y="
            r"\phi_i\sum_j p_j\bm\Sigma_{W,j}$, linear in $\bm p$, $\phi_i=a_i^2n_i$",
            "exact",
            "pass" if ok_sc and lin else "fail",
            r"Hessian at $\bm 0$ of $\sum_j m_jK_j(a\bm t)$ is $a^2n\{p\bm\Sigma_1+(1-p)\bm\Sigma_2\}$; "
            r"$\partial_p^2=\bm 0$; conditional on the counts $m_j$",
            489,
            detail_full=(
                "Hessian at 0 of K_y(t) = sum_j m_j K_j(a t) equals a^2 n {p Sigma_1 + (1-p) "
                "Sigma_2}; d2/dp2 = 0. Holds conditional on the composition m_j = n p_j"
            ),
        )
    )

    # the two cannot agree for every p
    k = sp.Symbol("k", positive=True)
    gap = (C_pop - k * (p * Sigs[0] + (1 - p) * Sigs[1])).applyfunc(sp.expand)
    lead = [sp.Poly(e, p).coeff_monomial(p**2) for e in gap]
    ok_diff = all(sp.simplify(lc - (S1 + S2)) == 0 for lc, S1, S2 in zip(lead, Sigs[0], Sigs[1]))
    out.append(
        _rec(
            "V6c",
            r"At $J=2$ the quadratic and the linear covariance differ for any common scale $k$",
            "exact",
            "pass" if ok_diff else "fail",
            r"the $p^2$ coefficient of $\sum_jp_j^2\bm\Sigma_j-k\sum_jp_j\bm\Sigma_j$ is "
            r"$\bm\Sigma_1+\bm\Sigma_2$, with a positive diagonal; equality only at isolated $p$",
            75,
            detail_full=(
                "p^2 coefficient of sum p_j^2 Sigma_j - k sum p_j Sigma_j is Sigma_1 + Sigma_2, "
                "whose diagonal is positive, so no k makes them equal for all p; they meet only "
                "at isolated p"
            ),
        )
    )

    # i.i.d. draws from the mixture: phi_cell(t)^n
    K_mix = n * sp.log(sum(w[j] * sp.exp(K(j, a * tv)) for j in range(2)))
    C_mix = _hess0(K_mix, ts)
    Dlt = mus[0] - mus[1]
    extra = (C_mix - C_sc).applyfunc(sp.factor)
    ok_mix = _zero(extra - a**2 * n * p * (1 - p) * Dlt * Dlt.T)
    out.append(
        _rec(
            "V6d",
            r"Cells drawn i.i.d. from the type mixture: $\operatorname{Cov}\bm y$ gains "
            r"$a^2np(1-p)\bm\Delta\bm\Delta^{\mathsf T}$, quadratic in $p$",
            "exact",
            "pass" if ok_mix else "fail",
            r"Hessian at $\bm 0$ of $n\log\sum_jp_je^{K_j(a\bm t)}$ minus the fixed-count covariance "
            r"equals $a^2np(1-p)(\bm\mu_1-\bm\mu_2)(\bm\mu_1-\bm\mu_2)^{\mathsf T}$",
            408,
            detail_full=(
                "Hessian at 0 of n log sum_j p_j exp K_j(a t) minus the fixed-composition "
                "covariance = a^2 n p (1-p) (mu_1 - mu_2)(mu_1 - mu_2)^T; linearity in p is "
                "verified only given the cell counts m_j"
            ),
        )
    )

    # V6e: the i.i.d.-draw model by a second route.
    # (i) characteristic function: given the counts m, the bulk c.f. is
    # prod_j phi_j^{m_j}; averaging over m ~ Multinomial(n, p) must give
    # phi_mix^n with phi_mix = sum_j p_j phi_j (one cell's c.f., by the law of
    # total expectation over its type).
    f = sp.symbols("f1:4")  # values of phi_j(t) at one t, left symbolic
    ok_cf = True
    for J_, nmax in ((2, 8), (3, 5)):
        pj = list(sp.symbols(f"q1:{J_}", positive=True))
        pj.append(1 - sum(pj))
        phi_mix = sum(pj[j] * f[j] for j in range(J_))
        for n_ in range(1, nmax + 1):
            avg = 0
            for m in _compositions(n_, J_):
                coef = sp.factorial(n_) / sp.prod([sp.factorial(mj) for mj in m])
                avg += coef * sp.prod([pj[j] ** m[j] * f[j] ** m[j] for j in range(J_)])
            ok_cf &= sp.expand(avg - phi_mix**n_) == 0
    # (ii) covariance by the law of total covariance, symbolic n. One cell's
    # type is a one-hot vector B with P(B = e_j) = p_j; m is the sum of n
    # independent copies, so E m = n E B and Cov m = n Cov B.
    e = [sp.Matrix([1, 0]), sp.Matrix([0, 1])]
    EB = sum((w[j] * e[j] for j in range(2)), sp.zeros(2, 1))
    EBB = sum((w[j] * e[j] * e[j].T for j in range(2)), sp.zeros(2, 2))
    Em, Cm = n * EB, n * (EBB - EB * EB.T)
    Mmat = sp.Matrix.hstack(mus[0], mus[1])  # E(y | m) = a M m
    # Cov(y | m) = a^2 sum_j m_j Sigma_j (V6b), averaged over m
    E_cond_cov = a**2 * (Em[0] * Sigs[0] + Em[1] * Sigs[1])
    cov_tot = (E_cond_cov + a**2 * Mmat * Cm * Mmat.T).applyfunc(sp.expand)
    printed = a**2 * n * (p * Sigs[0] + (1 - p) * Sigs[1]) + a**2 * n * p * (1 - p) * Dlt * Dlt.T
    ok_tot = _zero(cov_tot - printed)
    curv_mix = (cov_tot.diff(p, 2) + 2 * a**2 * n * Dlt * Dlt.T).applyfunc(sp.expand)
    ok_quad = _zero(curv_mix)
    out.append(
        _rec(
            "V6e",
            r"With i.i.d. draws, the bulk characteristic function is $\varphi_{\mathrm{mix}}(\bm t)^{n}$, "
            r"$\varphi_{\mathrm{mix}}=\sum_jp_j\varphi_j$, and the covariance is "
            r"$a^2n\sum_jp_j\bm\Sigma_{W,j}+a^2np(1-p)\bm\Delta\bm\Delta^{\mathsf T}$ at $J=2$",
            "exact",
            "pass" if ok_cf and ok_tot and ok_quad else "fail",
            r"averaging $\prod_j\varphi_j^{m_j}$ over $\bm m\sim\mathrm{Mult}(n,\bm p)$ gives "
            r"$\varphi_{\mathrm{mix}}^{\,n}$ ($n\le8$ at $J=2$, $n\le5$ at $J=3$); law of total covariance "
            r"with $\operatorname{Cov}\bm m=n\{\operatorname{diag}(\bm p)-\bm p\bm p^{\mathsf T}\}$ returns the printed "
            r"term; $\partial_p^2=-2a^2n\bm\Delta\bm\Delta^{\mathsf T}\neq\bm 0$",
            476,
            detail_full=(
                f"c.f. average over the multinomial equals phi_mix^n: {ok_cf} (J = 2, n = 1..8; "
                f"J = 3, n = 1..5). E m = {list(Em)}, Cov m = {sp.sstr(Cm.applyfunc(sp.factor))}. "
                f"Total covariance minus printed form = {sp.sstr((cov_tot - printed).applyfunc(sp.simplify))}; "
                f"d2/dp2 + 2 a^2 n Delta Delta^T = {sp.sstr(curv_mix)}"
            ),
        )
    )
    return out


def _compositions(n_: int, J_: int):
    """All count vectors (m_1, ..., m_J) of non-negative integers summing to n."""
    if J_ == 1:
        yield (n_,)
        return
    for first in range(n_ + 1):
        for rest in _compositions(n_ - first, J_ - 1):
            yield (first, *rest)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
