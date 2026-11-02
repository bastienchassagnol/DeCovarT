# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V10: rank of the empirical covariance and existence of the unconstrained
Gaussian maximum-likelihood estimate (identifiability box, part (a)).

The note says: with the mean estimated, V = X_c' X_c / C has rank equal to
rank(X_c) <= min(C - 1, G), with equality for data in general position; V is
positive definite only if C >= G + 1; an exact linear (or affine) relation
between genes lowers the rank further; when V is singular the Gaussian
log-likelihood is unbounded above, so the MLE does not exist.

Derivation routes:
- V10a: X_c = H X with H = I - 11'/C of rank C - 1, so rank V <= C - 1 is an
  algebraic identity. G = 3: with symbolic data, det V is identically 0 at
  C = 3 and is a non-zero polynomial at C = 4 (evaluated at a rational point);
  the exact rank at rational points is min(C - 1, G) for C = 2, ..., 6. A
  polynomial that is not identically zero vanishes on a Lebesgue-null set, so
  the rank is min(C - 1, G) almost surely for data with a density.
- V10b: exact affine collinearity x3 = 2 x1 - x2 + 5, and a constant gene,
  at C = 4 and C = 10 (rational data): the rank drops below min(C - 1, G), and
  the relation is a null vector of V.
- V10c: profile log-likelihood l(Sigma) = -(C/2){log det Sigma + tr(V Sigma^-1)}
  along Sigma_eps = V + eps u u' (u a unit null vector of a singular V):
  tr(V Sigma_eps^-1) stays equal to rank V, log det Sigma_eps = log eps + const,
  and l -> +oo as eps -> 0+.
- V10d (flag): Zwiernik, Uhler & Richards (2017, p. 1270 and Proposition 5,
  p. 1278) write that S_n, built with the sample mean, is positive definite
  with probability 1 "when n >= p". V10a gives rank p - 1 at n = p, so the
  condition is n >= p + 1.
"""

from __future__ import annotations

import json
import random

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


def _centred_cov(X: sp.Matrix) -> sp.Matrix:
    C = X.rows
    H = sp.eye(C) - sp.ones(C, C) / C
    Xc = H * X
    return (Xc.T * Xc) / C


def _rational_data(rng: random.Random, C: int, G: int) -> sp.Matrix:
    return sp.Matrix(C, G, lambda i, j: sp.Rational(rng.randint(-40, 40), rng.randint(1, 9)))


def run() -> list[dict]:
    out: list[dict] = []
    rng = random.Random(20261030)
    G = 3

    # ---------------- V10a: generic rank min(C - 1, G) ----------------
    X3 = sp.Matrix(3, G, sp.symbols("x1:10", real=True))
    det_c3 = sp.expand(_centred_cov(X3).det())
    zero_c3 = det_c3 == 0
    X4 = sp.Matrix(4, G, sp.symbols("y1:13", real=True))
    V4 = _centred_cov(X4)
    pt = {s: sp.Rational(rng.randint(-30, 30), rng.randint(1, 7)) for s in X4.free_symbols}
    det_c4_pt = sp.nsimplify(V4.subs(pt).det())
    nonzero_c4 = det_c4_pt != 0
    rank_H_ok = (sp.eye(5) - sp.ones(5, 5) / 5).rank() == 4
    ranks = {}
    ok_ranks = True
    for C in range(2, 7):
        X = _rational_data(rng, C, G)
        H = sp.eye(C) - sp.ones(C, C) / C
        r_v, r_xc = _centred_cov(X).rank(), (H * X).rank()
        ranks[C] = r_v
        ok_ranks &= (r_v == r_xc == min(C - 1, G))
    out.append(
        _rec(
            "V10a",
            r"With the mean estimated, $\operatorname{rank}\bm V=\operatorname{rank}\bm X_c\le\min(C-1,G)$, "
            r"with equality for data in general position; $\bm V\succ0$ needs $C\ge G+1$",
            "exact",
            "pass" if zero_c3 and nonzero_c4 and ok_ranks and rank_H_ok else "fail",
            r"$G=3$, symbolic data: $\det\bm V\equiv0$ at $C=3$; at $C=4$, $\det\bm V$ is a non-zero polynomial "
            r"(value $" + sp.latex(det_c4_pt) + r"$ at a rational point); exact ranks at rational points for "
            r"$C=2,\dots,6$: " + ", ".join(str(ranks[c]) for c in sorted(ranks)),
            0,
            detail_full=(
                f"det V (G=3, C=3, symbolic) == 0: {zero_c3}; det V (G=3, C=4) at a rational point = "
                f"{det_c4_pt}; rank(I - 11'/5) = 4: {rank_H_ok}; ranks of V (= ranks of X_c) for C=2..6: {ranks}"
            ),
        )
    )

    # ---------------- V10b: exact collinearity lowers the rank further ----------------
    drops = {}
    ok_drop = True
    for C in (4, 10):
        base = _rational_data(rng, C, 2)
        x3 = 2 * base[:, 0] - base[:, 1] + sp.Matrix([5] * C)  # affine relation
        V_aff = _centred_cov(base.row_join(x3))
        u = sp.Matrix([2, -1, -1])  # x3 - 2 x1 + x2 = 5, constant removed by centring
        null_ok = (V_aff * u).applyfunc(sp.simplify) == sp.zeros(3, 1)
        const = _centred_cov(base.row_join(sp.Matrix([7] * C)))  # a gene constant within the type
        drops[C] = (V_aff.rank(), const.rank())
        ok_drop &= null_ok and V_aff.rank() == 2 and const.rank() == 2 and min(C - 1, G) == 3
    out.append(
        _rec(
            "V10b",
            r"An exact linear or affine relation between genes, or a constant gene, lowers the rank below "
            r"$\min(C-1,G)$",
            "exact",
            "pass" if ok_drop else "fail",
            r"$G=3$, $x_3=2x_1-x_2+5$: $\bm V(2,-1,-1)^{\mathsf T}=\bm 0$; rank (affine, constant gene) $="
            + "(%d,%d)" % drops[4]
            + r"$ at $C=4$ and $"
            + "(%d,%d)" % drops[10]
            + r"$ at $C=10$, against $\min(C-1,G)=3$",
            0,
            detail_full=f"ranks (affine relation, constant gene) by C: {drops}; null vector checked: {ok_drop}",
        )
    )

    # ---------------- V10c: unbounded likelihood when V is singular ----------------
    X = _rational_data(rng, 3, G)  # C = 3 <= G: V singular
    V = _centred_cov(X)
    ns = V.nullspace()
    u = ns[0] / sp.sqrt((ns[0].T * ns[0])[0])
    eps = sp.Symbol("epsilon", positive=True)
    Ceps = sp.Symbol("C", positive=True)
    Sig = V + eps * u * u.T
    tr_term = sp.simplify((V * Sig.inv()).trace())
    logdet = sp.log(sp.factor(Sig.det()))
    loglik = -Ceps / 2 * (logdet + tr_term)
    lim = sp.limit(loglik.subs(Ceps, 3), eps, 0, "+")
    pd_eps = all(sp.simplify(m.subs(eps, sp.Rational(1, 10**6))) > 0 for m in
                 [Sig[:k, :k].det() for k in range(1, G + 1)])
    ok_unb = (tr_term == V.rank()) and lim == sp.oo and pd_eps and len(ns) == 1
    out.append(
        _rec(
            "V10c",
            r"If $\bm V$ is singular, the Gaussian log-likelihood is unbounded above: the MLE does not exist",
            "exact",
            "pass" if ok_unb else "fail",
            r"$C=G=3$; $\bm\Sigma_\varepsilon=\bm V+\varepsilon\bm u\bm u^{\mathsf T}$, $\bm V\bm u=\bm 0$: "
            r"$\operatorname{tr}(\bm V\bm\Sigma_\varepsilon^{-1})=" + sp.latex(tr_term)
            + r"$ for every $\varepsilon$, $\det\bm\Sigma_\varepsilon\propto\varepsilon$, and "
            r"$\ell(\bm\Sigma_\varepsilon)\to" + sp.latex(lim) + r"$ as $\varepsilon\to0^+$",
            0,
            detail_full=(
                f"V = {sp.sstr(V)}; nullity {len(ns)}; tr(V Sigma_eps^-1) = {tr_term}; det Sigma_eps = "
                f"{sp.sstr(sp.factor(Sig.det()))}; limit of l = {lim}; Sigma_eps PD at eps=1e-6: {pd_eps}"
            ),
        )
    )

    # ---------------- V10d: flag on Zwiernik, Uhler & Richards (2017) ----------------
    flag = zero_c3 and ranks[3] == 2  # n = p = 3 gives rank p - 1
    out.append(
        _rec(
            "V10d",
            r"Zwiernik et al. (2017, pp. 1270, 1278): with the sample mean, $S_n\succ0$ almost surely "
            r"``when $n\geqslant p$''; the condition is $n\geqslant p+1$",
            "exact",
            "flag" if flag else "fail",
            r"$n=p=3$: $\det S_n\equiv0$ and rank $" + str(ranks[3]) + r"=p-1$ (V10a); their own "
            r"$nS_n\sim W_p(n-1,\bm\Sigma)$ (p. 1274) is singular at $n=p$",
            0,
            detail_full=f"det V identically 0 at n = p = 3: {zero_c3}; rank at a rational point: {ranks[3]}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
