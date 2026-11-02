# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V1: the deconvolution model y = mu p + sum_j p_j eps_j + e as a linear mixed
model; its marginal covariance, and the BLUP of eps_j against the Wiener
filter / posterior mean.

Derivation routes (none re-types the note's result):
- covariance: every random term is written as a linear map of one standard
  normal vector w (eps_j = L_j w_j with L_j lower triangular, e = D^{1/2} w_e),
  so Var(y) = M M^T for the stacked map M;
- BLUP / MAP: Henderson's mixed-model equations
  (Z^T D^-1 Z + G^-1) u = Z^T D^-1 r, i.e. the stationarity condition of the
  log posterior of u; the note's predictor must satisfy them identically
  (r is parametrised as V s to avoid a symbolic 4 x 4 inverse);
- posterior mean: Gaussian conditioning E(u | y) = Cov(u, y) Var(y)^-1 r with
  Cov(u, y) computed from the same linear map.
The three results are compared with the note's sum_j p_j^2 Sigma_j + D and
p_j Sigma_j V^-1 (y - mu p) (Appendix A.3, which is also the Wiener filter
printed in the DeCovarT article, sec:wiener).
G = 2 genes and J = 2 types with fully symbolic 2 x 2 covariances.
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


def _is_zero_matrix(M: sp.Matrix) -> bool:
    return all(sp.cancel(sp.together(e)) == 0 for e in M)


def run() -> list[dict]:
    out: list[dict] = []
    G = 2
    p1, p2 = sp.symbols("p1 p2", positive=True)
    p = [p1, p2]
    # Cholesky factors of Sigma_j (positive diagonal)
    Ls = []
    for j in (1, 2):
        l11, l22 = sp.symbols(f"l{j}_11 l{j}_22", positive=True)
        l21 = sp.Symbol(f"l{j}_21", real=True)
        Ls.append(sp.Matrix([[l11, 0], [l21, l22]]))
    Sig = [L * L.T for L in Ls]
    dvals = sp.symbols("d1 d2", positive=True)
    Dm = sp.diag(*dvals)
    Dhalf = sp.diag(*[sp.sqrt(d) for d in dvals])

    # y - mu p = sum_j p_j L_j w_j + D^{1/2} w_e = M w,  w ~ N(0, I)
    M = sp.Matrix.hstack(p1 * Ls[0], p2 * Ls[1], Dhalf)
    V_engine = (M * M.T).applyfunc(sp.expand)
    V_note = (p1**2 * Sig[0] + p2**2 * Sig[1] + Dm).applyfunc(sp.expand)
    ok_v = _is_zero_matrix(V_engine - V_note)
    out.append(
        _rec(
            "V1a",
            r"Marginal covariance of the mixed model: $\operatorname{Var}\bm y=\sum_j p_j^2\bm\Sigma_j+\bm D$",
            "exact",
            "pass" if ok_v else "fail",
            r"$\bm y-\bm\mu\bm p=\bm M\bm w$ with $\bm w$ standard normal; $\bm M\bm M^{\mathsf T}$ "
            r"equals the printed form entrywise ($G=J=2$, symbolic Cholesky factors)",
            417,
            detail_full=(
                "y - mu p written as M w with w standard normal (eps_j = L_j w_j, "
                "e = D^1/2 w_e); M M^T equals the note's form entrywise (G = J = 2, "
                "symbolic Cholesky factors)"
            ),
        )
    )

    # random-effect design and covariance of u = (eps_1, eps_2)
    Z = sp.Matrix.hstack(p1 * sp.eye(G), p2 * sp.eye(G))
    Gm = sp.diag(Sig[0], Sig[1])
    # parametrise r = y - mu p as r = V s (V is invertible), which avoids a
    # symbolic 4 x 4 inverse: the note's predictor then reads p_j Sigma_j s
    s_vec = sp.Matrix(sp.symbols("s1 s2", real=True))
    r = V_note * s_vec
    u_note = sp.Matrix.vstack(p1 * Sig[0] * s_vec, p2 * Sig[1] * s_vec)

    # Henderson's mixed-model equations (= stationarity of the log posterior of
    # u); the left-hand matrix is positive definite, so the solution is unique
    Dinv = Dm.inv()
    Ginv = sp.diag(Sig[0].inv(), Sig[1].inv())
    lhs = Z.T * Dinv * Z + Ginv
    rhs = Z.T * Dinv * r
    ok_blup = _is_zero_matrix(lhs * u_note - rhs)
    out.append(
        _rec(
            "V1b",
            r"BLUP (= MAP of $\bm\varepsilon_j$ given $\bm p$) is $p_j\bm\Sigma_j\bm V^{-1}(\bm y-\bm\mu\bm p)$, the Wiener filter",
            "exact",
            "pass" if ok_blup else "fail",
            "the printed predictor solves Henderson's mixed-model equations identically; their "
            "matrix is positive definite, so it is the unique BLUP and MAP",
            496,
            detail_full=(
                "the note's predictor (with r = V s) satisfies Henderson's equations "
                "(Z'D^-1 Z + G^-1) u = Z'D^-1 r identically; the system matrix is positive "
                "definite, so it is the unique BLUP/MAP; same formula as the article's "
                "Wiener filter (sec:wiener)"
            ),
        )
    )

    # Gaussian conditioning: Cov(u, y) from the linear map, u = Mu w
    Mu = sp.Matrix.hstack(
        sp.Matrix.vstack(Ls[0], sp.zeros(2, 2)),
        sp.Matrix.vstack(sp.zeros(2, 2), Ls[1]),
        sp.zeros(4, 2),
    )
    C_uy = (Mu * M.T).applyfunc(sp.expand)
    post_mean = C_uy * s_vec  # Cov(u, y) V^-1 r with r = V s
    ok_post = _is_zero_matrix(post_mean - u_note)
    out.append(
        _rec(
            "V1c",
            r"Posterior mean $\mathbb E(\bm\varepsilon_j\mid\bm y,\bm p)$ equals the BLUP in the Gaussian model",
            "exact",
            "pass" if ok_post else "fail",
            r"$\operatorname{Cov}(\bm u,\bm y)$ from the linear map has blocks $p_j\bm\Sigma_j$; "
            r"Gaussian conditioning returns the predictor of V1b",
            496,
            detail_full=(
                "Cov(u, y) = Mu M^T derived from the linear map has blocks p_j Sigma_j; "
                "E(u | y) = Cov(u, y) V^-1 r therefore equals the predictor of V1b"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
