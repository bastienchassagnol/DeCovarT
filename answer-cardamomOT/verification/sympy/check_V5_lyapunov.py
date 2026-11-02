# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V5: stationary covariance and precision of a linearised (Ornstein-Uhlenbeck)
gene network, from the Lyapunov equation A S + S A^T + Q = 0.

Derivation route: the Lyapunov equation is solved as a linear system in the
free entries of S (nothing is typed from the note), S is inverted, and the
result is compared with the closed forms printed in the note (eq:lyap2 and
Appendix A.1). The static structural equation model is derived from
x = B x + e, so its covariance is (I - B)^{-1} D (I - B)^{-T}.
Decays a_g and noise intensities q_g are positive; edge strengths b, c are
real of either sign.
"""

from __future__ import annotations

import json

import sympy as sp

a1, a2, a3, q1, q2, q3 = sp.symbols("a1 a2 a3 q1 q2 q3", positive=True)
bb, cc = sp.symbols("b c", real=True)


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


def solve_lyapunov(A: sp.Matrix, Q: sp.Matrix) -> sp.Matrix:
    """Solve A S + S A^T + Q = 0 for symmetric S as a linear system."""
    n = A.shape[0]
    S = sp.zeros(n, n)
    syms = []
    for i in range(n):
        for j in range(i, n):
            s = sp.Symbol(f"s{i}{j}")
            syms.append(s)
            S[i, j] = s
            S[j, i] = s
    M = A * S + S * A.T + Q
    sol = sp.linsolve([M[i, j] for i in range(n) for j in range(i, n)], syms)
    assert len(sol) == 1
    (vals,) = sol
    return S.subs(dict(zip(syms, vals))).applyfunc(sp.cancel)


def _zero(expr) -> bool:
    return sp.cancel(sp.together(expr)) == 0


def _pos_coeffs(expr, gens) -> bool:
    """Sufficient condition for positivity on the positive orthant (with b, c
    entering only through even powers): every expanded coefficient is > 0."""
    p = sp.Poly(sp.expand(expr), *gens)
    ok_coeffs = all(co > 0 for co in p.coeffs())
    # b and c must appear with even exponents only (they may be negative)
    idx = [gens.index(v) for v in (bb, cc) if v in gens]
    even = all(all(m[k] % 2 == 0 for k in idx) for m in p.monoms())
    return ok_coeffs and even


def run() -> list[dict]:
    out: list[dict] = []

    # ---------------------------------------------------------- two genes
    A = sp.Matrix([[-a1, 0], [bb, -a2]])
    S = solve_lyapunov(A, sp.diag(q1, q2))
    Om = S.inv(method="LU").applyfunc(sp.factor)

    # forms printed in the note (eq:lyap2 and Appendix A.1)
    note_cov = {
        (0, 0): q1 / (2 * a1),
        (0, 1): bb * q1 / (2 * a1 * (a1 + a2)),
        (1, 1): (a1 * (a1 + a2) * q2 + bb**2 * q1) / (2 * a1 * a2 * (a1 + a2)),
    }
    den = (a1 + a2) ** 2 * q2 + bb**2 * q1
    note_prec = {
        (0, 0): 2 * (a1 * (a1 + a2) ** 2 * q2 + (a1 + a2) * bb**2 * q1) / (q1 * den),
        (0, 1): -2 * a2 * (a1 + a2) * bb / den,
        (1, 1): 2 * a2 * (a1 + a2) ** 2 / den,
    }
    ok_cov = all(_zero(S[k] - f) for k, f in note_cov.items())
    out.append(
        _rec(
            "V5a",
            r"Two-gene Lyapunov solution: $\sigma_{11}$, $\sigma_{12}$ (eq:lyap2) and $\sigma_{22}$ (A.1)",
            "exact",
            "pass" if ok_cov else "fail",
            r"$\bm A\bm S+\bm S\bm A^{\mathsf T}+\bm Q=\bm 0$ solved as a linear system: $\sigma_{12}=" + sp.latex(sp.factor(S[0, 1])) + r"$, $\sigma_{22}=" + sp.latex(sp.factor(S[1, 1])) + "$" + ("; equal to the printed forms" if ok_cov else "; MISMATCH"),
            361,
            detail_full=(
                "A S + S A^T + Q = 0 solved as a linear system in (s11, s12, s22); engine: "
                f"s12 = {sp.sstr(sp.factor(S[0, 1]))}, s22 = {sp.sstr(sp.factor(S[1, 1]))}; "
                + ("differences with the note cancel to 0" if ok_cov else "MISMATCH")
            ),
        )
    )

    ok_prec = all(_zero(Om[k] - f) for k, f in note_prec.items())
    out.append(
        _rec(
            "V5b",
            r"Two-gene precision: $\omega_{12}$ (eq:lyap2), $\omega_{11}$ and $\omega_{22}$ (A.1)",
            "exact",
            "pass" if ok_prec else "fail",
            r"$\bm S^{-1}$ from the engine's $\bm S$: $\omega_{12}=" + sp.latex(sp.factor(Om[0, 1])) + "$" + ("; all three equal the printed forms" if ok_prec else "; MISMATCH"),
            361,
            detail_full=(
                f"S^-1 from the engine's S; omega12 = {sp.sstr(Om[0, 1])}; "
                + (
                    "omega11, omega12, omega22 equal the note's forms"
                    if ok_prec
                    else f"MISMATCH: {sp.sstr(Om)}"
                )
            ),
        )
    )

    num12 = sp.numer(sp.together(Om[0, 1]))
    roots_b = sp.solve(num12, bb)
    ok_iff = roots_b == [0] and _pos_coeffs(sp.denom(sp.together(Om[0, 1])), (a1, a2, q1, q2, bb))
    out.append(
        _rec(
            "V5c",
            r"$\omega_{12}\neq0$ if and only if $b\neq0$",
            "exact",
            "pass" if ok_iff else "fail",
            r"numerator of $\omega_{12}$: $" + sp.latex(sp.factor(num12)) + r"$, only root in $b$: " + ", ".join("$" + sp.latex(r) + "$" for r in roots_b) + r"; denominator $(a_1+a_2)^2q_2+b^2q_1>0$",
            366,
            detail_full=(
                f"numerator of omega12 = {sp.sstr(sp.factor(num12))}, only root in b: {roots_b}; "
                "denominator (a1+a2)^2 q2 + b^2 q1 > 0"
            ),
        )
    )

    # Diagonal versus self-loop. a_g = d_g - varsigma_g, so self-activation
    # lowers omega_gg iff d omega_gg / d a_g > 0.
    gens2 = (a1, a2, q1, q2, bb)
    d11 = sp.factor(sp.diff(Om[0, 0], a1))
    d22 = sp.factor(sp.diff(Om[1, 1], a2))
    pos = all(
        _pos_coeffs(part, gens2)
        for d in (d11, d22)
        for part in sp.fraction(sp.together(d))
    )
    ds22 = sp.factor(-sp.diff(S[1, 1], a2))
    sig_ok = _zero(sp.diff(S[0, 0], a1) + q1 / (2 * a1**2)) and all(
        _pos_coeffs(part, gens2) for part in sp.fraction(sp.together(ds22))
    )
    d22_a1 = sp.factor(sp.diff(Om[1, 1], a1))
    mixes = all(sp.diff(Om[0, 0], v) != 0 for v in (bb, a2, q2)) and d22_a1 != 0
    out.append(
        _rec(
            "V5d",
            r"$\omega_{gg}$ falls as self-activation rises ($a_g$ falls), and depends on the edge and on the other gene",
            "exact",
            "pass" if (pos and sig_ok and mixes) else "fail",
            r"$\partial\omega_{11}/\partial a_1$ and $\partial\omega_{22}/\partial a_2$ have numerators and denominators with positive coefficients; $\partial\omega_{22}/\partial a_1=" + sp.latex(d22_a1) + "$",
            391,
            detail_full=(
                "d omega11/d a1 and d omega22/d a2 have numerator and denominator with positive "
                "coefficients (b in even powers), so both are > 0; sigma11 and sigma22 decrease "
                "in their own a_g; omega11 depends on b, a2, q2, and "
                f"d omega22/d a1 = {sp.sstr(d22_a1)}"
            ),
        )
    )

    # ---------------------------------------------------------- three-gene fork
    A3 = sp.Matrix([[-a1, 0, 0], [bb, -a2, 0], [cc, 0, -a3]])
    S3 = solve_lyapunov(A3, sp.diag(q1, q2, q3))
    # clear denominators: P = L S3 is polynomial, Omega = L adj(P) / det(P)
    L = 2 * a1 * a2 * a3 * (a1 + a2) * (a1 + a3) * (a2 + a3)
    P = (S3 * L).applyfunc(lambda e: sp.expand(sp.cancel(e)))
    assert all(sp.denom(sp.together(e)) == 1 for e in P)
    detP = sp.expand(P.det(method="bareiss"))
    adj12 = sp.expand(-P.minor_submatrix(2, 1).det())  # cofactor C_{3,2} = adj(P)_{2,3}
    w23 = sp.factor(L * adj12 / detP)
    note_num = -4 * a2 * a3 * bb * cc * q1 * (a1 + a2) * (a1 + a3) * (a2 + a3)
    D = sp.factor(sp.cancel(note_num / w23))
    Dn, Dd = sp.fraction(sp.together(D))
    gens3 = (a1, a2, a3, q1, q2, q3, bb, cc)
    D_poly = Dd == 1
    # D has one negative monomial (-2 a2 a3 b^2 c^2 q1^2); it combines into
    # b^2 c^2 q1^2 (a2 - a3)^2, and the remainder has positive coefficients.
    sq = bb**2 * cc**2 * q1**2 * (a2 - a3) ** 2
    D_pos = D_poly and _pos_coeffs(sp.expand(Dn) - sp.expand(sq), gens3)
    out.append(
        _rec(
            "V5e",
            r"Fork $1\to2$ ($b$), $1\to3$ ($c$): $\omega_{23}=-4a_2a_3bc\,q_1(a_1+a_2)(a_1+a_3)(a_2+a_3)/\mathcal D$ with $\mathcal D>0$",
            "exact",
            "pass" if (D_poly and D_pos) else "fail",
            r"six-unknown Lyapunov system solved exactly; $\mathcal D=b^2c^2q_1^2(a_2-a_3)^2$ plus a polynomial with positive coefficients, so $\mathcal D>0$",
            370,
            detail_full=(
                "6-unknown Lyapunov system solved exactly; omega23 = L adj(P)_23 / det P with "
                "P = L S polynomial; note numerator / omega23 is a polynomial D; "
                "D = b^2 c^2 q1^2 (a2-a3)^2 + (polynomial with positive coefficients), hence D > 0"
            ),
        )
    )

    # ---------------------------------------------------------- static SEM fork
    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    B = sp.Matrix([[0, 0, 0], [bb, 0, 0], [cc, 0, 0]])
    IB = sp.eye(3) - B
    S_sem = IB.inv() * sp.diag(d1, d2, d3) * IB.inv().T
    Om_sem = S_sem.inv(method="LU").applyfunc(sp.cancel)
    ok_sem = Om_sem[1, 2] == 0 and _zero(Om_sem[0, 1] + bb / d2)
    out.append(
        _rec(
            "V5f",
            r"Static structural equation model with the same fork: $\omega_{23}=0$",
            "exact",
            "pass" if ok_sem else "fail",
            r"$(\bm I-\bm B)^{-1}\bm D(\bm I-\bm B)^{-\mathsf T}$ inverted exactly: $\omega_{23}=" + sp.latex(Om_sem[1, 2]) + r"$, $\omega_{12}=" + sp.latex(Om_sem[0, 1]) + "$",
            370,
            detail_full=(
                "x = Bx + e, e ~ (0, diag(d1,d2,d3)); covariance (I-B)^-1 D (I-B)^-T inverted "
                f"exactly: omega23 = {sp.sstr(Om_sem[1, 2])}, omega12 = {sp.sstr(Om_sem[0, 1])}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
