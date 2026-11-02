# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V17: what the presence gate of vignette 03-1 (#sec-fs-presence) excludes.

The empirical within-type covariance of type j is
V_j = (C - 1)^{-1} sum_c (x_c - xbar)(x_c - xbar)', built here from symbolic
counts (nothing typed from the vignette).

Derivation route:
- V17a: a gene with x_gc = 0 in every cell of the type. Row and column g of
  V_j are identically zero, det V_j = 0 (rank at most G - 1), and the sample
  correlation v_gh / sqrt(v_gg v_hh) has numerator and denominator both zero
  (0/0, undefined), for generic counts of the other genes (C = 4, G = 3).
- V17b: a gene detected in one cell c only, with count a > 0. Its centred
  vector is a (e_c - 1/C); its correlation with any gene h is
  h~_c / (sqrt(1 - 1/C) ||h~||), h~ the centred vector of h (generic C = 4
  data). Two genes detected only in the same cell have correlation exactly 1,
  and their 2 x 2 covariance is singular, whatever their counts.
- V17c: with delta_min C_j >= 2 a gene passing the gate has at least two
  non-zero cells, which excludes both cases (linear certificate).
"""

from __future__ import annotations

import json

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line, detail_full=None):
    return {
        "id": id_,
        "claim": claim,
        "kind": kind,
        "status": status,
        "detail": detail,
        "detail_full": detail if detail_full is None else detail_full,
        "tex_line": tex_line,
    }


def sample_cov(X):
    """X: G x C matrix of counts; returns the G x G sample covariance."""
    G, C = X.shape
    xbar = X * sp.ones(C, 1) / C
    Xc = X - xbar * sp.ones(1, C)
    return (Xc * Xc.T) / (C - 1), Xc


def run() -> list[dict]:
    out = []
    C = 4

    # ---- V17a: all-zero gene ------------------------------------------------
    h = sp.symbols("h1:5", positive=True)
    k = sp.symbols("k1:5", positive=True)
    X = sp.Matrix([[0] * C, list(h), list(k)])
    V, _ = sample_cov(X)
    row_zero = all(sp.simplify(V[0, j]) == 0 for j in range(3))
    det0 = sp.simplify(V.det()) == 0
    num, den2 = V[0, 1], V[0, 0] * V[1, 1]
    corr_undefined = sp.simplify(num) == 0 and sp.simplify(den2) == 0
    rank_le = V.rank(simplify=True) <= 2
    ok_a = row_zero and det0 and corr_undefined and rank_le
    out.append(
        _rec(
            "V17a",
            r"A gene with zero counts in every cell of type $j$ gives a zero row and column of $\bm V_j$, "
            r"$\det\bm V_j=0$, and an undefined ($0/0$) correlation with every other gene",
            "exact",
            "pass" if ok_a else "fail",
            r"Symbolic counts, $C=4$, $G=3$: $v_{1h}=0$ for all $h$, $\det\bm V_j=0$, rank $\le2$; "
            r"$v_{12}=0$ and $v_{11}v_{22}=0$",
            "",
            detail_full=f"row zero {row_zero}, det {det0}, corr 0/0 {corr_undefined}, rank<=2 {rank_le}",
        )
    )

    # ---- V17b: gene detected in one cell only ------------------------------
    a, b = sp.symbols("a b", positive=True)
    Xb = sp.Matrix([[a, 0, 0, 0], list(h)])
    Vb, Xc = sample_cov(Xb)
    r = Vb[0, 1] / sp.sqrt(Vb[0, 0] * Vb[1, 1])
    htil = Xc[1, :]
    target = htil[0] / (sp.sqrt(1 - sp.Rational(1, C)) * sp.sqrt((htil * htil.T)[0, 0]))
    ok_formula = sp.simplify(sp.cancel(sp.together((r - target)))) == 0
    if not ok_formula:  # compare squares, signs agree by construction
        ok_formula = sp.simplify(sp.cancel(sp.together(r**2 - target**2))) == 0 and \
            sp.sign(r.subs({h[0]: 5, h[1]: 1, h[2]: 2, h[3]: 3, a: 2})) == \
            sp.sign(target.subs({h[0]: 5, h[1]: 1, h[2]: 2, h[3]: 3}))
    Xs = sp.Matrix([[a, 0, 0, 0], [b, 0, 0, 0]])
    Vs, _ = sample_cov(Xs)
    r_same = sp.simplify(Vs[0, 1] / sp.sqrt(Vs[0, 0] * Vs[1, 1]))
    det_same = sp.simplify(Vs.det())
    ok_b = ok_formula and r_same == 1 and det_same == 0
    out.append(
        _rec(
            "V17b",
            r"A gene detected in one cell $c$ only has correlation $\tilde h_c/(\sqrt{1-1/C_j}\,\lVert\tilde{\bm h}\rVert)$ "
            r"with gene $h$; two genes detected only in the same cell have correlation $1$ and a singular covariance",
            "exact",
            "pass" if ok_b else "fail",
            rf"Symbolic counts, $C=4$: formula {'holds' if ok_formula else 'fails'}; "
            rf"counts $(a,0,0,0)$, $(b,0,0,0)$: $r={sp.latex(r_same)}$, $\det={sp.latex(det_same)}$",
            "",
        )
    )

    # ---- V17c: the gate threshold -------------------------------------------
    # delta_gj = (number of non-zero cells)/C_j >= delta_min and delta_min*C_j >= 2
    # imply at least 2 non-zero cells: check over a grid of integers.
    # Certificate: nz - 2 = Cj (nz/Cj - dmin) + (dmin Cj - 2), both brackets >= 0.
    nz, Cj, dmin = sp.symbols("n_z C_j delta_min", positive=True)
    cert = Cj * (nz / Cj - dmin) + (dmin * Cj - 2)
    ok_c = sp.simplify(sp.expand(cert - (nz - 2))) == 0
    out.append(
        _rec(
            "V17c",
            r"If $\delta_{\min}C_j\ge2$, a gene that passes $\delta_{gj}\ge\delta_{\min}$ is non-zero in at least two "
            r"cells of type $j$, so the cases of V17a and V17b are excluded",
            "exact",
            "pass" if ok_c else "fail",
            r"Linear certificate $n_z-2=C_j(n_z/C_j-\delta_{\min})+(\delta_{\min}C_j-2)$, both terms $\ge0$ "
            r"($n_z$ non-zero cells)",
            "",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
