# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V18: D-, E- and A-optimal gene panels (vignette 03-1, #sec-optimal-design;
note sec:decovart-genes and app:fisher, which give the J = 2 case, V4a/V4d).

Model: single-cell bulk law y ~ N_G(c M p, phi S_W(p)) with
S_W(p) = sum_j p_j Sigma_{W,j}; M is the G x J matrix of type means of the
panel. With a common within-type covariance Sigma_{W,j} = Sigma for all j,
S_W(p) = Sigma because sum_j p_j = 1. The free coordinates are theta with
p = p0 + K theta, K a J x (J-1) matrix whose columns are orthogonal to 1...
(any basis of the contrast space).

Derivation route:
- V18a: the information for theta is minus the expected Hessian of the
  Gaussian log-likelihood, here computed symbolically (G = 2, J = 3, generic
  M, generic SPD Sigma, generic theta); with a covariance free of theta the
  Hessian does not depend on y. Result: (c^2/phi) K' M' Sigma^{-1} M K. With
  J = 2 and K = (1, -1)' it is (c^2/phi) Delta' Sigma^{-1} Delta,
  Delta = mu_1 - mu_2 (the note's V4d with c^2/phi = n_i).
- V18b: with Sigma = I the information is (c^2/phi) A'A, A = M K: its
  determinant is the Gram determinant |a1|^2 |a2|^2 - (a1.a2)^2 (squared
  area spanned by the contrasted mean columns; generic 3 x 2 A), and its
  eigenvalues are the squared singular values of A, so
  kappa(I) = kappa(A)^2 (checked exactly on the example panels).
- V18c: Cov(p-hat) = K I_K^{-1} K' does not depend on the basis K
  (K -> K B, generic invertible B); det I_{KB} = det(B)^2 det I_K, so the
  D-ranking of panels is basis-free; the E- and A-criteria are basis-free
  only when computed on Cov(p-hat) (equivalently with orthonormal K). With
  the raw coordinates (p1, p2) the E-ranking of two panels can reverse:
  exact counter-example.
- V18d: worked example of the vignette, J = 3, G = 3, n = c^2/phi = 1,
  orthonormal (Helmert) K. Panel A: M_A rows (4,2,2), (4,2,0), (2,1,4);
  panel B: M_B rows (2,4,2), (0,2,4), (0,2,3); unit variances. With
  Sigma = I panel A wins on D, E and A. With correlation rho = 4/5 between
  genes 1 and 3 of panel A, D still prefers A, but E and A prefer B.
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


def zero(expr) -> bool:
    return sp.simplify(sp.cancel(sp.together(sp.expand(expr)))) == 0


HELMERT = sp.Matrix(
    [
        [1 / sp.sqrt(2), 1 / sp.sqrt(6)],
        [-1 / sp.sqrt(2), 1 / sp.sqrt(6)],
        [0, -2 / sp.sqrt(6)],
    ]
)
RAW = sp.Matrix([[1, 0], [0, 1], [-1, -1]])  # (p1, p2) free, p3 = 1 - p1 - p2


def info(M, S, K, n=1):
    return sp.simplify(n * K.T * M.T * S.inv() * M * K)


def criteria(I):
    ev = [sp.nsimplify(e) for e in I.eigenvals(multiple=True)]
    ev = sorted(ev, key=lambda e: float(e))
    return {
        "D": sp.simplify(I.det()),
        "E": sp.simplify(ev[0]),  # lambda_min of I
        "A": sp.simplify((I.inv()).trace()),
        "kappa": sp.simplify(ev[-1] / ev[0]),
    }


def run() -> list[dict]:
    out = []

    # ---- V18a ----------------------------------------------------------
    G, J = 2, 3
    M = sp.Matrix(G, J, lambda g, j: sp.Symbol(f"m{g+1}{j+1}", real=True))
    s11, s22 = sp.symbols("s11 s22", positive=True)
    s12 = sp.Symbol("s12", real=True)
    Sig = sp.Matrix([[s11, s12], [s12, s22]])
    c, phi = sp.symbols("c phi", positive=True)
    th = sp.Matrix(sp.symbols("t1 t2", real=True))
    K = RAW
    p = sp.Matrix([0, 0, 1]) + K * th
    S_W = sum((p[j] * Sig for j in range(J)), sp.zeros(G, G))
    ok_SW = all(zero(e) for e in (S_W - Sig))
    y = sp.Matrix(sp.symbols("y1 y2", real=True))
    mean = c * M * p
    V = phi * Sig
    ll = -(y - mean).T * V.inv() * (y - mean) / 2
    ll = ll[0, 0] - sp.log(V.det()) / 2
    H = sp.Matrix(2, 2, lambda a, b: sp.diff(ll, th[a], th[b]))
    free_of_y = all(sp.diff(H[a, b], yy) == 0 for a in range(2) for b in range(2) for yy in y)
    I_derived = -H
    I_claim = (c**2 / phi) * K.T * M.T * Sig.inv() * M * K
    ok_I = all(zero(e) for e in (I_derived - I_claim))
    # J = 2 reduction
    M2 = sp.Matrix(G, 2, lambda g, j: sp.Symbol(f"q{g+1}{j+1}", real=True))
    K2 = sp.Matrix([1, -1])
    Delta = M2[:, 0] - M2[:, 1]
    ok_J2 = zero(((K2.T * M2.T * Sig.inv() * M2 * K2) - Delta.T * Sig.inv() * Delta)[0, 0])
    ok_a = ok_SW and free_of_y and ok_I and ok_J2
    out.append(
        _rec(
            "V18a",
            r"With a common within-type covariance, $\bm S_W(\bm p)=\bSigma_W$ and the information for the free "
            r"coordinates $\bm\theta$ ($\bm p=\bm p_0+\bm K\bm\theta$) is $n_i\bm K\T\bm M\T\bSigma_W^{-1}\bm M\bm K$; "
            r"for $J=2$ it is $n_i\bm\Delta\T\bSigma_W^{-1}\bm\Delta$",
            "exact",
            "pass" if ok_a else "fail",
            r"Expected negative Hessian of the Gaussian log-likelihood, $G=2$, $J=3$, generic $\bm M$ and $\bSigma_W$; "
            r"Hessian free of $\bm y$; $n_i=c^2/\phi$",
            "",
            detail_full=f"S_W=Sigma {ok_SW}; free of y {free_of_y}; info {ok_I}; J=2 {ok_J2}",
        )
    )

    # ---- V18b ----------------------------------------------------------
    A = sp.Matrix(3, 2, lambda i, j: sp.Symbol(f"a{i+1}{j+1}", real=True))
    a1, a2 = A[:, 0], A[:, 1]
    gram = (a1.T * a1)[0, 0] * (a2.T * a2)[0, 0] - ((a1.T * a2)[0, 0]) ** 2
    ok_gram = zero((A.T * A).det() - gram)
    MA = sp.Matrix([[4, 2, 2], [4, 2, 0], [2, 1, 4]])
    MB = sp.Matrix([[2, 4, 2], [0, 2, 4], [0, 2, 3]])
    ok_kappa = True
    kap = {}
    for name, Mx in (("A", MA), ("B", MB)):
        Ax = Mx * HELMERT
        ev = sorted(sp.Matrix(Ax.T * Ax).eigenvals(multiple=True), key=float)
        sv = sorted([sp.sqrt(e) for e in ev], key=float)
        kap[name] = sp.nsimplify(sp.simplify(sv[-1] / sv[0]))
        ok_kappa &= zero(ev[-1] / ev[0] - (sv[-1] / sv[0]) ** 2)
    ok_b = ok_gram and ok_kappa
    out.append(
        _rec(
            "V18b",
            r"With $\bSigma_W=\bm I$ the information is $n_i\bm A\T\bm A$, $\bm A=\bm M\bm K$: its determinant is the "
            r"squared Gram volume of the contrasted mean columns and $\kappa(\bm I)=\kappa(\bm A)^2$",
            "exact",
            "pass" if ok_b else "fail",
            r"Lagrange identity $\det\bm A\T\bm A=\lVert\bm a_1\rVert^2\lVert\bm a_2\rVert^2-(\bm a_1\T\bm a_2)^2$ for "
            rf"generic $3\times2$ $\bm A$; $\kappa(\bm A)$ = {float(kap['A']):.3f} (A), {float(kap['B']):.3f} (B)",
            "",
        )
    )

    # ---- V18c ----------------------------------------------------------
    Ik = sp.Matrix([[sp.Symbol("i11", positive=True), sp.Symbol("i12")], [sp.Symbol("i12"), sp.Symbol("i22", positive=True)]])
    B = sp.Matrix(2, 2, lambda i, j: sp.Symbol(f"b{i+1}{j+1}", real=True))
    Kg = sp.Matrix(3, 2, lambda i, j: sp.Symbol(f"k{i+1}{j+1}", real=True))
    I_KB = B.T * Ik * B  # information for theta' with p = p0 + K B theta'
    cov_K = Kg * Ik.inv() * Kg.T
    cov_KB = (Kg * B) * I_KB.inv() * (Kg * B).T
    ok_cov = all(zero(e) for e in (cov_K - cov_KB))
    ok_det = zero(I_KB.det() - B.det() ** 2 * Ik.det())
    # counter-example: raw (p1, p2) against orthonormal coordinates
    MA3 = sp.Matrix([[0, 4, 3], [1, 2, 4], [3, 2, 4]])
    MB3 = sp.Matrix([[2, 1, 0], [1, 4, 3], [0, 4, 1]])
    S3 = sp.eye(3)
    S3[0, 2] = S3[2, 0] = sp.Rational(-4, 5)
    E_raw = (criteria(info(MA3, S3, RAW))["E"], criteria(info(MB3, sp.eye(3), RAW))["E"])
    E_ort = (criteria(info(MA3, S3, HELMERT))["E"], criteria(info(MB3, sp.eye(3), HELMERT))["E"])
    reversed_ = (float(E_raw[0]) < float(E_raw[1])) and (float(E_ort[0]) > float(E_ort[1]))
    ok_c = ok_cov and ok_det and reversed_
    out.append(
        _rec(
            "V18c",
            r"$\operatorname{Cov}\hat{\bm p}=\bm K\bm I_{\bm K}^{-1}\bm K\T$ does not depend on the basis $\bm K$ and "
            r"$\det\bm I_{\bm K\bm B}=(\det\bm B)^2\det\bm I_{\bm K}$; E and A computed in the raw coordinates "
            r"$(p_1,p_2)$ are basis-dependent",
            "exact",
            "pass" if ok_c else "fail",
            rf"Generic $\bm K$, $\bm B$, $\bm I_{{\bm K}}$; counter-example: $\lambda_{{\min}}$ of two panels "
            rf"{float(E_raw[0]):.2f} vs {float(E_raw[1]):.2f} in $(p_1,p_2)$, {float(E_ort[0]):.2f} vs "
            rf"{float(E_ort[1]):.2f} with orthonormal $\bm K$",
            "",
        )
    )

    # ---- V18d ----------------------------------------------------------
    rho = sp.Rational(4, 5)
    SA = sp.eye(3)
    SA[0, 2] = SA[2, 0] = rho
    cA0 = criteria(info(MA, sp.eye(3), HELMERT))
    cB0 = criteria(info(MB, sp.eye(3), HELMERT))
    cA1 = criteria(info(MA, SA, HELMERT))
    base = cA0["D"] > cB0["D"] and cA0["E"] > cB0["E"] and cA0["A"] < cB0["A"]
    split = cA1["D"] > cB0["D"] and cA1["E"] < cB0["E"] and cA1["A"] > cB0["A"]
    ok_d = bool(base and split)

    def fmt(cv):
        return f"D ${sp.latex(cv['D'])}$, E ${sp.latex(cv['E'])}$, A ${sp.latex(cv['A'])}$"

    out.append(
        _rec(
            "V18d",
            r"Example ($J=3$, $G=3$): with $\bSigma_W=\bm I$ panel A beats panel B on D, E and A; with correlation "
            r"$4/5$ between genes 1 and 3 of panel A, D still prefers A while E and A prefer B",
            "exact",
            "pass" if ok_d else "fail",
            rf"Orthonormal $\bm K$, $n_i=1$. A, $\rho=0$: {fmt(cA0)}; B: {fmt(cB0)}; A, $\rho=4/5$: {fmt(cA1)}",
            "",
            detail_full=json.dumps({k: str(v) for k, v in {"A0": cA0, "B0": cB0, "A1": cA1}.items()}),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
