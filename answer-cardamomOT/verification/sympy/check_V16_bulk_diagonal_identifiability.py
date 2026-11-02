# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V16: identifiability of the bulk sampling diagonal against the variance
scale phi_i (article sec:cov-estimation and app:which-covariance; vignette
01-2, final choice; extension of V13).

Model, given the composition p, the reference means mu~ and the denoised
within-type covariance S = S^bio(p) (all treated as known), with
m = mu~ p (entries > 0) and D = diag(m):

  article mean model (eq:sc-model, eq:bulk-poisson-layer, V13b):
      y ~ N_G(c m, phi S + c D),  c = a n, phi = a c = a^2 n;
  covariance-only map (a free coefficient k on the diagonal, e.g. a
  read-count dispersion omega c, or a depth not tied to the mean):
      (phi, k) -> phi S + k D.

Derivation route:
- V16a: Jacobian of (phi, k) -> vec(phi S + k D) is [vec S, vec D]; its
  2 x 2 minors, for generic symmetric S at G = 2 and G = 3, vanish exactly
  when S = kappa D (S diagonal with S_gg / m_g constant). Counter-example
  with S = D / 2: (phi, k) = (4, 6) and (8, 4) give the same covariance, and
  the Fisher information of (phi, k), derived from the expected negative
  Hessian of the Gaussian log-likelihood (moments of y - mean substituted
  monomial by monomial), is singular there and regular for a non-proportional
  S.
- V16b: with the article's mean model, the Fisher information of (c, phi),
  derived the same way at G = 2, equals I_mean + I_cov with
  I_mean = diag(m' Sigma^{-1} m, 0); hence det I = (m' Sigma^{-1} m) I_phiphi +
  det I_cov, positive whenever S != 0: (c, phi) is identified even when
  S = kappa D. Explicit determinant in the degenerate case. The map
  (phi, a) -> (c, phi) = (phi / a, phi) has Jacobian determinant phi / a^2,
  so (phi, a) is identified under the same condition.
- V16c: cost. With phi = a c, Sigma = c M, M = a S + D. For fixed (p, a) the
  score in c times 2 c^2 is A - G c - C c^2 (A = y' M^{-1} y,
  C = m' M^{-1} m), so c-hat is the unique positive root of
  C c^2 + G c - A = 0 and a strict maximum; a is left to a one-dimensional
  numerical search. Without the diagonal the article's closed forms
  c-hat (GLS) and phi-hat = Q / G hold (checked as the D = 0 case).
"""

from __future__ import annotations

import itertools
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


def _sym(G, name):
    """Generic symmetric G x G matrix of real symbols."""
    ent = {}
    for g in range(G):
        for h in range(g, G):
            ent[(g, h)] = ent[(h, g)] = sp.Symbol(f"{name}{g + 1}{h + 1}", real=True)
    return sp.Matrix(G, G, lambda g, h: ent[(g, h)]), ent


def _fisher(params, mean, Sig):
    """Expected negative Hessian of the Gaussian log-likelihood, derived from
    first principles: differentiate l(theta; y) twice, write y = mean + r,
    then replace each monomial in r by its moment (E r = 0, E r r' = Sig)."""
    G = Sig.shape[0]
    y = sp.Matrix(sp.symbols(f"y1:{G + 1}", real=True))
    r = sp.symbols(f"r1:{G + 1}", real=True)
    res = y - mean
    ll = -sp.log(Sig.det()) / 2 - (res.T * Sig.inv() * res)[0] / 2
    k = len(params)
    info = sp.zeros(k, k)
    sub = {y[g]: mean[g] + r[g] for g in range(G)}
    for i in range(k):
        for j in range(i, k):
            h = -sp.diff(ll, params[i], params[j])
            h = sp.together(h.subs(sub))
            num, den = sp.fraction(h)
            poly = sp.Poly(sp.expand(num), *r)
            e = 0
            for mon, coef in poly.terms():
                deg = sum(mon)
                if deg == 0:
                    e += coef
                elif deg == 2:
                    idx = [g for g, d in enumerate(mon) for _ in range(d)]
                    e += coef * Sig[idx[0], idx[1]]
                elif deg != 1:
                    raise ValueError("Hessian is not quadratic in y")
            info[i, j] = info[j, i] = sp.simplify(e / den)
    return info


def _closed_fisher(params, mean, Sig):
    """Textbook form: dmu' Sig^-1 dmu + tr(Sig^-1 dSig Sig^-1 dSig) / 2."""
    Si = Sig.inv()
    k = len(params)
    return sp.Matrix(
        k,
        k,
        lambda i, j: sp.simplify(
            (sp.diff(mean, params[i]).T * Si * sp.diff(mean, params[j]))[0]
            + (Si * sp.diff(Sig, params[i]) * Si * sp.diff(Sig, params[j])).trace() / 2
        ),
    )


def run() -> list[dict]:
    out: list[dict] = []
    R = sp.Rational

    # ---------------- V16a: covariance-only map (phi, k) ----------------
    ok_minors = True
    minor_detail = []
    for G in (2, 3):
        S, ent = _sym(G, "s")
        m = sp.symbols(f"m1:{G + 1}", positive=True)
        vecS = [S[g, h] for g in range(G) for h in range(g, G)]
        vecD = [m[g] if g == h else 0 for g in range(G) for h in range(g, G)]
        J = sp.Matrix([vecS, vecD]).T
        minors = [sp.expand(J[[i, j], :].det()) for i, j in itertools.combinations(range(J.rows), 2)]
        unknowns = [ent[(g, h)] for g in range(G) for h in range(g, G) if (g, h) != (0, 0)]
        sols = sp.solve(minors, unknowns, dict=True)
        kappa = ent[(0, 0)] / m[0]
        target = {ent[(g, h)]: (kappa * m[g] if g == h else 0) for g in range(G) for h in range(g, G) if (g, h) != (0, 0)}
        ok_G = len(sols) == 1 and all(sp.simplify(sols[0][u] - target[u]) == 0 for u in unknowns)
        # and S = kappa D does make every minor vanish
        kap = sp.Symbol("kappa", positive=True)
        degen = {ent[(g, h)]: (kap * m[g] if g == h else 0) for g in range(G) for h in range(g, G)}
        ok_G = ok_G and all(sp.simplify(mi.subs(degen)) == 0 for mi in minors)
        ok_minors &= ok_G
        minor_detail.append(f"G={G}: minors vanish iff {sols}")

    phi, k = sp.symbols("phi k", positive=True)
    m2 = sp.Matrix([2, 3])
    D2 = sp.diag(*m2)
    S_deg = D2 / 2
    cov = lambda ph, kk, S: ph * S + kk * D2  # noqa: E731
    ok_counter = cov(4, 6, S_deg) == cov(8, 4, S_deg) == sp.diag(16, 24)
    zero = sp.zeros(2, 1)
    I_deg = _fisher([phi, k], zero, cov(phi, k, S_deg))
    I_deg_closed = _closed_fisher([phi, k], zero, cov(phi, k, S_deg))
    S_gen = sp.Matrix([[1, R(1, 2)], [R(1, 2), 2]])
    I_gen = _fisher([phi, k], zero, cov(phi, k, S_gen))
    det_deg = sp.simplify(I_deg.det())
    det_gen = sp.simplify(I_gen.subs({phi: 4, k: 6}).det())
    ok_fisher_form = all(sp.simplify(e) == 0 for e in I_deg - I_deg_closed)
    ok_a = ok_minors and ok_counter and det_deg == 0 and det_gen > 0 and ok_fisher_form
    out.append(
        _rec(
            "V16a",
            r"Covariance alone, given $\bm p$: $(\phi,k)\mapsto\phi\,\bm S^{\mathrm{bio}}(\bm p)+k\operatorname{diag}"
            r"(\tilde{\bm\mu}\bm p)$ is identifiable iff $\bm S^{\mathrm{bio}}(\bm p)$ is not a multiple of "
            r"$\operatorname{diag}(\tilde{\bm\mu}\bm p)$",
            "exact",
            "pass" if ok_a else "fail",
            r"Jacobian $[\operatorname{vec}\bm S,\operatorname{vec}\bm D]$: all $2\times2$ minors vanish iff "
            r"$\bm S=\kappa\bm D$ ($G=2,3$, generic symmetric $\bm S$); $\bm S=\bm D/2$, $\bm m=(2,3)$: "
            r"$(\phi,k)=(4,6)$ and $(8,4)$ give $\operatorname{diag}(16,24)$; Fisher determinant $0$ there, "
            r"$" + sp.latex(det_gen) + r"$ for a non-proportional $\bm S$",
            "",
            detail_full=(
                "; ".join(minor_detail)
                + f"; counter-example equal: {ok_counter}; Fisher (phi,k) at S = D/2: {sp.sstr(I_deg)}, "
                f"det {det_deg}; derived = textbook form: {ok_fisher_form}; det at S = [[1,1/2],[1/2,2]], "
                f"(phi,k) = (4,6): {det_gen}"
            ),
        )
    )

    # ---------------- V16b: the article's mean model (c, phi) -----------
    c = sp.Symbol("c", positive=True)
    S, _ = _sym(2, "s")
    msym = sp.Matrix(sp.symbols("m1 m2", positive=True))
    Dm = sp.diag(*msym)
    Sig = phi * S + c * Dm
    I_full = _fisher([c, phi], c * msym, Sig)
    I_closed = _closed_fisher([c, phi], c * msym, Sig)
    I_cov = _closed_fisher([c, phi], sp.zeros(2, 1), Sig)
    q = sp.simplify((msym.T * Sig.inv() * msym)[0])
    # cancel(together(.)) decides zero for rational functions and, unlike
    # simplify, returns quickly when the difference is not zero
    _zero = lambda e: sp.cancel(sp.together(e)) == 0  # noqa: E731
    ok_form = all(_zero(e) for e in I_full - I_closed)
    ok_split = all(_zero(e) for e in I_full - I_cov - sp.diag(q, 0))
    # given the split, det(I_cov + diag(q, 0)) = q I_phiphi + det I_cov is 2 x 2
    # algebra; checked on abstract entries (the expanded determinant is large)
    i11, i12, i22, qq = sp.symbols("i11 i12 i22 qq", real=True)
    Iab = sp.Matrix([[i11, i12], [i12, i22]])
    ok_detid = ok_split and sp.expand((Iab + sp.diag(qq, 0)).det() - (qq * i22 + Iab.det())) == 0
    # degenerate case S = kappa D: covariance part singular, full information regular
    kap = sp.Symbol("kappa", positive=True)
    degen = {S[0, 0]: kap * msym[0], S[1, 1]: kap * msym[1], S[0, 1]: 0}
    det_cov_deg = sp.simplify(I_cov.subs(degen).det())
    det_full_deg = sp.factor(sp.simplify(I_full.subs(degen).det()))
    ok_deg = det_cov_deg == 0 and det_full_deg.is_positive is True
    # (phi, a) parametrisation: c = phi / a
    a = sp.Symbol("a", positive=True)
    jac = sp.Matrix([phi / a, phi]).jacobian([phi, a])
    det_jac = sp.simplify(jac.det())
    ok_jac = sp.simplify(det_jac - phi / a**2) == 0
    ok_b = ok_form and ok_split and ok_detid and ok_deg and ok_jac
    out.append(
        _rec(
            "V16b",
            r"With the mean $c\,\tilde{\bm\mu}\bm p$ of the article and $\operatorname{Cov}\bm y=\phi\,\bm S^{\mathrm{bio}}"
            r"(\bm p)+c\operatorname{diag}(\tilde{\bm\mu}\bm p)$, $(c,\phi)$, and so $(\phi,a)$, is identified whenever "
            r"$\bm S^{\mathrm{bio}}(\bm p)\neq\bm 0$, also when it is a multiple of the diagonal",
            "exact",
            "pass" if ok_b else "fail",
            r"Fisher information of $(c,\phi)$ from the expected Hessian, $G=2$: $\bm I=\bm I_{\mathrm{cov}}+"
            r"\operatorname{diag}(\bm m^\top\bm\Sigma^{-1}\bm m,0)$, $\det\bm I=(\bm m^\top\bm\Sigma^{-1}\bm m)I_{\phi\phi}"
            r"+\det\bm I_{\mathrm{cov}}$; at $\bm S=\kappa\bm D$: $\det\bm I_{\mathrm{cov}}=0$, $\det\bm I="
            + sp.latex(det_full_deg)
            + r"$; $\partial(c,\phi)/\partial(\phi,a)$ has determinant $\phi/a^2$",
            "",
            detail_full=(
                f"derived = textbook form: {ok_form}; split I = I_cov + diag(m'Sigma^-1 m, 0): {ok_split}; "
                f"det identity: {ok_detid}; S = kappa D: det I_cov = {det_cov_deg}, det I = {det_full_deg}; "
                f"Jacobian determinant (phi, a) -> (c, phi): {det_jac}"
            ),
        )
    )

    # ---------------- V16c: profile of c for fixed a --------------------
    Mm, _ = _sym(2, "M")
    yv = sp.Matrix(sp.symbols("y1 y2", real=True))
    G = 2
    Mi = Mm.inv()
    A = (yv.T * Mi * yv)[0]
    C = (msym.T * Mi * msym)[0]
    ll = -G * sp.log(c) / 2 - sp.log(Mm.det()) / 2 - ((yv - c * msym).T * Mi * (yv - c * msym))[0] / (2 * c)
    score = sp.diff(ll, c)
    ok_score = sp.cancel(sp.together(2 * c**2 * score - (A - G * c - C * c**2))) == 0
    Asym, Csym, Gs = sp.symbols("A C G", positive=True)
    roots = sp.solve(sp.Eq(Csym * c**2 + Gs * c - Asym, 0), c)
    c_hat = (-Gs + sp.sqrt(Gs**2 + 4 * Asym * Csym)) / (2 * Csym)
    ok_root = len(roots) == 1 and sp.simplify(roots[0] - c_hat) == 0  # c positive: one root
    cc = sp.Symbol("cc", positive=True)
    ll_gen = -Gs * sp.log(cc) / 2 - Asym / (2 * cc) - Csym * cc / 2
    d2 = sp.diff(ll_gen, cc, 2)
    # at a root, A = C c^2 + G c; then d2 = -G / (2 c^2) - C / c < 0
    d2_at = sp.simplify(d2.subs(Asym, Csym * cc**2 + Gs * cc))
    ok_max = sp.simplify(d2_at - (-Gs / (2 * cc**2) - Csym / cc)) == 0
    # over all complex c the two roots multiply to -A/C < 0: exactly one is positive
    z = sp.Symbol("z")
    all_roots = sp.solve(Csym * z**2 + Gs * z - Asym, z)
    ok_prod = len(all_roots) == 2 and sp.simplify(all_roots[0] * all_roots[1] + Asym / Csym) == 0
    # D = 0: Sigma = phi S with phi free; closed forms of the article
    Sw, _ = _sym(2, "w")
    Wi = Sw.inv()
    ll0 = -G * sp.log(phi) / 2 - ((yv - c * msym).T * Wi * (yv - c * msym))[0] / (2 * phi)
    c_gls = sp.solve(sp.diff(ll0, c), c)
    ok_gls = len(c_gls) == 1 and sp.simplify(c_gls[0] - (msym.T * Wi * yv)[0] / (msym.T * Wi * msym)[0]) == 0
    Q = sp.Symbol("Q", positive=True)
    phi_hat = sp.solve(sp.diff(-G * sp.log(phi) / 2 - Q / (2 * phi), phi), phi)
    ok_q = phi_hat == [Q / G]
    ok_c = ok_score and ok_root and ok_max and ok_prod and ok_gls and ok_q
    out.append(
        _rec(
            "V16c",
            r"Cost: with $\phi=ac$ and $\bm M=a\bm S^{\mathrm{bio}}(\bm p)+\operatorname{diag}(\tilde{\bm\mu}\bm p)$, "
            r"$\hat c$ for fixed $(\bm p,a)$ is the positive root of $Cc^2+Gc-A=0$; $a$ needs a one-dimensional "
            r"numerical profile; without the diagonal, $\hat c$ is the GLS scale and $\hat\phi=Q/G$",
            "exact",
            "pass" if ok_c else "fail",
            r"$2c^2\,\partial_c\ell=A-Gc-Cc^2$, $A=\bm y^\top\bm M^{-1}\bm y$, $C=\bm m^\top\bm M^{-1}\bm m$ ($G=2$, "
            r"generic $\bm M$); one positive root; $\partial_c^2\ell=-G/(2c^2)-C/c<0$ there",
            "",
            detail_full=(
                f"score identity: {ok_score}; positive root {roots}: {ok_root}; second derivative at root "
                f"{d2_at}: {ok_max}; D = 0: GLS scale {ok_gls}, phi-hat {phi_hat}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
