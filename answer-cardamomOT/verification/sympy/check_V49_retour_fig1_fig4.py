# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0", "numpy==2.5.4"]
# ///
"""V49: claims of the note on figures 1 and 4
(review-decovart-bastien-elias-marielle-anais/retour_figure_1_and_4.qmd).

- V49a (exact): Gamma(alpha, beta) prior, Poisson counts with exposure L,
  zero count: posterior mean alpha / (beta + L).
- V49b (numeric, flag): alpha = beta = 1, L = 25 gives 1/26 = 0.03846, not
  the 0.03829 of the ChatGPT draft.
- V49c (exact): two perfect markers, reference 100, bulk 20 and 80: the
  factors 0.2/pA and 0.8/(1 - pA) reproduce the bulk for every pA; pA = 0.4
  gives 1/2 and 4/3.
- V49d (exact): Welch's G with J = 2 equals the square of Welch's t.
- V49e (exact): shrinkage keeps variances, Sigma(lambda) =
  Delta^1/2 R(lambda) Delta^1/2, and maps eigenvalues w -> (1 - lambda) w + lambda.
- V49f (exact): participation ratio and entropy rank equal G for equal
  eigenvalues and 1 for a single non-zero eigenvalue.
- V49g (exact): independent bivariate normal genes, n cells: E r^2 = 1/(n - 1)
  (n = 4..12), hence the null participation ratio
  G^2 / (G + G(G - 1)/(n - 1)) = G / (1 + (G - 1)/(n - 1)) (ratio of
  expectations).
- V49h (exact): X = s Y with s independent of Y:
  Cov(Xg, Xh) = E(s^2) Cov(Yg, Yh) + Var(s) mu_g mu_h.
- V49i (exact): two types, independent genes, y ~ N(c m(p), phi v(p)):
  I(p) = c^2 (muA - muB)^2/(phi v) + (1/2)((sA - sB)/v)^2; a perfect marker
  of A gives c^2 mu^2/(phi p s) + 1/(2 p^2).
- V49j (series): profile log-likelihood with a perfect marker of variance
  p s: l / log p -> (G - 1)/2 (so l -> -inf) when y > 0, and -> -1/2 (so
  l -> +inf) when y = 0.
- V49k (numeric): stageR level q R/G = 1.25 % for R = 2500, G = 10000;
  ORA example: E = 10, fold 3, hypergeometric tail P(X >= 30).
- V49l (numeric): 48 h RMSE of an all-epiblast prediction against the
  observed proportions is 0.02134.
"""

from __future__ import annotations

import sympy as sp


def _rec(id_, claim, kind, status, detail, tex_line=""):
    return {"id": id_, "claim": claim, "kind": kind, "status": status, "detail": detail,
            "detail_full": detail, "tex_line": tex_line}


def run() -> list[dict]:
    out = []

    # ---- V49a: Gamma-Poisson posterior after a zero count ---------------------
    a, b, L, r, x = sp.symbols("alpha beta L r x", positive=True)
    prior = b**a / sp.gamma(a) * r ** (a - 1) * sp.exp(-b * r)
    lik0 = sp.exp(-L * r)  # Poisson probability of zero counts at rate r L
    post_unnorm = prior * lik0
    z = sp.integrate(post_unnorm, (r, 0, sp.oo), conds="none")
    mean = sp.integrate(r * post_unnorm, (r, 0, sp.oo), conds="none") / z
    ok = sp.simplify(mean - a / (b + L)) == 0
    out.append(_rec("V49a", r"Gamma--Poisson, zero count: $\mathrm{E}(r\mid X=0)=\alpha/(\beta+L)$",
                    "exact", "pass" if ok else "fail", "integration of the posterior kernel, symbolic $\\alpha,\\beta,L$"))

    # ---- V49b: the draft's number ----------------------------------------------
    val = (a / (b + L)).subs({a: 1, b: 1, L: 25})
    draft = sp.Rational(3829, 100000)
    dev = abs(sp.N(val) - sp.N(draft)) > 1e-4
    out.append(_rec("V49b", r"$\alpha=\beta=1$, $L=25$: posterior mean $1/26$",
                    "numeric", "flag" if (val == sp.Rational(1, 26) and dev) else "fail",
                    f"1/26 = {sp.N(val, 6)}; the draft printed 0.03829"))

    # ---- V49c: identifiability counter-example ---------------------------------
    pA = sp.symbols("p_A", positive=True)
    bA, bB = sp.Rational(1, 5) / pA, sp.Rational(4, 5) / (1 - pA)
    ok = sp.simplify(100 * pA * bA - 20) == 0 and sp.simplify(100 * (1 - pA) * bB - 80) == 0
    ok = ok and bA.subs(pA, sp.Rational(2, 5)) == sp.Rational(1, 2) and bB.subs(pA, sp.Rational(2, 5)) == sp.Rational(4, 3)
    out.append(_rec("V49c", r"Free gene factors $b_A=0.2/p_A$, $b_B=0.8/(1-p_A)$ give bulk 20 and 80 for every $p_A$",
                    "exact", "pass" if ok else "fail", "symbolic $p_A$; $p_A=0.4$ gives $1/2$ and $4/3$"))

    # ---- V49d: Welch G for J = 2 ------------------------------------------------
    m1, m2, s1, s2, n1, n2 = sp.symbols("m1 m2 s1 s2 n1 n2", positive=True)
    w1, w2 = n1 / s1**2, n2 / s2**2
    W = w1 + w2
    mw = (w1 * m1 + w2 * m2) / W
    J = 2
    num = (w1 * (m1 - mw) ** 2 + w2 * (m2 - mw) ** 2) / (J - 1)
    den = 1 + sp.Rational(2 * (J - 2), J**2 - 1) * ((1 - w1 / W) ** 2 / (n1 - 1) + (1 - w2 / W) ** 2 / (n2 - 1))
    t2 = (m1 - m2) ** 2 / (s1**2 / n1 + s2**2 / n2)
    ok = sp.simplify(num / den - t2) == 0
    out.append(_rec("V49d", r"Welch's $G$ with $J=2$ equals Welch's $t^2$", "exact", "pass" if ok else "fail",
                    "symbolic means, variances and sizes"))

    # ---- V49e: shrinkage keeps variances and maps eigenvalues ------------------
    lam, rho, d1, d2 = sp.symbols("lambda rho d1 d2", positive=True)
    R = sp.Matrix([[1, rho], [rho, 1]])
    Dh = sp.diag(sp.sqrt(d1), sp.sqrt(d2))
    Sig = Dh * R * Dh
    Rl = (1 - lam) * R + lam * sp.eye(2)
    ok1 = sp.simplify(Dh * Rl * Dh - ((1 - lam) * Sig + lam * sp.diag(d1, d2))) == sp.zeros(2, 2)
    ev_R = sorted(R.eigenvals().keys(), key=str)
    ev_Rl = set(sp.simplify(e) for e in Rl.eigenvals().keys())
    ok2 = ev_Rl == set(sp.simplify((1 - lam) * e + lam) for e in ev_R)
    out.append(_rec("V49e", r"$\Delta^{1/2}\{(1-\lambda)R+\lambda I\}\Delta^{1/2}=(1-\lambda)\Sigma+\lambda\Delta$; eigenvalues $(1-\lambda)\omega+\lambda$",
                    "exact", "pass" if (ok1 and ok2) else "fail", "$2\\times2$, symbolic $\\rho$, $\\lambda$, variances"))

    # ---- V49f: limits of the two effective ranks --------------------------------
    def pr(w):
        return sum(w) ** 2 / sum(v**2 for v in w)

    def ent(w):
        t = sum(w)
        return sp.exp(-sum((v / t) * sp.log(v / t) for v in w if v != 0))

    G = 4
    ok = pr([1] * G) == G and sp.simplify(ent([1] * G) - G) == 0 and pr([G, 0, 0, 0]) == 1 and sp.simplify(ent([G, 0, 0, 0]) - 1) == 0
    out.append(_rec("V49f", r"$r_{\mathrm{PR}}$ and $r_{\mathrm H}$ equal $G$ for a flat spectrum and 1 for one direction",
                    "exact", "pass" if ok else "fail", "$G=4$"))

    # ---- V49g: E r^2 under independence, null participation ratio ---------------
    rr = sp.symbols("r", real=True)
    ok = True
    for n in range(4, 13):
        dens = (1 - rr**2) ** sp.Rational(n - 4, 2)
        er2 = sp.integrate(rr**2 * dens, (rr, -1, 1)) / sp.integrate(dens, (rr, -1, 1))
        ok = ok and sp.simplify(er2 - sp.Rational(1, n - 1)) == 0
    Gs, ns = sp.symbols("G n", positive=True)
    ok = ok and sp.simplify(Gs**2 / (Gs + Gs * (Gs - 1) / (ns - 1)) - Gs / (1 + (Gs - 1) / (ns - 1))) == 0
    out.append(_rec("V49g", r"Independent normal genes: $\mathrm E r^2=1/(n-1)$, so $r^0_{\mathrm{PR}}\approx G/\{1+(G-1)/(n-1)\}$",
                    "exact", "pass" if ok else "fail",
                    "null density of $r$, $n=4,\\ldots,12$; the rank is a ratio of expectations (approximation)"))

    # ---- V49h: capture-efficiency covariance -------------------------------------
    Es, Es2, mg, mh, cyy = sp.symbols("Es Es2 mu_g mu_h c_gh")
    exy = cyy + mg * mh  # E(Yg Yh)
    cov_x = Es2 * exy - (Es * mg) * (Es * mh)  # E(s^2) E(YgYh) - E(Xg) E(Xh)
    ok = sp.simplify(cov_x - (Es2 * cyy + (Es2 - Es**2) * mg * mh)) == 0
    out.append(_rec("V49h", r"$X=sY$, $s\perp Y$: $\mathrm{Cov}(X_g,X_h)=\mathrm E(s^2)\mathrm{Cov}(Y_g,Y_h)+\mathrm{Var}(s)\mu_g\mu_h$",
                    "exact", "pass" if ok else "fail", "moments of $s$ symbolic"))

    # ---- V49i: Fisher information, two types, independent genes --------------
    p, c, phi, muA, muB, sA, sB, yv = sp.symbols("p c phi mu_A mu_B s_A s_B y", positive=True)
    m = c * (p * muA + (1 - p) * muB)
    V = phi * (p * sA + (1 - p) * sB)
    loglik = -sp.log(V) / 2 - (yv - m) ** 2 / (2 * V)
    d2 = sp.diff(loglik, p, 2)
    # E over y ~ N(m, V): replace (y - m) and (y - m)^2 by their moments
    e = sp.symbols("e")
    d2e = sp.expand(d2.subs(yv, m + e))
    fisher = -sp.expand(d2e).subs(e**2, V).subs(e, 0)
    v = p * sA + (1 - p) * sB
    target = c**2 * (muA - muB) ** 2 / (phi * v) + sp.Rational(1, 2) * ((sA - sB) / v) ** 2
    ok1 = sp.simplify(fisher - target) == 0
    marker = target.subs({muB: 0, sB: 0})
    ok2 = sp.simplify(marker - (c**2 * muA**2 / (phi * p * sA) + 1 / (2 * p**2))) == 0
    out.append(_rec("V49i", r"$I(p_A)=c^2(\mu_A-\mu_B)^2/(\phi v)+\tfrac12\{(\sigma^2_A-\sigma^2_B)/v\}^2$; perfect marker: $c^2\mu^2/(\phi p\sigma^2)+1/(2p^2)$",
                    "exact", "pass" if (ok1 and ok2) else "fail", "$-\\mathrm E\\,\\partial^2_p\\ell$ for one gene"))

    # ---- V49j: barrier or unbounded supremum ------------------------------------
    eps_p, s, mu, Q0, Gn, yy = sp.symbols("p s mu Q0 G y", positive=True)
    def prof(yval, Gval):
        q = Q0 + (yval - eps_p * mu) ** 2 / (eps_p * s)
        return -sp.Rational(Gval, 2) * sp.log(q / Gval) - sp.log(eps_p * s) / 2
    lim_pos = sp.limit(prof(yy, 5) / sp.log(eps_p), eps_p, 0, "+")
    lim_zero = sp.limit(prof(0, 5) / sp.log(eps_p), eps_p, 0, "+")
    ok = sp.simplify(lim_pos - 2) == 0 and sp.simplify(lim_zero + sp.Rational(1, 2)) == 0
    out.append(_rec("V49j", r"Perfect marker of variance $p\sigma^2$: $\ell/\log p\to(G-1)/2$ if $y_g>0$, $\to-1/2$ if $y_g=0$",
                    "series", "pass" if ok else "fail",
                    f"$G=5$, other genes $Q_0>0$: limits {lim_pos} and {lim_zero}"))

    # ---- V49k: stageR level and ORA example -------------------------------------
    alpha_i = sp.Rational(5, 100) * 2500 / 10000
    N, K, nsel, k = 10000, 200, 500, 30
    tail = sum(sp.binomial(K, i) * sp.binomial(N - K, nsel - i) for i in range(k, min(nsel, K) + 1)) / sp.binomial(N, nsel)
    ok = alpha_i == sp.Rational(1, 80) and sp.Rational(nsel * K, N) == 10
    out.append(_rec("V49k", r"$\alpha_I=qR/G=1.25\,\%$; ORA $E=nK/N=10$, fold 3",
                    "numeric", "pass" if ok else "fail", f"hypergeometric tail $P(X\\ge 30)={sp.N(tail, 4)}$"))

    # ---- V49l: 48 h RMSE of an all-epiblast prediction --------------------------
    obs = [0.972362501147737, 0.0243320172619594, 0.00330548159030392]
    pred = [1, 0, 0]
    rmse = (sum((o - q) ** 2 for o, q in zip(obs, pred)) / 3) ** 0.5
    ok = abs(rmse - 0.021344825) < 1e-6
    out.append(_rec("V49l", r"48 h: predicting 100 \% epiblast gives RMSE 0.02134", "numeric",
                    "pass" if ok else "fail", f"RMSE = {rmse:.6f}; lambda\\_selection.csv minimum 0.021345"))
    return out


if __name__ == "__main__":
    for r_ in run():
        print(r_["id"], r_["status"], r_["detail"])
