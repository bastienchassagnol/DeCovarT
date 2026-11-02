# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V13: which within-type covariance enters the bulk, raw or Poisson-denoised
(vignette 02-1, sec-raw-vs-denoised; article app:which-covariance; note,
noise section).

Model. A cell of type j has a latent expression vector lambda_c on the scale
of the reference capture, with E lambda_c = mu_j and Cov lambda_c =
Sigma^bio_j; the reference counts are x_c | lambda_c ~ Poisson(lambda_c),
genes independent given lambda_c. The bulk pools m_j cells of type j
(independent cells) and has a depth a per reference unit. Three sampling
models for the bulk:

  (A) cell-depth sampling, then scaling: y = a * sum_c x_c;
  (B) bulk-depth sampling of the pooled transcripts: y | Lambda ~
      Poisson(a * Lambda), Lambda = sum_c lambda_c;
  (C) binomial thinning of the reference-depth counts (a <= 1): each counted
      molecule of each cell is kept with probability a, y = sum_c Bin(x_c, a).

Derivation route (G = 2, J = 2, no distributional form beyond a generic
cumulant generating function):
- the CGF of lambda for type j is written K_j(u) = mu_j'u + u'Sigma_j u / 2 +
  a generic cubic (symbolic third cumulants), standing for any law with
  finite third moments; independence of cells makes the CGF of the pool
  sum_j m_j K_j;
- the Poisson MGF E exp(s x) = exp(lambda (e^s - 1)) and the binomial MGF
  (1 - a + a e^t)^x are summed from their pmfs;
- model (A): K_y(t) = sum_j m_j K_j(e^{a t} - 1); model (B): K_y(t) =
  sum_j m_j K_j(a (e^t - 1)); model (C): the binomial MGF turns x into
  exp(x log(1 + a(e^t - 1))), and the Poisson MGF then gives the argument
  a(e^t - 1) of model (B), so (B) and (C) have the same CGF, not only the
  same covariance;
- mean and covariance are the gradient and Hessian of K_y at t = 0.
V13b: with m_j = n p_j, c = a n, phi = a^2 n, (B) is phi S^bio(p) +
c diag(mu p) = phi {S^bio(p) + a^{-1} diag(mu p)}; the raw plug-in phi S^raw(p)
matches it iff a = 1; the sampling term of gene g exceeds the biological one
iff mu_g(p) > a S^bio_gg(p).
V13c: the toy example of the vignette, exact rational numbers.
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


def _mgf_checks():
    """Poisson and binomial MGFs summed from the pmf."""
    k = sp.Symbol("k", integer=True, nonnegative=True)
    lam, s = sp.symbols("lam s", positive=True)
    # e^{sk} lam^k = z^k with z = lam e^s, so the series is an exponential series
    z = sp.Symbol("z", positive=True)
    pois = sp.summation(sp.exp(-lam) * z**k / sp.factorial(k), (k, 0, sp.oo)).subs(z, lam * sp.exp(s))
    ok_pois = sp.simplify(pois - sp.exp(lam * (sp.exp(s) - 1))) == 0
    a, t = sp.symbols("a t", positive=True)
    ok_bin = True
    for x in range(0, 7):
        mgf = sum(sp.binomial(x, kk) * a**kk * (1 - a) ** (x - kk) * sp.exp(t * kk) for kk in range(x + 1))
        ok_bin &= sp.expand(mgf - (1 - a + a * sp.exp(t)) ** x) == 0
    return ok_pois, ok_bin


def _cgf_type(j, u):
    """Generic CGF of lambda for type j up to third order (G = 2)."""
    mu = sp.symbols(f"mu{j}_1 mu{j}_2", positive=True)
    s11, s22 = sp.symbols(f"S{j}_11 S{j}_22", positive=True)
    s12 = sp.Symbol(f"S{j}_12", real=True)
    Sig = sp.Matrix([[s11, s12], [s12, s22]])
    kap = {
        idx: sp.Symbol(f"k{j}_" + "".join(map(str, idx)), real=True)
        for idx in itertools.combinations_with_replacement((1, 2), 3)
    }
    uv = sp.Matrix(u)
    K = (sp.Matrix(mu).T * uv)[0] + (uv.T * Sig * uv)[0] / 2
    K += sum(c * u[i - 1] * u[jj - 1] * u[kk - 1] for (i, jj, kk), c in kap.items()) / 6
    return K, sp.Matrix(mu), Sig


def run() -> list[dict]:
    out: list[dict] = []
    ok_pois, ok_bin = _mgf_checks()

    a = sp.Symbol("a", positive=True)
    t1, t2 = sp.symbols("t1 t2", real=True)
    t = (t1, t2)
    m1, m2 = sp.symbols("m1 m2", positive=True)

    def pool_cgf(argfun):
        tot = 0
        mus, sigs = [], []
        for j, m in ((1, m1), (2, m2)):
            u = [argfun(ti) for ti in t]
            K, mu, Sig = _cgf_type(j, u)
            tot += m * K
            mus.append(mu)
            sigs.append(Sig)
        return tot, mus, sigs

    KA, mus, sigs = pool_cgf(lambda ti: sp.exp(a * ti) - 1)
    KB, _, _ = pool_cgf(lambda ti: a * (sp.exp(ti) - 1))
    # (C): binomial thinning, then Poisson: argument exp(log(1 + a(e^t - 1))) - 1
    KC, _, _ = pool_cgf(lambda ti: sp.exp(sp.log(1 + a * (sp.exp(ti) - 1))) - 1)
    ok_bc = sp.simplify(KB - KC) == 0

    def moments(K):
        at0 = {t1: 0, t2: 0}
        mean = sp.Matrix([sp.diff(K, ti).subs(at0) for ti in t]).applyfunc(sp.simplify)
        cov = sp.Matrix(2, 2, lambda g, h: sp.simplify(sp.diff(K, t[g], t[h]).subs(at0)))
        return mean, cov

    meanA, covA = moments(KA)
    meanB, covB = moments(KB)
    mu_mix = m1 * mus[0] + m2 * mus[1]
    S_bio = m1 * sigs[0] + m2 * sigs[1]
    D = sp.diag(*mu_mix)
    ok_mean = all(sp.simplify(e) == 0 for e in list(meanA - a * mu_mix) + list(meanB - a * mu_mix))
    ok_A = all(sp.simplify(e) == 0 for e in covA - a**2 * (S_bio + D))
    ok_B = all(sp.simplify(e) == 0 for e in covB - (a**2 * S_bio + a * D))
    ok_a = ok_pois and ok_bin and ok_bc and ok_mean and ok_A and ok_B
    out.append(
        _rec(
            "V13a",
            r"Bulk covariance under cell-depth sampling then scaling, $a^2\sum_j m_j(\bm\Sigma^{\mathrm{bio}}_j"
            r"+\operatorname{diag}\bm\mu_j)$, against bulk-depth sampling, $a^2\sum_j m_j\bm\Sigma^{\mathrm{bio}}_j"
            r"+a\operatorname{diag}(\sum_j m_j\bm\mu_j)$; binomial thinning equals bulk-depth sampling",
            "exact",
            "pass" if ok_a else "fail",
            r"Poisson and binomial MGFs summed from the pmf; generic third-order CGF per type, $G=J=2$; "
            r"Hessian of the pooled CGF at $0$; thinning and bulk-depth sampling have the same CGF",
            "",
            detail_full=(
                f"Poisson MGF: {ok_pois}; binomial MGF x<=6: {ok_bin}; CGF(B) == CGF(C): {ok_bc}; "
                f"means a*sum m_j mu_j: {ok_mean}; Cov(A) = a^2 sum m_j (Sigma_j + diag mu_j): {ok_A}; "
                f"Cov(B) = a^2 sum m_j Sigma_j + a diag(sum m_j mu_j): {ok_B}; Cov(A) = {sp.sstr(covA)}; "
                f"Cov(B) = {sp.sstr(covB)}"
            ),
        )
    )

    # ---------------- V13b: per-cell scales and the plug-in -----------------
    n, p = sp.symbols("n p", positive=True)
    sub = {m1: n * p, m2: n * (1 - p)}
    c, phi = a * n, a**2 * n
    S_p = (p * sigs[0] + (1 - p) * sigs[1])
    mu_p = p * mus[0] + (1 - p) * mus[1]
    covB_p = covB.subs(sub)
    ok_scales = all(
        sp.simplify(e) == 0 for e in covB_p - (phi * S_p + c * sp.diag(*mu_p))
    ) and all(sp.simplify(e) == 0 for e in covB_p - phi * (S_p + sp.diag(*mu_p) / a))
    raw_plugin = phi * (S_p + sp.diag(*mu_p))
    gap = sp.factor(sp.simplify((raw_plugin - covB_p)[0, 0]))
    sol = sp.solve(sp.Eq(raw_plugin[0, 0], covB_p[0, 0]), a)
    ok_raw = sol == [1] and all(sp.simplify(e) == 0 for e in (raw_plugin - covB_p).subs(a, 1))
    ok_offdiag = sp.simplify((raw_plugin - covB_p)[0, 1]) == 0
    # denoised alone is the a -> oo limit of the per-phi covariance
    lims = [sp.limit(sp.simplify(e), a, sp.oo) for e in (covB_p / phi - S_p)]
    ok_lim = all(sp.simplify(v) == 0 for v in lims)
    # sampling > biology for gene g iff mu_g > a S_gg
    g_samp, g_bio = c * mu_p[0], phi * S_p[0, 0]
    ok_cross = sp.simplify(g_samp / g_bio - mu_p[0] / (a * S_p[0, 0])) == 0
    ok_gapform = sp.simplify(gap - a * n * (a - 1) * mu_p[0]) == 0
    ok_b = ok_scales and ok_raw and ok_offdiag and ok_lim and ok_cross and ok_gapform
    out.append(
        _rec(
            "V13b",
            r"With $m_j=np_j$: $\operatorname{Cov}\bm y=\phi\,\bm S^{\mathrm{bio}}(\bm p)+c\operatorname{diag}"
            r"(\tilde{\bm\mu}\bm p)=\phi\{\bm S^{\mathrm{bio}}(\bm p)+a^{-1}\operatorname{diag}(\tilde{\bm\mu}\bm p)\}$; "
            r"the raw plug-in $\phi\,\bm S^{\mathrm{raw}}(\bm p)$ is exact iff $a=1$; the denoised one only as "
            r"$a\to\infty$; sampling exceeds biology for gene $g$ iff $\tilde\mu_g(\bm p)>a S^{\mathrm{bio}}_{gg}(\bm p)$",
            "exact",
            "pass" if ok_b else "fail",
            r"raw minus bulk-depth covariance: $an(a-1)(\tilde{\bm\mu}\bm p)_g$ on the diagonal, $0$ off it; "
            r"its root in $a$ is $1$; limit of $\operatorname{Cov}\bm y/\phi-\bm S^{\mathrm{bio}}$ as $a\to\infty$ is $0$",
            "",
            detail_full=(
                f"scales: {ok_scales}; raw - bulk (1,1) = {gap}; roots in a: {sol}; off-diagonal equal: "
                f"{ok_offdiag}; a->oo limits {lims}: {ok_lim}; crossover: {ok_cross}"
            ),
        )
    )

    # ---------------- V13c: toy example, G = 2, one type -------------------
    R = sp.Rational
    mu = sp.Matrix([10, 1])
    Sb = sp.Matrix([[40, 2], [2, 1]])
    Sraw = Sb + sp.diag(*mu)
    rows = {}
    for aval in (R(1, 4), R(1), R(4)):
        true = Sb + sp.diag(*mu) / aval  # Cov y / phi
        corr = sp.nsimplify(true[0, 1] / sp.sqrt(true[0, 0] * true[1, 1]))
        rows[aval] = (true, corr, Sraw[1, 1] / true[1, 1], Sb[1, 1] / true[1, 1])
    expected = {
        R(1, 4): (sp.Matrix([[80, 2], [2, 5]]), R(1, 10), R(2, 5), R(1, 5)),
        R(1): (sp.Matrix([[50, 2], [2, 2]]), R(1, 5), R(1), R(1, 2)),
        R(4): (sp.Matrix([[R(85, 2), 2], [2, R(5, 4)]]), 2 / sp.sqrt(R(85, 2) * R(5, 4)), R(8, 5), R(4, 5)),
    }
    ok_c = all(
        rows[k][0] == v[0] and sp.simplify(rows[k][1] - v[1]) == 0 and rows[k][2] == v[2] and rows[k][3] == v[3]
        for k, v in expected.items()
    )
    cross = [mu[g] / Sb[g, g] for g in range(2)]
    ok_c = ok_c and cross == [R(1, 4), R(1)]
    corr4 = sp.N(rows[R(4)][1], 4)
    out.append(
        _rec(
            "V13c",
            r"Toy example: $\tilde{\bm\mu}=(10,1)$, $\bm\Sigma^{\mathrm{bio}}=\begin{pmatrix}40&2\\2&1\end{pmatrix}$; "
            r"$\operatorname{Cov}\bm y/\phi$ at $a=\tfrac14,1,4$; gene-2 variance ratios raw/true and "
            r"denoised/true; crossover depths $\tfrac14$ and $1$",
            "exact",
            "pass" if ok_c else "fail",
            r"$a=\tfrac14$: $\begin{pmatrix}80&2\\2&5\end{pmatrix}$, correlation $\tfrac1{10}$, ratios $\tfrac25$, $\tfrac15$; "
            r"$a=1$: raw, correlation $\tfrac15$, ratios $1$, $\tfrac12$; $a=4$: $\begin{pmatrix}42.5&2\\2&1.25\end{pmatrix}$, "
            r"correlation $" + str(corr4) + r"$, ratios $\tfrac85$, $\tfrac45$",
            "",
            detail_full=(
                "; ".join(
                    f"a={k}: Cov/phi={sp.sstr(v[0])}, corr={v[1]} ({sp.N(v[1], 6)}), raw/true={v[2]}, "
                    f"denoised/true={v[3]}"
                    for k, v in rows.items()
                )
                + f"; crossover depths mu_g/S_gg = {cross}"
            ),
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
