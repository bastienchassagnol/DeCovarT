# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V29: cumulant generating functions and Legendre transforms of the
bulk-versus-reference discrepancy D (vignette 03-1, #sec-fs-five-bounds,
@eq-fs-cgf-d, @eq-fs-legendre, @eq-fs-cgf-pois).

Setting (A1)-(A4): per-cell NB counts with mean mu_j and size psi,
pi_j = mu_j/(mu_j + psi); m_j = p_j n bulk cells of type j, C_j reference
cells; D = (1/n) sum_bulk X - sum_j (p_j/C_j) sum_ref_j X.

Derivation route:
- V29a: the per-cell MGF from the pmf weights: the coefficients of the series
  of (1 - z)^(-psi) equal C(x + psi - 1, x), so M(u) = (1 - pi)^psi
  (1 - pi e^u)^(-psi); its log is compared with the printed second form
  -psi log(1 - mu/psi (e^u - 1)), and the domain of that form with
  u < -log pi. The CGF of D is the log of the product of per-cell MGFs
  (independence); for a toy (J = 2, n = 4, m = (3, 1), C = (2, 3)) E e^{uD} is
  recomputed from the convolved NB pmfs (mpmath, 40 digits, truncated where
  the tail is below 1e-40) and compared with the printed K_D. Derivatives at
  0 with symbolic J = 3 parameters: K_D(0) = K_D'(0) = 0, K_D''(0) = Var D
  with sigma^2 = M''(0) - M'(0)^2, K_D'''(0) = kappa_3(D) with the per-cell
  third cumulant recomputed from M. K_gj'' > 0 (convexity). Skewness
  kappa_3 / Var^(3/2) at the parameters of fig_03_1_tail_bounds.R against
  the printed -0.08 and -0.06.
- V29b: Legendre transforms. u t - K_D(u) is concave; K_D'(u) -> +oo as
  u -> u_+ and -> -oo as u -> u_- (one-sided limits of each per-cell term at
  the printed endpoints), so for t > 0 the supremum over (0, u_+) is at the
  root of K_D'(u) = t, and over (u_-, 0) at the root of K_D'(u) = -t. For a
  single sum of n cells (no reference term) the stationary point and the
  value are solved symbolically and compared with n K*(m +- t) of
  @eq-fs-nb-rate; numeric root-finding (mpmath) at the spec grid confirms.
- V29c: the Poisson CGF mu (e^u - 1) by summing the Poisson series; the
  Poisson K_D^Pois has K(0) = K'(0) = 0; the single-sum rate is solved
  symbolically and compared with n {x log(x/m) - x + m}; the variance ratio
  NB / Poisson per cell is 1 + mu/psi; as psi -> oo the NB per-cell CGF tends
  to the Poisson CGF and the NB rate to the Poisson rate.
- V29d: module CGFs added over genes "as if independent": the CGF of
  Z_g = D_g / sd D_g is K_{D_g}(u / sd D_g), so the summed NB CGF has second
  derivative |M| at 0; the summed Poisson CGF with NB standard deviations has
  second derivative sum_g Var^Pois D_g / Var^NB D_g < |M| (difference
  sum_j w_j mu_j^2 / psi > 0).
"""

from __future__ import annotations

import json

import mpmath as mp
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


def zero(e) -> bool:
    if sp.cancel(sp.together(sp.expand(e))) == 0:
        return True
    return sp.simplify(e) == 0


u, v, t = sp.symbols("u v t", real=True)
psi = sp.Symbol("psi", positive=True)


def K_nb(arg, mu):
    """Printed per-cell CGF of @eq-fs-cgf-d."""
    return -psi * sp.log(1 - mu / psi * (sp.exp(arg) - 1))


# script parameters (scripts/figures/fig_03_1_tail_bounds.R)
P_SCRIPT = (0.6, 0.3, 0.1)
N_SCRIPT = 10_000
C_SCRIPT = (2000, 300, 30)
PSI_SCRIPT = 2
MU_SCRIPT = {"low": (0.2, 0.6, 0.1), "high": (20, 60, 10)}


def run() -> list[dict]:
    out = []
    J = 3
    mu = sp.symbols("mu1:4", positive=True)
    p = sp.symbols("p1:4", positive=True)
    Cs = sp.symbols("C1:4", positive=True)
    n = sp.Symbol("n", positive=True)
    pi_ = sp.Symbol("pi", positive=True)
    m1 = sp.Symbol("m", positive=True)

    # ---- V29a ---------------------------------------------------------------
    z = sp.Symbol("z")
    ser = sp.series((1 - z) ** (-psi), z, 0, 8).removeO()
    ok_ser = all(zero(ser.coeff(z, x) - sp.binomial(x + psi - 1, x).expand(func=True)) for x in range(8))
    M_pmf = ((1 - pi_) / (1 - pi_ * sp.exp(v))) ** psi  # sum_x C(x+psi-1,x) pi^x (1-pi)^psi e^{vx}
    M_mu = M_pmf.subs(pi_, m1 / (m1 + psi))
    # bases agree: (1-pi)/(1-pi e^v) = 1/(1 - m/psi (e^v - 1))
    base_ratio = ((1 - pi_) / (1 - pi_ * sp.exp(v))).subs(pi_, m1 / (m1 + psi)) * (1 - m1 / psi * (sp.exp(v) - 1))
    ok_form = zero(base_ratio - 1)
    # domain: 1 - m/psi (e^v - 1) > 0  <=>  v < log(1 + psi/m) = -log pi
    vb = sp.solve(sp.Eq(1 - m1 / psi * (sp.exp(v) - 1), 0), v)
    ok_dom = len(vb) == 1 and zero(sp.expand_log(vb[0] + sp.log(m1 / (m1 + psi)), force=True))
    # cumulants of one cell from M
    Mm = M_mu
    c1 = sp.diff(Mm, v).subs(v, 0)
    c2 = sp.diff(Mm, v, 2).subs(v, 0)
    c3 = sp.diff(Mm, v, 3).subs(v, 0)
    var_cell = sp.simplify(c2 - c1**2)
    k3_cell = sp.simplify(c3 - 3 * c2 * c1 + 2 * c1**3)
    ok_cell = zero(c1 - m1) and zero(var_cell - (m1 + m1**2 / psi)) and zero(
        k3_cell - (m1 + 3 * m1**2 / psi + 2 * m1**3 / psi**2)
    )
    # K_D with symbolic J = 3 and m_j = p_j n
    KD = sum(p[j] * n * K_nb(u / n, mu[j]) for j in range(J)) + sum(
        Cs[j] * K_nb(-u * p[j] / Cs[j], mu[j]) for j in range(J)
    )
    sig2 = [var_cell.subs(m1, mu[j]) for j in range(J)]
    k3 = [k3_cell.subs(m1, mu[j]) for j in range(J)]
    d0 = KD.subs(u, 0)
    d1 = sp.diff(KD, u).subs(u, 0)
    d2 = sp.diff(KD, u, 2).subs(u, 0)
    d3 = sp.diff(KD, u, 3).subs(u, 0)
    varD = sum(p[j] * sig2[j] / n for j in range(J)) + sum(p[j] ** 2 * sig2[j] / Cs[j] for j in range(J))
    k3D = sum(p[j] * k3[j] / n**2 for j in range(J)) - sum(p[j] ** 3 * k3[j] / Cs[j] ** 2 for j in range(J))
    ok_der = zero(d0) and zero(d1) and zero(d2 - varD) and zero(d3 - k3D)
    # convexity of the per-cell CGF: K'' = mu e^v (1 + mu/psi) / (1 - mu/psi (e^v - 1))^2 > 0
    K2 = sp.diff(K_nb(v, m1), v, 2)
    ok_cvx = zero(K2 - m1 * sp.exp(v) * (1 + m1 / psi) / (1 - m1 / psi * (sp.exp(v) - 1)) ** 2)
    # toy: exact E e^{uD} from convolved pmfs
    mp.mp.dps = 40
    toy = {"n": 4, "m": (3, 1), "C": (2, 3), "mu": (mp.mpf("0.7"), mp.mpf("2.5")), "psi": mp.mpf("1.5")}
    p_toy = [mp.mpf(mm) / toy["n"] for mm in toy["m"]]

    def pmf_cell(muv, xmax):
        pv = muv / (muv + toy["psi"])
        return [
            mp.exp(mp.loggamma(x + toy["psi"]) - mp.loggamma(toy["psi"]) - mp.loggamma(x + 1))
            * pv**x
            * (1 - pv) ** toy["psi"]
            for x in range(xmax)
        ]

    def conv_power(f, k):
        g = [mp.mpf(1)]
        for _ in range(k):
            h = [mp.mpf(0)] * min(len(g) + len(f) - 1, len(f))
            for i, gi in enumerate(g):
                for jj, fj in enumerate(f):
                    if i + jj < len(h):
                        h[i + jj] += gi * fj
            g = h
        return g

    XMAX = 300
    pmfs = [pmf_cell(mm, XMAX) for mm in toy["mu"]]
    # law of the bulk total (m_j cells) and of the reference total (C_j cells) of each type
    totals = {(j, k): conv_power(pmfs[j], k) for j in range(2) for k in (toy["m"][j], toy["C"][j])}
    tails_small = all(sum(d[-50:]) < mp.mpf("1e-40") for d in totals.values())

    def mgf_sum(j, k, a):
        return mp.fsum(pr * mp.e ** (a * x) for x, pr in enumerate(totals[(j, k)]))

    def K_toy(uu):
        val = mp.mpf(0)
        for j in range(2):
            val += toy["m"][j] * -toy["psi"] * mp.log(1 - toy["mu"][j] / toy["psi"] * mp.expm1(uu / toy["n"]))
            val += toy["C"][j] * -toy["psi"] * mp.log(
                1 - toy["mu"][j] / toy["psi"] * mp.expm1(-uu * p_toy[j] / toy["C"][j])
            )
        return val

    errs = []
    for uu in (mp.mpf("-0.8"), mp.mpf("-0.3"), mp.mpf("0.3"), mp.mpf("0.8")):
        logE = mp.fsum(
            [mp.log(mgf_sum(j, toy["m"][j], uu / toy["n"])) for j in range(2)]
            + [mp.log(mgf_sum(j, toy["C"][j], -uu * p_toy[j] / toy["C"][j])) for j in range(2)]
        )
        errs.append(abs(logE - K_toy(uu)))
    ok_toy = tails_small and max(errs) < mp.mpf("1e-25")
    # skewness at the script parameters
    skews = {}
    for name, mus in MU_SCRIPT.items():
        sub = {psi: PSI_SCRIPT, n: N_SCRIPT}
        sub.update({mu[j]: sp.nsimplify(mus[j]) for j in range(J)})
        sub.update({p[j]: sp.nsimplify(P_SCRIPT[j]) for j in range(J)})
        sub.update({Cs[j]: C_SCRIPT[j] for j in range(J)})
        skews[name] = float(k3D.subs(sub) / varD.subs(sub) ** sp.Rational(3, 2))
    ok_skew = round(skews["low"], 2) == -0.08 and round(skews["high"], 2) == -0.06
    ok_a = ok_ser and ok_form and ok_dom and ok_cell and ok_der and ok_cvx and ok_toy and ok_skew
    out.append(
        _rec(
            "V29a",
            r"$K_D(u)=\sum_jm_jK_{gj}(u/n)+\sum_jC_jK_{gj}(-up_j/C_j)$ with "
            r"$K_{gj}(u)=-\psi_g\log\{1-\tilde\mu_{gj}(e^u-1)/\psi_g\}$ for $u<-\log\pi_{gj}$; $K_D$ convex, "
            r"$K_D(0)=K_D'(0)=0$, $K_D''(0)=\operatorname{Var}D$, $K_D'''(0)=\kappa_3(D)$; $D$ left-skewed in "
            r"the simulation",
            "exact",
            "pass" if ok_a else "fail",
            rf"Per-cell MGF from the series of $(1-z)^{{-\psi}}$; printed form and domain recovered; symbolic $J=3$ "
            rf"derivatives match $\operatorname{{Var}}D$ and $\kappa_3(D)$ (per-cell cumulants from $M$); "
            rf"$K_{{gj}}''>0$. Numeric cross-check: toy $J=2$, $E e^{{uD}}$ from convolved pmfs, largest "
            rf"log error {float(max(errs)):.1e} (40 digits); skewness at the script parameters "
            rf"{skews['low']:.3f} and {skews['high']:.3f}",
            "",
            detail_full=(
                f"series {ok_ser}; form {ok_form}; domain {ok_dom}; cell cumulants {ok_cell}; derivatives {ok_der}; "
                f"convex {ok_cvx}; toy {ok_toy} errs {[mp.nstr(e, 5) for e in errs]}; skew {skews}"
            ),
        )
    )

    # ---- V29b ---------------------------------------------------------------
    # endpoints: K_gj'(v) -> +oo as v -> -log pi from below
    eps = sp.Symbol("epsilon", positive=True)
    K1 = sp.diff(K_nb(v, m1), v)
    vplus = sp.log(1 + psi / m1)  # = -log pi
    lim_plus = sp.limit(K1.subs(v, vplus - eps), eps, 0, "+")
    ok_lim = lim_plus == sp.oo
    # bulk term m_j K(u/n): u/n -> -log pi_j  <=> u -> n(-log pi_j); reference term C_j K(-u p_j/C_j):
    # -u p_j/C_j -> -log pi_j <=> u -> C_j log pi_j / p_j, where d/du = -p_j K'(.) -> -oo
    pj, Cj = sp.symbols("p_j C_j", positive=True)
    uref = sp.solve(sp.Eq(-u * pj / Cj, vplus), u)[0]
    ok_uminus = zero(sp.expand_log(uref - Cj * sp.log(m1 / (m1 + psi)) / pj, force=True))
    ubulk = sp.solve(sp.Eq(u / n, vplus), u)[0]
    ok_uplus = zero(sp.expand_log(ubulk - n * (-sp.log(m1 / (m1 + psi))), force=True))
    # concavity of u t - K_D: second derivative -K_D'' < 0 (K_gj'' > 0 from V29a)
    # single sum of n cells: K(u) = n K1(u/n) - m u (centred mean)
    Kss = n * K_nb(u / n, m1) - m1 * u
    a = sp.Symbol("a", positive=True)  # a = m + t or m - t
    root = sp.solve(sp.Eq(sp.diff(Kss, u), a - m1), u)
    rate = n * (a * sp.log(a * (m1 + psi) / (m1 * (a + psi))) + psi * sp.log((m1 + psi) / (a + psi)))
    ok_root = len(root) == 1
    val = (u * (a - m1) - Kss).subs(u, root[0]) if ok_root else sp.nan
    ok_val = ok_root and zero(sp.expand_log(sp.simplify(val - rate), force=True))
    # the root is positive iff a > m (upper tail) and negative iff a < m (lower tail)
    ok_sign = ok_root and zero(sp.expand_log(root[0] - n * sp.log(a * (m1 + psi) / (m1 * (a + psi))), force=True))
    # numeric confirmation on the spec grid
    mp.mp.dps = 40
    worst = mp.mpf(0)
    for mu0 in (mp.mpf("0.2"), mp.mpf(20)):
        for fac in ("0.3", "0.7", "1.4", "3"):
            av = mu0 * mp.mpf(fac)
            nn, ps = mp.mpf(N_SCRIPT), mp.mpf(PSI_SCRIPT)
            Kf = lambda uu: nn * -ps * mp.log(1 - mu0 / ps * mp.expm1(uu / nn)) - mu0 * uu
            dKf = lambda uu: mu0 * mp.e ** (uu / nn) / (1 - mu0 / ps * mp.expm1(uu / nn)) - mu0
            # bisection on (0, u_+) or (-200 n, 0); K' is increasing (convexity)
            lo, hi = (mp.mpf(0), nn * mp.log(1 + ps / mu0) * (1 - mp.mpf("1e-30"))) if av > mu0 else (-200 * nn, mp.mpf(0))
            for _ in range(400):
                mid = (lo + hi) / 2
                if dKf(mid) - (av - mu0) > 0:
                    hi = mid
                else:
                    lo = mid
            ustar = (lo + hi) / 2
            num = ustar * (av - mu0) - Kf(ustar)
            closed = nn * (av * mp.log(av * (mu0 + ps) / (mu0 * (av + ps))) + ps * mp.log((mu0 + ps) / (av + ps)))
            worst = max(worst, abs(num / closed - 1))
    ok_num = worst < mp.mpf("1e-25")
    ok_thr = sp.exp(-sp.log(100)) == sp.Rational(1, 100)
    ok_b = ok_lim and ok_uminus and ok_uplus and ok_val and ok_sign and ok_num and ok_thr
    out.append(
        _rec(
            "V29b",
            r"Legendre transforms $K_D^{\pm}$: supremum at the root of $K_D'(u)=\pm t$ inside "
            r"$(0,u_+)$ or $(u_-,0)$, $u_+=n\min_j(-\log\pi_{gj})$, $u_-=\max_jC_j\log\pi_{gj}/p_j$; for a single sum "
            r"of $n$ cells the rate is $nK^\ast(m\pm t)$ of the NB rate; $K_D^{\pm}=\log100$ gives the $0.01$ bound",
            "exact",
            "pass" if ok_b else "fail",
            rf"$K_{{gj}}'(v)\to\infty$ as $v\uparrow-\log\pi_{{gj}}$; endpoints solved from the arguments $u/n$ and "
            rf"$-up_j/C_j$; single sum: root $u^\ast=n\log\frac{{a(m+\psi)}}{{m(a+\psi)}}$ and value $nK^\ast(a)$ "
            rf"symbolic; mpmath root-finding on the grid $a\in m\{{0.3,0.7,1.4,3\}}$, $m\in\{{0.2,20\}}$, "
            rf"largest relative error {float(worst):.1e}",
            "",
            detail_full=(
                f"limit {ok_lim}; u- {ok_uminus}; u+ {ok_uplus}; value {ok_val}; root sign {ok_sign}; "
                f"numeric {ok_num} ({mp.nstr(worst, 5)}); threshold {ok_thr}"
            ),
        )
    )

    # ---- V29c ---------------------------------------------------------------
    x = sp.Symbol("x", integer=True, nonnegative=True)
    yy = sp.Symbol("y", positive=True)  # y = m e^v: m^x e^{vx} = y^x
    series_y = sp.summation(yy**x / sp.factorial(x), (x, 0, sp.oo))
    mgf_pois = sp.exp(-m1) * series_y.subs(yy, m1 * sp.exp(v))
    K_pois = lambda arg, mm: mm * (sp.exp(arg) - 1)
    ok_pmgf = zero(sp.expand_log(sp.log(mgf_pois), force=True) - K_pois(v, m1))
    KP = sum(p[j] * n * K_pois(u / n, mu[j]) for j in range(J)) + sum(
        Cs[j] * K_pois(-u * p[j] / Cs[j], mu[j]) for j in range(J)
    )
    ok_p0 = zero(KP.subs(u, 0)) and zero(sp.diff(KP, u).subs(u, 0))
    varP = sum(p[j] * mu[j] / n for j in range(J)) + sum(p[j] ** 2 * mu[j] / Cs[j] for j in range(J))
    ok_p2 = zero(sp.diff(KP, u, 2).subs(u, 0) - varP)
    Kps = n * K_pois(u / n, m1) - m1 * u
    rootp = sp.solve(sp.Eq(sp.diff(Kps, u), a - m1), u)
    valp = (u * (a - m1) - Kps).subs(u, rootp[0])
    ok_prate = len(rootp) == 1 and zero(sp.expand_log(sp.simplify(valp - n * (a * sp.log(a / m1) - a + m1)), force=True))
    ratio = sp.simplify(sp.diff(K_nb(v, m1), v, 2).subs(v, 0) / sp.diff(K_pois(v, m1), v, 2).subs(v, 0))
    ok_ratio = zero(ratio - (1 + m1 / psi))
    lim_cgf = sp.limit(K_nb(v, m1), psi, sp.oo)
    ok_lim_cgf = zero(lim_cgf - K_pois(v, m1))
    rate1 = a * sp.log(a * (m1 + psi) / (m1 * (a + psi))) + psi * sp.log((m1 + psi) / (a + psi))
    ok_lim_rate = zero(sp.limit(rate1, psi, sp.oo) - (a * sp.log(a / m1) - a + m1))
    ok_c = ok_pmgf and ok_p0 and ok_p2 and ok_prate and ok_ratio and ok_lim_cgf and ok_lim_rate
    out.append(
        _rec(
            "V29c",
            r"Poisson CGF $\tilde\mu_{gj}(e^u-1)$ and $K^{\mathrm{Pois}}_D$; single-sum rate "
            r"$n\{\bar x\log(\bar x/m)-\bar x+m\}$; the Poisson model understates the per-cell variance by "
            r"$1+\tilde\mu_{gj}/\psi_g$; NB CGF and rate tend to the Poisson ones as $\psi_g\to\infty$",
            "exact",
            "pass" if ok_c else "fail",
            rf"Poisson series summed to ${sp.latex(sp.simplify(mgf_pois))}$; $K^{{\mathrm{{Pois}}}}_D(0)=K'(0)=0$; "
            rf"rate solved at the root of $K'(u)=\bar x-m$; NB/Poisson variance ratio ${sp.latex(ratio)}$; limits "
            rf"$\psi\to\infty$ of the per-cell CGF and of the NB rate symbolic",
            "",
            detail_full=(
                f"mgf {ok_pmgf}; K(0),K'(0) {ok_p0}; K''(0) {ok_p2}; rate {ok_prate}; ratio {ok_ratio}; "
                f"cgf limit {ok_lim_cgf}; rate limit {ok_lim_rate}"
            ),
        )
    )

    # ---- V29d ---------------------------------------------------------------
    # two genes, independent: CGF of Z_1 + Z_2 with Z_g = D_g / sd_g
    s = sp.Symbol("s", real=True)
    ng = 2
    mug = [[sp.Symbol(f"mu_{g}{j}", positive=True) for j in range(1, J + 1)] for g in range(1, ng + 1)]

    def KD_gene(arg, mus, cgf):
        return sum(p[j] * n * cgf(arg / n, mus[j]) for j in range(J)) + sum(
            Cs[j] * cgf(-arg * p[j] / Cs[j], mus[j]) for j in range(J)
        )

    var_nb = [sp.diff(KD_gene(u, mug[g], K_nb), u, 2).subs(u, 0) for g in range(ng)]
    var_po = [sp.diff(KD_gene(u, mug[g], K_pois), u, 2).subs(u, 0) for g in range(ng)]
    sd = [sp.sqrt(vv) for vv in var_nb]
    KS_nb = sum(KD_gene(s / sd[g], mug[g], K_nb) for g in range(ng))
    KS_po = sum(KD_gene(s / sd[g], mug[g], K_pois) for g in range(ng))
    ok_nb2 = zero(sp.diff(KS_nb, s, 2).subs(s, 0) - ng)
    ok_po2 = zero(sp.diff(KS_po, s, 2).subs(s, 0) - sum(var_po[g] / var_nb[g] for g in range(ng)))
    w = [p[j] / n + p[j] ** 2 / Cs[j] for j in range(J)]
    gap = [sp.cancel(var_nb[g] - var_po[g]) for g in range(ng)]
    ok_gap = all(zero(gap[g] - sum(w[j] * mug[g][j] ** 2 / psi for j in range(J))) for g in range(ng))
    ok_d = ok_nb2 and ok_po2 and ok_gap
    out.append(
        _rec(
            "V29d",
            r"Module T3 and T4 add $K_{D_g}(u/\operatorname{sd}D_g)$ over genes as if independent; with NB standard "
            r"deviations the summed NB CGF has variance $\lvert M\rvert$ and the summed Poisson CGF less",
            "exact",
            "pass" if ok_d else "fail",
            r"Two genes, symbolic $J=3$: second derivative at 0 of the summed NB CGF is 2; of the summed Poisson CGF "
            r"$\sum_g\operatorname{Var}^{\mathrm{Pois}}D_g/\operatorname{Var}D_g$, with "
            r"$\operatorname{Var}D_g-\operatorname{Var}^{\mathrm{Pois}}D_g=\sum_j(p_j/n+p_j^2/C_j)\tilde\mu_{gj}^2/\psi>0$",
            "",
            detail_full=f"NB {ok_nb2}; Poisson {ok_po2}; gap {ok_gap}",
        )
    )
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
