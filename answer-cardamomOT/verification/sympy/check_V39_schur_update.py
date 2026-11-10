# /// script
# requires-python = ">=3.10"
# dependencies = ["sympy==1.14.0"]
# ///
"""V39: removing gene u from the panel updates the precision by a Schur
complement (vignette 03-2, eq-schur-update).

Derivation route: a generic symmetric 3x3 precision Omega over genes
(g, g', u); Sigma = Omega^{-1}; the panel precision is the inverse of the
{g, g'} block of Sigma, compared with
omega_{gg'} - omega_{gu} omega_{g'u} / omega_{uu}.
"""

from __future__ import annotations

import json

import sympy as sp


def run():
    w = sp.symbols("w11 w12 w13 w22 w23 w33", real=True)
    Om = sp.Matrix([[w[0], w[1], w[2]], [w[1], w[3], w[4]], [w[2], w[4], w[5]]])
    panel = Om.inv()[:2, :2].inv()
    claim = Om[:2, :2] - Om[:2, 2] * Om[2, :2] / Om[2, 2]
    ok = all(sp.simplify(e) == 0 for e in panel - claim)
    # gene g not adjacent to u (omega_gu = 0): row g of the panel precision is unchanged
    unchanged = sp.simplify((claim - Om[:2, :2])[0, :].subs(w[2], 0)) == sp.zeros(1, 2)
    return [
        {
            "id": "V39",
            "claim": r"Panel precision after removing gene $u$: $\omega'_{gg'}=\omega_{gg'}-\omega_{gu}\omega_{g'u}/\omega_{uu}$; entries of a non-neighbour of $u$ are unchanged",
            "kind": "exact",
            "status": "pass" if ok and unchanged else "fail",
            "detail": r"generic symmetric $3\times3$ precision; residual $\bm 0$; $\omega_{gu}=0$ leaves row $g$ unchanged",
            "tex_line": "",
        }
    ]


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
