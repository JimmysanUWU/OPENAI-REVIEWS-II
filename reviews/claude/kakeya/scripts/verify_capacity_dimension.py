"""The capacity-dimension principle (TOOLKIT Proposition K2): why the imaginary displacement must be
two-dimensional.

Parameter space R_T x R^k_tau, radial in s = |tau|.  For each k we build the analogue of Remark 6.2:

    p_k(s): radial harmonic in R^k with p(b) = 1, p(3) = 0
    g_k(s): radial, Delta_k g = p, g(b) = g'(b) = 0           (so g is increasing and subharmonic)
    L_k(s): radial harmonic with L = g at s = b and s = 3       (L = g(3)(1 - p))
    E_k(T, s) = (1 - cos T) p(s) + L(s) - g(s)

Checks (exact, sympy): Delta_{T,tau} E = (cos T - 1) p <= 0 (superharmonic), E(T, b) = 1 - cos T,
and E >= 0 on a grid.  Hence int_1^2 E(T,b) dT = 1 + sin 1 - sin 2 is a fixed constant while E(0, 2b)
is computed exactly.  Output: E(0,2b) for b = 1e-2 ... 1e-10 and fitted decay exponents.

Expected: k = 2 decays only like 1/log^2(1/b) (Lemma 6.1 shows this is attained); every k != 2 decays
like a power of b, so any scalar supersolution argument in a parameter space with k != 2 loses N^{-c}.
"""

import json
import math

import sympy as sp

s, T, b, u = sp.symbols("s T b u", positive=True)


def build(k):
    if k == 2:
        p = sp.log(3 / s) / sp.log(3 / b)
    else:
        p = (s ** (2 - k) - sp.Integer(3) ** (2 - k)) / (b ** (2 - k) - sp.Integer(3) ** (2 - k))
    inner = sp.integrate((u ** (k - 1) * p.subs(s, u)), (u, b, s))
    gprime = sp.simplify(inner / s ** (k - 1))
    g = sp.integrate(gprime.subs(s, u), (u, b, s))
    L = g.subs(s, 3) * (1 - p)
    E = (1 - sp.cos(T)) * p + L - g
    return p, g, L, E


def radial_laplacian(f, k):
    return sp.diff(f, T, 2) + sp.diff(f, s, 2) + (k - 1) / s * sp.diff(f, s)


def main():
    out = {}
    ok = True
    for k in (1, 2, 3, 4, 5):
        p, g, L, E = build(k)
        lap = sp.simplify(radial_laplacian(E, k) - (sp.cos(T) - 1) * p)
        trace = sp.simplify(E.subs(s, b) - (1 - sp.cos(T)))
        lap_ok = lap == 0
        trace_ok = trace == 0
        fE = sp.lambdify((T, s, b), E, "mpmath")
        nonneg = True
        rows = []
        for bv in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10):
            import mpmath as mp
            mp.mp.dps = 50
            bm = mp.mpf(bv)
            for Tv in (0, 0.5, 1.5, 3.0):
                for sv in (bm, 2 * bm, 10 * bm, mp.mpf("0.1"), mp.mpf(1), mp.mpf(2.9)):
                    nonneg &= fE(Tv, sv, bm) >= -mp.mpf(10) ** -40
            val = fE(0, 2 * bm, bm)
            rows.append({"b": bv, "E(0,2b)": float(val)})
        # fitted exponent between b=1e-6 and b=1e-10, and log-decay check for k = 2
        e1, e2 = rows[2]["E(0,2b)"], rows[4]["E(0,2b)"]
        exponent = math.log(e1 / e2) / math.log(1e-6 / 1e-10)
        logfit = [r["E(0,2b)"] * math.log(3 / r["b"]) ** 2 for r in rows]
        out[k] = {"laplacian_identity": lap_ok, "trace_identity": trace_ok, "nonnegative_on_grid": bool(nonneg),
                  "rows": rows, "power_exponent_fit": exponent,
                  "E*log^2(3/b)": logfit}
        ok &= lap_ok and trace_ok and nonneg
        print(k, lap_ok, trace_ok, nonneg, "exponent~%.3f" % exponent, [round(x, 4) for x in logfit])
    ok &= abs(out[2]["power_exponent_fit"]) < 0.2 and all(out[k]["power_exponent_fit"] > 0.9 for k in (1, 3, 4, 5))
    out["all_ok"] = bool(ok)
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/capacity_dimension.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
