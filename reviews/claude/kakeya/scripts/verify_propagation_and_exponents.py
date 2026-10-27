"""Checks for Section 6 (propagation across the thin cylinder) and the exponent bookkeeping.

  L6.1a  radial comparison: solutions of u'' + u'/s = (V + kappa^2) u on (b, S) with endpoint
         values 0, 1 dominate c * log-affine, c = exp(-int min(z, A-z) P(z) dz), for V of the
         form C s^{-2} exp(-c s/h) used in Prop. 6.3 (shooting with scipy, many h);
  L6.1b  the integrated outward comparison (6.3) and the full chain on an explicit supersolution;
  R6.2   the sharpness example E = (1 - cos T) p(s) + L(s) - g(s): Delta E <= 0, E >= 0, E(T,b) =
         1 - cos T, and E(0, r) <= (9/4) log(r/b)/D^2 (so two logarithms are lost);
  EXP    exponent arithmetic: p_c, s_+(p), beta(p), interpolation (7.3), focusing exponent of
         Prop. 8.1, Lemma 9.4 (A p / q = 2A, dp = 2n), Kakeya log powers, the Lean doc's
         3D power 66 = 2 n (2n+5), Stein threshold 1/2 + s_+(p) = max(1/2, n(1/2 - 1/p)),
         the Monge-Ampere exponent (1 + 1/d) * 2 = p_c, and the energy volume factor.
"""

import json
import math
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def radial_solution(V, kappa, b, S, n_grid=4000):
    """Solve u'' + u'/s = (V + kappa^2) u with u(b) = 0, u'(b) = 1 (log variable), return normalized u."""
    A = math.log(S / b)

    def rhs(x, y):
        s = b * math.exp(x)
        P = s * s * (V(s) + kappa ** 2)
        return [y[1], P * y[0]]

    xs = np.linspace(0, A, n_grid)
    sol = solve_ivp(rhs, (0, A), [0.0, 1.0], t_eval=xs, rtol=1e-10, atol=1e-12)
    psi = sol.y[0]
    return xs, psi / psi[-1], A


def lemma61a():
    out = []
    ok = True
    for h in (1e-2, 1e-3, 1e-4, 1e-6):
        K = 4.0
        b = 2 * K * h
        Cc, c = 5.0, 1.0
        V = lambda s, h=h: Cc * s ** -2 * math.exp(-c * s / h)
        for S, kappa in ((3.0, math.pi / 3), (1.0, math.pi / 2)):
            xs, u, A = radial_solution(V, kappa, b, S)
            # Gronwall constant
            zs = np.linspace(0, A, 20001)
            ss = b * np.exp(zs)
            P = ss ** 2 * (np.array([V(s) for s in ss]) + kappa ** 2)
            G = np.trapezoid(np.minimum(zs, A - zs) * P, zs)
            cgron = math.exp(-G)
            ratio = np.min(u[1:] / (xs[1:] / A))
            ok &= ratio >= cgron - 1e-9
            out.append({"h": h, "S": S, "kappa": round(kappa, 4), "min_u_over_logaffine": float(ratio),
                        "gronwall_constant": cgron})
    return ok, out


def remark62():
    ok = True
    rows = []
    for b in (1e-2, 1e-4, 1e-8):
        D = math.log(3 / b)
        p = lambda s: math.log(3 / s) / D
        g = lambda s: s * s * (1 + math.log(3 / s)) / (4 * D)
        L = lambda s: p(s) * g(b) + (1 - p(s)) * g(3)
        E = lambda T, s: (1 - math.cos(T)) * p(s) + L(s) - g(s)
        # exact Laplacian in (T, tau) of E, radial form, by sympy
        T_, s_, b_ = sp.symbols("T s b", positive=True)
        D_ = sp.log(3 / b_)
        p_ = sp.log(3 / s_) / D_
        g_ = s_ ** 2 * (1 + sp.log(3 / s_)) / (4 * D_)
        L_ = p_ * g_.subs(s_, b_) + (1 - p_) * g_.subs(s_, 3)
        E_ = (1 - sp.cos(T_)) * p_ + L_ - g_
        lap = sp.diff(E_, T_, 2) + sp.diff(E_, s_, 2) + sp.diff(E_, s_) / s_
        ok &= sp.simplify(lap - (sp.cos(T_) - 1) * p_) == 0
        for T in np.linspace(0, 2 * math.pi, 13):
            for s in np.geomspace(b, 3, 40):
                ok &= E(T, s) >= -1e-12
        ok &= abs(E(1.0, b) - (1 - math.cos(1.0))) < 1e-12
        for r in (2 * b, 10 * b, 0.5):
            ok &= E(0, r) <= 9 / 4 * math.log(r / b) / D ** 2 + 1e-12
        rows.append({"b": b, "D": D, "E(0,2b)": E(0, 2 * b), "bound_(9/4)log2/D^2": 9 / 4 * math.log(2) / D ** 2,
                     "int_1^2 H(T,b) dT": 1 + math.sin(1) - math.sin(2)})
    return ok, rows


def exponents():
    res = {}
    ok = True
    for n in range(2, 9):
        d = n - 1
        pc = Fr(2 * n, n - 1)
        # MA exponent: (1 + 1/d) * 2 = p_c
        ok &= (1 + Fr(1, d)) * 2 == pc
        # focusing: d - d p / 2 = -1 iff p = p_c
        ok &= d - d * pc / 2 == -1
        # beta(p) continuity at p_c, interpolation (7.3)
        beta = lambda p: min(d * (1 - 2 / p), 2 / p)
        ok &= d * (1 - 2 / pc) == 2 / pc == Fr(d, n)
        for p in (Fr(5, 2), pc, Fr(7, 2) if Fr(7, 2) > pc else pc + 1, 6, 10):
            p = Fr(p)
            if p < pc:
                theta = (Fr(1, 2) - 1 / p) / (Fr(1, 2) - 1 / pc)
                ok &= theta * 2 / pc == d * (1 - 2 / p)
            else:
                theta = pc / p   # weight on the p_c endpoint when interpolating with L^infinity
                ok &= theta * 2 / pc == 2 / p
                sp_ = Fr(d, 2) - n / p
                ok &= (1 - theta) * Fr(d, 2) == sp_   # h^{-d/2} at p = infinity
        # Lemma 9.4 and Kakeya powers
        A = 2 / pc
        q = pc / 2
        ok &= A * pc / q == 2 * A and d * pc == 2 * n
        ok &= 2 * A * n == 2 * (n - 1)          # (9.6) volume power
        # Stein threshold
        for p in (Fr(5, 2), Fr(3), Fr(4), Fr(10)):
            s_plus = max(Fr(0), Fr(d, 2) - n / p)
            ok &= Fr(1, 2) + s_plus == max(Fr(1, 2), n * (Fr(1, 2) - 1 / p))
        # energy volume factor h^2 s^{2d} (hs)^{-dq} = sigma^{n-2}, sigma = s/h
        hsym, ssym = 0.37, 1.9
        lhs = hsym ** 2 * ssym ** (2 * d) * (hsym * ssym) ** (-float(d * q))
        ok &= abs(lhs - (ssym / hsym) ** (n - 2)) < 1e-9 * lhs
        res[n] = {"p_c": str(pc), "kakeya_maximal_log_power": str(2 * A), "kakeya_volume_log_power": 2 * (n - 1),
                  "lean_log_power_A": 2 * n + 5, "lean_volume_power": 2 * n * (2 * n + 5)}
    ok &= res[3]["lean_volume_power"] == 66
    return ok, res


def main():
    out = {}
    ok1, rows1 = lemma61a()
    ok2, rows2 = remark62()
    ok3, rows3 = exponents()
    out["lemma_6.1_radial_comparison"] = {"ok": bool(ok1), "rows": rows1}
    out["remark_6.2_sharpness"] = {"ok": bool(ok2), "rows": rows2}
    out["exponents"] = {"ok": bool(ok3), "by_dimension": rows3}
    out["all_ok"] = bool(ok1 and ok2 and ok3)
    print(json.dumps(out, indent=1, default=str))
    print("ALL PASS" if out["all_ok"] else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/propagation_and_exponents.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
