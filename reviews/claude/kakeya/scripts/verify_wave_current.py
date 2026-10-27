"""Symbolic/numeric checks of Section 2 (the conserved current) for n = 2 and n = 3.

  S2.1  null coordinates: the plane-wave phase X.xi - T|xi| equals kappa(s0 + v.x - |v|^2 t/2);
  S2.2  the symbol identity (2.2): -C(p) = kappa [[|v|^2, -v^t], [-v, I]] and
        -w^t C(p) w = kappa |w' - v w0|^2; P f_sys = 0 for f_sys = sum kappa^{d/2} E(v) u;
  S2.4  p.q = -(kappa/2)|zeta - v R|^2 for z = w1 + k w2 = (R, zeta)    (identity (2.4));
  L2.1  B = Re sum_l Z_l (x) conj(Z_l) is PSD and every column has zero divergence in
        (X, w1, w2), evaluated at random points for a random finite superposition;
        negative control: replacing conj(f(X - iq)) by conj(f(X + iq)) breaks conservation;
  S2.3  the holomorphic flux map (2.9): Theta(R, zeta, eta) = (R^2, R zeta, eta - zeta.zeta/2),
        Jacobian 2R^{d+1}, M E(v) = 2 R V^v, (dPhi) V^v = 0, (dW) V^v = (v - W)/(2Y).
"""

import json
import random

import numpy as np
import sympy as sp


def null_phase_check(rng):
    ok = True
    for _ in range(20):
        d = rng.choice([1, 2])
        xi = [rng.uniform(-2, 2) for _ in range(d)] + [rng.uniform(0.1, 2)]
        x = [rng.uniform(-3, 3) for _ in range(d)]
        b, T = rng.uniform(-3, 3), rng.uniform(-3, 3)
        nrm = float(np.linalg.norm(xi))
        phys = sum(a * c for a, c in zip(x, xi[:d])) + b * xi[d] - T * nrm
        t, s0 = (T + b) / np.sqrt(2), (b - T) / np.sqrt(2)
        kappa = (xi[d] + nrm) / np.sqrt(2)
        v = [c / kappa for c in xi[:d]]
        null = kappa * (s0 + sum(a * c for a, c in zip(v, x)) - 0.5 * sum(c * c for c in v) * t)
        ok &= abs(phys - null) < 1e-10
    return ok


def matrices(d):
    n = d + 1
    E = lambda a, b: sp.Matrix(n, n, lambda i, j: 1 if (i, j) == (a, b) else 0)
    Ct = 2 * E(0, 0)
    Cx = [E(0, j) + E(j, 0) for j in range(1, d + 1)]
    Cs = -sum((E(j, j) for j in range(1, d + 1)), sp.zeros(n, n))
    return [Ct] + Cx + [Cs]   # coordinates (t, x_1..x_d, s0)


def symbol_check(d):
    kap = sp.symbols("kappa", positive=True)
    v = sp.symbols(f"v1:{d + 1}", real=True)
    w = sp.symbols(f"w0:{d + 1}", real=True)
    C = matrices(d)
    p = [kap * (-sum(c ** 2 for c in v) / 2)] + [kap * c for c in v] + [kap]
    Cp = sum((pi * Ci for pi, Ci in zip(p, C)), sp.zeros(d + 1, d + 1))
    vv = sp.Matrix(v)
    target = kap * sp.Matrix(sp.BlockMatrix([[sp.Matrix([[vv.dot(vv)]]), -vv.T], [-vv, sp.eye(d)]]))
    ok1 = sp.simplify(-Cp - target) == sp.zeros(d + 1, d + 1)
    W = sp.Matrix(w)
    quad = sp.expand(-(W.T * Cp * W)[0])
    ok2 = sp.simplify(quad - kap * sum((w[j + 1] - v[j] * w[0]) ** 2 for j in range(d))) == 0
    # P f = 0 for one mode
    X = sp.symbols(f"t x1:{d + 1} s0", real=True)
    u = sp.exp(sp.I * kap * (X[-1] + sum(v[j] * X[1 + j] for j in range(d)) - X[0] * vv.dot(vv) / 2))
    f = sp.Matrix([u] + [v[j] * u for j in range(d)])
    Pf = sum((C[m] * f.diff(X[m]) for m in range(d + 2)), sp.zeros(d + 1, 1))
    ok3 = all(sp.simplify(e) == 0 for e in Pf)
    return ok1, ok2, ok3


def gaussian_check(d, rng):
    ok = True
    for _ in range(20):
        kap = rng.uniform(0.2, 3)
        v = np.array([rng.uniform(-1, 1) for _ in range(d)])
        w1 = np.array([rng.uniform(-1, 1) for _ in range(d + 1)])
        w2 = np.array([rng.uniform(-1, 1) for _ in range(d + 1)])
        p = np.concatenate([[kap * (-v @ v / 2)], kap * v, [kap]])
        C = [np.array(M.tolist(), dtype=float) for M in matrices(d)]
        q = np.array([0.5 * (w1 @ Ci @ w1 + w2 @ Ci @ w2) for Ci in C])
        R = complex(w1[0], w2[0])  # C_k identified with C for the modulus
        zeta = w1[1:] + 1j * w2[1:]
        rhs = -kap / 2 * np.sum(np.abs(zeta - v * R) ** 2)
        ok &= abs(p @ q - rhs) < 1e-10
    return ok


def current_check(d, seed, n_modes=3, n_points=4, flip=False):
    """Lemma 2.1: build B symbolically for a random superposition and test PSD + divergence."""
    rng = random.Random(seed)
    n = d + 1
    X = list(sp.symbols(f"t x1:{d + 1} s0", real=True))
    W1 = list(sp.symbols(f"a0:{n}", real=True))
    W2 = list(sp.symbols(f"b0:{n}", real=True))
    coords = X + W1 + W2
    C = matrices(d)
    q = [sp.Rational(1, 2) * ((sp.Matrix([W1]) * Cm * sp.Matrix(W1))[0] + (sp.Matrix([W2]) * Cm * sp.Matrix(W2))[0])
         for Cm in C]
    sgn = 1 if flip else -1
    Xc = [X[m] + sgn * sp.I * q[m] for m in range(len(X))]   # X - iq
    modes = []
    for _ in range(n_modes):
        kap = sp.Rational(rng.randint(2, 6), 4)
        v = [sp.Rational(rng.randint(-3, 3), 8) for _ in range(d)]
        c = sp.Rational(rng.randint(-3, 3), 2) + sp.I * sp.Rational(rng.randint(-3, 3), 2)
        modes.append((kap, v, c))
    f = [0] * n
    for kap, v, c in modes:
        phase = sp.I * kap * (Xc[-1] + sum(v[j] * Xc[1 + j] for j in range(d)) - Xc[0] * sum(x ** 2 for x in v) / 2)
        u = c * kap ** sp.Rational(d, 2) * sp.exp(phase)
        Ev = [1] + v
        for a in range(n):
            f[a] += Ev[a] * u
    G = [sp.conjugate(e) for e in f]
    # Z_l components on coords
    Z = []
    for ell, Wl in enumerate((W1, W2)):
        comp = [0] * len(coords)
        Cw = [Cm * sp.Matrix(Wl) for Cm in C]
        for m in range(len(X)):
            comp[m] = sp.I * sum(G[a] * Cw[m][a] for a in range(n))
        off = len(X) + ell * n
        for a in range(n):
            comp[off + a] = G[a]
        Z.append(comp)
    N = len(coords)
    B = [[sp.re(sp.expand_complex(0)) for _ in range(N)] for _ in range(N)]
    Bexpr = [[sum(Z[l][I] * sp.conjugate(Z[l][J]) for l in range(2)) for J in range(N)] for I in range(N)]
    # divergence of columns: sum_I d_I Re(B^{IJ})
    div = [sum(sp.diff(Bexpr[I][J], coords[I]) for I in range(N)) for J in range(N)]
    f_div = sp.lambdify(coords, [sp.re(e) for e in div], "numpy")
    f_B = sp.lambdify(coords, [[sp.re(Bexpr[I][J]) for J in range(N)] for I in range(N)], "numpy")
    prng = np.random.default_rng(seed)
    max_div, min_eig, scale = 0.0, 1e9, 0.0
    for _ in range(n_points):
        pt = prng.uniform(-0.6, 0.6, size=N)
        dv = np.array(f_div(*pt), dtype=complex)
        Bm = np.array(f_B(*pt), dtype=complex).real
        scale = max(scale, np.abs(Bm).max())
        max_div = max(max_div, np.abs(dv).max())
        min_eig = min(min_eig, np.linalg.eigvalsh((Bm + Bm.T) / 2).min())
    return {"d": d, "seed": seed, "flip": flip, "max_abs_divergence": float(max_div),
            "max_abs_B": float(scale), "min_eigenvalue_B": float(min_eig)}


def flux_map_check(d):
    R = sp.symbols("R")
    zeta = sp.Matrix(sp.symbols(f"zeta1:{d + 1}"))
    eta = sp.symbols("eta")
    v = sp.Matrix(sp.symbols(f"v1:{d + 1}"))
    z = sp.Matrix([R] + list(zeta))
    C = matrices(d)
    c = [0] * (d + 1) + [1]
    omega = [sp.expand(sp.Rational(1, 2) * (z.T * Cm * z)[0] + c[m] * eta) for m, Cm in enumerate(C)]
    expect = [R ** 2] + [R * zeta[j] for j in range(d)] + [eta - (zeta.T * zeta)[0] / 2]
    ok_map = all(sp.simplify(a - b) == 0 for a, b in zip(omega, expect))
    Jac = sp.Matrix(omega).jacobian([R] + list(zeta) + [eta])
    ok_jac = sp.simplify(Jac.det() - 2 * R ** (d + 1)) == 0
    # M E(v) with M_{mu a} = (C^mu z)_a
    M = sp.Matrix([[(Cm * z)[a] for a in range(d + 1)] for Cm in C])
    Ev = sp.Matrix([1] + list(v))
    W = zeta / R
    Vv = sp.Matrix([1] + [(v[j] + W[j]) / 2 for j in range(d)] + [-(W.T * v)[0] / 2])
    ok_ME = all(sp.simplify(e) == 0 for e in (M * Ev - 2 * R * Vv))
    # target coordinates (Y, P, S1) = omega; Phi = S1 + P.P/(2Y), W = P/Y
    Y, S1 = sp.symbols("Y S1")
    P = sp.Matrix(sp.symbols(f"P1:{d + 1}"))
    Phi = S1 + (P.T * P)[0] / (2 * Y)
    Wt = P / Y
    vars_ = [Y] + list(P) + [S1]
    subs = {Y: R ** 2, S1: eta - (zeta.T * zeta)[0] / 2}
    subs.update({P[j]: R * zeta[j] for j in range(d)})
    dPhi = sp.Matrix([[sp.diff(Phi, x) for x in vars_]]).subs(subs)
    ok_phi = sp.simplify((dPhi * Vv)[0]) == 0 and sp.simplify(Phi.subs(subs) - eta) == 0
    ok_W = True
    for j in range(d):
        dW = sp.Matrix([[sp.diff(Wt[j], x) for x in vars_]]).subs(subs)
        ok_W &= sp.simplify((dW * Vv)[0] - (v[j] - W[j]) / (2 * R ** 2)) == 0
    return {"d": d, "theta_formula": ok_map, "jacobian_2R^(d+1)": bool(ok_jac), "ME=2RV": ok_ME,
            "dPhi.V=0_and_Phi=eta": bool(ok_phi), "dW.V=(v-W)/2Y": bool(ok_W)}


def main():
    rng = random.Random(7)
    out = {"null_phase": null_phase_check(rng)}
    out["symbol"] = {d: dict(zip(["-C(p)_matrix", "quadratic_form", "P_f_sys=0"], symbol_check(d))) for d in (1, 2)}
    out["gaussian_(2.4)"] = {d: gaussian_check(d, rng) for d in (1, 2)}
    out["flux_map"] = [flux_map_check(d) for d in (1, 2)]
    cur = [current_check(1, s) for s in (1, 2)] + [current_check(2, 3, n_modes=2, n_points=3)]
    neg = [current_check(1, 4, flip=True)]
    out["lemma_2.1"] = cur
    out["lemma_2.1_negative_control"] = neg
    ok = out["null_phase"] and all(all(v.values()) for v in out["symbol"].values())
    ok &= all(out["gaussian_(2.4)"].values()) and all(all(v for k, v in r.items() if k != "d") for r in out["flux_map"])
    ok &= all(r["max_abs_divergence"] < 1e-9 * max(1, r["max_abs_B"]) and r["min_eigenvalue_B"] > -1e-9 for r in cur)
    ok &= all(r["max_abs_divergence"] > 1e-3 for r in neg)
    out["all_ok"] = bool(ok)
    print(json.dumps(out, indent=1, default=str))
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/wave_current.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
