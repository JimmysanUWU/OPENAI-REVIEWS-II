"""Checks for 13_PLANAR_ENTRY_POINT.md (sharp L^4 local smoothing in R^{2+1} by a linear Dirichlet problem).

Part A  the Dirac-pencil fibre chart (exact, sympy; C_k is modelled by sympy's I):
        W = P1/(tau + sqrt(tau^2 - P1^2)), Phi = P2 - sqrt(tau^2 - P1^2) invert P = tau B(W) + Phi e_2 with
        B(W) = (2W, 1 - W^2)/(1 + W^2); dPhi = (-(1+W^2), 2W, 1-W^2)/(1-W^2) in (tau, P1, P2);
        e = (1, -2W/(1-W^2)) spans ker d_P Phi; Jcal = e . d_P W = (1+W^2)^2/(2 tau (1-W^2)).
Part L  the disk operator (Lemma 2.2), numerically at random points in 50-digit arithmetic:
        C(B_fib) : H_fib(u) = 4 sum Re(conj(e_a) e_b) d_{X_a X_b} u + 4 |Jcal|^2 Lap_W u for u = U(X, W, Q),
        with B_fib = H_fib(|Phi|^2) = 4 c^* c and C(B_fib) = 4 e e^*; the 3-slot Hessian by central differences.
Part E  the energy algebra (Lemma 3.2), numerically: D_B(H) = beta det(H|_E) for B = beta v v^*;
        with alpha = A|_E >= 0 and alpha_00 = beta b_11: A:H - D_B(H) = beta |b_01 + alpha_01/beta|^2 + det(alpha)/beta >= 0.
Part M  the mixed cofactor (Lemma 3.1): M_TT = C(B_fib) : H_fib, M : H = 2 D_B(H), nu^* M nu = 0 when H|nu^perp = 0.
Part P  polarized cofactor divergence (exact): Q(H(v), H(u)) = C(H(u)+H(v)) - C(H(u)) - C(H(v)) is divergence free
        for random polynomials u, v in three 3-slots; negative control with a non-Hessian field.
Writes ../results/planar_entry.json.
"""

import json
import os
import random
import sys

import mpmath as mp
import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_quaternion_calculus import (cofactor_numeric, from_complex, herm_inner, moore_det, qadd, qconj,  # noqa: E402
                                        qmul, rand_herm, random_poly, sym_cofactor, to_complex)
from verify_slot_calculus import frame_divergence, frame_hessian  # noqa: E402

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "planar_entry.json")
UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]


# ----------------------------------------------------------------- Part A (exact chart identities)

def part_A():
    tau, P1, P2 = sp.symbols("tau P1 P2")
    root = sp.sqrt(tau ** 2 - P1 ** 2)
    W = P1 / (tau + root)
    Phi = P2 - root
    B1 = 2 * W / (1 + W ** 2)
    B2 = (1 - W ** 2) / (1 + W ** 2)
    out = {}
    out["chart_inverts"] = bool(sp.simplify(tau * B1 - P1) == 0 and sp.simplify(tau * B2 + Phi - P2) == 0)
    dPhi = [sp.diff(Phi, v) for v in (tau, P1, P2)]
    target = [-(1 + W ** 2) / (1 - W ** 2), 2 * W / (1 - W ** 2), sp.Integer(1)]
    out["dPhi"] = bool(all(sp.simplify(a - b) == 0 for a, b in zip(dPhi, target)))
    e = [sp.Integer(1), -2 * W / (1 - W ** 2)]
    out["e_in_kernel"] = bool(sp.simplify(e[0] * dPhi[1] + e[1] * dPhi[2]) == 0)
    Jcal = e[0] * sp.diff(W, P1) + e[1] * sp.diff(W, P2)
    out["Jcal"] = bool(sp.simplify(Jcal - (1 + W ** 2) ** 2 / (2 * tau * (1 - W ** 2))) == 0)
    # the characteristic conormal at real W = v(xi): dPhi = (xi.dP - |xi| dtau)/xi_2
    x1, x2 = sp.symbols("x1 x2", positive=True)
    r = sp.sqrt(x1 ** 2 + x2 ** 2)
    v = x1 / (r + x2)
    sub = {W: v}
    conormal = [sp.simplify(sp.radsimp(t.subs(sub))) for t in target]
    out["conormal"] = bool(sp.simplify(conormal[0] + r / x2) == 0 and sp.simplify(conormal[1] - x1 / x2) == 0)
    return out, all(out.values())


# ----------------------------------------------------------------- 3-slot numerics (C_k numbers as mpc)

def kq(c):
    return (mp.re(c), mp.mpf(0), mp.mpf(0), mp.im(c))


def chart_inverse(tau, P1, P2):
    root = mp.sqrt(tau ** 2 - P1 ** 2)
    return P1 / (tau + root), P2 - root


def fib_coords(X, P):
    """real coordinates of the two spatial 3-slots: (X_a, Re P_a, Im P_a)."""
    return [[X[a], mp.re(P[a]), mp.im(P[a])] for a in range(2)]


def hessian_3slot(fn, point, h=mp.mpf("1e-12")):
    """H_ab = sum e_alpha conj(e_beta) d_{a alpha} d_{b beta} fn, by central differences, at point[a][alpha]."""
    H = [[None, None], [None, None]]

    def ev(delta):
        pt = [[point[a][al] + delta.get((a, al), 0) for al in range(3)] for a in range(2)]
        return fn(pt)

    D = {}
    keys = [(a, al) for a in range(2) for al in range(3)]
    f0 = ev({})
    for i, ki in enumerate(keys):
        for j, kj in enumerate(keys):
            if j < i:
                D[(ki, kj)] = D[(kj, ki)]
                continue
            if ki == kj:
                D[(ki, kj)] = (ev({ki: h}) - 2 * f0 + ev({ki: -h})) / h ** 2
            else:
                D[(ki, kj)] = (ev({ki: h, kj: h}) - ev({ki: h, kj: -h}) - ev({ki: -h, kj: h}) + ev({ki: -h, kj: -h})) / (4 * h ** 2)
    for a in range(2):
        for b in range(2):
            acc = (mp.mpf(0),) * 4
            for al in range(3):
                for be in range(3):
                    unit = qmul(UNITS[al], qconj(UNITS[be]))
                    acc = qadd(acc, tuple(D[((a, al), (b, be))] * x for x in unit))
            H[a][b] = acc
    return H


def herm2_to_complex(H):
    A = np.zeros((2, 2, 4))
    for a in range(2):
        for b in range(2):
            A[a, b] = [float(x) for x in H[a][b]]
    return A


def part_L(rng):
    mp.mp.dps = 50
    worst = mp.mpf(0)
    worst_B = mp.mpf(0)
    for trial in range(4):
        tau = mp.mpc(rng.uniform(0.7, 1.3), rng.uniform(-0.3, 0.3))
        a = mp.mpf("0.2")
        W0 = mp.mpc(rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1))
        Q0 = mp.mpc(rng.uniform(-0.05, 0.05), rng.uniform(-0.05, 0.05))
        Phi0 = a * Q0
        P = [tau * 2 * W0 / (1 + W0 ** 2), tau * (1 - W0 ** 2) / (1 + W0 ** 2) + Phi0]
        X = [mp.mpf(rng.uniform(-1, 1)), mp.mpf(rng.uniform(-1, 1))]
        co = [mp.mpf(rng.uniform(-1, 1)) for _ in range(12)]

        def U(X1, X2, w, q):
            wr, wi, qr, qi = mp.re(w), mp.im(w), mp.re(q), mp.im(q)
            return (co[0] * X1 ** 2 + co[1] * X1 * X2 + co[2] * X2 ** 2 + co[3] * mp.sin(X1 + 2 * X2)
                    + co[4] * wr ** 2 + co[5] * wr * wi + co[6] * wi ** 3 + co[7] * X1 * wr + co[8] * X2 * wi
                    + co[9] * qr * wr + co[10] * qi ** 2 + co[11] * X1 * qr)

        def u_of(pt):
            Pc = [mp.mpc(pt[0][1], pt[0][2]), mp.mpc(pt[1][1], pt[1][2])]
            w, ph = chart_inverse(tau, Pc[0], Pc[1])
            return U(pt[0][0], pt[1][0], w, ph / a)

        def phi2(pt):
            Pc = [mp.mpc(pt[0][1], pt[0][2]), mp.mpc(pt[1][1], pt[1][2])]
            return abs(chart_inverse(tau, Pc[0], Pc[1])[1]) ** 2

        pt = fib_coords(X, P)
        Hu = hessian_3slot(u_of, pt)
        Bf = hessian_3slot(phi2, pt)
        W, _ = chart_inverse(tau, P[0], P[1])
        c = [2 * W / (1 - W ** 2), mp.mpc(1)]
        e = [mp.mpc(1), -2 * W / (1 - W ** 2)]
        # B_fib = 4 c^* c, C(B_fib) = 4 e e^* (as quaternionic matrices with C_k entries)
        for i in range(2):
            for j in range(2):
                exp_B = kq(4 * mp.conj(c[i]) * c[j])
                worst_B = max(worst_B, max(abs(x - y) for x, y in zip(Bf[i][j], exp_B)))
        Cmat = [[kq(4 * e[i] * mp.conj(e[j])) for j in range(2)] for i in range(2)]
        lhs = sum(qmul(Cmat[i][j], Hu[j][i])[0] for i in range(2) for j in range(2))
        # right side: 4 sum Re(conj(e_a) e_b) U_{X_a X_b} + 4 |Jcal|^2 Lap_W U
        hh = mp.mpf("1e-12")

        def UX(da, db):
            Xp = list(X)
            def f(dx):
                Z = [X[0] + dx[0], X[1] + dx[1]]
                return U(Z[0], Z[1], W, Phi0 / a)
            if da == db:
                d = [0, 0]; d[da] = hh
                dm = [0, 0]; dm[da] = -hh
                return (f(d) - 2 * f([0, 0]) + f(dm)) / hh ** 2
            dpp = [0, 0]; dpp[da] += hh; dpp[db] += hh
            dpm = [0, 0]; dpm[da] += hh; dpm[db] -= hh
            dmp = [0, 0]; dmp[da] -= hh; dmp[db] += hh
            dmm = [0, 0]; dmm[da] -= hh; dmm[db] -= hh
            return (f(dpp) - f(dpm) - f(dmp) + f(dmm)) / (4 * hh ** 2)

        xpart = 4 * sum(mp.re(mp.conj(e[i]) * e[j]) * UX(i, j) for i in range(2) for j in range(2))
        lapW = sum((U(X[0], X[1], W + s * hh, Phi0 / a) - 2 * U(X[0], X[1], W, Phi0 / a) + U(X[0], X[1], W - s * hh, Phi0 / a)) / hh ** 2
                   for s in (1, mp.mpc(0, 1)))
        Jcal = (1 + W ** 2) ** 2 / (2 * tau * (1 - W ** 2))
        rhs = xpart + 4 * abs(Jcal) ** 2 * lapW
        worst = max(worst, abs(lhs - rhs) / max(abs(rhs), mp.mpf(1)))
    ok = worst < mp.mpf("1e-15") and worst_B < mp.mpf("1e-15")
    return {"disk_operator_relative_error": float(worst), "B_fib_error": float(worst_B)}, bool(ok)


# ----------------------------------------------------------------- Parts E, M (numeric quaternionic linear algebra)

def compress(A, basis):
    """A restricted to the span of a quaternionic orthonormal family (columns), in that basis."""
    m = len(basis)
    out = np.zeros((m, m, 4))
    n = A.shape[0]
    for i in range(m):
        for j in range(m):
            acc = (0.0, 0.0, 0.0, 0.0)
            for a in range(n):
                for b in range(n):
                    acc = qadd(acc, qmul(qmul(qconj(tuple(basis[i][a])), tuple(A[a, b])), tuple(basis[j][b])))
            out[i, j] = acc
    return out


def outer(v, w):
    n = len(v)
    out = np.zeros((n, n, 4))
    for a in range(n):
        for b in range(n):
            out[a, b] = qmul(tuple(v[a]), qconj(tuple(w[b])))
    return out


def unit_complement(v, rng):
    """an orthonormal pair spanning v^perp in H^3 (Gram-Schmidt with right scalars)."""
    from verify_quaternion_calculus import quat_orthonormal_complement
    basis = quat_orthonormal_complement(np.array(v))
    return [basis[:, k] for k in range(1, 3)]


def part_EM(seed):
    rng = np.random.default_rng(seed)
    worst = {"D_B=beta det b": 0.0, "completed_square": 0.0, "M_TT": 0.0, "M:H=2D": 0.0, "nu M nu": 0.0}
    min_gap = np.inf
    for _ in range(6):
        v = rng.normal(size=(3, 4))
        v /= np.sqrt(np.sum(v ** 2))
        beta = rng.uniform(0.5, 2.0)
        B = beta * outer(v, v)
        Ebasis = unit_complement(v, rng)
        H = rand_herm(3, rng)
        b = compress(H, Ebasis)
        detb = b[0, 0, 0] * b[1, 1, 0] - np.sum(b[0, 1] ** 2)
        # D_B(H) = coefficient of lambda in det(lambda B + H) = B : C(H)
        DB = (moore_det(H + B) - moore_det(H - B)) / 2       # det B = 0, so this is the lambda-coefficient
        worst["D_B=beta det b"] = max(worst["D_B=beta det b"], abs(DB - beta * detb) / max(1, abs(DB)))
        # A supported on E with alpha_00 = beta b_11 (need b_11 > 0)
        if b[1, 1, 0] <= 0.05:
            H = H + (0.05 - b[1, 1, 0] + 0.5) * outer(Ebasis[1], Ebasis[1])
            b = compress(H, Ebasis)
            detb = b[0, 0, 0] * b[1, 1, 0] - np.sum(b[0, 1] ** 2)
            DB = (moore_det(H + B) - moore_det(H - B)) / 2
        a01 = rng.normal(size=4)
        a00 = beta * b[1, 1, 0]
        a11 = np.sum(a01 ** 2) / a00 + rng.uniform(0.1, 1.0)
        alpha = np.zeros((2, 2, 4))
        alpha[0, 0, 0], alpha[1, 1, 0] = a00, a11
        alpha[0, 1] = a01
        alpha[1, 0] = qconj(tuple(a01))
        A = np.zeros((3, 3, 4))
        for i in range(2):
            for j in range(2):
                A = A + np.array(outer_scaled(Ebasis[i], alpha[i, j], Ebasis[j]))
        lhs = herm_inner(A, H) - DB
        b01 = b[0, 1]
        sq = np.array(b01) + np.array(a01) / beta
        deta = a00 * a11 - np.sum(a01 ** 2)
        rhs = beta * np.sum(sq ** 2) + deta / beta
        worst["completed_square"] = max(worst["completed_square"], abs(lhs - rhs) / max(1, abs(rhs)))
        min_gap = min(min_gap, lhs)
        # mixed cofactor: M : K = d/dt [B : C(H + tK)]
        def DBf(Hm):
            return (moore_det(Hm + B) - moore_det(Hm - B)) / 2
        eps = 1.0   # D_B is quadratic in H, so central differences with step 1 are exact
        ET = np.zeros((3, 3, 4)); ET[0, 0, 0] = 1
        MTT = (DBf(H + eps * ET) - DBf(H - eps * ET)) / (2 * eps)
        Bf = B[1:, 1:]
        Hf = H[1:, 1:]
        CBf = cofactor_numeric(Bf)
        worst["M_TT"] = max(worst["M_TT"], abs(MTT - herm_inner(CBf, Hf)) / max(1, abs(MTT)))
        MH = (DBf(H + eps * H) - DBf(H - eps * H)) / (2 * eps)
        worst["M:H=2D"] = max(worst["M:H=2D"], abs(MH - 2 * DB) / max(1, abs(DB)))
        # nu^* M nu = 0 when H|nu^perp = 0: take H = nu nu^* x + (nu w^* + w nu^*), w random
        nu = rng.normal(size=(3, 4)); nu /= np.sqrt(np.sum(nu ** 2))
        w = rng.normal(size=(3, 4))
        Hn = 1.7 * outer(nu, nu) + outer(nu, w) + outer(w, nu)
        perp = unit_complement(nu, rng)
        Hn_perp_norm = np.max(np.abs(compress(Hn, perp)))
        Nn = outer(nu, nu)
        nuMnu = (DBf(Hn + eps * Nn) - DBf(Hn - eps * Nn)) / (2 * eps)
        worst["nu M nu"] = max(worst["nu M nu"], abs(nuMnu) + 0 * Hn_perp_norm)
    ok = all(v_ < 1e-6 for v_ in worst.values()) and min_gap >= -1e-9
    worst["min_A:H-D_B"] = float(min_gap)
    return worst, bool(ok)


def outer_scaled(x, q, y):
    """x q y^* for quaternionic columns x, y and a quaternion q."""
    n = len(x)
    out = np.zeros((n, n, 4))
    for a in range(n):
        for b in range(n):
            out[a, b] = qmul(qmul(tuple(x[a]), tuple(q)), qconj(tuple(y[b])))
    return out


# ----------------------------------------------------------------- Part P (polarized cofactor divergence, exact)

def part_P(rng):
    frame = [UNITS[:3]] * 3
    ys = [[sp.Symbol(f"y{a}_{i}", real=True) for i in range(3)] for a in range(3)]
    flat = [y for row in ys for y in row]
    u = random_poly(flat, 4, rng)
    v = random_poly(flat, 3, rng)
    Hu = frame_hessian(u, ys, frame)
    Hv = frame_hessian(v, ys, frame)
    Hs = [[tuple(sp.expand(x + y) for x, y in zip(Hu[a][b], Hv[a][b])) for b in range(3)] for a in range(3)]
    Cs, Cu, Cv = sym_cofactor(Hs, 3), sym_cofactor(Hu, 3), sym_cofactor(Hv, 3)
    Qm = [[tuple(sp.expand(x - y - z) for x, y, z in zip(Cs[a][b], Cu[a][b], Cv[a][b])) for b in range(3)] for a in range(3)]
    div = frame_divergence(Qm, ys, frame)
    # control: replace H(v) by a Hermitian field that is not a Hessian
    Hbad = [row[:] for row in Hv]
    Hbad[1][1] = (sp.expand(Hbad[1][1][0] + flat[0] ** 2 * flat[4]), 0, 0, 0)
    Hs2 = [[tuple(sp.expand(x + y) for x, y in zip(Hu[a][b], Hbad[a][b])) for b in range(3)] for a in range(3)]
    Cs2, Cb = sym_cofactor(Hs2, 3), sym_cofactor(Hbad, 3)
    Qbad = [[tuple(sp.expand(x - y - z) for x, y, z in zip(Cs2[a][b], Cu[a][b], Cb[a][b])) for b in range(3)] for a in range(3)]
    div_bad = frame_divergence(Qbad, ys, frame)
    return {"mixed_cofactor_divergence_nonzero_components": div, "control_non_hessian": div_bad}, bool(div == 0 and div_bad > 0)


def main():
    rng = random.Random(2002)
    results = {}
    all_ok = True
    for name, fn in (("A", part_A), ("L", lambda: part_L(rng)), ("EM", lambda: part_EM(5)), ("P", lambda: part_P(rng))):
        res, ok = fn()
        results[name] = res
        all_ok = all_ok and ok
        print(name, "ok" if ok else "FAIL", json.dumps(res, default=str))
    results["all_ok"] = bool(all_ok)
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print("ALL OK" if all_ok else "FAILED")


if __name__ == "__main__":
    main()
