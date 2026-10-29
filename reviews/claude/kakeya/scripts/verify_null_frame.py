"""Checks for 16_NULL_FRAME_AND_SCHRODINGER.md: the null frame, the Husimi current, and a direct engine for
the Schrodinger flow (roadmap target T.KAK.12).

The current is the paper's (08 Prop. 3.3), kept in the null slots (t, x_1..x_d, s0) with auxiliary parts
(Y, P, S1), *without* rotating to the physical time T.  The parameter slot is t + Y i.

  G  the null chart, exactly (sympy, d = 1, 2, 3): the inverse chart P = YW, S1 = aQ - Y W.W/2; its Jacobian
     a Y^d; dPhi = l(W); Phi o Theta = eta; (dW_j) V^v = (v_j - W_j)/(2Y); tangency l(W) V^v = 0; and the
     imaginary displacement y (1, W, -W.W/2) is a null vector of 2 dt ds0 + |dx|^2.
  H  the time row is the analytically continued wave (50 digits, d = 1, 2, 3):
         |chi_t|^2 = y^{-d} |f0(t - i y, x - i y Re W, s0 + i y |W|^2/2)|^2,   f0 = sum kappa^{d/2} u_{kappa,v};
     for one kappa, |chi_t|^2 = (kappa/y)^d e^{-kappa y |W|^2} |psi(t - i y, x - i y Re W)|^2 (psi the
     Schrodinger solution).  Negative control: the opposite imaginary shift in x fails.
  A  the engine rows: |chi_t|^2 = |sum m u|^2 (one real Gaussian multiplier, so (A4) needs no splitting);
     |chi_mu|^2 = |sum Re_k(V_mu) m u|^2 + |sum Im_k(V_mu) m u|^2 for every row (A3, L = 2).
  C  conservation on one-kappa (Schrodinger) data, d = 1, 2, 3: the real divergence of (chi chi^*)^R vanishes
     (central differences, 60 digits) and the tensor does not depend on s0.  Controls: two different kappas make
     the tensor depend on s0; dropping the rotation iota_e breaks conservation.
  S  the null reduction, exactly: the null system on e^{i kappa s0} F is the Schrodinger equation
     d_t psi = (i/2kappa) Laplace psi; the Galilean identity; the exponent bookkeeping of Theorem 16.8.
  W  why a lift is needed (exact): for a Gaussian G = exp(-|a|^2/alpha - |b|^2/beta), b.grad_a G = -(beta/2) div_a grad_b G;
     hence the Husimi density H = |Gabor transform|^2 (window width w, beta = 1/w^2) of a Schrodinger solution satisfies
     d_t H + d_x( (xi/kappa) H + (1/(2 kappa w^2)) d_xi H ) = 0, whose flux has no sign.  Checked symbolically on a
     two-mode solution; control: without the d_xi H term the identity fails.
  F  focusing (numpy, d = 1, p = 4): R(kappa) = kappa int_{1/2kappa}^{1/kappa} ||e^{is Laplace} psi||_4^4 ds / ||psi||_4^4
     grows like log kappa for data focusing inside the window (one logarithm is necessary; Theorem 16.8 allows
     two).  Control: data focusing outside the window keep R bounded.

Writes ../results/null_frame.json.  Runtime about 15 s.
"""

import json
import os
import random
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import simpson

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------------------------- Part G

def part_G():
    out = {}
    ok = True
    for d in (1, 2, 3):
        Y, a, Q, S1, eta, R = sp.symbols("Y a Q S1 eta R")
        P = sp.symbols("P1:%d" % (d + 1))
        W = sp.symbols("W1:%d" % (d + 1))
        z = sp.symbols("z1:%d" % (d + 1))
        v = sp.symbols("v1:%d" % (d + 1), real=True)
        W_of = [p / Y for p in P]
        Phi_of = S1 + sum(p * p for p in P) / (2 * Y)
        P_inv = [Y * w for w in W]
        S1_inv = a * Q - Y * sum(w * w for w in W) / 2
        sub = dict(zip(P, P_inv))
        sub[S1] = S1_inv
        checks = {}
        checks["inverse chart"] = all(sp.simplify(W_of[j].subs(sub) - W[j]) == 0 for j in range(d)) and \
            sp.simplify(Phi_of.subs(sub) - a * Q) == 0
        Jm = sp.Matrix([[sp.diff(e, x) for x in list(W) + [Q]] for e in P_inv + [S1_inv]])
        checks["Jacobian a Y^d"] = sp.simplify(Jm.det() - a * Y ** d) == 0
        dPhi = [sp.diff(Phi_of, Y)] + [sp.diff(Phi_of, p) for p in P] + [sp.diff(Phi_of, S1)]
        ell = [-sum(w * w for w in W_of) / 2] + W_of + [sp.Integer(1)]
        checks["dPhi = l(W)"] = all(sp.simplify(x - y) == 0 for x, y in zip(dPhi, ell))
        th = {Y: R ** 2, S1: eta - sum(zz * zz for zz in z) / 2}
        th.update({P[j]: R * z[j] for j in range(d)})
        checks["Phi o Theta = eta"] = sp.simplify(Phi_of.subs(th) - eta) == 0
        checks["W o Theta = zeta/R"] = all(sp.simplify(W_of[j].subs(th) - z[j] / R) == 0 for j in range(d))
        Vv = [sp.Integer(1)] + [(v[j] + W_of[j]) / 2 for j in range(d)] + [-sum(W_of[j] * v[j] for j in range(d)) / 2]
        checks["tangency l(W) V^v = 0"] = sp.simplify(sum(x * y for x, y in zip(ell, Vv))) == 0
        dW_ok = True
        for j in range(d):
            dWj = [sp.diff(W_of[j], Y)] + [sp.diff(W_of[j], p) for p in P] + [sp.diff(W_of[j], S1)]
            dW_ok &= sp.simplify(sum(x * y for x, y in zip(dWj, Vv)) - (v[j] - W_of[j]) / (2 * Y)) == 0
        checks["(dW_j) V^v = (v_j - W_j)/(2Y)"] = dW_ok
        disp = [sp.Integer(1)] + list(W) + [-sum(w * w for w in W) / 2]
        checks["displacement is null"] = sp.expand(2 * disp[0] * disp[-1] + sum(w * w for w in W)) == 0
        out["d=%d" % d] = {k: bool(c) for k, c in checks.items()}
        ok &= all(checks.values())
    # negative control: the physical-frame chart (08 (4.2)) is not polynomial in (W, Q): its Jacobian has the
    # denominator (1 + W.W/2), so the simplest form above is special to the null frame.
    W1, tau, a, Q = sp.symbols("W1 tau a Q")
    Yp = (sp.sqrt(2) * tau + a * Q) / (1 + W1 * W1 / 2)
    out["control: physical-frame Y is not polynomial"] = bool(not sp.simplify(Yp).is_polynomial(W1, Q))
    ok &= out["control: physical-frame Y is not polynomial"]
    return ok, out


# ----------------------------------------------------------------------------------------------- quaternions

def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(p):
    return (p[0], -p[1], -p[2], -p[3])


def qadd(p, q):
    return tuple(x + y for x, y in zip(p, q))


def qscale(c, p):
    return tuple(c * x for x in p)


E3 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]          # slot units 1, i, j


def kdiv(x, y):
    den = y[0] ** 2 + y[1] ** 2
    return ((x[0] * y[0] + x[1] * y[1]) / den, (x[1] * y[0] - x[0] * y[1]) / den)


def kq(x):
    return (x[0], mp.mpf(0), mp.mpf(0), x[1])          # a + k b


def unpack(coords, d):
    nsl = d + 2
    X = coords[:nsl]
    om = [(coords[nsl + 2 * m], coords[nsl + 2 * m + 1]) for m in range(nsl)]
    Y, P = om[0], om[1:1 + d]
    y = mp.sqrt(Y[0] ** 2 + Y[1] ** 2)
    W = [kdiv(Pj, Y) for Pj in P]
    return X, Y, P, y, W


def chi(coords, d, modes, rotate=True):
    """Rows chi_mu, mu = (t, x_1..x_d, s0), of 08 Prop. 3.3 in the null slots; also the parts (Re_k V, Im_k V) m u."""
    nsl = d + 2
    X, Y, P, y, W = unpack(coords, d)
    e = (Y[0] / y, Y[1] / y) if rotate else (mp.mpf(1), mp.mpf(0))
    out = [(mp.mpf(0),) * 4 for _ in range(nsl)]
    parts = [[mp.mpc(0), mp.mpc(0)] for _ in range(nsl)]
    for kap, v, c in modes:
        Wv = (sum(W[j][0] * v[j] for j in range(d)), sum(W[j][1] * v[j] for j in range(d)))
        V = [(mp.mpf(1), mp.mpf(0))] + [((v[j] + W[j][0]) / 2, W[j][1] / 2) for j in range(d)] \
            + [(-Wv[0] / 2, -Wv[1] / 2)]
        dist2 = sum((W[j][0] - v[j]) ** 2 + W[j][1] ** 2 for j in range(d))
        m = (kap / y) ** (mp.mpf(d) / 2) * mp.exp(-kap * y * dist2 / 2)
        phase = kap * (X[-1] + sum(v[j] * X[1 + j] for j in range(d)) - X[0] * sum(x * x for x in v) / 2)
        u = c * mp.expj(phase)
        ubar = mp.conj(u)
        iota = (mp.re(ubar), mp.im(ubar) * e[0], mp.im(ubar) * e[1], mp.mpf(0))     # iota_e(a + i b) = a + e i b
        for mu in range(nsl):
            out[mu] = qadd(out[mu], qscale(m, qmul(kq(V[mu]), iota)))
            parts[mu][0] += V[mu][0] * m * u
            parts[mu][1] += V[mu][1] * m * u
    return out, parts


def qnorm2(p):
    return sum(x * x for x in p)


def random_modes(d, rng, n_modes, kappas=None):
    modes = []
    for i in range(n_modes):
        kap = mp.mpf(kappas[i % len(kappas)]) if kappas else mp.mpf(rng.uniform(0.6, 1.6))
        v = [mp.mpf(rng.uniform(-0.3, 0.3)) for _ in range(d)]
        c = mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1))
        modes.append((kap, v, c))
    return modes


def random_point(d, rng):
    nsl = d + 2
    X = [mp.mpf(rng.uniform(-1, 1)) for _ in range(nsl)]
    aux = []
    for mu in range(nsl):
        if mu == 0:
            aux += [mp.mpf(rng.uniform(0.6, 1.2)), mp.mpf(rng.uniform(-0.4, 0.4))]
        else:
            aux += [mp.mpf(rng.uniform(-0.3, 0.3)), mp.mpf(rng.uniform(-0.3, 0.3))]
    return X + aux


# ----------------------------------------------------------------------------------------------- Part H

def continued_wave(d, modes, t, x, s0):
    tot = mp.mpc(0)
    for kap, v, c in modes:
        ph = kap * (s0 + sum(v[j] * x[j] for j in range(d)) - t * sum(vv * vv for vv in v) / 2)
        tot += kap ** (mp.mpf(d) / 2) * c * mp.exp(1j * ph)
    return tot


def part_H(rng):
    mp.mp.dps = 50
    out = {"multi_kappa": [], "one_kappa": [], "control": []}
    ok = True
    for d in (1, 2, 3):
        for trial in range(2):
            modes = random_modes(d, rng, 4)
            pt = random_point(d, rng)
            ch, _ = chi(pt, d, modes)
            X, Y, P, y, W = unpack(pt, d)
            ReW = [w[0] for w in W]
            W2 = sum(w[0] ** 2 + w[1] ** 2 for w in W)
            lhs = qnorm2(ch[0])
            f0 = continued_wave(d, modes, X[0] - 1j * y, [X[1 + j] - 1j * y * ReW[j] for j in range(d)],
                                X[-1] + 1j * y * W2 / 2)
            rhs = abs(f0) ** 2 / y ** d
            err = abs(lhs - rhs) / rhs
            out["multi_kappa"].append({"d": d, "rel_err": float(err)})
            ok &= err < mp.mpf("1e-40")
            # one kappa: the Schrodinger solution psi continued to complex time and space
            kap = mp.mpf(rng.uniform(0.8, 1.4))
            modes1 = random_modes(d, rng, 4, kappas=[kap])
            ch1, _ = chi(pt, d, modes1)
            psi = mp.mpc(0)
            for _, v, c in modes1:
                tt = X[0] - 1j * y
                xx = [X[1 + j] - 1j * y * ReW[j] for j in range(d)]
                psi += c * mp.exp(1j * kap * (sum(v[j] * xx[j] for j in range(d)) - tt * sum(vv * vv for vv in v) / 2))
            rhs1 = (kap / y) ** d * mp.exp(-kap * y * W2) * abs(psi) ** 2
            err1 = abs(qnorm2(ch1[0]) - rhs1) / rhs1
            out["one_kappa"].append({"d": d, "rel_err": float(err1)})
            ok &= err1 < mp.mpf("1e-40")
        # control: shift x the other way
        modes = random_modes(d, rng, 4)
        pt = random_point(d, rng)
        ch, _ = chi(pt, d, modes)
        X, Y, P, y, W = unpack(pt, d)
        W2 = sum(w[0] ** 2 + w[1] ** 2 for w in W)
        f0 = continued_wave(d, modes, X[0] - 1j * y, [X[1 + j] + 1j * y * W[j][0] for j in range(d)],
                            X[-1] + 1j * y * W2 / 2)
        err = abs(qnorm2(ch[0]) - abs(f0) ** 2 / y ** d) / qnorm2(ch[0])
        out["control"].append({"d": d, "rel_err_wrong_sign": float(err)})
        ok &= err > mp.mpf("1e-3")
    return ok, out


# ----------------------------------------------------------------------------------------------- Part A

def part_A(rng):
    mp.mp.dps = 50
    out = []
    ok = True
    for d in (1, 2, 3):
        modes = random_modes(d, rng, 5)
        pt = random_point(d, rng)
        ch, parts = chi(pt, d, modes)
        errs = []
        for mu in range(d + 2):
            lhs = qnorm2(ch[mu])
            rhs = abs(parts[mu][0]) ** 2 + abs(parts[mu][1]) ** 2
            errs.append(abs(lhs - rhs) / max(rhs, mp.mpf("1e-30")))
        time_row_single = abs(parts[0][1]) == 0
        out.append({"d": d, "max_rel_err_rows": float(max(errs)), "time_row_is_one_multiplier": bool(time_row_single)})
        ok &= max(errs) < mp.mpf("1e-40") and time_row_single
    return ok, out


# ----------------------------------------------------------------------------------------------- Part C

def real_tensor(coords, d, modes, **kw):
    ch, _ = chi(coords, d, modes, **kw)
    nsl = d + 2
    Rt = [[None] * (3 * nsl) for _ in range(3 * nsl)]
    for mu in range(nsl):
        for nu in range(nsl):
            A = qmul(ch[mu], qconj(ch[nu]))
            for al in range(3):
                for be in range(3):
                    Rt[3 * mu + al][3 * nu + be] = qmul(qmul(qconj(E3[al]), A), E3[be])[0]
    return Rt


def coord_index(nu, beta, nsl):
    return nu if beta == 0 else nsl + 2 * nu + (beta - 1)


def divergence(coords, d, modes, h=mp.mpf("1e-22"), **kw):
    nsl = d + 2
    N = 3 * nsl
    div = [mp.mpf(0)] * N
    for nu in range(nsl):
        for beta in range(3):
            k = coord_index(nu, beta, nsl)
            cp, cm = list(coords), list(coords)
            cp[k] += h
            cm[k] -= h
            Rp, Rm = real_tensor(cp, d, modes, **kw), real_tensor(cm, d, modes, **kw)
            col = 3 * nu + beta
            for row in range(N):
                div[row] += (Rp[row][col] - Rm[row][col]) / (2 * h)
    R0 = real_tensor(coords, d, modes, **kw)
    scale = max(abs(x) for r in R0 for x in r)
    return max(abs(x) for x in div), scale


def s0_derivative(coords, d, modes, h=mp.mpf("1e-22")):
    nsl = d + 2
    cp, cm = list(coords), list(coords)
    cp[nsl - 1] += h
    cm[nsl - 1] -= h
    Rp, Rm = real_tensor(cp, d, modes), real_tensor(cm, d, modes)
    R0 = real_tensor(coords, d, modes)
    scale = max(abs(x) for r in R0 for x in r)
    return max(abs(Rp[i][j] - Rm[i][j]) / (2 * h) for i in range(3 * nsl) for j in range(3 * nsl)) / scale


def part_C(rng):
    mp.mp.dps = 60
    out = {"conservation_one_kappa": [], "s0_independence": [], "controls": []}
    ok = True
    for d in (1, 2, 3):
        kap = mp.mpf(rng.uniform(0.8, 1.4))
        modes = random_modes(d, rng, 3, kappas=[kap])
        pt = random_point(d, rng)
        dv, sc = divergence(pt, d, modes)
        rel = dv / sc
        out["conservation_one_kappa"].append({"d": d, "relative_divergence": float(rel)})
        ok &= rel < mp.mpf("1e-30")
        ds0 = s0_derivative(pt, d, modes)
        out["s0_independence"].append({"d": d, "relative_ds0": float(ds0)})
        ok &= ds0 < mp.mpf("1e-30")
    # controls (d = 2)
    d = 2
    modes2 = random_modes(d, rng, 4, kappas=[0.8, 1.3])
    pt = random_point(d, rng)
    ds0 = s0_derivative(pt, d, modes2)
    dv, sc = divergence(pt, d, modes2)
    out["controls"].append({"control": "two kappas: tensor depends on s0", "relative_ds0": float(ds0),
                            "still_conserved_relative_div": float(dv / sc)})
    ok &= ds0 > mp.mpf("1e-6") and dv / sc < mp.mpf("1e-30")
    kap = mp.mpf(1.1)
    modes1 = random_modes(d, rng, 3, kappas=[kap])
    dv, sc = divergence(pt, d, modes1, rotate=False)
    out["controls"].append({"control": "one kappa, no rotation iota_e: not conserved", "relative_div": float(dv / sc)})
    ok &= dv / sc > mp.mpf("1e-6")
    return ok, out


# ----------------------------------------------------------------------------------------------- Part S

def part_S():
    out = {}
    ok = True
    for d in (1, 2, 3):
        t, s0 = sp.symbols("t s0", real=True)
        kap = sp.symbols("kappa", positive=True)
        x = sp.symbols("x1:%d" % (d + 1), real=True)
        psi = sp.Function("psi")(t, *x)
        n = d + 1
        # null system: C^t = 2E00, C^{x_j} = E0j + Ej0, C^{s0} = -sum_j Ejj (08 Section 3.1)
        Ct = sp.zeros(n, n)
        Ct[0, 0] = 2
        Cx = []
        for j in range(d):
            M = sp.zeros(n, n)
            M[0, j + 1] = M[j + 1, 0] = 1
            Cx.append(M)
        Cs = sp.zeros(n, n)
        for j in range(d):
            Cs[j + 1, j + 1] = -1
        F = sp.Matrix([psi] + [sp.diff(psi, x[j]) / (sp.I * kap) for j in range(d)]) * sp.exp(sp.I * kap * s0)
        PF = Ct * sp.diff(F, t) + sum((Cx[j] * sp.diff(F, x[j]) for j in range(d)), sp.zeros(n, 1)) + Cs * sp.diff(F, s0)
        lap = sum(sp.diff(psi, xx, 2) for xx in x)
        rows_ok = all(sp.simplify(PF[j + 1]) == 0 for j in range(d))
        row0_ok = sp.simplify(PF[0] - sp.exp(sp.I * kap * s0) * (2 * sp.diff(psi, t) - sp.I / kap * lap)) == 0
        # plane wave of the null phase solves Schrodinger
        v = sp.symbols("v1:%d" % (d + 1), real=True)
        pw = sp.exp(sp.I * kap * (sum(v[j] * x[j] for j in range(d)) - t * sum(vv ** 2 for vv in v) / 2))
        schr_ok = sp.simplify(sp.diff(pw, t) - sp.I / (2 * kap) * sum(sp.diff(pw, xx, 2) for xx in x)) == 0
        # null symbol phi(xi) = |xi_x|^2/(2 xi_s0) at xi = (kappa v, kappa)
        phi_ok = sp.simplify((sum((kap * vv) ** 2 for vv in v) / (2 * kap)) - kap * sum(vv ** 2 for vv in v) / 2) == 0
        # Galilean identity on plane waves
        s = sp.symbols("s", real=True)
        k = sp.symbols("k1:%d" % (d + 1), real=True)
        xi = sp.symbols("xi1:%d" % (d + 1), real=True)
        lhs = sp.exp(sp.I * sum((k[j] + xi[j]) * x[j] for j in range(d)) - sp.I * s * sum((k[j] + xi[j]) ** 2 for j in range(d)))
        rhs = sp.exp(sp.I * sum(k[j] * x[j] for j in range(d)) - sp.I * s * sum(kk ** 2 for kk in k)) * \
            sp.exp(sp.I * sum(xi[j] * (x[j] - 2 * s * k[j]) for j in range(d)) - sp.I * s * sum(q ** 2 for q in xi))
        gal_ok = sp.simplify(sp.expand(sp.log(lhs) - sp.log(rhs), force=True)) == 0 or sp.simplify(lhs / rhs) == 1
        out["d=%d" % d] = {"system rows j>=1 vanish": bool(rows_ok), "row 0 = e^{i kappa s0}(2 d_t psi - (i/kappa) Lap psi)": bool(row0_ok),
                           "null plane wave solves Schrodinger": bool(schr_ok), "null symbol": bool(phi_ok),
                           "Galilean identity": bool(gal_ok)}
        ok &= rows_ok and row0_ok and schr_ok and phi_ok and gal_ok
    # exponent bookkeeping (Theorem 16.8): N^{d+2} * N^{-2} * N^{-d} = 1 (scaling (b) -> (c)); dyadic pieces carry a
    N, dd, a = sp.symbols("N d a", positive=True)
    out["scaling (b)->(c) exponent"] = bool(sp.simplify(sp.log(N ** (dd + 2) * N ** (-2) * N ** (-dd), N)) == 0)
    out["dyadic piece weight"] = bool(sp.simplify(a * a ** (dd / 2) * a ** (-dd / 2) - a) == 0)
    out["time change s = t/(2 kappa)"] = bool(sp.simplify(sp.integrate(1, (sp.Symbol("t"), 1, 2)) / (2 * sp.Symbol("k")) -
                                                          sp.integrate(1, (sp.Symbol("s"), 1 / (2 * sp.Symbol("k")), 1 / sp.Symbol("k")))) == 0)
    ok &= out["scaling (b)->(c) exponent"] and out["dyadic piece weight"] and out["time change s = t/(2 kappa)"]
    return bool(ok), out


# ----------------------------------------------------------------------------------------------- Part W

def part_W():
    out = {}
    ok = True
    # Gaussian identity b . grad_a G = -(beta/2) div_a grad_b G, in 1 and 2 dimensions
    alpha, beta = sp.symbols("alpha beta", positive=True)
    for dim in (1, 2):
        a = sp.symbols("a1:%d" % (dim + 1), real=True)
        b = sp.symbols("b1:%d" % (dim + 1), real=True)
        G = sp.exp(-sum(x * x for x in a) / alpha - sum(x * x for x in b) / beta)
        lhs = sum(b[i] * sp.diff(G, a[i]) for i in range(dim))
        rhs = -beta / 2 * sum(sp.diff(G, a[i], b[i]) for i in range(dim))
        good = sp.simplify(lhs - rhs) == 0
        out["gaussian identity dim=%d" % dim] = bool(good)
        ok &= good
    # the Husimi flux identity on a two-mode Schrodinger solution d_t psi = (i/2kappa) psi_xx
    x, xi, t = sp.symbols("x xi t", real=True)
    kap, w = sp.symbols("kappa w", positive=True)
    k1, k2, c1r, c1i, c2r, c2i = sp.symbols("k1 k2 c1r c1i c2r c2i", real=True)
    modes = ((c1r + sp.I * c1i, k1), (c2r + sp.I * c2i, k2))
    psi = sum(c * sp.exp(sp.I * k * x - sp.I * k ** 2 * t / (2 * kap)) for c, k in modes)
    schr = sp.simplify(sp.diff(psi, t) - sp.I / (2 * kap) * sp.diff(psi, x, 2)) == 0
    V = sum(c * sp.exp(sp.I * k * x - sp.I * k ** 2 * t / (2 * kap)) * sp.exp(-w ** 2 * (k - xi) ** 2 / 2) for c, k in modes)
    Vb = sum(sp.conjugate(c) * sp.exp(-sp.I * k * x + sp.I * k ** 2 * t / (2 * kap)) * sp.exp(-w ** 2 * (k - xi) ** 2 / 2)
             for c, k in modes)
    H = sp.expand(V * Vb)
    flux = xi / kap * H + 1 / (2 * kap * w ** 2) * sp.diff(H, xi)
    resid = sp.expand(sp.diff(H, t) + sp.diff(flux, x))
    resid_ctrl = sp.expand(sp.diff(H, t) + sp.diff(xi / kap * H, x))
    rng = random.Random(7)
    worst, worst_ctrl = 0.0, 0.0
    for _ in range(6):
        sub = {x: sp.Rational(rng.randint(-9, 9), 7), xi: sp.Rational(rng.randint(-9, 9), 5),
               t: sp.Rational(rng.randint(1, 9), 11), kap: sp.Rational(rng.randint(5, 15), 7),
               w: sp.Rational(rng.randint(3, 9), 6), k1: sp.Rational(rng.randint(-9, 9), 4),
               k2: sp.Rational(rng.randint(-9, 9), 3), c1r: sp.Rational(rng.randint(-5, 5), 3),
               c1i: sp.Rational(rng.randint(-5, 5), 2), c2r: sp.Rational(rng.randint(-5, 5), 4),
               c2i: sp.Rational(rng.randint(-5, 5), 5)}
        worst = max(worst, abs(complex(sp.N(resid.subs(sub), 60))))
        scale = abs(complex(sp.N((sp.diff(H, t)).subs(sub), 60))) + 1e-30
        worst_ctrl = max(worst_ctrl, abs(complex(sp.N(resid_ctrl.subs(sub), 60))) / scale)
    out["two-mode solution solves Schrodinger"] = bool(schr)
    out["Husimi flux identity, max |residual|"] = worst
    out["control: no d_xi H term, max relative residual"] = worst_ctrl
    ok &= schr and worst < 1e-40 and worst_ctrl > 1e-3
    return bool(ok), out


# ----------------------------------------------------------------------------------------------- Part F

def bump(r):
    out = np.zeros_like(r)
    m = np.abs(r) < 1
    out[m] = np.exp(-1.0 / (1.0 - r[m] ** 2))
    return out


def focusing_ratio(kap, s_star, L=20.0, nu=200):
    kmax = kap * L / (2 * np.pi)
    M = int(2 ** np.ceil(np.log2(8 * kmax + 64)))
    xi = np.fft.fftfreq(M, d=L / M) * 2 * np.pi
    chi_ = bump(xi / kap)

    def l4(fhat):
        f = np.fft.ifft(fhat) * M
        return np.sum(np.abs(f) ** 4) * (L / M)

    psi_hat = chi_ * np.exp(1j * s_star * xi ** 2)            # psi = e^{-i s* Lap} phi, phi_hat = bump
    n0 = l4(psi_hat)
    # s in [1/(2 kap), 1/kap]; sigma = s - s* parametrized by sinh to resolve the focal scale 1/kap^2
    lo, hi = 1 / (2 * kap) - s_star, 1 / kap - s_star
    ulo, uhi = np.arcsinh(lo * kap ** 2), np.arcsinh(hi * kap ** 2)
    us = np.linspace(ulo, uhi, 2 * nu + 1)
    vals = np.array([l4(chi_ * np.exp(-1j * (np.sinh(u) / kap ** 2) * xi ** 2)) * np.cosh(u) / kap ** 2 for u in us])
    return kap * simpson(vals, x=us) / n0


def part_F():
    kaps = [2.0 ** e for e in range(5, 12)]
    R_in = [focusing_ratio(k, 0.75 / k) for k in kaps]
    R_out = [focusing_ratio(k, -1.0 / k) for k in kaps]           # focus before the window
    inc_in = np.diff(R_in)
    inc_out = np.diff(R_out)
    last = inc_in[-3:]
    stable = (np.max(last) - np.min(last)) / np.mean(last) < 0.2
    grows = bool(np.all(inc_in > 0) and R_in[-1] > 1.5 * R_in[0])
    one_log = bool(stable)
    under_two_logs = bool(np.all(np.diff(np.array(R_in) / np.log(kaps) ** 2) < 0))
    control_bounded = bool(max(R_out) < 2 * min(R_out) and abs(inc_out[-1]) < 0.1 * inc_in[-1])
    out = {"kappa": kaps, "R_focus_inside": [float(r) for r in R_in],
           "increments_per_doubling": [float(x) for x in inc_in],
           "R_over_log_kappa": [float(r / np.log(k)) for r, k in zip(R_in, kaps)],
           "R_focus_outside_control": [float(r) for r in R_out],
           "grows": grows, "increments_stable(one log)": one_log, "R/log^2 decreasing": under_two_logs,
           "control_bounded": control_bounded}
    return grows and one_log and under_two_logs and control_bounded, out


# ----------------------------------------------------------------------------------------------- main

def main():
    rng = random.Random(16)
    t0 = time.time()
    results = {}
    ok_all = True
    for name, fn in (("G_null_chart", part_G), ("H_husimi_time_row", lambda: part_H(rng)),
                     ("A_engine_rows", lambda: part_A(rng)), ("C_conservation_one_kappa", lambda: part_C(rng)),
                     ("S_null_reduction", part_S), ("W_husimi_flux", part_W), ("F_focusing", part_F)):
        ok, out = fn()
        results[name] = {"ok": bool(ok), "detail": out}
        print("%-28s %s" % (name, "ok" if ok else "FAIL"))
        ok_all &= ok
    results["all_ok"] = bool(ok_all)
    results["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "..", "results", "null_frame.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print(json.dumps(results["F_focusing"]["detail"], indent=1, default=str))
    print("ALL OK" if ok_all else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
