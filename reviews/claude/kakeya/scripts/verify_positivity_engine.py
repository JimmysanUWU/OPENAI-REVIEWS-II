"""Checks for 11_FRAMEWORK_POSITIVITY_ENGINE.md: symmetric pencils give engine data; the wave is the arrow pencil.

Part G  eigenvector geometry of a pencil C(xi) = sum xi_j C_j near a direction omega (Lemmas 6.1-6.2):
        Hellmann-Feynman B(v(xi)) = grad phi(xi); the curvature formula D^2 phi(omega)[w,w] = 2 (Lw)^T S_perp^{-1} (Lw);
        Dv(omega) = S_perp^{-1} L; dB(0) = 2 L^T; ker L = R omega; det[dB(0), omega] != 0; stationarity of W -> xi.B(W)
        at W = v(xi); the characteristic conormal dPhi = (xi.dP - phi dtau)/(xi.omega) at W = v(xi).
Part C  the pencil current in the target coordinates (T, tau, X, P) (Lemma 6.4), for pencils with n = 2, 3, 4:
        every real column of (chi chi^*)^R is divergence free; (dPhi) chi = 0; J = |R|^{2n} j(W);
        the Gaussian identity xi.q - phi q^0 = -Psi; the row identity |chi_mu|^2 = |M'_mu g|^2 + |M''_mu g|^2.
        Negative controls: f(X + i q) instead of f(X - i q); dropping J^{-1/2}.
Part M  a middle eigenvalue: the current is still conserved, but exp(-w^* S w / 2) > 1 (the rows grow): only extremal
        sheets give a Gaussian window.
Part W  the wave is the arrow pencil (exact, sympy): the null-coordinate system of 08 §3.1, rotated to (T, b) and
        rescaled, is [I; arrow]; the arrow pencil has phi = |xi|, v(xi) = xi'/(|xi| + xi_n) (stereographic) and
        B(W) = (2W, 1 - W.W)/(1 + W.W) (inverse stereographic).

Pencils: the arrow pencils (n = 2, 3, 4), the cubic pencil of the paper's Example E.2, and a seeded random
perturbation of the arrow pencil (n = 3).  Derivatives are central differences in 60-digit arithmetic.
Writes ../results/positivity_engine.json.
"""

import json
import os
import random

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "positivity_engine.json")


# ----------------------------------------------------------------- quaternions (1, i, j, k)

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


def qnorm2(p):
    return sum(x * x for x in p)


UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]        # slot units e_0, e_1, e_2 = 1, i, j
QI = (0, 1, 0, 0)


def kq(c):
    """a + k b in C_k, stored as the mpc a + 1j b, as a quaternion."""
    return (mp.re(c), mp.mpf(0), mp.mpf(0), mp.im(c))


def iq(c):
    """a + i b in C_i (physical complex numbers) as a quaternion."""
    return (mp.re(c), mp.im(c), mp.mpf(0), mp.mpf(0))


# ----------------------------------------------------------------- pencils

def mat(rows):
    return mp.matrix([[mp.mpf(x) for x in r] for r in rows])


def arrow(n):
    """C_j = E_0j + E_j0 (j < n), C_n = diag(1, -1, ..., -1): top eigenvalue |xi|."""
    Cs = []
    for j in range(1, n):
        M = mp.zeros(n, n)
        M[0, j] = M[j, 0] = 1
        Cs.append(M)
    D = mp.zeros(n, n)
    D[0, 0] = 1
    for j in range(1, n):
        D[j, j] = -1
    Cs.append(D)
    return Cs


def cubic():
    """Example E.2: C(x, y, z) = [[z, x, y], [x, -z, 0], [y, 0, -2z]]."""
    return [mat([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), mat([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
            mat([[1, 0, 0], [0, -1, 0], [0, 0, -2]])]


def perturbed_arrow(n, seed, eps):
    rng = random.Random(seed)
    Cs = arrow(n)
    out = []
    for C in Cs:
        M = C.copy()
        for a in range(n):
            for b in range(a, n):
                x = mp.mpf(rng.uniform(-eps, eps))
                M[a, b] += x
                if a != b:
                    M[b, a] += x
        out.append(M)
    return out


def generic(n, seed):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        M = mp.zeros(n, n)
        for a in range(n):
            for b in range(a, n):
                M[a, b] = M[b, a] = mp.mpf(rng.uniform(-1, 1))
        out.append(M)
    return out


class Pencil:
    """A pencil in the adapted state basis (e_0 = chosen eigenvector of C(omega), then the others), with C^0 = I."""

    def __init__(self, Cs, omega, sheet=0):
        self.n = len(Cs)
        self.d = self.n - 1
        self.omega = [mp.mpf(x) for x in omega]
        nrm = mp.sqrt(sum(x * x for x in self.omega))
        self.omega = [x / nrm for x in self.omega]
        self.sheet = sheet                       # 0 = top eigenvalue, 1 = second, ...
        Cw = self.combo(Cs, self.omega)
        ev, Q = mp.eigsy(Cw)
        order = sorted(range(self.n), key=lambda k: -ev[k])
        first = order[sheet]
        basis = [first] + [k for k in order if k != first]
        O = mp.matrix(self.n, self.n)
        for col, k in enumerate(basis):
            for r in range(self.n):
                O[r, col] = Q[r, k]
        self.C = [O.T * C * O for C in Cs]       # C_j in the adapted basis
        self.gap = min(abs(ev[first] - ev[k]) for k in range(self.n) if k != first)

    @staticmethod
    def combo(Cs, xi):
        n = len(Cs)
        M = mp.zeros(n, n)
        for x, C in zip(xi, Cs):
            M += x * C
        return M

    def Cxi(self, xi):
        return self.combo(self.C, xi)

    def phi_E(self, xi):
        ev, Q = mp.eigsy(self.Cxi(xi))
        order = sorted(range(self.n), key=lambda k: -ev[k])
        k = order[self.sheet]
        vec = [Q[r, k] for r in range(self.n)]
        return ev[k], [x / vec[0] for x in vec]

    def phi(self, xi):
        return self.phi_E(xi)[0]

    def v(self, xi):
        return self.phi_E(xi)[1][1:]

    def S(self, xi):
        ph = self.phi(xi)
        return ph * mp.eye(self.n) - self.Cxi(xi)

    def L(self):
        """L : R^n -> e_0^perp, Lw = (C(w) e_0)_perp; column j is (C_j e_0)_perp."""
        Lm = mp.matrix(self.d, self.n)
        for j in range(self.n):
            for a in range(self.d):
                Lm[a, j] = self.C[j][a + 1, 0]
        return Lm

    # complexified group velocity map (C_k numbers stored as mpc)
    def B(self, W):
        e = [mp.mpc(1)] + list(W)
        den = 1 + sum(w * w for w in W)
        return [sum(e[a] * self.C[j][a, b] * e[b] for a in range(self.n) for b in range(self.n)) / den
                for j in range(self.n)]

    def dB(self, W):
        e = [mp.mpc(1)] + list(W)
        den = 1 + sum(w * w for w in W)
        out = mp.matrix(self.n, self.d)
        for j in range(self.n):
            Nj = sum(e[a] * self.C[j][a, b] * e[b] for a in range(self.n) for b in range(self.n))
            Ce = [sum(self.C[j][a, b] * e[b] for b in range(self.n)) for a in range(self.n)]
            for k in range(self.d):
                out[j, k] = (2 * Ce[k + 1] * den - Nj * 2 * W[k]) / den ** 2
        return out

    # fibre chart: P = tau B(W) + Phi omega
    def chart(self, tau, W, Phi):
        Bw = self.B(W)
        return [tau * Bw[j] + Phi * self.omega[j] for j in range(self.n)]

    def inverse_chart(self, tau, P, guess=None):
        W = list(guess[0]) if guess else [mp.mpc(0)] * self.d
        Phi = guess[1] if guess else sum(P[j] * self.omega[j] for j in range(self.n))
        for _ in range(60):
            F = [x - y for x, y in zip(self.chart(tau, W, Phi), P)]
            if max(abs(x) for x in F) < mp.mpf(10) ** (-(mp.mp.dps - 3)):
                break
            Jm = mp.matrix(self.n, self.n)
            dBw = self.dB(W)
            for j in range(self.n):
                for k in range(self.d):
                    Jm[j, k] = tau * dBw[j, k]
                Jm[j, self.d] = self.omega[j]
            step = mp.lu_solve(Jm, mp.matrix(F))
            W = [W[k] - step[k] for k in range(self.d)]
            Phi = Phi - step[self.d]
        return W, Phi


def maxabs(xs):
    return max(abs(x) for x in xs)


# ----------------------------------------------------------------- Part G

def part_G(name, pen, rng):
    n, d = pen.n, pen.d
    res = {"pencil": name, "n": n, "gap": float(pen.gap)}
    Lm = pen.L()
    S0 = pen.S(pen.omega)
    Sp = mp.matrix(d, d)
    for a in range(d):
        for b in range(d):
            Sp[a, b] = S0[a + 1, b + 1]
    Spinv = Sp ** -1
    errs = {"hellmann_feynman": [], "curvature": [], "Dv": [], "dB0": [], "stationarity": []}
    eps = mp.mpf("1e-20")
    for _ in range(3):
        xi = [o + mp.mpf(rng.uniform(-0.08, 0.08)) for o in pen.omega]
        # Hellmann-Feynman: B(v(xi)) = grad phi(xi), gradient by central differences
        grad = []
        for j in range(n):
            xp, xm = list(xi), list(xi)
            xp[j] += eps
            xm[j] -= eps
            grad.append((pen.phi(xp) - pen.phi(xm)) / (2 * eps))
        Bv = pen.B([mp.mpc(x) for x in pen.v(xi)])
        errs["hellmann_feynman"].append(maxabs([Bv[j] - grad[j] for j in range(n)]))
        # stationarity of W -> xi.B(W) at W = v(xi), in every direction of C^d
        v = pen.v(xi)
        for k in range(d):
            for direction in (mp.mpc(1), mp.mpc(0, 1)):
                Wp = [mp.mpc(x) for x in v]
                Wm = [mp.mpc(x) for x in v]
                Wp[k] += eps * direction
                Wm[k] -= eps * direction
                fp = sum(xi[j] * pen.B(Wp)[j] for j in range(n))
                fm = sum(xi[j] * pen.B(Wm)[j] for j in range(n))
                errs["stationarity"].append(abs((fp - fm) / (2 * eps)))
    # curvature formula and Dv(omega), for random w perpendicular to omega
    h2 = mp.mpf("1e-15")
    for _ in range(3):
        w = [mp.mpf(rng.uniform(-1, 1)) for _ in range(n)]
        dot = sum(a * b for a, b in zip(w, pen.omega))
        w = [a - dot * b for a, b in zip(w, pen.omega)]
        p0 = pen.phi(pen.omega)
        pp = pen.phi([o + h2 * a for o, a in zip(pen.omega, w)])
        pm = pen.phi([o - h2 * a for o, a in zip(pen.omega, w)])
        second = (pp - 2 * p0 + pm) / h2 ** 2
        Lw = Lm * mp.matrix(w)
        formula = 2 * (Lw.T * Spinv * Lw)[0, 0]
        errs["curvature"].append(abs(second - formula) / abs(formula))
        vp = pen.v([o + eps * a for o, a in zip(pen.omega, w)])
        vm = pen.v([o - eps * a for o, a in zip(pen.omega, w)])
        Dv = [(a - b) / (2 * eps) for a, b in zip(vp, vm)]
        target = Spinv * Lw
        errs["Dv"].append(maxabs([Dv[a] - target[a] for a in range(d)]))
    # dB(0) = 2 L^T, image perpendicular to omega, [dB(0), omega] invertible, ker L = R omega
    dB0 = pen.dB([mp.mpc(0)] * d)
    errs["dB0"].append(maxabs([dB0[j, k] - 2 * Lm[k, j] for j in range(n) for k in range(d)]))
    omega_dot_dB0 = maxabs([sum(pen.omega[j] * dB0[j, k] for j in range(n)) for k in range(d)])
    Lomega = maxabs(list(Lm * mp.matrix(pen.omega)))
    Mfull = mp.matrix(n, n)
    for j in range(n):
        for k in range(d):
            Mfull[j, k] = mp.re(dB0[j, k])
        Mfull[j, d] = pen.omega[j]
    det_chart = mp.det(Mfull)
    # transverse curvature is definite: smallest eigenvalue of 2 L^T S^-1 L on omega^perp
    Hess = 2 * Lm.T * Spinv * Lm
    hev = sorted(mp.eigsy(Hess)[0])
    res.update({k: float(max(v)) for k, v in errs.items()})
    res.update({"omega_dot_dB0": float(omega_dot_dB0), "L_omega": float(Lomega), "det_dB0_omega": float(det_chart),
                "transverse_hessian_eigs": [float(x) for x in hev]})
    ok = (all(max(v) < mp.mpf("1e-18") for k, v in errs.items() if k != "curvature")
          and max(errs["curvature"]) < mp.mpf("1e-12") and omega_dot_dB0 < mp.mpf("1e-40")
          and Lomega < mp.mpf("1e-40") and abs(det_chart) > mp.mpf("1e-3")
          and abs(hev[0]) < mp.mpf("1e-40") and hev[1] > mp.mpf("1e-3"))
    res["ok"] = bool(ok)
    return res


def conormal_check(pen, rng):
    """At W = v(xi): dPhi = (xi.dP - phi dtau)/(xi.omega) as holomorphic differentials in (tau, P)."""
    n = pen.n
    xi = [o + mp.mpf(rng.uniform(-0.05, 0.05)) for o in pen.omega]
    v = [mp.mpc(x) for x in pen.v(xi)]
    tau = mp.mpc(rng.uniform(0.6, 1.2), rng.uniform(-0.4, 0.4))
    Phi = mp.mpc(rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1))
    P = pen.chart(tau, v, Phi)
    eps = mp.mpf("1e-20")
    xo = sum(a * b for a, b in zip(xi, pen.omega))
    ph = pen.phi(xi)
    errs = []
    tp, tm = pen.inverse_chart(tau + eps, P, (v, Phi)), pen.inverse_chart(tau - eps, P, (v, Phi))
    errs.append(abs((tp[1] - tm[1]) / (2 * eps) - (-ph / xo)))
    for j in range(n):
        Pp, Pm = list(P), list(P)
        Pp[j] += eps
        Pm[j] -= eps
        dphi = (pen.inverse_chart(tau, Pp, (v, Phi))[1] - pen.inverse_chart(tau, Pm, (v, Phi))[1]) / (2 * eps)
        errs.append(abs(dphi - xi[j] / xo))
    return float(max(errs))


# ----------------------------------------------------------------- Part C: the current in target coordinates

def modes_for(pen, rng, count):
    out = []
    for _ in range(count):
        rho = mp.mpf(rng.uniform(0.8, 1.6))
        xi = [rho * (o + mp.mpf(rng.uniform(-0.06, 0.06))) for o in pen.omega]
        c = mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1))
        ph, E = pen.phi_E(xi)
        Sxi = pen.S(xi)
        out.append({"xi": xi, "c": c, "phi": ph, "E": E, "S": Sxi})
    return out


def unpack(coords, n):
    nsl = n + 1
    Xs = coords[:nsl]                                   # T, X_1..X_n
    imag = [mp.mpc(coords[nsl + 2 * m], coords[nsl + 2 * m + 1]) for m in range(nsl)]
    return Xs[0], Xs[1:], imag[0], imag[1:]


def current(pen, modes, coords, guess, sign=-1, volume=True):
    """chi_mu = J^{-1/2} sum_a M_{mu a} i G_a, with M_{mu a} = (C^mu z)_a and G = conj f(X - i q)."""
    n, d = pen.n, pen.d
    T, X, tau, P = unpack(coords, n)
    W, Phi = pen.inverse_chart(tau, P, guess)
    den = 1 + sum(w * w for w in W)
    R = mp.sqrt(2 * tau / den)
    e = [mp.mpc(1)] + W
    z = [R * x for x in e]
    Cmu = [mp.eye(n)] + pen.C
    q = [mp.re(sum(mp.conj(z[a]) * Cmu[m][a, b] * z[b] for a in range(n) for b in range(n))) / 2 for m in range(n + 1)]
    G = [mp.mpc(0)] * n
    for md in modes:
        xi, ph = md["xi"], md["phi"]
        amp = mp.sqrt(sum(x * x for x in xi)) ** (mp.mpf(d) / 2)
        expo = sum(xi[j] * q[j + 1] for j in range(n)) - ph * q[0]
        u = md["c"] * mp.expj(sum(xi[j] * X[j] for j in range(n)) - ph * T)
        fac = mp.exp(-sign * expo) if sign == 1 else mp.exp(expo)
        for a in range(n):
            G[a] += mp.conj(amp * md["E"][a] * u) * fac
    M = [[sum(Cmu[m][a, b] * z[b] for b in range(n)) for a in range(n)] for m in range(n + 1)]
    Mfull = mp.matrix(n + 1, n + 1)
    for m in range(n + 1):
        for a in range(n):
            Mfull[m, a] = M[m][a]
        Mfull[m, n] = 0 if m == 0 else pen.omega[m - 1]
    J = abs(mp.det(Mfull)) ** 2
    scale = 1 / mp.sqrt(J) if volume else mp.mpf(1)
    chi = []
    for m in range(n + 1):
        acc = (mp.mpf(0),) * 4
        for a in range(n):
            acc = qadd(acc, qmul(qmul(kq(M[m][a]), QI), iq(G[a])))
        chi.append(qscale(scale, acc))
    return chi, {"W": W, "Phi": Phi, "R": R, "z": z, "q": q, "J": J, "Mfull": Mfull, "e": e}


def real_tensor(chi):
    N = len(chi)
    Rt = [[None] * (3 * N) for _ in range(3 * N)]
    for m in range(N):
        for nu in range(N):
            A = qmul(chi[m], qconj(chi[nu]))
            for a in range(3):
                for b in range(3):
                    Rt[3 * m + a][3 * nu + b] = qmul(qmul(qconj(UNITS[a]), A), UNITS[b])[0]
    return Rt


def coord_index(nu, beta, nsl):
    return nu if beta == 0 else nsl + 2 * nu + (beta - 1)


def divergence(pen, modes, coords, guess, step=mp.mpf("1e-22"), **kw):
    nsl = pen.n + 1
    N = 3 * nsl
    div = [mp.mpf(0)] * N
    for nu in range(nsl):
        for beta in range(3):
            k = coord_index(nu, beta, nsl)
            cp, cm = list(coords), list(coords)
            cp[k] += step
            cm[k] -= step
            Rp = real_tensor(current(pen, modes, cp, guess, **kw)[0])
            Rm = real_tensor(current(pen, modes, cm, guess, **kw)[0])
            col = 3 * nu + beta
            for row in range(N):
                div[row] += (Rp[row][col] - Rm[row][col]) / (2 * step)
    R0 = real_tensor(current(pen, modes, coords, guess, **kw)[0])
    scale = max(abs(x) for r in R0 for x in r)
    return max(abs(x) for x in div), scale


def random_target(pen, rng):
    n = pen.n
    T = mp.mpf(rng.uniform(-1, 1))
    X = [mp.mpf(rng.uniform(-1, 1)) for _ in range(n)]
    tau = mp.mpc(rng.uniform(0.6, 1.2), rng.uniform(-0.4, 0.4))
    W = [mp.mpc(rng.uniform(-0.12, 0.12), rng.uniform(-0.12, 0.12)) for _ in range(pen.d)]
    Phi = mp.mpc(rng.uniform(-0.15, 0.15), rng.uniform(-0.15, 0.15))
    P = pen.chart(tau, W, Phi)
    coords = [T] + X + [mp.re(tau), mp.im(tau)]
    for p in P:
        coords += [mp.re(p), mp.im(p)]
    return coords, (W, Phi)


def tangency(pen, modes, coords, guess):
    """sum_mu (d Phi / d P_mu) chi_mu, with the holomorphic derivatives of the inverse chart."""
    n = pen.n
    chi, _ = current(pen, modes, coords, guess)
    T, X, tau, P = unpack(coords, n)
    eps = mp.mpf("1e-22")
    grads = [(pen.inverse_chart(tau + eps, P, guess)[1] - pen.inverse_chart(tau - eps, P, guess)[1]) / (2 * eps)]
    for j in range(n):
        Pp, Pm = list(P), list(P)
        Pp[j] += eps
        Pm[j] -= eps
        grads.append((pen.inverse_chart(tau, Pp, guess)[1] - pen.inverse_chart(tau, Pm, guess)[1]) / (2 * eps))
    tot = (mp.mpf(0),) * 4
    for m in range(n + 1):
        tot = qadd(tot, qmul(kq(grads[m]), chi[m]))
    size = mp.sqrt(sum(qnorm2(c) for c in chi))
    return max(abs(x) for x in tot) / size


def rows_and_identities(pen, modes, coords, guess):
    """J = |R|^{2n} j(W); xi.q - phi q^0 = -Psi; |chi_mu|^2 = |M' g|^2 + |M'' g|^2 with the explicit symbols."""
    n, d = pen.n, pen.d
    chi, info = current(pen, modes, coords, guess)
    T, X, tau, P = unpack(coords, n)
    W, R, q, e = info["W"], info["R"], info["q"], info["e"]
    y = abs(R) ** 2
    Cmu = [mp.eye(n)] + pen.C
    Dm = mp.matrix(n + 1, n + 1)
    for m in range(n + 1):
        Ce = [sum(Cmu[m][a, b] * e[b] for b in range(n)) for a in range(n)]
        for a in range(n):
            Dm[m, a] = Ce[a]
        Dm[m, n] = 0 if m == 0 else pen.omega[m - 1]
    jW = abs(mp.det(Dm)) ** 2
    err_J = abs(info["J"] - y ** n * jW) / info["J"]
    err_gauss = mp.mpf(0)
    err_rows = mp.mpf(0)
    S1 = [mp.mpc(0)] * (n + 1)
    S2 = [mp.mpc(0)] * (n + 1)
    for md in modes:
        xi, ph, E, Sx = md["xi"], md["phi"], md["E"], md["S"]
        v = E[1:]
        diff = [W[a] - v[a] for a in range(d)]
        Psi = mp.re(sum(mp.conj(diff[a]) * Sx[a + 1, b + 1] * diff[b] for a in range(d) for b in range(d))) * y / 2
        expo = sum(xi[j] * q[j + 1] for j in range(n)) - ph * q[0]
        err_gauss = max(err_gauss, abs(expo + Psi))
        rho = mp.sqrt(sum(x * x for x in xi)) ** (mp.mpf(d) / 2) * y ** (-mp.mpf(d) / 2) * mp.exp(-Psi)
        u = md["c"] * mp.expj(sum(xi[j] * X[j] for j in range(n)) - ph * T)
        for m in range(n + 1):
            cmu = sum(e[a] * Cmu[m][a, b] * E[b] for a in range(n) for b in range(n)) / mp.sqrt(jW)
            S1[m] += mp.re(cmu) * rho * u
            S2[m] += mp.im(cmu) * rho * u
    for m in range(n + 1):
        lhs = qnorm2(chi[m])
        rhs = abs(S1[m]) ** 2 + abs(S2[m]) ** 2
        err_rows = max(err_rows, abs(lhs - rhs) / max(rhs, mp.mpf("1e-30")))
    return float(err_J), float(err_gauss), float(err_rows)


def part_C(name, pen, rng, n_modes=3, trials=2, controls=True):
    out = []
    ok = True
    for trial in range(trials):
        modes = modes_for(pen, rng, n_modes)
        coords, guess = random_target(pen, rng)
        dv, sc = divergence(pen, modes, coords, guess)
        rel = dv / sc
        tg = tangency(pen, modes, coords, guess)
        eJ, eG, eR = rows_and_identities(pen, modes, coords, guess)
        rec = {"pencil": name, "n": pen.n, "trial": trial, "divergence_relative": float(rel),
               "tangency_relative": float(tg), "J_identity": eJ, "gaussian_identity": eG, "row_identity": eR}
        good = rel < mp.mpf("1e-30") and tg < mp.mpf("1e-30") and eJ < 1e-40 and eG < 1e-40 and eR < 1e-40
        if controls and trial == 0:
            d1, s1 = divergence(pen, modes, coords, guess, sign=1)
            d2, s2 = divergence(pen, modes, coords, guess, volume=False)
            rec["control_wrong_conjugation"] = float(d1 / s1)
            rec["control_no_volume_factor"] = float(d2 / s2)
            good = good and d1 / s1 > mp.mpf("1e-4") and d2 / s2 > mp.mpf("1e-4")
        rec["ok"] = bool(good)
        ok = ok and good
        out.append(rec)
    return out, ok


# ----------------------------------------------------------------- Part M: a middle eigenvalue

def part_M(rng):
    """Generic n = 3 pencil, middle sheet at a generic direction: conserved current, but the window grows."""
    for seed in range(1, 50):
        Cs = generic(3, seed)
        omega = [mp.mpf(0.3), mp.mpf(-0.2), mp.mpf(0.93)]
        pen = Pencil(Cs, omega, sheet=1)
        if pen.gap < 0.3:
            continue
        Lm = pen.L()
        S0 = pen.S(pen.omega)
        Sp = mp.matrix(2, 2)
        for a in range(2):
            for b in range(2):
                Sp[a, b] = S0[a + 1, b + 1]
        sev = sorted(mp.eigsy(Sp)[0])
        Hess = 2 * Lm.T * (Sp ** -1) * Lm
        hev = sorted(mp.eigsy(Hess)[0], key=abs)
        if abs(hev[1]) < 0.05:              # need an invertible chart: L injective on omega^perp
            continue
        res, ok_c = part_C("middle sheet (generic n=3, seed %d)" % seed, pen, rng, trials=1, controls=False)
        # the Gaussian factor exp(-w^* S w / 2) for w = R e(W): with W - v along the negative eigenvector of S_perp
        md = modes_for(pen, rng, 1)[0]
        Sx = md["S"]
        Sxp = mp.matrix(2, 2)
        for a in range(2):
            for b in range(2):
                Sxp[a, b] = Sx[a + 1, b + 1]
        ev, Q = mp.eigsy(Sxp)
        kneg = 0 if ev[0] < ev[1] else 1
        direction = [Q[0, kneg], Q[1, kneg]]
        v = md["E"][1:]
        W = [mp.mpc(v[a] + mp.mpf("0.1") * direction[a]) for a in range(2)]
        y = mp.mpf(1)
        w = [mp.sqrt(y)] + [mp.sqrt(y) * x for x in W]
        quad = mp.re(sum(mp.conj(w[a]) * Sx[a, b] * w[b] for a in range(3) for b in range(3)))
        growth = mp.exp(-quad / 2)
        ok = ok_c and sev[0] < 0 < sev[1] and growth > 1
        return {"seed": seed, "S_perp_eigs": [float(x) for x in sev], "conservation": res,
                "gaussian_factor_at_|W-v|=0.1_y=1": float(growth),
                "note": "y/h = sigma multiplies the exponent: the factor is exp(+c sigma), not a window"}, ok
    return {"error": "no admissible middle sheet found"}, False


# ----------------------------------------------------------------- Part W: the wave is the arrow pencil (exact)

def part_W():
    out = {}
    ok = True
    s2 = sp.sqrt(2)
    for n in (2, 3, 4):
        d = n - 1
        E = lambda a, b: sp.Matrix(n, n, lambda r, c: 1 if (r, c) == (a, b) else 0)
        Ct = 2 * E(0, 0)
        Cs0 = -sum((E(j, j) for j in range(1, n)), sp.zeros(n, n))
        Cx = [E(0, j) + E(j, 0) for j in range(1, n)]
        CT = (Ct - Cs0) / s2
        Cb = (Ct + Cs0) / s2
        Dm = sp.diag(*[1 / sp.sqrt(CT[k, k]) for k in range(n)])
        arrow_b = sp.diag(1, *([-1] * d))
        checks = [sp.simplify(Dm * CT * Dm - sp.eye(n)) == sp.zeros(n, n),
                  sp.simplify(Dm * Cb * Dm - arrow_b) == sp.zeros(n, n)]
        checks += [sp.simplify(Dm * C * Dm - C) == sp.zeros(n, n) for C in Cx]
        xs = sp.symbols("x1:%d" % (n + 1), real=True)
        t = sp.symbols("t")
        A = sum((xs[j - 1] * (E(0, j) + E(j, 0)) for j in range(1, n)), sp.zeros(n, n)) + xs[-1] * arrow_b
        r2 = sum(x ** 2 for x in xs)
        charpoly = sp.expand(sp.det(t * sp.eye(n) - A))
        target = sp.expand((t + xs[-1]) ** (n - 2) * (t ** 2 - r2))
        checks.append(sp.simplify(charpoly - target) == 0)
        rr = sp.sqrt(r2)
        Evec = sp.Matrix([1] + [xs[j] / (rr + xs[-1]) for j in range(d)])
        checks.append(sp.simplify(A * Evec - rr * Evec) == sp.zeros(n, 1))
        Ws = sp.symbols("W1:%d" % n)
        e = sp.Matrix([1] + list(Ws))
        WW = sum(w ** 2 for w in Ws)
        Bw = [sp.simplify((e.T * (E(0, j) + E(j, 0)) * e)[0] / (1 + WW)) for j in range(1, n)]
        Bw.append(sp.simplify((e.T * arrow_b * e)[0] / (1 + WW)))
        stereo = [2 * w / (1 + WW) for w in Ws] + [(1 - WW) / (1 + WW)]
        checks.append(all(sp.simplify(a - b) == 0 for a, b in zip(Bw, stereo)))
        sub = {Ws[j]: xs[j] / (rr + xs[-1]) for j in range(d)}
        checks.append(all(sp.simplify(sp.radsimp(Bw[j].subs(sub)) - xs[j] / rr) == 0 for j in range(n)))
        out["n=%d" % n] = [bool(c) for c in checks]
        ok = ok and all(checks)
    return out, ok


# ----------------------------------------------------------------- main

def main():
    rng = random.Random(20261028)
    results = {"G": [], "C": [], "M": None, "W": None}
    all_ok = True
    pencils = [("arrow n=2 (Dirac)", Pencil(arrow(2), [0, 1])),
               ("arrow n=3", Pencil(arrow(3), [0, 0, 1])),
               ("arrow n=4", Pencil(arrow(4), [0, 0, 0, 1])),
               ("cubic (Example E.2)", Pencil(cubic(), [0, 0, 1])),
               ("perturbed arrow n=3", Pencil(perturbed_arrow(3, 7, 0.15), [0, 0, 1]))]
    for name, pen in pencils:
        r = part_G(name, pen, rng)
        r["conormal"] = conormal_check(pen, rng)
        r["ok"] = r["ok"] and r["conormal"] < 1e-18
        results["G"].append(r)
        all_ok = all_ok and r["ok"]
        print("G %-24s HF %.1e curv %.1e Dv %.1e dB0 %.1e stat %.1e conormal %.1e  %s"
              % (name, r["hellmann_feynman"], r["curvature"], r["Dv"], r["dB0"], r["stationarity"], r["conormal"],
                 "ok" if r["ok"] else "FAIL"))
    cubic_pen = pencils[3][1]
    Hc = 2 * cubic_pen.L().T * (lambda S0: mp.matrix([[S0[1, 1], S0[1, 2]], [S0[2, 1], S0[2, 2]]]))(cubic_pen.S(cubic_pen.omega)) ** -1 * cubic_pen.L()
    cubic_hess = [[float(Hc[a, b]) for b in range(3)] for a in range(3)]
    cubic_ok = abs(Hc[0, 0] - 1) < 1e-40 and abs(Hc[1, 1] - mp.mpf(2) / 3) < 1e-40 and abs(Hc[0, 1]) < 1e-40
    results["cubic_transverse_hessian"] = cubic_hess
    print("G cubic transverse Hessian at (0,0,1): diag(%.6f, %.6f) (paper: diag(1, 2/3))  %s"
          % (cubic_hess[0][0], cubic_hess[1][1], "ok" if cubic_ok else "FAIL"))
    all_ok = all_ok and bool(cubic_ok)
    for name, pen in pencils[:4] + pencils[4:]:
        recs, ok = part_C(name, pen, rng, trials=2 if pen.n <= 3 else 1)
        results["C"] += recs
        all_ok = all_ok and ok
        for r in recs:
            extra = ""
            if "control_wrong_conjugation" in r:
                extra = "  controls: conj %.1e, no J %.1e" % (r["control_wrong_conjugation"], r["control_no_volume_factor"])
            print("C %-24s div %.1e tang %.1e J %.1e gauss %.1e rows %.1e%s  %s"
                  % (name, r["divergence_relative"], r["tangency_relative"], r["J_identity"], r["gaussian_identity"],
                     r["row_identity"], extra, "ok" if r["ok"] else "FAIL"))
    m, ok = part_M(rng)
    results["M"] = m
    all_ok = all_ok and ok
    if "error" in m:
        print("M", m["error"])
    else:
        print("M middle sheet: S_perp eigenvalues %s, conserved (div %.1e), Gaussian factor %.3f > 1  %s"
              % (["%.3f" % x for x in m["S_perp_eigs"]], m["conservation"][0]["divergence_relative"],
                 m["gaussian_factor_at_|W-v|=0.1_y=1"], "ok" if ok else "FAIL"))
    w, ok = part_W()
    results["W"] = w
    all_ok = all_ok and ok
    print("W wave = arrow pencil (exact): %s  %s" % (w, "ok" if ok else "FAIL"))
    results["all_ok"] = bool(all_ok)
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(results, fh, indent=1, default=float)
    print("ALL OK" if all_ok else "FAILED")


if __name__ == "__main__":
    main()
