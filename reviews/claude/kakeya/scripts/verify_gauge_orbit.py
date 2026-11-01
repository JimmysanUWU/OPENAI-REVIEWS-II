"""Checks for 20_GAUGE_ORBIT_OBSTRUCTION.md: GPT 18-19 re-derived, extended to the complex target and to every
Hermitian pencil.

Notation: auxiliary states z = w1 + k w2 in C_k^N are stored as complex numpy arrays (1j plays the role of k);
Q(z) = (z^T C^mu z / 2)_mu is the bilinear flux map, q(z) = ((w1^T C w1 + w2^T C w2)/2)_mu the real displacement.

  W  [Proposition 20.1] the realified Weyl pencil of GPT 18:
     - the (lambda, a, b, c) chart reproduces Q(z) = ((a+d)/2, (b+c)/2, (c-b)/(2k), (a-d)/2), d = bc/a;
     - the coarea density 1/(|a-d|^2 |lambda|^2) by finite-difference real Jacobians (10x10 and 8x8), not by
       the symbolic determinant used in GPT's script;
     - q^0 - q^3 = (|u2|^2 + |v2|^2)/2 and the time entry |z1|^2 exp(-2 kappa (q^0 - q^3)) from the matrices;
     - the closed complex profile F = pi |a|^2 J(kappa^2|b|^2|c|^2/|a|^2) (|b|^-2 + |c|^-2) / (4 kappa |a-d|^2),
       J(v) = 2 sqrt(v) K_1(2 sqrt(v)), against 2D quadrature of the fibre integral at random complex (a, b, c);
     - on real targets it reduces to GPT 19 (19.4); |b|^2 F -> pi |a|^2/(4 kappa |a-d|^2) as b -> 0 with c
       fixed (the singular set is a complex hypersurface, so F is not locally integrable on C_k^4);
     - with b = 0 or c = 0 one end of the radial integral is undamped: the truncated integrals grow.
  H  [Theorem 20.2] random complex Hermitian pencils m' x m', n = 2m' - 1 (m' = 2, 3, 4), realified:
     - J (multiplication by i) commutes with every C^mu, so g_theta = cos(theta) + sin(theta) J (theta complex)
       preserves Q and has det 1;
     - along the orbit g_theta z0 of a real unit top eigenvector: S g_theta z0 = 0 (S = phi I - C(xi0)), so the
       plane-wave damping exponent phi q^0 - xi0.q' vanishes identically, and the time entry is
       cos^2(a) cosh^2(t) + sin^2(a) sinh^2(t), theta = a + k t;
     - the coarea density of Theta(z, eta) = Q(z) + c eta in the orbit coordinate theta is the same constant at
       every orbit point (complex Jacobians, any complement); for the Weyl pencil with c = e_3 it equals 1,
       matching GPT 18's d^2 lambda/|lambda|^2 with lambda = e^{k theta};
     - the orbit contribution over |Im theta| <= T equals rho * pi * sinh(2T) (quadrature).
  Q  [Theorem 20.2, quaternionic case] random quaternionic Hermitian pencils 2 x 2 and 3 x 3, n = 4m'' - 3:
     - right multiplications I_1, I_2, I_3 by i, j, k commute with every C^mu and satisfy the quaternion rules;
     - for complex theta in C_k^3, g = exp(sum theta_a I_a) preserves Q, has det 1, keeps g z0 in E(xi0)_C (no
       damping), and <g z0, z0> = cos(s), s^2 = sum theta_a^2;
     - the coarea density in the frame (I_a x) is constant along the 3-dimensional orbit;
     - the average over real unit quaternions u of <u z0, z0>^2 is 1/4 (the polar-coordinate lower bound).
  N  controls: a random antisymmetric "gauge" matrix does not commute and does not preserve Q; a point pushed
     off the orbit transversally has positive damping; the real symmetric arrow pencil (N = n, m = 1) has no
     continuous fibre (Theta is a local diffeomorphism, kernel 0); at b, c != 0 the Weyl integral converges.

Requires numpy and mpmath.  Deterministic; writes ../results/gauge_orbit.json (about 10 s).
"""

import json
import math
import os
import random
import time

import mpmath as mp
import numpy as np

RESULTS = {}
RS = np.random.default_rng(20261101)


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def Qmap(Cs, z):
    return np.array([0.5 * z @ C @ z for C in Cs])


def qdisp(Cs, z):
    w1, w2 = z.real, z.imag
    return np.array([0.5 * (w1 @ C @ w1 + w2 @ C @ w2) for C in Cs])


def real_jac(f, x, h=1e-6):
    """Real Jacobian of f: R^p -> R^q by central differences."""
    x = np.asarray(x, dtype=float)
    cols = []
    for i in range(len(x)):
        e = np.zeros_like(x)
        e[i] = h
        cols.append((f(x + e) - f(x - e)) / (2 * h))
    return np.array(cols).T


def c2r(v):
    v = np.asarray(v, dtype=complex)
    return np.concatenate([v.real, v.imag])


def r2c(x):
    n = len(x) // 2
    return x[:n] + 1j * x[n:]


# =========================================================================== Part W: Weyl

def weyl():
    C0 = np.eye(4)
    C1 = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], float)
    C2 = np.array([[0, 0, 0, 1], [0, 0, -1, 0], [0, -1, 0, 0], [1, 0, 0, 0]], float)
    C3 = np.diag([1., 1., -1., -1.])
    return [C0, C1, C2, C3]


def z_from_chart(lam, a, b, c):
    u1, u2, v1, v2 = lam, c * lam / a, a / lam, b / lam
    return np.array([(u1 + v1) / 2, (u1 - v1) / (2j), (u2 + v2) / 2, (u2 - v2) / (2j)])


def J_profile(v):
    if v == 0:
        return 1.0
    s = 2 * math.sqrt(v)
    return float(s * mp.besselk(1, s))


def F_closed(a, b, c, kappa):
    d = b * c / a
    return (math.pi * abs(a) ** 2 * J_profile(kappa ** 2 * abs(b) ** 2 * abs(c) ** 2 / abs(a) ** 2)
            * (1 / abs(b) ** 2 + 1 / abs(c) ** 2) / (4 * kappa * abs(a - d) ** 2))


def F_quad(Cs, a, b, c, kappa, nr=1200, nth=64, rmin=1e-4, rmax=1e4):
    """Fibre integral of the time entry against d^2 lambda/(|a-d|^2 |lambda|^2), lambda = r e^{k th}."""
    d = b * c / a
    us = np.linspace(math.log(rmin), math.log(rmax), nr)
    ths = np.linspace(0, 2 * math.pi, nth, endpoint=False)
    total = 0.0
    for uu in us:
        r = math.exp(uu)
        acc = 0.0
        for th in ths:
            z = z_from_chart(r * np.exp(1j * th), a, b, c)
            q = qdisp(Cs, z)
            acc += abs(z[0]) ** 2 * math.exp(-2 * kappa * (q[0] - q[3]))
        total += acc * (2 * math.pi / nth)
    du = us[1] - us[0]
    return total * du / abs(a - d) ** 2           # d^2 lambda/|lambda|^2 = d(log r) d theta


def part_W():
    Cs = weyl()
    out = {}
    # Clifford relations
    for i in range(1, 4):
        for j in range(1, 4):
            check(np.allclose(Cs[i] @ Cs[j] + Cs[j] @ Cs[i], 2 * (i == j) * np.eye(4)), "Clifford")
    # chart reproduces Q
    for _ in range(20):
        lam, a, b, c = [complex(*RS.normal(size=2)) for _ in range(4)]
        z = z_from_chart(lam, a, b, c)
        d = b * c / a
        target = np.array([(a + d) / 2, (b + c) / 2, (c - b) / (2j), (a - d) / 2])
        check(np.allclose(Qmap(Cs, z), target), "chart reproduces Q(z)")
        q = qdisp(Cs, z)
        u2, v2 = c * lam / a, b / lam
        check(abs((q[0] - q[3]) - (abs(u2) ** 2 + abs(v2) ** 2) / 2) < 1e-9, "q0 - q3 = (|u2|^2+|v2|^2)/2")
    # coarea density by finite-difference real Jacobians
    dens = []
    for _ in range(8):
        lam, a, b, c = [complex(*RS.normal(size=2)) for _ in range(4)]
        a = a + 2.0
        eta = complex(*RS.normal(size=2))
        x = c2r([lam, a, b, c, eta])

        def chart(x):
            lam, a, b, c, eta = r2c(x)
            return c2r(np.concatenate([z_from_chart(lam, a, b, c), [eta]]))

        def target(x):
            a, b, c, eta = r2c(x)
            d = b * c / a
            return c2r([(a + d) / 2, (b + c) / 2, (c - b) / (2j), (a - d) / 2 + eta])

        J10 = abs(np.linalg.det(real_jac(chart, x)))
        J8 = abs(np.linalg.det(real_jac(target, c2r([a, b, c, eta]))))
        predicted = 1 / (abs(a - b * c / a) ** 2 * abs(lam) ** 2)
        check(abs(J10 / J8 / predicted - 1) < 1e-6, "coarea density 1/(|a-d|^2|lambda|^2)")
        dens.append(J10 / J8 / predicted)
    out["coarea_ratio_to_formula"] = [float(x) for x in dens]
    # closed profile against quadrature at complex points
    prof = []
    kappa = 1.3
    for _ in range(4):
        a = complex(1.0 + 0.3 * RS.normal(), 0.3 * RS.normal())
        b = complex(*(0.4 * RS.normal(size=2)))
        c = complex(*(0.4 * RS.normal(size=2)))
        Fq = F_quad(Cs, a, b, c, kappa)
        Fc = F_closed(a, b, c, kappa)
        check(abs(Fq / Fc - 1) < 2e-3, "complex profile vs quadrature")
        prof.append(dict(a=str(a), b=str(b), c=str(c), quad=Fq, closed=Fc))
    out["profile_complex"] = prof
    # reduction to GPT 19 on real targets
    for _ in range(6):
        w0 = 0.5 + 0.3 * RS.random()
        s = 0.2 * RS.random() + 1e-3
        ang = 2 * math.pi * RS.random()
        w1, w2 = s * math.cos(ang), s * math.sin(ang)
        a = w0 + math.sqrt(w0 ** 2 - s ** 2)
        d = w0 - math.sqrt(w0 ** 2 - s ** 2)
        b, c = complex(w1, -w2), complex(w1, w2)
        check(abs(b * c / a - d) < 1e-12, "real slice: d = bc/a")
        gpt = math.pi * a ** 2 * J_profile(kappa ** 2 * s ** 4 / a ** 2) / (2 * kappa * s ** 2 * (a - d) ** 2)
        check(abs(F_closed(a, b, c, kappa) / gpt - 1) < 1e-12, "complex profile restricts to GPT (19.4)")
    # transverse singularity in the complex target: |b|^2 F -> pi|a|^2/(4 kappa |a-d|^2) as b -> 0, c fixed
    a, c = complex(1.1, 0.2), complex(0.3, -0.2)
    lim = []
    for eps in [1e-2, 1e-3, 1e-4, 1e-5]:
        b = eps * complex(0.6, 0.8)
        val = abs(b) ** 2 * F_closed(a, b, c, kappa)
        lim.append(val)
    target = math.pi * abs(a) ** 2 / (4 * kappa * abs(a) ** 2)
    check(abs(lim[-1] / target - 1) < 1e-4, "|b|^2 F limit")
    out["b_to_0_limit"] = dict(values=lim, predicted=target)
    # undamped end when c = 0: truncated radial integrals grow like R^2
    a, b = 1.0 + 0j, 0.3 + 0.1j
    grow = []
    for R in [10, 30, 100]:
        F = F_quad(Cs, a, b, 0j, kappa, nr=600, nth=16, rmin=1e-3, rmax=R)
        grow.append(F)
    check(grow[2] > 50 * grow[0], "c = 0: integral grows with the cutoff")
    out["c_zero_growth"] = grow
    RESULTS["W"] = out
    print("Part W: Weyl coarea (finite differences), complex profile (quadrature), real-slice reduction, "
          "complex-hypersurface singularity: OK")


# =========================================================================== Part H: complex Hermitian

def realify_complex(H):
    A, B = H.real, H.imag
    return np.block([[A, -B], [B, A]])


def random_hermitian(m):
    X = RS.normal(size=(m, m)) + 1j * RS.normal(size=(m, m))
    return (X + X.conj().T) / 2


def top_eigen(Cs, xi):
    Cx = sum(x * C for x, C in zip(xi, Cs[1:]))
    w, U = np.linalg.eigh(Cx)
    return w, U, Cx


def coarea_density(Cs, cdir, x, tangents, E):
    """Fibre density of Theta(z, eta) = Q(z) + c eta at the source point (x, eta) in the coordinates whose
    tangent vectors are `tangents` (columns, source space C^{N+1}); E: complement (N+1) x (n+1)."""
    Mt = np.array([C @ x for C in Cs])                     # (n+1) x N complex Jacobian of Q
    DTheta = np.concatenate([Mt, cdir[:, None]], axis=1)    # (n+1) x (N+1)
    num = abs(np.linalg.det(np.concatenate([tangents, E], axis=1))) ** 2
    den = abs(np.linalg.det(DTheta @ E)) ** 2
    return num / den


def part_H():
    out = []
    for m in (2, 3, 4):
        n = 2 * m - 1
        if m == 2:
            Cs = weyl()                      # the Weyl pencil is the m' = 2 member, in GPT's real form
            J = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]], float)
            cdir = np.array([0, 0, 0, 1], complex)
            xi0 = np.array([0., 0., 1.])
        else:
            Hs = [random_hermitian(m) for _ in range(n)]
            Cs = [np.eye(2 * m)] + [realify_complex(H) for H in Hs]
            J = np.block([[np.zeros((m, m)), -np.eye(m)], [np.eye(m), np.zeros((m, m))]])
            cdir = RS.normal(size=n + 1) + 1j * RS.normal(size=n + 1)
            xi0 = RS.normal(size=n)
            xi0 /= np.linalg.norm(xi0)
        N = Cs[0].shape[0]
        check(np.allclose(J @ J, -np.eye(N)) and np.allclose(J.T, -J), "complex structure")
        for C in Cs:
            check(np.allclose(J @ C, C @ J), "J commutes with the pencil")
        w, U, Cx = top_eigen(Cs, xi0)
        phi = w[-1]
        check(w[-1] - w[-3] > 1e-6 and abs(w[-1] - w[-2]) < 1e-9, "top eigenvalue: real multiplicity exactly 2")
        z0 = U[:, -1]
        S = phi * np.eye(N) - Cx
        E = np.concatenate([np.eye(N + 1, dtype=complex)[:, 1:n + 2]], axis=1)   # a fixed complement
        dens = []
        for _ in range(12):
            al, t = 2 * math.pi * RS.random(), 3 * RS.normal()
            th = complex(al, t)
            g = math.cos(0) * 0 + np.cos(th) * np.eye(N) + np.sin(th) * J
            check(abs(np.linalg.det(g) - 1) < 1e-8 * max(1, abs(np.cos(th)) ** N), "det g = 1")
            z = g @ z0
            check(np.allclose(Qmap(Cs, z), Qmap(Cs, z0), atol=1e-8 * (1 + abs(np.cos(th)) ** 2)), "Q invariant")
            check(np.allclose(S @ z, 0, atol=1e-8 * (1 + abs(np.cos(th)))), "orbit stays in E(xi0)_C")
            q = qdisp(Cs, z)
            damp = phi * q[0] - xi0 @ q[1:]
            check(abs(damp) < 1e-7 * (1 + abs(z) @ abs(z)), "no damping on the orbit")
            time_entry = abs(z @ z0) ** 2
            check(abs(time_entry - (math.cos(al) ** 2 * math.cosh(t) ** 2 + math.sin(al) ** 2 * math.sinh(t) ** 2))
                  < 1e-8 * math.cosh(t) ** 2, "time entry along the orbit")
            tangent = np.concatenate([J @ z, [0]])[:, None]
            dens.append(coarea_density(Cs, cdir, z, tangent, E))
        dens = np.array(dens)
        check(np.max(np.abs(dens / dens[0] - 1)) < 1e-6, "coarea density constant along the orbit")
        if m == 2:
            check(abs(dens[0] - 1) < 1e-9, "Weyl: density 1 in theta (GPT's d^2 lambda/|lambda|^2)")
        # orbit contribution over |Im theta| <= T: rho * pi * sinh(2T) (quadrature in alpha, t)
        T = 2.0
        al = np.linspace(0, 2 * math.pi, 64, endpoint=False)
        ts = np.linspace(-T, T, 801)
        vals = []
        for t in ts:
            row = [abs((np.cos(complex(a_, t)) * z0 + np.sin(complex(a_, t)) * (J @ z0)) @ z0) ** 2 for a_ in al]
            vals.append(np.mean(row) * 2 * math.pi)
        integral = dens[0] * np.trapezoid(vals, ts)
        check(abs(integral / (dens[0] * math.pi * math.sinh(2 * T)) - 1) < 1e-4, "orbit integral rho pi sinh 2T")
        out.append(dict(mprime=m, n=n, N=N, density=float(dens[0]), orbit_integral_T2=float(integral)))
    RESULTS["H"] = out
    print("Part H: complex Hermitian pencils m' = 2, 3, 4: gauge orbit, no damping, constant coarea density, "
          "pi sinh(2T) growth: OK")


# =========================================================================== Part Q: quaternionic

def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def lmat(q):
    return np.array([qmul(q, e) for e in np.eye(4)]).T


def rmat(q):
    return np.array([qmul(e, q) for e in np.eye(4)]).T


def random_hermitian_quaternion(mq):
    M = np.zeros((4 * mq, 4 * mq))
    for i in range(mq):
        for j in range(i, mq):
            if i == j:
                q = np.array([RS.normal(), 0, 0, 0])
                M[4 * i:4 * i + 4, 4 * i:4 * i + 4] = lmat(q)
            else:
                q = RS.normal(size=4)
                qc = q * np.array([1, -1, -1, -1])
                M[4 * i:4 * i + 4, 4 * j:4 * j + 4] = lmat(q)
                M[4 * j:4 * j + 4, 4 * i:4 * i + 4] = lmat(qc)
    return M


def expm(A):
    w, V = np.linalg.eig(A)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


def part_Q():
    out = []
    units = [np.array(u, float) for u in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1])]
    for mq in (2, 3):
        n = 4 * mq - 3
        Cs = [np.eye(4 * mq)] + [random_hermitian_quaternion(mq) for _ in range(n)]
        N = 4 * mq
        for C in Cs:
            check(np.allclose(C, C.T), "realified quaternionic Hermitian is symmetric")
        Is = [np.kron(np.eye(mq), rmat(u)) for u in units]
        for I_ in Is:
            check(np.allclose(I_.T, -I_) and np.allclose(I_ @ I_, -np.eye(N)), "unit right multiplications")
            for C in Cs:
                check(np.allclose(I_ @ C, C @ I_), "right multiplication commutes with the pencil")
        check(np.allclose(Is[0] @ Is[1], -Is[2]) or np.allclose(Is[0] @ Is[1], Is[2]), "quaternion relation")
        xi0 = RS.normal(size=n)
        xi0 /= np.linalg.norm(xi0)
        w, U, Cx = top_eigen(Cs, xi0)
        phi = w[-1]
        check(abs(w[-1] - w[-4]) < 1e-9 and w[-1] - w[-5] > 1e-6, "top eigenvalue: real multiplicity exactly 4")
        z0 = U[:, -1]
        S = phi * np.eye(N) - Cx
        cdir = RS.normal(size=n + 1) + 1j * RS.normal(size=n + 1)
        E = np.eye(N + 1, dtype=complex)[:, 3:n + 4]
        dens = []
        for _ in range(10):
            th = RS.normal(size=3) + 1j * RS.normal(size=3)
            A = sum(t_ * I_ for t_, I_ in zip(th, Is))
            g = expm(A)
            check(abs(np.linalg.det(g) - 1) < 1e-6, "det g = 1")
            z = g @ z0
            check(np.allclose(Qmap(Cs, z), Qmap(Cs, z0), atol=1e-7), "Q invariant under the complexified Sp(1)")
            check(np.allclose(S @ z, 0, atol=1e-7), "orbit stays in E(xi0)_C")
            q = qdisp(Cs, z)
            check(abs(phi * q[0] - xi0 @ q[1:]) < 1e-7 * (1 + abs(z) @ abs(z)), "no damping")
            s = np.sqrt(np.sum(th ** 2))
            check(abs(z @ z0 - np.cos(s)) < 1e-8 * (1 + abs(np.cos(s))), "<g z0, z0> = cos s")
            tangents = np.array([np.concatenate([I_ @ z, [0]]) for I_ in Is]).T
            dens.append(coarea_density(Cs, cdir, z, tangents, E))
        dens = np.array(dens)
        check(np.max(np.abs(dens / dens[0] - 1)) < 1e-5, "coarea density constant along the Sp(1)_C orbit")
        # polar lower bound: mean over unit quaternions u of <u z0, z0>^2 = 1/4
        us = RS.normal(size=(20000, 4))
        us /= np.linalg.norm(us, axis=1)[:, None]
        mean_sq = float(np.mean(us[:, 0] ** 2))
        vals = []
        for u in us[:2000]:
            g = np.kron(np.eye(mq), rmat(u))
            vals.append((z0 @ g @ z0) ** 2)
        check(abs(np.mean(vals) - 0.25) < 0.02 and abs(mean_sq - 0.25) < 0.01, "E <u z0, z0>^2 = 1/4")
        out.append(dict(mq=mq, n=n, N=N, density=float(dens[0]), mean_sq=float(np.mean(vals))))
    RESULTS["Q"] = out
    print("Part Q: quaternionic Hermitian pencils 2x2, 3x3: Sp(1)_C orbit, no damping, cos s, constant density: OK")


# =========================================================================== Part N: controls

def part_N():
    out = {}
    Cs = weyl()
    A = RS.normal(size=(4, 4))
    A = A - A.T
    out["fake_gauge_commutes"] = bool(all(np.allclose(A @ C, C @ A) for C in Cs))
    check(not out["fake_gauge_commutes"], "control: a random antisymmetric matrix should not commute")
    th = 0.7 + 0.4j
    g = expm(th * A)
    z0 = np.array([1, 0, 0, 0], complex)
    out["fake_gauge_preserves_Q"] = bool(np.allclose(Qmap(Cs, g @ z0), Qmap(Cs, z0)))
    check(not out["fake_gauge_preserves_Q"], "control: a fake gauge should move Q")
    # off-orbit point: damping becomes positive
    z = z0 + 0.1 * np.array([0, 0, 1, 0.5j])
    q = qdisp(Cs, z)
    out["off_orbit_damping"] = float(q[0] - q[3])
    check(q[0] - q[3] > 1e-3, "control: off-orbit damping is positive")
    # arrow pencil (N = n, m = 1): Theta has trivial kernel, so no continuous fibre
    n = 3
    Ca = [np.eye(3)]
    for j in range(n - 1):
        C = np.zeros((3, 3))
        C[j, n - 1] = C[n - 1, j] = 1
        Ca.append(C)
    C = np.zeros((3, 3))
    C[n - 1, n - 1] = 1
    Ca.append(C)
    z = RS.normal(size=3) + 1j * RS.normal(size=3)
    Mt = np.array([C @ z for C in Ca])
    out["arrow_rank_DQ"] = int(np.linalg.matrix_rank(Mt))
    check(out["arrow_rank_DQ"] == 3, "control: arrow pencil flux map has full rank N (no fibre directions)")
    # Weyl off the axis (b, c != 0): finite and equal to the closed profile (checked in Part W); record one value
    out["weyl_offaxis_F"] = F_closed(1.0 + 0j, 0.2 + 0.1j, 0.1 - 0.3j, 1.0)
    RESULTS["N"] = out
    print("Part N: controls (fake gauge, off-orbit damping, arrow pencil has no fibre): OK")


def main():
    t0 = time.time()
    part_W()
    part_H()
    part_Q()
    part_N()
    RESULTS["seconds"] = round(time.time() - t0, 1)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "gauge_orbit.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK")


if __name__ == "__main__":
    main()
