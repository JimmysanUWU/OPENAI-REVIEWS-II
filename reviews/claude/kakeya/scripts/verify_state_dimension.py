"""Checks for 17_STATE_DIMENSION_RIGIDITY.md: which symmetric systems can carry a tangent lifted current.

For a pencil C^0 = I, C^1..C^n (real symmetric N x N) the lifted current pushed through the flux map
Q(z) = (z^T C^mu z / 2)_mu admits a passive coordinate iff M(z) = [C^0 z | ... | C^n z] has rank <= n for all
complex z (Theorem 17.2).  With an extremal sheet of constant multiplicity m this forces N' = n + m - 1
(Theorem 17.3), N' the dimension of the span of the extremal eigenvectors.

  R  rank criterion:
     - arrow pencils n = 2, 3, 4 (N = n, m = 1): the flux map lands on the light cone (exact polynomial identity),
       rank M(z) = n at random complex z;
     - random pencils with N = n + 1 and a simple top eigenvalue: det M(z) is a nonzero polynomial (exact) and
       rank M(z) = n + 1 at random z, so no passive coordinate exists;
     - the arrow pencil weakly coupled to one extra state (n = 3, N = 4): the top eigenvalue stays simple and its
       eigenvectors span R^4, and rank M(z) = 4;
     - the uncoupled direct sum arrow (+) [one scalar]: N = n + 1, but the top eigenvectors span only V = R^n
       (N' = n), and the pencil compressed to V has rank n: N > n is harmless when N' = n;
     - the fibre identity z^T S(xi0) z = 2(phi(xi0) Q_0(z) - xi0 . Q'(z)) (exact).
  H  Hopf (Clifford) pencils C(r, a) = [[r, L_a^T], [L_a, -r]] over R, C, H, O (k = 1, 2, 4, 8; n = k + 1, N = 2k):
     - anticommutation C_j C_l + C_l C_j = 2 delta_jl I (exact);
     - the Hopf identity sum_j (z^T C_j z)^2 = (z^T z)^2 (exact polynomial identity);
     - the top eigenvalue |xi| has multiplicity m = k, so N = n + m - 1 holds with equality;
     - rank M(z) = n at random complex z;
     - the passive coordinate Phi(tau, P) = w.P - sqrt(tau^2 - P.P + (w.P)^2) satisfies Phi o Theta = eta, and
       tangency sum_mu dPhi_mu C^mu z = 0 holds at Theta(z, eta) (complex points, 40 digits);
     - the polarization is the Hopf map: C(h(z)) z = z for real z, h(z) = (z^T C_j z)_j / |z|^2.
     Controls: four of the five quaternionic matrices (n = 4, N = 8: N > n + m - 1 = 7) fail the Hopf identity and
     have rank n + 1; a generic perturbation of the k = 2 pencil has rank n + 1; arrow (+) arrow (N = 2n, m = 2)
     has rank n + 1 for one lift.
  Q  division-algebra pencils (numerical): random complex Hermitian pencils (m' x m', realified, m = 2) with
     n = 2m' - 1 and random quaternionic Hermitian pencils (m'' x m'', m = 4) with n = 4m'' - 3 satisfy
     N' = n + m - 1, are transversally curved at a sample direction, and have rank M(z) = n, so they pass the
     tangency test; with one direction fewer the rank is n' + 1 (fail).  Maxwell's equations (N = 6, n = 3,
     |xi| of multiplicity 2, N' = 6 > n + m - 1 = 4) fail: rank M(z) = 4.
  W  Hurwitz-Radon bookkeeping: N = 2(n - 1) Clifford pencils need rho(2(n-1)) >= n, true exactly for
     n = 2, 3, 5, 9 among n <= 40.

Writes ../results/state_dimension.json.  Runtime about 15 s.
"""

import json
import os
import random
import time

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------------------------- pencils

def arrow(n):
    """A(xi) = [[xi_n, xi'^T], [xi', -xi_n I]]: C_j = E_0j + E_j0 (j < n), C_n = diag(1, -I)."""
    Cs = []
    for j in range(1, n):
        M = sp.zeros(n, n)
        M[0, j] = M[j, 0] = 1
        Cs.append(M)
    M = -sp.eye(n)
    M[0, 0] = 1
    Cs.append(M)
    return Cs


def cayley_dickson_mult(k):
    """Left-multiplication matrices L_e (k x k, integer) for the basis e of R, C, H, O (k = 1, 2, 4, 8)."""
    def mul(a, b):
        if len(a) == 1:
            return [a[0] * b[0]]
        h = len(a) // 2
        a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
        conj = lambda x: [x[0]] + [-y for y in x[1:]]
        # (a1, a2)(b1, b2) = (a1 b1 - conj(b2) a2, b2 a1 + a2 conj(b1))
        p = [u - v for u, v in zip(mul(a1, b1), mul(conj(b2), a2))]
        q = [u + v for u, v in zip(mul(b2, a1), mul(a2, conj(b1)))]
        return p + q
    basis = [[1 if i == j else 0 for i in range(k)] for j in range(k)]
    Ls = []
    for e in basis:
        L = sp.zeros(k, k)
        for j, b in enumerate(basis):
            col = mul(e, b)
            for i in range(k):
                L[i, j] = col[i]
        Ls.append(L)
    return Ls


def hopf_pencil(k):
    """C(r, a) = [[r, L_a^T], [L_a, -r]] on R^{2k}: returns [C_r, C_{e_1}, ..., C_{e_k}] (n = k + 1)."""
    Ls = cayley_dickson_mult(k)
    Cr = sp.diag(*([1] * k + [-1] * k))
    Cs = [Cr]
    for L in Ls:
        M = sp.zeros(2 * k, 2 * k)
        M[:k, k:] = L.T
        M[k:, :k] = L
        Cs.append(M)
    return Cs


def flux_rank(Cs, z, tol=1e-9):
    N = len(z)
    cols = [z] + [np.array(C, dtype=float) @ z for C in Cs]
    M = np.stack(cols, axis=1)
    s = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(s > tol * s[0]))
    return rank, s


def random_complex(N, rng):
    return np.array([complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(N)])


# ----------------------------------------------------------------------------------------------- Part R

def part_R(rng):
    out = {}
    ok = True
    # arrow pencils: cone identity and rank n
    for n in (2, 3, 4):
        Cs = arrow(n)
        z = sp.symbols("z0:%d" % n)
        zv = sp.Matrix(z)
        q0 = (zv.T * zv)[0]
        qs = [(zv.T * C * zv)[0] for C in Cs]
        cone = sp.expand(sum(q * q for q in qs) - q0 ** 2) == 0
        ranks = [flux_rank(Cs, random_complex(n, rng))[0] for _ in range(5)]
        out["arrow n=%d" % n] = {"flux image on the light cone": bool(cone), "rank M(z)": ranks}
        ok &= cone and all(r == n for r in ranks)
    # random pencils with N = n + 1
    for n in (2, 3, 4):
        N = n + 1
        rs = random.Random(100 + n)
        Cs = []
        for _ in range(n):
            A = sp.Matrix(N, N, lambda i, j: 0)
            for i in range(N):
                for j in range(i, N):
                    A[i, j] = A[j, i] = rs.randint(-3, 3)
            Cs.append(A)
        z = sp.symbols("z0:%d" % N)
        zv = sp.Matrix(z)
        M = sp.Matrix.hstack(zv, *[C * zv for C in Cs])
        detM = sp.expand(M.det())
        nonzero = detM != 0
        # simple top eigenvalue at a sample direction
        xi = np.array([rs.uniform(-1, 1) for _ in range(n)])
        Cxi = sum(float(x) * np.array(C, dtype=float) for x, C in zip(xi, Cs))
        ev = np.linalg.eigvalsh(Cxi)
        simple = ev[-1] - ev[-2] > 1e-6
        ranks = [flux_rank(Cs, random_complex(N, rng))[0] for _ in range(5)]
        out["random N=n+1, n=%d" % n] = {"det M(z) nonzero polynomial": bool(nonzero), "top eigenvalue simple": bool(simple),
                                         "rank M(z)": ranks}
        ok &= nonzero and simple and all(r == n + 1 for r in ranks)
    # coupled arrow, n = 3, N = 4
    n, eps = 3, sp.Rational(1, 5)
    A = arrow(n)
    Cs = []
    couple = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    diag = [sp.Rational(1, 3), -sp.Rational(1, 4), sp.Rational(1, 6)]
    for j in range(n):
        M = sp.zeros(4, 4)
        M[:3, :3] = A[j]
        M[3, 3] = diag[j]
        for i in range(3):
            M[i, 3] = M[3, i] = eps * couple[j][i]
        Cs.append(M)
    vecs = []
    simple_all = True
    for _ in range(12):
        xi = np.array([rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3), 1.0])
        Cxi = sum(float(x) * np.array(C, dtype=float) for x, C in zip(xi, Cs))
        w, U = np.linalg.eigh(Cxi)
        simple_all &= w[-1] - w[-2] > 1e-3
        vecs.append(U[:, -1])
    span = int(np.linalg.matrix_rank(np.stack(vecs, axis=1), tol=1e-8))
    ranks = [flux_rank(Cs, random_complex(4, rng))[0] for _ in range(5)]
    out["coupled arrow n=3, N=4"] = {"top simple on the cap": bool(simple_all), "span of top eigenvectors": span,
                                     "rank M(z)": ranks}
    ok &= simple_all and span == 4 and all(r == 4 for r in ranks)
    # uncoupled direct sum arrow (+) scalar: N = 4 > n, but the top eigenvectors span only V = R^3 (N' = n);
    # compressed to V the pencil is the arrow pencil, rank n, so the engine exists after compression
    Cs0 = []
    for j in range(n):
        M = sp.zeros(4, 4)
        M[:3, :3] = A[j]
        M[3, 3] = diag[j]
        Cs0.append(M)
    vecs0 = []
    for _ in range(12):
        xi = np.array([rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3), 1.0])
        Cxi = sum(float(x) * np.array(C, dtype=float) for x, C in zip(xi, Cs0))
        w, U = np.linalg.eigh(Cxi)
        vecs0.append(U[:, -1])
    span0 = int(np.linalg.matrix_rank(np.stack(vecs0, axis=1), tol=1e-8))
    restricted = [flux_rank([C[:3, :3] for C in Cs0], random_complex(3, rng))[0] for _ in range(5)]
    out["arrow (+) scalar (uncoupled)"] = {"span of top eigenvectors N'": span0, "rank of the compressed pencil": restricted}
    ok &= span0 == 3 and all(r == 3 for r in restricted)
    # fibre identity z^T S z = 2 (phi Q0 - xi . Q')
    xis = sp.symbols("x1:4")
    ph = sp.Symbol("phi")
    z = sp.Matrix(sp.symbols("z0:3"))
    Cxi = sum((x * C for x, C in zip(xis, A)), sp.zeros(3, 3))
    S = ph * sp.eye(3) - Cxi
    lhs = sp.expand((z.T * S * z)[0])
    rhs = sp.expand(2 * (ph * (z.T * z)[0] / 2 - sum(x * (z.T * C * z)[0] / 2 for x, C in zip(xis, A))))
    out["fibre identity"] = bool(sp.expand(lhs - rhs) == 0)
    ok &= out["fibre identity"]
    return bool(ok), out


# ----------------------------------------------------------------------------------------------- Part H

def part_H(rng):
    mp.mp.dps = 40
    out = {}
    ok = True
    for k in (1, 2, 4, 8):
        Cs = hopf_pencil(k)
        n, N = k + 1, 2 * k
        res = {}
        anti = True
        for i in range(n):
            for j in range(n):
                P = Cs[i] * Cs[j] + Cs[j] * Cs[i]
                anti &= P == (2 * sp.eye(N) if i == j else sp.zeros(N, N))
        res["anticommutation"] = bool(anti)
        z = sp.Matrix(sp.symbols("z0:%d" % N))
        q0 = (z.T * z)[0]
        hopf = sp.expand(sum(((z.T * C * z)[0]) ** 2 for C in Cs) - q0 ** 2) == 0
        res["Hopf identity"] = bool(hopf)
        xi = np.array([rng.uniform(-1, 1) for _ in range(n)])
        Cxi = sum(float(x) * np.array(C, dtype=float) for x, C in zip(xi, Cs))
        ev = np.linalg.eigvalsh(Cxi)
        m = int(np.sum(np.abs(ev - np.linalg.norm(xi)) < 1e-9))
        res["multiplicity m of |xi|"] = m
        res["N = n + m - 1"] = bool(N == n + m - 1)
        ranks = [flux_rank(Cs, random_complex(N, rng))[0] for _ in range(4)]
        res["rank M(z)"] = ranks
        # passive coordinate and tangency at complex points near a real state
        Cm = [mp.matrix(np.array(C, dtype=float).tolist()) for C in Cs]
        z0 = [mp.mpf(rng.gauss(0, 1)) for _ in range(N)]
        nz = mp.sqrt(sum(x * x for x in z0))
        z0 = [x / nz for x in z0]
        Qp0 = [sum(z0[a] * sum(Cm[j][a, b] * z0[b] for b in range(N)) for a in range(N)) / 2 for j in range(n)]
        nrm = mp.sqrt(sum(x * x for x in Qp0))
        w = [x / nrm for x in Qp0]                       # passive direction: the base point's own direction
        worst_phi, worst_tan = mp.mpf(0), mp.mpf(0)
        for _ in range(4):
            zc = [z0[a] + mp.mpc(rng.uniform(-0.05, 0.05), rng.uniform(-0.05, 0.05)) for a in range(N)]
            eta = mp.mpc(rng.uniform(-0.05, 0.05), rng.uniform(-0.05, 0.05))
            Q0 = sum(x * x for x in zc) / 2
            Qp = [sum(zc[a] * sum(Cm[j][a, b] * zc[b] for b in range(N)) for a in range(N)) / 2 for j in range(n)]
            P = [Qp[j] + w[j] * eta for j in range(n)]
            wP = sum(w[j] * P[j] for j in range(n))
            Phi = wP - mp.sqrt(Q0 ** 2 - sum(x * x for x in P) + wP ** 2)
            worst_phi = max(worst_phi, abs(Phi - eta))
            # gradient of Phi at (tau, P) = (Q0, P), then tangency sum_mu dPhi_mu C^mu z = 0
            root = mp.sqrt(Q0 ** 2 - sum(x * x for x in P) + wP ** 2)
            dtau = -Q0 / root
            dP = [w[j] - (-P[j] + wP * w[j]) / root for j in range(n)]
            vec = [dtau * zc[a] + sum(dP[j] * sum(Cm[j][a, b] * zc[b] for b in range(N)) for j in range(n)) for a in range(N)]
            worst_tan = max(worst_tan, max(abs(x) for x in vec))
        res["Phi o Theta = eta, max error"] = float(worst_phi)
        res["tangency, max |sum dPhi C z|"] = float(worst_tan)
        # polarization = Hopf map
        zr = np.array([rng.gauss(0, 1) for _ in range(N)])
        h = np.array([zr @ np.array(C, dtype=float) @ zr for C in Cs]) / (zr @ zr)
        Ch = sum(x * np.array(C, dtype=float) for x, C in zip(h, Cs))
        res["|h(z)| = 1"] = float(abs(np.linalg.norm(h) - 1))
        res["C(h(z)) z = z, error"] = float(np.linalg.norm(Ch @ zr - zr) / np.linalg.norm(zr))
        good = (anti and hopf and m == k and N == n + m - 1 and all(r == n for r in ranks) and worst_phi < 1e-30
                and worst_tan < 1e-30 and res["|h(z)| = 1"] < 1e-12 and res["C(h(z)) z = z, error"] < 1e-12)
        res["ok"] = bool(good)
        out["k=%d (n=%d, N=%d)" % (k, n, N)] = res
        ok &= good
    # controls
    ctrl = {}
    Cs = hopf_pencil(4)[:4]                               # n = 4 directions on R^8: m = 4, n + m - 1 = 7 < 8
    z = sp.Matrix(sp.symbols("z0:8"))
    q0 = (z.T * z)[0]
    ctrl["H, four matrices: Hopf identity fails"] = bool(sp.expand(sum(((z.T * C * z)[0]) ** 2 for C in Cs) - q0 ** 2) != 0)
    ctrl["H, four matrices: rank M(z)"] = [flux_rank(Cs, random_complex(8, rng))[0] for _ in range(4)]
    ok &= ctrl["H, four matrices: Hopf identity fails"] and all(r == 5 for r in ctrl["H, four matrices: rank M(z)"])
    Cs = hopf_pencil(2)
    rs = random.Random(5)
    pert = []
    for C in Cs:
        D = sp.Matrix(4, 4, lambda i, j: 0)
        for i in range(4):
            for j in range(i, 4):
                D[i, j] = D[j, i] = sp.Rational(rs.randint(-3, 3), 20)
        pert.append(C + D)
    ctrl["perturbed C pencil: rank M(z)"] = [flux_rank(pert, random_complex(4, rng))[0] for _ in range(4)]
    ok &= all(r == 4 for r in ctrl["perturbed C pencil: rank M(z)"])
    A = arrow(3)
    dsum = [sp.diag(C, C) for C in A]                    # arrow (+) arrow: N = 6, m = 2, n + m - 1 = 4 < 6
    ctrl["arrow (+) arrow, one lift: rank M(z)"] = [flux_rank(dsum, random_complex(6, rng))[0] for _ in range(4)]
    ok &= all(r == 4 for r in ctrl["arrow (+) arrow, one lift: rank M(z)"])
    out["controls"] = ctrl
    return bool(ok), out


# ----------------------------------------------------------------------------------------------- Part Q

def realify_complex(H):
    A, B = H.real, H.imag
    return np.block([[A, -B], [B, A]])


def qmat(q):
    """4 x 4 real matrix of left multiplication by the quaternion q = (a, b, c, d)."""
    a, b, c, d = q
    return np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]], dtype=float)


def realify_quaternion(entries, mq):
    M = np.zeros((4 * mq, 4 * mq))
    for i in range(mq):
        for j in range(mq):
            M[4 * i:4 * i + 4, 4 * j:4 * j + 4] = qmat(entries[i][j])
    return M


def random_hermitian_complex(mc, rs):
    X = np.array([[complex(rs.gauss(0, 1), rs.gauss(0, 1)) for _ in range(mc)] for _ in range(mc)])
    return (X + X.conj().T) / 2


def random_hermitian_quaternion(mq, rs):
    ent = [[None] * mq for _ in range(mq)]
    for i in range(mq):
        for j in range(i, mq):
            if i == j:
                ent[i][j] = (rs.gauss(0, 1), 0.0, 0.0, 0.0)
            else:
                q = tuple(rs.gauss(0, 1) for _ in range(4))
                ent[i][j] = q
                ent[j][i] = (q[0], -q[1], -q[2], -q[3])
    return realify_quaternion(ent, mq)


def pencil_profile(Cs, rs, n_dirs=40):
    """Multiplicity and span of the top eigenspaces near a random direction, and transversal curvature there."""
    n = len(Cs)
    xi0 = np.array([rs.gauss(0, 1) for _ in range(n)])
    xi0 /= np.linalg.norm(xi0)

    def top(xi):
        w, U = np.linalg.eigh(sum(x * C for x, C in zip(xi, Cs)))
        return w, U

    w, U = top(xi0)
    m = int(np.sum(w > w[-1] - 1e-9))
    vecs = []
    for _ in range(n_dirs):
        xi = xi0 + 0.05 * np.array([rs.gauss(0, 1) for _ in range(n)])
        w2, U2 = top(xi)
        vecs.append(U2[:, -m:])
    span = int(np.linalg.matrix_rank(np.concatenate(vecs, axis=1), tol=1e-8))
    # transversal curvature: Hessian of phi = lambda_max by central differences, restricted to xi0^perp
    h = 1e-4
    phi = lambda xi: top(xi)[0][-1]
    Hs = np.zeros((n, n))
    E = np.eye(n)
    for i in range(n):
        for j in range(n):
            Hs[i, j] = (phi(xi0 + h * E[i] + h * E[j]) - phi(xi0 + h * E[i] - h * E[j])
                        - phi(xi0 - h * E[i] + h * E[j]) + phi(xi0 - h * E[i] - h * E[j])) / (4 * h * h)
    P = np.eye(n) - np.outer(xi0, xi0)
    ev = np.linalg.eigvalsh(P @ Hs @ P)
    curved = bool(np.sort(ev)[1] > 1e-5)          # n - 1 positive eigenvalues (one zero, along xi0)
    gap = float(w[-1] - w[-m - 1]) if m < len(w) else float("inf")
    return m, span, curved, gap


def part_Q(rng):
    out = {}
    ok = True
    rs = random.Random(31)
    cases = []
    for mc in (2, 3, 4):
        n = 2 * mc - 1
        Cs = [realify_complex(random_hermitian_complex(mc, rs)) for _ in range(n)]
        cases.append(("complex Hermitian %dx%d, n=%d" % (mc, mc, n), Cs, 2, True))
        Cs_ctrl = Cs[:-1]
        cases.append(("control: complex %dx%d, n=%d" % (mc, mc, n - 1), Cs_ctrl, 2, False))
    for mq in (2, 3):
        n = 4 * mq - 3
        Cs = [random_hermitian_quaternion(mq, rs) for _ in range(n)]
        cases.append(("quaternionic Hermitian %dx%d, n=%d" % (mq, mq, n), Cs, 4, True))
        cases.append(("control: quaternionic %dx%d, n=%d" % (mq, mq, n - 1), Cs[:-1], 4, False))
    # Maxwell: (E, B) in R^6, C_j = [[0, -eps_j], [eps_j, 0]] with (eps_j)_{ab} = eps_{jab}
    eps = np.zeros((3, 3, 3))
    for (a, b, c), sgn in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
        eps[a, b, c] = sgn
    Cm = [np.block([[np.zeros((3, 3)), -eps[j]], [eps[j], np.zeros((3, 3))]]) for j in range(3)]
    cases.append(("Maxwell (isotropic), n=3", Cm, 2, False))
    for label, Cs, m_expected, should_pass in cases:
        n = len(Cs)
        N = Cs[0].shape[0]
        m, span, curved, gap = pencil_profile(Cs, rs)
        ranks = [flux_rank(Cs, random_complex(N, rng))[0] for _ in range(4)]
        meets = span == n + m - 1
        passes = all(r == n for r in ranks)
        res = {"N": N, "multiplicity m": m, "span N'": span, "n + m - 1": n + m - 1, "transversally curved": curved,
               "rank M(z)": ranks, "passes tangency test": passes}
        good = (m == m_expected) and (passes == should_pass) and (meets == should_pass)
        if should_pass:
            good &= curved
        res["consistent with Theorem 17.3"] = bool(good)
        out[label] = res
        ok &= good
    return bool(ok), out


# ----------------------------------------------------------------------------------------------- Part W

def rho(N):
    """Hurwitz-Radon number: N = 2^(4a+b) * odd, rho = 8a + 2^b."""
    c = 0
    while N % 2 == 0:
        N //= 2
        c += 1
    a, b = divmod(c, 4)
    return 8 * a + 2 ** b


def part_W():
    good_n = [n for n in range(2, 41) if rho(2 * (n - 1)) >= n]
    ok = good_n == [2, 3, 5, 9]
    return ok, {"n <= 40 with rho(2(n-1)) >= n": good_n}


def main():
    rng = random.Random(17)
    t0 = time.time()
    results = {}
    ok_all = True
    for name, fn in (("R_rank_criterion", lambda: part_R(rng)), ("H_hopf_pencils", lambda: part_H(rng)),
                     ("Q_division_algebra_pencils", lambda: part_Q(rng)), ("W_hurwitz_radon", part_W)):
        ok, out = fn()
        results[name] = {"ok": bool(ok), "detail": out}
        print("%-22s %s" % (name, "ok" if ok else "FAIL"))
        ok_all &= ok
    results["all_ok"] = bool(ok_all)
    results["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "..", "results", "state_dimension.json"), "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print(json.dumps(results, indent=1, default=str)[:6000])
    print("ALL OK" if ok_all else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
