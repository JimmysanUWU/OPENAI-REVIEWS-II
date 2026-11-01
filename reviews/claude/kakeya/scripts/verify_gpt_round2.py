"""Independent re-check of GPT's Kakeya 22 (curvature design, universal fibre divergence) and Kakeya 23
(global rigidity), with methods different from GPT's scripts: floating-point spectral computations, Hessians by
finite differences, and coarea densities from the full complex Jacobian with a random complement.

Pencils: Hopf pencils over C, H, O (Weyl n=3, n=5, n=9), a random complex Hermitian 3x3 pencil in 5 variables
(m=2) and a random quaternionic Hermitian 2x2 pencil in 5 variables (m=4); all have N = n + m - 1.

  K  [Kakeya 22 Lemma 22.1] R_a S^-1 R_a^T = (a.a) H / 2 for random real and complex a in E, with the Hessian H
     of the top eigenvalue computed by finite differences (independent of the identity); Corollary 22.2:
     L(a) = sqrt2 H^-1/2 R_a S^-1/2 satisfies L(a) L(a)^T = |a|^2 I.
  J  [Kakeya 22 Lemma 22.3, Theorem 22.4] at random points of the complex quadric {a in E_C, a.a = 1}, the
     coarea density in the coordinates a' = (a_2..a_m), computed from the full Jacobian of
     Theta(z, eta) = (Q_0, ..., Q_{n-1}, Q_n + eta) with a random complement, equals gamma/|a_1|^2 with
     gamma = 1/(det S det(H/2)); the plane-wave time entry |a_1|^2 makes the product the constant gamma.
     gamma = 1 for the three Hopf pencils (GPT's octonionic value included).
  G  [Kakeya 23] Hopf pencils satisfy the spin-factor relations D_i D_j + D_j D_i = 2 delta_ij I (trace-free);
     a random non-Clifford complex Hermitian 3x3 pencil violates GLOBAL constant multiplicity: a direction where
     the top two eigenvalues of H(xi) meet is found numerically, as Theorem 23.2 requires.
  R  [round-10 route log] real commutants of the Hopf pencils have dimensions 1, 2, 4, 1 (k = 1, 2, 4, 8).

Requires numpy and scipy.  Deterministic; writes ../results/gpt_round2.json (about 10 s).
"""

import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_state_dimension import hopf_pencil  # noqa: E402

RS = np.random.default_rng(2211)
RESULTS = {}


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)


def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def lmat(q):
    return np.array([qmul(q, e) for e in np.eye(4)]).T


def herm_c(m):
    X = RS.normal(size=(m, m)) + 1j * RS.normal(size=(m, m))
    H = (X + X.conj().T) / 2
    return np.block([[H.real, -H.imag], [H.imag, H.real]])


def herm_h(m):
    M = np.zeros((4 * m, 4 * m))
    for i in range(m):
        for j in range(i, m):
            if i == j:
                M[4 * i:4 * i + 4, 4 * i:4 * i + 4] = lmat(np.array([RS.normal(), 0, 0, 0]))
            else:
                q = RS.normal(size=4)
                M[4 * i:4 * i + 4, 4 * j:4 * j + 4] = lmat(q)
                M[4 * j:4 * j + 4, 4 * i:4 * i + 4] = lmat(q * np.array([1, -1, -1, -1]))
    return M


def pencils():
    out = []
    for k, name in ((2, "Weyl (C Hopf, n=3)"), (4, "H Hopf (n=5)"), (8, "O Hopf (n=9)")):
        Cs = [np.array(C, dtype=float) for C in hopf_pencil(k)]
        out.append((name, Cs, k, True))
    out.append(("random complex Hermitian 3x3, n=5", [herm_c(3) for _ in range(5)], 2, False))
    out.append(("random quaternionic Hermitian 2x2, n=5", [herm_h(2) for _ in range(5)], 4, False))
    return out


def top(Cs, xi):
    w, U = np.linalg.eigh(sum(x * C for x, C in zip(xi, Cs)))
    return w, U


def phi(Cs, xi):
    return top(Cs, xi)[0][-1]


def setup(Cs, m):
    """Normalize xi0 = e_n; return E basis, F basis, S, d, H (finite differences)."""
    n = len(Cs)
    xi0 = np.zeros(n)
    xi0[-1] = 1
    w, U = top(Cs, xi0)
    check(abs(w[-1] - w[-m]) < 1e-9 and w[-1] - w[-m - 1] > 1e-6, "top multiplicity m at xi0")
    E = U[:, -m:]
    Fb = U[:, :-m]
    phi0 = w[-1]
    S = Fb.T @ (phi0 * np.eye(len(w)) - Cs[-1]) @ Fb
    d = np.array([E[:, 0] @ C @ E[:, 0] for C in Cs])
    for C in Cs:      # Hellmann-Feynman: P C P = d P
        Pd = E.T @ C @ E
        check(np.allclose(Pd, Pd[0, 0] * np.eye(m), atol=1e-9), "P C^j P is scalar on E")
    k = n - 1

    def fd(h):
        H = np.zeros((k, k))
        for i in range(k):
            for j in range(k):
                def f(si, sj):
                    x = xi0.copy()
                    x[i] += si * h
                    x[j] += sj * h
                    return phi(Cs, x)
                H[i, j] = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) / (4 * h * h)
        return H
    H = (4 * fd(1e-3) - fd(2e-3)) / 3          # Richardson extrapolation: error O(h^4)
    return E, Fb, S, d, H, phi0


def R_of(Cs, E, Fb, d, a):
    k = len(Cs) - 1
    vec = E @ a                       # a in coordinates of E (complex allowed)
    return np.array([Fb.T @ ((Cs[j] - d[j] * np.eye(len(vec))) @ vec) for j in range(k)])


def part_K_J():
    rows = []
    for name, Cs, m, is_hopf in pencils():
        n = len(Cs)
        N = Cs[0].shape[0]
        check(N == n + m - 1, "critical state dimension N = n + m - 1")
        E, Fb, S, d, H, phi0 = setup(Cs, m)
        Sinv = np.linalg.inv(S)
        errs = []
        for complex_a in (False, True):
            for _ in range(4):
                a = RS.normal(size=m) + (1j * RS.normal(size=m) if complex_a else 0)
                R = R_of(Cs, E, Fb, d, a)
                lhs = R @ Sinv @ R.T
                rhs = 0.5 * (a @ a) * H
                errs.append(np.max(np.abs(lhs - rhs)) / np.max(np.abs(rhs)))
        check(max(errs) < 1e-7, "Lemma 22.1 (curvature identity) against finite-difference Hessian")
        # Corollary 22.2: orthogonal design
        wH, VH = np.linalg.eigh(H)
        Hm12 = VH @ np.diag(wH ** -0.5) @ VH.T
        wS, VS = np.linalg.eigh(S)
        Sm12 = VS @ np.diag(wS ** -0.5) @ VS.T
        for _ in range(3):
            a = RS.normal(size=m)
            L = np.sqrt(2) * Hm12 @ R_of(Cs, E, Fb, d, a) @ Sm12
            check(np.allclose(L @ L.T, (a @ a) * np.eye(n - 1), atol=1e-6 * (a @ a)), "Corollary 22.2: L L^T = |a|^2 I")
        gamma = 1 / (np.linalg.det(S) * np.linalg.det(H / 2))
        if is_hopf:
            check(abs(gamma - 1) < 1e-8, "gamma = 1 for Hopf pencils")
        # Lemma 22.3: coarea density from the full Jacobian with a random complement
        cdir = np.zeros(n + 1, dtype=complex)
        cdir[-1] = 1
        Call = [np.eye(N)] + Cs
        consts = []
        for _ in range(4):
            ap = (RS.normal(size=m - 1) + 1j * RS.normal(size=m - 1)) * 0.7
            a1 = np.sqrt(1 - ap @ ap + 0j)
            a = np.concatenate([[a1], ap])
            z = E @ a
            # tangent vectors of the quadric in coordinates a': dz/da'_alpha = E (e_alpha - (a'_alpha/a1) e_1)
            tangents = []
            for al in range(m - 1):
                da = np.zeros(m, dtype=complex)
                da[al + 1] = 1
                da[0] = -ap[al] / a1
                tangents.append(np.concatenate([E @ da, [0]]))
            T = np.array(tangents).T
            Mt = np.array([C @ z for C in Call])
            DTheta = np.concatenate([Mt, cdir[:, None]], axis=1)
            Ecomp = RS.normal(size=(N + 1, n + 1)) + 1j * RS.normal(size=(N + 1, n + 1))
            num = abs(np.linalg.det(np.concatenate([T, Ecomp], axis=1))) ** 2
            den = abs(np.linalg.det(DTheta @ Ecomp)) ** 2
            dens = num / den
            # Theorem 22.4's mechanism: time entry |z . e_1|^2 = |a_1|^2 times the density is constant
            check(abs(abs(z @ E[:, 0]) ** 2 - abs(a1) ** 2) < 1e-12 * (1 + abs(a1) ** 2), "time entry is |a_1|^2")
            consts.append(abs(z @ E[:, 0]) ** 2 * dens)
        consts = np.array(consts)
        flat = float(np.max(np.abs(consts / consts[0] - 1)))
        check(flat < 1e-9, "flat measure: |a_1|^2 x coarea density is constant on the quadric (exact algebra)")
        gap_gamma = float(abs(consts[0] / gamma - 1))
        check(gap_gamma < 1e-4, "Lemma 22.3: the constant equals gamma = 1/(det S det(H/2)) (H by finite differences)")
        rows.append(dict(pencil=name, n=n, N=N, m=m, gamma_from_hessian=float(gamma), gamma_from_coarea=float(consts[0]),
                         curvature_identity_err=float(max(errs)), flatness_err=flat, gamma_relative_gap=gap_gamma))
        print("   %-42s n=%d N=%2d m=%d  gamma=%.6g  identity err %.1e  flatness %.1e  gamma gap %.1e"
              % (name, n, N, m, consts[0], max(errs), flat, gap_gamma))
    RESULTS["K_J"] = rows
    print("Parts K, J: curvature identity, orthogonal design, coarea density gamma/|a_1|^2 and flat measure: OK")


def part_G():
    out = {}
    for k in (2, 4, 8):
        Cs = [np.array(C, dtype=float) for C in hopf_pencil(k)]
        N = Cs[0].shape[0]
        for i, Ci in enumerate(Cs):
            check(abs(np.trace(Ci)) < 1e-12, "trace-free: alpha = 0")
            for j, Cj in enumerate(Cs):
                check(np.allclose(Ci @ Cj + Cj @ Ci, 2 * (i == j) * np.eye(N)), "spin-factor relations")
    # non-Clifford pencil: find a crossing of the top two eigenvalues (global multiplicity fails)
    m = 3
    Hs = []
    for _ in range(5):
        X = RS.normal(size=(m, m)) + 1j * RS.normal(size=(m, m))
        Hs.append((X + X.conj().T) / 2)

    def gap(x):
        xi = x / np.linalg.norm(x)
        w = np.linalg.eigvalsh(sum(c * H for c, H in zip(xi, Hs)))
        return w[-1] - w[-2]
    best = min((minimize(gap, RS.normal(size=5), method="Nelder-Mead",
                         options=dict(xatol=1e-10, fatol=1e-12, maxiter=20000, maxfev=40000)) for _ in range(12)),
               key=lambda r: r.fun)
    xi = best.x / np.linalg.norm(best.x)
    w = np.linalg.eigvalsh(sum(c * H for c, H in zip(xi, Hs)))
    scale = np.max(np.abs(w))
    check((w[-1] - w[-2]) / scale < 1e-4, "a top-eigenvalue crossing exists for a non-Clifford pencil")
    out["noncliford_crossing_relative_gap"] = float((w[-1] - w[-2]) / scale)
    RESULTS["G"] = out
    print("Part G: Hopf pencils are spin factors; a random non-Clifford Hermitian 3x3 pencil has a top-eigenvalue "
          "crossing (relative gap %.1e), so the global hypotheses of Kakeya 23 fail for it: OK"
          % out["noncliford_crossing_relative_gap"])


def part_R():
    dims = {}
    for k in (1, 2, 4, 8):
        Cs = [np.array(C, dtype=float) for C in hopf_pencil(k)]
        N = Cs[0].shape[0]
        rows = []
        for C in Cs:
            # X C - C X = 0, vectorized: (C^T kron I - I kron C) vec(X) = 0
            rows.append(np.kron(C.T, np.eye(N)) - np.kron(np.eye(N), C))
        A = np.concatenate(rows, axis=0)
        s = np.linalg.svd(A, compute_uv=False)
        dims[k] = int(N * N - np.sum(s > 1e-9))
    check([dims[k] for k in (1, 2, 4, 8)] == [1, 2, 4, 1], "commutant dimensions 1, 2, 4, 1")
    RESULTS["R"] = dims
    print("Part R: real commutants of the Hopf pencils have dimensions 1, 2, 4, 1: OK")


def main():
    part_K_J()
    part_G()
    part_R()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "gpt_round2.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1)
    print("ALL OK")


if __name__ == "__main__":
    main()
