"""Checks for 21_GAUGE_QUOTIENT.md: the flux map of a pencil factors through a rank-one cone.

For the three division algebras F = R, C, H, a pencil of F-Hermitian matrices (realified to real symmetric
C^0 = I, C^1..C^n on R^N) commutes with the unit group K = {+-1}, U(1), Sp(1).  Claims checked:

  F  [Theorem 21.1, factorization]
     - R: Q^mu(z) = z^T C^mu z / 2 is linear in the rank-one symmetric matrix z z^T (Veronese cone, dim N).
     - C: on the +-k eigenspaces U_+, U_- of J, each C^mu pairs U_+ with U_- only (U_+- isotropic), so
       Q^mu = u^T C^mu v is linear in the rank-one matrix u v^T (Segre cone, dim 2m' - 1).
     - H: every C^mu commutes with I_1, I_2, I_3, so Q^mu lies in the space of K-invariant quadratic forms,
       which has dimension binom(2m'', 2) (quaternionic Hermitian forms); the invariant map
       pi(z) = (z^T A_k z)_k has generic rank 4m'' - 3 (the cone of rank-two skew matrices).
  D  [dimension count] generic rank of the quotient map pi equals N - dim K = N - (m - 1):
     N (R), 2m' - 1 (C), 4m'' - 3 (H).  So the tangency count N' = n + m - 1 says dim(quotient) = n.
  T  [tangency = generic finiteness] for pencils in the tangency class the flux map has generic rank exactly
     n = dim(quotient), i.e. the linear map L restricted to the cone is generically finite onto its image;
     with one direction removed (n - 1) the rank is n = dim(quotient) > n - 1, and tangency fails.
  N  controls: a generic real symmetric pencil with N = n + 1 has generic flux rank n + 1 (no gauge group, no
     drop); a random antisymmetric matrix does not commute with the Hermitian pencil.

Requires numpy.  Deterministic; writes ../results/gauge_quotient.json (about 2 s).
"""

import json
import os
from math import comb

import numpy as np

RS = np.random.default_rng(1101)
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


def rmat(q):
    return np.array([qmul(e, q) for e in np.eye(4)]).T


def realify_complex(H):
    return np.block([[H.real, -H.imag], [H.imag, H.real]])


def herm_c(m):
    X = RS.normal(size=(m, m)) + 1j * RS.normal(size=(m, m))
    return (X + X.conj().T) / 2


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


def jac_rank(f, N, trials=3, tol=1e-7):
    """Generic rank of the complex Jacobian of a polynomial map C^N -> C^K (finite differences, complex step)."""
    ranks = []
    for _ in range(trials):
        z = RS.normal(size=N) + 1j * RS.normal(size=N)
        h = 1e-6
        cols = [(f(z + h * e) - f(z - h * e)) / (2 * h) for e in np.eye(N)]
        Jm = np.array(cols).T
        s = np.linalg.svd(Jm, compute_uv=False)
        ranks.append(int(np.sum(s > tol * s[0])))
    return max(ranks)


def flux(Cs):
    return lambda z: np.array([0.5 * z @ C @ z for C in Cs])


def part_R():
    out = []
    for n in (2, 3, 4):
        # a pencil with N = n: arrow pencil plus a small random symmetric perturbation keeps the class generic
        Cs = [np.eye(n)] + [(lambda A: (A + A.T) / 2)(RS.normal(size=(n, n))) for _ in range(n)]
        pi = lambda z: np.array([z[a] * z[b] for a in range(n) for b in range(a, n)])
        r_pi = jac_rank(pi, n)
        r_Q = jac_rank(flux(Cs), n)
        check(r_pi == n, "Veronese cone has dimension N")
        # factorization: Q^mu = <C^mu, z z^T>/2 is linear in z z^T
        z = RS.normal(size=n) + 1j * RS.normal(size=n)
        check(np.allclose(flux(Cs)(z), [0.5 * np.sum(C * np.outer(z, z)) for C in Cs]), "R factorization")
        check(r_Q == n, "R: generic flux rank = n = dim quotient")
        out.append(dict(n=n, N=n, quotient_dim=r_pi, flux_rank=r_Q))
    RESULTS["R"] = out


def part_C():
    out = []
    for m in (2, 3, 4):
        n = 2 * m - 1
        Cs = [np.eye(2 * m)] + [realify_complex(herm_c(m)) for _ in range(n)]
        J = np.block([[np.zeros((m, m)), -np.eye(m)], [np.eye(m), np.zeros((m, m))]])
        for C in Cs:
            check(np.allclose(J @ C, C @ J), "J commutes")
        # eigenspaces of J over C_k (here 1j plays k): J e = +1j e on U_+, -1j e on U_-
        w, V = np.linalg.eig(J)
        Up = V[:, np.isclose(w, 1j)]
        Um = V[:, np.isclose(w, -1j)]
        check(Up.shape[1] == m and Um.shape[1] == m, "eigenspace dimensions")
        for C in Cs:
            check(np.allclose(Up.T @ C @ Up, 0) and np.allclose(Um.T @ C @ Um, 0), "U_+- isotropic for C^mu")
        # factorization through M = u v^T
        Lmu = [Up.T @ C @ Um for C in Cs]
        Bfull = np.concatenate([Up, Um], axis=1)
        for _ in range(5):
            u = RS.normal(size=m) + 1j * RS.normal(size=m)
            v = RS.normal(size=m) + 1j * RS.normal(size=m)
            z = Up @ u + Um @ v
            check(np.allclose(flux(Cs)(z), [u @ L @ v for L in Lmu]), "Q = u^T L v")
        Binv = np.linalg.inv(Bfull)

        def pi(z):
            c = Binv @ z
            return np.outer(c[:m], c[m:]).ravel()
        r_pi = jac_rank(pi, 2 * m)
        r_Q = jac_rank(flux(Cs), 2 * m)
        check(r_pi == 2 * m - 1, "Segre cone dimension 2m'-1")
        check(r_Q == n, "C: generic flux rank = n = dim quotient (tangency class)")
        # control: one direction fewer (n - 1): quotient dim n exceeds n - 1, flux rank n
        r_Q_ctrl = jac_rank(flux(Cs[:-1]), 2 * m)
        check(r_Q_ctrl == n, "control: n-1 directions, flux rank n > n-1 (tangency fails)")
        out.append(dict(mprime=m, n=n, N=2 * m, quotient_dim=r_pi, flux_rank=r_Q, flux_rank_one_fewer=r_Q_ctrl))
    RESULTS["C"] = out


def part_H():
    out = []
    units = [np.array(u, float) for u in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1])]
    for m in (2, 3):
        n = 4 * m - 3
        N = 4 * m
        Cs = [np.eye(N)] + [herm_h(m) for _ in range(n)]
        Is = [np.kron(np.eye(m), rmat(u)) for u in units]
        # basis of K-invariant quadratic forms: symmetric A commuting with I_1, I_2, I_3
        syms = []
        for a in range(N):
            for b in range(a, N):
                E = np.zeros((N, N))
                E[a, b] = E[b, a] = 1
                syms.append(E)
        rows = []
        for E in syms:
            rows.append(np.concatenate([(E @ I_ - I_ @ E).ravel() for I_ in Is]))
        Amat = np.array(rows).T          # columns: symmetric basis; kernel = invariant forms
        _, s, Vt = np.linalg.svd(Amat)
        null = Vt[np.sum(s > 1e-9):]
        dim_inv = null.shape[0]
        check(dim_inv == comb(2 * m, 2), "invariant quadratic forms: binom(2m'',2)")
        inv_forms = [sum(c * E for c, E in zip(vec, syms)) for vec in null]
        # each C^mu lies in their span
        span = np.array([F.ravel() for F in inv_forms]).T
        for C in Cs:
            coef, res, *_ = np.linalg.lstsq(span, C.ravel(), rcond=None)
            check(np.allclose(span @ coef, C.ravel(), atol=1e-8), "C^mu is a K-invariant form")
        pi = lambda z: np.array([z @ F @ z for F in inv_forms])
        r_pi = jac_rank(pi, N)
        r_Q = jac_rank(flux(Cs), N)
        check(r_pi == 4 * m - 3, "rank-two skew cone dimension 4m''-3")
        check(r_Q == n, "H: generic flux rank = n = dim quotient (tangency class)")
        r_Q_ctrl = jac_rank(flux(Cs[:-1]), N)
        check(r_Q_ctrl == n, "control: n-1 directions, flux rank n")
        out.append(dict(mpp=m, n=n, N=N, invariant_forms=dim_inv, quotient_dim=r_pi, flux_rank=r_Q,
                        flux_rank_one_fewer=r_Q_ctrl))
    RESULTS["H"] = out


def part_N():
    # a generic real symmetric pencil with N = n + 1: no gauge group, flux rank n + 1
    n = 3
    Cs = [np.eye(n + 1)] + [(lambda A: (A + A.T) / 2)(RS.normal(size=(n + 1, n + 1))) for _ in range(n)]
    r = jac_rank(flux(Cs), n + 1)
    check(r == n + 1, "control: generic N = n+1 real pencil has flux rank n+1")
    m = 3
    Hs = [realify_complex(herm_c(m)) for _ in range(5)]
    A = RS.normal(size=(2 * m, 2 * m))
    A = A - A.T
    check(not all(np.allclose(A @ C, C @ A) for C in Hs), "control: random antisymmetric does not commute")
    RESULTS["N"] = dict(real_N_eq_n_plus_1_flux_rank=r)


def main():
    part_R()
    part_C()
    part_H()
    part_N()
    print("R: Veronese cone, dim N = n, flux rank n: OK")
    print("C: Segre cone u v^T, dim 2m'-1 = n, Q = u^T L v, flux rank n (n-1 directions: tangency fails): OK")
    print("H: invariant forms binom(2m'',2), rank-two skew cone dim 4m''-3 = n, flux rank n: OK")
    print("Controls: generic N = n+1 real pencil has flux rank n+1; fake gauge does not commute: OK")
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "results", "gauge_quotient.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1)
    print("ALL OK")


if __name__ == "__main__":
    main()
