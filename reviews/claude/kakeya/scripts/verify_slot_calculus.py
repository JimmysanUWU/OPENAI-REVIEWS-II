"""Checks for 12_QUATERNIONIC_SLOT_CALCULUS.md (the slot calculus with arbitrary quaternionic frames).

A *frame* gives each slot a finite list of quaternions f_{a,1..k_a}; for a real function u of the variables
y_{a,i} put G_a(u) = sum_i f_{a,i} d_{a,i} u and H_f(u)_{ab} = sum_{i,j} f_{a,i} conj(f_{b,j}) d_{a,i} d_{b,j} u.

Part F  cofactor divergence for arbitrary frames (Theorem 12.1, new):
        sum_{a,i} conj(f_{a,i}) d_{a,i} C(H_f(u))_{ab} = 0, exact symbolic, for random rational frames of
        mixed sizes (1..5 vectors per slot, so including dependent frames), m = 2, 3 slots, random quartic u.
        Negative controls: the divergence taken with f instead of conj(f); a Hermitian field that is not a
        frame Hessian.
Part R  frame reduction (Lemma 12.2): H_f(u) o A = H^{(4)}(u o A), with A_{a, i alpha} = (f_{a,i})_alpha.
Part K  classical cases: one real unit gives the real Hessian and the real determinant; units (1, i) give
        H_ab = 4 d_{z-bar_a} d_{z_b} u and the Moore determinant equals the complex determinant.
Part D  the 3-slot decomposition (Proposition 12.5): for u(X, P), P = p + k q,
        H_ab = d_{X_a X_b} u + 4 dbar_a d_b u + 2 (dbar_a d_{X_b} u - d_{X_a} dbar_b u) i.
Part E  ellipticity (Proposition 12.4): the symbol g^* Q^{-1} g of Q^{-1} : H_f is positive definite iff every
        slot frame is R-linearly independent; a 5-vector frame has a null covector.
Part U  Moore determinant facts used by the energy (Theorem 12.3): Euler's identity C(A) : A = m det A, and
        C(A) > 0 for A > 0 (the sign that makes the Monge-Ampere energy concave).
Writes ../results/slot_calculus.json.
"""

import json
import os
import random
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_quaternion_calculus import (moore_det, qadd, qconj, qmul, qscale, rand_herm, random_poly,  # noqa: E402
                                        sym_cofactor, sym_moore_det, to_complex, herm_inner, cofactor_numeric)

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "slot_calculus.json")
UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]


def rand_quat(rng):
    return tuple(sp.Rational(rng.randint(-4, 4), rng.randint(1, 3)) for _ in range(4))


def frame_hessian(u, ys, frame):
    m = len(frame)
    H = [[None] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            acc = (0, 0, 0, 0)
            for i, fa in enumerate(frame[a]):
                for j, fb in enumerate(frame[b]):
                    c = sp.diff(u, ys[a][i], ys[b][j])
                    if c != 0:
                        acc = qadd(acc, qscale(c, qmul(fa, qconj(fb))))
            H[a][b] = tuple(sp.expand(x) for x in acc)
    return H


def frame_divergence(C, ys, frame, conj=True):
    m = len(frame)
    worst = 0
    for b in range(m):
        acc = (0, 0, 0, 0)
        for a in range(m):
            for i, fa in enumerate(frame[a]):
                dC = tuple(sp.diff(c, ys[a][i]) for c in C[a][b])
                acc = qadd(acc, qmul(qconj(fa) if conj else fa, dC))
        acc = tuple(sp.expand(x) for x in acc)
        worst = max(worst, sum(1 for x in acc if x != 0))
    return worst


def part_F(rng):
    cases = [("complex + 3-slot", [[UNITS[0], UNITS[1]], [UNITS[0], UNITS[1], UNITS[2]]]),
             ("random sizes (2, 3)", None), ("random sizes (5, 2), dependent", None),
             ("random sizes (1, 2, 4)", None), ("random sizes (3, 3, 3)", None)]
    sizes = {"random sizes (2, 3)": (2, 3), "random sizes (5, 2), dependent": (5, 2),
             "random sizes (1, 2, 4)": (1, 2, 4), "random sizes (3, 3, 3)": (3, 3, 3)}
    out = []
    ok = True
    for name, frame in cases:
        if frame is None:
            frame = [[rand_quat(rng) for _ in range(k)] for k in sizes[name]]
        m = len(frame)
        ys = [[sp.Symbol(f"y{a}_{i}", real=True) for i in range(len(frame[a]))] for a in range(m)]
        flat = [y for row in ys for y in row]
        u = random_poly(flat, 4, rng)
        H = frame_hessian(u, ys, frame)
        hermitian = all(sp.expand(H[a][b][c] - qconj(H[b][a])[c]) == 0 for a in range(m) for b in range(m) for c in range(4))
        C = sym_cofactor(H, m)
        div = frame_divergence(C, ys, frame)
        rec = {"case": name, "m": m, "frame_sizes": [len(f) for f in frame], "hermitian": hermitian,
               "nonzero_divergence_components": div}
        good = hermitian and div == 0
        if name in ("random sizes (2, 3)", "random sizes (1, 2, 4)"):
            rec["control_divergence_with_f"] = frame_divergence(C, ys, frame, conj=False)
            Hp = [row[:] for row in H]
            Hp[0][0] = (Hp[0][0][0] + flat[-1] ** 2 * flat[0], 0, 0, 0)
            rec["control_non_hessian_field"] = frame_divergence(sym_cofactor(Hp, m), ys, frame)
            good = good and rec["control_divergence_with_f"] > 0 and rec["control_non_hessian_field"] > 0
        rec["ok"] = bool(good)
        ok = ok and good
        out.append(rec)
    return out, ok


def part_R(rng):
    """H_f(u) o A = H^{(4)}(u o A) for a random frame (sizes 2, 3) and a random cubic u."""
    frame = [[rand_quat(rng) for _ in range(2)], [rand_quat(rng) for _ in range(3)]]
    m = 2
    ys = [[sp.Symbol(f"y{a}_{i}", real=True) for i in range(len(frame[a]))] for a in range(m)]
    xs = [[sp.Symbol(f"x{a}_{al}", real=True) for al in range(4)] for a in range(m)]
    flat = [y for row in ys for y in row]
    u = random_poly(flat, 3, rng)
    sub = {ys[a][i]: sum(frame[a][i][al] * xs[a][al] for al in range(4)) for a in range(m) for i in range(len(frame[a]))}
    U = sp.expand(u.subs(sub))
    H4 = frame_hessian(U, xs, [UNITS, UNITS])
    Hf = frame_hessian(u, ys, frame)
    ok = all(sp.expand(H4[a][b][c] - Hf[a][b][c].subs(sub)) == 0 for a in range(m) for b in range(m) for c in range(4))
    return {"frame_sizes": [2, 3], "identity_holds": bool(ok)}, bool(ok)


def part_K(rng):
    out = {}
    # one real unit: real Hessian
    xs = sp.symbols("x0:3", real=True)
    u = random_poly(list(xs), 4, rng)
    H = frame_hessian(u, [[x] for x in xs], [[UNITS[0]] for _ in xs])
    real_ok = all(sp.expand(H[a][b][0] - sp.diff(u, xs[a], xs[b])) == 0 and all(H[a][b][c] == 0 for c in (1, 2, 3))
                  for a in range(3) for b in range(3))
    detM = sym_moore_det(H, 3)
    detR = sp.expand(sp.Matrix(3, 3, lambda a, b: sp.diff(u, xs[a], xs[b])).det())
    out["real_hessian"] = bool(real_ok)
    out["real_determinant"] = bool(sp.expand(detM - detR) == 0)
    # units (1, i): H_ab = 4 d_{zbar_a} d_{z_b} u, Moore det = complex det
    ps = sp.symbols("p0:2", real=True)
    qs = sp.symbols("q0:2", real=True)
    u = random_poly(list(ps) + list(qs), 4, rng)
    H = frame_hessian(u, [[ps[a], qs[a]] for a in range(2)], [[UNITS[0], UNITS[1]] for _ in range(2)])
    dz = lambda f, a: (sp.diff(f, ps[a]) - sp.I * sp.diff(f, qs[a])) / 2
    dzb = lambda f, a: (sp.diff(f, ps[a]) + sp.I * sp.diff(f, qs[a])) / 2
    cx_ok = True
    Hc = sp.zeros(2, 2)
    for a in range(2):
        for b in range(2):
            val = sp.expand(4 * dzb(dz(u, b), a))
            Hc[a, b] = val
            cx_ok &= sp.expand(sp.re(val) - H[a][b][0]) == 0 and sp.expand(sp.im(val) - H[a][b][1]) == 0
            cx_ok &= H[a][b][2] == 0 and H[a][b][3] == 0
    out["complex_hessian"] = bool(cx_ok)
    out["complex_determinant"] = bool(sp.expand(sym_moore_det(H, 2) - sp.expand(Hc.det())) == 0)
    return out, all(out.values())


def part_D(rng):
    """3-slot decomposition for u(X, p, q), P = p + k q, Wirtinger dbar = (d_p + k d_q)/2, d = (d_p - k d_q)/2."""
    m = 2
    X = sp.symbols("X0:2", real=True)
    p = sp.symbols("p0:2", real=True)
    q = sp.symbols("q0:2", real=True)
    u = random_poly(list(X) + list(p) + list(q), 4, rng)
    H = frame_hessian(u, [[X[a], p[a], q[a]] for a in range(m)], [UNITS[:3] for _ in range(m)])

    def kc(re_, im_):        # C_k number re + k im as a quaternion
        return (re_, 0, 0, im_)

    def dbar(f, a):
        return (sp.diff(f, p[a]) / 2, sp.diff(f, q[a]) / 2)          # (real part, k part)

    def d(f, a):
        return (sp.diff(f, p[a]) / 2, -sp.diff(f, q[a]) / 2)

    def apply(op_pair_fn, pair, a):
        """apply a C_k-valued first-order operator to a C_k-valued function given as (re, im)."""
        r, i_ = pair
        A_, B_ = op_pair_fn(r, a)
        C_, D_ = op_pair_fn(i_, a)
        # (A + kB) applied to r + k i  = (A r) + k(B r) + k(A i) + k^2 (B i)
        return (A_ - D_, B_ + C_)

    ok = True
    for a in range(m):
        for b in range(m):
            real_part = (sp.diff(u, X[a], X[b]), 0, 0, 0)
            cpx = apply(dbar, d(u, b), a)                                  # dbar_a d_b u in C_k
            term2 = qscale(4, kc(*cpx))
            m1 = apply(dbar, (sp.diff(u, X[b]), 0), a)                     # dbar_a d_{X_b} u
            m2 = apply(lambda f, c: (sp.diff(f, X[c]), 0), dbar(u, b), a)  # d_{X_a} dbar_b u
            mixed = qmul(kc(2 * (m1[0] - m2[0]), 2 * (m1[1] - m2[1])), UNITS[1])
            rhs = qadd(qadd(real_part, term2), mixed)
            ok &= all(sp.expand(H[a][b][c] - rhs[c]) == 0 for c in range(4))
    return {"decomposition_holds": bool(ok)}, bool(ok)


def part_E(seed):
    rng = np.random.default_rng(seed)
    res = {}

    def symbol_min(frame_sizes, dependent_slot=None):
        m = len(frame_sizes)
        frames = [rng.normal(size=(k, 4)) for k in frame_sizes]
        A = rand_herm(m, rng, pd=True)
        Ainv_c = np.linalg.inv(to_complex(A))
        K = sum(frame_sizes)
        # real symbol matrix S_{(a,i),(b,j)} = Re( conj(f_{a,i}) Q^{-1}_{ab} f_{b,j} ), via the complex representation
        def qc(qv):
            return np.array([[qv[0] + 1j * qv[1], qv[2] + 1j * qv[3]], [-qv[2] + 1j * qv[3], qv[0] - 1j * qv[1]]])
        S = np.zeros((K, K))
        idx = [(a, i) for a in range(m) for i in range(frame_sizes[a])]
        for r, (a, i) in enumerate(idx):
            for c, (b, j) in enumerate(idx):
                blk = Ainv_c[np.ix_([a, a + m], [b, b + m])]
                val = qc(frames[a][i]).conj().T @ blk @ qc(frames[b][j])
                S[r, c] = np.real(np.trace(val)) / 2
        return float(np.linalg.eigvalsh((S + S.T) / 2).min())

    res["independent_frames_(3,4,2)_min_symbol"] = symbol_min((3, 4, 2))
    res["dependent_frame_(5,3)_min_symbol"] = symbol_min((5, 3))
    ok = res["independent_frames_(3,4,2)_min_symbol"] > 1e-6 and abs(res["dependent_frame_(5,3)_min_symbol"]) < 1e-10
    return res, bool(ok)


def part_U(seed):
    rng = np.random.default_rng(seed)
    worst_euler = 0.0
    min_cof_eig = np.inf
    for m in (2, 3, 4):
        for _ in range(4):
            A = rand_herm(m, rng, pd=True)
            C = cofactor_numeric(A)
            worst_euler = max(worst_euler, abs(herm_inner(C, A) - m * moore_det(A)) / abs(moore_det(A)))
            min_cof_eig = min(min_cof_eig, float(np.linalg.eigvalsh(to_complex(C)).min()))
    ok = worst_euler < 1e-6 and min_cof_eig > 0
    return {"euler_relative": worst_euler, "min_eigenvalue_of_C(A)_for_A>0": min_cof_eig}, bool(ok)


def main():
    rng = random.Random(1212)
    results = {}
    all_ok = True
    for name, fn in (("F", lambda: part_F(rng)), ("R", lambda: part_R(rng)), ("K", lambda: part_K(rng)),
                     ("D", lambda: part_D(rng)), ("E", lambda: part_E(7)), ("U", lambda: part_U(11))):
        res, ok = fn()
        results[name] = res
        all_ok = all_ok and ok
        print(name, "ok" if ok else "FAIL", json.dumps(res, default=str)[:400])
    results["all_ok"] = bool(all_ok)
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print("ALL OK" if all_ok else "FAILED")


if __name__ == "__main__":
    main()
