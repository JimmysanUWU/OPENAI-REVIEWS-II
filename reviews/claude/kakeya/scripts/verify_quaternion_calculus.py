"""Checks of Section 3 (quaternionic Hessian calculus) and Lemma 4.3.

Quaternions are 4-vectors (w, x, y, z) = w + xi + yj + zk; e_0..e_3 = 1, i, j, k.

  L3.1   H(xi xi^t) = g g^*, tr H(u) = Laplacian, product rule (3.2), and (3.4):
         H(sum |F_nu|^2) = 4 J^* J for F holomorphic in P in C_k^m (coordinates x_a + P_a i);
  L3.2   Moore determinant of Hermitian quaternionic matrices: (3.6) C(A) = det(A) A^{-1},
         (3.7) v^* C(A) v = det(A|v_perp) and det(A + t v v^*) = det A + t v^* C(A) v (also for
         indefinite A), (3.8) Schur formula, log det concavity along random segments;
  L3.3   cofactor divergence (3.12): sum_{a,alpha} conj(e_alpha) d_{a alpha} C(H(u))_{ab} = 0 for random
         polynomials u, in 3 slots (m = 2, 3) and 4 slots (m = 2); negative control with a
         Hermitian field that is not a Hessian;
  L4.3   det((1-theta) v v^* + theta D) = theta^n det D + (1-theta) theta^(n-1) v^* C(D) v and the
         resulting energy weights n/(n+1), 1/(n(n+1)).
"""

import json
import random
from itertools import product

import numpy as np
import sympy as sp

# ---------------------------------------------------------------- quaternion arithmetic (generic)


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


E = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
ZERO = (0, 0, 0, 0)

# ---------------------------------------------------------------- numeric Hermitian quaternionic matrices


def to_complex(A):
    """A: m x m array of quaternions (m, m, 4). q = a + b j, a = w + x i, b = y + z i."""
    m = A.shape[0]
    a = A[:, :, 0] + 1j * A[:, :, 1]
    b = A[:, :, 2] + 1j * A[:, :, 3]
    return np.block([[a, b], [-np.conj(b), np.conj(a)]])


def moore_det(A):
    ev = np.sort(np.linalg.eigvalsh(to_complex(A)))
    return float(np.prod(ev[::2]))


def rand_herm(m, rng, pd=False):
    A = np.zeros((m, m, 4))
    for a in range(m):
        A[a, a, 0] = rng.normal()
        for b in range(a + 1, m):
            q = rng.normal(size=4)
            A[a, b] = q
            A[b, a] = np.array(qconj(tuple(q)))
    if pd:
        Cm = to_complex(A)
        shift = -np.linalg.eigvalsh(Cm).min() + 0.5
        for a in range(m):
            A[a, a, 0] += shift
    return A


def qmat_mul(A, B):
    m, k, n = A.shape[0], A.shape[1], B.shape[1]
    C = np.zeros((m, n, 4))
    for i in range(m):
        for j in range(n):
            acc = (0.0, 0.0, 0.0, 0.0)
            for l in range(k):
                acc = qadd(acc, qmul(tuple(A[i, l]), tuple(B[l, j])))
            C[i, j] = acc
    return C


def qmat_star(A):
    m, n = A.shape[0], A.shape[1]
    S = np.zeros((n, m, 4))
    for i in range(m):
        for j in range(n):
            S[j, i] = qconj(tuple(A[i, j]))
    return S


def from_complex(Cm, m):
    a = Cm[:m, :m]
    b = Cm[:m, m:]
    A = np.zeros((m, m, 4))
    A[:, :, 0], A[:, :, 1], A[:, :, 2], A[:, :, 3] = a.real, a.imag, b.real, b.imag
    return A


def herm_inner(A, B):
    """A : B = Re tr(AB)."""
    m = A.shape[0]
    return sum(qmul(tuple(A[i, j]), tuple(B[j, i]))[0] for i in range(m) for j in range(m))


def cofactor_numeric(A, eps=1e-6):
    """C(A) via C(A):B = d/dt det(A+tB) (central differences over a Hermitian basis)."""
    m = A.shape[0]
    C = np.zeros((m, m, 4))
    for a in range(m):
        B = np.zeros((m, m, 4)); B[a, a, 0] = 1
        C[a, a, 0] = (moore_det(A + eps * B) - moore_det(A - eps * B)) / (2 * eps)
        for b in range(a + 1, m):
            for al in range(4):
                B = np.zeros((m, m, 4))
                B[a, b, al] = 1
                B[b, a] = np.array(qconj(tuple(B[a, b])))
                der = (moore_det(A + eps * B) - moore_det(A - eps * B)) / (2 * eps)
                C[a, b, al] = der / 2
            C[b, a] = np.array(qconj(tuple(C[a, b])))
    return C


def quat_orthonormal_complement(v):
    """Quaternionic orthonormal basis with first vector v/|v| (Gram-Schmidt, right scalars)."""
    m = v.shape[0]
    basis = [v / np.sqrt(np.sum(v ** 2))]
    for k in range(m):
        e = np.zeros((m, 4)); e[k, 0] = 1
        w = e.copy()
        for b in basis:
            # coefficient c = b^* w (quaternion), subtract b c
            c = (0.0, 0.0, 0.0, 0.0)
            for i in range(m):
                c = qadd(c, qmul(qconj(tuple(b[i])), tuple(w[i])))
            for i in range(m):
                w[i] = np.array(qadd(tuple(w[i]), qscale(-1, qmul(tuple(b[i]), c))))
        nrm = np.sqrt(np.sum(w ** 2))
        if nrm > 1e-8 and len(basis) < m:
            basis.append(w / nrm)
    U = np.zeros((m, m, 4))
    for j, b in enumerate(basis):
        U[:, j] = b
    return U


def moore_checks(rng):
    res = {"(3.6)": True, "(3.7)_compression": True, "(3.7)_rank_one_linear": True, "(3.8)_schur": True,
           "logdet_concave": True}
    for m in (2, 3, 4):
        for trial in range(6):
            A = rand_herm(m, rng, pd=True)
            C = cofactor_numeric(A)
            # (3.6): C = det(A) A^{-1}
            Ainv = from_complex(np.linalg.inv(to_complex(A)), m)
            res["(3.6)"] &= np.allclose(C, moore_det(A) * Ainv, atol=1e-5 * max(1, abs(moore_det(A))))
            # (3.7) for a random v and indefinite A too
            for Atest in (A, rand_herm(m, rng)):
                Ct = cofactor_numeric(Atest)
                v = rng.normal(size=(m, 4)); v /= np.sqrt(np.sum(v ** 2))
                vCv = sum(qmul(qconj(tuple(v[i])), qmul(tuple(Ct[i, j]), tuple(v[j])))[0] for i in range(m) for j in range(m))
                U = quat_orthonormal_complement(v)
                Ac = qmat_mul(qmat_star(U), qmat_mul(Atest, U))[1:, 1:]
                dcomp = moore_det(Ac) if m > 1 else 1.0
                res["(3.7)_compression"] &= abs(vCv - dcomp) < 1e-5 * max(1, abs(dcomp))
                vv = np.zeros((m, m, 4))
                for i in range(m):
                    for j in range(m):
                        vv[i, j] = qmul(tuple(v[i]), qconj(tuple(v[j])))
                for t in (-1.3, 0.4, 2.0):
                    lhs = moore_det(Atest + t * vv)
                    res["(3.7)_rank_one_linear"] &= abs(lhs - (moore_det(Atest) + t * vCv)) < 1e-5 * max(1, abs(lhs))
            # (3.8) Schur
            Dm = A[1:, 1:]
            c = A[1:, 0:1]
            Dinv = from_complex(np.linalg.inv(to_complex(Dm)), m - 1) if m > 1 else None
            cDc = qmat_mul(qmat_star(c), qmat_mul(Dinv, c))[0, 0, 0]
            sigma = A[0, 0, 0] - cDc
            res["(3.8)_schur"] &= abs(moore_det(A) - moore_det(Dm) * sigma) < 1e-6 * max(1, abs(moore_det(A)))
            # concavity of log det along a segment
            B = rand_herm(m, rng, pd=True)
            vals = [np.log(moore_det((1 - s) * A + s * B)) for s in np.linspace(0, 1, 9)]
            res["logdet_concave"] &= all(vals[i - 1] + vals[i + 1] - 2 * vals[i] <= 1e-9 for i in range(1, 8))
    return res


def segment_identity(rng):
    ok = True
    for m in (2, 3, 4):
        for _ in range(5):
            Dm = rand_herm(m, rng, pd=True)
            v = rng.normal(size=(m, 4)) * 3
            vv = np.zeros((m, m, 4))
            for i in range(m):
                for j in range(m):
                    vv[i, j] = qmul(tuple(v[i]), qconj(tuple(v[j])))
            C = cofactor_numeric(Dm)
            vCv = herm_inner(C, vv)
            for th in (0.1, 0.5, 0.9):
                lhs = moore_det((1 - th) * vv + th * Dm)
                rhs = th ** m * moore_det(Dm) + (1 - th) * th ** (m - 1) * vCv
                ok &= abs(lhs - rhs) < 1e-5 * max(1, abs(lhs))
            # integral weights
            xs = np.linspace(0, 1, 2001)
            vals = [moore_det((1 - th) * vv + th * Dm) for th in xs]
            integral = np.trapezoid(vals, xs)
            ok &= abs(integral - (moore_det(Dm) / (m + 1) + vCv / (m * (m + 1)))) < 1e-3 * max(1, abs(integral))
    return ok

# ---------------------------------------------------------------- symbolic Hessian calculus


def sym_hessian(u, xs, m, k):
    """H_ab(u) = sum_{alpha,beta} e_alpha conj(e_beta) d_{a alpha} d_{b beta} u, xs[a][alpha]."""
    H = [[None] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            acc = (0, 0, 0, 0)
            for al in range(k):
                for be in range(k):
                    coef = sp.diff(u, xs[a][al], xs[b][be])
                    if coef != 0:
                        acc = qadd(acc, qscale(coef, qmul(E[al], qconj(E[be]))))
            H[a][b] = tuple(sp.expand(x) for x in acc)
    return H


def sym_moore_det(A, m):
    if m == 1:
        return A[0][0][0]
    if m == 2:
        return sp.expand(A[0][0][0] * A[1][1][0] - sum(x ** 2 for x in A[0][1]))
    if m == 3:
        t = qmul(A[0][1], qmul(A[1][2], A[2][0]))[0]
        return sp.expand(A[0][0][0] * A[1][1][0] * A[2][2][0] + 2 * t
                         - A[0][0][0] * sum(x ** 2 for x in A[1][2])
                         - A[1][1][0] * sum(x ** 2 for x in A[0][2])
                         - A[2][2][0] * sum(x ** 2 for x in A[0][1]))
    raise ValueError


def sym_cofactor(A, m):
    """C(A) from the gradient of the Moore determinant with respect to A : B = Re tr(AB)."""
    syms = {}
    P = [[None] * m for _ in range(m)]
    for a in range(m):
        s = sp.Symbol(f"p{a}{a}_0")
        syms[(a, a, 0)] = s
        P[a][a] = (s, 0, 0, 0)
        for b in range(a + 1, m):
            q = tuple(sp.Symbol(f"p{a}{b}_{al}") for al in range(4))
            for al in range(4):
                syms[(a, b, al)] = q[al]
            P[a][b] = q
            P[b][a] = qconj(q)
    detP = sym_moore_det(P, m)
    subs = {}
    for (a, b, al), s in syms.items():
        subs[s] = A[a][b][al]
    C = [[None] * m for _ in range(m)]
    for a in range(m):
        C[a][a] = (sp.expand(sp.diff(detP, syms[(a, a, 0)]).subs(subs)), 0, 0, 0)
        for b in range(a + 1, m):
            comps = tuple(sp.expand(sp.Rational(1, 2) * sp.diff(detP, syms[(a, b, al)]).subs(subs)) for al in range(4))
            C[a][b] = comps
            C[b][a] = qconj(comps)
    return C


def check_moore_formula(rng):
    ok = True
    for _ in range(20):
        A = rand_herm(3, rng)
        Al = [[tuple(A[i, j]) for j in range(3)] for i in range(3)]
        ok &= abs(float(sym_moore_det(Al, 3)) - moore_det(A)) < 1e-9 * max(1, abs(moore_det(A)))
        A2 = rand_herm(2, rng)
        Al2 = [[tuple(A2[i, j]) for j in range(2)] for i in range(2)]
        ok &= abs(float(sym_moore_det(Al2, 2)) - moore_det(A2)) < 1e-9 * max(1, abs(moore_det(A2)))
    return ok


def random_poly(vars_, deg, rng):
    terms = 0
    for _ in range(14):
        mon = 1
        for _ in range(rng.randint(2, deg)):
            mon *= rng.choice(vars_)
        terms += rng.randint(-3, 3) * mon
    return terms


def cofactor_divergence(m, k, seed, perturb=False):
    rng = random.Random(seed)
    xs = [[sp.Symbol(f"x{a}_{al}", real=True) for al in range(k)] for a in range(m)]
    flat = [x for row in xs for x in row]
    u = random_poly(flat, 4, rng)
    H = sym_hessian(u, xs, m, k)
    if perturb:   # a Hermitian field that is not a Hessian
        H[0][0] = (H[0][0][0] + xs[1][0] ** 2 * xs[0][1], 0, 0, 0)
    C = sym_cofactor(H, m)
    worst = 0
    for b in range(m):
        acc = (0, 0, 0, 0)
        for a in range(m):
            for al in range(k):
                dC = tuple(sp.diff(c, xs[a][al]) for c in C[a][b])
                acc = qadd(acc, qmul(qconj(E[al]), dC))
        acc = tuple(sp.expand(x) for x in acc)
        worst = max(worst, sum(1 for x in acc if x != 0))
    return {"m": m, "slots": k, "seed": seed, "perturbed": perturb, "nonzero_components": worst}


def lemma31_checks(rng):
    out = {}
    # H(xi xi^t) = g g^*, trace = Laplacian, product rule
    for k in (3, 4):
        m = 2
        xs = [[sp.Symbol(f"y{a}_{al}", real=True) for al in range(k)] for a in range(m)]
        flat = [x for row in xs for x in row]
        xi = [rng.randint(-3, 3) for _ in flat]
        lin = sum(c * x for c, x in zip(xi, flat))
        H = sym_hessian(lin ** 2 / 2, xs, m, k)
        g = [tuple(sum(xi[a * k + al] * E[al][c] for al in range(k)) for c in range(4)) for a in range(m)]
        ok = all(tuple(sp.simplify(x) for x in H[a][b]) == tuple(sp.simplify(x) for x in qmul(g[a], qconj(g[b])))
                 for a in range(m) for b in range(m))
        u = random_poly(flat, 3, rng)
        v = random_poly(flat, 3, rng)
        Hu = sym_hessian(u, xs, m, k)
        tr_ok = sp.expand(sum(Hu[a][a][0] for a in range(m)) - sum(sp.diff(u, x, 2) for x in flat)) == 0
        Huv = sym_hessian(u * v, xs, m, k)
        Hv = sym_hessian(v, xs, m, k)
        Gu = [tuple(sum(sp.diff(u, xs[a][al]) * E[al][c] for al in range(k)) for c in range(4)) for a in range(m)]
        Gv = [tuple(sum(sp.diff(v, xs[a][al]) * E[al][c] for al in range(k)) for c in range(4)) for a in range(m)]
        prod_ok = True
        for a in range(m):
            for b in range(m):
                rhs = qadd(qadd(qscale(u, Hv[a][b]), qscale(v, Hu[a][b])),
                           qadd(qmul(Gu[a], qconj(Gv[b])), qmul(Gv[a], qconj(Gu[b]))))
                prod_ok &= all(sp.expand(x - y) == 0 for x, y in zip(Huv[a][b], rhs))
        out[f"slots={k}"] = {"H(xi xi^t)=gg*": ok, "trace=Laplacian": tr_ok, "product_rule": prod_ok}
    # (3.4) H(|F|^2) = 4 J^* J for F holomorphic in P in C_k^m, coordinates (x_a, p_a, q_a) in slots (1, i, j)
    m = 2
    xs = [[sp.Symbol(f"z{a}_{al}", real=True) for al in range(3)] for a in range(m)]
    Pc = [xs[a][1] + sp.I * xs[a][2] for a in range(m)]   # P_a = p_a + k q_a, k <-> complex unit
    Pa = sp.symbols("Pa0:2")
    Fpa = [Pa[0] ** 2 + 3 * Pa[0] * Pa[1] - Pa[1] + 1, Pa[1] ** 3 - 2 * Pa[0]]
    back = {Pa[0]: Pc[0], Pa[1]: Pc[1]}
    Fs = [F.subs(back) for F in Fpa]
    u = sp.expand(sum(sp.expand(F * sp.conjugate(F)) for F in Fs))
    u = sp.simplify(sp.re(sp.expand_complex(u)))
    H = sym_hessian(u, xs, m, 3)
    Jc = [[sp.diff(F, Pa[a]).subs(back) for a in range(m)] for F in Fpa]
    # C_k element c = re + k im  ->  quaternion (re, 0, 0, im)
    def ck(cx):
        cx = sp.expand_complex(cx)
        return (sp.re(cx), 0, 0, sp.im(cx))
    ok34 = True
    for a in range(m):
        for b in range(m):
            acc = (0, 0, 0, 0)
            for nu in range(len(Fs)):
                acc = qadd(acc, qmul(qconj(ck(Jc[nu][a])), ck(Jc[nu][b])))
            acc = qscale(4, acc)
            ok34 &= all(sp.simplify(sp.expand(x - y)) == 0 for x, y in zip(H[a][b], acc))
    out["(3.4)_H(|F|^2)=4J*J"] = ok34
    return out


def main():
    rng = np.random.default_rng(11)
    prng = random.Random(5)
    out = {"moore_formula_3x3_matches_eigen": check_moore_formula(rng)}
    out["lemma_3.2"] = moore_checks(rng)
    out["lemma_4.3_segment"] = segment_identity(rng)
    out["lemma_3.1"] = lemma31_checks(prng)
    div = [cofactor_divergence(2, 3, 1), cofactor_divergence(3, 3, 2), cofactor_divergence(3, 3, 3),
           cofactor_divergence(2, 4, 4), cofactor_divergence(3, 2, 5)]
    neg = [cofactor_divergence(3, 3, 6, perturb=True)]
    out["lemma_3.3_cofactor_divergence"] = div
    out["lemma_3.3_negative_control"] = neg
    ok = out["moore_formula_3x3_matches_eigen"] and all(out["lemma_3.2"].values()) and out["lemma_4.3_segment"]
    ok &= all(all(v.values()) if isinstance(v, dict) else v for v in out["lemma_3.1"].values())
    ok &= all(d["nonzero_components"] == 0 for d in div) and all(d["nonzero_components"] > 0 for d in neg)
    out["all_ok"] = bool(ok)
    print(json.dumps(out, indent=1, default=str))
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/quaternion_calculus.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
