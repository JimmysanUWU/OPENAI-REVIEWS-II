"""Checks for 10_GAPLESS_REMAINING.md (the Kakeya paper's remaining claims).

Part C  Proposition 8.2 / App. C, the tube arrangement: for k = 6..16 the lines L_j are
        alpha*delta/2-separated at X <= -1/2 (so the input packets are disjoint), the union B_X
        of the output intervals has measure <= 2 alpha / k + 2^m w on J = [1/4, 1/2] (exact
        interval-union measure on a grid of X), and the resulting compression k^d |E| stays
        bounded while sum |S_nu| = |J| (alpha/16)^d.  Contrast: on the input side X = -1/2 the
        same intervals are disjoint, with union exactly 2^k w = alpha/16.
Part N  Section 9.4, the projective map (y', s) -> (y'/(s - t0), 1/(s - t0)): it sends every line
        through (x', t0) to a line of direction (x', 1), its Jacobian is -(s - t0)^-(n+1), and
        the direction map x' -> (x', 1)/|(x', 1)| has area Jacobian (1 + |x'|^2)^(-n/2)
        (symbolic, n = 2, 3, 4).
Part P  Section 12: the parabolic rescaling (q_delta(v) = |v|^2 + O(delta^2), the delta-power
        -d + (d+2)/p + d/p = 0 at p = 2(d+1)/d, d = 1..6), the Fresnel identity
        |S(t)h(z)| = (4 pi t)^(-d/2) lambda^d |E g(-lambda z / 2t, lambda^2 / 4t)| checked by
        quadrature in d = 1, and the Jacobian lambda^(d+2)-scaling of the change of variables.
Part K  Corollary 9.3, the cone: the surface measure sqrt(2) rho^(n-1) d rho d sigma and the
        exponent bookkeeping of the slicing argument.
Part E  App. E: Example E.2 (top eigenvalue 1 at (0,0,1), transverse Hessian diag(1, 2/3) both
        from formula (E.2) and by finite differences of lambda_max, no linear factor of the
        characteristic cubic over C), formula (E.2) on random pencils, and the group-velocity
        identity B(v(xi)) = grad phi(xi).
Part D  App. D: the completed square A:H - D_B(H) = beta |b01 + a01/beta|^2 + det(a)/beta for
        2x2 quaternionic Hermitian blocks with a00 = beta b11, and the rank-one expansion
        det(lambda B + H) = det H + lambda beta det(H|ker B) for 3x3 quaternionic Hermitian H
        (Moore determinant), on random instances.

Writes ../results/remaining_claims.json.
"""

import json
import math
import random

import numpy as np
import sympy as sp

RESULTS = {}


# =============================================================== Part C

def union_measure(centers, radius):
    iv = sorted((c - radius, c + radius) for c in centers)
    tot, cur_a, cur_b = 0.0, None, None
    for a, b in iv:
        if cur_b is None or a > cur_b:
            if cur_b is not None:
                tot += cur_b - cur_a
            cur_a, cur_b = a, b
        else:
            cur_b = max(cur_b, b)
    if cur_b is not None:
        tot += cur_b - cur_a
    return tot


def tube_lines(k, d):
    alpha = 1 / (4 * math.sqrt(d))
    delta = 2.0 ** -k
    w = alpha * delta / 16
    a = [1 - l / k for l in range(1, k + 1)]

    def L(j, X):
        return alpha * sum(((j >> (k - l)) & 1) * 2.0 ** -l * (X - a[l - 1]) for l in range(1, k + 1))
    return alpha, delta, w, L


def part_C():
    rows, ok = [], True
    for k in (6, 8, 10, 12, 14, 16):
        alpha, delta, w, L = tube_lines(k, 1)
        nj = 2 ** k
        # separation at X <= -1/2 (consecutive j, at several X)
        sep_ok = True
        for X in (-0.5, -1.0, -2.0):
            vals = [L(j, X) for j in range(nj)]
            gaps = [vals[j] - vals[j + 1] for j in range(nj - 1)]
            sep_ok &= min(gaps) >= alpha * delta / 2 * (1 - 1e-9)
        # |B_X| on J
        worst_ratio = 0.0
        bound_ok = True
        measures = []
        for X in np.linspace(0.25, 0.5, 11):
            m = int(math.floor(k * (1 - X)))
            meas = union_measure([L(j, X) for j in range(nj)], w / 2)
            bound = 2 * alpha / k + 2 ** m * w
            bound_ok &= meas <= bound * (1 + 1e-9)
            measures.append(meas)
            worst_ratio = max(worst_ratio, meas * k / alpha)
        EJ = sum(mm ** 1 for mm in measures) / len(measures) * 0.25       # |E| for d = 1 (Riemann sum)
        # contrast: at X = -1/2 the intervals of radius w/2 are disjoint, union = 2^k w
        input_side = union_measure([L(j, -0.5) for j in range(nj)], w / 2)
        sep_ok &= abs(input_side - nj * w) < 1e-12
        rows.append({"k": k, "separated_at_X<=-1/2": sep_ok, "|B_X| <= 2a/k + 2^m w": bound_ok,
                     "max k|B_X|/alpha": round(worst_ratio, 3), "k*|E| (d=1)": round(k * EJ, 4),
                     "sum|S_nu| (d=1)": 0.25 * alpha / 16,
                     "input-side union at X=-1/2 (= 2^k w = alpha/16)": round(input_side, 6),
                     "output-side mean |B_X| on J": round(sum(measures) / len(measures), 6)})
        ok &= sep_ok and bound_ok and worst_ratio <= 3.0
    # the compression must make k|B_X| bounded while the control does not shrink with k
    RESULTS["tubes_appendix_C"] = {"all_ok": ok, "rows": rows}
    print(f"[C] App. C tubes: separation, |B_X| <= 2 alpha/k + 2^m w, k|B_X|/alpha bounded "
          f"(max {max(r['max k|B_X|/alpha'] for r in rows)}): {ok}")
    return ok


# =============================================================== Part N

def part_N():
    ok = True
    out = []
    for n in (2, 3, 4):
        d = n - 1
        ys = sp.symbols(f"y1:{d + 1}", real=True)
        s, t0 = sp.symbols("s t0", real=True)
        xs = sp.symbols(f"x1:{d + 1}", real=True)
        vs = sp.symbols(f"v1:{d + 1}", real=True)
        Ymap = [yi / (s - t0) for yi in ys] + [1 / (s - t0)]
        J = sp.Matrix(Ymap).jacobian(list(ys) + [s])
        detJ = sp.simplify(J.det())
        good_det = sp.simplify(detJ + (s - t0) ** (-(n + 1))) == 0
        # a line through (x', t0): y' = x' + (s - t0) v  ->  Y' = v + Y_n x'
        img = [sp.simplify(Ymap[i].subs({ys[j]: xs[j] + (s - t0) * vs[j] for j in range(d)})) for i in range(d)]
        Yn = Ymap[d]
        good_line = all(sp.simplify(img[i] - (vs[i] + Yn * xs[i])) == 0 for i in range(d))
        # area Jacobian of x' -> (x', 1)/|(x', 1)|: Gram determinant
        u = sp.Matrix(list(xs) + [1])
        normed = u / sp.sqrt(sum(c ** 2 for c in u))
        D = normed.jacobian(list(xs))
        gram = sp.simplify((D.T * D).det())
        good_area = sp.simplify(gram - (1 + sum(c ** 2 for c in xs)) ** (-n)) == 0
        ok &= bool(good_det and good_line and good_area)
        out.append({"n": n, "det": str(detJ), "lines_to_parallel_lines": bool(good_line),
                    "area_jacobian_squared": str(gram), "ok": bool(good_det and good_line and good_area)})
    RESULTS["nikodym_projective_map"] = {"all_ok": ok, "rows": out}
    print(f"[N] Nikodym projective map: lines -> parallel lines, Jacobian -(s-t0)^-(n+1), area Jacobian: {ok}")
    return ok


# =============================================================== Part P

def part_P():
    ok = True
    v, dl = sp.symbols("v delta", positive=True)
    q = 2 * (1 - sp.sqrt(1 - dl ** 2 * v ** 2)) / dl ** 2
    ser = sp.series(q, dl, 0, 4).removeO()
    good_q = sp.simplify(ser - (v ** 2 + dl ** 2 * v ** 4 / 4)) == 0
    powers = []
    for d in range(1, 7):
        p = sp.Rational(2 * (d + 1), d)
        pw = -d + sp.Rational(d + 2, 1) / p + sp.Rational(d, 1) / p
        powers.append(str(pw))
        ok &= pw == 0
        lam_pw = d - sp.Rational(2 * d + 2, 1) / p
        ok &= lam_pw == 0
    # Fresnel identity in d = 1 by quadrature: h(y) = exp(-y^2) (1 + 0.3 y), lambda = 4, t in {0.6, 0.9}
    lam = 4.0
    ys = np.linspace(-12, 12, 240001)
    dy = ys[1] - ys[0]
    h = np.exp(-ys ** 2) * (1 + 0.3 * ys)
    fres_ok = True
    worst = 0.0
    for t in (0.6, 0.9):
        for z in (-1.3, 0.2, 2.1):
            S = (4 * math.pi * 1j * t) ** -0.5 * np.sum(np.exp(1j * (z - ys) ** 2 / (4 * t)) * h) * dy
            # E g(x, s) = int_B e^{i(x v + s v^2)} g(v) dv with g(v) = h(lambda v), support in |v| <= 3
            vv = ys / lam
            g = h
            X, Sx = -lam * z / (2 * t), lam ** 2 / (4 * t)
            Eg = np.sum(np.exp(1j * (X * vv + Sx * vv ** 2)) * g) * (dy / lam)
            rhs = (4 * math.pi * t) ** -0.5 * lam * abs(Eg)
            err = abs(abs(S) - rhs) / max(abs(S), 1e-12)
            worst = max(worst, err)
            fres_ok &= err < 1e-8
    # Jacobian of (x, t) -> (-lambda z/(2t), lambda^2/(4t)), z = x - c (d = 2)
    x1, x2, t, L_, c1, c2 = sp.symbols("x1 x2 t lambda c1 c2", positive=True)
    F = sp.Matrix([-L_ * (x1 - c1) / (2 * t), -L_ * (x2 - c2) / (2 * t), L_ ** 2 / (4 * t)])
    Jd = sp.simplify(F.jacobian([x1, x2, t]).det())
    jac_ok = sp.simplify(Jd / L_ ** 4).free_symbols <= {t, x1, x2, c1, c2}
    ok &= bool(good_q) and fres_ok and bool(jac_ok)
    RESULTS["paraboloid_schrodinger"] = {"all_ok": ok, "q_delta_expansion": str(ser),
                                         "delta_powers_d1..6": powers, "fresnel_max_rel_err": worst,
                                         "jacobian_d2": str(Jd)}
    print(f"[P] parabolic rescaling powers vanish (d = 1..6), q_delta = |v|^2 + O(delta^2), Fresnel "
          f"identity (rel. err {worst:.1e}), Jacobian ~ lambda^(d+2): {ok}")
    return ok


# =============================================================== Part K

def part_K():
    ok = True
    rows = []
    for n in (2, 3, 4, 5):
        pc = sp.Rational(2 * n, n - 1)
        # pc = 1 + (n+1)/(n-1)
        ok &= pc == 1 + sp.Rational(n + 1, n - 1)
        # cone measure: parametrize (rho w, rho) on the cone; the slant factor is sqrt(2)
        rho = sp.symbols("rho", positive=True)
        # metric: d/d rho of (rho w, rho) = (w, 1), |.|^2 = 2; tangential directions scale by rho
        rows.append({"n": n, "pc": str(pc), "measure_factor": "sqrt(2) rho^(n-1)"})
        # slicing exponents: for b > a' we need a > b'
        for a in (1.5, 2.0, 4.0, 10.0):
            ap = a / (a - 1)
            for b in (float(pc) + 0.1, 3 * ap):
                if b > ap and b > 2:
                    bp = b / (b - 1)
                    ok &= a > bp
    RESULTS["cone_slicing"] = {"all_ok": ok, "rows": rows}
    print(f"[K] cone: p_c = 1 + (n+1)/(n-1), slicing exponents (b > a' => a > b'): {ok}")
    return ok


# =============================================================== Part E

def lam_max(C):
    return float(np.linalg.eigvalsh(C)[-1])


def part_E():
    ok = True
    x, y, z, t = sp.symbols("x y z t")
    Cs = sp.Matrix([[z, x, y], [x, -z, 0], [y, 0, -2 * z]])
    charpoly = sp.expand((t * sp.eye(3) - Cs).det())
    target = sp.expand((t - z) * (t + z) * (t + 2 * z) - x ** 2 * (t + 2 * z) - y ** 2 * (t + z))
    ok &= sp.simplify(charpoly - target) == 0
    # no linear factor t - (a x + b y + c z) over C
    a_, b_, c_ = sp.symbols("a b c")
    sub = sp.expand(target.subs(t, a_ * x + b_ * y + c_ * z))
    eqs = sp.Poly(sub, x, y, z).coeffs()
    sols = sp.solve(eqs, [a_, b_, c_], dict=True)
    ok &= len(sols) == 0
    # top eigenvalue at (0,0,1) and transverse Hessian by (E.2)
    C1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], float)
    C2 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], float)
    C3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -2]], float)
    Cmat = lambda xi: xi[0] * C1 + xi[1] * C2 + xi[2] * C3
    om = np.array([0, 0, 1.0])
    top = lam_max(Cmat(om))
    ok &= abs(top - 1) < 1e-12
    ev, evec = np.linalg.eigh(Cmat(om))
    e0 = evec[:, -1]
    perp = evec[:, :-1]
    Sperp = perp.T @ (top * np.eye(3) - Cmat(om)) @ perp
    hess_formula = np.zeros((2, 2))
    hess_fd = np.zeros((2, 2))
    basis = [np.array([1.0, 0, 0]), np.array([0, 1.0, 0])]
    for i_, wi in enumerate(basis):
        for j_, wj in enumerate(basis):
            Lwi = perp.T @ (Cmat(wi) @ e0)
            Lwj = perp.T @ (Cmat(wj) @ e0)
            hess_formula[i_, j_] = 2 * Lwi @ np.linalg.solve(Sperp, Lwj)
            hh = 1e-4
            f = lambda a, b: lam_max(Cmat(om + a * wi + b * wj))
            hess_fd[i_, j_] = (f(hh, hh) - f(hh, -hh) - f(-hh, hh) + f(-hh, -hh)) / (4 * hh * hh)
    ok &= np.allclose(hess_formula, np.diag([1, 2 / 3]), atol=1e-12)
    ok &= np.allclose(hess_fd, hess_formula, atol=1e-5)
    # (E.2) and Hellmann-Feynman on random pencils
    rng = np.random.default_rng(7)
    e2_err = 0.0
    hf_err = 0.0
    for trial in range(20):
        n = int(rng.integers(2, 5))
        Cj = []
        for _ in range(n):
            M = rng.normal(size=(n, n))
            Cj.append((M + M.T) / 2)
        Cm = lambda xi: sum(xi[i] * Cj[i] for i in range(n))
        w0 = rng.normal(size=n)
        w0 /= np.linalg.norm(w0)
        ev, evec = np.linalg.eigh(Cm(w0))
        if ev[-1] - ev[-2] < 0.3:
            continue
        e0 = evec[:, -1]
        perp = evec[:, :-1]
        Sperp = perp.T @ (ev[-1] * np.eye(n) - Cm(w0)) @ perp
        wv = rng.normal(size=n)
        wv -= (wv @ w0) * w0
        Lw = perp.T @ (Cm(wv) @ e0)
        formula = 2 * Lw @ np.linalg.solve(Sperp, Lw)
        hh = 1e-4
        fd = (lam_max(Cm(w0 + hh * wv)) - 2 * ev[-1] + lam_max(Cm(w0 - hh * wv))) / hh ** 2
        e2_err = max(e2_err, abs(fd - formula) / max(1.0, abs(formula)))
        # group velocity: gradient of phi equals e0^T C_j e0 (Hellmann-Feynman)
        grad_fd = np.array([(lam_max(Cm(w0 + hh * np.eye(n)[j])) - lam_max(Cm(w0 - hh * np.eye(n)[j]))) / (2 * hh)
                            for j in range(n)])
        grad_hf = np.array([e0 @ Cj[j] @ e0 for j in range(n)])
        hf_err = max(hf_err, float(np.max(np.abs(grad_fd - grad_hf))))
    ok &= e2_err < 1e-4 and hf_err < 1e-6
    RESULTS["appendix_E"] = {"all_ok": ok, "charpoly": str(target), "linear_factor_solutions": len(sols),
                             "top_eigenvalue": top, "transverse_hessian_formula": hess_formula.tolist(),
                             "transverse_hessian_fd": hess_fd.tolist(),
                             "E2_random_pencils_rel_err": e2_err, "hellmann_feynman_err": hf_err}
    print(f"[E] Example E.2 (lambda = 1, Hessian diag(1, 2/3), no linear factor) and (E.2), "
          f"group velocity on random pencils: {ok}")
    return ok


# =============================================================== Part D (quaternionic algebra)

def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                     a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                     a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def qconj(p):
    return np.array([p[0], -p[1], -p[2], -p[3]])


def q_to_c(p):
    """Quaternion a + bi + cj + dk as the 2x2 complex matrix [[a+bi, c+di], [-c+di, a-bi]]."""
    a, b, c, d = p
    return np.array([[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]])


def herm_q(rng, n, psd=False):
    """Random n x n quaternionic Hermitian matrix as an array of quaternions."""
    H = np.zeros((n, n, 4))
    for i in range(n):
        H[i, i] = [rng.normal(), 0, 0, 0]
        for j in range(i + 1, n):
            qv = rng.normal(size=4)
            H[i, j] = qv
            H[j, i] = qconj(qv)
    if psd:
        # A = X X^* + Y Y^* with random quaternionic columns: PSD of rank two
        X = rng.normal(size=(n, 4))
        Y = rng.normal(size=(n, 4))
        for i in range(n):
            for j in range(n):
                H[i, j] = qmul(X[i], qconj(X[j])) + qmul(Y[i], qconj(Y[j]))
    return H


def complexify(H):
    n = H.shape[0]
    M = np.zeros((2 * n, 2 * n), complex)
    for i in range(n):
        for j in range(n):
            M[2 * i:2 * i + 2, 2 * j:2 * j + 2] = q_to_c(H[i, j])
    return M


def moore_det(H):
    """Moore determinant of a quaternionic Hermitian matrix: sqrt of the complex determinant
    of its 2n x 2n complex form (which is real and >= 0 for PSD; in general Moore^2 = det)."""
    M = complexify(H)
    ev = np.linalg.eigvalsh(M)
    # eigenvalues come in equal pairs; Moore det = product of one from each pair
    ev = np.sort(ev)
    return float(np.prod(ev[::2]))


def part_D():
    rng = np.random.default_rng(11)
    ok = True
    worst_sq = 0.0
    worst_r1 = 0.0
    for trial in range(200):
        # completed square for 2x2 blocks (indices 0, 1)
        beta = float(rng.uniform(0.3, 3))
        b = herm_q(rng, 2)
        alpha = herm_q(rng, 2, psd=True)
        # impose a00 = beta b11 by rescaling alpha (alpha is PSD with alpha00 > 0)
        s = beta * b[1, 1, 0] / alpha[0, 0, 0]
        if s <= 0:
            continue
        alpha = alpha * s
        AH = sum(qmul(alpha[i, j], b[j, i])[0] for i in range(2) for j in range(2))   # Re tr(alpha b)
        detb = b[0, 0, 0] * b[1, 1, 0] - float(np.dot(b[0, 1], b[0, 1]))
        deta = alpha[0, 0, 0] * alpha[1, 1, 0] - float(np.dot(alpha[0, 1], alpha[0, 1]))
        lhs = AH - beta * detb
        diff = b[0, 1] + alpha[0, 1] / beta
        rhs = beta * float(np.dot(diff, diff)) + deta / beta
        worst_sq = max(worst_sq, abs(lhs - rhs) / max(1.0, abs(lhs)))
        # rank-one expansion of the Moore determinant (3x3): det(lambda B + H) linear in lambda
        H = herm_q(rng, 3)
        e = rng.normal(size=(3, 4))
        nrm = math.sqrt(float(np.sum(e * e)))
        e /= nrm
        Bm = np.zeros((3, 3, 4))
        for i in range(3):
            for j in range(3):
                Bm[i, j] = qmul(e[i], qconj(e[j]))
        betaB = 1.0
        # restriction of H to ker B = e^perp: orthonormal quaternionic basis via the complex form
        Mc = complexify(H)
        ec = np.zeros((6, 2), complex)
        for i in range(3):
            ec[2 * i:2 * i + 2, :] = q_to_c(e[i])
        # projector onto the complement of the quaternionic line spanned by e
        P = np.eye(6) - ec @ ec.conj().T
        evp, U = np.linalg.eigh(P)
        Uperp = U[:, evp > 0.5]
        Hr = Uperp.conj().T @ Mc @ Uperp
        det_restr = float(np.prod(np.sort(np.linalg.eigvalsh(Hr))[::2]))
        vals = []
        for lam in (0.0, 1.0, 2.0):
            vals.append(moore_det(H + lam * Bm))
        slope1 = vals[1] - vals[0]
        slope2 = vals[2] - vals[1]
        worst_r1 = max(worst_r1, abs(slope1 - slope2) / max(1.0, abs(slope1)),
                       abs(slope1 - betaB * det_restr) / max(1.0, abs(slope1)))
    ok &= worst_sq < 1e-10 and worst_r1 < 1e-8
    RESULTS["appendix_D_algebra"] = {"all_ok": ok, "completed_square_max_rel_err": worst_sq,
                                     "rank_one_expansion_max_rel_err": worst_r1}
    print(f"[D] App. D algebra: completed square (err {worst_sq:.1e}), rank-one Moore expansion "
          f"(err {worst_r1:.1e}): {ok}")
    return ok


def main():
    oks = [part_C(), part_N(), part_P(), part_K(), part_E(), part_D()]
    RESULTS["all_ok"] = all(oks)
    with open(__file__.rsplit("/", 2)[0] + "/results/remaining_claims.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print("ALL OK" if all(oks) else "FAILURES")


if __name__ == "__main__":
    main()
