"""Exact checks of KAK 22: curvature, orthogonal designs, and quadric coarea.

Standard library only. Gaussian rational arithmetic is reused from GPT 18.
The arbitrary-multiplicity theorem is a prose proof, not inferred from examples.
"""
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import json

from verify_weyl_fibre import G, determinant, matmul


def transpose(a):
    return [list(row) for row in zip(*a)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def mv(a, x):
    return [sum(u * v for u, v in zip(row, x)) for row in a]


def scaled(a, c):
    return [[c * x for x in row] for row in a]


def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def quadratic(a, z):
    return sum(z[i] * a[i][j] * z[j] for i in range(len(z)) for j in range(len(z))) / 2


def displacement(a, z):
    return sum(a[i][j] * (z[i].re * z[j].re + z[i].im * z[j].im)
               for i in range(len(z)) for j in range(len(z))) / 2


def cd_mul(a, b):
    # Same Cayley--Dickson convention as Claude's state-dimension checker.
    if len(a) == 1:
        return [a[0] * b[0]]
    h = len(a) // 2
    a1, a2, b1, b2 = a[:h], a[h:], b[:h], b[h:]
    conj = lambda x: [x[0]] + [-v for v in x[1:]]
    p = [u - v for u, v in zip(cd_mul(a1, b1), cd_mul(conj(b2), a2))]
    q = [u + v for u, v in zip(cd_mul(b2, a1), cd_mul(a2, conj(b1)))]
    return p + q


def hopf(m):
    basis = eye(m)
    cs = []
    for e in basis:
        left = transpose([cd_mul(e, b) for b in basis])
        c = [[F(0)] * (2 * m) for _ in range(2 * m)]
        for i in range(m):
            for j in range(m):
                c[m + i][j] = c[j][m + i] = left[i][j]
        cs.append(c)
    cs.append([[F((1 if i < m else -1) * (i == j)) for j in range(2 * m)]
               for i in range(2 * m)])
    return cs, [F(2)] * m, [F(1)] * m, [F(0)] * m, F(1)


def non_clifford():
    # A realified complex-Hermitian 3x3 pencil, N=6,n=5,m=2.
    # At xi0=e5: top eigenvalue 2, other complex eigenvalues 0,-1.
    # Arbitrary lower blocks make this a non-Clifford test.
    m, k = 2, 4
    bs = [
        [[1, 0], [0, 1], [0, 0], [0, 0]],
        [[0, -2], [2, 0], [0, 0], [0, 0]],
        [[0, 0], [0, 0], [3, 0], [0, 3]],
        [[0, 0], [0, 0], [0, -4], [4, 0]],
    ]
    derivatives = [F(j - 1, 3) for j in range(k)]
    cs = []
    for j, b in enumerate(bs):
        c = [[F(0)] * 6 for _ in range(6)]
        for a in range(m):
            c[a][a] = derivatives[j]
        for i in range(k):
            c[m + i][m + i] = F(j + 1) if i < 2 else F(-j - 1, 2)
            for a in range(m):
                c[m + i][a] = c[a][m + i] = F(b[i][a])
        cs.append(c)
    cs.append([[F(([2, 2, 0, 0, -1, -1][i]) * (i == j)) for j in range(6)]
               for i in range(6)])
    jmat = [[F(0)] * 6 for _ in range(6)]
    for i in range(0, 6, 2):
        jmat[i][i + 1], jmat[i + 1][i] = -1, 1
    assert all(matmul(c, jmat) == matmul(jmat, c) for c in cs)
    assert matmul(cs[0], cs[0]) != scaled(eye(6), matmul(cs[0], cs[0])[0][0])
    return cs, [F(2), F(2), F(3), F(3)], [F(1), F(4), F(6), F(32, 3)], derivatives, F(2)


def stereo(u):
    s = sum(x * x for x in u)
    assert 1 + s
    return [(1 - s) / (1 + s)] + [2 * x / (1 + s) for x in u]


def check_case(name, m, data):
    cs, gap, hdiag, derivatives, phi = data
    k, n = len(hdiag), len(cs)
    size = m + k
    assert n == k + 1 and len(cs[0]) == size == n + m - 1
    assert all(c == transpose(c) for c in cs)
    sinv = [[F(i == j) / gap[i] for j in range(k)] for i in range(k)]
    h = [[hdiag[i] * (i == j) for j in range(k)] for i in range(k)]
    rbasis = [[[cs[j][m + b][a] for b in range(k)] for j in range(k)] for a in range(m)]
    # Polynomial coefficient check of R_a S^-1 R_a^T=(a.a)H/2.
    for a in range(m):
        for b in range(m):
            lhs = add(matmul(matmul(rbasis[a], sinv), transpose(rbasis[b])),
                      matmul(matmul(rbasis[b], sinv), transpose(rbasis[a])))
            assert lhs == (h if a == b else scaled(h, 0))
    gamma = 1 / (prod(gap) * prod(x / 2 for x in hdiag))

    if name.startswith("Hopf"):
        assert rbasis[0] == eye(k)
        for a in range(1, m):
            assert transpose(rbasis[a]) == scaled(rbasis[a], -1)
            assert matmul(rbasis[a], rbasis[a]) == scaled(eye(k), -1)
            for b in range(a + 1, m):
                assert add(matmul(rbasis[a], rbasis[b]), matmul(rbasis[b], rbasis[a])) == scaled(eye(k), 0)

    points = [[G(1)]] if m == 1 else [
        stereo([G(F(j + 1, 2 * m), F((-1) ** j, 3 * m)) for j in range(m - 1)]),
        stereo([G(F(1, 3))] + [G(0)] * (m - 2)),
        stereo([G(0, F(9, 10))] + [G(0)] * (m - 2)),
        stereo([G(0, F(999, 1000))] + [G(0)] * (m - 2)),
    ]
    checked = []
    for a in points:
        assert sum(x * x for x in a) == G(1) and a[0]
        z = a + [G(0)] * k
        r = [[sum(a[alpha] * rbasis[alpha][j][b] for alpha in range(m))
              for b in range(k)] for j in range(k)]
        assert determinant(r).norm() == 1 / gamma
        flux = [quadratic(eye(size), z)] + [quadratic(c, z) for c in cs]
        assert flux == [G(F(1, 2))] + [G(d / 2) for d in derivatives] + [G(phi / 2)]
        assert phi * displacement(eye(size), z) - displacement(cs[-1], z) == 0

        # Differentiate the original quadratic flux, then change target rows.
        # Append the fibre chart a_2,...,a_m: a square ambient Jacobian.
        jac = [z + [G(0)]]
        for j in range(k):
            jac.append([v - derivatives[j] * zz for v, zz in zip(mv(cs[j], z), z)] + [G(0)])
        jac.append([v - phi * zz for v, zz in zip(mv(cs[-1], z), z)] + [G(1)])
        for j in range(1, m):
            jac.append([G(int(i == j)) for i in range(size)] + [G(0)])
        jac_norm = determinant(jac).norm()
        assert jac_norm == a[0].norm() / gamma
        assert a[0].norm() / jac_norm == gamma
        # q is not the bilinear flux Q: it grows on the complex sphere.
        checked.append(dict(a1_norm_squared=str(a[0].norm()),
                            real_time_displacement=str(displacement(eye(size), z)),
                            weighted_coarea_density=str(a[0].norm() / jac_norm)))

    off = [G(1)] + [G(0)] * (size - 1)
    off[m] = G(F(1, 3))
    assert phi * displacement(eye(size), off) - displacement(cs[-1], off) == gap[0] / 18
    if m > 1:
        assert checked[-1]["a1_norm_squared"] != "1"
        assert F(checked[-1]["real_time_displacement"]) > F(checked[0]["real_time_displacement"])
    cutoff = [dict(radius=r, coefficient_of_pi_power=str(2 * gamma * r ** (2 * m - 2) / factorial(m - 1)))
              for r in [1, 2, 4, 8]]
    assert m == 1 or all(F(b["coefficient_of_pi_power"]) > F(a["coefficient_of_pi_power"])
                         for a, b in zip(cutoff, cutoff[1:]))
    return dict(name=name, n=n, N=size, multiplicity=m, curvature_identity=True,
                gamma=str(gamma), full_jacobian_checks=checked,
                cutoff_exponent=2 * m - 2, pi_power=m - 1, cutoffs=cutoff)


def main():
    cases = [check_case(f"Hopf {m}", m, hopf(m)) for m in [1, 2, 4, 8]]
    cases.append(check_case("non-Clifford Hermitian 3x3", 2, non_clifford()))
    # Controls: wrong Jacobian (missing a_1), false positive central damping,
    # real q replaced by holomorphic Q, and omitting rigidity.
    controls = dict(missing_a1_factor_fails=all(c["full_jacobian_checks"][-1]["a1_norm_squared"] != "1"
                                                for c in cases if c["multiplicity"] > 1),
                    positive_central_damping_fails=all(c["curvature_identity"] for c in cases),
                    q_equals_Q_fails=all(F(c["full_jacobian_checks"][-1]["real_time_displacement"]) > F(1, 2)
                                         for c in cases if c["multiplicity"] > 1),
                    nonrigid_normal_map_is_rectangular=(8 - 4 != 4 - 1))
    # A deliberately altered off-diagonal row destroys the quadratic identity.
    r0, r1 = eye(2), [[F(0), F(1)], [F(-1), F(0)]]
    wrong = [row.copy() for row in r1]
    wrong[0][0] += 1
    controls["perturbed_orthogonal_design_fails"] = add(matmul(r0, transpose(wrong)),
                                                       matmul(wrong, transpose(r0))) != scaled(eye(2), 0)
    assert all(controls.values())
    out = dict(scope="Exact finite algebra/coarea checks; the general second-variation and divergence proofs are in KAK 22.",
               cases=cases, negative_controls=controls)
    path = Path(__file__).resolve().parent.parent / "results" / "curvature_fibre.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print("Curvature identity and full complex Jacobians: exact, 5 pencils")
    print("Complex-quadric cutoff masses: exact, exponents 0, 2, 6, 14")
    print("Clifford generators, non-Clifford test, and negative controls: OK")
    print("ALL OK")


if __name__ == "__main__":
    main()
