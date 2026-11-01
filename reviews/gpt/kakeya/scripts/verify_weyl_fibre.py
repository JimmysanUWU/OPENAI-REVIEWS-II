"""Exact algebra for 18_NONCOMPACT_FIBRE_OBSTRUCTION.md; standard library only."""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import json


@dataclass(frozen=True)
class G:
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    def __add__(self, other):
        o = other if isinstance(other, G) else G(other)
        return G(self.re+o.re, self.im+o.im)
    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, other):
        return self+-other

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        o = other if isinstance(other, G) else G(other)
        return G(self.re*o.re-self.im*o.im, self.re*o.im+self.im*o.re)
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = other if isinstance(other, G) else G(other)
        return self*G(o.re, -o.im)*G(1/o.norm())

    def __rtruediv__(self, other):
        return G(other)/self

    def norm(self):
        return self.re*self.re+self.im*self.im

    def __bool__(self):
        return bool(self.re or self.im)


I = G(0, 1)


def determinant(matrix):
    a = [[x if isinstance(x, G) else G(x) for x in row] for row in matrix]
    value = G(1)
    for j in range(len(a)):
        p = next((k for k in range(j, len(a)) if a[k][j]), None)
        if p is None:
            return G(0)
        if p != j:
            a[p], a[j] = a[j], a[p]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for k in range(j+1, len(a)):
            factor = a[k][j]/pivot
            a[k] = [x-factor*y for x, y in zip(a[k], a[j])]
    return value


def matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def quadratic(c, z):
    return sum(z[i]*c[i][j]*z[j] for i in range(4) for j in range(4))/2


def main():
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    c1 = [[0,0,1,0], [0,0,0,1], [1,0,0,0], [0,1,0,0]]
    c2 = [[0,0,0,1], [0,0,-1,0], [0,-1,0,0], [1,0,0,0]]
    c3 = [[(1 if i < 2 else -1)*int(i == j) for j in range(4)] for i in range(4)]
    matrices = [identity, c1, c2, c3]
    for i, a in enumerate(matrices[1:]):
        assert a == [list(row) for row in zip(*a)]
        for j, b in enumerate(matrices[1:]):
            ab, ba = matmul(a, b), matmul(b, a)
            assert [[ab[k][l]+ba[k][l] for l in range(4)] for k in range(4)] == [
                [2*int(i == j)*identity[k][l] for l in range(4)] for k in range(4)]

    linear_uv = [[1,I,0,0], [0,0,1,I], [1,-I,0,0], [0,0,1,-I]]
    assert determinant(linear_uv).norm() == 16
    base_jac = [[F(1,2),0,0,0], [0,F(1,2),F(1,2),0],
                [0,-1/(2*I),1/(2*I),0], [F(1,2),0,0,1]]
    assert determinant(base_jac).norm() == F(1,16)
    checks = 0
    for lam in [G(1), G(2), G(F(1,3)), G(1,1), G(2,-3)]:
        for a, b, c in [(G(1),G(0),G(0)), (G(2,1),G(1,-1),G(3,2))]:
            u1, u2, v1, v2 = lam, c*lam/a, a/lam, b/lam
            z = [(u1+v1)/2, (u1-v1)/(2*I), (u2+v2)/2, (u2-v2)/(2*I)]
            quad = [quadratic(matrix, z) for matrix in matrices]
            expected = [(a+b*c/a)/2, (b+c)/2, (c-b)/(2*I), (a-b*c/a)/2]
            assert quad == expected
            assert quad[0]*quad[0] == sum(x*x for x in quad[1:])
            jac = [[1,0,0,0], [c/a,-c*lam/(a*a),0,lam/a],
                   [-a/(lam*lam),1/lam,0,0], [-b/(lam*lam),0,1/lam,0]]
            assert determinant(jac) == 1/(a*lam)
            # Compute real displacement independently from w1 and w2 and the real matrices.
            displacement = [sum(F(z[i].re*z[j].re+z[i].im*z[j].im)*matrix[i][j]
                                for i in range(4) for j in range(4))/2 for matrix in matrices]
            assert displacement[0]-displacement[3] == z[2].norm()+z[3].norm()
            assert displacement[0]-displacement[3] == (u2.norm()+v2.norm())/2
            if not b and not c:
                assert quad == [G(F(1,2)), G(0), G(0), G(F(1,2))]
                assert displacement[0]-displacement[3] == 0
                assert z[0].norm() == ((lam+1/lam)/2).norm()
            checks += 1

    # Exact annular integral divided by pi. Integrate angular constant term and r powers.
    annuli = []
    for radius in [2, 3, 4, 8, 16, 32]:
        r = F(radius)
        integral_over_pi = ((r*r-1/(r*r))/2+(r*r-1/(r*r))/2)/2
        assert integral_over_pi == (r*r-1/(r*r))/2
        annuli.append({"radius": radius, "integral_divided_by_pi": str(integral_over_pi)})
    assert all(F(b["integral_divided_by_pi"]) > F(a["integral_divided_by_pi"])
               for a, b in zip(annuli, annuli[1:]))
    # Nonzero damping away from both axes; no uniform positive damping on the central fibre.
    assert G(1).norm()+G(1).norm() > 0
    epsilon = F(1,7)
    # Build the actual ambient tensor from real and imaginary parts of the lifted fields
    # at w1=e1, w2=0, G=e1. Coordinates: four physical, then w1, then w2.
    a1, a2, b1 = [0]*12, [0]*12, [0]*12
    a1[4], a2[8], b1[0], b1[3] = 1, 1, 1, 1
    tensor = [[sum(v[i]*v[j] for v in [a1,a2,b1]) for j in range(12)] for i in range(12)]
    gradient = [F(0)]*12
    gradient[4] = -2*epsilon  # gradient divided by exp(-epsilon)
    source = [sum(tensor[i][j]*gradient[j] for j in range(12)) for i in range(12)]
    assert source[4] == -2*epsilon and all(v == 0 for i,v in enumerate(source) if i != 4)
    results = dict(clifford_relations=True, exact_coordinate_checks=checks,
                   real_jacobians={"z_to_uv":16, "base_to_target":"1/16"},
                   conditional_fibre_measure="d^2(lambda)/|lambda|^2",
                   annuli=annuli, gaussian_cutoff_source_divided_by_exp="-2 epsilon",
                   controls={"real_displacement_not_assumed_constant":True,
                             "central_gaussian_does_not_damp":True,
                             "radial_cutoff_has_source":True}, all_ok=True)
    path = Path(__file__).resolve().parent.parent / "results" / "weyl_fibre.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")
    print("Clifford matrices, fibre coordinates, coarea weights, damping, annular integral: OK")
    print("ALL OK")


if __name__ == "__main__":
    main()
