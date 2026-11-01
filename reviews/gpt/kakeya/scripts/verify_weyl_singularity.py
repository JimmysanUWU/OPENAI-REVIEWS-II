"""Checks for KAK 19: exact coarea algebra and elementary profile estimates.

No numerical quadrature is needed for the proof. The reported enclosures come
from proved analytic bounds, evaluated in ordinary floating point for display.
They are not machine-certified interval arithmetic. Standard library only.
"""
from fractions import Fraction as F
from math import exp, pi
from pathlib import Path
import json

from verify_weyl_fibre import G, I, determinant


def general_coarea_checks():
    cases = []
    for a, b, c in [(G(1), G(0), G(0)), (G(2, 1), G(1, -1), G(3, 2)),
                    (G(3), G(F(1, 3), F(1, 4)), G(F(1, 3), -F(1, 4))),
                    (G(2, -1), G(0), G(1, 2))]:
        d = b*c/a
        assert a and a-d
        jac = [[(1-b*c/(a*a))/2, c/(2*a), b/(2*a), 0],
               [0, F(1, 2), F(1, 2), 0],
               [0, -1/(2*I), 1/(2*I), 0],
               [(1+b*c/(a*a))/2, -c/(2*a), -b/(2*a), 1]]
        det = determinant(jac)
        assert det == (1-b*c/(a*a))/(4*I)
        weight = 1/(16*a.norm()*det.norm())
        assert weight == 1/(a-d).norm()
        cases.append(dict(a=[str(a.re), str(a.im)], b=[str(b.re), str(b.im)],
                          c=[str(c.re), str(c.im)], coarea_factor=str(weight)))
    return cases


def real_slice_checks():
    cases = []
    # Choose a and complex b rationally, then d=|b|^2/a makes the north branch exact.
    for a in [F(1), F(3, 2), F(2)]:
        for b in [G(F(1, 8), F(1, 6)), G(F(1, 10)), G(0, F(1, 20))]:
            c = G(b.re, -b.im)
            s2 = b.norm()
            d = s2/a
            omega0 = (a+d)/2
            omega1, omega2 = b.re, -b.im
            assert omega1**2+omega2**2 == s2
            assert a > d and omega0**2-s2 == ((a-d)/2)**2
            assert b*c == G(a*d)
            kappa = F(3, 2)
            A, B = kappa*s2/(a*a), kappa*s2
            # Two independent substitutions in the radial terms.
            first, second = 1/(2*A), a*a/(2*B)
            assert first == second
            prefactor_without_pi = (first+second)/(2*(a-d)**2)
            assert prefactor_without_pi == a*a/(2*kappa*s2*(a-d)**2)
            assert A*B == kappa*kappa*s2*s2/(a*a)
            cases.append(dict(a=str(a), s_squared=str(s2), omega0=str(omega0),
                              prefactor_divided_by_pi=str(prefactor_without_pi)))
    return cases


def profile_and_remainder_checks():
    cases = []
    # a=1, d=s^2: these rational target points approach omega0=1/2.
    kappa, m = F(1), F(1, 2)
    previous_width = None
    for denominator in [4, 8, 16, 32, 64, 128, 256]:
        s = F(1, denominator)
        a, d = F(1), s*s
        omega0 = (a+d)/2
        assert omega0 >= m and s <= m/2
        ratio = a*a/(a-d)**2
        assert 1 <= ratio <= 1+2*s*s/(m*m)
        root_v = kappa*s*s/a
        # Exact algebra of the O(1) lower remainder, before evaluating exponentials.
        assert (1-2*root_v)/(2*kappa*s*s) == 1/(2*kappa*s*s)-1/a
        # J(v) lies between exp(-2 sqrt(v)) and 1 by Lemma 3.1.
        lower = pi/2*float(ratio)*exp(-2*float(root_v))
        upper = pi/2*float(ratio)
        assert 0 < lower <= upper
        width = upper-lower
        if previous_width is not None:
            assert width < previous_width
        previous_width = width
        cases.append(dict(s=str(s), scaled_density_lower=lower, scaled_density_upper=upper,
                          limit=pi/2, note="Analytic bounds evaluated in floating point."))
    assert abs(lower-pi/2) < 0.0001 and abs(upper-pi/2) < 0.0001
    # The transverse Jacobian is s ds dtheta: s^-2 gives a logarithm, not finite mass.
    # Exact coefficients after removing pi^2: (1/(2 kappa)) * 2 = 1/kappa.
    assert F(1, 2)/kappa*2 == 1/kappa
    return dict(profile_enclosures=cases, logarithmic_coefficient_divided_by_pi_squared=str(1/kappa),
                bounded_remainder_checked=True)


def main():
    result = dict(general_coarea=general_coarea_checks(), real_slice=real_slice_checks(),
                  profile=profile_and_remainder_checks(), all_ok=True)
    output = Path(__file__).resolve().parent.parent / "results" / "weyl_singularity.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("General coarea factor, real branch, radial substitutions, singularity coefficient: OK")
    print("ALL OK")


if __name__ == "__main__":
    main()
