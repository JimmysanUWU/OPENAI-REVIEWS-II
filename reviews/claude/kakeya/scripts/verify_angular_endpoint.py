"""Checks for 14_GAPLESS_ANGULAR_ENDPOINT.md (the paper's App. B, Theorem B.1, line by line).

Part X  exponent algebra (exact, sympy, symbolic n): p = 2n/d, alpha = d/2, alpha p = n, d - alpha p = -1,
        alpha - d/p = 1/p, alpha - 1/p = d/p, n - alpha = alpha + 1, the small-source exponent, the threshold
        B > A + 1/2 - 1/p = 1/2 + 1/p at A = 2/p, the shell-transplant error R^{2+n/p}/S = 1/(100 R).
Part F  Funk-Hecke (Lemma 14.2): int_{S^d} e^{i z th.om} Y_nu(om) dom = (2pi)^{n/2} i^nu z^{1-n/2} J_{nu+n/2-1}(z) Y_nu(th)
        on S^1 and S^2 by quadrature, and the two-term form with constants (2pi)^alpha e^{-+ i pi alpha/2}.
        Negative control: the same formula without i^nu fails.
Part H  Hankel's integral (Lemma 14.3): H^(1)_mu(x) = sqrt(2/(pi x)) e^{i(x - mu pi/2 - pi/4)} / Gamma(mu+1/2)
        * int_0^oo e^{-u} u^{mu-1/2} (1 + iu/(2x))^{mu-1/2} du, against mpmath's hankel1.
        Negative control: exponent mu + 1/2 in the bracket fails.
Part E  the explicit remainder bound (Lemma 14.3): for mu >= 0, x >= (mu+1)^2, k = 0..3,
        |(x d/dx)^k eps_mu(x)| <= e^{3/2} (2mu+2k+1)^{max(k-1,0)} (mu+1)^2 / (2x).
        Necessity of the degree threshold: x |eps_mu(x)| at x = (mu+1)^2 grows like mu^2 (no bound C/x uniform in
        degree), and |eps_mu(mu+1)| >= 1 below the threshold.
Part S  the paper's own route (Stein): lambda d/dlambda I_lambda(a) = -1/2 I_lambda(w . grad a), exact for a Gaussian
        amplitude in one variable; and the band facts: sum_{nu<=M} dim H_nu(S^d) <= 2(M+1)^d, sup|Z_nu| = Z_nu(1).
Part K  the four-kernel envelope in R^2 (Lemma 14.5, numerically): the degree-nu kernel
        int e^{-it rho} chi(rho/N) rho J_nu(r rho) J_nu(s rho) d rho minus its four-phase leading term is
        <= C (1+nu)^2/(N min(r,s)) (rs)^{-1/2} N sum <N(t +- r +- s)>^{-2}, uniformly in nu;
        and the kernel itself is <= (rs)^{-1/2} N sum <N(t +- r +- s)>^{-2} (four-line concentration);
        negative control: in the transition region N r ~ nu the leading term is off by O(1).
Part G  geometry and bookkeeping: 0 <= |s om - x| - (s - x.om) <= |x|^2/(2(s-|x|)) (random); the radial integrals of
        (14.9) exactly; int_0^4 (r+2^-j)^-n r^d dr <= 1/n + log(4 2^j); L(degree)/2^m bounded above and below on the
        support of Q_m (m <= 12); N_j R_k <= 2^{-(k-j-D)-1} for D = ceil(log2(2 H)), exhaustively.
Writes ../results/angular_endpoint.json.
"""

import json
import math
import os
import random

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.special import eval_gegenbauer, jv

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "angular_endpoint.json")


# ----------------------------------------------------------------- Part X (exact exponent algebra)

def part_X():
    n = sp.symbols("n", positive=True)
    d = n - 1
    p = 2 * n / d
    al = d / 2
    A = 2 / p
    N, R = sp.symbols("N R", positive=True)
    eta = al / (2 * n)
    out = {
        "alpha_p_eq_n": sp.simplify(al * p - n) == 0,
        "d_minus_alpha_p": sp.simplify(d - al * p + 1) == 0,
        "alpha_minus_d_over_p": sp.simplify(al - d / p - 1 / p) == 0,
        "alpha_minus_1_over_p": sp.simplify(al - 1 / p - d / p) == 0,
        "n_minus_alpha": sp.simplify(n - al - (al + 1)) == 0,
        "p_minus_2": sp.simplify(p - 2 - 2 / d) == 0,
        # (N^{alpha+1} b^n)^{1-2/p} at b = N^{eta-1} equals N^{-(alpha/2)(1-2/p)}
        "small_source": sp.simplify((al + 1 + n * (eta - 1)) + al / 2) == 0,
        "holder_b_exponent": sp.simplify(n * ((sp.Rational(1, 2) - 1 / p) * (2 / p) + (1 - 1 / p) * (1 - 2 / p))
                                         - n * (1 - 2 / p)) == 0,
        "threshold": sp.simplify(A + sp.Rational(1, 2) - 1 / p - (sp.Rational(1, 2) + 1 / p)) == 0,
        "transplant_error": sp.simplify(R ** (2 + n / p) / (100 * R ** (3 + n / p)) - 1 / (100 * R)) == 0,
        "weak_radial_power": sp.simplify(n / al - p) == 0,
    }
    # numerically for n = 2..8: delta > 0 and p_c values
    out["delta_positive_n2_8"] = all(((k - 1) / 4) * (1 - 2 * (k - 1) / (2 * k)) > 0 for k in range(2, 9))
    return out, all(out.values())


# ----------------------------------------------------------------- Part F (Funk-Hecke)

def part_F():
    mp.mp.dps = 25
    errs, ctrl = [], []
    # S^1: Y_nu = e^{i nu om}
    for nu in range(0, 6):
        for z in (0.7, 3.0, 11.5):
            th = mp.mpf("0.37")
            lhs = mp.quad(lambda om: mp.e ** (1j * z * mp.cos(th - om)) * mp.e ** (1j * nu * om), [0, mp.pi, 2 * mp.pi])
            rhs = 2 * mp.pi * (1j) ** nu * mp.besselj(nu, z) * mp.e ** (1j * nu * th)
            errs.append(abs(lhs - rhs))
            ctrl.append(abs(lhs - 2 * mp.pi * mp.besselj(nu, z) * mp.e ** (1j * nu * th)))
    # S^2: zonal Y_nu = P_nu(e3 . om), evaluate at th with cos(th_polar) = c
    for nu in range(0, 5):
        for z in (0.9, 4.0, 9.0):
            c = mp.mpf("0.3")
            st = mp.sqrt(1 - c ** 2)

            # the azimuthal integral is 2 pi J_0(z sin(th) sin(om)) in closed form
            def integrand(u):
                return mp.e ** (1j * z * c * u) * 2 * mp.pi * mp.besselj(0, z * st * mp.sqrt(1 - u ** 2)) * mp.legendre(nu, u)
            lhs = mp.quad(integrand, [-1, 0, 1])
            rhs = (2 * mp.pi) ** 1.5 * (1j) ** nu * z ** -0.5 * mp.besselj(nu + 0.5, z) * mp.legendre(nu, c)
            errs.append(abs(lhs - rhs))
            ctrl.append(abs(lhs - (2 * mp.pi) ** 1.5 * z ** -0.5 * mp.besselj(nu + 0.5, z) * mp.legendre(nu, c)))
    # two-term form for n = 3, nu = 0 (mu = 1/2, eps = 0): 4 pi sin z / z = (2pi) z^-1 sum_s e^{is(z - pi/2)}
    z = mp.mpf("7.3")
    two_term = abs(4 * mp.pi * mp.sin(z) / z - 2 * mp.pi / z * (mp.e ** (1j * (z - mp.pi / 2)) + mp.e ** (-1j * (z - mp.pi / 2))))
    nonzero_ctrl = [c for c, e in zip(ctrl, errs)]
    out = {"max_error": float(max(errs)), "two_term_n3": float(two_term),
           "control_without_i_nu_min_error_odd_degrees": float(min(c for c in nonzero_ctrl if c > 1e-12))}
    ok = max(errs) < 1e-12 and two_term < 1e-20 and sum(c > 1e-3 for c in ctrl) >= 10
    return out, ok


# ----------------------------------------------------------------- Parts H and E (Hankel's integral, the bound)

def hankel_rep(mu, x, shift=mp.mpf(-0.5)):
    integ = mp.quad(lambda u: mp.e ** (-u) * u ** (mu - 0.5) * (1 + 1j * u / (2 * x)) ** (mu + shift), [0, 1, 10, mp.inf])
    return mp.sqrt(2 / (mp.pi * x)) * mp.e ** (1j * (x - mu * mp.pi / 2 - mp.pi / 4)) / mp.gamma(mu + 0.5) * integ


def part_H():
    mp.mp.dps = 30
    errs, ctrl = [], []
    for mu in (0.25, 0.5, 1, 2.5, 7, 15):
        for x in (0.3, 2, 10, 60):
            h = mp.hankel1(mu, x)
            errs.append(abs(hankel_rep(mu, x) - h) / abs(h))
            ctrl.append(abs(hankel_rep(mu, x, mp.mpf(0.5)) - h) / abs(h))
    # mu = 0: the integrand has u^{-1/2}; check separately with a substitution u = v^2
    for x in (0.3, 2, 10, 60):
        integ = mp.quad(lambda v: 2 * mp.e ** (-v ** 2) * (1 + 1j * v ** 2 / (2 * x)) ** -0.5, [0, 1, 4, mp.inf])
        rep = mp.sqrt(2 / (mp.pi * x)) * mp.e ** (1j * (x - mp.pi / 4)) / mp.gamma(0.5) * integ
        errs.append(abs(rep - mp.hankel1(0, x)) / abs(mp.hankel1(0, x)))
    out = {"max_rel_error": float(max(errs)), "control_min_rel_error": float(min(c for c in ctrl if c > 0) if any(c > 0 for c in ctrl) else 0),
           "control_failures": sum(c > 1e-3 for c in ctrl), "control_cases": len(ctrl)}
    # mu = 1/2 makes the bracket exponent 0 in both forms only for the true formula; the control is wrong elsewhere
    return out, max(errs) < 1e-20 and out["control_failures"] >= len(ctrl) - 6


def eps(mu, x):
    lead = mp.sqrt(2 / (mp.pi * x)) * mp.e ** (1j * (x - mu * mp.pi / 2 - mp.pi / 4))
    return mp.hankel1(mu, x) / lead - 1


def bound(mu, x, k):
    return mp.e ** 1.5 * (2 * mu + 2 * k + 1) ** max(k - 1, 0) * (mu + 1) ** 2 / (2 * x)


def part_E():
    mp.mp.dps = 40
    worst = {}
    for mu in (0, 0.25, 0.5, 1, 1.5, 3, 5.5, 10, 20, 40):
        for f in (1, 1.5, 3, 10, 50):
            x = mp.mpf(f) * (mu + 1) ** 2
            for k in range(4):
                if k == 0:
                    v = abs(eps(mu, x))
                else:
                    v = abs(mp.diff(lambda t: eps(mu, mp.e ** t), mp.log(x), k))
                worst[k] = max(worst.get(k, 0), float(v / bound(mu, x, k)))
    growth = {mu: float((mu + 1) ** 2 * abs(eps(mu, (mu + 1) ** 2))) for mu in (1, 5, 10, 20, 40)}
    below = {mu: float(abs(eps(mu, mu + 1))) for mu in (10, 20, 40)}
    out = {"max_ratio_to_bound_by_k": worst, "x_eps_at_threshold": growth, "eps_below_threshold": below}
    ok = (max(worst.values()) < 1 and growth[40] > 100 * growth[1] and min(below.values()) > 1)
    return out, ok


# ----------------------------------------------------------------- Part S (Stein's route; band facts)

def part_S():
    lam, w = sp.symbols("lambda w", positive=True)
    c = 1 - sp.I * lam / 2  # e^{i lam w^2/2 - w^2} = e^{-c w^2}

    def moment(k):  # int_R w^k e^{-c w^2} dw
        if k % 2:
            return sp.Integer(0)
        return sp.gamma(sp.Rational(k + 1, 2)) / c ** sp.Rational(k + 1, 2)

    def I(poly):  # I_lambda(poly(w) e^{-w^2}) = lambda^{1/2} int e^{i lam w^2/2} poly e^{-w^2}
        P = sp.Poly(sp.expand(poly), w)
        return sp.sqrt(lam) * sum(coef * moment(m[0]) for m, coef in zip(P.monoms(), P.coeffs()))
    a = 1 + w + 3 * w ** 2 - w ** 3 + w ** 4
    lhs = lam * sp.diff(I(a), lam)
    wgrad = sp.expand(sp.exp(w ** 2) * w * sp.diff(a * sp.exp(-w ** 2), w))
    rhs = -sp.Rational(1, 2) * I(wgrad)
    stein = sp.simplify(sp.expand(lhs - rhs)) == 0
    bad = sp.simplify(sp.expand(lhs + rhs)) == 0
    # band counting on S^d: dim H_nu = C(nu+d, d) - C(nu+d-2, d)
    counts = True
    for dd in range(1, 6):
        tot = 0
        for M in range(0, 201):
            dim = math.comb(M + dd, dd) - (math.comb(M + dd - 2, dd) if M >= 2 else 0)
            tot += dim
            counts = counts and tot <= 2 * (M + 1) ** dd
    # sup |C^lam_nu(t)| on [-1,1] is attained at t = 1 (lam > 0): zonal kernels peak on the diagonal
    peak = True
    ts = np.linspace(-1, 1, 4001)
    for lam_g in (0.5, 1.0, 1.5):
        for nu in range(0, 25):
            vals = eval_gegenbauer(nu, lam_g, ts)
            peak = peak and np.max(np.abs(vals)) <= eval_gegenbauer(nu, lam_g, 1.0) * (1 + 1e-12)
    out = {"stein_log_derivative_identity": bool(stein), "control_wrong_sign": bool(bad),
           "band_dimension_bound": bool(counts), "zonal_peak_on_diagonal": bool(peak)}
    return out, stein and not bad and counts and peak


# ----------------------------------------------------------------- Part K (four-kernel envelope, n = 2)

def chi(x):
    y = np.zeros_like(x)
    m = (x > 0.5) & (x < 2)
    y[m] = np.exp(-1 / ((x[m] - 0.5) * (2 - x[m])))
    return y


def kern(nu, N, r, s, ts, lead=False):
    rho = np.linspace(0.5 * N, 2 * N, 40001)
    wgt = chi(rho / N)
    dr = rho[1] - rho[0]
    if lead:
        ph = nu * np.pi / 2 + np.pi / 4
        A = (2 / (np.pi * np.sqrt(r * s))) * np.cos(r * rho - ph) * np.cos(s * rho - ph)
    else:
        A = rho * jv(nu, r * rho) * jv(nu, s * rho)
    return np.array([np.sum(np.exp(-1j * t * rho) * wgt * A) * dr for t in ts])


def part_K():
    N = 200.0
    ts = np.linspace(1, 2, 401)
    ratios, conc = [], []
    for nu in (0, 2, 5, 10):
        for (r, s) in ((0.5, 0.9), (1.3, 0.4), (3.0, 1.6), (0.2, 1.5)):
            k = kern(nu, N, r, s, ts)
            l_ = kern(nu, N, r, s, ts, True)
            env = sum((1 + (N * (ts + a * r + b * s)) ** 2) ** -1 for a in (1, -1) for b in (1, -1))
            th = (1 + nu) ** 2 / (N * min(r, s))
            ratios.append(float(np.max(np.abs(k - l_) * np.sqrt(r * s) / N / th / env)))
            # the kernel itself obeys the envelope (rs)^{-1/2} N sum <N(t +- r +- s)>^{-2}
            conc.append(float(np.max(np.abs(k) * np.sqrt(r * s) / N / env)))
    neg = []
    for nu in (20, 40):
        r = s = (nu + 2) / N
        k = kern(nu, N, r, s, ts)
        l_ = kern(nu, N, r, s, ts, True)
        neg.append(float(np.max(np.abs(k - l_)) / np.max(np.abs(l_))))
    out = {"max_remainder_ratio": max(ratios), "max_kernel_envelope_ratio": max(conc),
           "transition_relative_error": neg}
    return out, max(ratios) < 1 and max(conc) < 1 and min(neg) > 0.3


# ----------------------------------------------------------------- Part G (geometry and bookkeeping)

def part_G(rng):
    geo = True
    for _ in range(20000):
        n = rng.randint(2, 6)
        om = np.array([rng.gauss(0, 1) for _ in range(n)])
        om /= np.linalg.norm(om)
        R = rng.uniform(0.1, 10)
        x = np.array([rng.gauss(0, 1) for _ in range(n)])
        x *= rng.uniform(0, R) / np.linalg.norm(x)
        s = rng.uniform(1.5 * R + 0.01, 1e4)
        diff = np.linalg.norm(s * om - x) - (s - x @ om)
        geo = geo and (-1e-9 <= diff <= np.dot(x, x) / (2 * (s - np.linalg.norm(x))) + 1e-9)
    r, Nn, lam, G = sp.symbols("r N lambda G", positive=True)
    radial = {}
    for n in range(2, 7):
        d = n - 1
        al = sp.Rational(d, 2)
        p = sp.Rational(2 * n, d)
        radial[n] = (sp.simplify(sp.integrate(r ** (d - al * p), (r, 1 / Nn, 4)) - sp.log(4 * Nn)) == 0
                     and sp.simplify(sp.integrate(r ** d, (r, 0, (G / lam) ** (1 / al))) - G ** p / (n * lam ** p)) == 0)
    jint = True
    for n in range(2, 7):
        for j in range(0, 41, 4):
            val = mp.quad(lambda t: (t + mp.mpf(2) ** -j) ** -n * t ** (n - 1), [0, mp.mpf(2) ** -j, 1, 4])
            jint = jint and val <= mp.mpf(1) / n + mp.log(4 * mp.mpf(2) ** j)
    # L(degree) ~ 2^m on supp Q_m = [M_{m-1}, 2 M_m], M_m = 2^{2^m}
    mp.mp.dps = 30
    lo, hi = mp.inf, 0
    for d in (1, 2, 5):
        for m in range(1, 13):
            for nu in (mp.mpf(2) ** (2 ** (m - 1)), 2 * mp.mpf(2) ** (2 ** m)):
                Lv = mp.log(2 + mp.sqrt(1 + nu * (nu + d - 1)))
                lo, hi = min(lo, Lv / 2 ** m), max(hi, Lv / 2 ** m)
    # shell bookkeeping: H = H_m, R_k = H 2^-k, D = ceil(log2(2H)): N_j R_k = H 2^{j-k} <= 2^{-(k-j-D)-1}
    book = True
    for H in [2, 3, 7, 16, 100, 4 ** 5, 10 ** 6]:
        D = math.ceil(math.log2(2 * H))
        for j in range(0, 30):
            for k in range(j + D + 1, j + D + 30):
                book = book and H * 2.0 ** (j - k) <= 2.0 ** (-(k - j - D) - 1) * (1 + 1e-12)
    young = all(sum(2.0 ** (-l * n / (2 * n / (n - 1)))) < 50 for n in range(2, 9) for l in [np.arange(1, 400)])
    out = {"transplant_geometry": bool(geo), "radial_integrals_exact": radial, "j_integral_bound": bool(jint),
           "L_over_2m_range": [float(lo), float(hi)], "shell_bookkeeping": bool(book), "discrete_young": bool(young)}
    ok = geo and all(radial.values()) and jint and 0.3 < lo and hi < 1.6 and book and young
    return out, ok


def main():
    rng = random.Random(1414)
    results = {}
    all_ok = True
    import time
    for name, fn in (("X", part_X), ("F", part_F), ("H", part_H), ("E", part_E), ("S", part_S),
                     ("K", part_K), ("G", lambda: part_G(rng))):
        t0 = time.time()
        res, ok = fn()
        res["seconds"] = round(time.time() - t0, 1)
        results[name] = res
        all_ok = all_ok and ok
        print(name, "ok" if ok else "FAIL", json.dumps(res, default=str))
    results["all_ok"] = bool(all_ok)
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print("ALL OK" if all_ok else "FAILED")


if __name__ == "__main__":
    main()
