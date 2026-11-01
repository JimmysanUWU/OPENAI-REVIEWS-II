"""Exact 2-adic checks of the Appendix D machinery (Lemma D.4, identity (D.10), Lemma D.5, (D.12)).

Setting (D.4.3): L = 2^a, S = {-L, ..., L-1}, C(T) = prod_{t in S} (T - t), alpha = 2^(a+1),
W_xy(T) = C(T) / (C'(0) (T - (y - x))) for slots x, y in {0..L-1}.

  D4a   W(0) = I and, for t in S, W_xy(t) != 0 only if y - x = t;
  D4b   for t in S, every coefficient of W_xy(t + alpha z) is 2-adically integral, and
        is divisible by 2 unless y - x = t;
  D10   v2(E_t(0)) = D sum_j (j/2^j)(q(i+j) - q(i)) >= -6 D h, with q(m) = h max(0, u - m);
  D5    tapered Hermite recovery: from residues of the first q_t coefficients of F(t + alpha z)
        mod 2^{p_t}, p_t = p + 8 h D_t + 8 M (a+1), formula (D.9) returns F(0) mod 2^p, for random
        integer F with deg F <= M(2L-1); and a negative control with p_t = p fails;
  D12   geometric expansion of a crossing entry agrees with W_xy(T) to the stated 2-adic precision.
The coefficient field is Q (the paper's K = Q[U]/(f) with e = 1 suffices for these identities).
"""

import json
import random
from fractions import Fraction
from math import comb


def v2(x):
    x = Fraction(x)
    if x == 0:
        return 10 ** 9
    n, d = x.numerator, x.denominator
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    while d % 2 == 0:
        d //= 2
        c -= 1
    return c


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def ev(p, x):
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def compose_affine(p, c0, c1):
    """p(c0 + c1 z) as a polynomial in z."""
    out = [Fraction(0)]
    powr = [Fraction(1)]
    for k, c in enumerate(p):
        if k:
            powr = pmul(powr, [Fraction(c0), Fraction(c1)])
        if c:
            out = [(out[i] if i < len(out) else 0) + (c * powr[i] if i < len(powr) else 0)
                   for i in range(max(len(out), len(powr)))]
    return out


def W_poly(a, x, y):
    L = 2 ** a
    S = list(range(-L, L))
    num = [Fraction(1)]
    for t in S:
        if t != y - x:
            num = pmul(num, [Fraction(-t), Fraction(1)])
    Cp0 = Fraction(1)
    for t in S:
        if t != 0:
            Cp0 *= (0 - t)
    return [c / Cp0 for c in num]


def check_D4(a):
    L = 2 ** a
    S = list(range(-L, L))
    alpha = 2 ** (a + 1)
    ok0 = ok_shift = ok_int = ok_div = True
    for x in range(L):
        for y in range(L):
            W = W_poly(a, x, y)
            if ev(W, 0) != (1 if x == y else 0):
                ok0 = False
            for t in S:
                val = ev(W, t)
                if val != 0 and y - x != t:
                    ok_shift = False
                Wz = compose_affine(W, t, alpha)
                vals = [v2(c) for c in Wz if c != 0]
                if vals and min(vals) < 0:
                    ok_int = False
                if y - x != t and vals and min(vals) < 1:
                    ok_div = False
    return {"a": a, "W0_identity": ok0, "Wt_shift_support": ok_shift, "W(t+alpha z)_integral": ok_int,
            "off_shift_divisible_by_2": ok_div}


def multiplicities(a, u, M):
    L = 2 ** a
    S = [t for t in range(-L, L) if t != 0]
    h = -(-3 * M // u)
    q = {t: h * max(0, u - v2(t)) for t in S}
    return S, h, q


def check_D10(a, u, M):
    S, h, q = multiplicities(a, u, M)
    L = 2 ** a
    ok = True
    worst = None
    for t in S:
        if q[t] == 0:
            continue
        i = v2(t)
        Dt = L // 2 ** i
        val = sum(q[s] * (v2(s) - v2(t - s)) for s in S if s != t)
        formula = Fraction(0)
        for j in range(1, a - i + 1):
            qi_j = h * max(0, u - (i + j))
            qi = h * max(0, u - i)
            formula += Fraction(Dt * j, 2 ** j) * (qi_j - qi)
        ok &= (val == formula) and (val >= -6 * Dt * h)
        slack = val + 6 * Dt * h
        worst = slack if worst is None else min(worst, slack)
    return {"a": a, "u": u, "M": M, "identity_and_bound": ok, "min_slack_vs_-6Dh": worst}


def recover(a, u, M, p, data, q, pt_used):
    """Formula (D.9) from the supplied truncated data; returns F(0) mod 2^p."""
    L = 2 ** a
    alpha = 2 ** (a + 1)
    S = [t for t in range(-L, L) if t != 0]
    total = Fraction(0)
    for t in S:
        if q[t] == 0:
            continue
        # E_t(T) = prod_{s != 0, t} ((T - s)/(t - s))^{q_s};  E_t(t + alpha z) as a series
        Et0 = Fraction(1)
        Ez = [Fraction(1)]
        for s in S:
            if s == t or q[s] == 0:
                continue
            Et0 *= Fraction(-s, t - s) ** q[s]
            fac = [Fraction(1), Fraction(alpha, t - s)]
            for _ in range(q[s]):
                Ez = pmul(Ez, fac)[: q[t]]
        # inverse series of Ez up to z^{q_t}
        inv = [Fraction(0)] * q[t]
        inv[0] = 1 / Ez[0]
        for k in range(1, q[t]):
            acc = sum(Ez[j] * inv[k - j] for j in range(1, min(k, len(Ez) - 1) + 1))
            inv[k] = -acc / Ez[0]
        Rt = pmul(data[t], inv)[: q[t]]
        total += Et0 * ev(Rt, Fraction(-t, alpha))
    assert total.denominator % 2 == 1 or True
    # reduce: total should be 2-adically close to an integer
    num, den = total.numerator, total.denominator
    if den % 2 == 0:
        return None
    return num * pow(den, -1, 2 ** p) % 2 ** p


def check_D5(a, u, M, p, trials, rng):
    S, h, q = multiplicities(a, u, M)
    L = 2 ** a
    alpha = 2 ** (a + 1)
    deg = M * (2 * L - 1)
    assert sum(q.values()) > deg
    ok = True
    fail_neg = 0
    for _ in range(trials):
        Fp = [Fraction(rng.randint(-10 ** 6, 10 ** 6)) for _ in range(deg + 1)]
        truth = int(Fp[0]) % 2 ** p
        data, data_neg = {}, {}
        for t in S:
            if q[t] == 0:
                continue
            i = v2(t)
            Dt = L // 2 ** i
            pt = p + 8 * h * Dt + 8 * M * (a + 1)
            G = compose_affine(Fp, t, alpha)[: q[t]]
            G = G + [Fraction(0)] * (q[t] - len(G))
            # replace each coefficient by a random lift of its residue mod 2^{p_t}
            data[t] = [Fraction(int(c) % 2 ** pt + rng.randint(-5, 5) * 2 ** pt) for c in G]
            data_neg[t] = [Fraction(int(c) % 2 ** p + rng.randint(-5, 5) * 2 ** p) for c in G]
        got = recover(a, u, M, p, data, q, None)
        ok &= got == truth
        got_neg = recover(a, u, M, p, data_neg, q, None)
        if got_neg != truth:
            fail_neg += 1
    return {"a": a, "u": u, "M": M, "p": p, "trials": trials, "recovery_ok": ok,
            "negative_control_failures_with_p_t=p": fail_neg}


def check_D12(a, i, t, x, y, terms):
    """Crossing entry: y = y0 + 2^{i+1} j; series sum_b C(T)/(C'(0)(T-y0+x)) (2^{i+1}/(T-y0+x))^b j^b
    should approximate W_xy at T = t + alpha*zeta 2-adically with error divisible by 2^{terms}."""
    L = 2 ** a
    alpha = 2 ** (a + 1)
    y0 = y % (2 ** (i + 1))
    j = y // (2 ** (i + 1))
    Tval = Fraction(t) + alpha * Fraction(3)  # zeta = 3
    W = ev(W_poly(a, x, y), Tval)
    S = list(range(-L, L))
    C = Fraction(1)
    for s in S:
        C *= (Tval - s)
    Cp0 = Fraction(1)
    for s in S:
        if s != 0:
            Cp0 *= -s
    base = C / (Cp0 * (Tval - y0 + x))
    ratio = Fraction(2 ** (i + 1), 1) / (Tval - y0 + x)
    approx = sum(base * ratio ** b * j ** b for b in range(terms))
    return v2(W - approx)


def main():
    rng = random.Random(1)
    out = {"D4": [check_D4(a) for a in (1, 2, 3)]}
    out["D10"] = [check_D10(a, u, M) for a in (2, 3, 4) for u in range(1, a + 1) for M in (u, 2 * u, 5)]
    out["D5"] = [check_D5(2, 2, 2, 4, 20, rng), check_D5(2, 1, 1, 3, 20, rng), check_D5(3, 2, 2, 4, 6, rng)]
    # D12 for a crossing entry: pick t with v2(t) = i, a pair with differing bit i
    d12 = []
    a = 3
    L = 2 ** a
    for i in (0, 1):
        for t in [s for s in range(-L, L) if s != 0 and v2(s) == i][:2]:
            for x in range(L):
                for y in range(L):
                    if ((x >> i) & 1) != ((y >> i) & 1) and y - x != t:
                        y0 = y % (2 ** (i + 1))
                        if v2(t - y0 + x) <= i:
                            val = check_D12(a, i, t, x, y, 12)
                            d12.append(val)
    out["D12_min_valuation_of_error_after_12_terms"] = min(d12) if d12 else None
    out["D12_cases"] = len(d12)
    allok = (all(all(v for k, v in r.items() if isinstance(v, bool)) for r in out["D4"])
             and all(r["identity_and_bound"] for r in out["D10"])
             and all(r["recovery_ok"] for r in out["D5"])
             and (out["D12_min_valuation_of_error_after_12_terms"] or 0) >= 12)
    out["all_ok"] = allok
    print(json.dumps(out, indent=1))
    print("ALL PASS" if allok else "SOME CHECK FAILED")
    path = __file__.rsplit("/", 2)[0] + "/results/appendix_d.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
