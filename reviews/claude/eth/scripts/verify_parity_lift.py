"""Exact checks of Lemma 5.1, Lemma 5.2 (cross-equation truncation) and Lemma B.1.

  L5.1  p_l(Y) = (1-Y)^l sum_{j<l} binom(l+j-1, j) Y^j satisfies p_l(y) = [2|y] mod 2^l,
        and division by h(Y) = Y(Y-1) gives p_l = sum_{j<l} e_j(Y) h(Y)^j with deg e_j <= 1;
  L5.2  for integer y_1..y_m: S_l(prod_i B_i(V) mod V^l) = prod_i [2|y_i] mod 2^l, where
        B_i(V) = sum_j e_j(y_i) h(y_i)^j V^j;
        and, symbolically, every retained V^s-coefficient has degree <= m + 2(l-1) in the y_i
        (d_* = 1), while the untruncated product has degree up to m(2l-1);
  ex.   p_2(Y) = 1 - 3Y^2 + 2Y^3 = (1-Y) + (2Y-1) h(Y)  (the worked example after Lemma 5.2);
  B.1   |a_j|, |b_j| for e_j = a_j + b_j Y, reported against 2^{C l^2}.
"""

import json
import random
import sys
from fractions import Fraction
from itertools import product
from math import comb

sys.path.insert(0, __file__.rsplit("/", 1)[0])


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def ev(p, y):
    acc = 0
    for c in reversed(p):
        acc = acc * y + c
    return acc


def p_ell(l):
    s = [comb(l + j - 1, j) for j in range(l)]
    one_minus = [1]
    for _ in range(l):
        one_minus = pmul(one_minus, [1, -1])
    return pmul(one_minus, s)


def divide_by_h(p):
    """Write p = sum_j e_j h^j with deg e_j <= 1, h = Y^2 - Y (monic division)."""
    h = [0, -1, 1]
    es = []
    cur = list(p)
    while len(cur) > 2 or (len(cur) == 2 and False):
        # quotient and remainder of cur by h
        quot = [0] * max(1, len(cur) - 2)
        rem = list(cur)
        for k in range(len(rem) - 1, 1, -1):
            c = rem[k]
            if c:
                quot[k - 2] += c
                for i, hv in enumerate(h):
                    rem[k - 2 + i] -= c * hv
        es.append(rem[:2] + [0] * (2 - len(rem[:2])))
        cur = quot
        while len(cur) > 1 and cur[-1] == 0:
            cur.pop()
    es.append(cur + [0] * (2 - len(cur)))
    return es


def main():
    out = {"L5.1": {}, "L5.2": {}, "B.1": {}}
    ok_all = True
    for l in range(1, 11):
        p = p_ell(l)
        mod = 2 ** l
        ok_cong = all((ev(p, y) - (1 if y % 2 == 0 else 0)) % mod == 0 for y in range(-200, 201))
        es = divide_by_h(p)
        es = es[:l] + [[0, 0]] * max(0, l - len(es))
        h = [0, -1, 1]
        recon = [0]
        hp = [1]
        for j, e in enumerate(es):
            recon = padd(recon, pmul(e, hp))
            hp = pmul(hp, h)
        while len(recon) > 1 and recon[-1] == 0:
            recon.pop()
        pp = list(p)
        while len(pp) > 1 and pp[-1] == 0:
            pp.pop()
        ok_div = recon == pp and len([e for e in divide_by_h(p) if any(e)]) <= l
        maxcoef = max(abs(c) for e in es for c in e)
        out["L5.1"][l] = {"congruence": ok_cong, "decomposition": ok_div, "deg_p": len(pp) - 1}
        out["B.1"][l] = {"max_abs_e_coeff": maxcoef, "log2": maxcoef.bit_length(), "l^2": l * l}
        ok_all &= ok_cong and ok_div
    # worked example
    p2 = p_ell(2)
    out["example_p2"] = p2 == [1, 0, -3, 2] and divide_by_h(p2)[:2] == [[1, -1], [-1, 2]]
    ok_all &= out["example_p2"]

    # Lemma 5.2 numerically for random integer tuples
    rng = random.Random(0)
    for l in range(1, 8):
        es = divide_by_h(p_ell(l))
        es = es[:l] + [[0, 0]] * max(0, l - len(es))
        mod = 2 ** l
        ok = True
        for _ in range(300):
            m = rng.randint(0, 6)
            ys = [rng.randint(-50, 50) for _ in range(m)]
            prodV = [1]
            for y in ys:
                hy = y * (y - 1)
                B = [ev(es[j], y) * hy ** j for j in range(l)]
                prodV = pmul(prodV, B)[:l]
            S = sum(prodV)
            target = 1
            for y in ys:
                target *= 1 if y % 2 == 0 else 0
            ok &= (S - target) % mod == 0
        out["L5.2"][l] = {"congruence_random_tuples": ok}
        ok_all &= ok

    # Lemma 5.2 degree bound: symbolic in y_1..y_m for small m, l (track total degree only)
    deg_ok = True
    deg_info = {}
    for l in range(1, 5):
        es = divide_by_h(p_ell(l))
        es = es[:l] + [[0, 0]] * max(0, l - len(es))
        for m in range(1, 4):
            # each B_i coefficient of V^j has degree deg(e_j) + 2j in y_i
            degs = [[(1 if es[j][1] != 0 else 0) + 2 * j for j in range(l)] for _ in range(m)]
            maxdeg = -1
            for js in product(range(l), repeat=m):
                if sum(js) < l:
                    maxdeg = max(maxdeg, sum(degs[i][js[i]] for i in range(m)))
            bound = m + 2 * (l - 1)
            deg_info[f"l={l},m={m}"] = {"retained_max_degree": maxdeg, "bound_m+2(l-1)": bound,
                                         "untruncated_bound": m * (2 * l - 1)}
            deg_ok &= maxdeg <= bound
    out["L5.2_degree"] = deg_info
    out["L5.2_degree_ok"] = deg_ok
    ok_all &= deg_ok
    out["all_ok"] = ok_all
    print(json.dumps(out, indent=1))
    print("ALL PASS" if ok_all else "SOME CHECK FAILED")
    path = __file__.rsplit("/", 2)[0] + "/results/parity_lift.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
