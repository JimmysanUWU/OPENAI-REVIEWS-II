"""Exact checks of the p-adic cross-equation truncation (TOOLKIT Theorem E2), a generalization of
Lemmas 5.1-5.2 from F_2 to F_p.

  F_l(U)   = (1-U)^l sum_{j<l} binom(l+j-1, j) U^j
  p_l(Y)   = F_l(Y^(p-1))                        claim: p_l(y) = [p | y]  (mod p^l) for all integers y
  h(Y)     = Y^p - Y                              h(y) = 0 (mod p) for all integers y (Fermat)
  p_l      = sum_j e_j(Y) h(Y)^j, deg e_j < p      (repeated division by the monic h)
  B_i(V)   = sum_{j<l} e_j(y_i) h(y_i)^j V^j
  claim:     S_l( prod_i B_i(V) mod V^l ) = prod_i [p | y_i]  (mod p^l)
  degree:    every retained coefficient has degree <= m(p-1) + p(l-1) in (y_1..y_m)
             (versus m (p-1)(2l-1) for the untruncated product of the p_l(y_i)).
Also checks the identity F_l(U) = 1 - U^l sum_{j<l} binom(l+j-1,j)(1-U)^j and p = 2 agrees with Lemma 5.1.
"""

import json
import random
from itertools import product
from math import comb


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
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


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def F(l):
    s = [comb(l + j - 1, j) for j in range(l)]
    om = [1]
    for _ in range(l):
        om = pmul(om, [1, -1])
    return pmul(om, s)


def compose_power(Fp, k):
    """F(Y^k) as a polynomial in Y."""
    out = [0] * ((len(Fp) - 1) * k + 1)
    for i, c in enumerate(Fp):
        out[i * k] += c
    return out


def p_l(p, l):
    return trim(compose_power(F(l), p - 1))


def divide_by_h(poly, p):
    """poly = sum_j e_j h^j with h = Y^p - Y, deg e_j < p."""
    h = [0, -1] + [0] * (p - 2) + [1]
    es = []
    cur = trim(poly)
    while True:
        if len(cur) <= p:
            es.append(cur + [0] * (p - len(cur)))
            break
        quot = [0] * (len(cur) - p)
        rem = list(cur)
        for k in range(len(rem) - 1, p - 1, -1):
            c = rem[k]
            if c:
                quot[k - p] += c
                for i, hv in enumerate(h):
                    rem[k - p + i] -= c * hv
        rem = trim(rem[:p])
        es.append(rem + [0] * (p - len(rem)))
        cur = trim(quot)
    return es


def check(p, l, rng, trials=300):
    mod = p ** l
    pl = p_l(p, l)
    ok_cong = all((ev(pl, y) - (1 if y % p == 0 else 0)) % mod == 0 for y in range(-400, 401))
    es = divide_by_h(pl, p)
    h = [0, -1] + [0] * (p - 2) + [1]
    recon = [0]
    hp = [1]
    for e in es:
        recon = padd(recon, pmul(e, hp))
        hp = pmul(hp, h)
    ok_decomp = trim(recon) == pl and all(len(trim(e)) <= p for e in es)
    es_l = (es + [[0] * p] * l)[:l]
    ok_trunc = True
    for _ in range(trials):
        m = rng.randint(0, 6)
        ys = [rng.randint(-60, 60) for _ in range(m)]
        prodV = [1]
        for y in ys:
            hy = y ** p - y
            B = [ev(es_l[j], y) * hy ** j for j in range(l)]
            prodV = pmul(prodV, B)[:l]
        target = 1
        for y in ys:
            target *= 1 if y % p == 0 else 0
        ok_trunc &= (sum(prodV) - target) % mod == 0
    # degree bound: coefficient of V^j in B_i has degree deg(e_j) + p j
    degs = [(len(trim(e)) - 1 if any(e) else -10 ** 9) + p * j for j, e in enumerate(es_l)]
    deg_rows = {}
    ok_deg = True
    for m in (1, 2, 3, 4):
        mx = -1
        for js in product(range(l), repeat=m):
            if sum(js) < l:
                mx = max(mx, sum(degs[j] for j in js))
        bound = m * (p - 1) + p * (l - 1)
        ok_deg &= mx <= bound
        deg_rows[m] = {"retained_max_degree": mx, "bound": bound, "untruncated": m * (len(pl) - 1)}
    return {"p": p, "l": l, "congruence": ok_cong, "decomposition": ok_decomp,
            "cross_equation_truncation": ok_trunc, "degree_bound": ok_deg, "degrees": deg_rows}


def main():
    rng = random.Random(2)
    # identity F_l(U) = 1 - U^l sum binom(l+j-1,j)(1-U)^j
    ok_identity = True
    for l in range(1, 9):
        acc = [0]
        om = [1]
        for j in range(l):
            acc = padd(acc, [comb(l + j - 1, j) * c for c in om])
            om = pmul(om, [1, -1])
        rhs = padd([1], [-c for c in pmul([0] * l + [1], acc)])
        ok_identity &= trim(rhs) == trim(F(l))
    rows = []
    for p in (2, 3, 5, 7):
        for l in range(1, 7 if p <= 3 else 5):
            rows.append(check(p, l, rng))
    ok = ok_identity and all(r["congruence"] and r["decomposition"] and r["cross_equation_truncation"] and r["degree_bound"]
                             for r in rows)
    out = {"identity_F": ok_identity, "rows": rows, "all_ok": ok}
    for r in rows:
        print(r["p"], r["l"], r["congruence"], r["decomposition"], r["cross_equation_truncation"], r["degree_bound"],
              r["degrees"][3])
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/padic_truncation.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
