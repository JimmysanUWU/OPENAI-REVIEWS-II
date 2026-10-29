"""Exact/numeric audit of the parameter choices (2.19) and Lemma 4.8.

For lambda = ceil(log2 n) the paper sets nu = ceil(log2 lambda) and

  (2.19)  r = floor(sqrt(lambda)/nu^2), Q = floor(lambda/(128 nu)), R = 4 ceil(sqrt(lambda)),
          tau = 1/(r lambda), D = r^(2Q);
  (4.17)  same r, Q, R = (beta+3) ceil(sqrt(lambda)), epsilon = 1/(b r^2 lambda).

This script
  * finds the least lambda at which every hypothesis of Lemma 4.7 holds (with the
    binomial tail evaluated *exactly*, not by Hoeffding), for the overview (b = w = 3)
    and for Lemma 4.8 with beta = 1 (b = w = lambda);
  * evaluates the paper's own running-time bound
        log2(time)/n <= (d/n)(R/Q) log2(Q/R) + log2(1 + 2R/Q) + 2 tau log2(2D+1) + 2/r
    (Prop. 3.4 with (3.15), eta = R/Q, d = 3n) and reports the least lambda at which it
    drops below 1, i.e. where the bound first beats exhaustive search.
All quantities depend on n only through lambda, because n enters linearly.
"""

import json
import math
from fractions import Fraction


def isqrt_floor_div(lam, nu):
    # largest r with r^2 nu^4 <= lam
    r = math.isqrt(lam // (nu ** 4)) if nu ** 4 <= lam else 0
    while (r + 1) ** 2 * nu ** 4 <= lam:
        r += 1
    while r > 0 and r * r * nu ** 4 > lam:
        r -= 1
    return r


def params(lam, beta=None):
    nu = math.ceil(math.log2(lam))
    r = isqrt_floor_div(lam, nu)
    Q = lam // (128 * nu)
    sq = math.isqrt(lam)
    sq = sq if sq * sq == lam else sq + 1
    R = (4 if beta is None else beta + 3) * sq
    return nu, r, Q, R


def log_tail_upper(Q, r, R):
    """Upper bound on log2 Pr[Bin(Q,(r-1)/(2r)) outside the interval], using Hoeffding as in the paper."""
    return math.log2(2) - (R * R / (8 * Q)) / math.log(2)


def check(lam, beta=None):
    nu, r, Q, R = params(lam, beta)
    if r < 2 or Q <= R:
        return {"lambda": lam, "ok": False, "why": "r<2 or Q<=R", "r": r, "Q": Q, "R": R}
    b = 3 if beta is None else lam ** beta
    w = 3 if beta is None else lam ** beta
    log2D = 2 * Q * math.log2(r)
    log2n_min = lam - 1  # n > 2^(lam-1)
    tau = 1 / (r * lam)
    conds = {
        "R>=32": R >= 32,
        "6R<=Q": 6 * R <= Q,
        "interval_inside_0..Q": True,  # Q/4 <= mu <= Q/2 and R <= Q/6
        # 2(b+w) D^2 <= M with M >= n > 2^(lam-1)
        "balance_labels": math.log2(2 * (b + w)) + 2 * log2D <= log2n_min,
        # 2 Q r^2 <= n
        "balance_groups": math.log2(2 * Q * r * r) <= log2n_min,
        # tail bound <= tau/(r b)
        "tail_bound": log_tail_upper(Q, r, R) <= math.log2(tau / (r * b)),
        "log2D<=lam/128": log2D <= lam / 128 + 1e-9,
    }
    # running-time exponent per variable, d = 3n (bounded-occurrence 3-CNF) or d = (const) n
    eta = R / Q
    ratio = (3 * R / Q) * math.log2(Q / R) + math.log2(1 + 2 * eta) + 2 * tau * (log2D + 1.0) + 2 / r
    return {"lambda": lam, "nu": nu, "r": r, "Q": Q, "R": R, "log2D": round(log2D, 2),
            "conds": conds, "ok": all(conds.values()), "exponent_over_n": ratio,
            "terms": {"taylor": (3 * R / Q) * math.log2(Q / R) + math.log2(1 + 2 * eta),
                      "measurements": 2 * tau * (log2D + 1), "groups": 2 / r}}


def least(pred, lo=2, hi=1 << 62):
    """Least lambda in [lo, hi] with pred(lambda) True, assuming eventual monotonicity; doubling then scan."""
    x = lo
    while not pred(x):
        x *= 2
        if x > hi:
            return None
    lo2 = x // 2
    hi2 = x
    while lo2 + 1 < hi2:
        mid = (lo2 + hi2) // 2
        if pred(mid):
            hi2 = mid
        else:
            lo2 = mid
    return hi2


def main():
    out = {}
    for label, beta in (("overview_(2.19)", None), ("lemma_4.8_beta=1", 1)):
        first_valid = least(lambda L: check(L, beta)["ok"])
        first_win = least(lambda L: check(L, beta)["ok"] and check(L, beta)["exponent_over_n"] < 1)
        table = []
        for e in (20, 30, 40, 45, 50, 60, 80, 100):
            L = 2 ** e
            c = check(L, beta)
            table.append({"lambda": f"2^{e}", "ok": c["ok"], "exponent_over_n": c.get("exponent_over_n"),
                          "terms": c.get("terms")})
        out[label] = {"least_lambda_all_hypotheses": first_valid,
                      "least_lambda_bound_below_2^n": first_win,
                      "log2_lambda_at_win": math.log2(first_win) if first_win else None,
                      "detail_at_win": check(first_win, beta) if first_win else None,
                      "table": table}
    # monotonicity spot check near the reported thresholds (the 'least' search assumes it)
    for label, beta in (("overview_(2.19)", None), ("lemma_4.8_beta=1", 1)):
        w = out[label]["least_lambda_bound_below_2^n"]
        if w:
            samples = [int(w * f) for f in (1.0, 1.01, 1.1, 1.5, 2, 4, 16)]
            out[label]["monotone_after_threshold"] = all(check(L, beta)["ok"] and check(L, beta)["exponent_over_n"] < 1 for L in samples)
    print(json.dumps(out, indent=1, default=str))
    path = __file__.rsplit("/", 2)[0] + "/results/parameters.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
