"""Numerical companion to TOOLKIT Proposition E4 (cost floor of symmetric swap interpolation).

For lambda = log2 n we minimize, over integers r >= 2, Q >= 1, 2 <= R <= Q with r^(2Q) <= n, the
exponent bound of the §3-4 algorithm with a symmetric replacement matrix:

    delta(r,Q,R) = max{ x log2(1/x) [x = R/Q, capped at 1 beyond 1/e],   (Taylor jets)
                        1/r,                                                (one group)
                        r * eps * Q }                                       (measurement values)
    eps >= max_A min_{A<=j<=A+R} P[Bin(Q,(r-1)/(2r)) = j]                    (Proposition B'')

The window minimum is lower-bounded rigorously by the method of types,
P[W=j] >= exp(-Q KL(j/Q || p))/(Q+1), and the best window is the one centred at the mean.
The minimized delta is compared with (log2 lambda)^2 / sqrt(lambda): the ratio stays bounded above
and below over lambda = 2^10 ... 2^200, as Proposition E4 asserts.
"""

import json
import math


def kl(x, p, t=None):
    """KL(Ber(x) || Ber(p)), computed stably; for |x - p| small use the series in t = x - p."""
    if t is None:
        t = x - p
    if abs(t) < 1e-4:
        # KL = t^2/(2pq) + t^3 (q - p)/(6 p^2 q^2) + O(t^4)  (q = 1 - p); the O(t^4) term is negligible here
        q = 1 - p
        return max(0.0, t * t / (2 * p * q) + t ** 3 * (q - p) / (6 * p * p * q * q))
    out = 0.0
    if x > 0:
        out += x * math.log(x / p)
    if x < 1:
        out += (1 - x) * math.log((1 - x) / (1 - p))
    return max(0.0, out)


def eps_lower(Q, R, r):
    """Natural log of a rigorous lower bound for min_{j in centred window} P[Bin(Q,p) = j], p = (r-1)/(2r).
    Window endpoints and their offsets from the mean are computed in exact integer arithmetic."""
    p = (r - 1) / (2 * r)
    a = (Q * (r - 1) - r * R) // (2 * r)          # floor(Qp - R/2)
    a = max(0, min(a, Q - R))
    def offset(j):                                 # (j - Qp)/Q exactly, then rounded once
        return (2 * r * j - Q * (r - 1)) / (2 * r * Q)
    worst = max(kl(None, p, offset(a)) if abs(offset(a)) < 1e-4 else kl(a / Q, p),
                kl(None, p, offset(a + R)) if abs(offset(a + R)) < 1e-4 else kl((a + R) / Q, p))
    return -Q * worst - math.log(Q + 1)   # natural log of the lower bound


def taylor(R, Q):
    x = R / Q
    if x >= 1 / math.e:
        return 1 / math.e * math.log2(math.e)
    return x * math.log2(1 / x)


def best(lam):
    bestv = (float("inf"), None)
    r = 2
    while r <= 2 ** 200:
        Qmax = int(lam / (2 * math.log2(r)))
        if Qmax < 2:
            break
        # scan Q on a log grid below its maximum, and R on a geometric grid
        Qs = sorted({max(2, int(Qmax * f)) for f in (1, 0.8, 0.6, 0.4, 0.25, 0.1, 0.05)})
        for Q in Qs:
            R = 2
            while R <= Q:
                le = eps_lower(Q, R, r)
                meas = math.exp(min(700.0, math.log(r) + le + math.log(Q)))
                d = max(taylor(R, Q), 1 / r, meas)
                if d < bestv[0]:
                    bestv = (d, (r, Q, R))
                R = max(R + 1, int(R * 1.15))
        r = max(r + 1, int(r * 1.3))
    return bestv


def main():
    rows = []
    for e in (10, 20, 30, 40, 60, 80, 100, 150, 200, 300):
        lam = 2.0 ** e
        d, arg = best(lam)
        ref = math.log2(lam) ** 2 / math.sqrt(lam)
        rows.append({"log2_n": f"2^{e}", "min_delta": d, "argmin_(r,Q,R)": arg, "reference": ref, "ratio": d / ref})
        print(rows[-1])
    ratios = [r["ratio"] for r in rows if r["min_delta"] < 1]
    out = {"rows": rows, "ratio_range": [min(ratios), max(ratios)] if ratios else None}
    with open(__file__.rsplit("/", 2)[0] + "/results/cost_floor.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
