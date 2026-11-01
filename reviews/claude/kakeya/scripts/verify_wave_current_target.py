"""Direct check of Proposition 2.2 in the target coordinates (companion to 08_GAPLESS_PROOF.md).

The paper obtains the wave current A^0 = chi chi^* by pushing the lifted current of Lemma 2.1
forward through the holomorphic flux map.  Here we bypass that derivation and verify its end
product directly, in null coordinates, for random finite superpositions (n = 2, 3):

  (C)  every real column of (chi chi^*)^R has zero divergence in the 3(n+1) real coordinates
       (X^mu, Re omega_mu, Im_k omega_mu), with omega = (Y, P, S1) and slot triples (1, i, j);
  (T)  (dPhi) chi = 0 for dPhi = l(W) = (-W.W/2, W, 1)                         [tangency];
  (S)  A^0_TT = |chi_T|^2 = |S_0 u|^2 + |S_1 u|^2, chi_T = (chi_t - chi_s0)/sqrt 2     [(2.15)];
  (N)  negative controls: dropping the rotation iota_e (e -> 1), or the volume factor y^{-d/2},
       destroys conservation.

chi_mu = sum_modes V^v_mu iota_e(conj u_{kappa,v}) m_{kappa,v},
V^v = (1, (v+W)/2, -W.v/2),  W = P/Y,  e = Y/|Y|,  iota_e(a + i b) = a + e i b,
m = (kappa/y)^{d/2} exp(-(kappa y/2)|W - v|^2),  u = c exp(i kappa (s0 + v.x - |v|^2 t/2)).

Derivatives are central differences in 60-digit arithmetic (step 1e-22), so the reported
divergences are exact to ~1e-40 relative to the tensor size.  Writes ../results/wave_current_target.json.
"""

import json
import random

import mpmath as mp

mp.mp.dps = 60


# ----------------------------------------------------------------- quaternions (1, i, j, k)

def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def qconj(p):
    return (p[0], -p[1], -p[2], -p[3])


def qadd(p, q):
    return tuple(x + y for x, y in zip(p, q))


def qscale(c, p):
    return tuple(c * x for x in p)


E = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]          # e_0, e_1, e_2 = 1, i, j


# ----------------------------------------------------------------- C_k numbers a + k b as (a, b)

def kmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def kdiv(x, y):
    den = y[0] ** 2 + y[1] ** 2
    return ((x[0] * y[0] + x[1] * y[1]) / den, (x[1] * y[0] - x[0] * y[1]) / den)


def kq(x):
    return (x[0], mp.mpf(0), mp.mpf(0), x[1])


def chi(coords, d, modes, rotate=True, volume=True):
    """chi_mu for mu in (t, x_1..x_d, s0); coords = (X^mu)_mu + (alpha_mu, beta_mu)_mu."""
    nsl = d + 2
    X = coords[:nsl]
    om = [(coords[nsl + 2 * m], coords[nsl + 2 * m + 1]) for m in range(nsl)]
    Y, P = om[0], om[1:1 + d]
    y = mp.sqrt(Y[0] ** 2 + Y[1] ** 2)
    e = (Y[0] / y, Y[1] / y) if rotate else (mp.mpf(1), mp.mpf(0))
    W = [kdiv(Pj, Y) for Pj in P]
    out = [(mp.mpf(0),) * 4 for _ in range(nsl)]
    for kap, v, c in modes:
        Wv = (sum(W[j][0] * v[j] for j in range(d)), sum(W[j][1] * v[j] for j in range(d)))
        V = [(mp.mpf(1), mp.mpf(0))] + [((v[j] + W[j][0]) / 2, W[j][1] / 2) for j in range(d)] \
            + [(-Wv[0] / 2, -Wv[1] / 2)]
        dist2 = sum((W[j][0] - v[j]) ** 2 + W[j][1] ** 2 for j in range(d))
        m = (kap / y) ** (mp.mpf(d) / 2) if volume else mp.mpf(1)
        m *= mp.exp(-kap * y * dist2 / 2)
        phase = kap * (X[-1] + sum(v[j] * X[1 + j] for j in range(d)) - X[0] * sum(x * x for x in v) / 2)
        ubar = mp.conj(c * mp.expj(phase))
        iota = (mp.re(ubar), mp.im(ubar) * e[0], mp.im(ubar) * e[1], mp.mpf(0))   # a + e i b
        for mu in range(nsl):
            out[mu] = qadd(out[mu], qscale(m, qmul(kq(V[mu]), iota)))
    return out


def real_tensor(coords, d, modes, **kw):
    ch = chi(coords, d, modes, **kw)
    nsl = d + 2
    R = [[None] * (3 * nsl) for _ in range(3 * nsl)]
    for mu in range(nsl):
        for nu in range(nsl):
            A = qmul(ch[mu], qconj(ch[nu]))
            for a in range(3):
                for b in range(3):
                    R[3 * mu + a][3 * nu + b] = qmul(qmul(qconj(E[a]), A), E[b])[0]
    return R


def coord_index(nu, beta, nsl):
    """position in coords of the beta-th real coordinate of slot nu (beta: 0 real, 1 alpha, 2 beta)."""
    return nu if beta == 0 else nsl + 2 * nu + (beta - 1)


def divergence(coords, d, modes, h=mp.mpf("1e-22"), **kw):
    nsl = d + 2
    N = 3 * nsl
    div = [mp.mpf(0)] * N
    for nu in range(nsl):
        for beta in range(3):
            k = coord_index(nu, beta, nsl)
            cp, cm = list(coords), list(coords)
            cp[k] += h
            cm[k] -= h
            Rp, Rm = real_tensor(cp, d, modes, **kw), real_tensor(cm, d, modes, **kw)
            col = 3 * nu + beta
            for row in range(N):
                div[row] += (Rp[row][col] - Rm[row][col]) / (2 * h)
    R0 = real_tensor(coords, d, modes, **kw)
    scale = max(abs(x) for r in R0 for x in r)
    return max(abs(x) for x in div), scale


def random_modes(d, rng, n_modes):
    modes = []
    for _ in range(n_modes):
        kap = mp.mpf(rng.uniform(0.6, 1.6))
        v = [mp.mpf(rng.uniform(-0.3, 0.3)) for _ in range(d)]
        c = mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1))
        modes.append((kap, v, c))
    return modes


def random_point(d, rng):
    nsl = d + 2
    X = [mp.mpf(rng.uniform(-1, 1)) for _ in range(nsl)]
    aux = []
    for mu in range(nsl):
        if mu == 0:   # Y away from 0
            aux += [mp.mpf(rng.uniform(0.6, 1.2)), mp.mpf(rng.uniform(-0.4, 0.4))]
        else:
            aux += [mp.mpf(rng.uniform(-0.3, 0.3)), mp.mpf(rng.uniform(-0.3, 0.3))]
    return X + aux


def tangency_and_time(coords, d, modes):
    nsl = d + 2
    ch = chi(coords, d, modes)
    om = [(coords[nsl + 2 * m], coords[nsl + 2 * m + 1]) for m in range(nsl)]
    Y, P = om[0], om[1:1 + d]
    W = [kdiv(Pj, Y) for Pj in P]
    WW = (sum(kmul(w, w)[0] for w in W), sum(kmul(w, w)[1] for w in W))
    ell = [(-WW[0] / 2, -WW[1] / 2)] + W + [(mp.mpf(1), mp.mpf(0))]
    tang = (mp.mpf(0),) * 4
    for mu in range(nsl):
        tang = qadd(tang, qmul(kq(ell[mu]), ch[mu]))
    tang_norm = max(abs(x) for x in tang)
    # time component in physical coordinates: T = (t - s0)/sqrt 2
    chiT = qscale(1 / mp.sqrt(2), qadd(ch[0], qscale(-1, ch[-1])))
    lhs = sum(x * x for x in chiT)
    # S_0, S_1: multipliers b_l m on each mode (b_0 + k b_1 = V_T^v)
    X = coords[:nsl]
    y = mp.sqrt(Y[0] ** 2 + Y[1] ** 2)
    S0 = mp.mpc(0)
    S1 = mp.mpc(0)
    for kap, v, c in modes:
        Wv = (sum(W[j][0] * v[j] for j in range(d)), sum(W[j][1] * v[j] for j in range(d)))
        b0 = (1 + Wv[0] / 2) / mp.sqrt(2)
        b1 = Wv[1] / (2 * mp.sqrt(2))
        dist2 = sum((W[j][0] - v[j]) ** 2 + W[j][1] ** 2 for j in range(d))
        m = (kap / y) ** (mp.mpf(d) / 2) * mp.exp(-kap * y * dist2 / 2)
        phase = kap * (X[-1] + sum(v[j] * X[1 + j] for j in range(d)) - X[0] * sum(x * x for x in v) / 2)
        u = c * mp.expj(phase)
        S0 += b0 * m * u
        S1 += b1 * m * u
    rhs = abs(S0) ** 2 + abs(S1) ** 2
    return tang_norm, abs(lhs - rhs) / max(rhs, mp.mpf("1e-30"))


def main():
    rng = random.Random(2026)
    out = {"conservation": [], "negative_controls": [], "tangency_and_time": []}
    ok = True
    for d in (1, 2):
        for trial in range(3):
            modes = random_modes(d, rng, 3)
            pt = random_point(d, rng)
            dv, sc = divergence(pt, d, modes)
            rel = float(dv / sc)
            out["conservation"].append({"n": d + 1, "trial": trial, "max_abs_div": float(dv),
                                        "tensor_scale": float(sc), "relative": rel})
            ok &= rel < 1e-30
            tn, terr = tangency_and_time(pt, d, modes)
            out["tangency_and_time"].append({"n": d + 1, "trial": trial, "dPhi_chi": float(tn),
                                             "time_component_rel_err": float(terr)})
            ok &= tn < mp.mpf("1e-40") and terr < mp.mpf("1e-40")
        modes = random_modes(d, rng, 3)
        pt = random_point(d, rng)
        for label, kw in (("no rotation iota_e -> iota_1", {"rotate": False}),
                          ("no volume factor (kappa/y)^{d/2}", {"volume": False})):
            dv, sc = divergence(pt, d, modes, **kw)
            rel = float(dv / sc)
            out["negative_controls"].append({"n": d + 1, "control": label, "relative": rel})
            ok &= rel > 1e-6
    out["all_ok"] = bool(ok)
    print(json.dumps(out, indent=1))
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    with open(__file__.rsplit("/", 2)[0] + "/results/wave_current_target.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
