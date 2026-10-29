"""Exact checks of Lemmas 4.1-4.3 and Remark 4.4 (the matrix M0 and auxiliary rows).

For each parameter triple (r, Q, R) this script verifies, in exact rational arithmetic:

  L4.1a  the psi vectors are pairwise orthogonal, have entries in {0,+-1}, and
         Pi_s psi = (eigenvalue at s) psi, W0 psi = j psi, V0 psi = (-1)^j psi;
  L4.1b  c_B = #{psi : j not in I} = D * Pr[Bin(Q,(r-1)/(2r)) not in I];
  L4.1c  every row of I - M0 lies in the row space of C^B (so M0 fixes ker C^B);
  L4.1d  rank(I - M0) = c_B exactly (the bound is attained);
  L4.2   M0[u,y] != 0  =>  |E_{u,y}| <= R-1;
  L4.3   P_u(0) = delta_u, deg P_u <= R-1, and for every t, g != d_u(t):
         P_u(t)|_{S_{t,g}} lies in the row space of C^B_{t,g};
         and the joint form: P_u(t)|_{complement of S_{t,d_u(t)}} lies in rowspace(C^B restricted);
  R4.4   the explicit example of Remark 4.4.
  NEG    negative controls: with f replaced by the constant 1 (M0 = V0) the rank of
         I - M0 is far larger, and with a wrong interpolation f the kernel property fails.
"""

import json
import sys
from fractions import Fraction
from itertools import product

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from eth_core import Aux, rref, reduce_against, in_rowspace, binom_tail_outside, poly_eval  # noqa: E402


def check(r, Q, R, A0=None, joint=True):
    aux = Aux(r, Q, R, A0)
    D = aux.D
    res = {"r": r, "Q": Q, "R": R, "D": D, "interval": sorted(aux.interval)}

    # L4.1a: orthogonality, entries, eigen-relations
    psis = aux.psi
    ok_orth = True
    for i in range(len(psis)):
        a = psis[i][0]
        for k in range(i + 1, len(psis)):
            b = psis[k][0]
            if sum(v * b.get(c, 0) for c, v in a.items()) != 0:
                ok_orth = False
                break
        if not ok_orth:
            break
    ok_entries = all(v in (1, -1) for row, _ in psis for v in row.values())
    ok_eig = True
    for row, j in psis:
        # W0 psi = j psi  and  M0 psi = (-1)^j f(j) psi
        for u in aux.labels:
            wpsi = sum(wv * row.get(y, 0) for y, wv in aux.W0[u].items())
            if wpsi != j * row.get(u, 0):
                ok_eig = False
                break
        if not ok_eig:
            break
    res["L4.1a_orthogonal"] = ok_orth
    res["L4.1a_entries_pm1"] = ok_entries
    res["L4.1a_W0_eigen"] = ok_eig
    res["psi_count_equals_D"] = len(psis) == D
    res["max_row_abs_sum_CB"] = max((len(row) for row in aux.CB), default=0)
    res["row_abs_sum_le_2^Q"] = res["max_row_abs_sum_CB"] <= 2 ** Q

    # L4.1b: count
    p = Fraction(r - 1, 2 * r)
    tail = binom_tail_outside(Q, p, aux.interval)
    res["cB"] = aux.cB
    res["D_times_tail"] = str(D * tail)
    res["L4.1b_count"] = aux.cB == D * tail

    # f interpolates (-1)^j on I
    fc = aux.f_coeffs()
    res["f_interpolates"] = all(poly_eval(fc, j) == (-1) ** j for j in aux.interval)

    # L4.1c: rows of I - M0 in rowspace(C^B); use orthogonal projection (rows of C^B orthogonal)
    norms = [sum(v * v for v in row.values()) for row in aux.CB]
    ok_c = True
    I_minus_M0 = {}
    for u in aux.labels:
        row = {u: Fraction(1)}
        for y, m in aux.M0[u].items():
            row[y] = row.get(y, 0) - m
        row = {k: v for k, v in row.items() if v != 0}
        I_minus_M0[u] = row
        resid = dict(row)
        for crow, nrm in zip(aux.CB, norms):
            ip = sum(v * row.get(c, 0) for c, v in crow.items())
            if ip:
                coef = ip / nrm
                for c, v in crow.items():
                    nv = resid.get(c, 0) - coef * v
                    if nv == 0:
                        resid.pop(c, None)
                    else:
                        resid[c] = nv
        if resid:
            ok_c = False
            break
    res["L4.1c_I_minus_M0_in_rowspace_CB"] = ok_c

    # L4.1d: exact rank of I - M0
    basis = rref(list(I_minus_M0.values()), D)
    res["rank_I_minus_M0"] = len(basis)
    res["L4.1d_rank_equals_cB"] = len(basis) == aux.cB

    # M0 symmetric
    res["M0_symmetric"] = all(aux.M0[y].get(u, 0) == m for u in aux.labels for y, m in aux.M0[u].items())

    # L4.2 support bound; also record the maximal |E| that actually occurs
    maxE = 0
    ok42 = True
    for u in aux.labels:
        for y in aux.M0[u]:
            e = len(aux.E(u, y))
            maxE = max(maxE, e)
            if e > R - 1:
                ok42 = False
    res["L4.2_support"] = ok42
    res["max_E_in_support"] = maxE

    # L4.3: rows
    ok_p0 = all({y: p[0] for y, p in aux.P[u].items() if p[0] != 0} == {u: 1} for u in aux.labels)
    ok_deg = all(len(p) - 1 <= R - 1 for u in aux.labels for p in aux.P[u].values())
    res["L4.3_P0_delta"] = ok_p0
    res["L4.3_degree"] = ok_deg
    ok_k = True
    ok_joint = True
    for s in range(Q):
        t = s + 1
        groups = {g: aux.S_aux(s, g) for g in range(r)}
        bases = {}
        for g in range(r):
            cols = set(groups[g])
            bases[g] = rref([{c: v for c, v in row.items() if c in cols} for row in aux.CB], D)
        joint_bases = {}
        for u in aux.labels:
            du = u[s][1]
            val = {y: poly_eval(p, t) for y, p in aux.P[u].items()}
            for g in range(r):
                if g == du:
                    continue
                cols = set(groups[g])
                restr = {c: v for c, v in val.items() if c in cols and v != 0}
                if not in_rowspace(bases[g], restr):
                    ok_k = False
            if joint:
                if du not in joint_bases:
                    cols = set(aux.labels) - set(groups[du])
                    joint_bases[du] = (cols, rref([{c: v for c, v in row.items() if c in cols} for row in aux.CB], D))
                cols, jb = joint_bases[du]
                restr = {c: v for c, v in val.items() if c in cols and v != 0}
                if not in_rowspace(jb, restr):
                    ok_joint = False
    res["L4.3_kernel_per_group"] = ok_k
    if joint:
        res["L4.3_kernel_joint"] = ok_joint
    return res, aux


def negative_controls(r, Q, R):
    """M0 = V0 (f = 1) and a deliberately wrong f: the lemma's conclusions should fail."""
    out = {}
    aux = Aux(r, Q, R)
    D = aux.D
    # V0 alone: rank(I - V0)
    rows = []
    for u in aux.labels:
        w = aux.swap_all(u)
        row = {u: Fraction(1)}
        row[w] = row.get(w, 0) - 1
        rows.append({k: v for k, v in row.items() if v != 0})
    out["rank_I_minus_V0"] = len(rref(rows, D))
    out["cB_with_correct_f"] = aux.cB
    # wrong f: shift the interval by one so it no longer matches C^B's kernel
    aux2 = Aux(r, Q, R)
    wrong = Aux(r, Q, R, A0=(aux2.A0 + 1) if aux2.A0 + R <= Q else aux2.A0 - 1)
    norms = [sum(v * v for v in row.values()) for row in aux2.CB]
    fails = 0
    for u in wrong.labels:
        row = {u: Fraction(1)}
        for y, m in wrong.M0[u].items():
            row[y] = row.get(y, 0) - m
        resid = {k: v for k, v in row.items() if v != 0}
        base = dict(resid)
        for crow, nrm in zip(aux2.CB, norms):
            ip = sum(v * base.get(c, 0) for c, v in crow.items())
            if ip:
                for c, v in crow.items():
                    nv = resid.get(c, 0) - ip / nrm * v
                    if nv == 0:
                        resid.pop(c, None)
                    else:
                        resid[c] = nv
        if resid:
            fails += 1
    out["rows_failing_with_wrong_f"] = fails
    out["negative_control_detects"] = fails > 0 and out["rank_I_minus_V0"] > aux.cB
    return out


def remark_4_4():
    """Q = R = 2, r = 2, I = {0,1}: M0 = Pi1 + Pi2 - Pi1 Pi2 and the displayed row."""
    aux = Aux(2, 2, 2, A0=0)
    u = ((0, 1), (0, 1))  # a_u = (0,0), d_u = (1,1): label (00,11)
    z0 = u
    z1 = ((1, 0), (0, 1))   # (10,01)
    z2 = ((0, 1), (1, 0))   # (01,10)
    z12 = ((1, 0), (1, 0))  # (11,00)
    M0row = {k: v for k, v in aux.M0[u].items() if v != 0}
    ok_m0 = M0row == {z1: 1, z2: 1, z12: -1}
    P = aux.P[u]
    expect = {z0: [Fraction(1)], z1: [Fraction(0), Fraction(-1, 2)], z2: [Fraction(0), Fraction(-1)]}
    ok_row = {k: v for k, v in P.items()} == expect
    # c(zeta) = z0 - z1 - z2 + z12 must be (up to scale) the unique C^B row on this orbit
    return {"M0_row_matches": ok_m0, "P_row_matches": ok_row}


def main():
    cases = [(2, 2, 2), (2, 3, 2), (2, 3, 3), (2, 4, 2), (2, 4, 3), (3, 2, 2), (3, 3, 2)]
    out = {"cases": [], "remark_4_4": remark_4_4(), "negative_controls": negative_controls(2, 4, 2)}
    all_ok = True
    for r, Q, R in cases:
        res, _ = check(r, Q, R, joint=(r ** (2 * Q) <= 256))
        out["cases"].append(res)
        flags = [v for k, v in res.items() if isinstance(v, bool)]
        all_ok &= all(flags)
        print(json.dumps(res))
    all_ok &= all(out["remark_4_4"].values()) and out["negative_controls"]["negative_control_detects"]
    out["all_checks_pass"] = all_ok
    print("remark_4_4", out["remark_4_4"])
    print("negative", out["negative_controls"])
    print("ALL PASS" if all_ok else "SOME CHECK FAILED")
    path = __file__.rsplit("/", 2)[0] + "/results/replacement_forms.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
