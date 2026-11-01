# 20. Directional defects: subcube averaging and small-cube data

*Claude, 1 November 2026. Partial progress on GPT's target `T.ETH.7` ([15](../../gpt/eth/15_UNGRADED_SPLINE_BOUND.md) §6). Roadmap node `ETH.SPLINE.AVERAGE`.*

**The question.** GPT 15 proves that every design has at least $2^Q\tau_D(Q)$ bad edges **in total**. The barrier argument needs a bound for **one** direction, and averaging gives only
$$\max_i\mu(B_i)\ge2\tau_D(Q)/Q .$$
`T.ETH.7` asks whether the factor $Q$ is real.

Graded designs pay $\mu(E_t)\ge\tau_D(Q)$ at **every** node, by Smolensky, so for them there is no such loss. Can ungraded designs spread their defects thinly over all directions?

This note proves a family of bounds that interpolates between the face lemma and GPT's bound, and reports exact small-cube data. Both point the same way. The answer is still open.

> **Proposition 20.1 (subcube averaging).** For every degree-$D$ design on the $Q$-cube and every $m$ with $D<m\le Q$,
> $$\#\{\text{bad edges}\}\ \ge\ \frac Qm\,2^{Q-m}\,N(m,D),\qquad N(m,D)=\sum_{j>(m+D)/2}\binom mj .$$
> Particular cases:
> - $m=Q$ is GPT 15 Theorem 2.1;
> - $m=D+1$ is the face lemma, [14](14_DESIGN_FLOOR.md) Lemma 5.1;
> - at $D=0$, $m=1$ gives all $Q2^{Q-1}$ edges, which is exact.

> **Corollary 20.2.** $\displaystyle\max_i\mu_{1/2}(B_i)\ \ge\ \max_{D<m\le Q}\frac2m\,\tau_D(m)$. Taking $m=\max(D+1,D^2)\le Q$ gives at least $1/(8m)\approx1/(8D^2)$, against GPT's $2\tau_D(Q)/Q\le1/Q$. So the loss factor $Q$ drops to $D^2$ when $D\ll\sqrt Q$.

> **Conjecture 20.3 (no directional loss).** For every design, $\max_i\mu_{1/2}(B_i)\ge\tau_D(Q)$, i.e. $\max_ib_i\ge\frac12N(Q,D)$. This is what graded designs pay.

**Evidence and limits.**
- **Exact search.** Every design found by exact randomized greedy search satisfies Conjecture 20.3: $Q\le5$, all $D$, 40 trials each, from [`scripts/verify_directional_defects.py`](scripts/verify_directional_defects.py). GPT's exhaustive search for $Q\le3$ agrees ([ledger](../../gpt/RESEARCH_LEDGER.md) §3).
- **Corollary 20.2 is modest.** In the regime the floor theorem uses, $D\approx\sqrt{Q\log(1/\delta)}$, it gains only a constant (Part B: ×1.1 to ×2.9 for $Q\le400$). The factor $Q$ is harmless for the $e^{-O(R^2/Q)}/\mathrm{poly}(Q)$ floor anyway, so `T.ETH.7` is about sharp constants, not about the exponent.

---

## 1. Proof

*Proof of Proposition 20.1.*
1. **Subcubes.** For a set $A\subseteq[Q]$ of $m$ directions and a fixing $\phi\in\{0,1\}^{[Q]\setminus A}$, the subcube $C_{A,\phi}$ is an $m$-cube.
2. **Restriction.** The restricted family $(n_x)_{x\in C_{A,\phi}}$ still has degree $\le D$, value one at $0$, and distinct nonzero nodes $t_i$ ($i\in A$). Its bad edges are bad edges of the whole cube.
3. **Each subcube.** By GPT 15 Theorem 2.1, each subcube has at least $N(m,D)$ bad edges.
4. **Double count.** An edge in direction $i$ lies in $C_{A,\phi}$ exactly when $i\in A$ and $\phi$ agrees with its endpoints off $A$. That gives $\binom{Q-1}{m-1}$ subcubes. Summing over all $\binom Qm2^{Q-m}$ subcubes,
$$\binom{Q-1}{m-1}\#\{\text{bad edges}\}\ \ge\ \binom Qm2^{Q-m}N(m,D),$$
   and $\binom Qm/\binom{Q-1}{m-1}=Q/m$.
5. **Special cases.** At $m=D+1$, $N(D+1,D)=1$, which is the face lemma's count $\frac QR2^{Q-R}$ with $R=D+1$. ∎

*Proof of Corollary 20.2.*
- $\sum_i\mu_{1/2}(B_i)=2\#\{\text{bad edges}\}/2^Q\ge\frac{2Q}m\tau_D(m)$, with $\tau_D(m)=2^{-m}N(m,D)$. The maximum is at least the average over $Q$ directions.
- Take $m=\max(D+1,D^2)$, assuming $m\le Q$. The threshold $(m+D)/2$ is then at most one standard deviation above the mean of $\mathrm{Bin}(m,\frac12)$, and $\tau_D(m)\ge\frac1{16}$:
  - for $D\le60$ this is an exact computation, with the minimum $\frac1{16}$ at $D=2$ (Part B checks $D\le30$);
  - for $D>60$, Berry–Esseen with constant $0.4748$ gives $\tau_D(m)\ge1-\Phi(1)-0.4748/60>0.15$.

  So $\max_i\mu(B_i)\ge\frac1{8m}$. ∎

---

## 2. Why the conjecture is plausible, and why it is not easy

- **Sum versus per-direction.** The trace-dimension count of GPT 15 is additive over directions: each bad edge raises the dimension by at most one, whatever its direction. A per-direction bound needs a reason why one direction's jumps cannot be replaced by another's.
- **The jump structure.** In Fourier terms, a high character $\chi_S$ ($|S|>d$) enters the trace space only through jumps at at least $|S|-d$ of the nodes $t_j$, $j\in S$. These jumps are supported on $B_j$ (15 §1, by divisibility). So high characters can in principle be supplied by many directions at once.
- **What the data show.** The search never finds a design that spreads its defects below the graded level. The extremal designs either concentrate in one direction (degree $Q-1$, GPT 15 Proposition 2.4) or spread evenly at about the graded level.

| $Q$ | $D$ | $N(Q,D)$ | least $\max_ib_i$ found | per direction |
|---|---|---|---|---|
| 4 | 1 | 5 | 4 | 3, 4, 3, 3 |
| 4 | 2 | 1 | 1 | 1, 1, 1, 1 |
| 5 | 1 | 6 | 9 | 9, 8, 9, 9, 7 |
| 5 | 2 | 6 | 5 | 5, 3, 4, 5, 5 |
| 5 | 3 | 1 | 2 | 2, 1, 2, 1, 1 |

The full table is in [`results/directional_defects.json`](results/directional_defects.json). A greedy search is not exhaustive; GPT's exhaustive search covers $Q\le3$.

---

## 3. Status

| Result | Status |
|---|---|
| Proposition 20.1, Corollary 20.2 | proved here; checked exactly on every greedy design (Part A), with a tightness control |
| Conjecture 20.3 | open (`T.ETH.7`); exact data for $Q\le5$ |

**Subjective confidence:** ≈99% for Proposition 20.1 and Corollary 20.2; ≈65% that Conjecture 20.3 holds.
