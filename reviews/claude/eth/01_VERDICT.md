# 1. Verdict: does "The Exponential Time Hypothesis Is False" answer the question?

Claude review, 27 October 2026. Source: `ETHDISPROVE.pdf`, 67 pp., OpenAI, dated 17 Oct 2026 and revised 25 Oct 2026 (SHA-256 in [`../../../source/eth/SOURCES.md`](../../../source/eth/SOURCES.md)). Lean statement: [`../../../source/eth/SubexponentialSAT.statement.lean`](../../../source/eth/SubexponentialSAT.statement.lean).

## Short answer

**Yes.** For every fixed $k$, the paper gives a deterministic algorithm that counts the satisfying assignments of a $k$-CNF in $2^{o(n)}\,\mathrm{poly}(n+m)$ time and polynomial space. If it is correct, ETH is false.

The part that refutes ETH is short and elementary. It also proves the result on its own, without the rest of the paper. I re-derived every step of it and found no gap. I also ran the complete algorithm in exact arithmetic on 17 random CNFs. It returned the exact count every time, and corrupting the design made it fail.

| Claim | Where proved | My status | Confidence |
|---|---|---|---|
| **ETH false** (and SETH, #ETH, NSETH, Gap-ETH): #$k$-SAT in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}$ | §2.8 + §3 + §4 + App. A + App. B.4 | Every lemma re-derived; core lemmas and full pipeline checked in exact arithmetic. **Also formally verified via Prop. 6.4 (below).** | **≈99%** |
| Prop. 6.4: #$k$-SAT in $2^{O_k(\rho(n))}$, $\rho(n)=n(\log\log n)^2/\sqrt{\log n}$ | §5, §6, App. A, App. B | Re-derived; Lemmas 5.1/5.2 checked exactly; **OpenAI's Lean proof built locally, only standard axioms, statement byte-identical** | **≈99%** |
| Thm 1.1 headline: $2^{O_k(n\log\log n/\sqrt{\log n})}$, and $2^{O(n/\sqrt{\log n})}$ at fixed density | App. D (tensor contraction) | **Re-derived completely (28 Oct, [10](10_GAPLESS_TENSOR_CONTRACTION.md))**, including the recursion's bookkeeping. The child identity (D.16), the tile expansion, one recursion level at the honest precision, and both encodings are checked exactly. Not formalized. | ≈95% |
| Cor. 1.2 (Max-$k$-SAT), Cor. 1.3 (colourings, matchings), Cor. C.2 (local maps) | standard reductions | Proved in [11](11_GAPLESS_COROLLARIES.md) | ≈98% |
| Cor. 1.4 / C.1 (E^NP circuit lower bounds) | Williams' witness-circuit argument, cited | Parameter bookkeeping checked in [11](11_GAPLESS_COROLLARIES.md). Cor. 1.4's $n\sqrt{\log n}$ needs App. D; the formally verified Prop. 6.4 alone gives $n\sqrt{\log n}/(\log\log n)^2$ | ≈93% (rests on cited theorems) |
| Lean statement `exact_count_and_sat` | — | Faithful to Prop. 6.4. **OpenAI's formal proof builds; only standard axioms; elaborated statement identical** ([06 §6.A](06_LEAN_STATEMENT_AUDIT.md)) | formally verified |

These are subjective confidences after verification, not measured error rates.

## The minimal ETH chain

1. **Disjoint sparsification** (App. A). This partitions the satisfying assignments over at most $\mathrm{poly}\cdot 2^{\sigma n}$ leaf formulas. In each leaf, every variable occurs at most $B=O_k((\sigma^{-1}\log(2/\sigma))^{k-1})$ times.
2. **Standard reductions** (§2.8). Binary OR gates and Tovey's implication cycles turn each leaf into a 3-CNF on $N=O_k(Bn)$ variables. Every variable occurs at most three times, and the count is preserved.
3. **The new counting algorithm** (§3–§4). It counts such a formula in $2^{O(N(\log\log N)^2/\sqrt{\log N})}$ time.
4. **Parameters.** Take $\sigma\asymp(\log n)^{-1/(2k)}$. Both $\sigma n$ and $BN(\log\log N)^2/\sqrt{\log N}$ are then $o(n)$.

The new mathematics is step 3. It needs only three external facts: Hoeffding's inequality, Hermite interpolation over $\mathbb Q$, and Möbius inversion on subsets.

## How step 3 works, in one paragraph

Write $\#F=\sum_x\prod_e\Phi_e(x)$, where $\Phi_e$ is the clause polynomial. Each clause input $x_v$ is replaced by a linear form $\langle P_{e,v}(T),x\rangle$. Each form has degree $R$ in a new parameter $T$, and equals $x_v$ at $T=0$.

At each node $T=t\in\{1,\dots,Q\}$, every form becomes "local modulo measurements". It reads variables of one group directly and other groups only through a few shared integer linear measurements. Fixing the measurement values therefore splits the $2^n$-term sum into independent sums over groups of size $2n/r$.

The deformed sum $H(T)$ has degree $\le 3nR$. Its Taylor jets of order $q\approx 3nR/Q$ at the $Q$ nodes determine $H$, hence $H(0)=\#F$. Each jet is extracted from Boolean substitutions on at most $\binom{n}{\le q}2^q=2^{O(n\frac RQ\log\frac QR)}$ subsets.

The construction needs $R\ll Q$: low degree, yet local at every node. This is the paper's key idea. Take a symmetric polynomial $M_0$ in commuting coordinate swaps. It is supported near the complete swap, which is perfectly local. It equals the identity except on a rank-$\varepsilon D$ subspace. This holds because a degree-$(R-1)$ polynomial matches $(-1)^j$ on the $R$ central values of a binomial distribution, which carry all but $e^{-\Omega(R^2/Q)}$ of its mass.

## What the result does and does not mean

- **It refutes ETH as stated**: there is no $\delta>0$ for which 3-SAT needs $2^{\delta n}$ time. It also refutes SETH, #ETH, NSETH and Gap-ETH, both randomized and deterministic, because the algorithm is deterministic and counts exactly.
- **The speed-up is barely subexponential.** The running time is $2^{n/(\log n)^{1/2-o(1)}}$. Hypotheses such as "3-SAT needs $2^{n/\log n}$", the Orthogonal Vectors Conjecture, and anything about $2^{n^{\Omega(1)}}$ are untouched. P vs NP is untouched.
- **It is galactic.** With the paper's own constants, all hypotheses of Lemma 4.7 first hold at $\log_2 n\approx1.1\cdot10^{10}$. The bound first beats $2^n$ at $\log_2 n\approx7.4\cdot10^{10}$ ([`results/parameters.json`](results/parameters.json)). Better constants are surely possible, but the method's $(\log\log n)^2/\sqrt{\log n}$ saving is intrinsic to it (see [04](04_MIRACLES_EXPLAINED.md)). This does not affect the refutation, which is asymptotic.
- Lower bounds conditioned on ETH or SETH lose that basis. Those that also follow from OVC or from weaker hypotheses survive, as the paper notes in §1.1.3.

## What remains open

- ~~**Appendix D.**~~ **Closed (28 Oct).** This separate mechanism (Cauchy deformation, 2-adic tapered Hermite recovery and mode expansion) supplies Theorem 1.1's sharper exponent and Corollary 1.4. Its recursion is now re-derived with explicit constants in [10](10_GAPLESS_TENSOR_CONTRACTION.md), and checked end to end on toy networks. I found no error. One stated inequality ($u\ge\lambda^{4/5}$) is stronger than the proof needs.
- ~~**The Lean proof.**~~ **Closed (27 Oct).** OpenAI's formal proof of exactly this statement builds locally (3,091 modules). It uses only `propext`, `Classical.choice` and `Quot.sound`. Its elaborated statement is byte-identical to the supplied file's. ETH's falsity is therefore formally verified, modulo Lean/Mathlib soundness. See [06 §6.A](06_LEAN_STATEMENT_AUDIT.md).
