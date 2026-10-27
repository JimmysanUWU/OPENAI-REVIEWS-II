# 3. Methods, leaps of insight, and how radical they are

## 3.1 Inventory of methods

Each row gives the method, where it is used, its known antecedent, and what is new.

| # | Method | Where | Known antecedent | What is new |
|---|---|---|---|---|
| M1 | **Deformed substitutions.** Replace each clause input $x_v$ by $\langle P_{e,v}(T),x\rangle$, with $P(0)=e_v$, so the count becomes the constant term of a polynomial $H(T)$. | §2.3, §3.1 | Interpolation in an auxiliary parameter: Valiant; Björklund 2016; Björklund–Kaski–Williams 2019 | The deformation targets a structural property (locality) at many nodes at once, at degree far below the number of nodes. |
| M2 | **Separation by shared integer measurements.** Fixing a few integer linear summaries splits the sum into independent group sums. | §2.4, Lemma 3.2 | Split-and-list (Williams 2005; Koivisto 2006); meet-in-the-middle on linear sums | The (r,τ,D)-design abstraction. Measurements are enumerated as values, with $\ell_1$-bounded rows controlling their range. |
| M3 | **Jets from Boolean substitutions.** Treat the shift $P(t+z)-P(t)$ as acting on an independent formal vector $X$. Then $[z^j]$ has $X$-degree $\le j$, so Möbius inversion over $\lvert E\rvert\le q$ recovers the diagonal sum. | Lemma 3.3 | Hasse-derivative/multiplicity evaluation: BGKM23, KSY14, Woodruff–Yekhanin | Combining jets with the design: the fixed substitution $X=\mathbf 1_U$ keeps separation intact. |
| M4 | **Swap interpolation.** $M_0=V_0f(W_0)$, a polynomial in commuting coordinate swaps, supported within Hamming distance $R-1$ of the complete swap and equal to $I$ off a rank-$\varepsilon D$ subspace. | §4.2, Lemmas 4.1–4.3 | Alman–Williams non-rigidity of Walsh–Hadamard (2017), which approximates parity on the binomial bulk; Dvir–Liu 2020 | Transplanting the parity-on-the-bulk phenomenon into the swap group algebra, so that it produces locality-preserving deformations. |
| M5 | **Universal two-schedule array.** Labels $u=(a_v,d_f)\in([r]^Q)^2$, so one construction serves every formula; incidences are placed into $h$ copies. | §4.1, §4.4–4.5 | Random schedule assignment with derandomization by conditional expectations | The device of putting every pair (variable group, clause group) on a common auxiliary array. |
| M6 | **Cross-equation 2-adic truncation.** $p_\ell=\sum e_jh^j$, with $h=Y(Y-1)$ even, truncated in total $h$-power across all equations; this keeps lifted $\mathbb F_2$ systems at linear degree. | Lemma 5.2 | Beigel–Tarui modulus amplification; Chan–Williams 2016 | The global truncation with a formal marker $V$. I have not seen it before; the degree saving ($m(2\ell-1)\to m+2\ell$) is exactly what the method needs. |
| M7 | **Tanner-code initial checks plus greedy rows chosen by exact counts.** | §6 | Tanner 1981; Sipser–Spielman 1996; $\varepsilon$-biased sets; conditional expectations | Standard; the counting oracle derandomizes the parity tests. |
| M8 | **Disjoint sparsification.** | App. A | IPZ01; CIP06; DHM+14 | A self-contained re-proof with explicit constants. |
| M9 | **Cauchy deformation with 2-adic tapering.** $W(T)=\mathcal C(T)/(\mathcal C'(0)(T-(y-x)))$, Hermite multiplicities tapered by $v_2(t)$, and "modes" of 2-adic weight replacing subset enumeration. | App. D | Lagrange/Cauchy kernels; 2-adic precision analysis in exact algebraic algorithms | Using 2-adic valuation, rather than low degree, to limit the interaction terms. This removes one $\log\log n$ factor. |

## 3.2 The leaps of insight

**Leap 1 (conceptual reformulation).** Counting #SAT faster is reformulated as finding a low-degree family of linear substitutions that is "local up to few measurements" at many nodes.

The three costs then appear transparently:
- the Taylor depth $q\sim nR/Q$;
- the measurement values $(2D+1)^{2\tau n}$;
- the group size $2^{2n/r}$.

Before this paper there was no known reason to expect all three could be made $2^{o(n)}$ at once. The paper's overview (§2.3–§2.7) explains this reformulation well.

**Leap 2 (the key lemma).** The algorithm needs degree $R\ll Q$, yet the condition at each node looks like it forces degree $\ge Q$.

A naive interpolation of the needed form at each of the $Q$ nodes has degree $Q$. The paper observes three things:
1. A clause's input can be split as (part determined by measurements) + (terms, each of which uses an out-of-group variable at only a few nodes).
2. The second part can be killed at those few nodes with the multipliers $\prod_{s\in E}(1-T/s)$.
3. The required matrix is a polynomial of degree $<R$ in commuting swaps that fixes most eigenvectors.

The third observation is where the analysis lives. It reduces to: *a degree-$(R-1)$ polynomial can equal $(-1)^j$ on the $R$ central values of $\mathrm{Bin}(Q,p)$*, which carry mass $1-e^{-\Omega(R^2/Q)}$. Hence $R\approx\sqrt{Q\log(1/\varepsilon)}\ll Q$.

This is the same phenomenon that makes Walsh–Hadamard matrices non-rigid (Alman–Williams). It has been moved to a new representation: the swap action on pairs of schedules.

**Leap 3 (universality).** Index auxiliary coordinates by *pairs of group schedules*. Then a single matrix $M_0$, independent of the formula, serves every incidence. Concentration of label multiplicities ($h\approx M/D$) converts the auxiliary rank $\varepsilon D$ into $\approx\varepsilon M$ actual measurements.

**Leap 4 (precision as degree).** Two places use 2-adic smallness as a substitute for low degree:
- the cross-equation truncation (M6);
- Appendix D's tapered recovery (M9).

## 3.3 How radical is it?

**In consequence: very.** ETH has been a working assumption of fine-grained complexity since 2001, and thousands of conditional lower bounds rest on it. A deterministic, polynomial-space $2^{o(n)}$ algorithm for #3-SAT overturns that foundation.

**In technique: a new architecture built from classical parts.** Every individual tool predates the paper:
- interpolation and polynomial-method thinking;
- split-and-list;
- the parity approximation behind matrix non-rigidity;
- sparsification and expander codes;
- derandomization by conditional expectations.

What is new is the architecture of M1 + M2 + M4 + M5. Particularly new is the realization that matrix non-rigidity, a phenomenon studied for circuit lower bounds, can be run *in reverse* as an algorithmic decoupling device for SAT.

I know of no prior work suggesting this route. It is not a tighter analysis of PPSZ or of local search, and it is not a variation on known $2^{n-n/\mathrm{polylog}}$ algebraic counting algorithms.

**In strength: modest.** The saving is from $2^n$ to $2^{n/(\log n)^{1/2-o(1)}}$. That is a "$\sqrt{\log}$-shaving" in the exponent, reminiscent of $n^3/2^{\Omega(\sqrt{\log n})}$ APSP. Both methods (§3–4 and App. D) stop at $\sqrt{\log n}$ for the reasons isolated in [04 §5](04_MIRACLES_EXPLAINED.md).

**Comparison with human prior expectations.** The community regarded ETH as plausible largely because no algorithmic paradigm came close. The paper does not find a weakness in SAT's combinatorics. It finds an *algebraic* slack: low-degree deformations exploit the fact that a clause polynomial cares only about a few linear forms at a time.

## 3.4 Is it "independently creative"?

Judged as mathematics, Leap 2 is the kind of idea one expects from a strong specialist in algebraic algorithms or matrix rigidity. It is non-obvious, simple once seen, and explained by a one-line spectral fact.

Leaps 1 and 3 are the decisive creative steps. They connect rigidity-style polynomial approximation to SAT algorithms through a deformation-and-locality framework, which, as far as I know, did not exist.

Appendix D shows further independent technique. It reaches the same exponent shape by an entirely different 2-adic mechanism.

I do not detect amalgamation. There is no lemma that looks transplanted from a specific source. The paper cites its antecedents accurately for the tools it uses, apart from the related-work omission noted in [02](02_VERIFICATION_LEDGER.md).
