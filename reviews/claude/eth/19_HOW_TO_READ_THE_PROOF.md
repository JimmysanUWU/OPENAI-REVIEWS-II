# 19. How to read the proof

*Claude, 1 November 2026. A reading guide to* The Exponential Time Hypothesis Is False *(OpenAI, revised 25 Oct 2026) and to this review's versions of its proof. It is the companion of [kakeya/15](../kakeya/15_HOW_TO_READ_THE_PROOF.md), with the same four levels.*

Each level stands on its own, so stop wherever you like. Every step names where it is proved, in the paper and in this repository, and which script checks it. A worked example runs the actual algorithm on a 4-variable formula and prints every intermediate number: [`scripts/verify_worked_example.py`](scripts/verify_worked_example.py).

| Level | For readers who know | What you get | Time |
|---|---|---|---|
| [0](#level-0-the-problem-and-the-idea) | nothing in particular | what the problem is, why the result is surprising, the idea of the proof | 15 minutes |
| [1](#level-1-the-proof-in-seven-moves) | high-school algebra, a little linear algebra | the proof as seven moves, one formula each, and a worked example | an hour |
| [2](#level-2-the-proof-on-one-page) | graduate discrete mathematics | the proof on one page; where every number comes from; the standard names of the objects | an afternoon |
| [3](#level-3-reading-routes) | the field, or Lean | reading routes through the paper, the gapless files and the formal proof; notation traps; checkpoints | days |

Then come [questions readers ask](#questions-readers-ask) and [what to reuse](#what-to-reuse).

---

## Level 0: the problem and the idea

**A puzzle with on/off switches.** A formula in *conjunctive normal form* is a list of requirements on $n$ switches, each on or off. Every requirement is small, for example "switch 1 is on, or switch 2 is off". A requirement is *satisfied* by a setting of the switches if at least one of its conditions holds. The question SAT asks is whether some setting satisfies every requirement at once. With at most $k$ conditions per requirement it is called $k$-SAT. The counting version, #SAT, asks *how many* settings do.

**The obvious method and the hypothesis.** Try all $2^n$ settings. Clever methods since the 1990s beat this, but only to $c^n$ for some $c<2$: still exponential.

In 1999 Impagliazzo and Paturi conjectured that this is unavoidable: 3-SAT needs time $2^{\delta n}$ for some fixed $\delta>0$. This is the **Exponential Time Hypothesis (ETH)**. Hundreds of "this problem needs exponential time" results assume it.

**The paper refutes ETH.** For each fixed $k$ it gives a deterministic algorithm that *counts* the satisfying settings in time
$$2^{O(n\log\log n/\sqrt{\log n})}.$$
The exponent is a vanishing fraction of $n$, so the time is "subexponential". Several things do **not** follow:
- **Not P = NP.** The bound is still far from polynomial.
- **Not practical.** It beats $2^n$ only for astronomically large $n$ (Level 2's ledger).

The stronger hypotheses SETH, #ETH and Gap-ETH fail with ETH.

**The idea: evaluate a polynomial where it is easy, and read off its value where it is hard.** A polynomial of low degree is pinned down by a little information at a few points. The count we want is the value at $T=0$ of a polynomial $H(T)$ built from the formula. At $T=0$, $H$ is the hard count. At special points $T=1,\dots,Q$, the *nodes*, the formula falls apart into a few independent pieces once a few numbers (*measurements*) are fixed. Each piece is small enough to brute-force. So:
1. **Simplify the formula** so that every switch appears in at most three requirements. This is standard, and keeps the count.
2. **Deform it.** Let each requirement read a mixture of switches that depends on a parameter $T$. At $T=0$ the mixture is the plain switch, so $H(0)$ is the count.
3. **Mix cleverly.** Give every (requirement, switch) pair an address, and let the mixtures be *swaps* of addresses. At each node the mixture only connects addresses inside one of $r$ groups, up to a small "defect".
4. **Pay for the defect with measurements.** The defect is small because of a fact about parity: a low-degree polynomial can agree with $\pm1$ parity on most of a binomial distribution. Fixing the defect's few linear readouts makes the groups independent.
5. **Collect derivatives at the nodes and interpolate back to $T=0$.** Each node needs only a few derivatives, because the degree of $H$ is spread over $Q$ nodes.

The time is a product of three costs: the size of the largest group, the number of measurement values, and the number of derivative terms. Balancing them gives $n/\sqrt{\log n}$ up to $\log\log$ factors.

**Where is the cleverness?** In step 3: commuting swaps, chosen so that the defect is a binomial tail. The rest is polynomial interpolation, inclusion–exclusion and Hoeffding's inequality. It is all undergraduate mathematics, assembled in a way nobody had tried.

**Is $\sqrt{\log n}$ an artefact?** No, for this design. This repository proves that no design inside the same algebra (graded or not) beats the paper's exponent:
- [14](14_DESIGN_FLOOR.md);
- GPT's [15](../../gpt/eth/15_UNGRADED_SPLINE_BOUND.md)–[16](../../gpt/eth/16_MEASUREMENT_TRANSFER_AND_FLOOR.md);
- [18](18_AGGREGATE_TRANSFER.md).

The square root is the agreement degree of parity, the same fact Smolensky used in 1987 for circuit lower bounds. Whether a different deformation could beat it is open (`T.ETH.5`).

---

## Level 1: the proof in seven moves

Notation: $\log=\log_2$. $F$ is a CNF with clauses $C_f$ and variables $x_1,\dots,x_n$. $\Phi_f$ is the clause polynomial: it equals $1$ exactly on assignments satisfying $C_f$.

```
 k-CNF ──(1)──▶ ≤3-occurrence 3-CNFs ──(2,3)──▶ H(T)=Σ_x Π_f Φ_f(P(T)x), H(0)=#F, deg H ≤ dR'
                                                          │
         (4) defect = binomial tail ──▶ measurements make each node split into r groups
                                                          │
         (5) jets of H at t=1..Q by group sums ──(6)──▶ H(0) by Chinese remainder ──▶ #F
                                                          │
         (7) time = (Taylor terms)·(measurement values)·2^(group size)  ⟹  2^{O(n(log log n)²/√log n)}
```

**Move 1. Make every variable rare.**
- *Disjoint sparsification* splits $F$ into at most $2^{\sigma n}\mathrm{poly}(n)$ formulas whose satisfying sets partition $\mathrm{Sat}(F)$, each with at most $B$ occurrences per variable.
- Binary OR gates and Tovey's implication cycles then give 3-CNFs where every variable occurs 1 to 3 times. The count is preserved exactly:
$$\#F=\sum_{\text{leaves }L}\#F'_L,\qquad F'_L:\ \text{width}\le3,\ \text{each variable in }1..3\text{ clauses}.$$

*Where:* paper App. A, §2.8; [08](08_GAPLESS_PROOF.md) §2–§3. *Check:* `verify_gapless_chain.py` Parts S, R.

**Move 2. Deform the readout.** Give each incidence (clause $f$, variable $v$) a row $P_{f,v}(T)$ with $P_{f,v}(0)=e_v$ and degree $\le R'$ in $T$. Then
$$H(T)=\sum_{x\in\{0,1\}^n}\prod_f\Phi_f\big((\langle P_{f,v}(T),x\rangle)_{v\in C_f}\big),\qquad H(0)=\#F,\qquad\deg H\le dR',$$
where $d=\sum_f|C_f|$.
*Idea:* only $H(0)$ matters, and a polynomial can be recovered from information elsewhere.
*Where:* paper §3; 08 §4 (Theorem 4.6, step 1).

**Move 3. Swap interpolation.**
- **Labels.** Each incidence gets a label $u=(a_v(s),d_f(s))_{s\le Q}\in([r]\times[r])^Q$ from a variable schedule $a_v$ and a clause schedule $d_f$. Swapping the pair at position $s$ is an involution $\Pi_s$ on labels.
- **The design.** Put $\varphi_S(T)=\prod_{s\in S}(1-T/s)$ and
$$P(T)=I-M_0+\sum_{k<R}\gamma_k\sum_{|S|=k}\varphi_S(T)\,\Pi_{S^c},\qquad M_0=\sum_{k<R}\gamma_k\sum_{|S|=k}\Pi_{S^c}.$$
- **Locality.** At a node $t$ every surviving term has $t\notin S$. It swaps position $t$, so it only connects labels whose group at $t$ agrees. So $P(t)$ is *local*: a row of group $g$ reads only variables of group $g$, except through the constant defect $I-M_0$.

*Where:* paper §4; 08 §5 (Theorem 5.4). *Check:* `verify_gapless_chain.py` Part K, `verify_distilled_swap_formula.py`.

**Move 4. The defect is a binomial tail.**
- The swaps commute, so they share an eigenbasis. An eigenvector $\psi$ with $j$ antisymmetric local factors gets the scalar $(-1)^jf(j)$ from $M_0$, where $f=\sum\gamma_kK_k$ (Krawtchouk polynomials).
- Choose $f$ to interpolate $(-1)^j$ on the $R$ levels around the mean $Qp_r$, $p_r=\frac{r-1}{2r}$. Then $I-M_0$ vanishes on those levels.
- Its rank is the binomial tail outside the window, by Hoeffding:
$$\frac{\operatorname{rank}(I-M_0)}{|\mathcal B|}=\Pr\big[W\notin\text{window}\big]\le2e^{-R^2/(8Q)},\qquad W\sim\mathrm{Bin}(Q,p_r).$$
- The measurements are the rows $\psi^{\mathsf T}$ of the measured eigenvectors, transferred to the variables.

*Idea:* parity agrees with a degree-$R$ polynomial on the bulk of a binomial distribution. That is the whole "miracle" ([04](04_MIRACLES_EXPLAINED.md)).
*Where:* paper §4; 08 Lemmas 5.1–5.3, Corollary 5.6, Lemma 6.1.

**Move 5. Jets by group sums.** Fix a node $t$.
- **Expand in the $T$-shift.** Write $T=t+z$ and expand to order $z^q$, $q=\lfloor dR'/Q\rfloor$. Boolean substitution (Möbius inversion on subsets) reduces the jet to $\sum_{|E|\le q}\sum_{U\subseteq E}\pm\sum_{x_E=1}p_x(\mathbf 1_U)$.
- **Split.** Each such sum splits over the $r$ groups once the measurement values $y$ are fixed:
$$\sum_{x}\prod_gA_g(x)=\sum_y\prod_g\ \sum_{x_g:\,C_gx_g=y_g}A_g\big(x_g,\ \text{fixed representatives elsewhere}\big).$$

*Where:* paper Lemmas 3.2–3.3; 08 Lemmas 4.4–4.5. *Check:* Part K, where separated jets equal direct jets.

**Move 6. Interpolate back.** $\deg H<QJ$ with $J=q+1$. The nodes $1,\dots,Q$ have unit differences modulo a prime $p>Q$, so the Chinese remainder theorem recovers $H(0)\bmod p^\mu$ from the $Q$ jets. Take $\mu=n+1$, so that $p^\mu>2^n\ge\#F$.
*Where:* paper §3, App. B; 08 Lemma 4.3.

**Move 7. Balance the three costs.** The time is
$$Q\cdot C_q(n)\cdot(2\beta+1)^{c_*}\cdot2r\,2^{s}\cdot\mathrm{poly},$$
with:
- $C_q(n)$ the Taylor terms;
- $(2\beta+1)^{c_*}$ the measurement values;
- $2^s$ the brute force inside the largest group.

Per variable, these are $x\log(e^2/x)$, $\approx24e^{-R^2/(8Q)}\cdot Q$ and $2/r$, with $x=3(R-1)/Q$. The choice
$$r=\lfloor\sqrt\lambda\rfloor,\qquad Q=\lfloor\lambda/(4\nu)\rfloor,\qquad R=2\lceil\sqrt\lambda\rceil$$
($\lambda=\lceil\log n\rceil$, $\nu=\lceil\log\lambda\rceil$) makes the total $\le13\nu^2/\sqrt\lambda$. The time is $2^{13n\nu^2/\sqrt\lambda}=2^{O(n(\log\log n)^2/\sqrt{\log n})}$ for the 3-CNFs of Move 1.

With the sparsification losses this refutes ETH for every $k$.
*Where:* paper §4, Lemma 4.8, (2.19); 08 §7 (Theorem 7.1, Lemma 7.2). *Check:* Part P; `verify_parameters.py`.

**Down the tower.**
- **The formal bound** $2^{O_k(n(\log\log n)^2/\sqrt{\log n})}$ for every $k$, without Tovey's loss: count common zeros of sparse $\mathbb F_2$ equations, lifted to $\mathbb Z/2^\ell$ by a closed-form parity polynomial (paper §5–6; [09](09_GAPLESS_SHARP_EXPONENT.md)). This is the statement of OpenAI's Lean proof ([06](06_LEAN_STATEMENT_AUDIT.md)).
- **The headline** $2^{O(n\log\log n/\sqrt{\log n})}$: contract a Boolean tensor network recursively, with a 2-adic Cauchy-kernel deformation and tapered Hermite recovery (App. D; [10](10_GAPLESS_TENSOR_CONTRACTION.md)).
- **Corollaries:** Max-$k$-SAT, colourings, matchings and circuit lower bounds ([11](11_GAPLESS_COROLLARIES.md)). Also counting zeros over any finite field, new here ([12](12_FP_COUNTING.md)).

### The worked example

[`scripts/verify_worked_example.py`](scripts/verify_worked_example.py) runs Moves 2–6 on
$$F=(x_1\vee\neg x_2)\wedge(x_2\vee x_3\vee\neg x_4)\wedge(\neg x_1\vee x_4)\wedge(x_3\vee x_4),$$
with $(r,Q,R)=(2,2,2)$. It uses the implementation checked in `verify_gapless_chain.py`.

| Step | What the script prints |
|---|---|
| Brute force | $\#F=5$ |
| Design (Moves 3–4) | 16 labels; levels $j=0,1,2$ hold $9,6,1$ eigenvectors ($=16\Pr[\mathrm{Bin}(2,\frac14)=j]$); window $\{0,1\}$; $\gamma=(-1,1)$; **1 measured eigenvector** (level 2) |
| Schedules | $a_v$ and $d_f$ from Lemma 6.2's derandomization; $M=9$ incidences, $h=1$ copy, potential $\Psi=0.073<2$ |
| Polynomial (Move 2) | $R'=1$, $d=9$, so $\deg H\le9$; $q=4$, $J=5$, $QJ=10>9$ |
| Arithmetic (Move 6) | modulo $3^5=243>2^4$ |
| Node $t=1$ (Move 5) | groups $\{x_1,x_3\},\{x_2,x_4\}$; one measurement row each; 2 and 3 measurement values; 81 correction terms; jet $(100,76,211,88,68)$, equal to the direct $2^n$-term jet |
| Node $t=2$ | groups $\{x_1,x_4\},\{x_2,x_3\}$; jet $(88,186,180,15,24)$ |
| Recovery | $H(0)\equiv5\pmod{243}$, which is $\#F$ |
| Controls | dropping one group's measurement row changes the jet; recovering from $J-1=4$ jets gives 182 |

At this size nothing is faster than brute force. The point is to see every object of the proof as a concrete number.

---

## Level 2: the proof on one page

**Theorem (08 Theorem 7.1).** Let $F$ be a CNF on $n$ variables with width $\le3$, every variable in 1 to 3 clauses, and $\lambda=\lceil\log n\rceil\ge2^{17}+1$. Then $\#F$ is computed in $\mathrm{poly}(n)\,2^{13n\nu^2/\sqrt\lambda}$ time and polynomial space.

*Proof on one page.*
1. **Design.** Take $(r,Q,R)$ as in Move 7, $p$ the least prime above $Q$, and $\mu=n+1$.
   - Theorem 5.4 gives $P(T)$ with $P(0)=I$ and degree $R-1$, local at $t=1..Q$ modulo the measured eigenvectors.
   - Those number $|\mathcal B|\Pr[W\notin\mathcal I]\le2|\mathcal B|e^{-R^2/(8Q)}$ (Corollary 5.6).
2. **Transfer.** Lemma 6.2 computes schedules with at most $h\le2M/|\mathcal B|$ copies of each label and groups of size $\le2n/r$. Lemma 6.1 turns $P$ into rows $P_{f,v}(T)$ and measurement rows of norm $\le2^Q$, at most $2hc$ per node.
3. **Count.** Theorem 4.6 computes $\#F\bmod p^\mu=\#F$. The cost per variable is
$$\underbrace{x\log(e^2/x)}_{\text{Taylor}}+\underbrace{24e^{-2\nu}(Q+2)}_{\text{measurements}}+\underbrace{2/r}_{\text{groups}}+O\Big(\frac{\log n}n\Big)\ \le\ \frac{13\nu^2}{\sqrt\lambda},$$
   by (P5)–(P8). ∎

**Every number.**

| Number | Value | Why this value | Where |
|---|---|---|---|
| $\sigma$ | $\asymp(\log n)^{-1/(2k)}$ | sparsification leaves $2^{\sigma n}$; occurrence bound $B=O_k((\sigma^{-1}\log\frac2\sigma)^{k-1})$ | 08 Prop. 2.4, §8 |
| 3 | occurrences and width after Tovey | the counting algorithm needs bounded $b,w$ | 08 §3 |
| $\lambda,\nu$ | $\lceil\log n\rceil$, $\lceil\log\lambda\rceil$ | scale parameters | 08 §7 |
| $r$ | $\lfloor\sqrt\lambda\rfloor$ groups | group cost $2/r\le3/\sqrt\lambda$ | (P8) |
| $Q$ | $\lfloor\lambda/(4\nu)\rfloor$ nodes | keeps $\lvert\mathcal B\rvert=r^{2Q}\le2^{\lambda/4}$, so $12\lvert\mathcal B\rvert^2<n$ | (P3) |
| $R$ | $2\lceil\sqrt\lambda\rceil$ (even) | $R^2/(8Q)\ge2\nu$: the defect is $\le e^{-2\nu}\approx\lambda^{-2.9}$ | (P2), (P5) |
| $p_r$ | $\frac{r-1}{2r}\in[\frac14,\frac12)$ | antisymmetric fraction of local vectors: $\binom r2/r^2$ | 08 Lemma 5.1(3) |
| window | $R$ levels from $\lceil Qp_r\rceil-R/2$ | contains $\lvert j-Qp_r\rvert\le R/4$ | 08 Cor. 5.6 |
| $\gamma_k$ | in $\mathbb Z[\frac12]$ | Krawtchouk expansion of the interpolant; units mod odd $p$ | 08 Lemma 5.3 |
| $p$ | least prime $>Q$, $\le2Q$ | node differences are units; $p\nmid2\cdot Q!$ | (P0) |
| $\mu$ | $n+1$ | $p^\mu\ge3^{n+1}>2^n\ge\#F$ | 08 Theorem 7.1 |
| $q$ | $\lfloor dR'/Q\rfloor\le xn$, $x=3(R-1)/Q\le25\nu/\sqrt\lambda$ | Taylor depth per node; $\deg H<QJ$ | (P6) |
| $c_*$ | $\le24ne^{-2\nu}$ measurement rows | $2hc$ with $h\le2M/\lvert\mathcal B\rvert$, $c\le2\lvert\mathcal B\rvert e^{-2\nu}$, $M\le3n$ | 08 Theorem 7.1 |
| 13 | constant in $13n\nu^2/\sqrt\lambda$ | $12.5\nu^2$ (Taylor) $+\,4$ (the rest) $\le13\nu^2$ | 08 Theorem 7.1 |
| $2^{17}+1$ | least $\lambda$ for the stated inequalities | only (P2) fails just below | 08 Lemma 7.2 |
| $9{,}297{,}119$ | least $\lambda$ where the bound beats $2^n$ | the exact cost drops below 1 | 08 Remark after Thm 7.1 |
| $\approx7\cdot10^{10}$ | the same for the paper's parameters | cruder constants | [`results/parameters.json`](results/parameters.json) |
| $\rho(n)$ | $n(\log\log n)^2/\sqrt{\log n}$ | Part II: no Tovey loss in $k$; the Lean `runtimeScale` | 09 Theorem 9.1; 06 |
| $64$, $s^{9/10}$, $u\asymp\log s/2^d$ | App. D's groups, degree cap, depth profile | geometric decay $(\sqrt2/8)^d$ of the recursion | 10 |

**Dictionary.**

| This proof | Standard name |
|---|---|
| deformation counting | polynomial interpolation of a parametrized sum: Lagrange, Hermite, and the Chinese remainder theorem for $\prod(T-t)^J$ |
| Boolean substitution | Möbius inversion on the subset lattice, i.e. inclusion–exclusion |
| separation by measurements | split-and-list: enumerate the joint values of a few linear forms, then sum each part independently |
| swap algebra | the group algebra of $(\mathbb Z/2)^Q$, i.e. Fourier analysis on the Boolean cube |
| Krawtchouk polynomials | the zonal spherical functions of the cube |
| window design | binomial concentration (Hoeffding); approximate degree of parity |
| design floor | Smolensky's dimension argument (1987) |
| parity lift (Part II) | Beigel–Tarui / Toda modulus amplification, here as an explicit idempotent |
| App. D | tensor-network contraction with $\varpi$-adic precision; Hermite interpolation with valuation loss |

---

## Level 3: reading routes

### Route A: the shortest complete refutation of ETH
- **Proof:** [08](08_GAPLESS_PROOF.md), self-contained apart from Bertrand's postulate, the Chinese remainder theorem, and RAM-to-Turing simulation.
- **Reading path:** `P.ETH.I` in [`PROOF_MAP.md`](../../../PROOF_MAP.md), seven node cards.
- **Script:** `verify_gapless_chain.py`, which runs the algorithm of §4–§7 on 16 formulas against brute force.

### Route B: the formally verified bound
- **Proof:** [09](09_GAPLESS_SHARP_EXPONENT.md), reusing 08's deformation counting with an algebra of coefficients.
- **Audit:** [06](06_LEAN_STATEMENT_AUDIT.md) §6.A. The Lean proof (3,091 modules) builds with only the standard axioms, and its statement is byte-identical to the challenge file.

### Route C: the headline exponent
- **Proof:** [10](10_GAPLESS_TENSOR_CONTRACTION.md), every step of App. D re-derived. One constant (D.5) holds only from $\log_2\lambda\approx500$, while the proof needs it only from $\approx26$.
- **Script:** `verify_tensor_contraction.py`.

### Route D: the distilled version and its limits
- **One theorem:** [13](13_UNIFIED_DEFORMATION.md) covers Parts I–III. Interpolation through local nodes has two instances (formal or 2-adic neighbourhoods).
- **The floors of both architectures:** [07](07_TOOLKIT.md) E4; [14](14_DESIGN_FLOOR.md) for graded designs; GPT [15](../../gpt/eth/15_UNGRADED_SPLINE_BOUND.md)–[16](../../gpt/eth/16_MEASUREMENT_TRANSFER_AND_FLOOR.md) for all group-algebra designs; [18](18_AGGREGATE_TRANSFER.md) for the measurement transfer under the algorithm's own schedules; [13](13_UNIFIED_DEFORMATION.md) Theorem 5.4 for App. D.

### Route E: the paper's own order

| Paper | Content | Here |
|---|---|---|
| §1 | Theorem 1.1, Corollaries 1.2–1.4 | [01](01_VERDICT.md), [11](11_GAPLESS_COROLLARIES.md) |
| §2 (and §2.8) | overview; OR gates and Tovey cycles | 08 §3 |
| §3 | deformation counting: measured partitions (Def. 3.1), Lemmas 3.2–3.3, the cost (3.15) | 08 §4 |
| §4 | swap interpolation (Lemmas 4.1–4.3, rows (4.10)), transfer (§4.5), balanced schedules (Lemma 4.5), parameters (Lemma 4.8) | 08 §5–§7 |
| §5–6 | sparse $\mathbb F_2$ systems, Prop. 6.4 | 09 |
| App. A | disjoint sparsification | 08 §2 |
| App. B | bit lengths; the parity lift | arithmetic mod $p^\mu$ (08 §7); closed-form lift (09 Lemma 2.1) |
| App. C | Corollaries C.1–C.2 | 11 |
| App. D | Boolean tensor contraction (Theorem D.1) | 10 |

### Notation traps

1. **$D$ has two meanings.** In 08 §5–§6 and 07, $D=r^{2Q}$ is the number of labels. In 14, GPT 15–17 and 18, $D$ is a polynomial degree, and labels are $\lvert\mathcal B\rvert$ or $L$.
2. **Three degrees.**
   - $R$: the design has degree $<R$.
   - $R'=R-1$: the actual deformation degree.
   - $q=\lfloor dR'/Q\rfloor$: the Taylor depth.
3. **$M$ has three meanings:** the number of incidences (08 §6); the number of wires (10, 13); and $M_0$, a matrix.
4. **$h$ has three meanings:** the number of copies (08 §6); $h=Y(Y-1)$ in the parity lift (09); and the frequency scale in the Kakeya files.
5. **$p$ has three meanings:** the prime above $Q$; the bias $p_r$ (or $p$ in 14); and nothing to do with $L^p$.
6. **$n$ versus $N$.** After sparsification and Tovey, the counting algorithm runs on $N=O(Bn)$ variables. 08 §7 calls that $n$ again.
7. **$\log$ is $\log_2$** throughout ETH; $\ln$ is natural.
8. **Symmetric versus graded versus group-algebra designs** are nested classes (07 E1, 14 Definition 1.1, GPT 15 §4). The floors hold for the largest of them.
9. **The cost model is an enumeration bound.** The floors bound the expression of 08 Theorem 4.6, not every possible algorithm (GPT 16 §2).

### Checkpoints

| Move | Question to check yourself | Script |
|---|---|---|
| 1 | Why do the leaves' counts add up to $\#F$ exactly? | `verify_gapless_chain.py` S, R |
| 2–3 | Why is $P(t)$ local at the node $t$, but not at other nodes? | Part K (`locality_ok`), `verify_worked_example.py` |
| 4 | Why is the rank of $I-M_0$ a binomial tail? | `verify_distilled_swap_formula.py`, `verify_design_floor.py` |
| 5 | Why do separated jets equal direct jets? | Part K, `verify_worked_example.py` |
| 6 | Why does $J-1$ jets per node fail? | `verify_worked_example.py` (control) |
| 7 | Where does each of the three costs come from, and why $13$? | Part P, `verify_parameters.py`, `verify_cost_floor.py` |
| floor | Why can no group-algebra design do better? | `verify_design_floor.py`, GPT `verify_spline_bound.py`, `verify_aggregate_transfer.py` |

---

## Questions readers ask

**Is it really correct?**
- **Formally.** The formally verified bound has a machine-checked proof that builds with the standard axioms only, and its statement matches the challenge file byte for byte ([06](06_LEAN_STATEMENT_AUDIT.md)).
- **By hand.** Every lemma of the ETH-refuting chain is re-derived in [08](08_GAPLESS_PROOF.md), and the algorithm returns exact counts in tests.
- **Subjective confidences:** ≈99% for ETH false; ≈95% for the headline exponent ([01](01_VERDICT.md)).

**Why is the result so slow to kick in?** The constants were never optimized. Our parameters beat $2^n$ from $\log_2n\approx9.3\cdot10^6$, the paper's from $\approx7\cdot10^{10}$, while the cost floor itself drops below 1 near $\log_2n\approx2^{10}$ (07 E4 remark). The mechanism is asymptotic.

**Does it say anything about P versus NP?** No. A $2^{n/\sqrt{\log n}}$ algorithm is far from polynomial.

**What exactly is "local"?** At a node $t$, every row of group $g$ reads only variables whose schedule puts them in group $g$ at $t$. The only exception is the defect, whose finitely many linear readouts are enumerated as measurements.

**Why $\sqrt{\log n}$ and not better?** Three costs pull against each other:
- small groups need many labels, hence few nodes;
- few nodes need high degree;
- high degree needs a small defect, which costs measurements unless the degree reaches the parity-agreement degree $\sqrt{Q\log(1/\varepsilon)}$.

The floors prove that this balance is optimal inside the architecture ([18](18_AGGREGATE_TRANSFER.md) Corollary 18.3).

**How was it found?** The paper is an OpenAI release; no account of its discovery has been published.

---

## What to reuse

| ID | Reusable statement | Where |
|---|---|---|
| E1 | Commuting-involution interpolation: a path of degree $<R$ from $I$ to elements local at $Q$ nodes, with binomial-tail defect | 07 E1, 08 §5 |
| E2 | $p$-adic cross-equation truncation, with a closed-form idempotent | 07 E2, 09 Lemma 2.1, 12 |
| E3 | Deformation counting: exact sums of products of local factors in $2^{o(n)}$ time | 07 E3; partition functions (07) |
| E4 | The cost floor of the method | 07 E4, 14, 18 |
| E6 | Interpolation through local nodes (one theorem for Parts I–III) | 13 |
| E8 | Parity completion in cube splines | GPT 15 |
| E9 | Measurement transfer: support injection, under-filled labels | GPT 16, 18 |
| E10 | Trace capacity and approximate locality | GPT 17 |

The full table, with node IDs and statuses, is [`ROADMAP.md`](../../../ROADMAP.md) §4.
