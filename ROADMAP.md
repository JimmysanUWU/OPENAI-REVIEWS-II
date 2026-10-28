# Mathematical roadmap

This is a map of every result reviewed in this repository, written for AI agents (GPT, Claude) and for people. Its machine-readable twin is [`roadmap.json`](roadmap.json). The tables and graphs between `GENERATED` markers are produced from that file by [`tools/check_roadmap.py`](tools/check_roadmap.py), which also validates it. Updated 28 October 2026.

## 0. Quick start for AI agents

**New here or after a context reset?** Read [`START_HERE.md`](START_HERE.md) first, then [`PROOF_MAP.md`](PROOF_MAP.md), which is this DAG as readable dependency-ordered paths with a proof sketch for every node. Verify everything with `python3 tools/verify_all.py` (about 3 minutes).

1. **One node per result.** Each result has a stable ID such as `ETH.II.LIFT` or `KAK.E.BAL`. The node records:
   - a precise statement (ASCII math) and a status;
   - the file and label where it is proved;
   - a **proof sketch** from the nodes it uses, and any alternative proofs;
   - the nodes its proof uses, and the scripts that check it;
   - its Lean name, if that exact statement is formally verified.
2. **Checking a claim.**
   - Open the node's `where` file at the given label.
   - Run its `checks` scripts. They are deterministic, print `ALL OK`, and write JSON to `results/`.
3. **Seeing what a result rests on.** Run `python3 tools/check_roadmap.py --node KAK.T12`. It prints the node, every node it transitively depends on, and its weakest prose link.
4. **Contributing.** Pick a target in §6 and follow the protocol in §7. Never renumber or rename a node.
5. **Standards.** See [`REVIEW_STANDARD.md`](REVIEW_STANDARD.md): definition → theorem → proof; narrowest verdicts; subjective confidences labelled as such.

**Reading orders.**

| Goal | ETH | Kakeya |
|---|---|---|
| Verdict in 5 minutes | [eth/01](reviews/claude/eth/01_VERDICT.md) | [kakeya/01](reviews/claude/kakeya/01_VERDICT.md) |
| Shortest complete proof | path *ETH is false* in [PROOF_MAP](PROOF_MAP.md#p-eth-i), full text [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) | path *Shortest complete proof* in [PROOF_MAP](PROOF_MAP.md#p-kak-planar), full text [kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md) |
| The distilled general proof | — | path *The distilled route* in [PROOF_MAP](PROOF_MAP.md#p-kak-engine), full text [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) |
| Everything the paper claims | eth/08 → [09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) → [10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) → [11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) | kakeya/08 → [09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) → [10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) |
| Reusable tools | [eth/07](reviews/claude/eth/07_TOOLKIT.md), [eth/12](reviews/claude/eth/12_FP_COUNTING.md) | [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) |
| Formal proofs | [eth/06](reviews/claude/eth/06_LEAN_STATEMENT_AUDIT.md) | [kakeya/06](reviews/claude/kakeya/06_LEAN_AUDIT.md) |

## 1. Status codes

| Code | Meaning |
|---|---|
| `FORMAL` | The exact statement is a theorem of OpenAI's Lean proof. That proof was built locally, uses only `propext`, `Classical.choice` and `Quot.sound`, and its elaborated statement was audited. Internal lemmas are never marked `FORMAL`; the optional field `lean_area` names the closest Lean module instead. |
| `PROVED` | A complete prose proof is in this repository. Textbook theorems (spectral theorem, Hoeffding, Riesz–Thorin, Mihlin, …) are cited by name. |
| `PROVED_CITED` | A complete prose proof, except for named non-elementary cited theorems whose hypotheses are checked (Evans–Krylov, Tao's ε-removal, Williams's algorithm, …). |
| `CHECKED` | Key identities, or the algorithm itself, are checked exactly by a script with negative controls. Finite checks confirm identities; they do not prove asymptotic claims. |
| `PARTIAL` | Exponents, new identities and the structure are checked. Some analytic estimates were followed but not rewritten line by line. |
| `HEURISTIC` | A blueprint, not yet a theorem. |
| `NEW` | Not in the source paper; proved here. |

The **prose chain** column gives the weakest prose link among a node and everything it depends on:
- `closed` means every link is `PROVED`;
- `cites X, Y` names the `PROVED_CITED` links.

The checker enforces that a node proved in prose depends only on nodes proved in prose.

<!-- BEGIN GENERATED: summary -->
| Project | Results | FORMAL | PROVED | PROVED_CITED | PARTIAL | HEURISTIC | NEW | CHECKED |
|---|---|---|---|---|---|---|---|---|
| ETH | 30 | 1 | 28 | 2 | 0 | 0 | 4 | 22 |
| KAK | 34 | 6 | 28 | 6 | 0 | 0 | 8 | 18 |

Targets: 14 open, 0 in progress, 4 done.
<!-- END GENERATED: summary -->

## 2. ETH: *The Exponential Time Hypothesis Is False*

**The question.** Is #$k$-SAT solvable in $2^{o(n)}$ time? Yes. There are three nested routes, and each is complete here.

| Route | Bound | Key idea | Where |
|---|---|---|---|
| I | $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}$ | Deform the substitution $x\mapsto P(T)x$ so the count becomes a polynomial in $T$. At integer nodes the sum splits into independent groups. Recover $T=0$ by interpolation. | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) |
| II | $2^{O_k(n(\log\log n)^2/\sqrt{\log n})}$ (formally verified) | Reduce to sparse $\mathbb F_2$ systems with a Tanner code and greedy rows. Lift parity to $\mathbb Z/2^\ell$ at linear degree. Count with algebra coefficients. | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) |
| III | $2^{O_k(n\log\log n/\sqrt{\log n})}$; $2^{O(n/\sqrt{\log n})}$ at fixed density | Contract a Boolean tensor network by a 2-adic Cauchy-kernel deformation, recursing on modes. | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) |

All three routes share one deformation principle (`ETH.I.DEFCOUNT`): find a polynomial deformation that is local at many nodes, then pay its degree in interpolation. Route III replaces degree with 2-adic precision. Target `T.ETH.2` asks for a single theorem covering both.

### 2.1 Dependency graph

An arrow $A\to B$ means that the proof of $B$ uses $A$. Colours:
- green: formally verified statement;
- blue: proved;
- yellow: proved modulo a named cited theorem;
- orange: partially verified;
- purple dashed: heuristic.

<!-- BEGIN GENERATED: eth-graph -->
```mermaid
flowchart TD
  subgraph ETH_g0["Part I: subexponential #35;k-SAT"]
    ETH_I_SPARS["ETH.I.SPARS<br/>Disjoint sparsification"]
    ETH_I_TOVEY["ETH.I.TOVEY<br/>OR gates and Tovey cycles"]
    ETH_I_DEFCOUNT["ETH.I.DEFCOUNT<br/>Deformation counting"]
    ETH_I_SWAP["ETH.I.SWAP<br/>Swap interpolation in closed form"]
    ETH_I_TRANSFER["ETH.I.TRANSFER<br/>Transfer to a formula; balanced schedules"]
    ETH_I_PARAM["ETH.I.PARAM<br/>Counting at #60;=3 occurrences, explicit parameters"]
    ETH_I_MAIN["ETH.I.MAIN<br/>ETH is false (Part I)"]
  end
  subgraph ETH_g1["Part II: the formally verified exponent rho(n)"]
    ETH_II_LIFT["ETH.II.LIFT<br/>Parity lift in closed form; cross-equation truncation"]
    ETH_II_ZEROS["ETH.II.ZEROS<br/>Counting common zeros of sparse F_2 systems"]
    ETH_II_TANNER["ETH.II.TANNER<br/>Initial tests from a Tanner code"]
    ETH_II_GREEDY["ETH.II.GREEDY<br/>Greedy rows by exact counts"]
    ETH_II_MAIN["ETH.II.MAIN<br/>#35;k-SAT in 2^O(rho(n)) (formally verified bound)"]
  end
  subgraph ETH_g2["Part III: App. D, tensor contraction"]
    ETH_III_RING["ETH.III.RING<br/>Coefficient ring and interpolation points"]
    ETH_III_CAUCHY["ETH.III.CAUCHY<br/>Cauchy kernel"]
    ETH_III_HERMITE["ETH.III.HERMITE<br/>Tapered 2-adic Hermite recovery"]
    ETH_III_PACK["ETH.III.PACK<br/>Packing edges into packets"]
    ETH_III_CHILD["ETH.III.CHILD<br/>Modes and the child identity"]
    ETH_III_TENSOR["ETH.III.TENSOR<br/>Boolean tensor contraction"]
    ETH_III_FIXED["ETH.III.FIXED<br/>Fixed clause density"]
    ETH_III_SPARSEZ["ETH.III.SPARSEZ<br/>Sparse systems by contraction"]
    ETH_III_MAIN["ETH.III.MAIN<br/>Headline Theorem 1.1"]
  end
  subgraph ETH_g3["Corollaries"]
    ETH_C_MAXSAT["ETH.C.MAXSAT<br/>Max-k-SAT approximation"]
    ETH_C_GRAPH["ETH.C.GRAPH<br/>Independent sets, colourings, perfect matchings"]
    ETH_C_CIRC1["ETH.C.CIRC1<br/>E^NP not in SIZE(O(n))"]
    ETH_C_CIRC2["ETH.C.CIRC2<br/>E^NP circuit lower bound c n sqrt(log n)"]
    ETH_C_LOCAL["ETH.C.LOCAL<br/>Inverting maps of constant output locality"]
  end
  subgraph ETH_g4["New results"]
    ETH_N_LIFTP["ETH.N.LIFTP<br/>p-adic lift for every prime (E2) with a closed form (new)"]
    ETH_N_FQ["ETH.N.FQ<br/>Counting zeros of sparse systems over any finite field (new) (new)"]
    ETH_K_E3["ETH.K.E3<br/>Exact partition functions of bounded-degree factor graphs (new)"]
    ETH_K_E4["ETH.K.E4<br/>Cost floor of the method (new)"]
  end
  ETH_I_SWAP --> ETH_I_TRANSFER
  ETH_I_DEFCOUNT --> ETH_I_PARAM
  ETH_I_SWAP --> ETH_I_PARAM
  ETH_I_TRANSFER --> ETH_I_PARAM
  ETH_I_SPARS --> ETH_I_MAIN
  ETH_I_TOVEY --> ETH_I_MAIN
  ETH_I_PARAM --> ETH_I_MAIN
  ETH_I_DEFCOUNT --> ETH_II_ZEROS
  ETH_I_SWAP --> ETH_II_ZEROS
  ETH_I_TRANSFER --> ETH_II_ZEROS
  ETH_II_LIFT --> ETH_II_ZEROS
  ETH_II_TANNER --> ETH_II_GREEDY
  ETH_I_SPARS --> ETH_II_MAIN
  ETH_II_ZEROS --> ETH_II_MAIN
  ETH_II_TANNER --> ETH_II_MAIN
  ETH_II_GREEDY --> ETH_II_MAIN
  ETH_III_CAUCHY --> ETH_III_CHILD
  ETH_III_PACK --> ETH_III_CHILD
  ETH_III_RING --> ETH_III_TENSOR
  ETH_III_CAUCHY --> ETH_III_TENSOR
  ETH_III_HERMITE --> ETH_III_TENSOR
  ETH_III_PACK --> ETH_III_TENSOR
  ETH_III_CHILD --> ETH_III_TENSOR
  ETH_III_TENSOR --> ETH_III_FIXED
  ETH_III_TENSOR --> ETH_III_SPARSEZ
  ETH_II_LIFT --> ETH_III_SPARSEZ
  ETH_III_FIXED --> ETH_III_MAIN
  ETH_III_SPARSEZ --> ETH_III_MAIN
  ETH_I_SPARS --> ETH_III_MAIN
  ETH_II_TANNER --> ETH_III_MAIN
  ETH_II_GREEDY --> ETH_III_MAIN
  ETH_II_MAIN --> ETH_C_MAXSAT
  ETH_I_MAIN --> ETH_C_GRAPH
  ETH_II_MAIN --> ETH_C_CIRC1
  ETH_III_FIXED --> ETH_C_CIRC2
  ETH_II_MAIN --> ETH_C_LOCAL
  ETH_II_LIFT --> ETH_N_LIFTP
  ETH_N_LIFTP --> ETH_N_FQ
  ETH_II_ZEROS --> ETH_N_FQ
  ETH_I_DEFCOUNT --> ETH_N_FQ
  ETH_I_DEFCOUNT --> ETH_K_E3
  ETH_I_SWAP --> ETH_K_E3
  ETH_I_TRANSFER --> ETH_K_E3
  ETH_I_SWAP --> ETH_K_E4
  class ETH_II_MAIN formal
  class ETH_I_SPARS,ETH_I_TOVEY,ETH_I_DEFCOUNT,ETH_I_SWAP,ETH_I_TRANSFER,ETH_I_PARAM,ETH_I_MAIN,ETH_II_LIFT,ETH_II_ZEROS,ETH_II_TANNER,ETH_II_GREEDY,ETH_III_RING,ETH_III_CAUCHY,ETH_III_HERMITE,ETH_III_PACK,ETH_III_CHILD,ETH_III_TENSOR,ETH_III_FIXED,ETH_III_SPARSEZ,ETH_III_MAIN,ETH_C_MAXSAT,ETH_C_GRAPH,ETH_C_LOCAL,ETH_N_LIFTP,ETH_N_FQ,ETH_K_E3,ETH_K_E4 proved
  class ETH_C_CIRC1,ETH_C_CIRC2 cited
  classDef formal fill:#c8e6c9,stroke:#2e7d32,color:#000
  classDef proved fill:#bbdefb,stroke:#1565c0,color:#000
  classDef cited fill:#fff9c4,stroke:#f9a825,color:#000
  classDef partial fill:#ffccbc,stroke:#d84315,color:#000
  classDef heuristic fill:#e1bee7,stroke:#6a1b9a,color:#000,stroke-dasharray:4 3
```
<!-- END GENERATED: eth-graph -->

### 2.2 Results

<!-- BEGIN GENERATED: eth-nodes -->
| ID | Result | Status | Prose chain | Depends on | Where | Checks | Lean |
|---|---|---|---|---|---|---|---|
| `ETH.I.SPARS` | **Disjoint sparsification.** For 0&lt;sigma&lt;=1 a k-CNF on n variables is split by a branching tree into at most (C_M n+1)^2 2^{sigma n} normalized CNFs of width &lt;=k, whose satisfying sets partition Sat(F), each variable occurring in at most B=O_k((sigma^{-1}log(2/sigma))^{k-1}) clauses. | PROVED, CHECKED | closed | — | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §2, Prop. 2.4, Lemma 2.5 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py) | — |
| `ETH.I.TOVEY` | **OR gates and Tovey cycles.** A CNF of width &lt;=k with occurrence bound B gives, in polynomial time, a CNF L' of width &lt;=3 with every variable in 1..3 clauses and N&lt;=6kBn variables, and z with #L=2^z #L'. | PROVED, CHECKED | closed | — | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §3, Prop. 3.3 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py) | — |
| `ETH.I.DEFCOUNT` | **Deformation counting.** Given a deformation of degree &lt;=R' with Q&lt;p nodes whose rows are local for measured partitions, the sum over {0,1}^n of a product of factors (clauses, or any factors with values in a free Z/p^mu-algebra, degree measured in the algebra) is computed with Q\*C_q(n)\*(2beta+1)^{c_\*}\*2r\*2^s\*poly operations, q=floor(dR'/Q). | PROVED, CHECKED | closed | — | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4, Thm 4.6; [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §1, Thm 1.1 (general factors) | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) | — |
| `ETH.I.SWAP` | **Swap interpolation in closed form.** P(T)=I-M_0+sum_k gamma_k sum_{\|S\|=k} prod_{s in S}(1-T/s) Pi_{S^c} has P(0)=I, degree &lt;=R-1, coefficients in Z[1/(2Q!)], and P(t) is local on ker C at every node t; c/D &lt;= 2exp(-R^2/(8Q)). | PROVED, CHECKED | closed | — | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §5, Thm 5.4, Cor. 5.6; [eth/07](reviews/claude/eth/07_TOOLKIT.md) E1 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_distilled_swap_formula](reviews/claude/eth/scripts/verify_distilled_swap_formula.py), [verify_replacement_forms](reviews/claude/eth/scripts/verify_replacement_forms.py) | — |
| `ETH.I.TRANSFER` | **Transfer to a formula; balanced schedules.** Labels from derandomized schedules turn the swap rows into a deformation with measurement rows of norm &lt;=2^Q, at most 2hc rows per node, h&lt;=2M/D, groups &lt;=2n/r. | PROVED, CHECKED | closed | `ETH.I.SWAP` | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §6, Lemmas 6.1-6.2 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py) | — |
| `ETH.I.PARAM` | **Counting at &lt;=3 occurrences, explicit parameters.** A CNF of width &lt;=3 with occurrences in 1..3 is counted in poly(n)\*2^{13 n nu^2/sqrt(lambda)} for n&gt;2^{2^17}; the exact bound is nontrivial for lambda&gt;9,297,119. | PROVED, CHECKED | closed | `ETH.I.DEFCOUNT`, `ETH.I.SWAP`, `ETH.I.TRANSFER` | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §7, Lemma 7.2, Thm 7.1 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py) | — |
| `ETH.I.MAIN` | **ETH is false (Part I).** #k-SAT is computable in 2^{O_k(n(log log n)^{k+1}/(log n)^{1/(2k)})} poly(n+m) time and polynomial space; hence ETH, SETH, #ETH fail. *Note:* Also implied by the formally verified ETH.II.MAIN. | PROVED, CHECKED | closed | `ETH.I.SPARS`, `ETH.I.TOVEY`, `ETH.I.PARAM` | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §8 | [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_counting_pipeline](reviews/claude/eth/scripts/verify_counting_pipeline.py) | — |
| `ETH.II.LIFT` | **Parity lift in closed form; cross-equation truncation.** e_0=1-Y, e_j=(-1)^{j+1}C(2j-1,j-1)(2Y-1); p_ell=sum_{j&lt;ell} e_j h^j (h=Y^2-Y) satisfies p_ell(y)=[2\|y] mod 2^ell; S_ell(prod_i B_i mod V^ell)=prod_i [2\|y_i] mod 2^ell with retained degree d_\*(m+2(ell-1)); \|\|e_j\|\|_1&lt;=2\*4^j. *Note:* The closed form is new; it replaces App. B.1's 2^{C ell^2}. | PROVED, CHECKED | closed | — | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §2, Lemmas 2.1-2.2 | [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py), [verify_parity_lift](reviews/claude/eth/scripts/verify_parity_lift.py) | — |
| `ETH.II.ZEROS` | **Counting common zeros of sparse F_2 systems.** m&lt;=cn equations over F_2, each with at most (log n)^c monomials of degree &lt;=k: common zeros counted in 2^{O_{k,c}(rho(n))} time, poly space. | PROVED, CHECKED | closed | `ETH.I.DEFCOUNT`, `ETH.I.SWAP`, `ETH.I.TRANSFER`, `ETH.II.LIFT` | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §3, Thm 3.1 | [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) | — |
| `ETH.II.TANNER` | **Initial tests from a Tanner code.** For m&lt;=Bn there is a deterministic A^(0) with &lt;=3n rows of &lt;=D_0 ones such that every nonzero kernel vector has weight &gt;=m/P, P=20000 B^2 mu^2. | PROVED, CHECKED | closed | — | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §4, Lemma 4.1 | [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) | — |
| `ETH.II.GREEDY` | **Greedy rows by exact counts.** Appending 3n rows, each chosen among (2m)^P candidates to minimize the exact common-zero count, makes N(A)=#F. | PROVED, CHECKED | closed | `ETH.II.TANNER` | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §4, Lemma 4.2 | [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) | — |
| `ETH.II.MAIN` | **#k-SAT in 2^{O(rho(n))} (formally verified bound).** #k-SAT in (n+m)^{O_k(1)} 2^{O_k(rho(n))} time and (n+m)^{O_k(1)} space, rho(n)=n(log log n)^2/sqrt(log n). | FORMAL, PROVED, CHECKED | closed | `ETH.I.SPARS`, `ETH.II.ZEROS`, `ETH.II.TANNER`, `ETH.II.GREEDY` | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) Thm 9.1; [eth/06](reviews/claude/eth/06_LEAN_STATEMENT_AUDIT.md) §6.A | [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) | `exact_count_and_sat` |
| `ETH.III.RING` | **Coefficient ring and interpolation points.** O=Z_(2)[U]/(f), f irreducible mod 2 of degree e: v_2 is multiplicative, units are the elements nonzero mod 2, and the 2^e binary-coordinate points have unit differences. | PROVED | closed | — | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §2, Lemma 2.1 | — | — |
| `ETH.III.CAUCHY` | **Cauchy kernel.** W_d(T)=C(T)/(C'(0)(T-d)), S={-L..L-1}: W(0)=I, W(t) is the shift by t, W_d(t+alpha z) is integral and divisible by 2^{a+1-v_2(t-d)} for d != t. | PROVED, CHECKED | closed | — | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §3, Lemma 3.1, Cor. 3.2 | [verify_appendix_d](reviews/claude/eth/scripts/verify_appendix_d.py) | — |
| `ETH.III.HERMITE` | **Tapered 2-adic Hermite recovery.** With multiplicities q_t=h max(0,u-v_2(t)), F(0) mod 2^p is recovered from F(t+alpha z) mod 2^{p_t}, p_t=p+8hD_t+8M(a+1); loss v_2(E_t(0))&gt;=-6hD_t. | PROVED, CHECKED | closed | — | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §4, Lemma 4.1; [eth/07](reviews/claude/eth/07_TOOLKIT.md) E5 | [verify_appendix_d](reviews/claude/eth/scripts/verify_appendix_d.py), [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.III.PACK` | **Packing edges into packets.** Derandomized labels in [r]^a and packets of L=2^a slots with orientation bits; &lt;=3M/L packets; &lt;=3M/r ports per symbol per side. | PROVED | closed | — | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §5, Lemma 5.1 | — | — |
| `ETH.III.CHILD` | **Modes and the child identity.** At T=t+alpha zeta, crossing entries expand in tile moments sum_j j^b X_y with 2-adic weight b+1, and F(T)=sum_{beta: sum w beta&lt;P} prod_g ct(N_{g,beta}) mod 2^P, with children of &lt;=6M/r edges and exp(2 pi sqrt(MP/D_t)) mode tuples. | PROVED, CHECKED | closed | `ETH.III.CAUCHY`, `ETH.III.PACK` | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §6, Lemmas 6.1-6.4 | [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.III.TENSOR` | **Boolean tensor contraction.** A bipartite Boolean tensor network with &lt;=Cs edges, degree &lt;=s^{9/10}, polynomially many vertices and suitable entry oracles is contracted in (s+L)^{O(1)} 2^{O(s/sqrt(log s))} time, poly space. *Note:* Paper's floor u&gt;=lambda^{4/5} holds only from log2(lambda)~500; absorption inequalities (A1)-(A2) suffice from ~26. | PROVED, CHECKED | closed | `ETH.III.RING`, `ETH.III.CAUCHY`, `ETH.III.HERMITE`, `ETH.III.PACK`, `ETH.III.CHILD` | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §7, Lemmas 7.1-7.3, Thm 10.2 | [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.III.FIXED` | **Fixed clause density.** #k-SAT with m&lt;=Dn in (n+m)^{O(1)} 2^{O_{k,D}(n/sqrt(log n))} via the clause network. | PROVED, CHECKED | closed | `ETH.III.TENSOR` | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §8, Prop. 8.1 | [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.III.SPARSEZ` | **Sparse systems by contraction.** Sparse F_2 systems are counted in n^{O(1)} 2^{O(n log log n/sqrt(log n))} via incidence-subset networks with O(n log log n) wires. | PROVED, CHECKED | closed | `ETH.III.TENSOR`, `ETH.II.LIFT` | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §9, Lemma 9.1 | [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.III.MAIN` | **Headline Theorem 1.1.** #k-SAT in (n+m)^{O_k(1)} 2^{O_k(n log log n/sqrt(log n))}, and 2^{O_{k,D}(n/sqrt(log n))} when m&lt;=Dn; polynomial space. *Note:* Not formalized. | PROVED, CHECKED | closed | `ETH.III.FIXED`, `ETH.III.SPARSEZ`, `ETH.I.SPARS`, `ETH.II.TANNER`, `ETH.II.GREEDY` | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) Thm 10.1 | [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) | — |
| `ETH.C.MAXSAT` | **Max-k-SAT approximation.** OPT(F)-eps m in 2^{O_{k,eps}(rho(n))} randomized time. | PROVED | closed | `ETH.II.MAIN` | [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §1 | — | — |
| `ETH.C.GRAPH` | **Independent sets, colourings, perfect matchings.** Counts in 2^{o(v)} (independent sets, q-colourings) and 2^{o(v+e)} (perfect matchings), polynomial space. | PROVED | closed | `ETH.I.MAIN` | [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §2 | — | — |
| `ETH.C.CIRC1` | **E^NP not in SIZE(O(n)).** Williams's theorem applied to 3-SAT in 2^{o(n)}. | PROVED_CITED | cites ETH.C.CIRC1 | `ETH.II.MAIN` | [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §3, Cor. 3.1 | — | — |
| `ETH.C.CIRC2` | **E^NP circuit lower bound c n sqrt(log n).** Some L in E^NP has C_L(n)&gt;c n sqrt(log n) infinitely often (needs the fixed-density bound; the formal bound alone gives n sqrt(log n)/(log log n)^2). | PROVED_CITED | cites ETH.C.CIRC2 | `ETH.III.FIXED` | [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §3, Cor. 3.2, Remark 3.2 | — | — |
| `ETH.C.LOCAL` | **Inverting maps of constant output locality.** \|G^{-1}(b)\| and a preimage in (n+m)^{O_d(1)} 2^{O_d(rho(n))}. | PROVED | closed | `ETH.II.MAIN` | [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §4 | — | — |
| `ETH.N.LIFTP` | **p-adic lift for every prime (E2) with a closed form.** p_ell=F_ell(Y^{p-1}) and its h-adic digits (h=Y^p-Y) give S_ell(prod B_i mod V^ell)=prod [p\|y_i] mod p^ell with degree d_\*(m(p-1)+p(ell-1)); closed form (1-(2z-1)(1+4 delta)^{-1/2})/2, z=Y^{p-1}, delta=Y^{p-2}h. | NEW, PROVED, CHECKED | closed | `ETH.II.LIFT` | [eth/07](reviews/claude/eth/07_TOOLKIT.md) E2; [eth/12](reviews/claude/eth/12_FP_COUNTING.md) §1, Lemma 1.1, Remark 1.2 | [verify_padic_truncation](reviews/claude/eth/scripts/verify_padic_truncation.py), [verify_fp_counting](reviews/claude/eth/scripts/verify_fp_counting.py) | — |
| `ETH.N.FQ` | **Counting zeros of sparse systems over any finite field (new).** Fixed q: m&lt;=cn polynomials over F_q in n variables, each with at most (log n)^c terms of degree &lt;=k: common zeros counted in 2^{O(rho(n))} time, poly space. | NEW, PROVED, CHECKED | closed | `ETH.N.LIFTP`, `ETH.II.ZEROS`, `ETH.I.DEFCOUNT` | [eth/12](reviews/claude/eth/12_FP_COUNTING.md) Thm 12.1, Cor. 12.2 | [verify_fp_counting](reviews/claude/eth/scripts/verify_fp_counting.py) | — |
| `ETH.K.E3` | **Exact partition functions of bounded-degree factor graphs.** Weighted factor graphs with rational weights: exact partition function in 2^{O(n(log log n)^2/sqrt(log n))}. | NEW, PROVED | closed | `ETH.I.DEFCOUNT`, `ETH.I.SWAP`, `ETH.I.TRANSFER` | [eth/07](reviews/claude/eth/07_TOOLKIT.md) E3 | — | — |
| `ETH.K.E4` | **Cost floor of the method.** For symmetric designs, the method's exponent cannot beat (log log n)^2/sqrt(log n). | NEW, PROVED, CHECKED | closed | `ETH.I.SWAP` | [eth/07](reviews/claude/eth/07_TOOLKIT.md) E4 | [verify_cost_floor](reviews/claude/eth/scripts/verify_cost_floor.py) | — |
<!-- END GENERATED: eth-nodes -->

## 3. Kakeya: *Local smoothing for the Euclidean wave equation*

**The question.** Is Sogge's local smoothing conjecture true in every dimension? Yes. The critical estimate (Theorem 1.2) loses only $\log^{2/p_c}$. It implies the restriction (interior range), Kakeya maximal, Kakeya set, Nikodym, Bochner–Riesz and Stein square-function conjectures.

**The proof as a pipeline.** Each arrow is one node.

```
quaternionic slot calculus (KAK.Q.CALC)
  → positive conserved current: lifted (KAK.C.LIFT), then for the wave (KAK.C.WAVE)
  → fibre geometry (KAK.F.GEOM) + periodic Dirichlet problem (KAK.DIR)
  → determinant energy balance: ΔF ≤ wall flux, F ≍ X_0 (KAK.E.BAL), density to depth (KAK.E.DEPTH)
  → estimates on one fibre: data circle, observation circle, ΔF ≤ V F (KAK.FIB)
  → two-circle propagation, log² loss (KAK.PROP); only 2-dimensional auxiliary time gives logs (KAK.CAP)
  → Theorem 1.2 (KAK.T12) → Theorem 1.1 (KAK.T11) → consequences
```

### 3.1 Dependency graph

<!-- BEGIN GENERATED: kak-graph -->
```mermaid
flowchart TD
  subgraph KAK_g0["Algebra and currents"]
    KAK_Q_CALC["KAK.Q.CALC<br/>Quaternionic slot calculus"]
    KAK_C_LIFT["KAK.C.LIFT<br/>Lifted current and push-forward"]
    KAK_C_WAVE["KAK.C.WAVE<br/>The positive wave current"]
    KAK_FRAME["KAK.FRAME<br/>Cofactor divergence for arbitrary quaternionic frames (Theorem 12.1) (new)"]
  end
  subgraph KAK_g1["Fibres and energy"]
    KAK_F_GEOM["KAK.F.GEOM<br/>Fibre geometry"]
    KAK_DIR["KAK.DIR<br/>Periodic Dirichlet problem (3-slot quaternionic Monge-Ampere)"]
    KAK_E_BAL["KAK.E.BAL<br/>Determinant energy balance"]
    KAK_E_DEPTH["KAK.E.DEPTH<br/>Density to depth (rank-r background) (new)"]
    KAK_FIB["KAK.FIB<br/>Energy estimates on a fibre"]
    KAK_MAENERGY["KAK.MAENERGY<br/>The Monge-Ampere operator is the gradient of a concave energy (Theorem 12.3) (new)"]
  end
  subgraph KAK_g2["Propagation"]
    KAK_PROP["KAK.PROP<br/>Two-circle propagation"]
    KAK_CAP["KAK.CAP<br/>Capacity-dimension principle (new) (new)"]
  end
  subgraph KAK_g3["Main theorems"]
    KAK_T12["KAK.T12<br/>Critical annular estimate (Theorem 1.2)"]
    KAK_T11["KAK.T11<br/>Sogge's local smoothing conjecture (Theorem 1.1)"]
  end
  subgraph KAK_g4["Consequences"]
    KAK_R_SPHERE["KAK.R.SPHERE<br/>Spherical restriction, interior range"]
    KAK_R_CONE["KAK.R.CONE<br/>Truncated cone restriction"]
    KAK_R_PARAB["KAK.R.PARAB<br/>Paraboloid restriction and Schrodinger local smoothing"]
    KAK_KM["KAK.KM<br/>Kakeya maximal conjecture"]
    KAK_KD["KAK.KD<br/>Kakeya sets have full Hausdorff dimension"]
    KAK_NIK["KAK.NIK<br/>Nikodym maximal conjecture"]
    KAK_BR["KAK.BR<br/>Bochner-Riesz conjecture (with maximal bounds for p#62;=2)"]
    KAK_SQ["KAK.SQ<br/>Stein's square function conjecture"]
  end
  subgraph KAK_g5["Sharpness"]
    KAK_OB_FOCUS["KAK.OB.FOCUS<br/>Radial focusing (sharpness of s_+(p) and of the power 1/p_c)"]
    KAK_OB_WEAK["KAK.OB.WEAK<br/>Weak-norm obstruction"]
  end
  subgraph KAK_g6["Refinements and frameworks"]
    KAK_ANG_HANKEL["KAK.ANG.HANKEL<br/>Uniform Hankel expansion on an angular band (new)"]
    KAK_ANG_BALL["KAK.ANG.BALL<br/>Observation on a ball (shell transplant)"]
    KAK_ANG_ENV["KAK.ANG.ENV<br/>Four-kernel envelope and small source"]
    KAK_ANG_BAND["KAK.ANG.BAND<br/>One annulus, bounded angular degree"]
    KAK_APP_B["KAK.APP.B<br/>Angular endpoint"]
    KAK_APP_D["KAK.APP.D<br/>Planar case by a linear Dirichlet problem (Theorem 13.1)"]
    KAK_APP_E["KAK.APP.E<br/>Curved characteristic sheets of symmetric systems"]
    KAK_ENGINE["KAK.ENGINE<br/>The positivity engine (Theorem 11.1) (new)"]
    KAK_PENCIL["KAK.PENCIL<br/>Symmetric pencils carry engine data (Theorem 11.2) (new)"]
    KAK_ARROW["KAK.ARROW<br/>The wave is the arrow pencil (Corollary 11.4) (new)"]
  end
  KAK_Q_CALC --> KAK_C_LIFT
  KAK_C_LIFT --> KAK_C_WAVE
  KAK_Q_CALC --> KAK_F_GEOM
  KAK_Q_CALC --> KAK_DIR
  KAK_Q_CALC --> KAK_E_BAL
  KAK_C_WAVE --> KAK_E_BAL
  KAK_F_GEOM --> KAK_E_BAL
  KAK_DIR --> KAK_E_BAL
  KAK_DIR --> KAK_E_DEPTH
  KAK_E_BAL --> KAK_FIB
  KAK_E_DEPTH --> KAK_FIB
  KAK_FIB --> KAK_T12
  KAK_PROP --> KAK_T12
  KAK_C_WAVE --> KAK_T12
  KAK_ARROW -.->|alt| KAK_T12
  KAK_T12 --> KAK_T11
  KAK_T12 --> KAK_R_SPHERE
  KAK_R_SPHERE --> KAK_R_CONE
  KAK_R_SPHERE --> KAK_R_PARAB
  KAK_R_SPHERE --> KAK_KM
  KAK_KM --> KAK_KD
  KAK_KM --> KAK_NIK
  KAK_T11 --> KAK_BR
  KAK_T11 --> KAK_SQ
  KAK_T12 --> KAK_ANG_BALL
  KAK_ANG_HANKEL --> KAK_ANG_BALL
  KAK_ANG_HANKEL --> KAK_ANG_ENV
  KAK_T12 --> KAK_ANG_BAND
  KAK_ANG_BALL --> KAK_ANG_BAND
  KAK_ANG_ENV --> KAK_ANG_BAND
  KAK_T12 --> KAK_APP_B
  KAK_ANG_BALL --> KAK_APP_B
  KAK_ANG_ENV --> KAK_APP_B
  KAK_PENCIL --> KAK_APP_D
  KAK_FRAME --> KAK_APP_D
  KAK_PROP --> KAK_APP_D
  KAK_ENGINE --> KAK_APP_E
  KAK_PENCIL --> KAK_APP_E
  KAK_Q_CALC --> KAK_ENGINE
  KAK_DIR --> KAK_ENGINE
  KAK_PROP --> KAK_ENGINE
  KAK_CAP --> KAK_ENGINE
  KAK_C_LIFT --> KAK_PENCIL
  KAK_Q_CALC --> KAK_PENCIL
  KAK_APP_E --> KAK_ARROW
  KAK_Q_CALC --> KAK_FRAME
  KAK_FRAME --> KAK_MAENERGY
  class KAK_T11,KAK_R_SPHERE,KAK_KM,KAK_KD,KAK_BR,KAK_SQ formal
  class KAK_Q_CALC,KAK_C_LIFT,KAK_C_WAVE,KAK_F_GEOM,KAK_E_BAL,KAK_E_DEPTH,KAK_FIB,KAK_PROP,KAK_CAP,KAK_T12,KAK_R_CONE,KAK_NIK,KAK_OB_FOCUS,KAK_OB_WEAK,KAK_ANG_HANKEL,KAK_ANG_BALL,KAK_ANG_ENV,KAK_ANG_BAND,KAK_APP_D,KAK_APP_E,KAK_ENGINE,KAK_PENCIL,KAK_ARROW,KAK_FRAME,KAK_MAENERGY proved
  class KAK_DIR,KAK_R_PARAB,KAK_APP_B cited
  classDef formal fill:#c8e6c9,stroke:#2e7d32,color:#000
  classDef proved fill:#bbdefb,stroke:#1565c0,color:#000
  classDef cited fill:#fff9c4,stroke:#f9a825,color:#000
  classDef partial fill:#ffccbc,stroke:#d84315,color:#000
  classDef heuristic fill:#e1bee7,stroke:#6a1b9a,color:#000,stroke-dasharray:4 3
```
<!-- END GENERATED: kak-graph -->

### 3.2 Results

<!-- BEGIN GENERATED: kak-nodes -->
| ID | Result | Status | Prose chain | Depends on | Where | Checks | Lean |
|---|---|---|---|---|---|---|---|
| `KAK.Q.CALC` | **Quaternionic slot calculus.** For the 3-slot operators G, H: Hessian calculus (H(sum\|F\|^2)=4J\*J), Moore determinant and cofactor identities, cofactor divergence (sum ebar_alpha d_{a alpha} C(H(u))_{ab}=0), Green's identity and compression. | PROVED, CHECKED | closed | — | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §2, Lemmas 2.1-2.4; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K6; [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) Lemma 1.2, Lemma 2.1, Prop. 12.5 | [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py), [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py) | — |
| `KAK.C.LIFT` | **Lifted current and push-forward.** For constant-coefficient symmetric systems, B=Re sum_l Z_l (x) conj(Z_l) is PSD and conserved on (X,w_1,w_2); its density push-forward by a holomorphic flux map is conserved and quaternionic rank one. | PROVED, CHECKED | closed | `KAK.Q.CALC` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Lemmas 3.1-3.2; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm I | [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py) | — |
| `KAK.C.WAVE` | **The positive wave current.** chi=sum V^v iota_e(conj u_{kappa,v}) m_{kappa,v}, A^0=chi chi^\*: on Y!=0, A^0 is smooth, PSD and conserved, annihilated by the passive differential (dPhi)A^0=0, and A^0_TT=\|S_0 u\|^2+\|S_1 u\|^2 with explicit cap multipliers S_0,S_1 (3.5). | PROVED, CHECKED | closed | `KAK.C.LIFT` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Prop. 3.3; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K1 | [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py), [verify_wave_current_target](reviews/claude/kakeya/scripts/verify_wave_current_target.py) | — |
| `KAK.F.GEOM` | **Fibre geometry.** On the closed fibres, uniformly in a&lt;s&lt;4: H_fib(r) &gt;= c s^-2 I, G_fib G_fib^\* &lt;= H_fib, volume element ~ a^2 s^{2d}, det H_fib ~ a^-2 s^-2d; the full H(r)&gt;=0 has a one-dimensional quaternionic kernel. | PROVED | closed | `KAK.Q.CALC` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §4, Lemma 4.1 | — | — |
| `KAK.DIR` | **Periodic Dirichlet problem (3-slot quaternionic Monge-Ampere).** det H(B_0+u)=g, H(B_0+u)&gt;0, u=0 on the wall: unique smooth solution with br&lt;=u&lt;=b'r. *Note:* Comparison, barriers and the boundary C^2 estimate are proved; existence uses cited Schauder theory and the Evans-Krylov/Krylov estimates (GT 6.14, 6.19, 17.14, 17.26'). OpenAI's Lean proof contains a formal existence theory in lean_area. | PROVED_CITED | cites KAK.DIR | `KAK.Q.CALC` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §5, Lemma 5.1, Thm 5.2 | — | — |
| `KAK.E.BAL` | **Determinant energy balance.** Positive conserved current + rank-one background + designed fibres: Delta_z F &lt;= wall flux and (n-1)/(n+1) X_0 &lt;= F &lt;= n/(n+1) X_0. | PROVED, CHECKED | cites KAK.DIR | `KAK.Q.CALC`, `KAK.C.WAVE`, `KAK.F.GEOM`, `KAK.DIR` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §6, Prop. 6.2, Lemma 6.3; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm II; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K2 | [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py) | — |
| `KAK.E.DEPTH` | **Density to depth (rank-r background).** Background of rank r&lt;n with eigenvalues ~lambda, f=n-r, det(B_fib+H_fib(u))=lambda^r mu, u=0 on the wall: interior -u(x) &gt;= eps rho^2 m^{1/f} if mu&gt;=m on B_rho(x) and m^{1/f}&lt;=lambda; barrier 0&lt;=u/r&lt;=C_2 P if mu&lt;=(c_1 kappa P)^f. The rank-one case (r=1) is the paper's. | NEW, PROVED | cites KAK.DIR | `KAK.DIR` | [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K3; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm III | — | — |
| `KAK.FIB` | **Energy estimates on a fibre.** For cap data at frequency ~1/h and a=Kh: F(0,tau) &lt;= C a^2 \|\|f\|\|_{p_c}^{p_c} at \|tau\|=4a; F(T,tau) &gt;= c a^2 \|\|g_T\|\|_{p_c}^{p_c} at \|tau\|=2a; Delta_z F &lt;= C s^-2 e^{-cs/h} F for 2a&lt;=s&lt;=3. | PROVED | cites KAK.DIR | `KAK.E.BAL`, `KAK.E.DEPTH` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §7, Prop. 7.1, Lemma 7.2 | — | — |
| `KAK.PROP` | **Two-circle propagation.** Delta E &lt;= V(\|tau\|) E on R x {b&lt;=\|tau\|&lt;=3}: H(0,r) &gt;= c log(r/b)/log^2(3/b) \* int_1^2 H(T,b) dT; log^{-2} sharp for scalar supersolutions. *Note:* Not used by the Lean proof, which takes a cruder route (corrected energy, power 2n+5). | PROVED, CHECKED | closed | — | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8, Lemma 8.1, Prop. 8.2; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K5 | [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py) | — |
| `KAK.CAP` | **Capacity-dimension principle (new).** Superharmonic comparison between two spheres in an auxiliary time of dimension k costs log^{-2} only for k=2; for k!=2 it loses a power. | NEW, PROVED, CHECKED | closed | — | [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K4 | [verify_capacity_dimension](reviews/claude/kakeya/scripts/verify_capacity_dimension.py) | — |
| `KAK.T12` | **Critical annular estimate (Theorem 1.2).** \|\|U(T)P_N f\|\|_{L^{p_c}(R^n x [1,2])} &lt;= C [log(2+N)]^{2/p_c} \|\|f\|\|_{p_c}. *Note:* Formalized only with the weaker power 2n+5, which suffices for every consequence. The prose chain cites the Dirichlet theory of KAK.DIR and textbook Riesz-Thorin, Mihlin, Littlewood-Paley. | PROVED | cites KAK.DIR | `KAK.FIB`, `KAK.PROP`, `KAK.C.WAVE` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 | — | — |
| `KAK.T11` | **Sogge's local smoothing conjecture (Theorem 1.1).** \|\|U(T)f\|\|_{L^p(R^n x [1,2])} &lt;= C \|\|J^{s_+(p)+eps} f\|\|_p for n&gt;=2, 2&lt;=p&lt;inf, eps&gt;0. | FORMAL, PROVED | cites KAK.DIR | `KAK.T12` | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 | — | `euclidean_local_smoothing` |
| `KAK.R.SPHERE` | **Spherical restriction, interior range.** E: L^a(S^{n-1}) -&gt; L^b(R^n) for b&gt;p_c, b&gt;(n+1)a'/(n-1). *Note:* Tao's epsilon-removal cited. | FORMAL, PROVED_CITED | cites KAK.DIR, KAK.R.SPHERE | `KAK.T12` | [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §1-2, Lemma 1.1, Cor. 2.1 | — | `spherical_restriction` |
| `KAK.R.CONE` | **Truncated cone restriction.** Same range for the truncated cone, by Nicola's slicing. | PROVED | cites KAK.DIR, KAK.R.SPHERE | `KAK.R.SPHERE` | [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §2, Cor. 2.1 | [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.R.PARAB` | **Paraboloid restriction and Schrodinger local smoothing.** Paraboloid extension in the interior range; \|\|e^{it Delta} f\|\|_{L^{p_c}} &lt;= C log^{2/p_c} \|\|f\|\|; mixed-norm bounds via Lee-Rogers-Seeger. | PROVED_CITED | cites KAK.DIR, KAK.R.PARAB, KAK.R.SPHERE | `KAK.R.SPHERE` | [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §4, Cor. 4.1 | [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.KM` | **Kakeya maximal conjecture.** \|\|K_delta f\|\|_{L^n(S^{n-1})} &lt;= C log(1/delta)^{2A} \|\|f\|\|_n. *Note:* The Lean statement is the weaker C_eps delta^-eps form; the paper and 09 prove a power of log(1/delta). | FORMAL, PROVED | cites KAK.DIR, KAK.R.SPHERE | `KAK.R.SPHERE` | [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §3-4, Lemma 3.1, Cor. 4.1 | — | `kakeya_maximal` |
| `KAK.KD` | **Kakeya sets have full Hausdorff dimension.** Every set containing a unit segment in every direction has dimension n. | FORMAL, PROVED | cites KAK.DIR, KAK.R.SPHERE | `KAK.KM` | [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §5, Thm 5.1 | — | `kakeya` |
| `KAK.NIK` | **Nikodym maximal conjecture.** \|\|N_delta f\|\|_n &lt;= C log(1/delta)^{2A} \|\|f\|\|_n, via a projective map. | PROVED | cites KAK.DIR, KAK.R.SPHERE | `KAK.KM` | [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §3, Thm 3.1 | [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.BR` | **Bochner-Riesz conjecture (with maximal bounds for p&gt;=2).** Uniform L^p bounds for alpha&gt;delta(p). | FORMAL, PROVED_CITED | cites KAK.BR, KAK.DIR | `KAK.T11` | [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §6, Thm 6.3 | — | `bochner_riesz` |
| `KAK.SQ` | **Stein's square function conjecture.** \|\|G_alpha f\|\|_p &lt;= C\|\|f\|\|_p for alpha&gt;max(1/2, n(1/2-1/p)). | FORMAL, PROVED_CITED | cites KAK.DIR, KAK.SQ | `KAK.T11` | [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §7, Thm 7.3 | — | `stein_square` |
| `KAK.OB.FOCUS` | **Radial focusing (sharpness of s_+(p) and of the power 1/p_c).** Lower bounds (8.1). | PROVED, CHECKED | closed | — | [kakeya/02](reviews/claude/kakeya/02_VERIFICATION_LEDGER.md) row for paper Prop. 8.1 | [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py) | — |
| `KAK.OB.WEAK` | **Weak-norm obstruction.** \|\|U(t)g_N\|\|_{L^{p,inf}} &gt;~ (log N)^{d(1/2-1/p)} \|\|g_N\|\|_p for 2&lt;p&lt;=p_c. | PROVED, CHECKED | closed | — | [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §1, Prop. 1.1 | [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.ANG.HANKEL` | **Uniform Hankel expansion on an angular band.** For degree &lt;= M and z &gt;= (2+M)^K: S_z F = (2pi)^alpha z^-alpha sum_{+-} e^{+-i(z - pi alpha/2)}(P_+- + E_{+-,z})F with P_+ = Id, P_- the antipodal map, and E (with log-derivatives) bounded L^1 -&gt; L^inf uniformly in M. Explicitly \|(x d_x)^k eps_mu(x)\| &lt;= e^{3/2}(2mu+2k+1)^{max(k-1,0)}(mu+1)^2/(2x) for x &gt;= (mu+1)^2. *Note:* Replaces the paper's pole-coordinate stationary phase (whose key identity lambda d_lambda I = -1/2 I(w.grad a) is also checked) by an explicit Bessel bound. | NEW, PROVED, CHECKED | closed | — | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 1.2, Proposition 1.3, Cor. 1.4 | [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) | — |
| `KAK.ANG.BALL` | **Observation on a ball (shell transplant).** For N &gt;= 2, 0 &lt; r &lt;= 1: \|\|U(t)Pi_N f\|\|_{L^p(I x B(x_0,r))} &lt;~ min{1,(Nr)^{n/p}} L(Nr)^A \|\|f\|\|_p, given Theorem 1.2 with power A. | PROVED, CHECKED | cites KAK.DIR | `KAK.T12`, `KAK.ANG.HANKEL` | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 3.2, Lemma 3.1 | [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) | — |
| `KAK.ANG.ENV` | **Four-kernel envelope and small source.** If g has angular degree &lt;= M, support in \|y\| &gt;= s_0, and N min(r,s_0) &gt;= (2+M)^K, then \|U(t)Pi_N g(r theta)\| &lt;~ sum_+- G_g(t +- r, theta) r^{-alpha} (r &lt;= 4; r^{-d/p} for r &gt; 4), with \|\|G_g\|\|_{L^p(dv dtheta)} &lt;~ \|\|g\|\|_p, uniformly in M, N. The input from \|y\| &lt; N^{eta-1} costs N^{-delta} in strong L^p. | PROVED, CHECKED | closed | `KAK.ANG.HANKEL` | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 4.2, Lemma 4.3, Definition 4.1 | [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) | — |
| `KAK.ANG.BAND` | **One annulus, bounded angular degree.** For angular degree &lt;= M: \|\|U(t)Pi_N f\|\|_p &lt;~ [L(N)^{1/p} + L(M)^A]\|\|f\|\|_p and \|\|U(t)Pi_N f\|\|_{p,inf} &lt;~ L(M)^A \|\|f\|\|_p, uniformly in M, N. | PROVED, CHECKED | cites KAK.DIR | `KAK.T12`, `KAK.ANG.BALL`, `KAK.ANG.ENV` | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Proposition 5.1 | [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) | — |
| `KAK.APP.B` | **Angular endpoint.** At p_c, for B &gt; 1/2 + 1/p: \|\|U(t)f\|\|_{L^p(R^n x I)} &lt;~ \|\|L(J)^{1/p} f\|\|_p + \|\|L(Lambda_omega)^B f\|\|_p, and \|\|U(t)f\|\|_{L^{p,inf}} &lt;~ \|\|L(Lambda_omega)^B f\|\|_p; Lambda_omega^sigma (sigma &gt; 0) may replace L(Lambda_omega)^B. *Note:* Cited: Littlewood-Paley and synthesis (Grafakos), Fefferman-Stein, vector-valued maximal on S^d (Grafakos-Liu-Yang 2009), localized spherical cut-off kernels (Narcowich-Petrushev-Ward 2006), spherical multipliers (Bonami-Clerc 1973), real interpolation; hypotheses checked in 14 §0.3. Whether B &gt; 1/2 + 1/p is sharp is open (T.KAK.13). | PROVED_CITED, CHECKED | cites KAK.APP.B, KAK.DIR | `KAK.T12`, `KAK.ANG.BALL`, `KAK.ANG.ENV` | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Theorem 14.1, §6; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 | [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) | — |
| `KAK.APP.D` | **Planar case by a linear Dirichlet problem (Theorem 13.1).** n=2: \|\|e^{-iT\|D\|} phi(D/N) f\|\|_{L^4(R^2 x [1,2])} &lt;= C log(2+N)^{1/2} \|\|f\|\|_4, hence local smoothing in R^{2+1} for all 2&lt;=p&lt;inf with eps loss. The fibre problem is linear: C(B_fib):H_fib(u)=mu, i.e. 4 sum Re(conj e_a e_b) d_{X_a X_b} u + 4\|J\|^2 Lap_W u = mu on a disk, solved by Fourier series and the Dirichlet Laplacian; the energy (1/2) int(-u)mu is controlled by a completed-square inequality A:H &gt;= D_B(H). No Monge-Ampere theory. *Note:* Uses the 2D Dirac pencil, whose chart is explicit: W is the stereographic coordinate of (P - Phi e_2)/tau on the complexified unit circle and Phi = P_2 - sqrt(tau^2 - P_1^2). For n=2 the eps-loss form is due to Guth-Wang-Zhang (2020) by a different method. | PROVED, CHECKED | closed | `KAK.PENCIL`, `KAK.FRAME`, `KAK.PROP` | [kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md) Thm 13.1, Cor. 13.2, §1-5; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 | [verify_planar_entry](reviews/claude/kakeya/scripts/verify_planar_entry.py), [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py), [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.APP.E` | **Curved characteristic sheets of symmetric systems.** For a simple extremal eigenvalue phi of a symmetric pencil with definite transverse Hessian: \|\|e^{-iT phi(D)} b(hD) f\|\|_{L^{p_c}} &lt;= C log(2/h)^{2/p_c} \|\|f\|\|_{p_c}. *Note:* Completed in 11 (chart and current quantified, wave proof reused as one abstract argument). The single-patch lower bound of App. E is not re-derived. | PROVED, CHECKED | cites KAK.DIR | `KAK.ENGINE`, `KAK.PENCIL` | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Cor. 11.3; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 | [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py), [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) | — |
| `KAK.ENGINE` | **The positivity engine (Theorem 11.1).** Engine data on a cap (phi, K', v) -- a holomorphic fibre chart (W, Phi) with the scales (G1)-(G2), a positive conserved current chi chi^\* tangent to Phi whose rows are Fourier multipliers with Gaussian symbols C (hs)^{-d/2} sigma^m exp(-c_0 sigma \|W - v(eta)\|^2) (A1)-(A3), and an invertible time row (A4) -- imply, for some aperture gamma and all 0&lt;h&lt;=1, int_1^2 \|\|e^{-iT phi(D)} f\|\|_{L^p(T^n)}^p dT &lt;= C log^2(2/h) \|\|f\|\|_p^p, p=2n/(n-1), for f with Fourier support in the gamma-cap at scale h. *Note:* Re-runs the proofs of KAK.E.BAL, KAK.E.DEPTH and KAK.FIB under the axioms; the equation enters only through (A1)-(A4). Hypothesis-usage matrix and failure modes in 11 §7. | NEW, PROVED | cites KAK.DIR | `KAK.Q.CALC`, `KAK.DIR`, `KAK.PROP`, `KAK.CAP` | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §1-5, Thm 11.1, Prop. 4.1; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K7 | — | — |
| `KAK.PENCIL` | **Symmetric pencils carry engine data (Theorem 11.2).** C(xi)=sum xi_j C_j real symmetric n x n, phi=lambda_max simple at omega with D^2 phi(omega)&gt;0 on omega^perp: with E(xi)=(1,v(xi)), S=phi I-C, L w=(C(w)e_0)_perp, B(W)=e(W)^T C e(W)/(1+W.W), the chart P=tau B(W)+Phi omega and App. E's pushed current chi_mu=J^{-1/2} sum_a (C^mu z)_a i G_a satisfy (G1)-(G2) and (A1)-(A4); D^2 phi(omega)=2 L^T S_perp^{-1} L, Dv(omega)=S_perp^{-1} L, dB(0)=2 L^T, and the window is exp(-(y/2h)(W-v)^\* S_perp (W-v)). *Note:* Conservation of the pencil current checked directly in target coordinates for 5 pencils (incl. the cubic of Ex. E.2), to 1e-38. | NEW, PROVED, CHECKED | closed | `KAK.C.LIFT`, `KAK.Q.CALC` | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Lemmas 6.1-6.4, Thm 11.2 | [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py) | — |
| `KAK.ARROW` | **The wave is the arrow pencil (Corollary 11.4).** A(xi)=[[xi_n, xi'^T],[xi', -xi_n I]] has lambda_max=\|xi\| (simple off the ray xi'=0, xi_n&lt;0), polarization v(xi)=xi'/(\|xi\|+xi_n) (stereographic) and group-velocity map B(W)=(2W,1-W.W)/(1+W.W) (inverse stereographic); the paper's null-coordinate wave system is A after the frame change (t,s_0)-&gt;(T,b) and the state scaling (C^T)^{-1/2}. Hence Theorem 1.2's cap estimate is Corollary 11.3 for two rotated arrow pencils (a second derivation, with a spatial passive coordinate). | NEW, PROVED, CHECKED | cites KAK.DIR | `KAK.APP.E` | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Cor. 11.4 | [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py) | — |
| `KAK.FRAME` | **Cofactor divergence for arbitrary quaternionic frames (Theorem 12.1).** For any quaternions f_{a,i} (any number per slot, possibly dependent) and H_f(u)_{ab}=sum f_{a,i} conj(f_{b,j}) d_{a,i} d_{b,j} u: sum_{a,i} conj(f_{a,i}) d_{a,i} C(H_f(u))_{ab}=0. Every frame calculus is the restriction H_f(u) o A = H^{(4)}(u o A) of the 4-slot calculus (Lemma 12.2). Q^{-1}:H_f is elliptic iff every slot frame is R-linearly independent, hence at most 4 real coordinates per slot (Prop. 12.4). | NEW, PROVED, CHECKED | closed | `KAK.Q.CALC` | [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §3-5, Thm 12.1, Lemma 12.2, Prop. 12.4 | [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py) | — |
| `KAK.MAENERGY` | **The Monge-Ampere operator is the gradient of a concave energy (Theorem 12.3).** Independent frames, H_f(B_0)&gt;=0, P={u: u=0 on the boundary, H_f(B_0+u)&gt;0}: E(u)=int u int_0^1 det H_f(B_0+theta u) dtheta has dE(u)[v]=int v det H_f(B_0+u) and d^2E(u)[v,v]=-int Re(G(v)^\* C(H_f(B_0+u)) G(v))&lt;0; so det H_f(B_0+u)=g is the Euler-Lagrange equation of the strictly concave E(u)-int g u, and the energy of 08 §6 is its Legendre-type transform. | NEW, PROVED | closed | `KAK.FRAME` | [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §6, Thm 12.3 | [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py) | — |
<!-- END GENERATED: kak-nodes -->

## 4. Frameworks extracted so far

These are statements that make sense outside the papers. The ranking is by expected reach; see also [`reviews/claude/SYNTHESIS.md`](reviews/claude/SYNTHESIS.md).

| Framework | One-line content | Nodes | Status |
|---|---|---|---|
| K1. Positive conserved phase-space current | Constant-coefficient symmetric systems carry an exactly conserved, positive semidefinite density on a lifted phase space. Neither Husimi nor Wigner densities have both properties. | `KAK.C.LIFT`, `KAK.C.WAVE` | proved, checked |
| K2. Determinant energy balance | A positive conserved current, a rank-one background and designed fibres give an energy with $\Delta F\le$ wall flux and $F\asymp$ the Monge–Ampère potential's integral. | `KAK.E.BAL` | proved, checked |
| E3. Deformation counting | Exact sums over $\{0,1\}^n$ of products of local factors in $2^{o(n)}$ time, whenever a low-degree deformation makes the factors local at many nodes. | `ETH.I.DEFCOUNT`, `ETH.K.E3` | proved, checked |
| E1. Commuting-involution interpolation | A path of degree $<R$ from $I$ to elements that are local at $Q$ nodes, with an exact binomial-tail defect $\le2e^{-R^2/(8Q)}$. | `ETH.I.SWAP` | proved, checked |
| E2. $p$-adic cross-equation truncation | Lift a system over $\mathbb F_p$ to $\mathbb Z/p^\ell$ at linear degree, with a closed-form idempotent. | `ETH.II.LIFT`, `ETH.N.LIFTP` | new, proved, checked |
| K4. Capacity-dimension principle | Superharmonic transport between two spheres costs $\log^{-2}$ only when the auxiliary time is 2-dimensional. | `KAK.CAP` | new, proved, checked |
| K3. Density to depth, any rank | A Monge–Ampère potential is deep where the density is large, at the exponent $1/f$ with $f$ free directions. | `KAK.E.DEPTH` | new (rank $r$), proved |
| E4. Cost floor | Symmetric swap designs cannot beat the exponent $(\log\log n)^2/\sqrt{\log n}$. | `ETH.K.E4` | new, proved, checked |
| E5. Tapered Hermite recovery | Exact interpolation from many nodes, where the precision needed tapers with the 2-adic depth of the node. | `ETH.III.HERMITE` | proved, checked |
| K5. Two-circle propagation | $\Delta E\le VE$ on a thin cylinder transports averages between two circles at cost $\log(r/b)/\log^2(3/b)$, sharp for supersolutions. | `KAK.PROP` | proved, checked |
| K6. Cofactor conservation, any frame | The Moore cofactor of a frame Hessian $\sum f_{a,i}\bar f_{b,j}\partial\partial u$ is divergence-free, for any quaternion frames. The Monge–Ampère operator is the gradient of a concave energy. | `KAK.Q.CALC`, `KAK.FRAME`, `KAK.MAENERGY` | proved, checked; generalized in [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) |
| K7. Positivity engine | Axioms on a fibre chart and a current imply the critical cap estimate with $\log^2$ loss. The axioms say nothing about the equation. Symmetric pencils satisfy them, and the wave is the arrow pencil. | `KAK.ENGINE`, `KAK.PENCIL`, `KAK.ARROW` | new, proved, checked ([kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md)) |

**One pattern behind both papers.** Each paper moves information along an auxiliary parameter and pays a capacity for the transport:
- ETH: interpolation in $T$, at the cost of the binomial bulk;
- Kakeya: superharmonic comparison in $\tau\in\mathbb R^2$, at the cost of the planar capacity.

Each paper's barrier is a sharp statement about that capacity (`ETH.K.E4`, `KAK.CAP`). Target `T.X.1` asks whether this analogy is a theorem.

## 5. The distillation programme

Verification is complete for both papers (see the verdict files). The remaining work is distillation: rebuilding the mathematics as reusable definitions and theorems with short proofs. It is ordered as follows.

1. **Kakeya, `T.KAK.1`: the positivity engine as a theorem. Done** ([kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md)).
   - Theorem 11.1 proves the critical cap estimate from axioms (G1)–(G2) and (A1)–(A4).
   - Theorem 11.2 shows that symmetric pencils satisfy them, which completes App. E.
   - Corollary 11.4 shows the wave is the arrow pencil.
   - The paraboloid follows only through the wave, by rescaling. A direct engine for it is open (`T.KAK.12`).
   - §7 is the hypothesis-usage matrix.
2. **Kakeya, `T.KAK.3`: the quaternionic slot calculus. Done** ([kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md)). Cofactor divergence holds for arbitrary quaternion frames (Theorem 12.1). The Monge–Ampère operator is the gradient of a concave energy (Theorem 12.3). Ellipticity holds iff the frames are independent (Proposition 12.4).
3. **Kakeya, `T.KAK.2`: the planar entry point. Done** ([kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md)). It gives a complete proof of sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ with $\log^{1/2}$ loss, using only linear PDE: a Dirichlet Laplacian on a disk, Fourier series, and the maximum principle. The prose chain is closed: no cited non-elementary theorem.
4. **ETH, `T.ETH.2`: one deformation theorem** covering the degree-controlled deformations of Parts I–II and the 2-adic one of Part III. Then `T.ETH.4`: explicit thresholds.

Research targets (`T.ETH.1`, `T.KAK.4`–`T.KAK.9`, `T.KAK.12`, `T.X.1`) are open problems made precise by the distillation. They are listed for anyone who wants them.

## 6. Targets

<!-- BEGIN GENERATED: targets -->
| ID | Kind | Title | Status | Owner | Builds on | Goal | Acceptance |
|---|---|---|---|---|---|---|---|
| `T.ETH.1` | research | Non-symmetric swap designs | OPEN | — | `ETH.I.SWAP`, `ETH.K.E4` | Find the least p-biased measure of the set where a real polynomial of degree &lt;R differs from parity on the cube; either prove e^{-O(R^2/Q)} for all polynomials (making the cost floor a barrier for the whole group algebra) or construct a non-symmetric design that beats the exponent (log log n)^2/sqrt(log n). | a proof or construction with a precise statement; an exact check script with a negative control; roadmap updated |
| `T.ETH.2` | distillation | A unified deformation theorem | OPEN | — | `ETH.I.DEFCOUNT`, `ETH.III.CHILD`, `ETH.III.HERMITE`, `ETH.K.E4` | One abstract statement covering the degree-controlled deformations (Parts I-II) and App. D's 2-adic deformation, with a cost functional that explains both sqrt(log) barriers. | abstract theorem with hypotheses; both instances derived as corollaries; barrier computation for each instance |
| `T.ETH.3` | formal | Formalize App. D (Theorem 1.1) in Lean | OPEN | — | `ETH.III.MAIN` | Extend OpenAI's Lean development from Prop. 6.4 to Theorem 1.1. | a Lean proof building against the pinned Mathlib with standard axioms only |
| `T.ETH.4` | distillation | Explicit thresholds for Parts II and III | OPEN | — | `ETH.II.MAIN`, `ETH.III.TENSOR` | Explicit n_0 in Theorem 9.1 and an explicit s_0 in Theorem 10.2, as was done for Part I (lambda &gt;= 2^17+1). | every inequality checked by script on its full range; thresholds stated in the docs |
| `T.KAK.1` | distillation | The positivity engine as a theorem | DONE | claude | `KAK.ENGINE`, `KAK.C.LIFT`, `KAK.E.BAL`, `KAK.E.DEPTH`, `KAK.FIB`, `KAK.PROP`, `KAK.APP.E` | State axioms on (symbol, current, passive coordinate, fibres) under which the critical estimate follows, prove the abstract theorem, and derive the wave (Thm 1.2), symmetric pencils (App. E) and the paraboloid as instances; include a hypothesis-usage matrix. Deliverable: [`11_FRAMEWORK_POSITIVITY_ENGINE.md`](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md). | axioms stated precisely; abstract proof referencing only the axioms; each instance verifies the axioms; checks for the new identities |
| `T.KAK.2` | distillation | The planar entry point | DONE | claude | `KAK.APP.D`, `KAK.C.WAVE`, `KAK.PROP` | A self-contained complete proof of sharp L^4 local smoothing in R^{2+1} with log^{1/2}, using the linear Dirichlet problem of App. D (no Monge-Ampere regularity). Deliverable: [`13_PLANAR_ENTRY_POINT.md`](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md). | every step proved or a cited textbook theorem; App. D's remaining analytic steps re-derived; checks for the planar identities |
| `T.KAK.3` | distillation | Quaternionic slot calculus as a standalone toolkit | DONE | claude | `KAK.Q.CALC` | A short self-contained account of the 3-slot calculus (G, H, Moore determinant, cofactor divergence, Green, Schur, rank-one expansion) usable outside the paper. Deliverable: [`12_QUATERNIONIC_SLOT_CALCULUS.md`](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md). | definitions and proofs; exact symbolic checks |
| `T.KAK.4` | research | A coordinate-free conserved current | OPEN | — | `KAK.C.LIFT` | Describe the lifted current on the complexified characteristic variety as a CR-geometric object. | invariant definition reproducing Prop. 3.3 |
| `T.KAK.5` | distillation | Avoid the fourth slot in the boundary C^2 estimate | OPEN | — | `KAK.DIR` | Prove the boundary second-derivative estimate of App. A with 3-slot rigid motions only, or explain why a fourth slot is needed. | proof or obstruction |
| `T.KAK.6` | research | One logarithm instead of two | OPEN | — | `KAK.PROP`, `KAK.OB.FOCUS`, `KAK.CAP` | Reach log^{1/p_c} (the lower bound) by using the energy's structure beyond superharmonicity. | improved propagation lemma with proof |
| `T.KAK.7` | research | Pencils whose size differs from the spatial dimension | OPEN | — | `KAK.PENCIL`, `KAK.ENGINE` | Extend Theorem 11.2 beyond 'state dimension = spatial dimension'. For N&gt;n the flux map is a submersion; push the lifted current forward by fibre integration (still positive, conserved and tangent) and re-derive (A3)-(A4). Connects to determinantal representations of hyperbolic polynomials. | a theorem with a new class of symbols |
| `T.KAK.8` | research | Variable coefficients: which hypothesis fails | OPEN | — | `KAK.ENGINE` | Explain the Minicozzi-Sogge counterexamples on curved metrics. By 11 §8.4 the obstruction must lie in the current axiom (A1): (A3)-(A4) are stable symbol estimates. Find the precise failure. | a precise failing hypothesis with an example |
| `T.KAK.9` | research | Degenerate symbols with a rank-r background | OPEN | — | `KAK.E.DEPTH`, `KAK.E.BAL` | Find currents tangent to r&gt;1 passive coordinates and prove the rank-r energy comparison, giving critical estimates at 2(f+1)/f. | rank-r analogue of Lemma 6.3 with proof |
| `T.KAK.10` | verification | Line-by-line verification of App. B | DONE | claude | `KAK.APP.B` | Re-derive the analytic estimates of App. B (Lemmas B.3, B.4, Prop. B.2, the summation). App. E was completed in 11 (Corollary 11.3). Deliverable: [`14_GAPLESS_ANGULAR_ENDPOINT.md`](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md). | complete derivations; KAK.APP.B upgraded to PROVED |
| `T.KAK.11` | formal | The sharp power 2/p_c in Lean | OPEN | — | `KAK.T12`, `KAK.PROP` | Formalize Theorem 1.2 with the sharp power. | Lean proof with standard axioms only |
| `T.X.1` | research | Capacity of transport as a theorem | OPEN | — | `ETH.K.E4`, `KAK.CAP` | Turn the cross-paper analogy (move information along an auxiliary parameter and pay a capacity) into a precise statement covering both papers' barriers. | a precise common statement or a reasoned negative answer |
| `T.KAK.12` | research | A direct engine for the paraboloid | OPEN | — | `KAK.ENGINE`, `KAK.ARROW`, `KAK.R.PARAB` | e^{it Delta} is not a first-order symmetric system, so it is not an instance of Theorem 11.2; its critical estimate currently comes from the wave by parabolic rescaling. Find an exactly conserved, positive current tangent to a passive coordinate for the Schrodinger flow satisfying (A1)-(A4), e.g. as a parabolic limit of arrow-pencil currents. | engine data for the paraboloid with proofs; a direct-coordinate conservation check |
| `T.KAK.13` | research | The angular threshold | OPEN | — | `KAK.APP.B`, `KAK.ANG.BALL` | Decide whether B &gt; 1/2 + 1/p in the angular endpoint is sharp. The 1/2 - 1/p comes from l^2-vs-l^p over the ~log H_m frequencies near the stationary-phase threshold, and A = 2/p from the ball estimate: a square-function form of the shell-transplant lemma, or A = 1/p (T.KAK.6), would lower it. | a lower threshold with proof, or a counterexample family |
<!-- END GENERATED: targets -->

## 7. Protocol for AI contributors

### 7.1 Claiming and finishing a target

Every new node needs a `sketch` (how it follows from the nodes it uses) and must be placed on a reading path in dependency order; the checker enforces both.
1. Set the target's `owner` to `gpt` or `claude` and its `status` to `IN_PROGRESS` in `roadmap.json`. Run `python3 tools/check_roadmap.py --write` and commit.
2. When the acceptance criteria are met:
   - set `status` to `DONE` and `deliverable` to the file;
   - add or upgrade the nodes the work proves (new IDs, `where`, `checks`, `status`);
   - run the checker again.

### 7.2 Writing a mathematical document
- **Location.** Claude writes under `reviews/claude/<project>/`, GPT under `reviews/gpt/<project>/`. Number files after the existing ones.
- **Order.** Each section goes definition → statement → proof → remarks. No procedural narrative.
- **Statements.** Number every statement (`Lemma 3.2`) and write it in bold (`**Lemma 3.2 (name).**`) so the checker can resolve labels.
- **Status.** Give every claim one verdict: derivation, accepted source, exact computation, or formal build. Name every cited theorem.
- **Cross-references.** Refer to other results by node ID (`KAK.E.BAL`) as well as by file and label.

### 7.3 Writing a check script
- Use the standard library only, unless the project's scripts already use `numpy`, `sympy` or `scipy`.
- Use fixed seeds and exact arithmetic wherever possible.
- Give every check a negative control: a deliberately wrong variant that must fail.
- Print `ALL OK` and write a JSON summary to the project's `results/`. List the script in the project's `scripts/README.md`.

### 7.4 Recording agreement or a dispute
A node may carry a list `reviews` of entries `{"by": "gpt", "verdict": "agree" | "dispute" | "partial", "file": "<path to the review>"}`.
- The rendered status column shows each entry.
- A dispute should name the step and the failing inequality or identity, ideally with a script.
- The other model answers in its own review file. Neither model deletes the other's entries.

### 7.5 Commits
- Use one topic per commit, with a descriptive message.
- Run the checker and the affected scripts before pushing.

## 8. Changelog
- **2026-10-28.** Roadmap created: 56 results (30 ETH, 26 Kakeya) and 16 targets. Distillation target `T.KAK.1` started.
- **2026-10-28 (later).** `T.KAK.1` done: [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md).
  - `KAK.ENGINE` is now a theorem (was heuristic).
  - New nodes `KAK.PENCIL` and `KAK.ARROW`.
  - `KAK.APP.E` is upgraded from PARTIAL to PROVED.
  - New target `T.KAK.12`; `T.KAK.7`, `T.KAK.8` and `T.KAK.10` refined.
- **2026-10-28 (later).** `T.KAK.3` done: [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md). New nodes `KAK.FRAME` and `KAK.MAENERGY`.
- **2026-10-28 (later).** `T.KAK.2` done: [kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md). `KAK.APP.D` is upgraded from PARTIAL to PROVED, and its prose chain is closed.
- **2026-10-28 (later).** Cold-start layer:
  - [`START_HERE.md`](START_HERE.md) and [`PROOF_MAP.md`](PROOF_MAP.md);
  - proof sketches on all 60 nodes, and an alternative proof for `KAK.T12`;
  - nine reading paths, checked for dependency order and coverage;
  - a 28-term glossary;
  - `tools/verify_all.py`, which runs every script.
- **2026-10-28 (later).** `T.KAK.10` done: [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md).
  - `KAK.APP.B` is upgraded from PARTIAL to PROVED_CITED; no node is PARTIAL any more.
  - New nodes `KAK.ANG.HANKEL`, `KAK.ANG.BALL`, `KAK.ANG.ENV` and `KAK.ANG.BAND`, on the new reading path `P.KAK.ANGULAR`; `P.KAK.LIMITS` now holds only the sharpness results.
  - New target `T.KAK.13` (the angular threshold).
  - Six glossary terms added.
