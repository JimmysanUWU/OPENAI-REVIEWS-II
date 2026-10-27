# Handoff to GPT (Claude, 27 October 2026)

The user tells me you are verifying and distilling the prose proof of *The Exponential Time Hypothesis Is False* (OpenAI, revised 25 Oct 2026). The user's priorities, in order:
1. Does it answer the question?
2. Verify rigorously, and identify the new methods and leaps.
3. Does it explain its ideas, or rely on miracles?
4. Distill the mathematics into better structure.

## Update (27 Oct): formal proof verified

OpenAI's Lean proof of the exact statement (`openai/math`, solution `OAI.Algorithms.SubexponentialSAT.Main`) builds locally against pinned Mathlib: 3,091 modules, 0 `sorry`. `exact_count_and_sat` uses only `propext`, `Classical.choice` and `Quot.sound`, and its `pp.all` statement is byte-identical to the challenge file's ([06 §6.A](06_LEAN_STATEMENT_AUDIT.md)). Prop. 6.4, and hence ETH's falsity, is formally verified. App. D (Thm 1.1's sharper exponent) remains prose-only.

New toolkit results are in [07](07_TOOLKIT.md): the $p$-adic truncation for all primes, weighted partition functions, and a proved cost floor.

## Update (27 Oct, later): gapless proof of the ETH-refuting chain

[08](08_GAPLESS_PROOF.md) is a complete proof, with no gaps, that #$k$-SAT runs in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}\mathrm{poly}(n+m)$ time and polynomial space. It refutes ETH and #ETH. Four simplifications over the paper:
- arithmetic mod $p^\mu$ removes App. B;
- closed-form swap rows;
- each measurement row meets at most 2 groups per node;
- explicit parameters with all inequalities proved for $\lambda\ge2^{17}+1$, giving exponent $13n\nu^2/\sqrt\lambda$, nontrivial from $\log n\approx9.3\cdot10^6$ (the paper's: $7\cdot10^{10}$).

Every lemma is checked by [`scripts/verify_gapless_chain.py`](scripts/verify_gapless_chain.py). §5–6 (the sharper exponent) is covered by the Lean proof, not re-proved in prose. **An independent read of 08 would be the most useful check.**

## Where I am

- **Verdict.** Yes, it refutes ETH. The minimal chain is App. A (disjoint sparsification), §2.8 (OR gates + Tovey), §3 (Lemmas 3.2, 3.3, Prop. 3.4) and §4 (Lemmas 4.1–4.8). I re-derived it completely and found no error. Confidence 97%. See [01](01_VERDICT.md) and [02](02_VERIFICATION_LEDGER.md).
- **Exact computations** ([`scripts/`](scripts/), outputs in [`results/`](results/)):
  - Lemmas 4.1–4.3 and Remark 4.4 for 7 parameter triples. $\operatorname{rank}(I-M_0)=c_{\mathcal B}$ exactly.
  - The full §3–§4 pipeline, including Lemma 3.2 per pair $(E,U)$, (3.13) and Hermite interpolation, on 17 CNFs: always $H(0)=\#F$. Negative controls fail as they should.
  - Lemmas 5.1 and 5.2.
  - Lemma D.4, (D.10), Lemma D.5 and (D.12).
  - The parameter thresholds of (2.19) and Lemma 4.8.
- **New structure** ([05](05_DISTILLATION.md)):
  - The replacement rows have a closed form: $P(T)=I-M_0+V_0\,e_f\big((1-T/s)\Pi_s\big)_s$, with $f=\sum\gamma_kK_k$ in the Krawtchouk basis.
  - The paper's $M_0$ acts by $(-1)^jf(j)$ on the $j$-th isotypic piece, which has dimension $D\Pr[W=j]$.
  - Proposition B″: the interval choice is optimal among $\mathrm{Sym}(Q)$-invariant constructions, up to $\mathrm{poly}(Q)$.
  - Heuristic barrier: the method cannot beat $(\log\log n)^2/\sqrt{\log n}$ ([04 §4.5](04_MIRACLES_EXPLAINED.md)).
- **Sharpening S1.** Measurement rows have $\ell_1\le2^Q$, not merely $\le r^{2Q}$. This changes no exponent.

## Suggested division (to avoid duplicate work)

| Target | Status here | Suggestion |
|---|---|---|
| Appendix D recursion: Lemma D.3 packing (D.4)–(D.6), precision (D.11), child sizes and enumeration (D.13)–(D.16), Lemma D.2 wire encoding | partially checked (local lemmas exact) | **Best use of an independent pass.** It carries Thm 1.1's headline bound and Cor. 1.4. |
| Cor. 1.4 and C.1: [Wil13, Thm 6.1] and [JMV13] interfaces | not checked | Check the cited statements against their sources. |
| §6 constants (Lemma 6.1, 6.2) | re-derived | A second reading is welcome but low priority. |
| Minimal chain (§3–§4, App. A) | complete | Cross-check [05](05_DISTILLATION.md), Theorems A–E, if you produce your own blueprint. |
| Lean | statement audited ([06](06_LEAN_STATEMENT_AUDIT.md)); proof not seen | Whoever obtains the Lean project should build it and run the comparator. |

## Open questions I would value your view on

1. Is there a reference for the agreement degree of parity under a $p$-biased measure, i.e. that degree $\Omega(\sqrt{Q\log(1/\varepsilon)})$ is necessary? It would make the optimality of swap interpolation unconditional over the whole group algebra.
2. Can Appendix D's Cauchy deformation and §4's swap interpolation be stated as one "local-at-nodes deformation" theorem with a single cost functional? See [05 §5.8](05_DISTILLATION.md).
3. Missing related work: Lokshtanov–Paturi–Tamaki–Williams–Yu (SODA'17), Björklund–Kaski–Williams (ICALP'19), Dinur (SODA'21). Do you agree they should be discussed next to §5?
