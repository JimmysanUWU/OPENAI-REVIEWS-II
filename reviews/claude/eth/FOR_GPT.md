# Handoff to GPT (Claude, 27 October 2026)

The user tells me you are verifying and distilling the prose proof of *The Exponential Time Hypothesis Is False* (OpenAI, revised 25 Oct 2026). The user's priorities, in order:
1. Does it answer the question?
2. Verify rigorously, and identify the new methods and leaps.
3. Does it explain its ideas, or rely on miracles?
4. Distill the mathematics into better structure.

## Update (27 Oct): formal proof verified

OpenAI's Lean proof of the exact statement (`openai/math`, solution `OAI.Algorithms.SubexponentialSAT.Main`) builds locally against pinned Mathlib: 3,091 modules, 0 `sorry`. `exact_count_and_sat` uses only `propext`, `Classical.choice` and `Quot.sound`, and its `pp.all` statement is byte-identical to the challenge file's ([06 §6.A](06_LEAN_STATEMENT_AUDIT.md)). Prop. 6.4, and hence ETH's falsity, is formally verified. App. D (Thm 1.1's sharper exponent) is not formalized; see the 28 Oct update below.

New toolkit results are in [07](07_TOOLKIT.md): the $p$-adic truncation for all primes, weighted partition functions, and a proved cost floor.

## Update (27 Oct, later): gapless proof of the ETH-refuting chain

[08](08_GAPLESS_PROOF.md) is a complete proof, with no gaps, that #$k$-SAT runs in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}\mathrm{poly}(n+m)$ time and polynomial space. It refutes ETH and #ETH. Four simplifications over the paper:
- arithmetic mod $p^\mu$ removes App. B;
- closed-form swap rows;
- each measurement row meets at most 2 groups per node;
- explicit parameters with all inequalities proved for $\lambda\ge2^{17}+1$, giving exponent $13n\nu^2/\sqrt\lambda$, nontrivial from $\log n\approx9.3\cdot10^6$ (the paper's: $7\cdot10^{10}$).

Every lemma is checked by [`scripts/verify_gapless_chain.py`](scripts/verify_gapless_chain.py). **An independent read of 08 would be the most useful check.**

## Update (28 Oct): gapless proof of the formally verified exponent

[09](09_GAPLESS_SHARP_EXPONENT.md) proves Prop. 6.4, #$k$-SAT in $2^{O_k(\rho(n))}$ with $\rho(n)=n(\log\log n)^2/\sqrt{\log n}$ and polynomial space. This is the bound of the Lean statement. It follows §5–6 and reuses 08's machinery. What it adds:
- **The parity lift in closed form.** The $h$-adic digits of the Beigel–Tarui polynomial are $e_j=(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)$, because the lift is the idempotent $\frac12\big(1-(2Y-1)(1+4h)^{-1/2}\big)$ of $\mathbb Q[Y]/(h^\ell)$. This replaces Lemma B.1's $2^{C\ell^2}$ by $2\cdot4^j$ and explains the truncation: digit $j$ carries $h^j$, which is divisible by $2^j$ at integers.
- Deformation counting with coefficients in $(\mathbb Z/p^\mu)[V]/(V^\ell)$, with exact recovery of the $V$-coefficients, in place of App. B's rational bookkeeping.
- 08's closed-form design in place of §4's construction.

Checked by [`scripts/verify_sharp_route.py`](scripts/verify_sharp_route.py). With 08, every bound that is both proved in the paper and formalized now has a complete prose proof here.

## Update (28 Oct, later): App. D verified; gapless proof of Theorem 1.1

[10](10_GAPLESS_TENSOR_CONTRACTION.md) proves the headline Theorem 1.1: #$k$-SAT in $2^{O_k(n\log\log n/\sqrt{\log n})}$, and $2^{O(n/\sqrt{\log n})}$ at fixed density. It goes through Theorem D.1 (Boolean tensor contraction in $2^{O(s/\sqrt{\log s})}$). **Verdict: App. D is correct.** What 10 adds:
- The recursion's bookkeeping is re-derived with explicit constants: input precision $\le41Cs(d+1)2^{a_{d-1}}$, $p_t/D_t\le34M/u$, at most $e^{37M/\sqrt u}$ mode tuples.
- The paper's floor $u\ge\lambda^{4/5}$ in (D.5) holds only from $\log_2\lambda\approx500$. Two absorption inequalities suffice, and they hold from $\log_2\lambda\approx26$.
- [`scripts/verify_tensor_contraction.py`](scripts/verify_tensor_contraction.py) runs the first end-to-end checks:
  - the child identity (D.16), with brute-force children on 66 cases;
  - one full recursion level at the honest precision (Lagrange over $\mathbb Z_{(2)}[U]/(f)$ plus tapered recovery), exact in 8/8;
  - Lemma D.2's encoding, and the clause network.

With 08–10, every complexity bound in the paper has a complete prose proof here. **An independent read of 10 §6–§7 would be the most useful check.**

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
