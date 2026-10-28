# OPENAI-REVIEWS-II

Joint Claude–GPT review of OpenAI's mathematical releases. The workflow is the one used in [Claude-GPT-reviews](https://github.com/JimmysanUWU/Claude-GPT-reviews).

**Start here: [`ROADMAP.md`](ROADMAP.md).** It maps every result reviewed here, with stable IDs, precise statements, statuses, dependency graphs, check scripts and open targets. Its machine-readable twin [`roadmap.json`](roadmap.json) is validated by `python3 tools/check_roadmap.py`. AI agents (GPT, Claude) should read it first and follow its protocol (§7).

## User priorities (each target)

1. Does the paper answer the question?
2. Verify rigorously, and identify methods, leaps of insight and new structure; assess how radical they are.
3. Does the paper explain its ideas, or rely on unexplained "miracles"?
4. Distill the mathematics into better structure.

## Targets

| Target | Claim | Claude verdict | GPT |
|---|---|---|---|
| [**ETH**](reviews/claude/eth/01_VERDICT.md): *The Exponential Time Hypothesis Is False* (OpenAI, Oct 2026) | #$k$-SAT in $2^{O_k(n\log\log n/\sqrt{\log n})}$ time, polynomial space | ETH refutation correct. The minimal chain was re-derived, and the algorithm was run in exact arithmetic. **OpenAI's formal proof of the exact statement builds locally** (3,091 modules, 0 sorry, only standard axioms, elaborated statement byte-identical) ([06](reviews/claude/eth/06_LEAN_STATEMENT_AUDIT.md)). | pending |
| [**Kakeya / local smoothing**](reviews/claude/kakeya/01_VERDICT.md): *Local smoothing for the Euclidean wave equation* (OpenAI, 25 Oct 2026) | Sogge's local smoothing conjecture in all dimensions, which implies Bochner–Riesz, restriction, Stein square function, Kakeya maximal and Kakeya set conjectures | Answers the question. The core identities were verified. **OpenAI's formal proof builds locally** (2,261 modules, 0 errors, 0 sorry; only standard axioms; elaborated statements byte-identical to the challenge) ([06](reviews/claude/kakeya/06_LEAN_AUDIT.md)). A new method: positive currents plus quaternionic Monge–Ampère energies. | pending |

## Gapless proofs (distillation round 3)

- **ETH** ([`08_GAPLESS_PROOF.md`](reviews/claude/eth/08_GAPLESS_PROOF.md)). A complete, self-contained proof that #$k$-SAT runs in $2^{O_k(n(\log\log n)^{k+1}/(\log n)^{1/(2k)})}$ time and polynomial space, refuting ETH. The arithmetic is mod $p^\mu$, the swap rows are in closed form, and the parameters are explicit (exponent $13n\nu^2/\sqrt\lambda$). Every lemma is checked by [`verify_gapless_chain.py`](reviews/claude/eth/scripts/verify_gapless_chain.py).
- **ETH, formally verified exponent** ([`09_GAPLESS_SHARP_EXPONENT.md`](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md)). #$k$-SAT in $2^{O_k(n(\log\log n)^2/\sqrt{\log n})}$ time and polynomial space: the paper's Prop. 6.4 and the Lean statement's bound. New: the parity lift's $h$-adic digits in closed form, $e_j=(-1)^{j+1}\binom{2j-1}{j-1}(2Y-1)$, replacing App. B.1's $2^{C\ell^2}$ bound by $2\cdot4^j$. Checked by [`verify_sharp_route.py`](reviews/claude/eth/scripts/verify_sharp_route.py).
- **ETH, corollaries and a new theorem** ([`11_GAPLESS_COROLLARIES.md`](reviews/claude/eth/11_GAPLESS_COROLLARIES.md), [`12_FP_COUNTING.md`](reviews/claude/eth/12_FP_COUNTING.md)).
  - Proofs of all of the paper's corollaries. Cor. 1.4 is found to need App. D.
  - **New:** counting zeros of sparse polynomial systems over any fixed finite field in $2^{O(n(\log\log n)^2/\sqrt{\log n})}$. This is the $\mathbb F_p$ analogue of the paper's Theorem 2.2, with a closed-form $p$-adic lift. Checked by [`verify_fp_counting.py`](reviews/claude/eth/scripts/verify_fp_counting.py).
- **ETH, the headline exponent** ([`10_GAPLESS_TENSOR_CONTRACTION.md`](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md)). Theorem 1.1, #$k$-SAT in $2^{O_k(n\log\log n/\sqrt{\log n})}$ ($2^{O(n/\sqrt{\log n})}$ at fixed density), via App. D's Boolean tensor contraction. App. D is verified completely, with explicit constants for its recursion. Its child identity and one recursion level at the honest precision are checked end to end by [`verify_tensor_contraction.py`](reviews/claude/eth/scripts/verify_tensor_contraction.py).
- **Kakeya** ([`08_GAPLESS_PROOF.md`](reviews/claude/kakeya/08_GAPLESS_PROOF.md), [`09_GAPLESS_CONSEQUENCES.md`](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md)). Theorems 1.1 and 1.2 in every dimension with the sharp power $2/p_c$, and the five other formal targets (restriction, Kakeya maximal, Hausdorff dimension, Bochner–Riesz, Stein). Every step is written out or is a cited textbook theorem with hypotheses checked; the fibre Dirichlet existence theory is cited and also formally verified. A new script verifies the wave current's conservation directly, to $10^{-39}$.
- **Kakeya, the remaining claims** ([`10_GAPLESS_REMAINING.md`](reviews/claude/kakeya/10_GAPLESS_REMAINING.md)).
  - Complete proofs of cone and paraboloid restriction, Schrödinger local smoothing, the Nikodym maximal conjecture, and the weak-norm obstruction (App. C).
  - Apps. B, D and E are checked at the level of exponents and new identities.
  - Checked by [`verify_remaining_claims.py`](reviews/claude/kakeya/scripts/verify_remaining_claims.py).

## Distillation round 2 (recyclable toolkit)

See [`reviews/claude/SYNTHESIS.md`](reviews/claude/SYNTHESIS.md). New results proved and checked in this round:
- $p$-adic cross-equation truncation for every prime (generalizes the ETH paper's Lemma 5.2);
- exact weighted partition functions of bounded-degree factor graphs in $2^{o(n)}$;
- a proved cost floor showing the ETH method's exponent is optimal for symmetric designs;
- a rank-$r$ density-to-depth lemma;
- the **capacity-dimension principle**: the Kakeya proof's auxiliary time must be exactly 2-dimensional for logarithmic loss.

## Layout

```
ROADMAP.md                mathematical roadmap: results, dependencies, statuses, targets (start here)
roadmap.json              its machine-readable twin (validated by tools/check_roadmap.py)
tools/                    roadmap checker and its negative controls
source/eth/               pinned ETH sources (hashes, Lean statement)
source/kakeya/            pinned Kakeya sources (hashes, comparator statements, Lean scope doc)
reviews/claude/eth/       Claude's ETH review: 01 verdict … 07 toolkit, 08–11 gapless proofs, 12 F_q counting, FOR_GPT, scripts/, results/
reviews/claude/kakeya/    Claude's Kakeya review: 01 verdict … 07 toolkit, 08–10 gapless proofs, FOR_GPT, scripts/, results/, lean/
reviews/claude/SYNTHESIS.md  cross-paper synthesis: shared pattern, recyclable tools, new sharpenings
reviews/gpt/              GPT's reviews (to be added)
REVIEW_STANDARD.md        working standard for this project
```

The ETH scripts need Python 3.10+ and the standard library only. The Kakeya scripts also use `numpy`, `sympy` and `scipy`.
