# OPENAI-REVIEWS-II

Joint Claude–GPT review of OpenAI's mathematical releases. The workflow is the one used in [Claude-GPT-reviews](https://github.com/JimmysanUWU/Claude-GPT-reviews).

## User priorities (each target)

1. Does the paper answer the question?
2. Verify rigorously, and identify methods, leaps of insight and new structure; assess how radical they are.
3. Does the paper explain its ideas, or rely on unexplained "miracles"?
4. Distill the mathematics into better structure.

## Targets

| Target | Claim | Claude verdict | GPT |
|---|---|---|---|
| [**ETH**](reviews/claude/eth/01_VERDICT.md): *The Exponential Time Hypothesis Is False* (OpenAI, Oct 2026) | #$k$-SAT in $2^{O_k(n\log\log n/\sqrt{\log n})}$ time, polynomial space | ETH refutation correct (97%). The minimal chain was re-derived, and the algorithm was run in exact arithmetic. | pending |
| [**Kakeya / local smoothing**](reviews/claude/kakeya/01_VERDICT.md): *Local smoothing for the Euclidean wave equation* (OpenAI, 25 Oct 2026) | Sogge's local smoothing conjecture in all dimensions, which implies Bochner–Riesz, restriction, Stein square function, Kakeya maximal and Kakeya set conjectures | Answers the question. The core identities were verified. **OpenAI's formal proof builds locally** (2,261 modules, 0 errors, 0 sorry; only standard axioms; elaborated statements byte-identical to the challenge) ([06](reviews/claude/kakeya/06_LEAN_AUDIT.md)). A new method: positive currents plus quaternionic Monge–Ampère energies. | pending |

## Distillation round 2 (recyclable toolkit)

See [`reviews/claude/SYNTHESIS.md`](reviews/claude/SYNTHESIS.md). New results proved and checked in this round:
- $p$-adic cross-equation truncation for every prime (generalizes the ETH paper's Lemma 5.2);
- exact weighted partition functions of bounded-degree factor graphs in $2^{o(n)}$;
- a proved cost floor showing the ETH method's exponent is optimal for symmetric designs;
- a rank-$r$ density-to-depth lemma;
- the **capacity-dimension principle**: the Kakeya proof's auxiliary time must be exactly 2-dimensional for logarithmic loss.

## Layout

```
source/eth/               pinned ETH sources (hashes, Lean statement)
source/kakeya/            pinned Kakeya sources (hashes, comparator statements, Lean scope doc)
reviews/claude/eth/       Claude's ETH review: 01 verdict … 07 toolkit, FOR_GPT, scripts/, results/
reviews/claude/kakeya/    Claude's Kakeya review: 01 verdict … 07 toolkit, FOR_GPT, scripts/, results/, lean/
reviews/claude/SYNTHESIS.md  cross-paper synthesis: shared pattern, recyclable tools, new sharpenings
reviews/gpt/              GPT's reviews (to be added)
REVIEW_STANDARD.md        working standard for this project
```

The ETH scripts need Python 3.10+ and the standard library only. The Kakeya scripts also use `numpy`, `sympy` and `scipy`.
