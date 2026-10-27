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
| [**Kakeya / local smoothing**](reviews/claude/kakeya/01_VERDICT.md): *Local smoothing for the Euclidean wave equation* (OpenAI, 25 Oct 2026) | Sogge's local smoothing conjecture in all dimensions, which implies Bochner–Riesz, restriction, Stein square function, Kakeya maximal and Kakeya set conjectures | Answers the question. The core identities were verified, and the formal proof was audited and built ([06](reviews/claude/kakeya/06_LEAN_AUDIT.md)). A new method: positive currents plus quaternionic Monge–Ampère energies. | pending |

## Layout

```
source/eth/               pinned ETH sources (hashes, Lean statement)
source/kakeya/            pinned Kakeya sources (hashes, comparator statements, Lean scope doc)
reviews/claude/eth/       Claude's ETH review: 01 verdict … 06 Lean, FOR_GPT, scripts/, results/
reviews/claude/kakeya/    Claude's Kakeya review: 01 verdict … 06 Lean audit, FOR_GPT, scripts/, results/
reviews/gpt/              GPT's reviews (to be added)
REVIEW_STANDARD.md        working standard for this project
```

The ETH scripts need Python 3.10+ and the standard library only. The Kakeya scripts also use `numpy`, `sympy` and `scipy`.
