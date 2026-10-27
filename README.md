# OPENAI-REVIEWS-II

Joint Claude–GPT review of OpenAI's mathematical releases. The workflow is the one used in [Claude-GPT-reviews](https://github.com/JimmysanUWU/Claude-GPT-reviews).

## Current target: *The Exponential Time Hypothesis Is False* (OpenAI, Oct 2026)

| | |
|---|---|
| Claim | A deterministic algorithm counts the satisfying assignments of any $k$-CNF in $(n+m)^{O_k(1)}2^{O_k(n\log\log n/\sqrt{\log n})}$ time and polynomial space. This refutes ETH, SETH, #ETH, NSETH and Gap-ETH. |
| Claude verdict | **The ETH refutation is correct** (confidence 97%). The minimal chain is self-contained and was re-derived completely, and the core algebra and full algorithm were checked in exact arithmetic. Appendix D (sharper exponent) is partly checked. The Lean statement is faithful, but its proof has not been seen. |
| GPT verdict | *pending* |

## User priorities

1. Does the paper answer the question?
2. Verify rigorously, and identify methods, leaps of insight and new structure; assess how radical they are.
3. Does the paper explain its ideas, or rely on unexplained "miracles"?
4. Distill the mathematics into better structure.

## Layout

```
source/                 pinned sources (hashes; Lean statement)
reviews/claude/         Claude's review
  01_VERDICT.md                    priority 1
  02_VERIFICATION_LEDGER.md        priority 2: lemma-by-lemma status
  03_METHODS_AND_NOVELTY.md        priority 2: methods, leaps, radicalness
  04_MIRACLES_EXPLAINED.md         priority 3: spectral explanation, optimality, barrier
  05_DISTILLATION.md               priority 4: blueprint (Theorems A–E)
  06_LEAN_STATEMENT_AUDIT.md       formal statement audit
  FOR_GPT.md                       handoff and suggested division of work
  scripts/  results/               exact-arithmetic checks and their outputs
reviews/gpt/            GPT's review (to be added)
REVIEW_STANDARD.md      working standard for this project
```

Reproduce every computation with Python 3.10+ and the standard library only (see [`reviews/claude/scripts/README.md`](reviews/claude/scripts/README.md)).
