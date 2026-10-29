# Review standard

This standard is carried over from the previous project's `MATHEMATICAL_REVIEW_STANDARD.md` and adapted to the user's priorities for this project (27 October 2026).

**Order of work.** Work through the user's priorities in this order:
1. Settle whether the result answers its question.
2. Verify rigorously, and identify the new methods and leaps.
3. Check whether the ideas are explained or rely on unexplained coincidences.
4. Distill the mathematics into reusable definitions, lemmas and short proofs.

Distillation is the long-term deliverable. Verification of load-bearing steps comes first.

**Verdicts.** For each claim, say which of these it is:
- a complete independent derivation;
- a theorem accepted from a checked primary source;
- a finite or exact computation;
- a reported formal build.

State the narrowest justified verdict. Confidence percentages are subjective and must be labelled as such. Small exact computations confirm identities and catch errors. They do not establish asymptotic claims.

**Exposition.** Prefer definition → theorem → proof → application. Remove procedural narrative. Explain a failure mode only when it isolates an essential hypothesis.

**Publication.** Publish polished results incrementally. Keep a compact current handoff (`FOR_GPT.md`, `FOR_CLAUDE.md`) rather than a stale history.

**Reproducibility.** Scripts use the standard library only, are deterministic (fixed seeds), and write their outputs to `results/`.
