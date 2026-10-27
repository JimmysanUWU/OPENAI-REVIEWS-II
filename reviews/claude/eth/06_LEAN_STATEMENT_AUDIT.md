# 6. Audit of the supplied Lean statement

File: [`../../../source/eth/SubexponentialSAT.statement.lean`](../../../source/eth/SubexponentialSAT.statement.lean)
(SHA-256 `4183b205…6685`, 221 lines, `import Mathlib`). The main theorem is `OAI.Algorithms.SubexponentialSAT.SqrtLogSpace.exact_count_and_sat`, and **its body is `sorry`**. This is a statement file of the kind used for comparator challenges, not a proof.

I read it by hand. I have not compiled it, because the session has no Mathlib toolchain and the file does not pin a Mathlib revision.

## What it states

For every $k\ge2$ there exist:
- a fixed `PrimitiveMachine` $\mathcal M$;
- distinct roles for a count tape and an answer tape;
- constants $C,A,S>0$ and $d,e>0$;

such that for every $n$ and every CNF $F$ over `Fin n` with all clauses of size $\le k$, the following hold.

1. Started on `encodeCNF F`, $\mathcal M$ reaches a halted configuration `final` in exactly `time` steps of the deterministic `locatedStep`.
2. The count tape is `Tape.mk₁ (countWord (count F))`. This is the little-endian binary count, with the head on its first bit and nothing else on the tape.
3. The answer tape is `Tape.mk₁ [some (decide (0 < count F))]`, and the halting state is `Sum.inl (decide (0 < count F))`.
4. $\mathrm{time}\le C\,(\lvert\mathrm{enc}\rvert+1)^d\cdot2^{A\,n(\log\log(n+16))^2/\sqrt{\log(n+16)}}$.
5. Every prefix of the run visits at most $S(\lvert\mathrm{enc}\rvert+1)^e$ tape cells in total.

The last two conjuncts about `countWord` (its value and its length $\le n+1$) are facts about `Nat.bits`. They do not constrain the machine.

## Fidelity

| Aspect | Check | Verdict |
|---|---|---|
| Machine model | Finite states and tapes over `Option Bool` with blank `none`. One scanned tape per state. Each transition halts with a Boolean or does one move/write and changes state. | A standard multitape TM, up to a constant-factor simulation. Deterministic. |
| Input | `encodeCNF`: unary $n$, then sorted clauses of sorted unary literal codes $2v+[\text{pos}]$, framed bitwise. | Injective given $n$, so the count is a function of the input. Length is $\Theta(n+\sum_C\sum_{\ell\in C}(2v_\ell+1))=\mathrm{poly}(n,m)$. Unary codes inflate only the polynomial factor. |
| Semantics | `Literal.eval a l = (a l.variable == l.positive)`; clause = ∃, CNF = ∀; `count` over `Fin n → Bool`. | Correct. The empty clause gives count $0$; the empty CNF gives $2^n$. |
| Width | `WidthAtMost k F`: every clause has card $\le k$ (as a `Finset` of literals). | Correct. Tautological clauses are allowed, which is harmless. |
| Quantifiers | The machine and constants are chosen after $k$ and before $n,F$. | Uniform algorithm for fixed $k$, as in Prop. 6.4. |
| Time | Exponent $A\cdot\rho(n)$ with $\rho$ evaluated at $n+16$; polynomial in encoding length. | Matches **Prop. 6.4**, not Thm 1.1's sharper bound. It is still $2^{o(n)}\mathrm{poly}(n,m)$, so it implies ETH false. |
| Space | `visitedSpace` sums `upper − lower + 1` over tapes; `locatedInit` counts the input cells; the monitor updates only on moves. | Correct polynomial-space formalization. |
| Output hygiene | Exact `Tape.mk₁` equality forces the machine to clean its output tapes. | Strong but fair. |
| Vacuity | If `countTape = answerTape`, conjunct 2 would force `countWord (count F) = [b]`, which fails for count $\ge2$. So a proof must use distinct tapes. No hypothesis is contradictory. | Not vacuous. |

**Conclusion.** The statement is a faithful and nontrivial formalization of Proposition 6.4. Taking $k=3$, it implies a deterministic $2^{o(n)}\mathrm{poly}(n+m)$-time, polynomial-space decision procedure for 3-SAT on a standard machine model. That is, it implies that ETH is false.

## What is needed to close the formal side

1. The full Lean project: proof files, `lakefile`, `lean-toolchain`, and the Mathlib revision. I need this from the user's link.
2. A build check that `exact_count_and_sat` has no `sorry` or extra axioms (`#print axioms` should show only `propext`, `Classical.choice`, `Quot.sound`).
3. A comparator run showing the proved theorem is definitionally this statement. The previous project used Nanoda and Lean replay for this.

Until then, the formal status is **statement audited, proof not seen**.
