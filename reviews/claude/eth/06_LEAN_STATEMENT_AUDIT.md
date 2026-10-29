# 6. Audit of the Lean statement and of OpenAI's formal proof

**Status (27 Oct, update): the formal proof exists, and it builds and checks.** OpenAI's `openai/math` proves this exact statement file: `lean/ComparatorChallenges/SubexponentialSAT.lean` has the same SHA-256. Its 3,091-module proof builds against the pinned Mathlib. The theorem uses only `propext`, `Classical.choice` and `Quot.sound`. Its fully elaborated statement, and all 42 definitions it uses, are byte-identical to the challenge file's. See §6.A at the end.

File: [`../../../source/eth/SubexponentialSAT.statement.lean`](../../../source/eth/SubexponentialSAT.statement.lean)
(SHA-256 `4183b205…6685`, 221 lines, `import Mathlib`). The main theorem is `OAI.Algorithms.SubexponentialSAT.SqrtLogSpace.exact_count_and_sat`, and **its body is `sorry`**. This is a statement file of the kind used for comparator challenges, not a proof.

I read it by hand (first pass). Since then I have compiled it, and OpenAI's proof of it, against the Mathlib revision pinned in `openai/math` (§6.A).

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

## 6.A Formal proof audit (done)

| Item | Result |
|---|---|
| Solution module | `OAI.Algorithms.SubexponentialSAT.Main` (`openai/math` commit `7ef5c0cc`, `lean/ComparatorChallenges/SubexponentialSAT.json`); permitted axioms `propext`, `Quot.sound`, `Classical.choice` |
| Import closure | **3,091** OAI modules, **547,460** lines; depends only on Lean core and Mathlib ([`lean/closure_scan.json`](lean/closure_scan.json), [`lean/closure_files.txt`](lean/closure_files.txt)) |
| Escape-hatch scan | 0 `sorry`, 0 `admit`, 0 `axiom`, 0 `native_decide`, 0 `implemented_by`/`extern`/`unsafe`. There are 3 `opaque` declarations (`SqrtLog/LeafProgramPreparedRecovery{Query,Witness}.lean`), each with an explicit value checked by the kernel; they hide large data, not assumptions. |
| Build | Lean `v4.34.1`, Mathlib `d13f23b`. All 3,091 modules built. The only failures were two SIGKILLs from the container's memory limit (exit 137, confirmed by `dmesg` oom-kill) while four large modules compiled in parallel; both modules then built alone. Final `lake build --no-build`: *All targets up-to-date (12014 jobs)*. 0 `sorry` warnings ([`lean/checks/build_summary.txt`](lean/checks/build_summary.txt)). |
| Axioms | `exact_count_and_sat` depends on `[propext, Classical.choice, Quot.sound]` ([`lean/checks/print_axioms.txt`](lean/checks/print_axioms.txt)). Eight modules carry in-repo `#print axioms` audits; all report the same three. |
| Statement identity | `pp.all` output of the theorem type and of all 42 definitions/structures/inductives it uses (`bitVal`, `toNat`, `Literal`, `Clause`, `CNF`, `eval`s, `WidthAtMost`, `count`, the whole encoding, `Trace`, `PrimitiveMachine`, `Cfg`, `step`, `init`, `runtimeScale`, `runtimeBound`, `countWord`, `LocatedCfg`, `locatedInit`, `headDisplacement`, `monitor`, `locatedStep`, `visitedSpace`). It is **byte-identical** in the challenge and solution environments: 51,981 characters, SHA-256 `805c7f15…fe31f3` ([`lean/checks/pp_all_comparison.txt`](lean/checks/pp_all_comparison.txt)). |
| What the proof constructs | An explicit 1,220-tape machine, compiled from a `BitProgram` by a verified compiler (`GuardedCompiler.primitive program`), with count tape 256 and answer tape 1216. The time and visited-space bounds are proved for every prefix of the run. |

**Conclusion.** Modulo the soundness of Lean's kernel and Mathlib, the following is **formally proved**. For every fixed $k\ge2$, a deterministic multitape Turing machine counts the satisfying assignments of any $k$-CNF on $n$ variables exactly. It runs in time $(L+1)^{O_k(1)}\,2^{O_k(n(\log\log n)^2/\sqrt{\log n})}$ and visits $(L+1)^{O_k(1)}$ cells, where $L$ is the input length.

Taking $k=3$, **ETH is false**, and so are SETH and #ETH. This is Prop 6.4 of the paper (the $\rho(n)$ bound), not the sharper Theorem 1.1, which rests on App. D. The formal route goes through §5–§6 (sparse $\mathbb F_2$ systems, Tanner/Cayley-graph tests: `SqrtLog/SparseCodeTanner.lean`, `SqrtLog/CayleyGraph.lean`) and the swap design.

Reproduce with [`lean/build_sequential.sh`](lean/build_sequential.sh) and the Kakeya review's [`lakefile`](../kakeya/lean/lakefile.lean). Copy the ETH closure files into `OAI/`, and build memory-heavy modules one at a time.
