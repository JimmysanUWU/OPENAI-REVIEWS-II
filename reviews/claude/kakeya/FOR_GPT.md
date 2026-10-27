# Handoff to GPT: the Kakeya / local smoothing paper (Claude, 27 October 2026)

The user's four priorities are the same as for ETH. The source is OpenAI's *Local smoothing for the Euclidean wave equation* (25 Oct 2026). It is byte-identical to `openai/math` `preprints/Local-smoothing-…-October-25-2026/paper.pdf`, and the TeX is in that directory.

## Update (27 Oct, later): gapless proofs

- [08](08_GAPLESS_PROOF.md) proves Theorems 1.1 and 1.2 (sharp power $2/p_c$) in every dimension. Every step is either written out or a cited textbook theorem with hypotheses checked.
  - The one substantial cited block is the existence and regularity half of the fibre Dirichlet theorem (Schauder, Evans–Krylov/Krylov). The comparison principle is proved, and the whole theorem is also in the Lean proof.
  - New check: [`scripts/verify_wave_current_target.py`](scripts/verify_wave_current_target.py) verifies conservation of the wave current $\chi\chi^*$ directly in the target coordinates, to $10^{-39}$, independently of the push-forward derivation.
- [09](09_GAPLESS_CONSEQUENCES.md) derives the other five formal targets: restriction, Kakeya maximal, Hausdorff dimension, Bochner–Riesz, Stein. Tao's $\varepsilon$-removal and Fefferman–Stein are cited.
- **The most useful independent check now** is a read of 08 §5 (Theorem 5.2, especially the 4-slot boundary $C^2$ step) and §7 (constants).

## Where I am

- **Verdict** ([01](01_VERDICT.md)). It answers the question. Theorem 1.2 (critical annular local smoothing at $p_c$ with a log loss) implies local smoothing, Bochner–Riesz, the Stein square function, restriction (strict range), Kakeya maximal, Kakeya sets of full dimension, and Nikodym, in all dimensions.
- **Verified** ([02](02_VERIFICATION_LEDGER.md)). §2 (current), §3 (quaternionic calculus), §4 (energy balance, Lemma 4.3), §6 (propagation, sharpness example), and the exponent bookkeeping in §7–§12. These are symbolic or numeric with negative controls ([`scripts/`](scripts/)). §5 logic is re-derived and its constants accepted. App. A is read but not re-derived.
- **Lean** ([06](06_LEAN_AUDIT.md)):
  - **Built locally** against pinned Mathlib: 2,261 modules, 0 errors, 0 `sorry`.
  - All six comparator theorems use only `propext`, `Classical.choice` and `Quot.sound`.
  - Their `pp.all` statements are byte-identical to the challenge's.
  - The formal log power is $2n+5$, not $2/p_c$; that is enough for every headline result.
  - The formal propagation step uses a corrected superharmonic energy and cubic logarithms instead of Lemma 6.1.
- **Distillation** ([05](05_DISTILLATION.md)): four modular theorems.
  - I. Positive wave current.
  - II. Determinant energy balance under abstract hypotheses (A), (B), (W), (D).
  - III. Density to depth at $p_c$.
  - IV. Two-circle propagation.
  - The $n=2$ case is linear and makes a five-page entry point.

## Suggested division (to avoid duplicate work)

| Target | Status here | Suggestion for GPT |
|---|---|---|
| **App. A**: 3-slot quaternionic MA Dirichlet problem on periodic fibers, especially the boundary $C^2$ estimate with the added fourth slot and smooth parameter dependence for moving domains | read only | **Highest-value independent prose check.** |
| **§5 constants**: kernel bounds (5.5)–(5.6), the inversion cost $e^{C\gamma^2\sigma}$ versus the wall gain $e^{-c_0\sigma}$, ball inclusions | logic only | Recompute. |
| Lemma 4.1 (fiber geometry bounds (4.3)–(4.4)) | structure only | Recompute. |
| App. B (angular endpoint), App. C (weak-norm obstruction), App. E (pencils) | not checked | Lower priority: not used by the formal targets. |
| Prior OpenAI Sept preprints [62]–[64] and the 3D companion | not reviewed | Compare methods and check whether the October paper supersedes them. |

## Open questions I would value your view on

1. Is Theorem II (energy balance) best seen as a *linear* analogue of Serre's compensated integrability? Is there a cleaner abstract inequality behind it?
2. Can the holomorphic flux map (2.9) be described intrinsically, e.g. on the complexified characteristic variety?
3. Remark 6.2 shows two logarithms are sharp for scalar supersolutions. Is one logarithm reachable with the actual energy profile? The known lower bound is $1/p_c$.

## Also relevant: the ETH paper's formal proof exists

`openai/math` `lean/ComparatorChallenges/SubexponentialSAT.lean` has the same SHA-256 as the user's statement file. Its solution module is `OAI.Algorithms.SubexponentialSAT.Main`. See [`../eth/06_LEAN_STATEMENT_AUDIT.md`](../eth/06_LEAN_STATEMENT_AUDIT.md) for my follow-up.
