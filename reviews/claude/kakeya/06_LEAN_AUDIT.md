# 6. Audit of OpenAI's Lean formalization (local smoothing ⇒ Kakeya)

**Status: complete.** OpenAI's formal proof builds locally against the pinned Mathlib with no errors and no `sorry`. Its six comparator theorems depend only on `propext`, `Classical.choice` and `Quot.sound`. Their fully elaborated statements are byte-identical to the challenge file's. See §6.5–6.6.

## 6.1 What is formalized

`openai/math` `lean/docs/079.md` (copied to [`../../../source/kakeya/lean-docs-079.md`](../../../source/kakeya/lean-docs-079.md)) lists six comparator targets. They are in [`../../../source/kakeya/LocalSmoothing.challenge.lean`](../../../source/kakeya/LocalSmoothing.challenge.lean), and each challenge body is `sorry`. The solution module is `OAI.Analysis.Fourier.LocalSmoothing.Comparator`, and the permitted axioms are `propext`, `Quot.sound` and `Classical.choice`.

| Comparator theorem | Paper result | Fidelity |
|---|---|---|
| `euclidean_local_smoothing` | Thm 1.1: for $n\ge2$, $2\le p<\infty$, and $\varepsilon>0$ (or $p=2,\varepsilon=0$): $\|U(T)f\|_{L^p(\mathbb R^n\times[1,2])}\le C\|J^{s_+(p)+\varepsilon}f\|_p$ for Schwartz $f$ | Faithful. The half-wave and $J^s$ are defined by explicit Fourier integrals with angular frequency, as in the paper. |
| `bochner_riesz` | Thm 10.1, (10.1): uniform in $R$, for $1\le p\le\infty$ and $\delta>\max(n|1/p-1/2|-1/2,0)$ | Faithful, including $p=\infty$. |
| `spherical_restriction` | Cor 9.3 (9.3): $E:L^a(S^{n-1})\to L^b(\mathbb R^n)$ for $a>1$, $b>2n/(n-1)$, $b>\frac{n+1}{n-1}a'$ | Faithful (with $a'=1$ at $a=\infty$). |
| `kakeya_maximal` | Cor 9.1: $\|K_\delta f\|_{L^n(S^{n-1})}\le C_\varepsilon\delta^{-\varepsilon}\|f\|_{L^n}$, $0<\delta\le1$ | Faithful. The tube is the open $\delta$-neighbourhood of a unit segment; the maximal function is the supremum over centres. (Weaker than the paper's log bound, as expected.) |
| `kakeya` | Cor 9.1: every $F\subseteq\mathbb R^n$ containing a unit segment in every direction has $\dim_HF=n$ | Faithful, and **no measurability assumed**. |
| `stein_square` | Thm 11.1: $\alpha>\max(1/2,n(1/2-1/p))$, $2\le p<\infty$, with the literal radial derivative | Faithful. |

**Junk values.** Bochner integrals of non-integrable functions are 0 in Mathlib. Every integrand here is integrable for the stated inputs: Schwartz data, an $L^a$ function on a finite-measure sphere, and Bochner–Riesz symbols that are bounded with compact support. So no statement is vacuous or weakened by junk values. `deriv` is used in the Stein target, where the radial derivative exists for Schwartz data and $\alpha>1/2$.

## 6.2 Statement comparison

Stripping comments and proofs, the definitions and statements in the solution module are **textually identical** to the challenge file. The script is in this review's history; the only differences are proof bodies.

The solution file writes `open _root_.MeasureTheory`, because it lives inside the `OAI` namespace. This resolves to the same Mathlib namespace as the challenge's `open MeasureTheory`.

A textual match is not the same as an elaborated-term match: an imported instance could, in principle, change what `volume` or `inner` elaborate to. §6.6 therefore compares the **fully elaborated** terms.

## 6.3 Import closure and escape-hatch scan

Computed by [`lean/import_closure.py`](lean/import_closure.py) on `openai/math` at commit `7ef5c0cc` (26 Oct 2026). Data: [`lean/closure_files.txt`](lean/closure_files.txt), [`lean/closure_scan.json`](lean/closure_scan.json).

| Quantity | Value |
|---|---|
| OAI modules in the closure of `…LocalSmoothing.Comparator` | **2,261** |
| Lines | **190,467** |
| External dependencies | **Mathlib only** (none of the repo's other git dependencies) |
| `sorry` / `admit` / `axiom` declarations | **0 / 0 / 0** |
| `native_decide`, `implemented_by`, `@[extern]`, `unsafe`, `opaque`, `ofReduceBool`, `debug.skipKernelTC` | **0** |

Distribution by area (largest first):

| Area | Modules |
|---|---|
| PDE/Dirichlet | 398 |
| PDE/Schauder | 274 |
| PDE/Concentration | 243 |
| PDE/Elliptic | 223 |
| Fourier/WavePackets | 137 |
| PDE/Weak | 134 |
| Fourier/SquareFunction | 132 |
| Fourier/Restriction | 124 |
| Fourier/BochnerRiesz | 91 |
| PDE/Green | 84 |
| Calculus/Current | 58 |
| Calculus/Quaternion | 16 |
| MeasureTheory/Kakeya | ~20 |

**The Dirichlet and Schauder theory of App. A is formalized**, not assumed.

**Top of the formal chain** (`OAI/Analysis/Fourier/LocalSmoothing/Main.lean`):
```
logarithmic_refinement   := logarithmicRefinement_of_packet_solvability positive_packet_solvability
euclidean_local_smoothing := euclideanLocalSmoothing_of_logarithmicRefinement logarithmic_refinement
```
`positive_packet_solvability` (existence of the fiber Monge–Ampère/Dirichlet solutions) is itself a proved theorem built from the Dirichlet family, not a hypothesis.

## 6.4 Formal versus paper statement of the key estimate

`LogarithmicRefinement` (`LocalSmoothing/Definitions.lean`) is Theorem 1.2 with $C(1+\log(1/h))^{2n+5}$ in place of $C[\log(2/h)]^{2/p_c}$. The doc's 3D tube-volume exponent $66=2n(2n+5)$ at $n=3$ is consistent with this.

All six targets need only *some* fixed power, so they are unaffected. The sharp $2/p_c$, Theorem B.1, and the sharpness propositions of §8 are **not** among the comparator targets. The repository contains a `Sharpness/` directory (40 files), which is not in the comparator's closure and was not audited.

## 6.5 Local build (done)

Reproduction:
- [`lean/lakefile.lean`](lean/lakefile.lean) and [`lean/lean-toolchain`](lean/lean-toolchain): Lean `v4.34.1`, Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612`, with `autoImplicit false` as in OpenAI's lakefile.
- Copy the closure files from `openai/math/lean/` into `OAI/`.
- `lake exe cache get`, then `lake build OAI.Analysis.Fourier.LocalSmoothing.Comparator`.

Practical note: Lake writes an ~8 MB `*.setup.json` per module under `.lake/build/ir`, about 18 GB for this closure. Deleting them is safe; `lake build --no-build` still reports every target up to date.

| Run | Modules built | Outcome |
|---|---|---|
| 1 (background) | 2,024 | stopped at the 2-hour session limit; 0 errors |
| 2 | 169 | stopped at the chunk time limit; 0 errors |
| 3 | 66 | stopped by **disk full** (an environment failure, not a Lean error) |
| 4 | 2 | `Build completed successfully (11184 jobs)` |
| **Total** | **2,261** | equals the closure size; **0 errors, 0 `sorry` warnings** |

Details: [`lean/checks/build_summary.txt`](lean/checks/build_summary.txt).

## 6.6 Axioms and statement identity

[`lean/checks/Solution.lean`](lean/checks/Solution.lean) prints every comparator theorem and definition with `pp.all`, then runs `#print axioms`. [`lean/checks/Challenge.lean`](lean/checks/Challenge.lean) prints the same in the challenge environment (`import Mathlib` plus the challenge file).

**Axioms.** All six theorems (`euclidean_local_smoothing`, `bochner_riesz`, `spherical_restriction`, `kakeya_maximal`, `kakeya`, `stein_square`) depend only on `[propext, Classical.choice, Quot.sound]`. These are exactly the permitted axioms ([`lean/checks/print_axioms.txt`](lean/checks/print_axioms.txt)).

**Statement identity.** The `pp.all` output covers the 16 definitions (`Space`, the Fourier transforms, the half-wave, the Bessel potential, the observation measure, sphere measure, extension, Bochner–Riesz, tubes, the Kakeya maximal function, `IsKakeyaSet`, the Stein square function) and the 6 theorem types. It is **byte-identical** in the two environments: 411,445 characters, SHA-256 `db5adcf8…dce9dce` ([`lean/checks/pp_all_comparison.txt`](lean/checks/pp_all_comparison.txt)). Every instance, coercion and auxiliary constant elaborates to the same term. This is the property OpenAI's Comparator tool certifies; I checked it independently.

## 6.7 Conclusion

Modulo the soundness of Lean 4.34.1's kernel and of Mathlib at the pinned revision, the following are **formally proved**, for every dimension $n\ge2$:
- Sogge's local smoothing conjecture (with $\varepsilon$ loss);
- uniform Bochner–Riesz bounds in the full conjectured range;
- spherical restriction in the strict conjectured range;
- the Kakeya maximal conjecture with $\delta^{-\varepsilon}$ loss;
- **full Hausdorff dimension of every Kakeya set**;
- Stein's square function conjecture.

Each is in exactly the form stated in [6.1](#61-what-is-formalized).

What the formal proof does **not** certify:
- the sharp logarithmic power $2/p_c$ of Theorem 1.2 (formalized with $2n+5$);
- Theorem B.1;
- the sharpness propositions of §8;
- that the prose is the same proof. The formal propagation step differs: it uses a *corrected superharmonic* energy (`Concentration/CorrectedEnergy.lean`, correction $16-z_i^2-z_j^2$) and a cubic-logarithm cylinder lemma (`Concentration/Cylinder/Estimate.lean`) in place of Lemma 6.1's supersolution with potential.
