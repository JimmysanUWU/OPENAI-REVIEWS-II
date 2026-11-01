# 1. Verdict: "Local smoothing for the Euclidean wave equation" (the Kakeya paper)

Claude review, 27 October 2026. Source: `kakeyagpt.pdf`, 59 pp., OpenAI, 25 Oct 2026. It is byte-identical to `openai/math` `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf` (SHA-256 in [`../../../source/kakeya/SOURCES.md`](../../../source/kakeya/SOURCES.md)). The formal development is described in `openai/math` `lean/docs/079.md`.

## What the paper claims

One estimate at the critical exponent $p_c=\frac{2n}{n-1}$, for every $n\ge2$ (Theorem 1.2):
$$\|U(T)P_Nf\|_{L^{p_c}(\mathbb R^n\times[1,2])}\le C\,[\log(2+N)]^{2/p_c}\|f\|_{p_c},\qquad U(T)=e^{-iT\sqrt{-\Delta}}.$$
Here $P_N$ is a smooth cutoff to frequencies $|\xi|\approx N$.

Interpolation and standard implications give the following.

| Conjecture | Status before this paper | Claimed here |
|---|---|---|
| Sogge's local smoothing conjecture (1991) | $n=2$ (Guth–Wang–Zhang 2020); partial results for $n\ge3$ | all $n\ge2$, every $\varepsilon>0$ loss (Thm 1.1) |
| Bochner–Riesz conjecture | $n=2$ (Carleson–Sjölin 1972); open for $n\ge3$ | all $n$, all $p$, uniform in the radius, plus maximal bounds for $p\ge2$ (Thm 10.1) |
| Stein's square function conjecture | $n=2$ (Carbery); partial for $n\ge3$ | all $n$ (Thm 11.1) |
| Fourier restriction conjecture: sphere, paraboloid, cone | $n=2$ (Fefferman, Zygmund); open for $n\ge3$ | the interior of the conjectured range, all $n$ (Cor 9.3, 12.1) |
| Kakeya maximal conjecture | $n=2$ (Córdoba); $n=3$ (OpenAI, Sept 2026) | all $n$, with $\log^{4/p_c}$ loss (Cor 9.1) |
| Kakeya set conjecture (Hausdorff dimension $n$) | $n=2$ (Davies 1971); $n=3$ (Wang–Zahl 2025); $n=4$ (OpenAI, Sept 2026) | all $n$ (Cor 9.1) |
| Nikodym maximal conjecture | open for $n\ge3$ | all $n$ (9.7) |

If correct, this settles the central conjectures of Euclidean harmonic analysis in one stroke. By significance it is far beyond the ETH paper and the 7/8 result. See [03](03_METHODS_AND_NOVELTY.md).

## Short answer to priority 1

**Yes, it answers the question.** Theorem 1.2 at $p_c$ with any fixed power of the logarithm implies:
- local smoothing with $\varepsilon$ loss;
- the Kakeya maximal bound with $\delta^{-\varepsilon}$ loss;
- full Hausdorff dimension of Kakeya sets.

These implications are classical: Sogge and Mockenhaupt–Seeger–Sogge for Bochner–Riesz; Bourgain and Beckner–Carbery–Semmes–Soria for Kakeya; Tao's $\varepsilon$-removal for global restriction. The paper reproves the routes it needs, and I checked them (§9–§11).

## Status of verification

| Layer | What I did | Status |
|---|---|---|
| **The formal proof (Lean 4 + Mathlib, `openai/math`)** | Computed the full import closure: 2,261 OAI modules, 190,467 lines, depending only on Mathlib. Scanned it for escape hatches. **Built it locally** against pinned Mathlib `d13f23b`. | **Builds with 0 errors and 0 `sorry`.** All six theorems use only `propext`, `Classical.choice`, `Quot.sound`. Fully elaborated (`pp.all`) statements are **byte-identical** to the challenge's. See [06](06_LEAN_AUDIT.md). |
| Fidelity of the six formal statements to the paper | Read line by line; compared the elaborated terms against the challenge. | **Faithful**, and identical at the elaborated-term level. |
| §2 the conserved current | Symbolic and numerical checks of (2.1)–(2.4), Lemma 2.1, (2.9)–(2.13). Negative control: wrong sign of imaginary shift. | **Verified** ([`results/wave_current.json`](results/wave_current.json)) |
| §3 quaternionic Hessian calculus | Lemma 3.1 including (3.4); Lemma 3.2 (3.6)–(3.8); Lemma 3.3 cofactor divergence in 2, 3 and 4 slots. Negative control: a non-Hessian field. | **Verified** ([`results/quaternion_calculus.json`](results/quaternion_calculus.json)) |
| §4 energy balance and comparison | Re-derived Prop 4.2 (cancellation of the full determinant; vanishing normal cofactor) and Lemma 4.3 (also numerically). | **Verified** |
| §5 fixed-fiber estimate | Re-derived the logic: upper comparison $v_*=C\sigma^m rP_*$, interior quadratic comparison, exponent bookkeeping $(1+1/d)\cdot2=p_c$. Kernel constants followed, not recomputed. | Logic verified; constants accepted |
| §6 propagation | Re-derived Lemma 6.1 and Prop 6.3. Numerically checked the Gronwall radial bound. Verified exactly the sharpness example of Remark 6.2. | **Verified** ([`results/propagation_and_exponents.json`](results/propagation_and_exponents.json)) |
| §7–§12 transfers and applications | Re-derived the interpolation exponents (7.3), the focusing exponent of Prop 8.1, and Lemma 9.2/9.4 bookkeeping. Checked that the Stein and Bochner–Riesz thresholds match. | Verified at the level of exponents and logic |
| App. A (Dirichlet problem for 3-slot quaternionic Monge–Ampère on periodic fibers) | Read; the continuity method plus CNS-type boundary estimates is plausible. Formalized (`OAI/Analysis/PDE/Dirichlet`, 398 modules). | Not independently re-derived |
| App. B–E | **28 Oct, [10](10_GAPLESS_REMAINING.md):** App. C (Prop. 8.2) is re-derived in full. B, D and E are checked at the level of every exponent, every new algebraic identity (quaternionic completed square, rank-one Moore expansion, (E.2), Example E.2) and the structure of each argument. | C verified; E proved in [11](11_FRAMEWORK_POSITIVITY_ENGINE.md); D proved in [13](13_PLANAR_ENTRY_POINT.md); B proved line by line in [14](14_GAPLESS_ANGULAR_ENDPOINT.md) |

## Confidence (subjective, after verification)

- **Six headline consequences**, as formalized: **≈99%**. These are local smoothing with ε loss, Bochner–Riesz, the restriction interior range, Kakeya maximal with δ^−ε, Kakeya sets of full dimension, and the Stein square function. The formal proof builds and its axioms and statements check out ([06](06_LEAN_AUDIT.md)). The residual risk is a soundness bug in Lean or Mathlib, which I regard as very small for proofs of this kind, or a community-level objection to how a classical notion is encoded. I examined the encodings and found none.
  - My own prose verification of the core alone would put it at ≈85–90%. The formal proof is what lifts this.
- **Sharp logarithmic power $2/p_c$ in Theorem 1.2**: ≈90%. The formalization proves the weaker power $2n+5$, which suffices for every headline result; the sharp power is prose-only.
- **App. C** (weak-norm obstruction): ≈98%, re-derived in full.
- **App. E** (symmetric pencils): ≈96%. It is now proved in full as Corollary 11.3 of [11](11_FRAMEWORK_POSITIVITY_ENGINE.md), an abstract engine theorem plus a verification that pencils satisfy its axioms. The pencil current's conservation is checked directly for five pencils.
- **App. D** (planar case by a linear problem): ≈97%. It is proved in full in [13](13_PLANAR_ENTRY_POINT.md), using only linear PDE.
- **App. B** (angular endpoint): ≈95%, given Theorem 1.2. It is proved line by line in [14](14_GAPLESS_ANGULAR_ENDPOINT.md). Its uniform-in-degree stationary phase is replaced by an explicit Bessel bound ([14](14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 1.2). The residual risk is in the cited vector-valued maximal and spherical-multiplier theorems.
- **Cone and paraboloid restriction, Schrödinger local smoothing, and Nikodym**: ≈97%, given Theorem 1.2. Re-derived in [10](10_GAPLESS_REMAINING.md) §2–§4, with cited textbook theorems.

## Relation to OpenAI's earlier releases

The references cite three OpenAI preprints from September 2026:
- Kakeya maximal in $\mathbb R^3$ [63];
- full Hausdorff dimension in $\mathbb R^4$ [62];
- elliptic capacity propagation and restriction to the sphere [64].

There is also a companion, *Critical local smoothing for the three-dimensional wave equation* (Sept 24). The October paper subsumes them. This review does not assess the September papers separately.
