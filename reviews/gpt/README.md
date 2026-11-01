# GPT research continuation — 1 November 2026

**Round 10 (abstraction search, latest):** [ETH 22](eth/22_CIRCUITS_AND_ERASURE.md) strips the erasure mechanism to minimal inconsistent circuits. Sign-compatible splitting of an optimal residual gives an exact splitting of its loss. If contradictions have at most c equations, inverse-weighted retention and refitting leave at least a theta^(c-1) fraction of the original loss, for every power s>1. The exponent and c-wise independence threshold are sharp. Exact locality becomes a circuit-hitting problem; no sharper asymptotic directional bound is claimed. The [route log](ABSTRACTION_SEARCH_ROUND10.md) preserves the tested failures, nonlinear-model counterexample and remaining algebraic-degree candidate.

Targeted checks pass; see the [round-10 validation report](check-logs/round10_validation.json). The roadmap has 106 nodes and 23 targets, with all 24 validator controls passing. The unchanged full suite was not repeated.

**Round 9:** [Route log](ABSTRACTION_SEARCH_ROUND9.md) retains the tested failures and the lead developed in Kakeya 23: global curvature plus critical state dimension forces a sphere bundle, then two spectral levels and a classical spin-factor algebra. The global dispersion is necessarily a drifted anisotropic wave cone. This uses explicitly cited topology and Clifford theory; local nonquadratic sheets remain possible.


**Round 8 (phone continuation):** [Product-bias energy](eth/21_PRODUCT_BIAS_ENERGY.md) solves **T.ETH.8** by diagonalizing squared edge readouts. It gives an explicit optimal family and a sharp degree-dependent comparison: for the actual swap bias, the relative difference from the uniform optimum is at most D/r^2. These are prose proofs with independent exact rational QPs; no new Lean build or SAT speedup is claimed.

The task has moved from auditing the original results to extracting and extending their mathematics. At the user's direction, **Claude's reported Lean runs are accepted as a starting point**. No new Lean build is claimed. The earlier [assessment](ASSESSMENT.md) is retained as a dated audit, including its narrower formal-scope observations and its measurement-transfer witness.

**Round 8, analytic continuation:** [Curvature and quadric fibres](kakeya/22_CURVATURE_AND_QUADRIC_FIBRES.md) proves that EVERY rigid curved multiple sheet has infinite central time density for the unweighted full-fibre single lift. The second eigenvalue variation gives an orthogonal design and computes the conditional measure directly. This includes the octonionic case left open in Claude 21. T.KAK.7 remains open for different currents and quotient reductions.

## Earlier round (reviewed and merged)

| Result | Status | Read |
|---|---|---|
| Arbitrary, including ungraded, swap designs have a directional defect of measure at least `exp(-O(R²/Q))/poly(Q)` | Complete prose proof; answers `T.ETH.6`; exact finite checks | [ETH 15](eth/15_UNGRADED_SPLINE_BOUND.md) |
| Parity completion in a cube spline algebra; unequal jet-order extension | Complete reusable algebraic derivation | [ETH 15, §§1–2 and 5](eth/15_UNGRADED_SPLINE_BOUND.md) |
| Aggregate measurement transfer from complete copies, and the ungraded enumeration floor under that condition | Complete prose proof with explicit filling and cost-model hypotheses | [ETH 16](eth/16_MEASUREMENT_TRANSFER_AND_FLOOR.md) |
| Multiplication-rank certificates, joint trace completion, odd-residue unit defects, and sharp uniform defect energy | Complete reusable proofs; classical orthogonal-polynomial ingredient derived explicitly | [ETH 17](eth/17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md) |
| Full complex-fibre integration diverges for the realified Weyl single lift at a plane wave's polarization centre | Complete explicit obstruction; leaves other reductions open | [Kakeya 18](kakeya/18_NONCOMPACT_FIBRE_OBSTRUCTION.md) |
| The nearby Weyl density has an inverse-square transverse singularity and is not locally integrable | Complete proof with an explicit logarithmic divergence coefficient | [Kakeya 19](kakeya/19_TRANSVERSE_SINGULARITY.md) |

All five research notes include definitions, proofs, applications, limitations, and deterministic checks. They are not formalized in Lean. “New” means new to this repository; worldwide novelty has not been established. The [research ledger](RESEARCH_LEDGER.md) records the abstract routes, small exhaustive experiments, failed shortcuts, and unfinished interfaces.

## The main mechanism, in plain language

An ETH swap design assigns a small polynomial to each corner of a cube. Across an edge, locality asks two polynomial values to cancel. Multiplying by parity changes cancellation into agreement. This turns the design into a graph spline: a polynomial-valued function whose endpoints agree at the number attached to an edge.

On the intact cube, small-degree splines can only produce a limited set of values at zero. A spline that produces parity lets multiplication generate every missing character in about half the total degree. That creates a dimension deficit. Deleting one agreement condition can repair at most one dimension, so many edges must fail. Random restriction transfers the conclusion to the biased spectrum used by the counting algorithm.

This combines an established graph-spline framework with Smolensky's parity-complement idea. The [primary references](eth/15_UNGRADED_SPLINE_BOUND.md#6-attribution-checks-and-scope) are credited in the proof. The concrete contribution is the parity-trace completion theorem and its application to the previously open ungraded-design question.

The next abstraction is multiplication of traces. Several low-degree traces sharing the same good edges can multiply to create additional independent characters; a weighted Hamming-ball union quantifies the shared constraint cost. Separately, the exact uniform quadratic defect minimum is `rho/(1-rho)`, where `rho = product_{j=0}^D (Q-j)/(Q+j+1)`. It has a sharp degree-of-order-`sqrt(Q)` transition. These give distinct tools for exact vanishing and approximate locality, with the arithmetic and measure assumptions kept explicit.

## Why this direction was chosen

The roadmaps, handoffs, the previous GPT assessment, and the principal deformation/current frameworks identify several viable directions. The pinned upstream materials remain OpenAI's October 25 preprints, with the [October 26 announcement](https://openai.com/index/building-scientific-understanding-with-academic-communities/) as provenance. This continuation does not claim a fresh line-by-line audit of every source file.

- **Ungraded designs:** the obstruction was a sharply defined finite-dimensional problem. The spline interpretation removes the exact assumption that stopped Claude's argument and yields a standalone theorem.
- **Measurement transfer:** the previous audit found a real missing interface. Complete copies and support injection supply a clean sufficient repair. The condition is kept visible instead of upgrading every historical floor claim indiscriminately.
- **Hermitian engines:** the roadmap's next analytic estimates presupposed a finite fibre integral. The Weyl calculation isolates an explicit failure before undertaking those estimates.
- **One logarithm:** still attractive, but the available scalar propagation bound is sharp for its stated class. No stronger propagation theorem was obtained here.
- **The broad square-root law:** still needs a precise worst-case class. Already-local systems are trivial counterexamples to an unrestricted universal formulation; the swap theorem does not settle the broad question.

## What this does and does not accomplish

The ETH result closes an abstract design question and strengthens a conditional barrier for a specified enumeration architecture. **It does not improve the SAT running time.** The useful design lesson is that removing grading inside the swap group algebra cannot evade this defect scale. Measurement compression and genuinely different deformations remain possible routes.

The Kakeya results rule out one natural integral, not Hermitian engines in general. The realified Weyl pencil has the wave's dispersion, so it supplies a controlled test case rather than a new dispersive estimate. The second calculation shows that ignoring the exceptional axis cannot help: nearby positive mass is already infinite. A quotient or altered current must account for the noncompact gauge direction and prove conservation after reduction.

## Reproduction

From the repository root, with Python 3.10 or later:

```text
python reviews/gpt/eth/scripts/verify_spline_bound.py
python reviews/gpt/eth/scripts/verify_complete_copies.py
python reviews/gpt/eth/scripts/verify_defect_capacity.py
python reviews/gpt/eth/scripts/verify_circuit_erasure.py
python reviews/gpt/kakeya/scripts/verify_weyl_fibre.py
python reviews/gpt/kakeya/scripts/verify_weyl_singularity.py
python tools/check_roadmap.py
python tools/test_check_roadmap.py
```

The five new check scripts use only the standard library. Their result files are committed beside the notes. The full `tools/verify_all.py` now discovers both Claude's and GPT's scripts and normalizes script paths for Windows. The old Kakeya scripts additionally require NumPy, SciPy, SymPy, and mpmath. The separately runnable `explore_directional_extremizers.py` records exhaustive tiny-instance experiments; it is not part of the verification suite.

The earlier 26-script rerun is recorded in [check-logs](check-logs/). The new spline tests include 42 feasible design cases, 30 genuinely ungraded, along with trace-space dimensions, biased coupling, weighted jets, and negative controls. Finite tests support the derivations; they are not substitutes for the proofs.

## Next useful work

1. Independently inspect the short parity-completion proof, especially the quotient-to-trace dimension step. Then formalize this small algebraic result if desired.
2. Determine whether the factor `Q` in the maximum-direction uniform bound is necessary (`T.ETH.7`). The exponential scale is settled; sharp constants and extremizers are not.
   ETH 22 supplies an exact circuit-hitting formulation. A direction-sensitive packing certificate would be more useful than brute-force circuit enumeration.
3. Bound the actual image size of correlated measurement maps. The present cost floor deliberately addresses the Cartesian enumeration expression.
4. Apply joint trace completion to a construction with shared exceptional edges (`T.ETH.9`); the biased energy problem (`T.ETH.8`) is now solved in ETH 21.
5. Construct a different finite conserved reduction (`T.KAK.7`); Kakeya 22 now rules out the unweighted full-fibre single lift for every rigid multiple sheet.

The compact continuation instructions are in [FOR_CLAUDE.md](FOR_CLAUDE.md).


## Round 8 validation

All 38 verification scripts passed (21 ETH, 17 Kakeya), together with the 101-node roadmap validator and its 24 negative controls. See [validation report](check-logs/round8_validation.json) and [analytic rerun](check-logs/round8_kakeya.txt). The original Lean builds are accepted as directed by the user and were not repeated.
