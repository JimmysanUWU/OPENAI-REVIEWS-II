# Start here

**Claude's review of GPT rounds 8–10 (1 November, latest):**
- **Verdict.** ETH 21–22 and Kakeya 22–23 are verified, with no gap found, and merged ([`reviews/claude/GPT_ROUND2_REVIEW.md`](reviews/claude/GPT_ROUND2_REVIEW.md)).
- **Independent checks.** [`eth/scripts/verify_gpt_round2.py`](reviews/claude/eth/scripts/verify_gpt_round2.py) and [`kakeya/scripts/verify_gpt_round2.py`](reviews/claude/kakeya/scripts/verify_gpt_round2.py).
- **Relation to Claude's notes.** Kakeya 22 generalizes Claude's kakeya/20 Theorem 20.2, with no gauge group and octonions included.

**Round 10 (latest):** [Minimal contradictions](reviews/gpt/eth/22_CIRCUITS_AND_ERASURE.md) explain the degree exponent without polynomials or squared loss. A sign-compatible splitting of the residual splits the full loss into circuit losses. Contradictions of size at most c give a sharp theta^(c-1) bound after inverse-weighted deletion and refitting; c-wise independence suffices and is sharp. Exact locality is a circuit-hitting problem, with no improved asymptotic bound yet. [Search log](reviews/gpt/ABSTRACTION_SEARCH_ROUND10.md).

**Round 9:** [Route log](reviews/gpt/ABSTRACTION_SEARCH_ROUND9.md) retains the tested failures and the lead developed in Kakeya 23: global curvature plus critical state dimension forces a sphere bundle, then two spectral levels and a classical spin-factor algebra. The global dispersion is necessarily a drifted anisotropic wave cone. This uses explicitly cited topology and Clifford theory; local nonquadratic sheets remain possible.


**Round 8 (phone continuation):** [Product-bias energy](reviews/gpt/eth/21_PRODUCT_BIAS_ENERGY.md) solves **T.ETH.8** by diagonalizing squared edge readouts. It gives an explicit optimal family and a sharp degree-dependent comparison: for the actual swap bias, the relative difference from the uniform optimum is at most D/r^2. These are prose proofs with independent exact rational QPs; no new Lean build or SAT speedup is claimed.

This page is for an AI agent (GPT, Claude, or this same Claude after its context was cleared) or a person picking up the project cold. Read it once, then follow §2. Everything it says is checkable in the repository.

**Current continuation (1 November):** Read [GPT / Astra research notes](reviews/gpt/README.md) first for the new work. The user accepts Claude's Lean runs and has shifted the task to mathematical extensions. `T.ETH.6` is now proved via cube splines; the enumeration floor has an explicit complete-copy hypothesis. ETH 17 adds joint trace completion, residual obstructions, and sharp uniform defect energy. The naive Weyl full-fibre push-forward is not locally integrable on the real target slice (Kakeya 19), sharpening `T.KAK.7`. New results are prose proofs with finite checks, not new Lean builds.

**Claude's review of that round (1 November, later):**
- **Verdict.** Every GPT proof checks out, and all scripts pass; see [`reviews/claude/GPT_ROUND_REVIEW.md`](reviews/claude/GPT_ROUND_REVIEW.md), item by item.
- **The measurement-transfer gap is closed** under the algorithm's own balanced schedules ([eth/18](reviews/claude/eth/18_AGGREGATE_TRANSFER.md)):
  - at most three labels are under-filled;
  - the floors of 07 E4 and 14 hold as stated;
  - GPT's ungraded floor no longer needs complete copies.
- **The Weyl divergence is universal** ([kakeya/20](reviews/claude/kakeya/20_GAUGE_ORBIT_OBSTRUCTION.md)):
  - every complex and quaternionic Hermitian single lift has infinite time density at every polarization centre, because the complexified gauge group preserves the flux map but is invisible to the damping;
  - the complex-target profile also closes a step in GPT 19.
  - `T.KAK.7` now asks for a current on the gauge quotient.

**Round 7 (1 November, latest):**
- **The ETH reading guide** [eth/19](reviews/claude/eth/19_HOW_TO_READ_THE_PROOF.md), with a worked example ([`verify_worked_example.py`](reviews/claude/eth/scripts/verify_worked_example.py)).
- **[kakeya/21](reviews/claude/kakeya/21_GAUGE_QUOTIENT.md).** The flux map is a linear projection of a classical rank-one cone: Veronese for $\mathbb R$, Segre for $\mathbb C$, Plücker for $\mathbb H$. This explains 17's counts $n=2m'-1$ and $4m''-3$ as cone dimensions, and puts tangency and divergence into one factorization.
- **[eth/20](reviews/claude/eth/20_DIRECTIONAL_DEFECTS.md).** Subcube averaging interpolates between the face lemma and GPT's bound (progress on `T.ETH.7`). Exact small-cube data support the conjecture that there is no directional loss.

**Round 8, analytic continuation:** [Curvature and quadric fibres](reviews/gpt/kakeya/22_CURVATURE_AND_QUADRIC_FIBRES.md) proves that EVERY rigid curved multiple sheet has infinite central time density for the unweighted full-fibre single lift. The second eigenvalue variation gives an orthogonal design and computes the conditional measure directly. This includes the octonionic case left open in Claude 21. T.KAK.7 remains open for different currents and quotient reductions.

## 1. The project in one minute

The user asked for a joint Claude–GPT review of two OpenAI mathematical releases, with four priorities for each:
1. Does it answer the question?
2. Verify rigorously; name the methods and the leaps.
3. Explain any "miracles".
4. **Distill** the mathematics into better structure, gaplessly.

| Paper | Claim | Verdict (subjective confidence) | Formal status |
|---|---|---|---|
| *The Exponential Time Hypothesis Is False* (OpenAI, Oct 2026) | #k-SAT in $2^{O_k(n\log\log n/\sqrt{\log n})}$ time, polynomial space | **Correct.** ETH is false (≈99%). The headline exponent via App. D is ≈95%. | OpenAI's Lean proof of the weaker $2^{O(n(\log\log n)^2/\sqrt{\log n})}$ bound builds here, with only the standard axioms |
| *Local smoothing for the Euclidean wave equation* (OpenAI, 25 Oct 2026) | Sogge's local smoothing conjecture in all dimensions, which implies Bochner–Riesz, restriction, Kakeya maximal and sets, Stein, Nikodym | **Correct** (≈99% for the six formal targets; ≈90% for the sharp power $2/p_c$, which is prose-only) | OpenAI's Lean proof of all six targets builds here, with only the standard axioms |

**Earlier Claude review (29 Oct 2026), with the qualification above.**
- **Claude reports completed verification of the original papers.** The present research accepts the recorded Lean runs; added framework claims are tracked individually, including the measurement-transfer qualification in GPT 16. The last one, Kakeya App. B, is proved line by line in [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md).
- **Distillation.** The Kakeya frameworks are written: the engine theorem, the slot calculus, and the planar entry point. For ETH, one deformation theorem covers Parts I–III ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md)), and the design floor ([eth/14](reviews/claude/eth/14_DESIGN_FLOOR.md)) shows that no graded swap design beats Part I's exponent. The open targets are listed in [`ROADMAP.md`](ROADMAP.md) §6.
- **Latest round (Kakeya).**
  - A reading guide at four levels: [kakeya/15](reviews/claude/kakeya/15_HOW_TO_READ_THE_PROOF.md).
  - The engine in light-cone coordinates is its simplest instance, and its one-mass sector is a direct engine for Schrödinger (`T.KAK.12` done): [kakeya/16](reviews/claude/kakeya/16_NULL_FRAME_AND_SCHRODINGER.md).
  - A linear-algebra test decides which symmetric systems can carry the engine's current. It corrects a claim in 11 §7.3: [kakeya/17](reviews/claude/kakeya/17_STATE_DIMENSION_RIGIDITY.md).

## 2. Read in this order

0. **New to either proof?** Start with a reading guide.
   - **Kakeya:** [kakeya/15](reviews/claude/kakeya/15_HOW_TO_READ_THE_PROOF.md).
   - **ETH:** [eth/19](reviews/claude/eth/19_HOW_TO_READ_THE_PROOF.md), with a worked example that runs the algorithm on a 4-variable formula and prints every intermediate number.

   Each has four levels, from no formulas to the Lean statements, with a ledger of every number in the proof and checkpoints tied to the scripts.
1. **[`PROOF_MAP.md`](PROOF_MAP.md).** Every result as a step in a dependency-ordered path, each with a statement and a proof sketch.
   - Start with the path *Shortest complete proof* (Kakeya, 6 steps, linear PDE only) or *ETH is false* (7 steps).
   - About 20 minutes per project.
2. **[`ROADMAP.md`](ROADMAP.md).** The dependency graphs, the per-node tables, the frameworks extracted so far, the open targets, and the contribution protocol (§7).
3. **The full proofs**, reached from each node's `Full proof` link. The core files:
   - **Kakeya:**
     - [`08`](reviews/claude/kakeya/08_GAPLESS_PROOF.md): the paper's proof, gapless;
     - [`11`](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md): the engine theorem, pencils, and the wave as the arrow pencil;
     - [`12`](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md): the calculus;
     - [`13`](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md): the planar proof;
     - [`14`](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md): the angular endpoint (App. B);
     - [`16`](reviews/claude/kakeya/16_NULL_FRAME_AND_SCHRODINGER.md): the null frame, the Husimi current, and a direct engine for Schrödinger;
     - [`17`](reviews/claude/kakeya/17_STATE_DIMENSION_RIGIDITY.md): which symmetric systems carry the engine (the tangency test);
     - GPT's [`18`](reviews/gpt/kakeya/18_NONCOMPACT_FIBRE_OBSTRUCTION.md)–[`19`](reviews/gpt/kakeya/19_TRANSVERSE_SINGULARITY.md) and [`20`](reviews/claude/kakeya/20_GAUGE_ORBIT_OBSTRUCTION.md): tangency is not enough; the fibre integral diverges for every Hermitian pencil;
     - [`09`](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) and [`10`](reviews/claude/kakeya/10_GAPLESS_REMAINING.md): the consequences.
   - **ETH:**
     - [`08`](reviews/claude/eth/08_GAPLESS_PROOF.md): Part I;
     - [`09`](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md): Part II;
     - [`10`](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md): Part III;
     - [`11`](reviews/claude/eth/11_GAPLESS_COROLLARIES.md): the corollaries;
     - [`12`](reviews/claude/eth/12_FP_COUNTING.md): the new $\mathbb F_q$ theorem;
     - [`13`](reviews/claude/eth/13_UNIFIED_DEFORMATION.md): one deformation theorem for Parts I–III, and both cost floors;
     - [`14`](reviews/claude/eth/14_DESIGN_FLOOR.md): the design floor: every graded swap design pays $(\log\log n)^2/\sqrt{\log n}$, because of Smolensky's bound for parity;
     - GPT's [`15`](reviews/gpt/eth/15_UNGRADED_SPLINE_BOUND.md)–[`17`](reviews/gpt/eth/17_TRACE_CAPACITY_AND_APPROXIMATE_LOCALITY.md): parity completion in cube splines (every design, graded or not), complete-copy transfer, trace capacity;
     - [`18`](reviews/claude/eth/18_AGGREGATE_TRANSFER.md): measurement transfer under the algorithm's own schedules; the floors without extra hypotheses.
4. **Background, only if needed:**
   - the verdict, ledger, methods and "miracles" files (`01`–`04`) in each review folder;
   - the toolkits (`07`);
   - the Lean audits (`06`);
   - [`reviews/claude/SYNTHESIS.md`](reviews/claude/SYNTHESIS.md) (cross-paper).

## 3. Verify everything in three minutes

```
python3 tools/verify_all.py            # 2 roadmap checks + 42 check scripts (Claude's and GPT's), about 3 min with 4 jobs; prints ALL OK
python3 tools/verify_all.py --quick    # skips the five slowest scripts
python3 tools/check_roadmap.py --path P.KAK.PLANAR   # read a path in the terminal
python3 tools/check_roadmap.py --node KAK.T12        # a node, its sketch, and everything beneath it
```

**Requirements.** Python ≥ 3.10. The ETH scripts use only the standard library; the Kakeya scripts need `numpy`, `sympy`, `scipy` and `mpmath`.

**Side effect.** The scripts rewrite their JSON in `results/`. Only timing fields such as `"seconds"` change between runs; restore them with `git checkout reviews/claude/*/results` before committing unless results really changed.

**The Lean builds.** These are not rerun by these commands (each takes hours). Their logs and audits are in `reviews/claude/*/06_*` and `reviews/claude/*/lean/`.

## 4. Repository map

```
START_HERE.md      this page
PROOF_MAP.md       the theorem DAG as readable paths (generated from roadmap.json)
ROADMAP.md         graphs, node tables, frameworks, targets, protocol (partly generated)
roadmap.json       the source of truth: nodes, sketches, dependencies, paths, glossary, targets
tools/             check_roadmap.py (validator and renderer), test_check_roadmap.py, verify_all.py
REVIEW_STANDARD.md the working standard
source/            pinned sources: hashes, Lean statements, where the papers live
reviews/claude/GPT_ROUND_REVIEW.md  Claude's verdicts on GPT's round (ETH 15-17, Kakeya 18-19, the assessment)
reviews/claude/eth/     ETH review (01-14, 18-20, FOR_GPT, scripts/, results/, lean/; 19 is the reading guide)
reviews/claude/kakeya/  Kakeya review (01-17, 20-21, FOR_GPT, scripts/, results/, lean/; 15 is the reading guide)
reviews/gpt/            research continuation: ETH 15-17, Kakeya 18-19, ledger, checks, earlier audit
```

**The papers themselves are not committed.**
- The Kakeya TeX and PDF are public in `openai/math` (see [`source/kakeya/SOURCES.md`](source/kakeya/SOURCES.md)).
- The ETH PDF was supplied by the user; its hash is in [`source/eth/SOURCES.md`](source/eth/SOURCES.md).

## 5. Rules that keep the repository trustworthy

- **Verdicts.**
  - Give the narrowest justified verdict: complete derivation, accepted source, exact computation, or formal build.
  - Confidences are subjective and labelled as such.
- **Node discipline.**
  - Never rename or renumber a node ID. Add new nodes instead.
  - `FORMAL` marks only statements whose Lean form was audited. Internal lemmas get `lean_area`, never `FORMAL`.
- **Workflow.**
  1. Edit `roadmap.json`.
  2. Run `python3 tools/check_roadmap.py --write`.
  3. Run `python3 tools/check_roadmap.py`. It must print `OK`.
  
  Never hand-edit text between `GENERATED` markers.
- **New mathematics.**
  - Write it in a numbered file in the project folder, in definition → theorem → proof order.
  - Write statements as `**Lemma 3.2 (name).**` so the checker can resolve labels.
  - Every new identity gets a deterministic check script with a negative control that must fail. The script prints `ALL OK` and writes JSON to `results/`.
- **Ownership.**
  - GPT writes under `reviews/gpt/<project>/`.
  - Agreement or dispute goes in a node's `reviews` list (ROADMAP §7.4). Neither model deletes the other's entries.

## 6. Facts that are easy to get wrong

1. **The ETH formal statement** is Prop. 6.4's $\rho(n)=n(\log\log n)^2/\sqrt{\log n}$ bound, not the headline Theorem 1.1. The headline exponent rests on App. D, which is proved in prose ([eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md)) and not formalized.
2. **ETH Cor. 1.4** ($E^{NP}$ circuits of size $n\sqrt{\log n}$) needs App. D. The formal bound alone gives $n\sqrt{\log n}/(\log\log n)^2$.
3. **ETH App. D's floor** $u\ge\lambda^{4/5}$ holds only from $\log_2\lambda\approx500$. The proof needs only two absorption inequalities, which hold from ≈26. This is not an error.
4. **The Kakeya formal proof** proves Theorem 1.2 with the power $2n+5$, not the sharp $2/p_c$. The sharp power is prose ([kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8). The formal Kakeya maximal statement is the $\delta^{-\varepsilon}$ form.
5. **The wave is the arrow pencil.** Its passive coordinate in 11 is spatial, while 08 (the paper) uses a null direction. Both satisfy the engine axioms, so they give two proofs of the cap estimate (`KAK.T12` has an alternative proof).
6. **The paraboloid has its own engine.** Schrödinger is the one-mass sector of the wave in light-cone coordinates ([kakeya/16](reviews/claude/kakeya/16_NULL_FRAME_AND_SCHRODINGER.md) Theorem 16.8), which proves the paper's (12.2) directly. (12.1) and (12.3) still go through the sphere ([kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §4).
7. **Ungraded defect bounds are now proved.** [GPT 15](reviews/gpt/eth/15_UNGRADED_SPLINE_BOUND.md) answers `T.ETH.6` without grading. [GPT 16](reviews/gpt/eth/16_MEASUREMENT_TRANSFER_AND_FLOOR.md) proves the enumeration floor with complete-copy transfer, and [eth/18](reviews/claude/eth/18_AGGREGATE_TRANSFER.md) removes that hypothesis: the algorithm's own balanced schedules suffice. This is a bound on that cost expression, not every faster implementation. The original justification of E4's measurement cost ("all $h$ copies") was wrong for partly filled copies; it is repaired, not deleted.
8. **App. B's threshold.** $B>\frac12+\frac1p$ is what the proof gives ($A+\frac12-\frac1p$ with $A=2/p$), not a proved optimum (`T.KAK.13`). Its uniform stationary phase holds only beyond $(\text{degree}+1)^2$ wavelengths. This is necessary, not an artefact (kakeya/14 Lemma 1.2 and its remark).
9. **Larger systems.** For a pencil with $N>n$ states, the lifted current is tangent to a passive coordinate only if $N'=n+m-1$, where $N'$ is the dimension of the span of the extremal eigenvectors and $m$ their multiplicity ([kakeya/17](reviews/claude/kakeya/17_STATE_DIMENSION_RIGIDITY.md) Theorem 17.4). An earlier version of 11 §7.3 said the push-forward stays tangent; it is corrected. Tangency does not imply integrability: [GPT 18](reviews/gpt/kakeya/18_NONCOMPACT_FIBRE_OBSTRUCTION.md) proves divergence of the naive full-fibre Weyl current, and [kakeya/20](reviews/claude/kakeya/20_GAUGE_ORBIT_OBSTRUCTION.md) proves it for every complex and quaternionic Hermitian pencil.

## 7. If you are continuing the work

- **Git.**
  - Working branch: `claude/pensive-franklin-4oovjw`, with draft PR [JimmysanUWU/OPENAI-REVIEWS-II#1](https://github.com/JimmysanUWU/OPENAI-REVIEWS-II/pull/1).
  - GPT's round 1 (`gpt/astra-spline-extensions`, [JimmysanUWU/OPENAI-REVIEWS-II#2](https://github.com/JimmysanUWU/OPENAI-REVIEWS-II/pull/2)) is merged into it after review.
  - New GPT work: branch from `claude/pensive-franklin-4oovjw` and open a PR into it.
  - Keep the PR description current.
- **Next targets** in the roadmap's order:
  - `T.ETH.4`: explicit thresholds;
  - the research targets (`T.KAK.4`, `T.KAK.6`–`T.KAK.9`, `T.KAK.13`, `T.ETH.5`, `T.ETH.7`–`T.ETH.9`, `T.X.1`).
- **Current handoffs:**
  - Claude → GPT: "Update (1 Nov)" at the top of [eth/FOR_GPT.md](reviews/claude/eth/FOR_GPT.md) and [kakeya/FOR_GPT.md](reviews/claude/kakeya/FOR_GPT.md).
  - GPT → Claude: [reviews/gpt/FOR_CLAUDE.md](reviews/gpt/FOR_CLAUDE.md), answered in [GPT_ROUND_REVIEW.md](reviews/claude/GPT_ROUND_REVIEW.md).

