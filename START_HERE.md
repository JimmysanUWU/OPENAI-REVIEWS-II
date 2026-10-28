# Start here

This page is for an AI agent (GPT, Claude, or this same Claude after its context was cleared) or a person picking up the project cold. Read it once, then follow §2. Everything it says is checkable in the repository.

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

**Where things stand (28 Oct 2026).**
- **Verification is complete for both papers.** Every claim has a prose proof in this repository; no node is `PARTIAL`. The last one, Kakeya App. B, is proved line by line in [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md).
- **Distillation.** The Kakeya frameworks are written: the engine theorem, the slot calculus, and the planar entry point. The open targets are listed in [`ROADMAP.md`](ROADMAP.md) §6.

## 2. Read in this order

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
     - [`09`](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) and [`10`](reviews/claude/kakeya/10_GAPLESS_REMAINING.md): the consequences.
   - **ETH:**
     - [`08`](reviews/claude/eth/08_GAPLESS_PROOF.md): Part I;
     - [`09`](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md): Part II;
     - [`10`](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md): Part III;
     - [`11`](reviews/claude/eth/11_GAPLESS_COROLLARIES.md): the corollaries;
     - [`12`](reviews/claude/eth/12_FP_COUNTING.md): the new $\mathbb F_q$ theorem;
     - [`13`](reviews/claude/eth/13_UNIFIED_DEFORMATION.md): one deformation theorem for Parts I–III, and both cost floors.
4. **Background, only if needed:**
   - the verdict, ledger, methods and "miracles" files (`01`–`04`) in each review folder;
   - the toolkits (`07`);
   - the Lean audits (`06`);
   - [`reviews/claude/SYNTHESIS.md`](reviews/claude/SYNTHESIS.md) (cross-paper).

## 3. Verify everything in three minutes

```
python3 tools/verify_all.py            # 2 roadmap checks + 23 check scripts, about 3 min with 4 jobs; prints ALL OK
python3 tools/verify_all.py --quick    # skips the four slowest scripts
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
reviews/claude/eth/     ETH review (01-13, FOR_GPT, scripts/, results/, lean/)
reviews/claude/kakeya/  Kakeya review (01-14, FOR_GPT, scripts/, results/, lean/)
reviews/gpt/            reserved for GPT's reviews (empty so far)
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
6. **The paraboloid** is reached only through the wave, by rescaling. A direct engine for it is open (`T.KAK.12`).
7. **App. B's threshold.** $B>\frac12+\frac1p$ is what the proof gives ($A+\frac12-\frac1p$ with $A=2/p$), not a proved optimum (`T.KAK.13`). Its uniform stationary phase holds only beyond $(\text{degree}+1)^2$ wavelengths. This is necessary, not an artefact (kakeya/14 Lemma 1.2 and its remark).

## 7. If you are continuing the work

- **Git.**
  - Branch: `claude/pensive-franklin-4oovjw`.
  - Draft PR: [JimmysanUWU/OPENAI-REVIEWS-II#1](https://github.com/JimmysanUWU/OPENAI-REVIEWS-II/pull/1).
  - Keep the PR description current.
- **Next targets** in the roadmap's order:
  - `T.ETH.4`: explicit thresholds;
  - the research targets (`T.KAK.4`–`T.KAK.9`, `T.KAK.12`, `T.KAK.13`, `T.ETH.1`, `T.ETH.5`, `T.X.1`).
- **Handoff notes** for the other model are in each project's `FOR_GPT.md`.
