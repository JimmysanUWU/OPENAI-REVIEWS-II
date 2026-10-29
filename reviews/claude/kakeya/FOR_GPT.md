# Handoff to GPT: the Kakeya / local smoothing paper (Claude, 27 October 2026)

The user's four priorities are the same as for ETH. The source is OpenAI's *Local smoothing for the Euclidean wave equation* (25 Oct 2026). It is byte-identical to `openai/math` `preprints/Local-smoothing-…-October-25-2026/paper.pdf`, and the TeX is in that directory.


**Cold start:** read [`START_HERE.md`](../../../START_HERE.md) and [`PROOF_MAP.md`](../../../PROOF_MAP.md) first. Run `python3 tools/verify_all.py` from the repository root to re-check everything.
## Update (29 Oct, newest): reading guide, the null frame, and which systems carry the engine

- **[15](15_HOW_TO_READ_THE_PROOF.md), a reading guide.** The proof at four levels, a ledger of every number in it, a dictionary to standard mathematics, notation traps, and checkpoints tied to scripts. Useful for checking that our account of the mechanism matches yours.
- **[16](16_NULL_FRAME_AND_SCHRODINGER.md), `T.KAK.12` done.**
  - The paper's current kept in light-cone coordinates satisfies the engine axioms with a polynomial chart (Jacobian $aY^d$) and a one-multiplier time row (Theorem 16.3).
  - Its time component is the wave continued along an imaginary light ray; for Schrödinger, the Husimi density (Prop. 16.5, Cor. 16.6).
  - Schrödinger is the one-mass sector, so the paper's (12.2) follows directly (Theorem 16.8).
  - Worth a second look: the Galilean covering and the two-lattice-mass step in Theorem 16.8(1).
- **[17](17_STATE_DIMENSION_RIGIDITY.md), `T.KAK.7` sharpened and corrected.**
  - A single lift is tangent iff $N'=n+m-1$ (Theorem 17.4). Our earlier 11 §7.3 claimed the $N>n$ push-forward stays tangent; that was wrong and is corrected there.
  - Please check the "only if" argument: constant rank, then the fibre identity $z^{\mathsf T}Sz=2(\varphi Q_0-\xi\cdot Q')$, then the fibre is an $(m-1)$-sphere.
- Scripts: [`verify_null_frame.py`](scripts/verify_null_frame.py) and [`verify_state_dimension.py`](scripts/verify_state_dimension.py).

## Update (28 Oct): App. B line by line, so no Kakeya node is `PARTIAL`

- **[14](14_GAPLESS_ANGULAR_ENDPOINT.md)** proves the angular endpoint (Theorem B.1) step by step. It adds four sub-nodes, `KAK.ANG.HANKEL`, `KAK.ANG.BALL`, `KAK.ANG.ENV` and `KAK.ANG.BAND`, on reading path `P.KAK.ANGULAR`.
- **The paper's one sketched step is replaced.** That step is the stationary phase uniform in the angular degree. In its place, Funk–Hecke plus Hankel's integral give $|(x\partial_x)^k\varepsilon_\mu(x)|\le e^{3/2}(2\mu+2k+1)^{\max(k-1,0)}(\mu+1)^2/(2x)$ for $x\ge(\mu+1)^2$, with explicit constants. A check shows the threshold $(\mu+1)^2$ cannot be dropped.
- **§7 explains the exponents.** Observation on a ball costs the logarithm of its radius in wavelengths (shell transplant). $d-\alpha p=-1$ makes the strong norm pay one logarithm and the weak norm none. The threshold is $B>A+\frac12-\frac1p$, where $A$ is the ball cost and $\frac12-\frac1p$ is the $\ell^2$-versus-$\ell^p$ price.
- **Open, `T.KAK.13`:** is $B>\frac12+\frac1p$ sharp?
- Script: [`scripts/verify_angular_endpoint.py`](scripts/verify_angular_endpoint.py).

## Update (28 Oct, earlier): roadmap and distillation

- The repository root now has a **mathematical roadmap** ([`../../../ROADMAP.md`](../../../ROADMAP.md)), with machine-readable twin `roadmap.json`. Each result has a stable ID (`KAK.E.BAL`, `KAK.PROP`, …), a precise statement, a status, its dependencies and its checks. Start there, and record agreement or disputes in a node's `reviews` list (ROADMAP §7.4).
- **Distillation, `T.KAK.1` done: [11](11_FRAMEWORK_POSITIVITY_ENGINE.md), the positivity engine as a theorem.**
  - Theorem 11.1 states axioms on a fibre chart (G1)–(G2) and a current (A1)–(A4) that never mention the equation, and proves the critical cap estimate from them alone.
  - Theorem 11.2 shows that every symmetric pencil with a simple, extremal, transversally curved eigenvalue satisfies the axioms. This completes App. E (Corollary 11.3).
  - Corollary 11.4: **the wave is the arrow pencil** $\begin{pmatrix}\xi_n&\xi'^{\mathsf T}\\\xi'&-\xi_nI\end{pmatrix}$. Its polarization is stereographic projection, and the paper's null-coordinate system is this pencil after a frame change and a rescaling. Theorem 1.2's cap estimate thus gets a second derivation.
  - §7 is a hypothesis-usage matrix with what breaks without each hypothesis. For example, a middle eigenvalue still gives a conserved current, but its Gaussian window grows (checked).
  - Script: [`scripts/verify_positivity_engine.py`](scripts/verify_positivity_engine.py). The pencil current's conservation is checked directly in target coordinates for five pencils, including the non-quadric cubic of Ex. E.2.
- **Distillation, `T.KAK.3` done: [12](12_QUATERNIONIC_SLOT_CALCULUS.md), the slot calculus as a standalone toolkit.**
  - The Moore cofactor of a frame Hessian is divergence-free for **arbitrary quaternion frames** (Theorem 12.1, checked symbolically, including dependent 5-vector frames).
  - Every frame calculus is a restriction of the 4-slot one (Lemma 12.2).
  - The Monge–Ampère operator is the gradient of a strictly concave energy (Theorem 12.3).
  - Ellipticity holds iff the frames are independent, so at most 4 coordinates per slot (Proposition 12.4).
  - The 3-slot Hessian is the real Hessian plus the complex Hessian plus an $\mathrm i$-coupling (Proposition 12.5).
- **Distillation, `T.KAK.2` done: [13](13_PLANAR_ENTRY_POINT.md), the planar entry point.**
  - It is a complete proof of sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ with $\log^{1/2}$ loss, using **only linear PDE**: a Dirichlet Laplacian with potential on a disk, Fourier series in $X$, and the maximum principle.
  - The chart is explicit: $W$ is the stereographic coordinate of $(P-\Phi e_2)/\tau$ on a complexified unit circle, and $\Phi=P_2-\sqrt{\tau^2-P_1^2}$.
  - The energy step is a completed-square inequality.
  - This is App. D, now with every step written out. It is my recommended entry point to the paper.
- **Most useful independent check now:** read 11 §1 (axioms) and §4 (Proposition 4.1) against 08 §7, and confirm that nothing wave-specific is used. Read 13 §3 (the planar energy) as a short, fully linear instance.

## Update (27 Oct, later): gapless proofs

- [08](08_GAPLESS_PROOF.md) proves Theorems 1.1 and 1.2 (sharp power $2/p_c$) in every dimension. Every step is either written out or a cited textbook theorem with hypotheses checked.
  - The one substantial cited block is the existence and regularity half of the fibre Dirichlet theorem (Schauder, Evans–Krylov/Krylov). The comparison principle is proved, and the whole theorem is also in the Lean proof.
  - New check: [`scripts/verify_wave_current_target.py`](scripts/verify_wave_current_target.py) verifies conservation of the wave current $\chi\chi^*$ directly in the target coordinates, to $10^{-39}$, independently of the push-forward derivation.
- [09](09_GAPLESS_CONSEQUENCES.md) derives the other five formal targets: restriction, Kakeya maximal, Hausdorff dimension, Bochner–Riesz, Stein. Tao's $\varepsilon$-removal and Fefferman–Stein are cited.
- **The most useful independent check now** is a read of 08 §5 (Theorem 5.2, especially the 4-slot boundary $C^2$ step) and §7 (constants).

## Update (28 Oct): the remaining claims

[10](10_GAPLESS_REMAINING.md) covers everything outside the formal targets.
- **Complete proofs:**
  - cone restriction (Nicola's slicing);
  - the Nikodym maximal conjecture (the projective map, checked symbolically);
  - paraboloid restriction and Schrödinger local smoothing (parabolic rescaling, Fresnel identity; Kim and Lee–Rogers–Seeger cited);
  - the weak-norm obstruction Prop. 8.2 / App. C (the tube arrangement, checked exactly).
- **Apps. B, D and E** (refinements) are checked at the level of every exponent and every new identity:
  - the quaternionic completed square and the rank-one Moore expansion in D;
  - formula (E.2) and Example E.2 in E.
  
  Their analytic estimates re-run §4–§7 and are not rewritten line by line.
- Script: [`scripts/verify_remaining_claims.py`](scripts/verify_remaining_claims.py). No errors found.

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
| App. B (angular endpoint), App. C (weak-norm obstruction), App. E (pencils) | **all proved here:** B in [14](14_GAPLESS_ANGULAR_ENDPOINT.md), C in [10](10_GAPLESS_REMAINING.md) §1, E in [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) | Lower priority: not used by the formal targets. For B, check the cited harmonic analysis (14 §0.3). |
| Prior OpenAI Sept preprints [62]–[64] and the 3D companion | not reviewed | Compare methods and check whether the October paper supersedes them. |

## Open questions I would value your view on

1. Is Theorem II (energy balance) best seen as a *linear* analogue of Serre's compensated integrability? Is there a cleaner abstract inequality behind it?
2. Can the holomorphic flux map (2.9) be described intrinsically, e.g. on the complexified characteristic variety?
3. Remark 6.2 shows two logarithms are sharp for scalar supersolutions. Is one logarithm reachable with the actual energy profile? The known lower bound is $1/p_c$.

## Also relevant: the ETH paper's formal proof exists

`openai/math` `lean/ComparatorChallenges/SubexponentialSAT.lean` has the same SHA-256 as the user's statement file. Its solution module is `OAI.Algorithms.SubexponentialSAT.Main`. See [`../eth/06_LEAN_STATEMENT_AUDIT.md`](../eth/06_LEAN_STATEMENT_AUDIT.md) for my follow-up.
