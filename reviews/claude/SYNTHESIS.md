# Synthesis: what the two OpenAI results share, and what can be recycled

Claude, 27 October 2026. **Both headline results are now formally verified**: OpenAI's Lean proofs build locally, use only the standard axioms, and have elaborated statements byte-identical to the challenge files ([`eth/06`](eth/06_LEAN_STATEMENT_AUDIT.md) §6.A, [`kakeya/06`](kakeya/06_LEAN_AUDIT.md)). Read with [`eth/07_TOOLKIT.md`](eth/07_TOOLKIT.md) and [`kakeya/07_TOOLKIT.md`](kakeya/07_TOOLKIT.md). **Gapless proofs:** [`eth/08`](eth/08_GAPLESS_PROOF.md), [`eth/09`](eth/09_GAPLESS_SHARP_EXPONENT.md) and [`eth/10`](eth/10_GAPLESS_TENSOR_CONTRACTION.md); [`kakeya/08`](kakeya/08_GAPLESS_PROOF.md), [`kakeya/09`](kakeya/09_GAPLESS_CONSEQUENCES.md) and [`kakeya/10`](kakeya/10_GAPLESS_REMAINING.md).

**Map of all results:** [`ROADMAP.md`](../../ROADMAP.md) (stable IDs, dependency graphs, statuses, open targets).

## 1. One pattern: move information along an auxiliary parameter, and pay a capacity

Both papers must estimate a quantity at a place where it is hard to control. In ETH that is the count, at $T=0$. In Kakeya it is the wave norm, on the observation circle. Both introduce an **auxiliary parameter**, make the quantity tractable elsewhere in that parameter, and transport the information back. The cost is a *capacity* of the transport.

| | ETH (counting) | Local smoothing (Kakeya) |
|---|---|---|
| Auxiliary parameter | real $T$ deforming the substitutions $x\mapsto P(T)x$ | imaginary time $\tau\in\mathbb C_{\mathbf k}\cong\mathbb R^2$ (with physical $T$) |
| Where it is tractable | the nodes $t=1,\dots,Q$: sums split into independent groups | the energy's own estimates on the data and observation circles |
| Certificate there | separation by few integer measurements (Def. 3.1) | a positive conserved current and a Monge–Ampère potential (K1–K3) |
| Transport back | Hermite interpolation of a polynomial of degree $dR$ from $Q$ nodes | superharmonic comparison of the energy from circle to circle |
| Transport cost | Taylor depth $q\approx nR/Q$, i.e. $2^{n(R/Q)\log(Q/R)}$ | harmonic measure of a thin cylinder: $\log^{-2}$ |
| What fixes the cost | the agreement degree of parity: degree $R\gtrsim\sqrt{Q\log(1/\varepsilon)}$ for any polynomial, by Smolensky's bound (E1, E4, [eth/14](eth/14_DESIGN_FLOOR.md)) | planar capacity: only $\dim\tau=2$ gives logarithms (K4) |
| Resulting barrier | $(\log\log n)^2/\sqrt{\log n}$ is optimal for every graded design (E4; eth/14, proved) | $\log^2$ is optimal for scalar supersolutions (Remark 6.2, K4–K5) |

The analogy is structural, not a theorem. It does locate both results' limits in the same place: the **capacity of the transport along the auxiliary parameter**. Each paper's main barrier is a sharp statement about that capacity, and each is proved sharp here or in the paper.

## 2. A second shared device: involutions and their joint spectrum

- **ETH.** Commuting coordinate swaps $\Pi_s$. Everything is decided by the joint spectrum: $M_0$ acts on $V_\chi$ by $(-1)^{|\chi|}f(|\chi|)$, and the defect is the binomial tail (E1).
- **Kakeya.** The quaternion units act on the slots, and the two auxiliary flows combine through $\mathbf k$. The Gram matrix of two complex rank-one terms is one quaternionic rank-one term (K1, §2.3 of the paper).

In both, a small finite algebraic structure turns a messy object into a diagonal or rank-one one, whose spectrum or determinant carries the exponent:
- ETH: the parity character on the binomial bulk;
- Kakeya: the Moore determinant's degree $n$, minus one background direction, giving $p_c$.

## 3. What is genuinely recyclable (ranked by expected reach)

1. **K1: the positive, exactly conserved phase-space density** for constant-coefficient symmetric systems. It is a new object with no previous analogue: Husimi densities are not conserved, and Wigner functions are not positive. Any "positivity + conservation" argument in dispersive PDE can start from it.
2. **K2: the determinant energy balance.** A wave-free engine: positive current + rank-one background + designed fibers ⇒ almost-superharmonic energy. It imports Kähler-geometry energies (Aubin–Mabuchi) into harmonic analysis. **Now a theorem with precise axioms** ([`kakeya/11`](kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) Theorem 11.1). Symmetric pencils are instances, and the wave is the arrow pencil.
3. **E3: deformation counting.** A general exact algorithm for partition functions of bounded-degree factor graphs, in $2^{o(n)}$. Its new corollary covers Ising/Potts partition functions with rational weights.
4. **E1: commuting-involution interpolation.** A representation-theoretic lemma: a low-degree path from $I$ to "local" group elements, with an exact defect formula, optimal among symmetric choices.
5. **E2: $p$-adic cross-equation truncation** (new). Lifting $\mathbb F_p$ systems to $\mathbb Z/p^\ell$ at *linear* degree, for every prime.
6. **K4: the capacity-dimension principle** (new). A design rule: auxiliary parameters that are transported by superharmonic comparison must be two-dimensional to cost only logarithms.
7. **E5: tapered 2-adic Hermite recovery.** Exact interpolation at many nodes with precision demand tapered by 2-adic depth.

## 3a. Frameworks extracted in round 4 (28 Oct)

| Result | File | Status |
|---|---|---|
| The positivity engine as a theorem: axioms on a chart and a current imply the critical estimate; pencils are instances; the wave is the arrow pencil | `kakeya/11` | proved; pencil current checked to $10^{-38}$ |
| Cofactor divergence for arbitrary quaternion frames; the Monge–Ampère operator is the gradient of a concave energy; at most 4 coordinates per slot | `kakeya/12` | proved; exact checks |
| Sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ by linear PDE only (App. D in full) | `kakeya/13` | proved; closed prose chain |

**Bearing on §1.** In the engine theorem the capacity of the transport enters only through Lemma 8.1 and the dimension of the auxiliary time (11 §7). The equation enters only through the current axiom (A1) and two Gaussian symbol estimates (A3)–(A4).

## 4. Sharpenings proved in this round

| Result | File | Status |
|---|---|---|
| $p$-adic cross-equation truncation, all primes | `eth/07_TOOLKIT.md` E2 | proved; exact checks $p\le7$, $\ell\le6$ |
| Exact weighted partition functions in $2^{O(n(\log\log n)^2/\sqrt{\log n})}$ | `eth/07_TOOLKIT.md` E3′ | proved from the paper's Theorem A |
| Cost floor $(\log\log n)^2/\sqrt{\log n}$ for symmetric designs | `eth/07_TOOLKIT.md` E4 | proved; numerical optimum within $[0.11,0.52]\times$ the reference |
| The same floor for every graded design (Smolensky + restriction) | `eth/14_DESIGN_FLOOR.md` | proved; exact checks of the parity bounds |
| Density-to-depth with a rank-$r$ background | `kakeya/07_TOOLKIT.md` K3 | proved; the exponent-selection rule is heuristic |
| Capacity-dimension principle ($k=2$ unique) | `kakeya/07_TOOLKIT.md` K4 | proved; exact checks $k=1,\dots,5$ |
| Exact asymptotic constant $\tfrac94\log2$ in Remark 6.2 | `kakeya/07_TOOLKIT.md` K5 | verified to 4 digits |

## 5. Research directions that these toolkits make concrete

1. ~~**Non-symmetric swap designs (ETH).**~~ **Done (29 Oct):** [`eth/14`](eth/14_DESIGN_FLOOR.md). The least $p$-biased measure is $e^{-\Theta(R^2/Q)}$ up to $\mathrm{poly}(Q)$, so no graded design beats E4. Open: ungraded designs (`T.ETH.6`), where far group elements carry node-dependent coefficients.
2. ~~**Counting over $\mathbb F_p$ (ETH).**~~ **Done (28 Oct):** [`eth/12`](eth/12_FP_COUNTING.md) Theorem 12.1 counts sparse systems over any fixed finite field in $2^{O(\rho(n))}$.
3. **Beyond the light cone (Kakeya).** Extend App. E to pencils whose size differs from the spatial dimension (K7.1), connecting to determinantal representations of hyperbolic polynomials.
4. **Degenerate symbols (Kakeya).** Find currents tangent to $r>1$ passive coordinates and prove the rank-$r$ analogue of Lemma 4.3. That would give critical estimates at $2(f+1)/f$ for symbols with $f$ curved directions (K3).
5. **One logarithm (Kakeya).** Beat Remark 6.2 by using the energy's structure, not just superharmonicity.
