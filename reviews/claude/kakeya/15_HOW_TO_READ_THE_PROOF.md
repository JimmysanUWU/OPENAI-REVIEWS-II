# 15. How to read the proof

*Claude, 29 October 2026. A reading guide to* Local smoothing for the Euclidean wave equation *(OpenAI, 25 Oct 2026) and to this review's versions of its proof.*

This guide goes from no mathematics at all to the full proof, in four levels. Each level stands on its own, so stop wherever you like. Every step names where it is proved, in the paper and in this repository, and which script checks it.

| Level | For readers who know | What you get | Time |
|---|---|---|---|
| [0](#level-0-the-problem-and-the-idea) | nothing in particular | what the problem is, why it was hard, the idea of the proof | 15 minutes |
| [1](#level-1-the-proof-in-seven-moves) | calculus, complex numbers, a little linear algebra | the proof as seven moves, one formula each | an hour |
| [2](#level-2-the-proof-on-one-page) | graduate analysis | the proof on one page; where every number comes from; the standard names of the objects | an afternoon |
| [3](#level-3-reading-routes) | the field, or Lean | reading routes through the paper, the gapless files and the formal proof; notation traps; checkpoints | days |

Then come [questions readers ask](#questions-readers-ask) and [what to reuse](#what-to-reuse).

---

## Level 0: the problem and the idea

**Waves pile up.** Think of a wave (sound, light, ripples on a pond) as a crowd of small packets. Each packet moves in its own direction, like the beam of a flashlight. Packets can cross the same spot at the same moment, and there they add up. That is *focusing*. A lens does it on purpose; a whispering gallery does it by accident.

**Measuring the piling-up.** The plain average of a wave's energy never changes, so it says nothing about focusing. Mathematicians therefore use averages that weigh tall peaks more than gentle swells: the $L^p$ norms, where a larger $p$ means more weight on peaks. The question is how much more peaked the wave can be, averaged over a short stretch of time, than it was at the start.

**The critical setting.** One choice of $p$ is special: $p=3$ in three-dimensional space ($2n/(n-1)$ in $n$ dimensions). There, every known bad example sits exactly on the edge:
- a wave converging to a single point;
- packets lined up along needles pointing in every direction.

Christopher Sogge conjectured in 1991 that at this setting the time-averaged wave is essentially no more peaked than the starting wave. **The paper proves this**, with only a tiny extra factor that grows like a logarithm. Through implications known since the 1990s, this also settles the Bochner–Riesz, Fourier restriction and Kakeya conjectures in every dimension.

**The old road (1991–2025).** The standard plan was:
1. chop the wave into packets, each living on a thin tube;
2. understand how tubes pointing in different directions can overlap. This is the Kakeya problem in disguise;
3. repeat at every scale from coarse to fine.

Each scale costs a little, and step 2 needed exactly the Kakeya information nobody had. Even Kakeya in three dimensions took until 2025.

**The new road: bookkeeping with a quantity that never goes negative.** Suppose a quantity is never negative and is *conserved*: it moves around but never appears from nowhere. Then it cannot pile up beyond what flows in. Water in pipes is the picture: a junction cannot hold more than reached it.

Waves seemed to forbid this. Their energy is conserved, but in more than one dimension its flow is not a positive count; contributions cancel. The paper gets around this with five ideas:
1. **Directional sunglasses.** Evaluate the wave at *complex* points, shifted along light rays. This keeps only the part moving in a chosen direction, at a chosen resolution.
2. **Square after summing.** Combine the filtered views, then square. The result is never negative, and in a space with a few extra variables it is exactly conserved.
3. **Quaternions as packaging.** Pairs of complex numbers become single quaternions, which plug into a rigid branch of geometry.
4. **A curvature meter.** For each position and time, solve a geometric equation (Monge–Ampère) that builds a bowl whose curvature matches the density. The bowl's depth grows like a specific power of the density, and adding up depths gives exactly the critical measurement $p=2n/(n-1)$. Nobody chooses this exponent; the geometry produces it.
5. **A temperature comparison.** Assemble the bowls into one "energy" depending on three extra parameters. It obeys a rule like temperature in a room, so its value where the wave is observed is controlled by its value where the starting data sit. The two places are joined through a long thin passage, and squeezing through costs exactly two logarithms. That is the tiny extra factor.

**Where Kakeya went.** No tubes, needles or overlap counting appear anywhere in the argument. Kakeya comes out at the end, by known implications. (An analogy: Dvir's 2008 proof of the finite-field Kakeya problem also bypassed the combinatorics entirely.)

**What is new.** None of the ingredients is new on its own:
- Monge–Ampère equations and their energies come from complex geometry;
- complex points as filters come from signal analysis;
- positive conserved tensors come from gas dynamics;
- capacities come from potential theory.

What is new is the bridges between them, each built to remove one specific obstacle. No one had used any of these on these problems before (03 §3.3–3.5).

---

## Level 1: the proof in seven moves

Notation: $n$ is the space dimension, $d=n-1$, and $p=p_c=2n/(n-1)$. The wave is $u(T)=e^{-\mathrm iT|D|}f$. The tower of consequences sits under Move 7.

```
 data f ──(2,3)──▶ positive conserved current A ──(4,5)──▶ fibre potentials u ──(6)──▶ energy F(T,τ)
                                                                                         │
   F/a² ≈ ‖f‖_p^p on the data circle       F/a² ≈ ‖u(T)‖_p^p on the observation circles
                     └──────────(7: ΔF ≤ tiny·F, two passages)──────────┘
                    ⟹  ∫₁² ‖u(T)‖_p^p dT ≤ C log²(1/h) ‖f‖_p^p
```

**Move 1. Localize.** Restrict to one frequency scale $|\xi|\approx1/h$ and a narrow cone of directions, and put the problem on a large torus. The target is the **cap estimate**
$$\int_1^2\|e^{-\mathrm iT|D|}f\|_{L^p}^p\,dT\le C\log^2(2/h)\,\|f\|_{L^p}^p .$$
Covering the sphere of directions by caps and the frequencies by dyadic shells is standard.
*Where:* paper §6–§7; [08](08_GAPLESS_PROOF.md) §8–§9.

**Move 2. Complex points are direction filters.** A plane wave at a complex point gains a real factor:
$$\big|e^{\mathrm i\xi\cdot(X-\mathrm iq)}\big|=e^{\xi\cdot q}.$$
In light-cone coordinates, with the imaginary shift $q=y\,(1,W,-\frac12|W|^2)$ along a light ray, a wave of direction $v$ is multiplied by $e^{-\frac{\kappa y}2|W-v|^2}$. That is a Gaussian, largest when the wave's direction $v$ matches the chosen $W$. The size $y$ of the shift sets the resolution.
*Idea:* Kaiser's "physical wavelets"; for the Schrödinger equation this is exactly the Husimi (Gabor) density.
*Where:* paper §2.1; 08 (3.3); [16](16_NULL_FRAME_AND_SCHRODINGER.md) Proposition 16.5. *Check:* `verify_null_frame.py` Part H.

**Move 3. Square after summing.** Add the filtered plane waves first, then square. In lifted variables (the shift becomes the flux of an auxiliary state $w$) the result is a matrix field $A=\chi\chi^*$ with
$$A\ge0,\qquad\operatorname{div}A=0,\qquad A_{tt}=y^{-d}\,|u_0(X-\mathrm iq)|^2 ,$$
where $u_0$ is the wave with a harmless weight ($|D|^{d/2}$ in light-cone form).
*Idea:* summing before squaring keeps the interference, which is the focusing we want to measure. Squaring makes the quantity positive, and the lift makes it conserved. The wave's ordinary energy tensor is not positive in dimension $\ge2$.
*Where:* paper §2.2–2.3; 08 §3; 16 §4. *Check:* `verify_wave_current.py`, `verify_wave_current_target.py`.

**Move 4. Quaternionic bookkeeping.** Each spacetime coordinate $X^\mu$ gets a *slot* $X^\mu+\omega_\mu\mathrm i$, with $\omega_\mu$ a complex number: three real coordinates. On slots there is a Hessian $H(u)$, a determinant (Moore's) and a cofactor $C(H(u))$, and
$$\operatorname{div}C(H(u))=0\quad\text{for every }u .$$
*Idea:* this is the slot version of "$(dd^cu)^{n-1}$ is closed" in complex geometry. It makes the Monge–Ampère operator a divergence, so energies can be integrated by parts.
*Where:* paper §3; 08 §2; [12](12_QUATERNIONIC_SLOT_CALCULUS.md). *Check:* `verify_quaternion_calculus.py`, `verify_slot_calculus.py`.

**Move 5. A Monge–Ampère meter.** For each value of the extra parameters, solve on a "fibre"
$$\det H_{\rm fib}\big(\lambda|\Phi|^2+u\big)=\lambda\,\mu,\qquad u=0\ \text{on the wall},\qquad\mu=A_{tt}.$$
The background $\lambda|\Phi|^2$ is huge in one direction, so the determinant must come from the remaining $n-1$ directions. A convex bowl with curvature density $\mu$ in $n-1$ directions has depth $\approx\ell^2\mu^{1/(n-1)}$ on a ball of radius $\ell$, so
$$\int(-u)\,\mu\;\approx\;\int\mu^{1+\frac1{n-1}}\;\approx\;\int|\text{wave}|^{\frac{2n}{n-1}} .$$
*Idea:* the exponent is the determinant's degree; no $p$ is chosen anywhere.
*Where:* paper §4.2, §5, App. A; 08 §5–§7; 11 §4 (Steps 3, 5); toolkit K3. *Check:* `verify_propagation_and_exponents.py`.

**Move 6. An energy that is almost harmonic in three parameters.** Let the parameters be the time $T$ and an imaginary time $\tau\in\mathbb R^2$. Put
$$\mathcal F(T,\tau)=\int(-u)\,\mu+\frac1\lambda\int u\int_0^1\det H_{\rm fib}\big(\lambda|\Phi|^2+\theta u\big)\,d\theta .$$
Then $\Delta_{T,\tau}\mathcal F\le$ (flux through the fibre wall) $\le Cs^{-2}e^{-cs/h}\mathcal F$, and $\mathcal F\approx\int(-u)\mu$.
*Idea:*
- conservation turns parameter derivatives into fibre divergences;
- a Schur complement shows the current sees at least $\lambda^{-1}$ times the full determinant;
- the second term, the Aubin–Mabuchi energy of Kähler geometry, has exactly that determinant as its Laplacian, so the two cancel;
- what is left is a flux through the fibre wall, which the Gaussian of Move 2 makes exponentially small.

*Where:* paper §4.3–4.4, §5; 08 §6; 11 §3. *Check:* `verify_quaternion_calculus.py` (the energy identities).

**Move 7. Two circles and a thin cylinder.** In the three parameters $(T,\tau)$:
- the data sit on a small circle $|\tau|=4a$ at $T=0$;
- the wave at times $T\in[1,2]$ sits on small circles $|\tau|=2a$.

A nearly superharmonic positive function relates them through harmonic measure:
$$\mathcal F(0,\cdot)\ \gtrsim\ \frac1{\log^2(1/h)}\int_1^2\mathcal F(T,\cdot)\,dT .$$
*Idea:* Brownian motion must escape a thin cylinder and then hit another; each passage costs a logarithm, which is the capacity of a thin cylinder. Only a two-dimensional imaginary time gives logarithms (toolkit K4); in other dimensions the cost is a power.
*Where:* paper §6; 08 §8; toolkit K4–K5. *Check:* `verify_propagation_and_exponents.py`, `verify_capacity_dimension.py`.

**Down the tower.**
- Caps and shells assemble into the estimate on $\mathbb R^n$ (Theorem 1.2).
- Interpolation gives Sogge's conjecture with $\varepsilon$ derivatives (Theorem 1.1).
- Known implications give Bochner–Riesz, restriction, the Kakeya maximal conjecture, Kakeya sets, Nikodym and Stein's square function.

*Where:* paper §7–§12; 08 §9; [09](09_GAPLESS_CONSEQUENCES.md); [10](10_GAPLESS_REMAINING.md).

---

## Level 2: the proof on one page

**The chain.** Write $g_T=e^{-\mathrm iT|D|}f$ for cap data at scale $h$, $a=Kh$ with $K$ large, $b=2a$, and $\|f\|_2=1$. Proposition 7.1 of 08 (11 Proposition 4.1 in abstract form) gives three estimates:
$$\text{(i) }\ \mathcal F(0,\tau)\le Ca^2\|f\|_p^p\ \ (|\tau|=4a),\qquad\text{(ii) }\ \mathcal F(T,\tau)\ge ca^2\|g_T\|_p^p\ \ (|\tau|=2a),\qquad\text{(iii) }\ \Delta\mathcal F\le Cs^{-2}e^{-cs/h}\mathcal F .$$
- (i) is the barrier $C\sigma^{m}\,\mathrm r\,P_*$ in the comparison principle.
- (ii) inverts the Gaussian at the fibre's centre and uses a quadratic barrier on a ball.
- (iii) is the energy balance plus Gaussian decay on the wall.

The two-circle lemma (08 Lemma 8.1) turns (iii) into
$$\mathcal F(0,2b)\ge c\,\frac{\log2}{\log^2(3/b)}\int_1^2\mathcal F(T,b)\,dT .$$
Chaining with (i)–(ii), $a^2$ cancels:
$$\int_1^2\|g_T\|_p^p\,dT\le C\log^2(2/h)\|f\|_p^p.$$

**Every number, and where it comes from.**

| Number | Value | Why | Where |
|---|---|---|---|
| critical exponent | $p_c=2n/(n-1)$ | Monge–Ampère in $n-1$ free directions; the background takes one | 11 §8.1; toolkit K3 |
| auxiliary copies | 2, combined as $z=w_1+\mathrm kw_2$ | a full Gaussian window in $W$; a *holomorphic*, hence harmonic, flux map, so pushing forward keeps conservation | 04 §4.3; 08 Lemma 3.2 |
| real coordinates per slot | 3: $(1,\mathrm i,\mathrm j)$ | one real coordinate plus one complex auxiliary coordinate | 08 §3.4 |
| dimension of the imaginary time $\tau$ | 2 | the auxiliary part of the time slot, and the only dimension where a thin cylinder's capacity is logarithmic | toolkit K4 |
| background $\lambda$ | $h^{-30n}$ (any large power) | must dominate $\mu^{1/d}$; cancels exactly | 08 Lemma 6.1(4) |
| energy constants | $\frac{n-1}{n+1}\le\mathcal F/X_0\le\frac n{n+1}$ | rank-one expansion of the Moore determinant along $\theta\mapsto B+\theta H(u)$ | 08 Lemma 6.3 |
| logarithms | $\log^2$ in the $p$-th power, i.e. $\log^{2/p}$ in the norm | two passages through a thin cylinder; sharp for scalar supersolutions | 08 Lemma 8.1; paper Remark 6.2 |
| lower bound | one logarithm | the focusing example | paper Prop. 8.1; 16 §6 |
| aperture $\gamma$ | a small constant | wall decay $e^{-c_0\sigma}$ must beat the inversion cost $e^{C_2\gamma^2\sigma}$ | 11 (4.0) |
| regularizer | $h^{50n}C(H(\mathrm r))$ | keeps the Dirichlet problem non-degenerate without touching the balance | 08 (6.1) |
| formal power | $\log^{2n+5}$ | the Lean proof uses a cruder propagation step | 06 |

**The standard names of the objects.**

| In the paper | Standard name |
|---|---|
| conserved Hermitian field $A\ge0$ on slots | a closed positive current. In complex geometry, a divergence-free positive Hermitian matrix field is a closed positive $(n-1,n-1)$-form. |
| cofactor $C(H(u))$ is divergence-free | $(dd^cu)^{m-1}$ is closed (Alesker, for quaternionic slots) |
| Green's identity for conserved fields | Stokes' theorem for closed currents |
| $E_z(u)=\frac1\lambda\int u\int_0^1\det D_{\theta u}$ | the Aubin–Mabuchi (Monge–Ampère) energy |
| $\Delta_zE_z(u)=\frac1\lambda\int\det H(B_0+u)$ | "$dd^c$ of the energy of a family = fibre integral of the total Monge–Ampère measure" |
| $A_{tt}=y^{-d}\lvert u_0(X-\mathrm iq)\rvert^2$ | boundary values of the holomorphic extension to the forward tube; for Schrödinger, the Husimi density (16) |
| two-circle lemma | harmonic measure and the capacity of a thin cylinder |
| depth $\propto\mu^{1/d}$ | an Alexandrov-type estimate for Monge–Ampère |
| the method as a whole | a *linear* compensated integrability. Serre's inequality is multilinear in $A$ and vanishes for the rank-one wave current. |

**What each hypothesis does.** See the hypothesis-usage matrix, 11 §7. In one line each:
- **positivity** gives a solvable Monge–Ampère equation;
- **conservation** gives Green's identity;
- **tangency** kills the huge background in the balance;
- the **rank-one background** fixes the exponent;
- the **Gaussian window** makes the wall flux tiny;
- the **two-dimensional imaginary time** makes the loss logarithmic.

**The simplest instance.** In light-cone coordinates ([16](16_NULL_FRAME_AND_SCHRODINGER.md)):
- the chart is $P=YW$, $S_1=aQ-\frac12YW\cdot W$, with Jacobian exactly $aY^d$;
- the time row is one Gaussian multiplier;
- the density is literally the wave continued to a complex point.

Read 16 §1–§4 next to 08 §3–§7: everything in the paper's physical frame is this after a rotation of two coordinates. On waves of one "mass" the same engine is a direct proof for the Schrödinger equation (16 Theorem 16.8).

**Which equations the method handles.** [17](17_STATE_DIMENSION_RIGIDITY.md) gives a linear-algebra test: a system can carry a tangent lifted current only if $z,C_1z,\dots,C_nz$ are linearly dependent for every complex state $z$. With an extremal sheet of multiplicity $m$, this is exactly the condition $N'=n+m-1$. The wave passes; Maxwell's equations do not.

---

## Level 3: reading routes

### Route A: the shortest complete proof (planar, linear PDE only)

Sharp $L^4$ local smoothing in $\mathbb R^{2+1}$ with $\log^{1/2}$ loss. There is no Monge–Ampère theory: the fibre problem is a linear Dirichlet problem on a disk.

1. [13](13_PLANAR_ENTRY_POINT.md) §0–§1: the set-up and the Dirac-pencil chart.
2. [12](12_QUATERNIONIC_SLOT_CALCULUS.md) §1–§4: the slot calculus.
3. [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Lemma 6.4: the current.
4. 13 §2–§4: the disk operator, the linear potential, the planar energy.
5. 08 Lemmas 7.2 and 8.1: kernels and the two circles.
6. 08 §9: the transfer to $\mathbb R^2$.

This is the reading path `P.KAK.PLANAR` in [PROOF_MAP.md](../../../PROOF_MAP.md); `python3 tools/check_roadmap.py --path P.KAK.PLANAR` prints it.

### Route B: every dimension, distilled

1. 12: the calculus.
2. 11 §1–§5: the engine theorem from axioms (G1)–(G2) and (A1)–(A4).
3. 16 §1–§3: the null frame, the cleanest instance of the axioms.
4. 11 §6: pencils, and the wave as the arrow pencil.
5. 08 §9, then 09: the transfers and consequences.

This is `P.KAK.ENGINE`, with 16 inserted as a worked instance. 17 then says which other systems fit.

### Route C: the paper's own order

| Paper | Content | Gapless version here | Watch for |
|---|---|---|---|
| §1.4 | the idea | 03, 04, this guide | |
| §2 | the positive conserved current | 08 §3; 16 §4 | the conjugation in $G$; the branch of $\sqrt Y$ |
| §3 | the quaternionic calculus | 08 §2; 12 | slots have three real coordinates |
| §4 | fibres and the determinant energy | 08 §4, §6 | the common null line on the wall |
| App. A | the Dirichlet problem | 08 §5 (Schauder and Evans–Krylov cited; also formally verified) | the fourth slot appears only here |
| §5 | estimates on one fibre | 08 §7; 11 §4 | the order of constants: $\gamma$, then $L$, then $K$, then $h$ |
| §6 | two-circle propagation | 08 §8; toolkit K4–K5 | Remark 6.2: two logarithms are sharp for this step |
| §7 | the transfer to $\mathbb R^n$ and interpolation | 08 §9 | |
| §8 | obstructions | 10 §1; 16 §6 | |
| §9–§11 | restriction, Kakeya, Nikodym, Bochner–Riesz, Stein | 09; 10 §2–§3 | the Kakeya maximal power $4/p_c$ |
| §12 | the paraboloid and Schrödinger | 10 §4 (via the sphere); 16 (directly) | |
| App. B | the angular endpoint | [14](14_GAPLESS_ANGULAR_ENDPOINT.md) | the threshold $B>\frac12+\frac1p$ |
| App. C | the weak-norm obstruction | 10 §1 | |
| App. D | the planar case | 13 | |
| App. E | symmetric pencils | 11 §6; 17 for the scope | |

### Route D: the formal proof

1. [`source/kakeya/LocalSmoothing.challenge.lean`](../../../source/kakeya/LocalSmoothing.challenge.lean): the six statements, each `sorry` in the challenge. They are:
   - `euclidean_local_smoothing`;
   - `bochner_riesz`;
   - `spherical_restriction`;
   - `kakeya_maximal`;
   - `kakeya`;
   - `stein_square`.
2. [06](06_LEAN_AUDIT.md):
   - OpenAI's proof builds here, from 2,261 modules;
   - it uses only `propext`, `Classical.choice` and `Quot.sound`;
   - its elaborated statements are byte-identical to the challenge's;
   - the Dirichlet and Schauder theory is formalized, not assumed.
3. What is *not* formal:
   - the sharp power $2/p_c$ (the formal power is $2n+5$);
   - App. B;
   - everything in 16–17.

### Notation traps

1. **The bar in $G=\overline{f(\mathsf X-\mathrm iq)}$.** The other sign is not conserved (08 Lemma 3.1, remark; a negative control in `verify_wave_current.py`).
2. **Two complex structures.** $\mathbb C_{\mathrm i}$ carries Fourier phases and $\mathbb C_{\mathrm k}$ carries auxiliary coordinates. Quaternions combine them, with $\mathrm{ki}=\mathrm j$.
3. **Three real coordinates per slot.** A slot is $(1,\mathrm i,\mathrm j)$, not four coordinates; a fourth appears only in App. A, for boundary estimates.
4. **Bilinear versus Hermitian.** In holomorphic formulas dots are bilinear ($W\cdot W$). The Gaussian uses the Hermitian $|W-v|^2$.
5. **Two times.** The physical time is $T$, the null time is $t$, and the paper rotates $(t,s_0)\to(T,b)$. The unrotated frame is simpler (16).
6. **"Imaginary time".** The parameter is $(T,\tau)\in\mathbb R^3$. The circles $|\tau|=2a$ and $4a$ are at the scale of the uncertainty principle, with $\sigma=s/h$ and $a=Kh$.
7. **$\lambda=h^{-30n}$** is enormous, but it cancels exactly. It only has to dominate $\mu^{1/d}$.
8. **$\chi$ does not depend on the branch of $\sqrt Y$** (08 Prop. 3.3, step 4).
9. **The Monge–Ampère determinant** is Moore's determinant of a quaternionic Hermitian matrix, the product of its real eigenvalues.
10. **The formal theorem has $\log^{2n+5}$**, not the sharp $\log^{2/p_c}$ (06; START_HERE fact 4).

### Checkpoints

You have understood a step if you can answer its question. The script checks the identity behind the answer.

| Step | Question | Script |
|---|---|---|
| Move 2 | Why does a complex shift along a light ray produce a Gaussian in the direction? | `verify_null_frame.py` H |
| Move 3 | Why is $\operatorname{div}\mathcal Z_\ell=2\mathrm iw_\ell^{\mathsf T}\mathscr PG=0$, and why does the "acceleration" drop out of the real part? | `verify_wave_current.py` |
| Move 3 | Why is the pushed current conserved, independently of the push-forward argument? | `verify_wave_current_target.py` |
| Move 4 | Why is the cofactor of every Hessian divergence-free? | `verify_quaternion_calculus.py`, `verify_slot_calculus.py` |
| Move 5 | Why is the exponent exactly $2n/(n-1)$? | `verify_propagation_and_exponents.py` |
| Move 6 | Why does $\Delta_zE_z(u)$ equal $\lambda^{-1}\int\det H(B_0+u)$, with no boundary term? | `verify_quaternion_calculus.py` |
| Move 7 | Why two logarithms, and why only in imaginary time of dimension 2? | `verify_propagation_and_exponents.py`, `verify_capacity_dimension.py` |
| engine | Which hypotheses does each step use? | 11 §7; `verify_positivity_engine.py` |
| scope | Can Maxwell's equations carry the current? (No: 17.) | `verify_state_dimension.py` |

`python3 tools/verify_all.py --project KAK` runs all of them.

---

## Questions readers ask

**Why quaternions?**
- Two auxiliary flows give two complex rank-one terms. Together they form one quaternionic rank-one term (04 §4.3).
- Each slot then has three real coordinates, so the time slot carries a two-dimensional imaginary time.
- Two is exactly the dimension in which a thin cylinder has logarithmic capacity (K4).

So the number of auxiliary flows and the size of the loss are tied together.

**Why is $\lambda$ so large?** It makes the equation see only the $n-1$ directions transverse to the background, which fixes the exponent. It cancels exactly in the energy balance, so its size never enters the result.

**Could the method prove something false?** Not above $p_c$. The exponent is the determinant's degree in $n-1$ directions; there is no free parameter to push. Above $p_c$ the estimate without derivative loss is false (the focusing example), and the method never claims it.

**Why two logarithms, when one is the truth for the focusing example?** Two harmonic-measure passages cost one logarithm each, and the paper's Remark 6.2 shows this is sharp for any positive superharmonic function. Doing better needs more information about the energy than superharmonicity (`T.KAK.6`). 16 §6 shows numerically, for Schrödinger, that focusing produces exactly one logarithm.

**Where does the equation enter?** Only through the current: positivity, conservation, tangency and Gaussian rows, the axioms (A1)–(A4) of 11. The rest of the proof never mentions the equation (11 §8.4). 16 runs the same engine for Schrödinger, and 17 says which systems can carry it.

**Why not on curved manifolds?** The current needs exact holomorphic extensions of exact solutions, which constant coefficients provide. On some curved metrics the critical estimate is false (Bourgain; Minicozzi–Sogge), so something must fail there. By 11 §8.4 it is the current (`T.KAK.8`).

**How sure is this?**
- **The six formal targets:** OpenAI's Lean proof builds with standard axioms only (06). About 99%, subjective.
- **The sharp power $2/p_c$:** proved in prose here (08 §8), not formalized. About 90%, subjective.
- **The consequences:** derived in 09–10.

See [01](01_VERDICT.md).

**How was it found?** The paper is an OpenAI release in the `openai/math` repository, authored as "OpenAI". No account of how the proof was discovered, such as a reasoning trace, is published, so this review can only assess the mathematics.

---

## What to reuse

| Tool | One line | Where |
|---|---|---|
| K1. Positive conserved phase-space current | Constant-coefficient symmetric systems carry a positive, exactly conserved current on a lifted space. Its time component is the continued wave, which for Schrödinger is the Husimi density. | 08 §3; 16 §4 |
| K2. Determinant energy balance | A positive conserved current, a rank-one background and designed fibres give $\Delta\mathcal F\le$ wall flux. | 08 §6 |
| K3. Density to depth | Monge–Ampère depth $\propto\mu^{1/f}$ with $f$ free directions. | 07 K3 |
| K4. Capacity–dimension | Transport between two spheres costs $\log^{-2}$ only in two-dimensional auxiliary time. | 07 K4 |
| K5. Two-circle propagation | $\Delta E\le VE$ on a thin cylinder transports averages at cost $\log(r/b)/\log^2(3/b)$. | 08 §8 |
| K6. Cofactor conservation | Moore cofactors of frame Hessians are divergence-free for any quaternion frames. | 12 |
| K7. Positivity engine | Axioms (G1)–(G2), (A1)–(A4) imply the critical estimate. | 11 |
| K8. The null frame and Schrödinger | The engine in light-cone coordinates has a polynomial chart, and its single-mass sector is the Schrödinger flow. | 16 |
| K9. The tangency test | A lift is tangent iff $\operatorname{rank}[z\,\lvert\,C_1z\,\lvert\cdots]\le n$; equivalently $N'=n+m-1$. | 17 |

The roadmap ([ROADMAP.md](../../../ROADMAP.md) §4, §6) lists these frameworks and the open targets they lead to.
