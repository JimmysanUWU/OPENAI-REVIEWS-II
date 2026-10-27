# 3. Methods, leaps of insight, and how radical they are

## 3.1 What the field was doing

Since Bourgain's 1991 work, progress on restriction, Kakeya, Bochner–Riesz and local smoothing has come from **multiscale combinatorial geometry**:
- wave-packet decompositions and induction on scales (Wolff, Tao);
- bilinear and multilinear restriction and Kakeya (Tao–Vargas–Vega, Bennett–Carbery–Tao, Guth's endpoint);
- decoupling (Bourgain–Demeter, Bourgain–Demeter–Guth);
- polynomial partitioning (Guth, Hickman–Rogers–Zhang);
- the square-function method for the cone (Guth–Wang–Zhang);
- the Kakeya-specific structure theory of Katz–Tao and Wang–Zahl (stickiness, grains, Furstenberg estimates).

Each of these loses $R^\varepsilon$ or worse. None has reached the conjectures for $n\ge3$, apart from Kakeya sets in $\mathbb R^3$ (Wang–Zahl 2025, a very long argument).

## 3.2 What this paper does instead

The paper uses **no wave packets, no induction on scales, no decoupling, no polynomial partitioning, and no tube combinatorics**. It works at a single frequency scale and proves one inequality by **positivity and elliptic PDE in auxiliary variables**.

| # | Method | Where | Antecedent | What is new |
|---|---|---|---|---|
| M1 | **Imaginary displacement along light rays** as a frequency-direction filter: $e^{-\ell|\xi||\omega-\hat\xi|^2/2}$ | §1.4, §2.1 | Kaiser's analytic signals and Kaiser–Streater windowed Radon transform; FBI/Bargmann transforms; Husimi densities | Choosing the displacement to be the *flux of an auxiliary state of the symmetric system* (2.3), so the Gaussian width is set by the system. |
| M2 | **Lifted conserved current**: $B=\operatorname{Re}\sum_\ell\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ in $(\mathsf X,w_1,w_2)$, PSD and divergence-free (Lemma 2.1) | §2.2 | Friedrichs' bilinear flux for symmetric systems; Kaiser's conserved vector at fixed imaginary displacement; CR vector fields (Baouendi–Trèves) | A **positive** conserved *tensor* for a superposition of waves in every dimension. The ordinary stress tensor is not positive for $n\ge2$ (Serre). Positivity is obtained by squaring *after* summing the modes, in lifted variables. |
| M3 | **Holomorphic flux map to quaternionic slots**: two complex rank-one terms become one quaternionic rank-one Gram matrix on triples $X^\mu+\omega_\mu\mathbf i$ | §2.3 | none known | Producing a *quaternionic* Hermitian positive current, with a time component equal to a Husimi density of the wave (2.15). |
| M4 | **Quaternionic Monge–Ampère potential with a rank-one background** on each spatial fiber: $\det H_{\rm fib}(\lambda|\Phi|^2+u)=\lambda\mu$ | §4.2, App. A | Alesker's quaternionic MA and Moore determinant; Zhu's Dirichlet problem; CNS boundary estimates | Using MA as a **density-to-depth converter**: the background's single huge eigenvalue leaves $d=n-1$ free directions, so depth scales like $\mu^{1/d}$ and $\int(-u)\mu\asymp\int|{\rm wave}|^{2n/(n-1)}$. This is where $p_c$ comes from. |
| M5 | **Aubin–Mabuchi / J-functional energy** $\mathcal F=X_0+E_z(u)$ whose parameter Laplacian cancels the full determinant (Prop 4.2) | §4.3–4.4 | Kähler geometry: Aubin, Mabuchi, Berman–Boucksom–Guedj–Zeriahi; Åhag et al.'s slice energy; Serre's compensated integrability | Transplanting variational pluripotential energies to *moving* fibers in a Fourier-analytic problem, with an exact boundary cancellation from fiber design. |
| M6 | **Scalar supersolution and two harmonic-measure passages** in the 3-dimensional parameter space $(T,\tau)$; each costs one logarithm (Lemma 6.1) | §5–6 | radial harmonic functions and Brownian exit probabilities; weak Harnack | Converting the whole problem to a 3D elliptic comparison. The logarithmic loss is a capacity of a thin cylinder, not a count of scales. |
| M7 | Standard transfers: torus to $\mathbb R^n$, interpolation, Bourgain/BCSS randomization to Kakeya, Tao's ε-removal, Nicola's slicing, Seeger's dilate combination | §7, §9–12 | classical | Careful log bookkeeping (e.g. Hausdorff measure is infinite for the gauge $r^n\log^{2An}(\cdots)$). |

## 3.3 The leaps of insight

**Leap 1: a positive conserved tensor for waves.** Positivity plus conservation is the kind of structure that turns an $L^p$ estimate into a calculus inequality, as Serre's compensated integrability does for gas dynamics. The obstacle is that the wave's stress tensor is indefinite for $n\ge2$. The paper's resolution has three steps:
1. Analytically continue the wave in an imaginary direction chosen *by an auxiliary state*.
2. Use two auxiliary states.
3. Square after summing.

This gives a positive tensor on a $3(n+1)$-dimensional lifted space whose time component is the wave seen through a Gaussian direction window. I know of nothing like this in the restriction literature.

**Leap 2: Monge–Ampère as an $L^{p_c}$ meter.** Serre-type inequalities control $(\det A)^{1/(n-1)}$, a *multilinear* quantity, and the lifted current has rank one, so its determinant vanishes. The paper does not use $\det A$ at all. It uses $A_{TT}$ as the *right-hand side* of a fiberwise Monge–Ampère equation, with a rank-one background $\lambda|\Phi|^2$ supplying the missing direction.

The exponent $n/(n-1)$ then comes from the degree of the determinant in the remaining $d$ directions. This is the conceptual heart. It is *why* the method lands exactly on $p_c$, and why it cannot prove anything false above $p_c$.

**Leap 3: a parameter-space energy that is almost harmonic.** The Kähler-geometry formula "Laplacian of the energy of a family = fiber integral of the total Monge–Ampère measure" is adapted to moving quaternionic fibers. Combined with conservation of the current and of the cofactor, every interior term becomes a nonnegative $A:\Sigma$ or cancels. What survives is a boundary flux through the fiber wall, where the Gaussian window makes the wave exponentially small.

So the energy satisfies $\Delta\mathcal F\le Cs^{-2}e^{-cs/h}\mathcal F$: almost a positive superharmonic function in three variables.

**Leap 4: losses as capacities.** Comparing the energy at the data circle $|\tau|=4a$ and the observation circle $|\tau|=2a$ through the unit circle costs exactly two logarithms of the cylinder's aspect ratio. Remark 6.2 proves this is sharp for scalar supersolutions. It explains the $[\log]^{2/p_c}$ and why it sits between the lower bound $1/p_c$ and $2/p_c$.

## 3.4 How radical is it?

**In consequence: the largest result of the three we have examined, by a wide margin.** If correct, it resolves in all dimensions:
- the local smoothing conjecture;
- the Bochner–Riesz conjecture;
- Stein's square function conjecture;
- the Fourier restriction conjecture (strict range);
- the Kakeya maximal and Kakeya set conjectures;
- the Nikodym maximal conjecture.

These are the organizing problems of Euclidean harmonic analysis since the 1970s. The Kakeya set conjecture in $\mathbb R^3$ alone was a celebrated 2025 breakthrough.

**In technique: a paradigm change, not an assembly of known tools.** The ingredients each have a history in some field:
- quaternionic Monge–Ampère theory (Alesker 2003–2013);
- Kähler energy functionals (Aubin 1984, Mabuchi 1986);
- symmetric hyperbolic systems (Friedrichs 1954);
- analytic signals (Kaiser 1990);
- compensated integrability (Serre 2018).

None had been brought to bear on restriction or Kakeya, as far as I know. The combination is not suggested by any prior roadmap in harmonic analysis.

The ETH paper transplanted one known phenomenon (parity approximation) into a new architecture. This paper builds a new *bridge between whole fields*: pluripotential theory and Fourier restriction. The bridge carries the critical exponent exactly.

**Kakeya geometry never appears.** Tubes enter only in §9, through Bourgain's 1991 randomization from a restriction estimate. The paper proves the Kakeya conjecture without any Kakeya combinatorics. This is analogous to how Dvir's polynomial method bypassed the combinatorics of finite-field Kakeya, but in the real, analytic setting where no such shortcut was believed to exist.

## 3.5 Prior work the paper builds on (all its own references)

- OpenAI Sept 2026: Kakeya maximal in $\mathbb R^3$ [63]; full dimension in $\mathbb R^4$ [62]; "elliptic capacity propagation" for restriction [64]; critical 3D local smoothing.
  - The title of [64] suggests that the capacity and propagation idea (M6) predates the October paper.
  - The October paper is the general-dimension culmination.
- Serre [69–71]: compensated integrability. This is the nearest conceptual ancestor of "positive divergence-free tensor ⇒ integrability of a determinant power".

## 3.6 Is it independently creative?

Yes, at the level of the field's best work, in my judgement. The decisive ideas (Leaps 1–3) are not recombinations of lemmas from the restriction literature. They import structures from complex and quaternionic geometry and from hyperbolic PDE to create a new mechanism. The paper explains it, and I verified its algebraic core.

The 7/8 result extended a known analytic framework, and the ETH paper re-purposed a known phenomenon. This paper creates a new method. If its analytic estimates hold (and the formal proof indicates they do; see [06](06_LEAN_AUDIT.md)), it is the strongest evidence yet for the user's thesis.

## 3.7 Stitching or creation? Why it was not found earlier

**A test for stitching.** In a stitched proof, most of the content is cited theorems and the glue is routine. Here the cited inputs are a minority of the argument:
- Alesker's cofactor identity;
- Friedrichs' flux;
- the Aubin–Mabuchi first variation;
- CNS boundary estimates;
- Bourgain's randomization.

The rest consists of constructions invented to remove specific obstructions:

| Obstruction | New construction |
|---|---|
| The stress tensor is indefinite | Lifted current with two auxiliary flows |
| Two complex squares | One quaternionic rank-one object, via the holomorphic flux map |
| The current has rank one, so $\det=0$ | Monge–Ampère with a rank-one background, used as an $L^{p_c}$ meter |
| Boundary terms have no sign | Fiber design making the normal cofactor vanish |
| Control in parameters | An almost-harmonic energy in $(T,\tau)$ |

**A plausible discovery path.** This is speculation; no reasoning trace is published for this paper (`openai/math/reasoning_traces` has 11 traces, none for it).
1. Serre's compensated integrability controls $(\det A)^{1/(n-1)}$, with the exponent $n/(n-1)=p_c/2$. That numerical match invites the attempt.
2. Each obstruction then points to the next tool: positivity → phase-space (Husimi) densities; conservation → Friedrichs and analytic signals; rank one → Monge–Ampère with a background; parameter control → Aubin–Mabuchi energies.
3. The September releases ([62]–[64] and 3D critical local smoothing, including "elliptic capacity propagation") suggest an iterated research program, not a single leap.
4. A full formal proof makes a long, fragile chain worth attempting.

**Why it was not found earlier.**
1. **Belief.** Since Wolff and Katz–Tao, and reinforced by Wang–Zahl, Kakeya was believed to be fundamentally combinatorial.
2. **Distance between fields.** The ingredients are spread across quaternionic pluripotential theory, Kähler energy functionals, and gas-dynamics compensated integrability. Few people work across all of these and restriction theory at once.
3. **Cost of a trial.** One must solve a degenerate 3-slot Monge–Ampère problem on moving periodic fibers before learning whether the exponents close.
4. **Counterintuition.** A huge background term ($\lambda=h^{-30n}$) is precisely what linearizes the problem at the right exponent.

**Human precedent.** Bennett–Carbery–Tao proved *multilinear* Kakeya by heat-flow monotonicity, a PDE/positivity method. The missing step was the passage to the *linear* estimate. The rank-one background supplies it.

**Assessment.** This is new mathematics, not a recombination. The decisive constructions are invented, internally motivated, and explained well enough to be learned. Among the three results reviewed so far (7/8 quasi-RH, ETH, this one), it is the clearest instance of AI producing a *method*.
