# 7. Recyclable toolkit: what in the local smoothing / Kakeya paper is reusable, and how far it sharpens

This file extracts the parts of *Local smoothing for the Euclidean wave equation* that are useful beyond the paper. Each is stated at the greatest generality its proof supports. Status labels:
- **P** — proved in the paper;
- **L** — covered by OpenAI's formal proof, which I built ([06](06_LEAN_AUDIT.md));
- **R** — proved here;
- **X** — checked by a script in [`scripts/`](scripts/);
- **H** — heuristic, clearly marked.

| # | Tool | General form | Status | What is new here |
|---|---|---|---|---|
| K1 | Positive conserved phase-space current | any constant-coefficient real symmetric hyperbolic system | P, X, L | Isolated as a standalone object: an *exactly conserved, positive* Husimi-type density |
| K2 | Determinant energy balance | abstract hypotheses (A), (B), (W), (D) | P, R | Abstract theorem separating geometry from the wave |
| K3 | Density-to-depth | background of any rank $r$ | R | **New:** rank-$r$ version; heuristic exponent-selection rule $p=2(f+1)/f$ |
| K4 | Capacity-dimension principle | auxiliary time of any dimension $k$ | R, X | **New:** only $k=2$ gives logarithmic loss, so two auxiliary flows are *necessary* |
| K5 | Two-circle propagation | supersolutions with exponentially small potential | P, R, X, L | Sharpness constant identified exactly |
| K6 | 3-slot cofactor conservation | real $C^3$ functions, 2/3/4 slots | P, X, L | |
| K7 | The positivity engine | axioms (G1)–(G2), (A1)–(A4) | R, X | **Now a theorem:** [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Theorem 11.1, with pencils (Theorem 11.2) and the wave as the arrow pencil (Corollary 11.4) |

---

## K1. A positive, exactly conserved phase-space density

**Theorem K1** (Lemma 2.1, Prop 2.2, App. E). Let $\mathscr P=\sum_\mu C^\mu\partial_\mu$ have constant real symmetric coefficients. Let $f$ be a finite superposition of holomorphic plane-wave solutions, and $q^\mu=\frac12\sum_{\ell=1,2}w_\ell^tC^\mu w_\ell$.

Then $B=\operatorname{Re}\sum_\ell\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ is positive semidefinite with zero column divergence on $(\mathsf X,w_1,w_2)$. Here $\mathcal Z_\ell=\sum_aG_aD_{\ell a}$, $G=\overline{f(\mathsf X-\mathrm iq)}$, and $D_{\ell a}=\partial_{w_{\ell a}}+\mathrm i\sum_\mu(C^\mu w_\ell)_a\partial_\mu$.

For the wave, after the holomorphic flux map its time component is $|\mathsf S_0u|^2+|\mathsf S_1u|^2$: the wave seen through a Gaussian direction window $e^{-\frac{\kappa y}2|W-v|^2}$.

**Why it is worth recycling.** Husimi densities are the standard positive phase-space densities, but they are not conserved; Wigner functions are conserved but not positive. K1 is a lifted density that is **both positive and exactly conserved**, for every constant-coefficient symmetric system (wave, Maxwell, Dirac, and the pencils of App. E). Any argument that wants "positivity + conservation ⇒ integrability" (compensated integrability, entropy or Lyapunov arguments, phase-space concentration) can start from it.

(Verified: [`scripts/verify_wave_current.py`](scripts/verify_wave_current.py), $n=2,3$, with a negative control.)

---

## K2. Determinant energy balance (abstract engine)

The precise statement and proof are [05](05_DISTILLATION.md) Theorem II.

**Hypotheses:**
- **(A)** a positive conserved 3-slot Hermitian field $A$ with $\mu=A_{00}>0$;
- **(B)** a background $B_0$ with $B=H(B_0)\ge0$, $\det B=0$, $A:B=0$, and $B_{\rm fib}$ of rank one;
- **(W)** on the fiber wall, a common null vector of $H(r)$ and $B$ in $G(r)^\perp$;
- **(D)** fiberwise Dirichlet solutions $\det H_{\rm fib}(B_0+u)=\lambda\mu$.

**Conclusion:** $\Delta_z\mathcal F\le$ wall flux, and $\frac{n-1}{n+1}X_0\le\mathcal F\le\frac n{n+1}X_0$.

**What recycles.** None of (A)–(D) mentions waves. The theorem is a general mechanism: *positive conserved current + rank-one background + designed fibers ⇒ an almost-superharmonic energy in the parameters*. The wave enters only through the size of the wall flux and the interior density (K3).

---

## K3. Density-to-depth, with a background of any rank (new generalization)

**Lemma K3.** Let $B_{\rm fib}\ge0$ be an $n\times n$ quaternionic Hermitian matrix field of rank $r<n$, with nonzero eigenvalues in $[c\lambda,C\lambda]$. Let $u\le0$ solve $\det(B_{\rm fib}+H_{\rm fib}(u))=\lambda^r\mu$ on a 3-slot fiber, with $u=0$ on the wall. Put $f=n-r$.

1. *(Interior.)* Suppose $\mu\ge m$ on a ball $B_\rho(x)$ inside the fiber and $m^{1/f}\le\lambda$. Then $-u(x)\ge\epsilon\rho^2m^{1/f}$, with $\epsilon=\epsilon(n,r,C)>0$.
2. *(Barrier.)* Suppose $H_{\rm fib}(\mathrm r)\ge\kappa I$, $G G^*\le C'H_{\rm fib}(\mathrm r)$ for the defining function $\mathrm r$, $P>0$ with $|\nabla^jP|\le\epsilon_0\kappa^{j/2}P$, and $\mu\le(c_1\kappa P)^f$. Then $0\le u/\mathrm r\le C_2P$.

*Proof.*
1. Let $v=\epsilon m^{1/f}(|y-x|^2-\rho^2)$. Then $H_{\rm fib}(v)=6\epsilon m^{1/f}I=:\kappa I$, and
$$\det(B_{\rm fib}+\kappa I)=\prod_{i\le r}(\beta_i+\kappa)\,\kappa^f\le(2C\lambda)^r(6\epsilon)^fm\le\lambda^rm,$$
for $\epsilon$ small (using $\kappa\le C\lambda$). Comparison on the ball gives $u\le v$.
2. $H(\mathrm rP)\ge(1-C\epsilon_0)P\,H(\mathrm r)$, so $\det(B_{\rm fib}+H(C_2\mathrm rP))\ge(c\lambda)^r(C_2\kappa P/2)^f\ge\lambda^r\mu$. Comparison gives $u\ge C_2\mathrm rP$. ∎

For $r=1$ this is [05](05_DISTILLATION.md) Theorem III, i.e. §5 of the paper.

**Exponent-selection rule (H).** Under Lemma K3, $(-u)\mu\approx\rho^2\mu^{1+1/f}$, so the method measures $\|\text{wave}\|_{L^p}$ with
$$p=2\Big(1+\frac1f\Big)=\frac{2(f+1)}f,\qquad f=\#\text{free fiber directions}.$$
- $f=n-1$ gives $p_c$.
- In Prop 8.1's focusing computation, replacing the $d$ curved directions by $f$ makes the critical exponent ($f-fp/2=-1$) exactly $2(f+1)/f$.

So the rule matches the numerology of characteristic varieties with $f$ curved directions. Whether a rank-$r$ background is *available* depends on the current being tangent to $r$ independent passive coordinates. That requires more degenerate symbols, and the analogue of Lemma 4.3 (energy comparability) is unproved for $r>1$. This is a research direction, not a theorem.

---

## K4. The capacity-dimension principle (new): why there must be exactly two auxiliary flows

The auxiliary time variable $\tau$ lives in $\mathbb C_{\mathbf k}\cong\mathbb R^2$ because there are two auxiliary flows. Suppose it lived in $\mathbb R^k$.

**Proposition K4.** For $k\ge1$ and $0<b<1/4$, consider $E\ge0$ on $\mathbb R_T\times\{b\le|\tau|\le3\}$, $\tau\in\mathbb R^k$, with $\Delta_{T,\tau}E\le0$ and spherical means $H(T,s)$.

1. **$k=2$.** $H(0,2b)\ge c\,\log^{-2}(3/b)\int_1^2H(T,b)\,dT$ (Lemma 6.1 with $V=0$). The order is attained: the example below has $E(0,2b)\log^2(3/b)\to\frac94\log2$.
2. **$k\ne2$.** There are such $E$ with $E(T,b)=1-\cos T$ and
$$E(0,2b)\le C_k\,b\log(3/b).$$
More precisely, the bound is $\asymp b$ for $k=1,3$, $\asymp b^2\log(1/b)$ for $k=4$, and $\asymp b^2$ for $k\ge5$.

*Proof of 2.* Let:
- $p$ be the radial harmonic function in $\mathbb R^k$ with $p(b)=1$, $p(3)=0$;
- $g$ be radial with $\Delta g=p$ and $g(b)=g'(b)=0$, so $g$ is increasing and subharmonic;
- $L=g(3)(1-p)$, the harmonic interpolant of $g$;
- $E=(1-\cos T)p+L-g$.

Then:
- $\Delta E=(\cos T-1)p\le0$;
- $E\ge0$, because $g\le L$ by the maximum principle;
- $E(T,b)=1-\cos T$;
- $E(0,2b)\le g(3)(1-p(2b))$.

For $k=1$, $1-p(2b)=b/(3-b)$. For $k\ge3$, $p(u)\le(b/u)^{k-2}$ gives $s^{k-1}g'(s)\le b^{k-2}s^2/2$, so $g(3)\le\frac{b^{k-2}}2\int_b^3s^{3-k}ds$. ∎

Verified exactly: the Laplacian and trace identities are symbolic, nonnegativity is checked on a grid in 50-digit arithmetic, and the fitted decay exponents are $1.000,\ 0.10\,(\text{i.e. logarithmic}),\ 1.000,\ 1.944,\ 2.000$ for $k=1,\dots,5$ ([`scripts/verify_capacity_dimension.py`](scripts/verify_capacity_dimension.py), [`results/capacity_dimension.json`](results/capacity_dimension.json)).

**Consequence.** In the paper $b\asymp h=N^{-1}$. For $k\ne2$, any argument that uses only "nonnegative supersolution compared on two spheres" loses a **power** $N^{-c}$, which destroys the $\varepsilon$-loss theorem. Two auxiliary flows are therefore forced twice over:
- they make the Gaussian window isotropic, (2.4);
- they put the auxiliary time in the unique dimension of **logarithmic capacity**.

Two flows produce quaternions. So the quaternionic structure, the $\log$ (not power) loss, and the planar capacity of a point are one fact seen three ways.

---

## K5. Two-circle propagation with potential (Lemma 6.1)

**Statement** ([05](05_DISTILLATION.md) Theorem IV). Let $E\ge0$ satisfy $\Delta E\le V(|\tau|)E$ on $\mathbb R\times\{b\le|\tau|\le3\}\subset\mathbb R^3$, with $\int_b^3s\log(s/b)V\le M$ and $\sup_{[1/2,3/2]}V\le M$. Then
$$H(0,r)\ge c(M)\,\frac{\log(r/b)}{\log^2(3/b)}\int_1^2H(T,b)\,dT .$$

**Sharp constant.** With $V=0$ and $r=2b$, Remark 6.2's example gives
$$E(0,2b)\log^2(3/b)\to\tfrac94\log2\approx1.5596,$$
which the K4 script confirms to 4 digits. The proof's Gronwall factor $\exp(-\int\min(z,A-z)P)$ is verified numerically for $V=Cs^{-2}e^{-cs/h}$ down to $h=10^{-6}$.

**Formal variant (L).** OpenAI's Lean proof does not use the potential version. It proves a $V=0$ cylinder lemma with cubic logarithms (`Concentration/Cylinder/Estimate.lean`) and applies it to a *corrected* energy (`Concentration/CorrectedEnergy.lean`, correction $16-z_i^2-z_j^2$, Laplacian $-4$). This is a second, cruder route, worth comparing in a textbook treatment.

---

## K6. 3-slot cofactor conservation

**Update (28 Oct).** Generalized in [12](12_QUATERNIONIC_SLOT_CALCULUS.md) to arbitrary quaternion frames (Theorem 12.1), with the variational structure of the Monge–Ampère energy (Theorem 12.3).

For every real $C^3$ function $u$, in 2, 3 or 4 slots, $\sum_{a,\alpha}\bar e_\alpha\partial_{a\alpha}C(H(u))_{ab}=0$, with no definiteness assumed (Lemma 3.3). This is the quaternionic null-Lagrangian identity that makes every integration by parts in §4 possible.
- Exact symbolic checks: [`scripts/verify_quaternion_calculus.py`](scripts/verify_quaternion_calculus.py).
- Formal statement: `OAI/Analysis/Calculus/Quaternion/CofactorDivergence.lean`, `mooreCofactor_quaternionHessian_divergence`, for arbitrary slot systems.

---

## K7. The positivity engine: a master blueprint

**Update (28 Oct).** This blueprint is now a theorem: [11](11_FRAMEWORK_POSITIVITY_ENGINE.md) Theorem 11.1, with precise axioms and a hypothesis-usage matrix (§7). H1 becomes (A1)–(A2) and (G1)–(G2). H2 becomes (G2). H3 and H4 become (A3) and (A4). H5 is built into the slot structure. Symmetric pencils satisfy the axioms (Theorem 11.2), and the wave is the arrow pencil (Corollary 11.4).

**Hypotheses for a dispersive family $e^{-iT\phi(D)}$ on a frequency cap:**
- **H1** (current) A constant-coefficient symmetric system realizes $\phi$ as a simple extremal eigenvalue sheet, with a holomorphic flux map making the lifted current (K1) quaternionic rank one and tangent to a passive coordinate $\Phi$.
- **H2** (fibers) Fiber domains satisfying (W) of K2, with $H_{\rm fib}(r)\gtrsim s^{-2}$ and volume $\asymp a^2s^{2d}$.
- **H3** (window) Gaussian decay $e^{-c_0s/h}$ of the direction window on the fiber wall.
- **H4** (inversion) Inverting the time-component multiplier near the fiber center costs at most $e^{C\gamma^2s/h}$ with $C\gamma^2<c_0/2$.
- **H5** (capacity) The auxiliary time is two-dimensional (K4).

**Conclusion:** the critical estimate $\int_1^2\|e^{-iT\phi(D)}f\|_{p}^{p}\lesssim\log^2(2/h)\|f\|_p^p$ on the cap, with $p=2(f+1)/f$ ($f=d$ under rank-one background). Assembly: K1 + K2 + K3 + K5 + §7.

**Instances:**
- the Euclidean wave (§2–7, formally verified);
- extremal sheets of symmetric pencils with definite transverse Hessian (App. E, including the irreducible cubic cone of Ex. E.2);
- the paraboloid and Schrödinger, by rescaling caps (§12).

**Open instances (H):**
1. **Pencils whose size differs from the spatial dimension.** App. E uses "state dimension = spatial dimension" to make the polarization a local angular coordinate. Removing that restriction would bring in determinantal representations of hyperbolic polynomials (Helton–Vinnikov in three variables). It would then ask which hyperbolicity cones satisfy critical local smoothing by this route.
2. **Variable coefficients.** The Minicozzi–Sogge counterexamples show the conclusion can fail on curved metrics, so some hypothesis must fail there. H1, which needs exact holomorphic extension and a polynomial flux map, is the natural candidate. Pinpointing the failure would give a structural explanation of those counterexamples.
3. **One logarithm.** H5 plus scalar supersolutions force two logarithms (K4–K5). Reaching the lower bound $1/p_c$ needs information about the energy beyond superharmonicity.
