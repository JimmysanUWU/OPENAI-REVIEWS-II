# Proof map: every result, in dependency order, with a proof sketch

This file is the **readable form of the theorem DAG** in [`roadmap.json`](roadmap.json). Read it top to bottom along a path: each step states one result, sketches its proof from the earlier steps, and links to the full proof and to the script that checks it.

**Guarantees.** All sections below between `GENERATED` markers are rendered from `roadmap.json` by `python3 tools/check_roadmap.py --write`. The checker refuses to render unless:
- every path lists its nodes in dependency order: each node's dependencies, for its main proof or for one listed alternative, appear earlier on the path or are declared as assumed;
- every node lies on at least one path;
- every file and every structured label (§, Lemma, Theorem, …) cited here exists in the repository;
- nothing proved in prose rests on something only partially verified or heuristic.

**Status words.**
- `FORMAL`: the exact statement is a theorem of OpenAI's Lean proof, built and audited here.
- `PROVED`: a complete prose proof is in this repository.
- `PROVED_CITED`: complete except for named cited theorems.
- `CHECKED`: an exact or high-precision script checks its identities.
- `PARTIAL`: exponents and identities checked, not every estimate rewritten.
- `NEW`: not in the source papers.

**Prose chain.**
- `closed`: every result it rests on is proved here with textbook citations only.
- `cites X`: it inherits the named cited theorems of `X`.

**Terminal use.**
- `python3 tools/check_roadmap.py --path P.KAK.PLANAR` prints a path.
- `--node KAK.T12` prints a node and everything beneath it.

## Paths

<!-- BEGIN GENERATED: path-index -->
| Path | Read it for | Steps | Assumes |
|---|---|---|---|
| [Shortest complete proof: sharp L^4 local smoothing in R^{2+1}](#p-kak-planar) | Read this first. Every idea of the paper with only linear PDE; the prose chain is closed (textbook citations only). | 6 | — |
| [The distilled route to the critical estimate in every dimension](#p-kak-engine) | Axioms -> engine theorem -> pencils -> the wave as the arrow pencil -> Theorems 1.2 and 1.1. | 13 | — |
| [The paper's own route (08)](#p-kak-paper) | The wave-specific proof as written in the paper, made gapless in 08. | 11 | — |
| [Consequences: restriction, Kakeya, Nikodym, Bochner-Riesz, Stein](#p-kak-conseq) | Everything downstream of Theorems 1.2 and 1.1. | 8 | `KAK.T12`, `KAK.T11` |
| [Sharpness](#p-kak-limits) | Why the losses are what they are: the focusing lower bound log^{1/p} and the weak-norm obstruction. The angular refinement that beats both is on P.KAK.ANGULAR. | 2 | `KAK.T12` |
| [The angular endpoint (App. B, made gapless in 14)](#p-kak-angular) | How a little angular regularity removes one logarithm (strong norm) or all of them (weak norm): uniform Hankel -> shell transplant -> four-kernel envelope -> one annulus -> summation. | 5 | `KAK.T12` |
| [ETH is false (Part I)](#p-eth-i) | The shortest complete refutation of ETH: deformation counting with swap interpolation. | 7 | — |
| [The formally verified exponent rho(n) (Part II)](#p-eth-ii) | Reduction to sparse F_2 systems; parity lift; the bound OpenAI proved in Lean. | 5 | `ETH.I.SPARS`, `ETH.I.DEFCOUNT`, `ETH.I.SWAP`, `ETH.I.TRANSFER` |
| [The headline exponent (Part III, App. D)](#p-eth-iii) | Boolean tensor contraction by a 2-adic Cauchy-kernel deformation. | 9 | `ETH.I.SPARS`, `ETH.II.LIFT`, `ETH.II.TANNER`, `ETH.II.GREEDY` |
| [Corollaries, new theorems, and the method's limits](#p-eth-more) | What the counting algorithms give, the F_q extension, and the cost floor. | 9 | `ETH.I.MAIN`, `ETH.II.MAIN`, `ETH.III.FIXED`, `ETH.II.LIFT`, `ETH.II.ZEROS`, `ETH.I.DEFCOUNT`, `ETH.I.SWAP`, `ETH.I.TRANSFER` |
| [One deformation theorem for Parts I-III (13)](#p-eth-unified) | The common skeleton of swap interpolation and App. D: readout systems, the filtered node lemma, Hermite recovery with loss, correction counts; both instances; both square-root floors. | 9 | `ETH.I.SWAP`, `ETH.III.CAUCHY`, `ETH.III.PACK` |
<!-- END GENERATED: path-index -->

## Kakeya / local smoothing

Notation used throughout the Kakeya paths:

<!-- BEGIN GENERATED: glossary-kak -->
| Term | Meaning | Defined in |
|---|---|---|
| p_c, d, q | n &gt;= 2 spatial dimensions, d = n-1, p_c = 2n/d (critical exponent), q = p_c/2 = n/d. | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §0 |
| U(T), s_+(p) | Half-wave e^{-iT\|D\|}; s_+(p) = max(0, d/2 - n/p), the conjectured derivative loss. | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §0 |
| slot, frame | A slot is a quaternionic coordinate sum_alpha e_alpha z_{a alpha}. The paper uses units (1,i,j) (3-slot) or (1,i,j,k) (4-slot); 12 allows any list of quaternions (a frame). | [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §1 |
| G(u), H(u) | Quaternionic gradient G_a = sum e_alpha d_{a alpha}u and Hessian H_ab = sum e_alpha conj(e_beta) d d u. H(xi xi^T) = g g^\*. | [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §1 |
| Moore det, C(A), A:B | Determinant of a quaternionic Hermitian matrix (product of real eigenvalues); its cofactor C(A), defined by C(A):B = d/dt det(A+tB); A:B = Re tr(AB). | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §2 |
| conserved field | A Hermitian field A whose real tensor A^R_{a alpha,b beta} = Re(conj(e_alpha) A_ab e_beta) has divergence-free columns. Then A:H(u) = div(A^R grad u). C(H(u)) is always conserved. | [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §1 |
| lifted current | B = Re sum_{l=1,2} Z_l (x) conj(Z_l) on (X, w_1, w_2), built from a constant-coefficient symmetric system and its complexified solution. It is positive and conserved. | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3 |
| flux map Theta, passive coordinate Phi | Theta(X,z,eta) = (X, z^T C^mu z/2 + c^mu eta) pushes the lifted current to the target slots. Phi is the coordinate with Phi o Theta = eta; tangency means (dPhi) chi = 0. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6 |
| chi, A^0, mu | Pushed current A^0 = chi chi^\* (quaternionic rank one); regularized A = g(Q)(J_h \*_X A^0) (+ h^{50n}C(H(r)) in n &gt;= 3); the density mu = A_TT. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §3 |
| polarization v, window | v(xi): eigenvector coordinate, E(xi) = (1, v(xi)). Rows of chi carry the Gaussian window exp(-(y/2h)(W-v)^\* S_perp (W-v)), S = phi I - C &gt;= 0. For the wave, v is stereographic projection. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6 |
| fibre, chart (W,Q), wall | At parameter z = (T,tau), s = \|tau\|, the fibre is T^n x {(W,Q): \|W\|^2/w_0^2 + 100\|Q\|^2 &lt; 1}, Q = Phi/a, a = Kh. The wall is r = 0. In 13 (n = 2) it is a product: disk in W times a disk in Q. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §1 |
| background B_0, lambda | B_0 = lambda\|Phi\|^2 with lambda = h^{-Lambda}. B = H(B_0) has rank one; its fibre block removes one direction, leaving d free directions. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §3 |
| energy F, X_0 | X_0 = int(-u)mu over the fibre; F = X_0 + lambda^{-1} int u int_0^1 det D_{theta u}, comparable to X_0. In 13, F = X_0/2. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §3 |
| cap, gamma, sigma | Frequencies with h xi in a compact annular cap K and \|v(h xi)\| &lt;= gamma; sigma = s/h. The energy is read on the circles \|tau\| = 4a (data, T = 0) and \|tau\| = 2a (observation). | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §4 |
| two-circle propagation | Delta_{T,tau}E &lt;= V(\|tau\|)E implies H(0,r) &gt;~ log(r/b)/log^2(3/b) int_1^2 H(T,b), where H is the circular mean. This is the source of the log^2 loss; only a 2-dim auxiliary time gives logarithms (K4). | [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8 |
| engine data (G1)-(G2), (A1)-(A4) | The axioms of Theorem 11.1: chart and scales; positive conserved current; tangency; Gaussian rows; invertible time row. | [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §1 |
| arrow pencil, Dirac pencil | A(xi) = [[xi_n, xi'^T],[xi', -xi_n I]], with lambda_max = \|xi\|; it is the paper's wave system. For n = 2 it is the Dirac pencil [[xi_2, xi_1],[xi_1, -xi_2]], whose chart is explicit: Phi = P_2 - sqrt(tau^2 - P_1^2). | [kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md) §1 |
| Lambda_omega, L(v) | Lambda_omega = (1 - Delta_sphere)^{1/2} on each sphere about 0 (the scalar sqrt(1+nu(nu+d-1)) on degree-nu harmonics); L(v) = log(2+v). | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) §0.1 |
| angular degree &lt;= M | Restricted to almost every sphere about 0, a sum of spherical harmonics of degree &lt;= M. | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) §0.1 |
| uniform Hankel expansion | S_z F = two point evaluations (theta, -theta) plus an L^1 -&gt; L^inf error, valid on degree &lt;= M once z &gt;= (2+M)^K; explicit via Hankel's integral. | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 1.2 |
| shell transplant | Seen from a ball of R wavelengths, an input shell at radius Q &gt;&gt; R^3 gives a wave independent of Q, so it can be moved to radius ~R^{3+n/p}: observation on a ball costs L(R)^A. | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 3.2 |
| radial profile h_g, G_g | h_g(v,theta) = \|v\|^{d/p}(\|g(\|v\|theta)\| + \|g(-\|v\|theta)\| + sphere average); G_g its maximal function in v; \|\|G_g\|\|_{L^p(dv dtheta)} &lt;~ \|\|g\|\|_p. | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Definition 4.1 |
| angular bands, selection m(j,r) | M_m = 2^{2^m}, Q_m = S_m - S_{m-1}; m(j,r) = the largest band whose stationary-phase threshold H_m fits below min(N_j^eta, N_j r). | [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) §6.3 |
<!-- END GENERATED: glossary-kak -->

### Shortest complete proof: sharp L^4 local smoothing in R^{2+1}
<a id="p-kak-planar"></a>

<!-- BEGIN GENERATED: path-p-kak-planar -->
*Read this first. Every idea of the paper with only linear PDE; the prose chain is closed (textbook citations only).*

#### <a id="p-kak-planar--kak-q-calc"></a>1. `KAK.Q.CALC` — Quaternionic slot calculus

*PROVED, CHECKED* · prose chain: closed

**Statement.** For the 3-slot operators G, H: Hessian calculus (H(sum|F|^2)=4J\*J), Moore determinant and cofactor identities, cofactor divergence (sum ebar_alpha d_{a alpha} C(H(u))_{ab}=0), Green's identity and compression.

**Proof sketch.** H(u) is the image of D^2u under the real-linear map that sends xi xi^T to g g^\*, g_a = sum_alpha e_alpha xi_{a alpha}; this gives Hermitian, positivity and the product rule. The Moore determinant is the product of the real eigenvalues (quaternionic spectral theorem); cofactor identities follow from the 4m x 4m real representation, whose determinant is det^4. Cofactor divergence: on exponentials e^{xi.z} the divergence of the polarized cofactor is a multiple of (sum g_j)^\* Q(g_1g_1^\*,...) = 0, because the cofactor of a rank < m PSD matrix kills its image; polynomial identities then extend this to all u.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §2, Lemmas 2.1-2.4; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K6; [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) Lemma 1.2, Lemma 2.1, Prop. 12.5 · **checks:** [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py), [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-planar--kak-frame"></a>2. `KAK.FRAME` — Cofactor divergence for arbitrary quaternionic frames (Theorem 12.1)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-planar--kak-q-calc)

**Statement.** For any quaternions f_{a,i} (any number per slot, possibly dependent) and H_f(u)_{ab}=sum f_{a,i} conj(f_{b,j}) d_{a,i} d_{b,j} u: sum_{a,i} conj(f_{a,i}) d_{a,i} C(H_f(u))_{ab}=0. Every frame calculus is the restriction H_f(u) o A = H^{(4)}(u o A) of the 4-slot calculus (Lemma 12.2). Q^{-1}:H_f is elliptic iff every slot frame is R-linearly independent, hence at most 4 real coordinates per slot (Prop. 12.4).

**Proof sketch.** The exponential proof of cofactor divergence uses only the rank-one map xi xi^T -> g g^\* and the cofactor annihilating images; both hold for any quaternion frame. Alternatively, with A_{a,i alpha} the components of f_{a,i}, H_f(u) o A = H^{(4)}(u o A) and the frame divergence pulls back to the 4-slot divergence. Ellipticity: the symbol g^\*Q^{-1}g vanishes iff g(xi) = sum f_i xi_i = 0, which has nonzero solutions iff a slot frame is dependent.

**Full proof:** [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §3-5, Thm 12.1, Lemma 12.2, Prop. 12.4 · **checks:** [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-planar--kak-c-lift"></a>3. `KAK.C.LIFT` — Lifted current and push-forward

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-planar--kak-q-calc)

**Statement.** For constant-coefficient symmetric systems, B=Re sum_l Z_l (x) conj(Z_l) is PSD and conserved on (X,w_1,w_2); its density push-forward by a holomorphic flux map is conserved and quaternionic rank one.

**Proof sketch.** The fields D_{la} = d/dw_{la} + i V_{la} annihilate every holomorphic function of X - iq, so Z_l = sum G_a D_{la} is divergence free (it reduces to w_l^T P G = 0). Its self-acceleration is purely imaginary, so B = Re sum Z_l (x) conj(Z_l) is PSD and conserved. A holomorphic flux map Theta with (d^2_{w1} + d^2_{w2}) Theta = 0 pushes conserved densities to conserved densities. The Gram rows of the pushed tensor are (Q, -iQ, -jQ), which makes it a quaternionic rank-one chi chi^\*.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Lemmas 3.1-3.2; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm I · **checks:** [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py)

#### <a id="p-kak-planar--kak-pencil"></a>4. `KAK.PENCIL` — Symmetric pencils carry engine data (Theorem 11.2)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`KAK.C.LIFT`](#p-kak-planar--kak-c-lift), [`KAK.Q.CALC`](#p-kak-planar--kak-q-calc)

**Statement.** C(xi)=sum xi_j C_j real symmetric n x n, phi=lambda_max simple at omega with D^2 phi(omega)>0 on omega^perp: with E(xi)=(1,v(xi)), S=phi I-C, L w=(C(w)e_0)_perp, B(W)=e(W)^T C e(W)/(1+W.W), the chart P=tau B(W)+Phi omega and App. E's pushed current chi_mu=J^{-1/2} sum_a (C^mu z)_a i G_a satisfy (G1)-(G2) and (A1)-(A4); D^2 phi(omega)=2 L^T S_perp^{-1} L, Dv(omega)=S_perp^{-1} L, dB(0)=2 L^T, and the window is exp(-(y/2h)(W-v)^\* S_perp (W-v)).

**Proof sketch.** Perturbation theory for the simple top eigenvalue gives a smooth polarization E = (1,v), S = phi I - C >= 0 with a spectral gap, D^2 phi = 2L^T S_perp^{-1} L and Dv = S_perp^{-1}L. Definiteness makes L an isomorphism on omega^perp, so dB(0) = 2L^T is invertible and the chart P = tau B(W) + Phi omega is a biholomorphism with the required scales. The lifted current pushed by Theta(z,eta) = (z^T C^mu z/2) + (0,omega)eta has rows c_mu(W,xi)|xi|^{d/2}y^{-d/2}e^{-Psi} with Psi = (y/2)(W-v)^\* S_perp (W-v) >= c sigma |W-v|^2: the Gaussian window. The time row 1 + W.v is invertible near the centre.

**Full proof:** [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Lemmas 6.1-6.4, Thm 11.2 · **checks:** [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py) · *Note:* Conservation of the pencil current checked directly in target coordinates for 5 pencils (incl. the cubic of Ex. E.2), to 1e-38.

#### <a id="p-kak-planar--kak-prop"></a>5. `KAK.PROP` — Two-circle propagation

*PROVED, CHECKED* · prose chain: closed

**Statement.** Delta E <= V(|tau|) E on R x {b<=|tau|<=3}: H(0,r) >= c log(r/b)/log^2(3/b) \* int_1^2 H(T,b) dT; log^{-2} sharp for scalar supersolutions.

**Proof sketch.** Average over rotations. Radial supersolution comparisons (Gronwall in x = log(s/b), using int s log(s/b) V < inf) carry the sin-weighted T-average outward from |tau| = b, losing 1/D. Weak Harnack on a fixed region chains to |tau| = 1. A return passage with the cos(pi T/2) radial solution brings it to |tau| = r, losing log(r/b)/log(1/b).

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8, Lemma 8.1, Prop. 8.2; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K5 · **checks:** [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py) · *Note:* Not used by the Lean proof, which takes a cruder route (corrected energy, power 2n+5).

#### <a id="p-kak-planar--kak-app-d"></a>6. `KAK.APP.D` — Planar case by a linear Dirichlet problem (Theorem 13.1)

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.PENCIL`](#p-kak-planar--kak-pencil), [`KAK.FRAME`](#p-kak-planar--kak-frame), [`KAK.PROP`](#p-kak-planar--kak-prop)

**Statement.** n=2: ||e^{-iT|D|} phi(D/N) f||_{L^4(R^2 x [1,2])} <= C log(2+N)^{1/2} ||f||_4, hence local smoothing in R^{2+1} for all 2<=p&lt;inf with eps loss. The fibre problem is linear: C(B_fib):H_fib(u)=mu, i.e. 4 sum Re(conj e_a e_b) d_{X_a X_b} u + 4|J|^2 Lap_W u = mu on a disk, solved by Fourier series and the Dirichlet Laplacian; the energy (1/2) int(-u)mu is controlled by a completed-square inequality A:H >= D_B(H). No Monge-Ampere theory.

**Proof sketch.** For n = 2 the lambda-coefficient of det(lambda B_fib + H_fib) is C(B_fib):H_fib, which is linear. With the Dirac pencil chart it equals 4 sum Re(conj(e_a)e_b) d_{X_aX_b} + 4|J|^2 Lap_W, solved mode by mode by the Dirichlet Laplacian on a disk. The energy (1/2)int(-u)mu: Green's identity with the conserved mixed cofactor gives -Delta X_0 = 2 int D_B, and with A it gives -Delta X_0 = int A:H - flux. The completed square A:H >= D_B then gives Delta F <= flux. Window estimates and propagation follow as in KAK.ENGINE, with linear comparisons.

**Full proof:** [kakeya/13](reviews/claude/kakeya/13_PLANAR_ENTRY_POINT.md) Thm 13.1, Cor. 13.2, §1-5; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 · **checks:** [verify_planar_entry](reviews/claude/kakeya/scripts/verify_planar_entry.py), [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py), [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) · *Note:* Uses the 2D Dirac pencil, whose chart is explicit: W is the stereographic coordinate of (P - Phi e_2)/tau on the complexified unit circle and Phi = P_2 - sqrt(tau^2 - P_1^2). For n=2 the eps-loss form is due to Guth-Wang-Zhang (2020) by a different method.
<!-- END GENERATED: path-p-kak-planar -->

### The distilled route to the critical estimate in every dimension
<a id="p-kak-engine"></a>

<!-- BEGIN GENERATED: path-p-kak-engine -->
*Axioms -> engine theorem -> pencils -> the wave as the arrow pencil -> Theorems 1.2 and 1.1.*

#### <a id="p-kak-engine--kak-q-calc"></a>1. `KAK.Q.CALC` — Quaternionic slot calculus

*PROVED, CHECKED* · prose chain: closed

**Statement.** For the 3-slot operators G, H: Hessian calculus (H(sum|F|^2)=4J\*J), Moore determinant and cofactor identities, cofactor divergence (sum ebar_alpha d_{a alpha} C(H(u))_{ab}=0), Green's identity and compression.

**Proof sketch.** H(u) is the image of D^2u under the real-linear map that sends xi xi^T to g g^\*, g_a = sum_alpha e_alpha xi_{a alpha}; this gives Hermitian, positivity and the product rule. The Moore determinant is the product of the real eigenvalues (quaternionic spectral theorem); cofactor identities follow from the 4m x 4m real representation, whose determinant is det^4. Cofactor divergence: on exponentials e^{xi.z} the divergence of the polarized cofactor is a multiple of (sum g_j)^\* Q(g_1g_1^\*,...) = 0, because the cofactor of a rank < m PSD matrix kills its image; polynomial identities then extend this to all u.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §2, Lemmas 2.1-2.4; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K6; [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) Lemma 1.2, Lemma 2.1, Prop. 12.5 · **checks:** [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py), [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-engine--kak-frame"></a>2. `KAK.FRAME` — Cofactor divergence for arbitrary quaternionic frames (Theorem 12.1)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-engine--kak-q-calc)

**Statement.** For any quaternions f_{a,i} (any number per slot, possibly dependent) and H_f(u)_{ab}=sum f_{a,i} conj(f_{b,j}) d_{a,i} d_{b,j} u: sum_{a,i} conj(f_{a,i}) d_{a,i} C(H_f(u))_{ab}=0. Every frame calculus is the restriction H_f(u) o A = H^{(4)}(u o A) of the 4-slot calculus (Lemma 12.2). Q^{-1}:H_f is elliptic iff every slot frame is R-linearly independent, hence at most 4 real coordinates per slot (Prop. 12.4).

**Proof sketch.** The exponential proof of cofactor divergence uses only the rank-one map xi xi^T -> g g^\* and the cofactor annihilating images; both hold for any quaternion frame. Alternatively, with A_{a,i alpha} the components of f_{a,i}, H_f(u) o A = H^{(4)}(u o A) and the frame divergence pulls back to the 4-slot divergence. Ellipticity: the symbol g^\*Q^{-1}g vanishes iff g(xi) = sum f_i xi_i = 0, which has nonzero solutions iff a slot frame is dependent.

**Full proof:** [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §3-5, Thm 12.1, Lemma 12.2, Prop. 12.4 · **checks:** [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-engine--kak-maenergy"></a>3. `KAK.MAENERGY` — The Monge-Ampere operator is the gradient of a concave energy (Theorem 12.3)

*NEW, PROVED* · prose chain: closed · uses [`KAK.FRAME`](#p-kak-engine--kak-frame)

**Statement.** Independent frames, H_f(B_0)>=0, P={u: u=0 on the boundary, H_f(B_0+u)>0}: E(u)=int u int_0^1 det H_f(B_0+theta u) dtheta has dE(u)[v]=int v det H_f(B_0+u) and d^2E(u)[v,v]=-int Re(G(v)^\* C(H_f(B_0+u)) G(v))<0; so det H_f(B_0+u)=g is the Euler-Lagrange equation of the strictly concave E(u)-int g u, and the energy of 08 §6 is its Legendre-type transform.

**Proof sketch.** First variation: differentiate under the integral, move H_f from v to u by Green's identity (cofactor conserved, u = v = 0 on the boundary), and integrate d/dtheta(theta det H_f(B_0 + theta u)). Second variation: int v C:H(v) = -int grad v^T C^R grad v < 0, since C(H_f(B_0+u)) > 0 and the frames are independent.

**Full proof:** [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) §6, Thm 12.3 · **checks:** [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-engine--kak-c-lift"></a>4. `KAK.C.LIFT` — Lifted current and push-forward

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-engine--kak-q-calc)

**Statement.** For constant-coefficient symmetric systems, B=Re sum_l Z_l (x) conj(Z_l) is PSD and conserved on (X,w_1,w_2); its density push-forward by a holomorphic flux map is conserved and quaternionic rank one.

**Proof sketch.** The fields D_{la} = d/dw_{la} + i V_{la} annihilate every holomorphic function of X - iq, so Z_l = sum G_a D_{la} is divergence free (it reduces to w_l^T P G = 0). Its self-acceleration is purely imaginary, so B = Re sum Z_l (x) conj(Z_l) is PSD and conserved. A holomorphic flux map Theta with (d^2_{w1} + d^2_{w2}) Theta = 0 pushes conserved densities to conserved densities. The Gram rows of the pushed tensor are (Q, -iQ, -jQ), which makes it a quaternionic rank-one chi chi^\*.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Lemmas 3.1-3.2; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm I · **checks:** [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py)

#### <a id="p-kak-engine--kak-pencil"></a>5. `KAK.PENCIL` — Symmetric pencils carry engine data (Theorem 11.2)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`KAK.C.LIFT`](#p-kak-engine--kak-c-lift), [`KAK.Q.CALC`](#p-kak-engine--kak-q-calc)

**Statement.** C(xi)=sum xi_j C_j real symmetric n x n, phi=lambda_max simple at omega with D^2 phi(omega)>0 on omega^perp: with E(xi)=(1,v(xi)), S=phi I-C, L w=(C(w)e_0)_perp, B(W)=e(W)^T C e(W)/(1+W.W), the chart P=tau B(W)+Phi omega and App. E's pushed current chi_mu=J^{-1/2} sum_a (C^mu z)_a i G_a satisfy (G1)-(G2) and (A1)-(A4); D^2 phi(omega)=2 L^T S_perp^{-1} L, Dv(omega)=S_perp^{-1} L, dB(0)=2 L^T, and the window is exp(-(y/2h)(W-v)^\* S_perp (W-v)).

**Proof sketch.** Perturbation theory for the simple top eigenvalue gives a smooth polarization E = (1,v), S = phi I - C >= 0 with a spectral gap, D^2 phi = 2L^T S_perp^{-1} L and Dv = S_perp^{-1}L. Definiteness makes L an isomorphism on omega^perp, so dB(0) = 2L^T is invertible and the chart P = tau B(W) + Phi omega is a biholomorphism with the required scales. The lifted current pushed by Theta(z,eta) = (z^T C^mu z/2) + (0,omega)eta has rows c_mu(W,xi)|xi|^{d/2}y^{-d/2}e^{-Psi} with Psi = (y/2)(W-v)^\* S_perp (W-v) >= c sigma |W-v|^2: the Gaussian window. The time row 1 + W.v is invertible near the centre.

**Full proof:** [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Lemmas 6.1-6.4, Thm 11.2 · **checks:** [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py) · *Note:* Conservation of the pencil current checked directly in target coordinates for 5 pencils (incl. the cubic of Ex. E.2), to 1e-38.

#### <a id="p-kak-engine--kak-dir"></a>6. `KAK.DIR` — Periodic Dirichlet problem (3-slot quaternionic Monge-Ampere)

*PROVED_CITED* · prose chain: cites KAK.DIR · uses [`KAK.Q.CALC`](#p-kak-engine--kak-q-calc)

**Statement.** det H(B_0+u)=g, H(B_0+u)>0, u=0 on the wall: unique smooth solution with br<=u<=b'r.

**Proof sketch.** Uniqueness and the barriers b r <= u <= b' r come from the comparison principle (concavity of log det). Existence uses the method of continuity: openness from Schauder theory for the linearization Q^{-1}:H; C^1 and interior C^2 bounds by the maximum principle; the boundary C^2 bound by CNS rigid-motion barriers in 4 slots (quaternionic covariance); C^{2,alpha} from Evans-Krylov/Krylov (cited).

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §5, Lemma 5.1, Thm 5.2 · *Note:* Comparison, barriers and the boundary C^2 estimate are proved; existence uses cited Schauder theory and the Evans-Krylov/Krylov estimates (GT 6.14, 6.19, 17.14, 17.26'). OpenAI's Lean proof contains a formal existence theory in lean_area.

#### <a id="p-kak-engine--kak-prop"></a>7. `KAK.PROP` — Two-circle propagation

*PROVED, CHECKED* · prose chain: closed

**Statement.** Delta E <= V(|tau|) E on R x {b<=|tau|<=3}: H(0,r) >= c log(r/b)/log^2(3/b) \* int_1^2 H(T,b) dT; log^{-2} sharp for scalar supersolutions.

**Proof sketch.** Average over rotations. Radial supersolution comparisons (Gronwall in x = log(s/b), using int s log(s/b) V < inf) carry the sin-weighted T-average outward from |tau| = b, losing 1/D. Weak Harnack on a fixed region chains to |tau| = 1. A return passage with the cos(pi T/2) radial solution brings it to |tau| = r, losing log(r/b)/log(1/b).

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8, Lemma 8.1, Prop. 8.2; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K5 · **checks:** [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py) · *Note:* Not used by the Lean proof, which takes a cruder route (corrected energy, power 2n+5).

#### <a id="p-kak-engine--kak-cap"></a>8. `KAK.CAP` — Capacity-dimension principle (new)

*NEW, PROVED, CHECKED* · prose chain: closed

**Statement.** Superharmonic comparison between two spheres in an auxiliary time of dimension k costs log^{-2} only for k=2; for k!=2 it loses a power.

**Proof sketch.** For k = 2 this is KAK.PROP with V = 0, and Remark 6.2's example is sharp. For k != 2 there is an explicit superharmonic E = (1 - cos T)p + L - g, built from the radial harmonic p and a subharmonic g, with E(T,b) = 1 - cos T but E(0,2b) <~ b log(3/b), a power loss.

**Full proof:** [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K4 · **checks:** [verify_capacity_dimension](reviews/claude/kakeya/scripts/verify_capacity_dimension.py)

#### <a id="p-kak-engine--kak-engine"></a>9. `KAK.ENGINE` — The positivity engine (Theorem 11.1)

*NEW, PROVED* · prose chain: cites KAK.DIR · uses [`KAK.Q.CALC`](#p-kak-engine--kak-q-calc), [`KAK.DIR`](#p-kak-engine--kak-dir), [`KAK.PROP`](#p-kak-engine--kak-prop), [`KAK.CAP`](#p-kak-engine--kak-cap)

**Statement.** Engine data on a cap (phi, K', v) -- a holomorphic fibre chart (W, Phi) with the scales (G1)-(G2), a positive conserved current chi chi^\* tangent to Phi whose rows are Fourier multipliers with Gaussian symbols C (hs)^{-d/2} sigma^m exp(-c_0 sigma |W - v(eta)|^2) (A1)-(A3), and an invertible time row (A4) -- imply, for some aperture gamma and all 0&lt;h<=1, int_1^2 ||e^{-iT phi(D)} f||_{L^p(T^n)}^p dT <= C log^2(2/h) ||f||_p^p, p=2n/(n-1), for f with Fourier support in the gamma-cap at scale h.

**Proof sketch.** Re-run 08 §4-§8 using only the axioms. (G1)-(G2) give the fibre geometry. (A1)-(A2) give the conserved, tangent regularized current with A:B = 0, so the energy balance and comparison hold verbatim. (A3) gives mu, ||A|| <~ s^{-2d} sigma^m P_\*^d, with e^{-c_0' sigma} on the wall. (A4) gives mu >= m_X at the centre. Choosing gamma with 2qC_2gamma^2 <= c_0'/2 makes Delta F/F exponentially small, and Lemma 8.1 concludes.

**Full proof:** [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §1-5, Thm 11.1, Prop. 4.1; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K7 · *Note:* Re-runs the proofs of KAK.E.BAL, KAK.E.DEPTH and KAK.FIB under the axioms; the equation enters only through (A1)-(A4). Hypothesis-usage matrix and failure modes in 11 §7.

#### <a id="p-kak-engine--kak-app-e"></a>10. `KAK.APP.E` — Curved characteristic sheets of symmetric systems

*PROVED, CHECKED* · prose chain: cites KAK.DIR · uses [`KAK.ENGINE`](#p-kak-engine--kak-engine), [`KAK.PENCIL`](#p-kak-engine--kak-pencil)

**Statement.** For a simple extremal eigenvalue phi of a symmetric pencil with definite transverse Hessian: ||e^{-iT phi(D)} b(hD) f||_{L^{p_c}} <= C log(2/h)^{2/p_c} ||f||_{p_c}.

**Proof sketch.** Theorem 11.2 (KAK.PENCIL) supplies engine data on caps around every direction where lambda_max is simple and curved, and Theorem 11.1 (KAK.ENGINE) gives the torus cap estimate. Rescale C -> C/M so that |grad phi| <= 1, which allows the cube transfer to R^n. The smallest eigenvalue follows from -C and f -> conj f(-x).

**Full proof:** [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Cor. 11.3; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 · **checks:** [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py), [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py) · *Note:* Completed in 11 (chart and current quantified, wave proof reused as one abstract argument). The single-patch lower bound of App. E is not re-derived.

#### <a id="p-kak-engine--kak-arrow"></a>11. `KAK.ARROW` — The wave is the arrow pencil (Corollary 11.4)

*NEW, PROVED, CHECKED* · prose chain: cites KAK.DIR · uses [`KAK.APP.E`](#p-kak-engine--kak-app-e)

**Statement.** A(xi)=[[xi_n, xi'^T],[xi', -xi_n I]] has lambda_max=|xi| (simple off the ray xi'=0, xi_n<0), polarization v(xi)=xi'/(|xi|+xi_n) (stereographic) and group-velocity map B(W)=(2W,1-W.W)/(1+W.W) (inverse stereographic); the paper's null-coordinate wave system is A after the frame change (t,s_0)->(T,b) and the state scaling (C^T)^{-1/2}. Hence Theorem 1.2's cap estimate is Corollary 11.3 for two rotated arrow pencils (a second derivation, with a spatial passive coordinate).

**Proof sketch.** Factor the characteristic polynomial: (t + xi_n)^{n-2}(t^2 - |xi|^2). The top eigenvector (1, xi'/(|xi|+xi_n)) is stereographic projection, and e(W)^T C_j e(W)/(1+W.W) is the inverse stereographic map. Rotating (t,s_0) -> (T,b) and scaling the states by (C^T)^{-1/2} turns 08's null-coordinate system into this pencil (exact check). Two rotated arrow pencils cover S^{n-1}.

**Full proof:** [kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) §6, Cor. 11.4 · **checks:** [verify_positivity_engine](reviews/claude/kakeya/scripts/verify_positivity_engine.py)

#### <a id="p-kak-engine--kak-t12"></a>12. `KAK.T12` — Critical annular estimate (Theorem 1.2)

*PROVED* · prose chain: cites KAK.DIR · uses [`KAK.FIB`](#p-kak-paper--kak-fib), [`KAK.PROP`](#p-kak-engine--kak-prop), [`KAK.C.WAVE`](#p-kak-paper--kak-c-wave)

**Statement.** ||U(T)P_N f||_{L^{p_c}(R^n x [1,2])} <= C [log(2+N)]^{2/p_c} ||f||_{p_c}.

**Proof sketch.** On one cap on the torus: Proposition 7.1 (KAK.FIB) gives F at |tau| = 4a <~ a^2||f||^p, F at |tau| = 2a >~ a^2||g_T||^p, and Delta F <= V F with V exponentially small. Two-circle propagation (KAK.PROP, b = 2a, r = 4a) gives int_1^2 ||g_T||^p <~ log^2(2/h)||f||^p. Finitely many caps cover the annulus; kernel decay and the cube argument transfer to R^n.

**On this path the result is reached by the alternative proof below**; the sketch above is the main proof.

**Alternative proof.** Uses [`KAK.ARROW`](#p-kak-engine--kak-arrow). Second proof: the engine theorem applied to two rotated arrow pencils (spatial passive coordinate), then the transfer of 08 §9. ([kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) Cor. 11.4)

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 · *Note:* Formalized only with the weaker power 2n+5, which suffices for every consequence. The prose chain cites the Dirichlet theory of KAK.DIR and textbook Riesz-Thorin, Mihlin, Littlewood-Paley.

#### <a id="p-kak-engine--kak-t11"></a>13. `KAK.T11` — Sogge's local smoothing conjecture (Theorem 1.1)

*FORMAL, PROVED* · prose chain: cites KAK.DIR · uses [`KAK.T12`](#p-kak-engine--kak-t12)

**Statement.** ||U(T)f||_{L^p(R^n x [1,2])} <= C ||J^{s_+(p)+eps} f||_p for n>=2, 2<=p&lt;inf, eps>0.

**Proof sketch.** Interpolate Theorem 1.2 at p_c with L^2 (Plancherel) and L^inf (kernel L^1 norm h^{-d/2}) to get h^{-s_+(p)} log^{beta(p)} on each dyadic annulus. The factor 2^{-j eps} sums the annuli; low frequencies by Mihlin.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.euclidean_local_smoothing`
<!-- END GENERATED: path-p-kak-engine -->

### The paper's own route (08)
<a id="p-kak-paper"></a>

<!-- BEGIN GENERATED: path-p-kak-paper -->
*The wave-specific proof as written in the paper, made gapless in 08.*

#### <a id="p-kak-paper--kak-q-calc"></a>1. `KAK.Q.CALC` — Quaternionic slot calculus

*PROVED, CHECKED* · prose chain: closed

**Statement.** For the 3-slot operators G, H: Hessian calculus (H(sum|F|^2)=4J\*J), Moore determinant and cofactor identities, cofactor divergence (sum ebar_alpha d_{a alpha} C(H(u))_{ab}=0), Green's identity and compression.

**Proof sketch.** H(u) is the image of D^2u under the real-linear map that sends xi xi^T to g g^\*, g_a = sum_alpha e_alpha xi_{a alpha}; this gives Hermitian, positivity and the product rule. The Moore determinant is the product of the real eigenvalues (quaternionic spectral theorem); cofactor identities follow from the 4m x 4m real representation, whose determinant is det^4. Cofactor divergence: on exponentials e^{xi.z} the divergence of the polarized cofactor is a multiple of (sum g_j)^\* Q(g_1g_1^\*,...) = 0, because the cofactor of a rank < m PSD matrix kills its image; polynomial identities then extend this to all u.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §2, Lemmas 2.1-2.4; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K6; [kakeya/12](reviews/claude/kakeya/12_QUATERNIONIC_SLOT_CALCULUS.md) Lemma 1.2, Lemma 2.1, Prop. 12.5 · **checks:** [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py), [verify_slot_calculus](reviews/claude/kakeya/scripts/verify_slot_calculus.py)

#### <a id="p-kak-paper--kak-c-lift"></a>2. `KAK.C.LIFT` — Lifted current and push-forward

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-paper--kak-q-calc)

**Statement.** For constant-coefficient symmetric systems, B=Re sum_l Z_l (x) conj(Z_l) is PSD and conserved on (X,w_1,w_2); its density push-forward by a holomorphic flux map is conserved and quaternionic rank one.

**Proof sketch.** The fields D_{la} = d/dw_{la} + i V_{la} annihilate every holomorphic function of X - iq, so Z_l = sum G_a D_{la} is divergence free (it reduces to w_l^T P G = 0). Its self-acceleration is purely imaginary, so B = Re sum Z_l (x) conj(Z_l) is PSD and conserved. A holomorphic flux map Theta with (d^2_{w1} + d^2_{w2}) Theta = 0 pushes conserved densities to conserved densities. The Gram rows of the pushed tensor are (Q, -iQ, -jQ), which makes it a quaternionic rank-one chi chi^\*.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Lemmas 3.1-3.2; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm I · **checks:** [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py)

#### <a id="p-kak-paper--kak-c-wave"></a>3. `KAK.C.WAVE` — The positive wave current

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.C.LIFT`](#p-kak-paper--kak-c-lift)

**Statement.** chi=sum V^v iota_e(conj u_{kappa,v}) m_{kappa,v}, A^0=chi chi^\*: on Y!=0, A^0 is smooth, PSD and conserved, annihilated by the passive differential (dPhi)A^0=0, and A^0_TT=|S_0 u|^2+|S_1 u|^2 with explicit cap multipliers S_0,S_1 (3.5).

**Proof sketch.** Apply KAK.C.LIFT to the wave system in null coordinates with the flux map omega_mu = z^T C^mu z/2 + passive eta. Evaluate: the Gaussian factor exp(-(kappa y/2)|W - v|^2) comes from Im q, and the row vectors V^v from M E(v). Tangency (dPhi)chi = 0 holds because Phi o Theta = eta.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §3, Prop. 3.3; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K1 · **checks:** [verify_wave_current](reviews/claude/kakeya/scripts/verify_wave_current.py), [verify_wave_current_target](reviews/claude/kakeya/scripts/verify_wave_current_target.py)

#### <a id="p-kak-paper--kak-f-geom"></a>4. `KAK.F.GEOM` — Fibre geometry

*PROVED* · prose chain: closed · uses [`KAK.Q.CALC`](#p-kak-paper--kak-q-calc)

**Statement.** On the closed fibres, uniformly in a&lt;s<4: H_fib(r) >= c s^-2 I, G_fib G_fib^\* <= H_fib, volume element ~ a^2 s^{2d}, det H_fib ~ a^-2 s^-2d; the full H(r)>=0 has a one-dimensional quaternionic kernel.

**Proof sketch.** With F = (W/w_0, 10Q) holomorphic, H(r) = 4J^\*J and G(r) = 2J^\* F i. The chart inverse bounds |J_fib^{-1}| by s, which gives H_fib >= c s^-2; |F| <= 1 gives G G^\* <= H; the Jacobian of (W,Q) -> (P,Z) gives the volume. The kernel of J is the tau-direction at fixed (W,Q), which lies in ker dPhi.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §4, Lemma 4.1

#### <a id="p-kak-paper--kak-dir"></a>5. `KAK.DIR` — Periodic Dirichlet problem (3-slot quaternionic Monge-Ampere)

*PROVED_CITED* · prose chain: cites KAK.DIR · uses [`KAK.Q.CALC`](#p-kak-paper--kak-q-calc)

**Statement.** det H(B_0+u)=g, H(B_0+u)>0, u=0 on the wall: unique smooth solution with br<=u<=b'r.

**Proof sketch.** Uniqueness and the barriers b r <= u <= b' r come from the comparison principle (concavity of log det). Existence uses the method of continuity: openness from Schauder theory for the linearization Q^{-1}:H; C^1 and interior C^2 bounds by the maximum principle; the boundary C^2 bound by CNS rigid-motion barriers in 4 slots (quaternionic covariance); C^{2,alpha} from Evans-Krylov/Krylov (cited).

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §5, Lemma 5.1, Thm 5.2 · *Note:* Comparison, barriers and the boundary C^2 estimate are proved; existence uses cited Schauder theory and the Evans-Krylov/Krylov estimates (GT 6.14, 6.19, 17.14, 17.26'). OpenAI's Lean proof contains a formal existence theory in lean_area.

#### <a id="p-kak-paper--kak-e-bal"></a>6. `KAK.E.BAL` — Determinant energy balance

*PROVED, CHECKED* · prose chain: cites KAK.DIR · uses [`KAK.Q.CALC`](#p-kak-paper--kak-q-calc), [`KAK.C.WAVE`](#p-kak-paper--kak-c-wave), [`KAK.F.GEOM`](#p-kak-paper--kak-f-geom), [`KAK.DIR`](#p-kak-paper--kak-dir)

**Statement.** Positive conserved current + rank-one background + designed fibres: Delta_z F <= wall flux and (n-1)/(n+1) X_0 <= F <= n/(n+1) X_0.

**Proof sketch.** Green's identity in the parameters z with two conserved fields. (i) M = A: the Schur split gives A:H(u) = lambda^{-1} det H(B_0+u) + A:Sigma, using A:B = 0 (tangency), and only the wave part reaches the wall. (ii) M = C(H(B_0 + theta u)): it vanishes in the normal direction on the wall (common null line of B and H(r)), and integrating in theta gives lambda^{-1} int det H(B_0+u). The determinant terms cancel, leaving Delta F + int A:Sigma = wall flux. The comparison comes from the rank-one segment expansion and Euler's identity.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §6, Prop. 6.2, Lemma 6.3; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm II; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K2 · **checks:** [verify_quaternion_calculus](reviews/claude/kakeya/scripts/verify_quaternion_calculus.py)

#### <a id="p-kak-paper--kak-e-depth"></a>7. `KAK.E.DEPTH` — Density to depth (rank-r background)

*NEW, PROVED* · prose chain: cites KAK.DIR · uses [`KAK.DIR`](#p-kak-paper--kak-dir)

**Statement.** Background of rank r&lt;n with eigenvalues ~lambda, f=n-r, det(B_fib+H_fib(u))=lambda^r mu, u=0 on the wall: interior -u(x) >= eps rho^2 m^{1/f} if mu>=m on B_rho(x) and m^{1/f}<=lambda; barrier 0<=u/r<=C_2 P if mu<=(c_1 kappa P)^f. The rank-one case (r=1) is the paper's.

**Proof sketch.** Compare u with quadratic and barrier functions through the Monge-Ampere comparison principle. A background of rank r leaves f = n - r free directions, so det(B_fib + kappa I) ~ lambda^r kappa^f, and depth scales like mu^{1/f}.

**Full proof:** [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K3; [kakeya/05](reviews/claude/kakeya/05_DISTILLATION.md) Thm III

#### <a id="p-kak-paper--kak-fib"></a>8. `KAK.FIB` — Energy estimates on a fibre

*PROVED* · prose chain: cites KAK.DIR · uses [`KAK.E.BAL`](#p-kak-paper--kak-e-bal), [`KAK.E.DEPTH`](#p-kak-paper--kak-e-depth)

**Statement.** For cap data at frequency ~1/h and a=Kh: F(0,tau) <= C a^2 ||f||_{p_c}^{p_c} at |tau|=4a; F(T,tau) >= c a^2 ||g_T||_{p_c}^{p_c} at |tau|=2a; Delta_z F <= C s^-2 e^{-cs/h} F for 2a<=s<=3.

**Proof sketch.** Kernel bounds for the Gaussian row multipliers give mu, ||A|| <~ s^{-2d} sigma^m P_\*^d, with an extra e^{-c sigma} on the wall. The upper barrier C sigma^m r P_\* gives F <~ a^2 ||g||^p and a small flux. Inverting the time-row window near the centre gives mu >= m_X, and quadratic comparison on a ball of radius h gives -u >~ h^2 m_X^{1/d}. Integrating gives F >~ a^2 ||g||^p, since (1+1/d)2 = p_c. The aperture gamma is chosen so that decay beats the inversion cost.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §7, Prop. 7.1, Lemma 7.2

#### <a id="p-kak-paper--kak-prop"></a>9. `KAK.PROP` — Two-circle propagation

*PROVED, CHECKED* · prose chain: closed

**Statement.** Delta E <= V(|tau|) E on R x {b<=|tau|<=3}: H(0,r) >= c log(r/b)/log^2(3/b) \* int_1^2 H(T,b) dT; log^{-2} sharp for scalar supersolutions.

**Proof sketch.** Average over rotations. Radial supersolution comparisons (Gronwall in x = log(s/b), using int s log(s/b) V < inf) carry the sin-weighted T-average outward from |tau| = b, losing 1/D. Weak Harnack on a fixed region chains to |tau| = 1. A return passage with the cos(pi T/2) radial solution brings it to |tau| = r, losing log(r/b)/log(1/b).

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §8, Lemma 8.1, Prop. 8.2; [kakeya/07](reviews/claude/kakeya/07_TOOLKIT.md) K5 · **checks:** [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py) · *Note:* Not used by the Lean proof, which takes a cruder route (corrected energy, power 2n+5).

#### <a id="p-kak-paper--kak-t12"></a>10. `KAK.T12` — Critical annular estimate (Theorem 1.2)

*PROVED* · prose chain: cites KAK.DIR · uses [`KAK.FIB`](#p-kak-paper--kak-fib), [`KAK.PROP`](#p-kak-paper--kak-prop), [`KAK.C.WAVE`](#p-kak-paper--kak-c-wave)

**Statement.** ||U(T)P_N f||_{L^{p_c}(R^n x [1,2])} <= C [log(2+N)]^{2/p_c} ||f||_{p_c}.

**Proof sketch.** On one cap on the torus: Proposition 7.1 (KAK.FIB) gives F at |tau| = 4a <~ a^2||f||^p, F at |tau| = 2a >~ a^2||g_T||^p, and Delta F <= V F with V exponentially small. Two-circle propagation (KAK.PROP, b = 2a, r = 4a) gives int_1^2 ||g_T||^p <~ log^2(2/h)||f||^p. Finitely many caps cover the annulus; kernel decay and the cube argument transfer to R^n.

**Alternative proof.** Uses [`KAK.ARROW`](#p-kak-engine--kak-arrow). Second proof: the engine theorem applied to two rotated arrow pencils (spatial passive coordinate), then the transfer of 08 §9. ([kakeya/11](reviews/claude/kakeya/11_FRAMEWORK_POSITIVITY_ENGINE.md) Cor. 11.4)

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 · *Note:* Formalized only with the weaker power 2n+5, which suffices for every consequence. The prose chain cites the Dirichlet theory of KAK.DIR and textbook Riesz-Thorin, Mihlin, Littlewood-Paley.

#### <a id="p-kak-paper--kak-t11"></a>11. `KAK.T11` — Sogge's local smoothing conjecture (Theorem 1.1)

*FORMAL, PROVED* · prose chain: cites KAK.DIR · uses [`KAK.T12`](#p-kak-paper--kak-t12)

**Statement.** ||U(T)f||_{L^p(R^n x [1,2])} <= C ||J^{s_+(p)+eps} f||_p for n>=2, 2<=p&lt;inf, eps>0.

**Proof sketch.** Interpolate Theorem 1.2 at p_c with L^2 (Plancherel) and L^inf (kernel L^1 norm h^{-d/2}) to get h^{-s_+(p)} log^{beta(p)} on each dyadic annulus. The factor 2^{-j eps} sums the annuli; low frequencies by Mihlin.

**Full proof:** [kakeya/08](reviews/claude/kakeya/08_GAPLESS_PROOF.md) §9 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.euclidean_local_smoothing`
<!-- END GENERATED: path-p-kak-paper -->

### Consequences: restriction, Kakeya, Nikodym, Bochner–Riesz, Stein
<a id="p-kak-conseq"></a>

<!-- BEGIN GENERATED: path-p-kak-conseq -->
*Everything downstream of Theorems 1.2 and 1.1.*

Assumes: [`KAK.T12`](#p-kak-engine--kak-t12), [`KAK.T11`](#p-kak-engine--kak-t11).

#### <a id="p-kak-conseq--kak-r-sphere"></a>1. `KAK.R.SPHERE` — Spherical restriction, interior range

*FORMAL, PROVED_CITED* · prose chain: cites KAK.DIR, KAK.R.SPHERE · uses [`KAK.T12`](#p-kak-engine--kak-t12)

**Statement.** E: L^a(S^{n-1}) -> L^b(R^n) for b>p_c, b>(n+1)a'/(n-1).

**Proof sketch.** Average the wave against e^{iNT}psi(T) to obtain a radial multiplier concentrated on |xi| = N. Its kernel is c N^{d/2}|Z|^{-d/2}e^{iN|Z|}psi(|Z|) on 1/2 <= |Z| <= 3, so testing on f_N(r omega) = e^{-iNr}psi(r)g(omega) reproduces Eg on B(0,R) (N = R^{O(1)}). The R^eps local bound becomes global by Tao's eps-removal, and interpolation with L^1 -> L^inf gives the off-diagonal range.

**Full proof:** [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §1-2, Lemma 1.1, Cor. 2.1 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.spherical_restriction` · *Note:* Tao's epsilon-removal cited.

#### <a id="p-kak-conseq--kak-r-cone"></a>2. `KAK.R.CONE` — Truncated cone restriction

*PROVED* · prose chain: cites KAK.DIR, KAK.R.SPHERE · uses [`KAK.R.SPHERE`](#p-kak-conseq--kak-r-sphere)

**Statement.** Same range for the truncated cone, by Nicola's slicing.

**Proof sketch.** Slice the cone by rho: E_Sigma G(x,t) = sqrt2 int e^{it rho} rho^{n-1} E(G_rho)(rho x) d rho. Apply Hausdorff-Young in t, Minkowski, dilation, the sphere bound, and Hoelder in rho.

**Full proof:** [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §2, Cor. 2.1 · **checks:** [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py)

#### <a id="p-kak-conseq--kak-r-parab"></a>3. `KAK.R.PARAB` — Paraboloid restriction and Schrodinger local smoothing

*PROVED_CITED* · prose chain: cites KAK.DIR, KAK.R.PARAB, KAK.R.SPHERE · uses [`KAK.R.SPHERE`](#p-kak-conseq--kak-r-sphere)

**Statement.** Paraboloid extension in the interior range; ||e^{it Delta} f||_{L^{p_c}} <= C log^{2/p_c} ||f||; mixed-norm bounds via Lee-Rogers-Seeger.

**Proof sketch.** Rescale a cap of S^d of size delta parabolically to the paraboloid graph q_delta = |v|^2 + O(delta^2); the powers of delta cancel at p_c. Replace q_delta by |v|^2 at cost O(delta^2 R). For Schroedinger, the Fresnel identity writes e^{it Delta}h as a rescaled paraboloid extension evaluated at (-lambda z/2t, lambda^2/4t).

**Full proof:** [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §4, Cor. 4.1 · **checks:** [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py)

#### <a id="p-kak-conseq--kak-km"></a>4. `KAK.KM` — Kakeya maximal conjecture

*FORMAL, PROVED* · prose chain: cites KAK.DIR, KAK.R.SPHERE · uses [`KAK.R.SPHERE`](#p-kak-conseq--kak-r-sphere)

**Statement.** ||K_delta f||_{L^n(S^{n-1})} <= C log(1/delta)^{2A} ||f||_n.

**Proof sketch.** Each tube is lit by one cap: modulated cap indicators h_j have |Eh_j| >~ delta^d on T_j after rescaling by delta^{-2}. Random signs and Jensen give ||sum a_j 1_{T_j}||_q <~ L^{2A}(delta^d sum a_j^q)^{1/q}, using 2dq = dp. Kakeya maximal follows by discretization and duality.

**Full proof:** [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §3-4, Lemma 3.1, Cor. 4.1 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.kakeya_maximal` · *Note:* The Lean statement is the weaker C_eps delta^-eps form; the paper and 09 prove a power of log(1/delta).

#### <a id="p-kak-conseq--kak-kd"></a>5. `KAK.KD` — Kakeya sets have full Hausdorff dimension

*FORMAL, PROVED* · prose chain: cites KAK.DIR, KAK.R.SPHERE · uses [`KAK.KM`](#p-kak-conseq--kak-km)

**Statement.** Every set containing a unit segment in every direction has dimension n.

**Proof sketch.** Group the covering balls by log(1/r) dyadically (Keich). Each unit segment meets the union of the k-th group in length <~ K_{delta_k}1_{G_k}(omega). Sum, take L^n(S^{n-1}) norms, and apply KAK.KM with the gauge r^n log^{2An} (log log)^gamma. Hoelder forces the gauge measure of the cover to be infinite, so the dimension is n.

**Full proof:** [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §5, Thm 5.1 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.kakeya`

#### <a id="p-kak-conseq--kak-nik"></a>6. `KAK.NIK` — Nikodym maximal conjecture

*PROVED* · prose chain: cites KAK.DIR, KAK.R.SPHERE · uses [`KAK.KM`](#p-kak-conseq--kak-km)

**Statement.** ||N_delta f||_n <= C log(1/delta)^{2A} ||f||_n, via a projective map.

**Proof sketch.** Restrict to half segments at distance ~1 and use the projective map (y',s) -> (y'/(s-t_0), 1/(s-t_0)). It sends all segments through x to parallel segments of direction (x',1), so a Nikodym average at x becomes a Kakeya average of a transported f. Integrate in x', then sum dyadic distances.

**Full proof:** [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §3, Thm 3.1 · **checks:** [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py)

#### <a id="p-kak-conseq--kak-br"></a>7. `KAK.BR` — Bochner-Riesz conjecture (with maximal bounds for p>=2)

*FORMAL, PROVED_CITED* · prose chain: cites KAK.BR, KAK.DIR · uses [`KAK.T11`](#p-kak-engine--kak-t11)

**Statement.** Uniform L^p bounds for alpha>delta(p).

**Proof sketch.** Decompose (1 - |xi|^2/r^2)^alpha into shells of width 2^{-j} around |xi| = r. The shell maximal function is a time average of the wave over |t| <~ 2^j: rescaled local smoothing gives 2^{j sigma}, sigma > s_+(p), summable against 2^{-j alpha}. Littlewood-Paley handles the octaves. Duality gives p < 2; the Bessel kernel gives p = 1, inf.

**Full proof:** [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §6, Thm 6.3 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.bochner_riesz`

#### <a id="p-kak-conseq--kak-sq"></a>8. `KAK.SQ` — Stein's square function conjecture

*FORMAL, PROVED_CITED* · prose chain: cites KAK.DIR, KAK.SQ · uses [`KAK.T11`](#p-kak-engine--kak-t11)

**Statement.** ||G_alpha f||_p <= C||f||_p for alpha>max(1/2, n(1/2-1/p)).

**Proof sketch.** Split the differentiated multiplier into shells delta = 2^{-j} at the boundary. Each shell is a vector-valued multiplier with L^2 bound delta^{1/2} and L^p bound delta^{1/2-sigma} (local smoothing in L^p_x L^2_t). Combine dilates by a sharp-function/Fefferman-Stein argument losing only log^{1/2-1/p}. The sum converges for alpha > 1/2 + s_+(p).

**Full proof:** [kakeya/09](reviews/claude/kakeya/09_GAPLESS_CONSEQUENCES.md) §7, Thm 7.3 · **Lean:** `OAI.Analysis.Fourier.LocalSmoothing.Comparator.stein_square`
<!-- END GENERATED: path-p-kak-conseq -->

### Sharpness
<a id="p-kak-limits"></a>

<!-- BEGIN GENERATED: path-p-kak-limits -->
*Why the losses are what they are: the focusing lower bound log^{1/p} and the weak-norm obstruction. The angular refinement that beats both is on P.KAK.ANGULAR.*

Assumes: [`KAK.T12`](#p-kak-engine--kak-t12).

#### <a id="p-kak-limits--kak-ob-focus"></a>1. `KAK.OB.FOCUS` — Radial focusing (sharpness of s_+(p) and of the power 1/p_c)

*PROVED, CHECKED* · prose chain: closed

**Statement.** Lower bounds (8.1).

**Proof sketch.** Radially focusing data (phase e^{(3/2)i|xi|} on the annulus |xi| ~ N) concentrate the wave near the origin at T ~ 3/2. Stationary phase gives output modulus ~ N^{n-d/2} delta^{-d/2} on a set of volume ~ delta^d/N at time 3/2 + delta. Integrating gives int delta^{d - dp/2} d delta, which is ~ log N exactly at p = p_c (d - dp/2 = -1): a log^{1/p_c} loss is necessary, and s_+(p) is sharp.

**Full proof:** [kakeya/02](reviews/claude/kakeya/02_VERIFICATION_LEDGER.md) row for paper Prop. 8.1 · **checks:** [verify_propagation_and_exponents](reviews/claude/kakeya/scripts/verify_propagation_and_exponents.py)

#### <a id="p-kak-limits--kak-ob-weak"></a>2. `KAK.OB.WEAK` — Weak-norm obstruction

*PROVED, CHECKED* · prose chain: closed

**Statement.** ||U(t)g_N||_{L^{p,inf}} >~ (log N)^{d(1/2-1/p)} ||g_N||_p for 2&lt;p<=p_c.

**Proof sketch.** Binary-digit families of 2^{kd} lines (k ~ (1/4)log N) are separated at the input side but compressed to measure ~ k^{-d} at the output side. Wave packets along them, with random signs, put mass ~1 on a set of measure k^{-d}, so the weak L^p norm grows like k^{d(1/2-1/p)}.

**Full proof:** [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §1, Prop. 1.1 · **checks:** [verify_remaining_claims](reviews/claude/kakeya/scripts/verify_remaining_claims.py)
<!-- END GENERATED: path-p-kak-limits -->

### The angular endpoint (App. B, made gapless in 14)
<a id="p-kak-angular"></a>

<!-- BEGIN GENERATED: path-p-kak-angular -->
*How a little angular regularity removes one logarithm (strong norm) or all of them (weak norm): uniform Hankel -> shell transplant -> four-kernel envelope -> one annulus -> summation.*

Assumes: [`KAK.T12`](#p-kak-engine--kak-t12).

#### <a id="p-kak-angular--kak-ang-hankel"></a>1. `KAK.ANG.HANKEL` — Uniform Hankel expansion on an angular band

*NEW, PROVED, CHECKED* · prose chain: closed

**Statement.** For degree <= M and z >= (2+M)^K: S_z F = (2pi)^alpha z^-alpha sum_{+-} e^{+-i(z - pi alpha/2)}(P_+- + E_{+-,z})F with P_+ = Id, P_- the antipodal map, and E (with log-derivatives) bounded L^1 -> L^inf uniformly in M. Explicitly |(x d_x)^k eps_mu(x)| <= e^{3/2}(2mu+2k+1)^{max(k-1,0)}(mu+1)^2/(2x) for x >= (mu+1)^2.

**Proof sketch.** Funk-Hecke makes S_z a multiplier by (2pi)^{n/2} i^nu z^{1-n/2} J_{nu+n/2-1}(z) on degree-nu harmonics. Write J = Re H^(1) and use Hankel's integral: eps_mu(x) is an average of (1 + iu/2x)^{mu-1/2} - 1 against e^{-u}u^{mu-1/2}. With w = (1+iy)^-1 every log-derivative carries a factor (1-w) = O(y), which gives the explicit bound for x >= (mu+1)^2. Summing over a band uses |Z_nu| <= Z_nu(1) = dim H_nu/|S^d| and sum dim H_nu <= 2(M+1)^d. The threshold (mu+1)^2 is necessary: x|eps| ~ mu^2/2 there.

**Full proof:** [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 1.2, Proposition 1.3, Cor. 1.4 · **checks:** [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) · *Note:* Replaces the paper's pole-coordinate stationary phase (whose key identity lambda d_lambda I = -1/2 I(w.grad a) is also checked) by an explicit Bessel bound.

#### <a id="p-kak-angular--kak-ang-ball"></a>2. `KAK.ANG.BALL` — Observation on a ball (shell transplant)

*PROVED, CHECKED* · prose chain: cites KAK.DIR · uses [`KAK.T12`](#p-kak-engine--kak-t12), [`KAK.ANG.HANKEL`](#p-kak-angular--kak-ang-hankel)

**Statement.** For N >= 2, 0 < r <= 1: ||U(t)Pi_N f||_{L^p(I x B(x_0,r))} <~ min{1,(Nr)^{n/p}} L(Nr)^A ||f||_p, given Theorem 1.2 with power A.

**Proof sketch.** Rescale to frequency 1, T = N, R = Nr. Small R: Young in time on the unit kernel. R <= T <= 10S with S = 100R^{3+n/p}: the global estimate, since L(T) <~ L(R). T > 10S: an input shell at radius Q produces on the ball a wave T h(b,x) = int k(b-u+x.omega)h independent of Q up to R^{2+n/p}/Q, so move the shell to radius S where the global estimate costs L(S)^A ~ L(R)^A; the factor S^{1/p} cancels against ||f_S||_p ~ S^{-1/p}||h||_p because d - alpha p = -1. Cut [T,2T] into intervals of length ~R; the remaining input is off the light cone.

**Full proof:** [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 3.2, Lemma 3.1 · **checks:** [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py)

#### <a id="p-kak-angular--kak-ang-env"></a>3. `KAK.ANG.ENV` — Four-kernel envelope and small source

*PROVED, CHECKED* · prose chain: closed · uses [`KAK.ANG.HANKEL`](#p-kak-angular--kak-ang-hankel)

**Statement.** If g has angular degree <= M, support in |y| >= s_0, and N min(r,s_0) >= (2+M)^K, then |U(t)Pi_N g(r theta)| <~ sum_+- G_g(t +- r, theta) r^{-alpha} (r <= 4; r^{-d/p} for r > 4), with ||G_g||_{L^p(dv dtheta)} <~ ||g||_p, uniformly in M, N. The input from |y| < N^{eta-1} costs N^{-delta} in strong L^p.

**Proof sketch.** Polar coordinates give s^d int e^{-it rho} chi(rho/N) rho^d S_{r rho} S_{-s rho} g_s d rho. Expand both spherical operators by KAK.ANG.HANKEL into four phases t +- r +- s with amplitudes bounded at theta by |g_s(+-theta)| + ||g_s||_1. Integrate by parts in rho to get N&lt;N(t +- r +- s)>^{-D}; the weight s^alpha = s^{1/p} s^{d/p} is absorbed, and evenness in s turns the two s-phases into one convolution, bounded by the maximal function. The small source: interpolate L^2 and the N^{alpha+1} dispersive bound, then Holder on the ball.

**Full proof:** [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Lemma 4.2, Lemma 4.3, Definition 4.1 · **checks:** [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py)

#### <a id="p-kak-angular--kak-ang-band"></a>4. `KAK.ANG.BAND` — One annulus, bounded angular degree

*PROVED, CHECKED* · prose chain: cites KAK.DIR · uses [`KAK.T12`](#p-kak-engine--kak-t12), [`KAK.ANG.BALL`](#p-kak-angular--kak-ang-ball), [`KAK.ANG.ENV`](#p-kak-angular--kak-ang-env)

**Statement.** For angular degree <= M: ||U(t)Pi_N f||_p <~ [L(N)^{1/p} + L(M)^A]||f||_p and ||U(t)Pi_N f||_{p,inf} <~ L(M)^A ||f||_p, uniformly in M, N.

**Proof sketch.** H = (2+M)^K. If H > N^eta, the global estimate. Otherwise: the inner ball r <= H/N costs L(H)^A by KAK.ANG.BALL; outside it, drop the small source and apply the envelope. Change variables to v = t +- r: int_{1/N}^4 (r^{-alpha}G)^p r^d dr = G^p log 4N (strong: L(N)^{1/p}) while |{r^{-alpha}G > lambda}|_{r^d dr} = G^p/(n lambda^p) (weak: no logarithm), both because d - alpha p = -1.

**Full proof:** [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Proposition 5.1 · **checks:** [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py)

#### <a id="p-kak-angular--kak-app-b"></a>5. `KAK.APP.B` — Angular endpoint

*PROVED_CITED, CHECKED* · prose chain: cites KAK.APP.B, KAK.DIR · uses [`KAK.T12`](#p-kak-engine--kak-t12), [`KAK.ANG.BALL`](#p-kak-angular--kak-ang-ball), [`KAK.ANG.ENV`](#p-kak-angular--kak-ang-env)

**Statement.** At p_c, for B > 1/2 + 1/p: ||U(t)f||_{L^p(R^n x I)} <~ ||L(J)^{1/p} f||_p + ||L(Lambda_omega)^B f||_p, and ||U(t)f||_{L^{p,inf}} <~ ||L(Lambda_omega)^B f||_p; Lambda_omega^sigma (sigma > 0) may replace L(Lambda_omega)^B.

**Proof sketch.** Littlewood-Paley in space; group angular degrees in double-exponential bands M_m = 2^{2^m} so that L(degree) ~ 2^m. At output radius r and frequency 2^j, select the largest band m(j,r) whose threshold H_m = (2+2M_m)^K is <= min(N_j^eta, N_j r): the selected low-pass wave obeys the four-kernel envelope with a majorant (via the spherical maximal function) independent of r, so the radial integral pays L(N_j)^{1/p} in strong norm and nothing in weak norm (vector-valued maximal inequalities). Unselected band m lives on balls of H_m wavelengths and costs L(H_m)^{A + 1/2 - 1/p} by the shell transplant and l^2-vs-l^p over ~log H_m frequencies; summing against ||Q_m f||_p <~ 2^{-mB}||L(Lambda_omega)^B f||_p needs B > A + 1/2 - 1/p = 1/2 + 1/p. Reassemble with the synthesis operator (real interpolation for the weak norm).

**Full proof:** [kakeya/14](reviews/claude/kakeya/14_GAPLESS_ANGULAR_ENDPOINT.md) Theorem 14.1, §6; [kakeya/10](reviews/claude/kakeya/10_GAPLESS_REMAINING.md) §5 · **checks:** [verify_angular_endpoint](reviews/claude/kakeya/scripts/verify_angular_endpoint.py) · *Note:* Cited: Littlewood-Paley and synthesis (Grafakos), Fefferman-Stein, vector-valued maximal on S^d (Grafakos-Liu-Yang 2009), localized spherical cut-off kernels (Narcowich-Petrushev-Ward 2006), spherical multipliers (Bonami-Clerc 1973), real interpolation; hypotheses checked in 14 §0.3. Whether B > 1/2 + 1/p is sharp is open (T.KAK.13).
<!-- END GENERATED: path-p-kak-angular -->

## ETH

Notation used throughout the ETH paths:

<!-- BEGIN GENERATED: glossary-eth -->
| Term | Meaning | Defined in |
|---|---|---|
| rho(n) | n (log log n)^2/sqrt(log n): the formally verified exponent (Prop. 6.4). | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §1 |
| deformation, nodes | Rows P_{f,v}(T) with P(0) = identity, of degree &lt;= R'. At each node t in [Q], the rows are local for a measured partition (groups S_g, measurement matrices C_g). | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4 |
| measured partition, local row | A partition [n] = S_1 u ... u S_r with integer matrices C_g. A row is g-local if its restriction to every other group lies in the row space of C_h, so the sum factorizes once the measurements y = C_h x_h are fixed. | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4 |
| swap design Pi_S, M_0, Krawtchouk | Commuting involutions on labels ([r]x[r])^Q. M_0 = sum gamma_k sum_{\|S\|=k} Pi_{S^c} acts by f(j) on level j; P(T) puts prod_{s in S}(1 - T/s) in front of Pi_{S^c}. | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §5 |
| Taylor depth q, jets | q = floor(dR'/Q); the jet of H at a node mod z^{q+1} costs C_q(n) Boolean substitutions. Its log, x log(e^2/x), is the dominant cost term. | [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4 |
| parity lift p_ell, h, digits e_j | h = Y^2 - Y; p_ell = sum_{j&lt;ell} e_j h^j is congruent to [2\|y] mod 2^ell; e_j = (-1)^{j+1}C(2j-1,j-1)(2Y-1). | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §2 |
| V-truncation | Multiply the per-equation series B_i(V) and keep V^{&lt;ell}: the degree stays d_\*(m + 2(ell-1)) instead of m \* ell. | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §2 |
| Tanner tests, greedy rows | A^(0): deterministic parity checks catching every nonzero violation vector of weight &lt; m/P. Greedy rows: append exact-count-minimizing parity rows until N(A) = #F. | [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §4 |
| Cauchy kernel W_d, nodes S | W_d(T) = prod_{q in S \ {d}} (T-q)/C'(0), S = {-L..L-1}, L = 2^a. W(t) is the shift by t; it is 2-adically small off the diagonal at t + 2^{a+1}z. | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §3 |
| packets, tiles, modes, children | Edges are packed into slots of L; crossing tiles expand into modes L_Z = sum j^b X_y of 2-adic weight b+1. Fixing a mode tuple splits the network into children of &lt;= 6M/r edges. | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §6 |
| ct(N) | The contraction (sum over Boolean edge assignments of the product of vertex tensors) of a tensor network N. | [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §1 |
| readout system | Boolean variables, factors Phi_f with ports, and a readout rho; Z(rho) = sum_x prod_f Phi_f(rho(x)). CNFs, weighted factor graphs and tensor networks are readout systems; the deformation acts on the readout. | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) §0 |
| local node, measured locality | A parameter value where every port of a group-g factor depends only on x_g and the measurements y = (C_h x_h)_h; the sum then factorizes over groups (separation). | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Definition 1.1 |
| mode, weight, source, target | Near a node the readout is local part + sum kappa_Z L_Z with kappa_Z of filtration weight w_Z; L_Z reads its source group; targeted modes are read only by one group. | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Definition 2.2 |
| Hermite loss l_t | -v(E_t(0)) = sum_q m_q (v(t-q) - v(q)): how much closer the other nodes are to t than to 0. Zero for unit-distance nodes; D_i Delta(i) on residue systems. | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Lemma 3.1 |
| depth profile, slope Delta(i) | Multiplicities m(v_2(t)) on the residue system mod 2^{a+1}; Delta(i) = sum_j (j/2^j)(m(i) - m(i+j)) sets both the loss and the node cost. | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 3.3 |
| square-root law | Both architectures pay M/sqrt(number of node coordinates) per call: Part I via the central limit window, Part III via Hardy-Ramanujan; open in general (T.ETH.5). | [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) §7 |
<!-- END GENERATED: glossary-eth -->

### ETH is false (Part I)
<a id="p-eth-i"></a>

<!-- BEGIN GENERATED: path-p-eth-i -->
*The shortest complete refutation of ETH: deformation counting with swap interpolation.*

#### <a id="p-eth-i--eth-i-spars"></a>1. `ETH.I.SPARS` — Disjoint sparsification

*PROVED, CHECKED* · prose chain: closed

**Statement.** For 0&lt;sigma<=1 a k-CNF on n variables is split by a branching tree into at most (C_M n+1)^2 2^{sigma n} normalized CNFs of width <=k, whose satisfying sets partition Sat(F), each variable occurring in at most B=O_k((sigma^{-1}log(2/sigma))^{k-1}) clauses.

**Proof sketch.** Branch on an eligible pair (c,H): H is a proper subset of at least M^{c-|H|} width-c clauses H u P_a. Child A adds the clause H (some literal of H is true); child B adds every P_a and the units not-l for l in H (all of H false). Satisfying sets split disjointly. Insertion accounting (the least-width, largest-|H| priority keeps fewer than 2M^{j-i} marked width-j clauses through any i-set) bounds width-j insertions on a path by 4(3M)^{j-1}n, which bounds the depth and the number of B-steps. Leaves have no eligible pair, hence occurrence <= B. Taking M ~ sigma^{-1}log(2/sigma) gives (C_M n+1)^2 2^{sigma n} nodes.

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §2, Prop. 2.4, Lemma 2.5 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py)

#### <a id="p-eth-i--eth-i-tovey"></a>2. `ETH.I.TOVEY` — OR gates and Tovey cycles

*PROVED, CHECKED* · prose chain: closed

**Statement.** A CNF of width <=k with occurrence bound B gives, in polynomial time, a CNF L' of width <=3 with every variable in 1..3 clauses and N<=6kBn variables, and z with #L=2^z #L'.

**Proof sketch.** Drop free variables (factor 2^z). Shorten clauses with OR gates z = l1 v l2 (four clauses that determine z, so counts are preserved), each adding 6 occurrences. Then replace every variable occurring p>=2 times by a cycle of copies x_1 -> x_2 -> ... -> x_p -> x_1, which forces equality. Every copy then occurs exactly 3 times, and N is the total number of occurrences, at most 6kBn.

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §3, Prop. 3.3 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py)

#### <a id="p-eth-i--eth-i-defcount"></a>3. `ETH.I.DEFCOUNT` — Deformation counting

*PROVED, CHECKED* · prose chain: closed

**Statement.** Given a deformation of degree <=R' with Q&lt;p nodes whose rows are local for measured partitions, the sum over {0,1}^n of a product of factors (clauses, or any factors with values in a free Z/p^mu-algebra, degree measured in the algebra) is computed with Q\*C_q(n)\*(2beta+1)^{c_\*}\*2r\*2^s\*poly operations, q=floor(dR'/Q).

**Proof sketch.** H(T) = sum_x prod_f Phi_{C_f}(&lt;P_{f,v}(T),x>) has degree <= dR' and H(0) = #F. At each node t, the jet of H(t+z) mod z^J is a sum over x of polynomials of degree <= q in x. Moebius inversion on the Boolean lattice reduces it to Boolean substitutions over sets E with |E| <= q (C_q(n) of them). Locality of the rows makes each sum factor over the groups once the measurement values y in [-beta,beta]^{c_\*} are fixed (separation by measurements). The nodes 1..Q have unit differences mod p > Q, so the Chinese remainder theorem recovers H(0) mod p^mu.

**Alternative proof.** Uses [`ETH.U.MAIN`](#p-eth-unified--eth-u-main). As the formal-neighbourhood instance of the unified deformation theorem: unit nodes, Boolean coordinate modes, measurements. ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 4.3)

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4, Thm 4.6; [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §1, Thm 1.1 (general factors) · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py)

#### <a id="p-eth-i--eth-i-swap"></a>4. `ETH.I.SWAP` — Swap interpolation in closed form

*PROVED, CHECKED* · prose chain: closed

**Statement.** P(T)=I-M_0+sum_k gamma_k sum_{|S|=k} prod_{s in S}(1-T/s) Pi_{S^c} has P(0)=I, degree <=R-1, coefficients in Z[1/(2Q!)], and P(t) is local on ker C at every node t; c/D <= 2exp(-R^2/(8Q)).

**Proof sketch.** The swaps Pi_s are commuting involutions on the labels ([r]x[r])^Q. Their joint eigenbasis is the tensor products of symmetric and antisymmetric local vectors, and sum_{|S|=k} Pi_{S^c} acts by (-1)^j K_k(j) (Krawtchouk) on vectors with j antisymmetric factors. Take f of degree < R with f(j) = (-1)^j on R consecutive levels around Q p_r, written as f = sum gamma_k K_k with gamma_k in Z[1/2]. Then M_0 = sum gamma_k sum Pi_{S^c} is the identity on ker C. P(T) replaces Pi_{S^c} by prod_{s in S}(1-T/s) Pi_{S^c}, which vanishes at the nodes t in S, so P(t) is local. The defect is the binomial tail Pr[W not in window] <= 2e^{-R^2/(8Q)} (Hoeffding).

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §5, Thm 5.4, Cor. 5.6; [eth/07](reviews/claude/eth/07_TOOLKIT.md) E1 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_distilled_swap_formula](reviews/claude/eth/scripts/verify_distilled_swap_formula.py), [verify_replacement_forms](reviews/claude/eth/scripts/verify_replacement_forms.py)

#### <a id="p-eth-i--eth-i-transfer"></a>5. `ETH.I.TRANSFER` — Transfer to a formula; balanced schedules

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.SWAP`](#p-eth-i--eth-i-swap)

**Statement.** Labels from derandomized schedules turn the swap rows into a deformation with measurement rows of norm <=2^Q, at most 2hc rows per node, h<=2M/D, groups <=2n/r.

**Proof sketch.** Give every (clause, variable) incidence a label u in B through schedules a_v(t), d_f(t) in [r]. The rows P_{f,v}(T) = P(T)_u pi iota are then local at each node, and the measurement rows are the rows of C restricted to the groups S_{t,g}. Random schedules have small variance in loads and group sizes; derandomize by conditional expectations of a quadratic potential.

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §6, Lemmas 6.1-6.2 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py)

#### <a id="p-eth-i--eth-i-param"></a>6. `ETH.I.PARAM` — Counting at <=3 occurrences, explicit parameters

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount), [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer)

**Statement.** A CNF of width <=3 with occurrences in 1..3 is counted in poly(n)\*2^{13 n nu^2/sqrt(lambda)} for n>2^{2^17}; the exact bound is nontrivial for lambda>9,297,119.

**Proof sketch.** Set lambda = ceil(log n), nu = ceil(log lambda), r = floor(sqrt lambda), Q = floor(lambda/(4 nu)), R = 2 ceil(sqrt lambda), p the least prime > Q and mu = n+1. The time exponent per variable is x log(e^2/x) (Taylor jets, x = 3(R-1)/Q) + 24e^{-2nu}(Q+2) (measurements) + 2/r (groups). Under the inequalities (P0)-(P8), valid for lambda > 2^17, this is at most 13 nu^2/sqrt(lambda).

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §7, Lemma 7.2, Thm 7.1 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py)

#### <a id="p-eth-i--eth-i-main"></a>7. `ETH.I.MAIN` — ETH is false (Part I)

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.SPARS`](#p-eth-i--eth-i-spars), [`ETH.I.TOVEY`](#p-eth-i--eth-i-tovey), [`ETH.I.PARAM`](#p-eth-i--eth-i-param)

**Statement.** #k-SAT is computable in 2^{O_k(n(log log n)^{k+1}/(log n)^{1/(2k)})} poly(n+m) time and polynomial space; hence ETH, SETH, #ETH fail.

**Proof sketch.** Normalize; sparsify with sigma = 1/ceil((log n)^{1/(2k)}); reduce each leaf to width 3 and occurrence 3 on N <= 6kBn variables; count each leaf by ETH.I.PARAM; add the counts. B = O_k((log n)^{(k-1)/(2k)}(log log n)^{k-1}) inflates N, giving the exponent n(log log n)^{k+1}/(log n)^{1/(2k)}. Polynomial space comes from depth-first traversal.

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §8 · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_counting_pipeline](reviews/claude/eth/scripts/verify_counting_pipeline.py) · *Note:* Also implied by the formally verified ETH.II.MAIN.
<!-- END GENERATED: path-p-eth-i -->

### The formally verified exponent rho(n) (Part II)
<a id="p-eth-ii"></a>

<!-- BEGIN GENERATED: path-p-eth-ii -->
*Reduction to sparse F_2 systems; parity lift; the bound OpenAI proved in Lean.*

Assumes: [`ETH.I.SPARS`](#p-eth-i--eth-i-spars), [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount), [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer).

#### <a id="p-eth-ii--eth-ii-lift"></a>1. `ETH.II.LIFT` — Parity lift in closed form; cross-equation truncation

*PROVED, CHECKED* · prose chain: closed

**Statement.** e_0=1-Y, e_j=(-1)^{j+1}C(2j-1,j-1)(2Y-1); p_ell=sum_{j&lt;ell} e_j h^j (h=Y^2-Y) satisfies p_ell(y)=[2|y] mod 2^ell; S_ell(prod_i B_i mod V^ell)=prod_i [2|y_i] mod 2^ell with retained degree d_\*(m+2(ell-1)); ||e_j||_1<=2\*4^j.

**Proof sketch.** In Q[Y]/(h^ell) with h = Y^2-Y, (2Y-1)^2 = 1+4h, so eps = (2Y-1)(1+4h)^{-1/2} (binomial series) satisfies eps^2 = 1, and (1-eps)/2 is an idempotent. Expanding it gives the digits e_j. It equals 1 mod Y^ell and 0 mod (Y-1)^ell, hence agrees with the Beigel-Tarui polynomial and is congruent to [2|y] mod 2^ell. For the truncation, multiply the series B_i(V) = sum e_j(y_i)h(y_i)^j V^j and keep V^{&lt;ell}. The coefficient of V^s has degree <= d_\*(m+2s), because digit j has degree 1 and h^j has degree 2j.

**Full proof:** [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §2, Lemmas 2.1-2.2 · **checks:** [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py), [verify_parity_lift](reviews/claude/eth/scripts/verify_parity_lift.py) · *Note:* The closed form is new; it replaces App. B.1's 2^{C ell^2}.

#### <a id="p-eth-ii--eth-ii-zeros"></a>2. `ETH.II.ZEROS` — Counting common zeros of sparse F_2 systems

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount), [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer), [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift)

**Statement.** m<=cn equations over F_2, each with at most (log n)^c monomials of degree <=k: common zeros counted in 2^{O_{k,c}(rho(n))} time, poly space.

**Proof sketch.** P_i(x) = 0 iff the integer count L_i(x) of satisfied monomials is even. Branch on the < n/sqrt(lambda) variables of high degree. On each branch, take the factors B_i(L_i(X),V) in the algebra (Z/p^mu)[V]/(V^ell) with ell = n+1. The V-truncation keeps the degree k(m+2n), linear in n. Count with general-factor deformation counting (Part I design, R = 2(a+2)ceil(sqrt lambda)), recover the integer coefficients b_s with p^mu > 2Z, and read N_branch = sum_s b_s mod 2^ell.

**Full proof:** [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §3, Thm 3.1 · **checks:** [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py)

#### <a id="p-eth-ii--eth-ii-tanner"></a>3. `ETH.II.TANNER` — Initial tests from a Tanner code

*PROVED, CHECKED* · prose chain: closed

**Statement.** For m<=Bn there is a deterministic A^(0) with <=3n rows of <=D_0 ones such that every nonzero kernel vector has weight >=m/P, P=20000 B^2 mu^2.

**Proof sketch.** Tanner code on a Cayley-type bipartite graph on {0,1}^{l_G} with eps-biased shifts (found by exhaustive search; they exist by Hoeffding). The local code at each vertex is a random parity-check matrix with kernel of minimum weight >= delta d_G. Expander mixing forces every nonzero accepted word to have weight >= delta(delta-eps)q_G d_G >= m/P. Shortening keeps <= 3n rows with <= D_0 ones each.

**Full proof:** [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §4, Lemma 4.1 · **checks:** [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py)

#### <a id="p-eth-ii--eth-ii-greedy"></a>4. `ETH.II.GREEDY` — Greedy rows by exact counts

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.II.TANNER`](#p-eth-ii--eth-ii-tanner)

**Statement.** Appending 3n rows, each chosen among (2m)^P candidates to minimize the exact common-zero count, makes N(A)=#F.

**Proof sketch.** N(A) = #F + D(A), where D counts unsatisfying x with A v(x) = 0. For such x, v(x) has weight >= m/P by ETH.II.TANNER, so a random P-term row rejects x with probability >= 1/4. Choosing, among all (2m)^P candidate rows, the one minimizing the exactly computed N gives D <= (3/4)D. After 3n rounds D < 1, so D = 0, without ever knowing #F.

**Full proof:** [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §4, Lemma 4.2 · **checks:** [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py)

#### <a id="p-eth-ii--eth-ii-main"></a>5. `ETH.II.MAIN` — #k-SAT in 2^{O(rho(n))} (formally verified bound)

*FORMAL, PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.SPARS`](#p-eth-i--eth-i-spars), [`ETH.II.ZEROS`](#p-eth-ii--eth-ii-zeros), [`ETH.II.TANNER`](#p-eth-ii--eth-ii-tanner), [`ETH.II.GREEDY`](#p-eth-ii--eth-ii-greedy)

**Statement.** #k-SAT in (n+m)^{O_k(1)} 2^{O_k(rho(n))} time and (n+m)^{O_k(1)} space, rho(n)=n(log log n)^2/sqrt(log n).

**Proof sketch.** Sparsify at sigma = nu^2/sqrt(lambda) (poly(n)2^{O(rho)} leaves with polylog occurrences). For each leaf, build the Tanner tests and run the greedy procedure. Each of its 2^{polylog} calls counts common zeros of <= 6n sparse F_2 equations by ETH.II.ZEROS, at cost 2^{O(rho(n))}. Sum the leaf counts.

**Full proof:** [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) Thm 9.1; [eth/06](reviews/claude/eth/06_LEAN_STATEMENT_AUDIT.md) §6.A · **checks:** [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py) · **Lean:** `OAI.Algorithms.SubexponentialSAT.SqrtLogSpace.exact_count_and_sat`
<!-- END GENERATED: path-p-eth-ii -->

### The headline exponent (Part III, App. D)
<a id="p-eth-iii"></a>

<!-- BEGIN GENERATED: path-p-eth-iii -->
*Boolean tensor contraction by a 2-adic Cauchy-kernel deformation.*

Assumes: [`ETH.I.SPARS`](#p-eth-i--eth-i-spars), [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift), [`ETH.II.TANNER`](#p-eth-ii--eth-ii-tanner), [`ETH.II.GREEDY`](#p-eth-ii--eth-ii-greedy).

#### <a id="p-eth-iii--eth-iii-ring"></a>1. `ETH.III.RING` — Coefficient ring and interpolation points

*PROVED* · prose chain: closed

**Statement.** O=Z_(2)[U]/(f), f irreducible mod 2 of degree e: v_2 is multiplicative, units are the elements nonzero mod 2, and the 2^e binary-coordinate points have unit differences.

**Proof sketch.** O/2O = F_2[U]/(f) = F_{2^e} is a field. So an element is a unit iff it is nonzero mod 2, v_2 is multiplicative, and distinct binary-coordinate points differ by units. Vandermonde determinants are therefore units, and Lagrange interpolation works over A_p.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §2, Lemma 2.1

#### <a id="p-eth-iii--eth-iii-cauchy"></a>2. `ETH.III.CAUCHY` — Cauchy kernel

*PROVED, CHECKED* · prose chain: closed

**Statement.** W_d(T)=C(T)/(C'(0)(T-d)), S={-L..L-1}: W(0)=I, W(t) is the shift by t, W_d(t+alpha z) is integral and divisible by 2^{a+1-v_2(t-d)} for d != t.

**Proof sketch.** S = {-L..L-1} is a complete residue system mod 2^{a+1}. W_d vanishes on S \ {d} and equals [d=0] at T=0, so W(t) is the shift by t on slots. At T = t + alpha z each factor (t-q)+alpha z is 2^{v_2(t-q)} times a 2-adic unit, so W_d(t+alpha z) is integral and divisible by 2^{a+1-v_2(t-d)}.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §3, Lemma 3.1, Cor. 3.2 · **checks:** [verify_appendix_d](reviews/claude/eth/scripts/verify_appendix_d.py)

#### <a id="p-eth-iii--eth-iii-hermite"></a>3. `ETH.III.HERMITE` — Tapered 2-adic Hermite recovery

*PROVED, CHECKED* · prose chain: closed

**Statement.** With multiplicities q_t=h max(0,u-v_2(t)), F(0) mod 2^p is recovered from F(t+alpha z) mod 2^{p_t}, p_t=p+8hD_t+8M(a+1); loss v_2(E_t(0))>=-6hD_t.

**Proof sketch.** Hermite interpolation F(0) = sum_t E_t(0) P_t(-t/alpha), with node multiplicities q_t = h(u - v_2(t))^+ tapered by 2-adic depth (sum q_t > deg F). The 2-adic valuation of E_t(0) is >= -6hD_t (a geometric sum over depth levels), so knowing each jet mod 2^{p_t} with p_t = p + 8hD_t + 8M(a+1) determines F(0) mod 2^p.

**Alternative proof.** Uses [`ETH.U.HERMITE`](#p-eth-unified--eth-u-hermite). As Hermite recovery with valuation loss on the residue system mod 2^{a+1}: loss D_i Delta(i) <= 6 h D_i for the taper. ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 3.3)

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §4, Lemma 4.1; [eth/07](reviews/claude/eth/07_TOOLKIT.md) E5 · **checks:** [verify_appendix_d](reviews/claude/eth/scripts/verify_appendix_d.py), [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-iii--eth-iii-pack"></a>4. `ETH.III.PACK` — Packing edges into packets

*PROVED* · prose chain: closed

**Statement.** Derandomized labels in [r]^a and packets of L=2^a slots with orientation bits; <=3M/L packets; <=3M/r ports per symbol per side.

**Proof sketch.** Uniform random vertex labels in [r]^a balance the counts of edges per (row label, column label) pair; derandomize. Edges sharing a label pair fill packets of L = 2^a slots, with bit c of the slot recording which symbol of the coordinate-c pair is on the row side. Each symbol then meets <= 3M/r edges per side.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §5, Lemma 5.1

#### <a id="p-eth-iii--eth-iii-child"></a>5. `ETH.III.CHILD` — Modes and the child identity

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack)

**Statement.** At T=t+alpha zeta, crossing entries expand in tile moments sum_j j^b X_y with 2-adic weight b+1, and F(T)=sum_{beta: sum w beta&lt;P} prod_g ct(N_{g,beta}) mod 2^P, with children of <=6M/r edges and exp(2 pi sqrt(MP/D_t)) mode tuples.

**Proof sketch.** At a node t of depth i, an entry W_{xy} for a crossing tile y_0 is a geometric series in j (y = y_0 + 2^{i+1}j) with ratio of valuation >= 1. So rows become lambda_x(X) + sum_Z kappa_{x,Z} L_Z(X), with modes L_Z = sum_j j^b X_y of weight b+1. Expanding in formal modes, the coefficient of Z^beta has valuation >= sum w beta, so only tuples with sum w beta < P survive. For each tuple the summand factors over the groups r into children with <= 6M/r edges. Rankin's bound counts the tuples by exp(2 pi sqrt(MP/D_t)).

**Alternative proof.** Uses [`ETH.U.NODE`](#p-eth-unified--eth-u-node), [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack). As the targeted node lemma with tile power moments as modes (weights b+1, source = tile group, target = crossing rows). ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 2.7)

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §6, Lemmas 6.1-6.4 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-iii--eth-iii-tensor"></a>6. `ETH.III.TENSOR` — Boolean tensor contraction

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.RING`](#p-eth-iii--eth-iii-ring), [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.HERMITE`](#p-eth-iii--eth-iii-hermite), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack), [`ETH.III.CHILD`](#p-eth-iii--eth-iii-child)

**Statement.** A bipartite Boolean tensor network with <=Cs edges, degree <=s^{9/10}, polynomially many vertices and suitable entry oracles is contracted in (s+L)^{O(1)} 2^{O(s/sqrt(log s))} time, poly space.

**Proof sketch.** Recursive algorithm CONTRACT(N,p,d): if N is small, sum directly. Otherwise pack, evaluate F at each node t + alpha z by the child identity (the children are contracted recursively at precision p_t), and recover F(0) by tapered Hermite. Explicit parameters make the recursion depth O(log log) and the exponent O(s/sqrt(log s)), via the absorption inequalities (A1)-(A2).

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §7, Lemmas 7.1-7.3, Thm 10.2 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) · *Note:* Paper's floor u>=lambda^{4/5} holds only from log2(lambda)~500; absorption inequalities (A1)-(A2) suffice from ~26.

#### <a id="p-eth-iii--eth-iii-fixed"></a>7. `ETH.III.FIXED` — Fixed clause density

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.TENSOR`](#p-eth-iii--eth-iii-tensor)

**Statement.** #k-SAT with m<=Dn in (n+m)^{O(1)} 2^{O_{k,D}(n/sqrt(log n))} via the clause network.

**Proof sketch.** The clause network has variables and clauses as vertices, literal tables as entries, and O_k(n+m) edges, so its contraction is #F. Apply ETH.III.TENSOR with s = n+m <= (1+D)n.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §8, Prop. 8.1 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-iii--eth-iii-sparsez"></a>8. `ETH.III.SPARSEZ` — Sparse systems by contraction

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.TENSOR`](#p-eth-iii--eth-iii-tensor), [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift)

**Statement.** Sparse F_2 systems are counted in n^{O(1)} 2^{O(n log log n/sqrt(log n))} via incidence-subset networks with O(n log log n) wires.

**Proof sketch.** Branch on high-degree variables (at most 2^{n/sqrt(lambda)} branches). Give each incidence its own indeterminate, so the parity-lifted product F = prod_i B_i (ETH.II.LIFT) is multilinear of degree K = O(n). Group tests and variables into n^{1/4} groups; for each array of per-block monomial degrees (2^{O(sqrt n log n)} arrays), a network whose wire bits index subsets of incidences (O(n nu) wires) contracts to the corresponding part of sum_x F. Contract each by ETH.III.TENSOR with s = n nu.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §9, Lemma 9.1 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-iii--eth-iii-main"></a>9. `ETH.III.MAIN` — Headline Theorem 1.1

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.FIXED`](#p-eth-iii--eth-iii-fixed), [`ETH.III.SPARSEZ`](#p-eth-iii--eth-iii-sparsez), [`ETH.I.SPARS`](#p-eth-i--eth-i-spars), [`ETH.II.TANNER`](#p-eth-ii--eth-ii-tanner), [`ETH.II.GREEDY`](#p-eth-ii--eth-ii-greedy)

**Statement.** #k-SAT in (n+m)^{O_k(1)} 2^{O_k(n log log n/sqrt(log n))}, and 2^{O_{k,D}(n/sqrt(log n))} when m<=Dn; polynomial space.

**Proof sketch.** Fixed density is ETH.III.FIXED. In general: sparsify (ETH.I.SPARS) to polylog occurrence, then use Tanner tests and greedy rows (ETH.II.TANNER/GREEDY) with ETH.III.SPARSEZ as the zero counter.

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) Thm 10.1 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py) · *Note:* Not formalized.
<!-- END GENERATED: path-p-eth-iii -->

### Corollaries, new theorems, and the method's limits
<a id="p-eth-more"></a>

<!-- BEGIN GENERATED: path-p-eth-more -->
*What the counting algorithms give, the F_q extension, and the cost floor.*

Assumes: [`ETH.I.MAIN`](#p-eth-i--eth-i-main), [`ETH.II.MAIN`](#p-eth-ii--eth-ii-main), [`ETH.III.FIXED`](#p-eth-iii--eth-iii-fixed), [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift), [`ETH.II.ZEROS`](#p-eth-ii--eth-ii-zeros), [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount), [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer).

#### <a id="p-eth-more--eth-c-maxsat"></a>1. `ETH.C.MAXSAT` — Max-k-SAT approximation

*PROVED* · prose chain: closed · uses [`ETH.II.MAIN`](#p-eth-ii--eth-ii-main)

**Statement.** OPT(F)-eps m in 2^{O_{k,eps}(rho(n))} randomized time.

**Proof sketch.** Sample s = O(n/eps^2) clauses (Hoeffding plus a union bound over 2^n assignments). Maximize the sampled fraction exactly: a threshold circuit of size O(n), Tseitin, binary search on the threshold, and self-reduction, each step a 3-SAT decision in 2^{O(rho)}.

**Full proof:** [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §1

#### <a id="p-eth-more--eth-c-graph"></a>2. `ETH.C.GRAPH` — Independent sets, colourings, perfect matchings

*PROVED* · prose chain: closed · uses [`ETH.I.MAIN`](#p-eth-i--eth-i-main)

**Statement.** Counts in 2^{o(v)} (independent sets, q-colourings) and 2^{o(v+e)} (perfect matchings), polynomial space.

**Proof sketch.** Each count is #F for an explicit CNF: (not x_u v not x_w) per edge for independent sets; one-hot colour clauses plus edge-colour clauses (q-CNF on qv variables); for matchings, a Tseitin-encoded exactly-one circuit per vertex on edge variables, which preserves counts.

**Full proof:** [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §2

#### <a id="p-eth-more--eth-c-circ1"></a>3. `ETH.C.CIRC1` — E^NP not in SIZE(O(n))

*PROVED_CITED* · prose chain: cites ETH.C.CIRC1 · uses [`ETH.II.MAIN`](#p-eth-ii--eth-ii-main)

**Statement.** Williams's theorem applied to 3-SAT in 2^{o(n)}.

**Proof sketch.** Williams's connection: a 2^{o(n)} algorithm for 3-SAT (here ETH.II.MAIN) implies E^NP has no linear-size circuits. Cited, with hypotheses checked.

**Full proof:** [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §3, Cor. 3.1

#### <a id="p-eth-more--eth-c-circ2"></a>4. `ETH.C.CIRC2` — E^NP circuit lower bound c n sqrt(log n)

*PROVED_CITED* · prose chain: cites ETH.C.CIRC2 · uses [`ETH.III.FIXED`](#p-eth-iii--eth-iii-fixed)

**Statement.** Some L in E^NP has C_L(n)>c n sqrt(log n) infinitely often (needs the fixed-density bound; the formal bound alone gives n sqrt(log n)/(log log n)^2).

**Proof sketch.** Williams's witness argument with Zak's hierarchy and the Jahanjou-Miles-Viola local reduction. If the witness language had circuits of size eta n sqrt(log n), the test circuit (decoder + 3 copies) would have U_n = O(eta n sqrt(log n)) variables and <= 5U_n clauses. The fixed-density bound 2^{AU/sqrt(log U)} then decides H in NTIME(2^{n/2}), a contradiction. The rho(n) bound alone loses (log log n)^2.

**Full proof:** [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §3, Cor. 3.2, Remark 3.2

#### <a id="p-eth-more--eth-c-local"></a>5. `ETH.C.LOCAL` — Inverting maps of constant output locality

*PROVED* · prose chain: closed · uses [`ETH.II.MAIN`](#p-eth-ii--eth-ii-main)

**Statement.** |G^{-1}(b)| and a preimage in (n+m)^{O_d(1)} 2^{O_d(rho(n))}.

**Proof sketch.** Forbid each bad local assignment of each output bit: a width-d CNF with <= 2^d m clauses whose models are G^{-1}(b). Count it; find a preimage by self-reduction.

**Full proof:** [eth/11](reviews/claude/eth/11_GAPLESS_COROLLARIES.md) §4

#### <a id="p-eth-more--eth-n-liftp"></a>6. `ETH.N.LIFTP` — p-adic lift for every prime (E2) with a closed form

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift)

**Statement.** p_ell=F_ell(Y^{p-1}) and its h-adic digits (h=Y^p-Y) give S_ell(prod B_i mod V^ell)=prod [p|y_i] mod p^ell with degree d_\*(m(p-1)+p(ell-1)); closed form (1-(2z-1)(1+4 delta)^{-1/2})/2, z=Y^{p-1}, delta=Y^{p-2}h.

**Proof sketch.** z = Y^{p-1} is 0 or 1 on F_p, and delta = z^2 - z = Y^{p-2}h with h = Y^p - Y. The same square-root idempotent (1 - (2z-1)(1+4delta)^{-1/2})/2, truncated mod h^ell, is congruent to [p|y] mod p^ell (Hensel-type congruences at y = 0 and y a unit). Its h-adic digits by division by Y^p - Y give the cross-equation truncation with degree d_\*(m(p-1)+p(ell-1)), and the coefficients are bounded by (5p ell)^ell 8^ell.

**Full proof:** [eth/07](reviews/claude/eth/07_TOOLKIT.md) E2; [eth/12](reviews/claude/eth/12_FP_COUNTING.md) §1, Lemma 1.1, Remark 1.2 · **checks:** [verify_padic_truncation](reviews/claude/eth/scripts/verify_padic_truncation.py), [verify_fp_counting](reviews/claude/eth/scripts/verify_fp_counting.py)

#### <a id="p-eth-more--eth-n-fq"></a>7. `ETH.N.FQ` — Counting zeros of sparse systems over any finite field (new)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.N.LIFTP`](#p-eth-more--eth-n-liftp), [`ETH.II.ZEROS`](#p-eth-ii--eth-ii-zeros), [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount)

**Statement.** Fixed q: m<=cn polynomials over F_q in n variables, each with at most (log n)^c terms of degree <=k: common zeros counted in 2^{O(rho(n))} time, poly space.

**Proof sketch.** Encode each F_p variable in t = ceil(log p) bits, with validity factors [sum 2^b y_b < p]. P_i = 0 iff p | L_i for an integer lift. Lift with ETH.N.LIFTP, multiply in the V-truncated algebra (degree stays linear), and count by general-factor deformation counting with an auxiliary prime P > Q for the Chinese remainder step. F_q reduces to F_p by restriction of scalars.

**Full proof:** [eth/12](reviews/claude/eth/12_FP_COUNTING.md) Thm 12.1, Cor. 12.2 · **checks:** [verify_fp_counting](reviews/claude/eth/scripts/verify_fp_counting.py)

#### <a id="p-eth-more--eth-k-e3"></a>8. `ETH.K.E3` — Exact partition functions of bounded-degree factor graphs

*NEW, PROVED* · prose chain: closed · uses [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount), [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer)

**Statement.** Weighted factor graphs with rational weights: exact partition function in 2^{O(n(log log n)^2/sqrt(log n))}.

**Proof sketch.** Deformation counting never uses that the factors are clauses. Take multilinear interpolants of rational weight tables, run Part I with Part II parameters, and clear denominators.

**Full proof:** [eth/07](reviews/claude/eth/07_TOOLKIT.md) E3

#### <a id="p-eth-more--eth-k-e4"></a>9. `ETH.K.E4` — Cost floor of the method

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.SWAP`](#p-eth-i--eth-i-swap)

**Statement.** For symmetric designs, the method's exponent cannot beat (log log n)^2/sqrt(log n).

**Proof sketch.** The three costs are: delta >= 1/r (largest group); delta >= Q eps_\*/32 (measurement rows at the worst level inside the window); delta >= x log(1/x) with x = R/Q (Taylor jets). The method of types gives eps_\* >= e^{-O(R^2/Q)}/(Q+1), forcing R^2 >~ Q log(1/delta). Together with Q <~ log n/log r and log r >~ log(1/delta), this forces x >~ log log n/sqrt(log n), hence delta >~ (log log n)^2/sqrt(log n).

**Full proof:** [eth/07](reviews/claude/eth/07_TOOLKIT.md) E4 · **checks:** [verify_cost_floor](reviews/claude/eth/scripts/verify_cost_floor.py)
<!-- END GENERATED: path-p-eth-more -->

### One deformation theorem for Parts I–III (13)
<a id="p-eth-unified"></a>

<!-- BEGIN GENERATED: path-p-eth-unified -->
*The common skeleton of swap interpolation and App. D: readout systems, the filtered node lemma, Hermite recovery with loss, correction counts; both instances; both square-root floors.*

Assumes: [`ETH.I.SWAP`](#p-eth-i--eth-i-swap), [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack).

#### <a id="p-eth-unified--eth-u-node"></a>1. `ETH.U.NODE` — Node lemma: filtered mode expansion

*NEW, PROVED, CHECKED* · prose chain: closed

**Statement.** Readout system Z(rho) = sum_x prod_f Phi_f(rho(x)). If rho = lambda + sum_Z kappa_Z L_Z mod a_P with lambda local (modulo measurements), kappa_Z in a_{w_Z} and L_Z reading only its source group, then Z(rho) = sum_y sum over correction indices of weight < P of products of group sums: (targeted modes) beta with sum w beta < P and symbolic coefficient extraction; (Boolean coordinate modes) pinned sets E with sum w < P and Moebius points U subset E.

**Proof sketch.** Replace rho by lambda + sum kappa L (Phi_f has coefficients in the ring). Expand prod_f Phi_f(lambda + sum kappa zeta) in formal modes zeta: the coefficient of zeta^beta lies in a_{sum w beta}, so only light monomials survive. Targeted: each group's factors see only modes targeted at it, so coefficients factor by group and L^beta by source. Boolean: multilinearize, extract [zeta^E] by Moebius inversion at the points 1_U, and the pinning x_E = 1 replaces the monomial. In both cases separation by measurements factorizes the sum over groups.

**Full proof:** [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Lemma 2.3, Lemma 2.4, Lemma 1.2, Definition 2.2 · **checks:** [verify_unified_deformation](reviews/claude/eth/scripts/verify_unified_deformation.py) · *Note:* Part I uses the formal (z-adic) neighbourhood with all coordinates as Boolean modes; Part III the 2-adic neighbourhood with structured tile modes. The script runs one generic implementation on genuine instances of both.

#### <a id="p-eth-unified--eth-u-hermite"></a>2. `ETH.U.HERMITE` — Hermite recovery with valuation loss

*NEW, PROVED, CHECKED* · prose chain: closed

**Statement.** Nodes t in a discretely valued field with multiplicities m_t and scalings alpha_t (v(alpha_t) >= v(t-q)): F(0) = sum_t E_t(0) sum_{j&lt;m_t} [z^j](F(t+alpha_t z)/E_t(t+alpha_t z)) (-t/alpha_t)^j, and knowing the jets mod varpi^{P_t}, P_t = P + l_t + (m_t-1) max(0,-v(t/alpha_t)), l_t = -v(E_t(0)), determines F(0) mod varpi^P; one digit less fails. For the residue system mod 2^{a+1} with a depth profile m: l_t = D_i Delta(i), Delta(i) = sum_j (j/2^j)(m(i)-m(i+j)).

**Proof sketch.** Hermite: G(T) = sum_t E_t(T) P_t((T-t)/alpha_t) agrees with F to order m_t at every node and has degree < sum m_t, so G = F. E_t(t+alpha_t z) = prod (1 + alpha_t z/(t-q))^{m_q} is integral with constant term 1. A change delta at coefficient j propagates with valuation >= P_t - l_t - (m_t-1)max(0,-v(t/alpha_t)); at the top coefficient the bound is attained. The loss l_t = sum_q m_q (v(t-q) - v(q)) counts how much closer the other nodes are to t than to 0: zero for unit-distance nodes (Part I), D_i Delta(i) for residue systems (Part III).

**Full proof:** [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Lemma 3.1, Cor. 3.2, Cor. 3.3 · **checks:** [verify_unified_deformation](reviews/claude/eth/scripts/verify_unified_deformation.py)

#### <a id="p-eth-unified--eth-u-count"></a>3. `ETH.U.COUNT` — Correction counts, two-sided

*NEW, PROVED, CHECKED* · prose chain: closed

**Statement.** With c_w modes of weight w: #{beta : sum w beta < P} <= inf_theta e^{theta(P-1)} prod_w (1-e^{-theta w})^{-c_w} <= exp(pi sqrt(2cP/3)) if c_w <= c; and log_2 # >= (1/8) min(P-1, sqrt(c(P-1))) if c_w >= c for all w < P. Boolean: #{U subset E : weight(E) < P} <= inf_theta e^{theta(P-1)} prod_Z (1 + 2 e^{-theta w_Z}) (= C_{P-1}(n) for n weight-1 modes).

**Proof sketch.** Rankin's trick with sum_w -log(1-e^{-theta w}) <= pi^2/(6 theta). Lower bound: if c >= (P-1)/8 take subsets of weight-1 modes; otherwise all subsets of the modes of weight <= sqrt(2(P-1)/c) - 1 have total weight < P, giving 2^{c W} >= 2^{0.7 sqrt(c(P-1))}.

**Full proof:** [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Lemma 4.1 · **checks:** [verify_unified_deformation](reviews/claude/eth/scripts/verify_unified_deformation.py)

#### <a id="p-eth-unified--eth-u-main"></a>4. `ETH.U.MAIN` — Interpolation through local nodes (the unified deformation theorem)

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.U.NODE`](#p-eth-unified--eth-u-node), [`ETH.U.HERMITE`](#p-eth-unified--eth-u-hermite), [`ETH.U.COUNT`](#p-eth-unified--eth-u-count)

**Statement.** Readout system, deformation rho(T) with rho(0) = incidence readout (degree R', deformed-port degree d), recovery data with sum m_t > dR', and at each node a mode expansion (formal, or varpi-adic at more than dR' unit-separated points): Z is computed mod varpi^P by node evaluation, interpolation and Hermite recovery, with work sum_t |Xi_t| Y_t N_t sum_g W_{t,g} + poly. Parts I-II (formal, unit nodes, Boolean modes, brute-force groups) and Part III (2-adic, residue system, tile modes, recursive groups) are corollaries.

**Proof sketch.** H(T) = Z(rho(T)) has degree <= dR' and H(0) = Z (substitution homomorphism, degree measured in the coefficient algebra). The node lemma gives H near each node: directly as a jet (formal case) or at points t + alpha zeta, interpolated in zeta over a Vandermonde unit (varpi-adic case). Hermite recovery with loss returns H(0). The cost is the number of correction indices times measurement values times group costs, per node and point.

**Full proof:** [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Theorem 4.2, Cor. 4.3, Cor. 4.4, Definition 0.3 · **checks:** [verify_unified_deformation](reviews/claude/eth/scripts/verify_unified_deformation.py)

#### <a id="p-eth-unified--eth-i-defcount"></a>5. `ETH.I.DEFCOUNT` — Deformation counting

*PROVED, CHECKED* · prose chain: closed

**Statement.** Given a deformation of degree <=R' with Q&lt;p nodes whose rows are local for measured partitions, the sum over {0,1}^n of a product of factors (clauses, or any factors with values in a free Z/p^mu-algebra, degree measured in the algebra) is computed with Q\*C_q(n)\*(2beta+1)^{c_\*}\*2r\*2^s\*poly operations, q=floor(dR'/Q).

**Proof sketch.** H(T) = sum_x prod_f Phi_{C_f}(&lt;P_{f,v}(T),x>) has degree <= dR' and H(0) = #F. At each node t, the jet of H(t+z) mod z^J is a sum over x of polynomials of degree <= q in x. Moebius inversion on the Boolean lattice reduces it to Boolean substitutions over sets E with |E| <= q (C_q(n) of them). Locality of the rows makes each sum factor over the groups once the measurement values y in [-beta,beta]^{c_\*} are fixed (separation by measurements). The nodes 1..Q have unit differences mod p > Q, so the Chinese remainder theorem recovers H(0) mod p^mu.

**Alternative proof.** Uses [`ETH.U.MAIN`](#p-eth-unified--eth-u-main). As the formal-neighbourhood instance of the unified deformation theorem: unit nodes, Boolean coordinate modes, measurements. ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 4.3)

**Full proof:** [eth/08](reviews/claude/eth/08_GAPLESS_PROOF.md) §4, Thm 4.6; [eth/09](reviews/claude/eth/09_GAPLESS_SHARP_EXPONENT.md) §1, Thm 1.1 (general factors) · **checks:** [verify_gapless_chain](reviews/claude/eth/scripts/verify_gapless_chain.py), [verify_sharp_route](reviews/claude/eth/scripts/verify_sharp_route.py)

#### <a id="p-eth-unified--eth-iii-child"></a>6. `ETH.III.CHILD` — Modes and the child identity

*PROVED, CHECKED* · prose chain: closed · uses [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack)

**Statement.** At T=t+alpha zeta, crossing entries expand in tile moments sum_j j^b X_y with 2-adic weight b+1, and F(T)=sum_{beta: sum w beta&lt;P} prod_g ct(N_{g,beta}) mod 2^P, with children of <=6M/r edges and exp(2 pi sqrt(MP/D_t)) mode tuples.

**Proof sketch.** At a node t of depth i, an entry W_{xy} for a crossing tile y_0 is a geometric series in j (y = y_0 + 2^{i+1}j) with ratio of valuation >= 1. So rows become lambda_x(X) + sum_Z kappa_{x,Z} L_Z(X), with modes L_Z = sum_j j^b X_y of weight b+1. Expanding in formal modes, the coefficient of Z^beta has valuation >= sum w beta, so only tuples with sum w beta < P survive. For each tuple the summand factors over the groups r into children with <= 6M/r edges. Rankin's bound counts the tuples by exp(2 pi sqrt(MP/D_t)).

**Alternative proof.** Uses [`ETH.U.NODE`](#p-eth-unified--eth-u-node), [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy), [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack). As the targeted node lemma with tile power moments as modes (weights b+1, source = tile group, target = crossing rows). ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 2.7)

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §6, Lemmas 6.1-6.4 · **checks:** [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-unified--eth-iii-hermite"></a>7. `ETH.III.HERMITE` — Tapered 2-adic Hermite recovery

*PROVED, CHECKED* · prose chain: closed

**Statement.** With multiplicities q_t=h max(0,u-v_2(t)), F(0) mod 2^p is recovered from F(t+alpha z) mod 2^{p_t}, p_t=p+8hD_t+8M(a+1); loss v_2(E_t(0))>=-6hD_t.

**Proof sketch.** Hermite interpolation F(0) = sum_t E_t(0) P_t(-t/alpha), with node multiplicities q_t = h(u - v_2(t))^+ tapered by 2-adic depth (sum q_t > deg F). The 2-adic valuation of E_t(0) is >= -6hD_t (a geometric sum over depth levels), so knowing each jet mod 2^{p_t} with p_t = p + 8hD_t + 8M(a+1) determines F(0) mod 2^p.

**Alternative proof.** Uses [`ETH.U.HERMITE`](#p-eth-unified--eth-u-hermite). As Hermite recovery with valuation loss on the residue system mod 2^{a+1}: loss D_i Delta(i) <= 6 h D_i for the taper. ([eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Cor. 3.3)

**Full proof:** [eth/10](reviews/claude/eth/10_GAPLESS_TENSOR_CONTRACTION.md) §4, Lemma 4.1; [eth/07](reviews/claude/eth/07_TOOLKIT.md) E5 · **checks:** [verify_appendix_d](reviews/claude/eth/scripts/verify_appendix_d.py), [verify_tensor_contraction](reviews/claude/eth/scripts/verify_tensor_contraction.py)

#### <a id="p-eth-unified--eth-k-e4"></a>8. `ETH.K.E4` — Cost floor of the method

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.I.SWAP`](#p-eth-i--eth-i-swap)

**Statement.** For symmetric designs, the method's exponent cannot beat (log log n)^2/sqrt(log n).

**Proof sketch.** The three costs are: delta >= 1/r (largest group); delta >= Q eps_\*/32 (measurement rows at the worst level inside the window); delta >= x log(1/x) with x = R/Q (Taylor jets). The method of types gives eps_\* >= e^{-O(R^2/Q)}/(Q+1), forcing R^2 >~ Q log(1/delta). Together with Q <~ log n/log r and log r >~ log(1/delta), this forces x >~ log log n/sqrt(log n), hence delta >~ (log log n)^2/sqrt(log n).

**Full proof:** [eth/07](reviews/claude/eth/07_TOOLKIT.md) E4 · **checks:** [verify_cost_floor](reviews/claude/eth/scripts/verify_cost_floor.py)

#### <a id="p-eth-unified--eth-u-floor3"></a>9. `ETH.U.FLOOR3` — Cost floor of the 2-adic architecture

*NEW, PROVED, CHECKED* · prose chain: closed · uses [`ETH.U.HERMITE`](#p-eth-unified--eth-u-hermite), [`ETH.U.COUNT`](#p-eth-unified--eth-u-count)

**Statement.** In App. D's architecture (Cauchy-kernel deformation, nodes S\{0}, any depth profile on depths < a, formula-based Hermite recovery at output precision p, tile modes with >= M/2 crossing wires), some node needs >= 2^{(1/8) min(p-1, sqrt(3/16) M/sqrt(a))} correction tuples; with a <= log_2 M and p >= M/C at the top level this is 2^{Omega(M/sqrt(log M))}, matching Thm D.1.

**Proof sketch.** Depth-free node cost: a depth-i node has >= M/D_i modes per weight and needs precision >= p + D_i Delta(i), so modes x precision >= M Delta(i). Profile lemma: Hermite needs sum 2^{a-i} m(i) > M(2L-1), and a downward induction gives m(i) <= 2 Delta\*(a-i), so Delta\* > 3M/(8a). The lower bound of ETH.U.COUNT gives the count. The exact optimum of the profile problem and the paper's linear taper both tend to 6M/a.

**Full proof:** [eth/13](reviews/claude/eth/13_UNIFIED_DEFORMATION.md) Theorem 5.4, Lemma 5.2, Lemma 5.3 · **checks:** [verify_unified_deformation](reviews/claude/eth/scripts/verify_unified_deformation.py) · *Note:* Together with ETH.K.E4 (Part I's floor), both sqrt(log) barriers are theorems about their architectures. Whether every local-node scheme pays M/sqrt(log M) is target T.ETH.5.
<!-- END GENERATED: path-p-eth-unified -->

## Index of nodes

Each node links to its first appearance on a path.

<!-- BEGIN GENERATED: node-index -->
| Node | Title | Status | On paths |
|---|---|---|---|
| [`ETH.I.SPARS`](#p-eth-i--eth-i-spars) | Disjoint sparsification | PROVED, CHECKED | P.ETH.I |
| [`ETH.I.TOVEY`](#p-eth-i--eth-i-tovey) | OR gates and Tovey cycles | PROVED, CHECKED | P.ETH.I |
| [`ETH.I.DEFCOUNT`](#p-eth-i--eth-i-defcount) | Deformation counting | PROVED, CHECKED | P.ETH.I, P.ETH.UNIFIED |
| [`ETH.I.SWAP`](#p-eth-i--eth-i-swap) | Swap interpolation in closed form | PROVED, CHECKED | P.ETH.I |
| [`ETH.I.TRANSFER`](#p-eth-i--eth-i-transfer) | Transfer to a formula; balanced schedules | PROVED, CHECKED | P.ETH.I |
| [`ETH.I.PARAM`](#p-eth-i--eth-i-param) | Counting at &lt;=3 occurrences, explicit parameters | PROVED, CHECKED | P.ETH.I |
| [`ETH.I.MAIN`](#p-eth-i--eth-i-main) | ETH is false (Part I) | PROVED, CHECKED | P.ETH.I |
| [`ETH.II.LIFT`](#p-eth-ii--eth-ii-lift) | Parity lift in closed form; cross-equation truncation | PROVED, CHECKED | P.ETH.II |
| [`ETH.II.ZEROS`](#p-eth-ii--eth-ii-zeros) | Counting common zeros of sparse F_2 systems | PROVED, CHECKED | P.ETH.II |
| [`ETH.II.TANNER`](#p-eth-ii--eth-ii-tanner) | Initial tests from a Tanner code | PROVED, CHECKED | P.ETH.II |
| [`ETH.II.GREEDY`](#p-eth-ii--eth-ii-greedy) | Greedy rows by exact counts | PROVED, CHECKED | P.ETH.II |
| [`ETH.II.MAIN`](#p-eth-ii--eth-ii-main) | #k-SAT in 2^{O(rho(n))} (formally verified bound) | FORMAL, PROVED, CHECKED | P.ETH.II |
| [`ETH.III.RING`](#p-eth-iii--eth-iii-ring) | Coefficient ring and interpolation points | PROVED | P.ETH.III |
| [`ETH.III.CAUCHY`](#p-eth-iii--eth-iii-cauchy) | Cauchy kernel | PROVED, CHECKED | P.ETH.III |
| [`ETH.III.HERMITE`](#p-eth-iii--eth-iii-hermite) | Tapered 2-adic Hermite recovery | PROVED, CHECKED | P.ETH.III, P.ETH.UNIFIED |
| [`ETH.III.PACK`](#p-eth-iii--eth-iii-pack) | Packing edges into packets | PROVED | P.ETH.III |
| [`ETH.III.CHILD`](#p-eth-iii--eth-iii-child) | Modes and the child identity | PROVED, CHECKED | P.ETH.III, P.ETH.UNIFIED |
| [`ETH.III.TENSOR`](#p-eth-iii--eth-iii-tensor) | Boolean tensor contraction | PROVED, CHECKED | P.ETH.III |
| [`ETH.III.FIXED`](#p-eth-iii--eth-iii-fixed) | Fixed clause density | PROVED, CHECKED | P.ETH.III |
| [`ETH.III.SPARSEZ`](#p-eth-iii--eth-iii-sparsez) | Sparse systems by contraction | PROVED, CHECKED | P.ETH.III |
| [`ETH.III.MAIN`](#p-eth-iii--eth-iii-main) | Headline Theorem 1.1 | PROVED, CHECKED | P.ETH.III |
| [`ETH.C.MAXSAT`](#p-eth-more--eth-c-maxsat) | Max-k-SAT approximation | PROVED | P.ETH.MORE |
| [`ETH.C.GRAPH`](#p-eth-more--eth-c-graph) | Independent sets, colourings, perfect matchings | PROVED | P.ETH.MORE |
| [`ETH.C.CIRC1`](#p-eth-more--eth-c-circ1) | E^NP not in SIZE(O(n)) | PROVED_CITED | P.ETH.MORE |
| [`ETH.C.CIRC2`](#p-eth-more--eth-c-circ2) | E^NP circuit lower bound c n sqrt(log n) | PROVED_CITED | P.ETH.MORE |
| [`ETH.C.LOCAL`](#p-eth-more--eth-c-local) | Inverting maps of constant output locality | PROVED | P.ETH.MORE |
| [`ETH.N.LIFTP`](#p-eth-more--eth-n-liftp) | p-adic lift for every prime (E2) with a closed form | NEW, PROVED, CHECKED | P.ETH.MORE |
| [`ETH.N.FQ`](#p-eth-more--eth-n-fq) | Counting zeros of sparse systems over any finite field (new) | NEW, PROVED, CHECKED | P.ETH.MORE |
| [`ETH.K.E3`](#p-eth-more--eth-k-e3) | Exact partition functions of bounded-degree factor graphs | NEW, PROVED | P.ETH.MORE |
| [`ETH.K.E4`](#p-eth-more--eth-k-e4) | Cost floor of the method | NEW, PROVED, CHECKED | P.ETH.MORE, P.ETH.UNIFIED |
| [`ETH.U.NODE`](#p-eth-unified--eth-u-node) | Node lemma: filtered mode expansion | NEW, PROVED, CHECKED | P.ETH.UNIFIED |
| [`ETH.U.HERMITE`](#p-eth-unified--eth-u-hermite) | Hermite recovery with valuation loss | NEW, PROVED, CHECKED | P.ETH.UNIFIED |
| [`ETH.U.COUNT`](#p-eth-unified--eth-u-count) | Correction counts, two-sided | NEW, PROVED, CHECKED | P.ETH.UNIFIED |
| [`ETH.U.MAIN`](#p-eth-unified--eth-u-main) | Interpolation through local nodes (the unified deformation theorem) | NEW, PROVED, CHECKED | P.ETH.UNIFIED |
| [`ETH.U.FLOOR3`](#p-eth-unified--eth-u-floor3) | Cost floor of the 2-adic architecture | NEW, PROVED, CHECKED | P.ETH.UNIFIED |
| [`KAK.Q.CALC`](#p-kak-planar--kak-q-calc) | Quaternionic slot calculus | PROVED, CHECKED | P.KAK.PLANAR, P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.C.LIFT`](#p-kak-planar--kak-c-lift) | Lifted current and push-forward | PROVED, CHECKED | P.KAK.PLANAR, P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.C.WAVE`](#p-kak-paper--kak-c-wave) | The positive wave current | PROVED, CHECKED | P.KAK.PAPER |
| [`KAK.F.GEOM`](#p-kak-paper--kak-f-geom) | Fibre geometry | PROVED | P.KAK.PAPER |
| [`KAK.DIR`](#p-kak-engine--kak-dir) | Periodic Dirichlet problem (3-slot quaternionic Monge-Ampere) | PROVED_CITED | P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.E.BAL`](#p-kak-paper--kak-e-bal) | Determinant energy balance | PROVED, CHECKED | P.KAK.PAPER |
| [`KAK.E.DEPTH`](#p-kak-paper--kak-e-depth) | Density to depth (rank-r background) | NEW, PROVED | P.KAK.PAPER |
| [`KAK.FIB`](#p-kak-paper--kak-fib) | Energy estimates on a fibre | PROVED | P.KAK.PAPER |
| [`KAK.PROP`](#p-kak-planar--kak-prop) | Two-circle propagation | PROVED, CHECKED | P.KAK.PLANAR, P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.CAP`](#p-kak-engine--kak-cap) | Capacity-dimension principle (new) | NEW, PROVED, CHECKED | P.KAK.ENGINE |
| [`KAK.T12`](#p-kak-engine--kak-t12) | Critical annular estimate (Theorem 1.2) | PROVED | P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.T11`](#p-kak-engine--kak-t11) | Sogge's local smoothing conjecture (Theorem 1.1) | FORMAL, PROVED | P.KAK.ENGINE, P.KAK.PAPER |
| [`KAK.R.SPHERE`](#p-kak-conseq--kak-r-sphere) | Spherical restriction, interior range | FORMAL, PROVED_CITED | P.KAK.CONSEQ |
| [`KAK.R.CONE`](#p-kak-conseq--kak-r-cone) | Truncated cone restriction | PROVED | P.KAK.CONSEQ |
| [`KAK.R.PARAB`](#p-kak-conseq--kak-r-parab) | Paraboloid restriction and Schrodinger local smoothing | PROVED_CITED | P.KAK.CONSEQ |
| [`KAK.KM`](#p-kak-conseq--kak-km) | Kakeya maximal conjecture | FORMAL, PROVED | P.KAK.CONSEQ |
| [`KAK.KD`](#p-kak-conseq--kak-kd) | Kakeya sets have full Hausdorff dimension | FORMAL, PROVED | P.KAK.CONSEQ |
| [`KAK.NIK`](#p-kak-conseq--kak-nik) | Nikodym maximal conjecture | PROVED | P.KAK.CONSEQ |
| [`KAK.BR`](#p-kak-conseq--kak-br) | Bochner-Riesz conjecture (with maximal bounds for p&gt;=2) | FORMAL, PROVED_CITED | P.KAK.CONSEQ |
| [`KAK.SQ`](#p-kak-conseq--kak-sq) | Stein's square function conjecture | FORMAL, PROVED_CITED | P.KAK.CONSEQ |
| [`KAK.OB.FOCUS`](#p-kak-limits--kak-ob-focus) | Radial focusing (sharpness of s_+(p) and of the power 1/p_c) | PROVED, CHECKED | P.KAK.LIMITS |
| [`KAK.OB.WEAK`](#p-kak-limits--kak-ob-weak) | Weak-norm obstruction | PROVED, CHECKED | P.KAK.LIMITS |
| [`KAK.ANG.HANKEL`](#p-kak-angular--kak-ang-hankel) | Uniform Hankel expansion on an angular band | NEW, PROVED, CHECKED | P.KAK.ANGULAR |
| [`KAK.ANG.BALL`](#p-kak-angular--kak-ang-ball) | Observation on a ball (shell transplant) | PROVED, CHECKED | P.KAK.ANGULAR |
| [`KAK.ANG.ENV`](#p-kak-angular--kak-ang-env) | Four-kernel envelope and small source | PROVED, CHECKED | P.KAK.ANGULAR |
| [`KAK.ANG.BAND`](#p-kak-angular--kak-ang-band) | One annulus, bounded angular degree | PROVED, CHECKED | P.KAK.ANGULAR |
| [`KAK.APP.B`](#p-kak-angular--kak-app-b) | Angular endpoint | PROVED_CITED, CHECKED | P.KAK.ANGULAR |
| [`KAK.APP.D`](#p-kak-planar--kak-app-d) | Planar case by a linear Dirichlet problem (Theorem 13.1) | PROVED, CHECKED | P.KAK.PLANAR |
| [`KAK.APP.E`](#p-kak-engine--kak-app-e) | Curved characteristic sheets of symmetric systems | PROVED, CHECKED | P.KAK.ENGINE |
| [`KAK.ENGINE`](#p-kak-engine--kak-engine) | The positivity engine (Theorem 11.1) | NEW, PROVED | P.KAK.ENGINE |
| [`KAK.PENCIL`](#p-kak-planar--kak-pencil) | Symmetric pencils carry engine data (Theorem 11.2) | NEW, PROVED, CHECKED | P.KAK.PLANAR, P.KAK.ENGINE |
| [`KAK.ARROW`](#p-kak-engine--kak-arrow) | The wave is the arrow pencil (Corollary 11.4) | NEW, PROVED, CHECKED | P.KAK.ENGINE |
| [`KAK.FRAME`](#p-kak-planar--kak-frame) | Cofactor divergence for arbitrary quaternionic frames (Theorem 12.1) | NEW, PROVED, CHECKED | P.KAK.PLANAR, P.KAK.ENGINE |
| [`KAK.MAENERGY`](#p-kak-engine--kak-maenergy) | The Monge-Ampere operator is the gradient of a concave energy (Theorem 12.3) | NEW, PROVED | P.KAK.ENGINE |
<!-- END GENERATED: node-index -->
