# 4. Does the paper rely on miracles? What makes it work

## 4.1 Summary

There is **no unexplained coincidence in the load-bearing chain.** Every identity in §2–§4 has a short algebraic reason, and I verified each symbolically or numerically ([02](02_VERIFICATION_LEDGER.md)). The paper's own explanations are compressed but adequate. The design choices whose motivation it leaves implicit are explained below.

| Apparent miracle | Real reason | Paper explains? |
|---|---|---|
| A *positive* conserved tensor exists for waves | Square after summing, in lifted variables; CR vector fields annihilate the analytic continuation (§4.2 below) | Yes (§1.4, §2.2), tersely |
| Quaternions appear | Two auxiliary flows → two complex rank-one terms → one quaternionic rank-one term (§4.3) | Partly: the Gram computation is given, the reason is not |
| The exponent is exactly $p_c=2n/(n-1)$ | The Monge–Ampère determinant has degree $n$, and one direction is taken by the background (§4.4) | Yes (§1.4), well |
| The full determinant cancels in $\Delta_z\mathcal F$ | Aubin–Mabuchi first variation plus cofactor conservation (§4.5) | Yes, by analogy with Kähler geometry |
| The cofactor boundary term vanishes without a sign | The null line of $H(r)$ lies in $\ker B\cap\nu^\perp$, by fiber design (§4.5) | Yes (Lemma 4.1, Prop 4.2) |
| The potential $V$ is harmless | Gaussian decay at the wall beats Gaussian inversion at the center (§4.6) | Yes (Fig. 1) |
| Exactly two logarithms | Two harmonic-measure passages through a thin cylinder; sharp for scalar supersolutions | Yes (Fig. 2, Remark 6.2) |
| The holomorphic flux map $\Theta$ and the passive variable $\eta$ | $\Theta$ is the complexified flux $z\mapsto\frac12z^tC^\mu z$; $\eta$ completes the dimension count | Stated, motivation implicit (§4.3) |

## 4.2 Why a positive conserved tensor exists

**Setup.**
- $\mathscr P=\sum_\mu C^\mu\partial_\mu$ is a constant-coefficient real symmetric system.
- $f$ is a holomorphic solution of $\mathscr Pf=0$, the analytic continuation of a finite superposition of plane waves.
- For auxiliary real states $w_1,w_2\in\mathbb R^n$, $q^\mu=\frac12\sum_\ell w_\ell^tC^\mu w_\ell$ is their combined flux.
- $G=\overline{f(\mathsf X-\mathrm iq)}$.

**The CR fields.** The fields $D_{\ell a}=\partial_{w_{\ell a}}+\mathrm iV_{\ell a}$, with $V_{\ell a}=\sum_\mu(C^\mu w_\ell)_a\partial_\mu$, annihilate every holomorphic function of $\mathsf X-\mathrm iq$. Hence:
- $\mathcal Z_\ell=\sum_aG_aD_{\ell a}$ has divergence $2\mathrm iw_\ell^t\mathscr PG=0$;
- the "acceleration" $\nabla_{\mathcal Z_\ell}\overline{\mathcal Z_m}$ is $\mathrm i$ times a real field.

So $\operatorname{Re}\sum_\ell\mathcal Z_\ell\otimes\overline{\mathcal Z_\ell}$ is a sum of real squares with zero column divergence.

**Why squaring after summing matters.** Summing modes before squaring keeps the interference, i.e. the focusing. Squaring each mode would destroy it. That is why $A_{TT}$ detects $\|U(T)f\|_{p_c}$ and not merely $\|f\|_2$.

## 4.3 Why quaternions

**Two flows are needed.** One auxiliary state would give a Gaussian in only one real direction of the transverse frequency. The weight (2.4) needs the full complex modulus $|\zeta-vR|^2$, so it needs *two* real states $w_1,w_2$.

**Two flows give a quaternion.** The two complex rank-one squares $\mathcal Z_1\otimes\bar{\mathcal Z}_1+\mathcal Z_2\otimes\bar{\mathcal Z}_2$ combine into $Q=\mathcal Z_1+\mathbf k\mathcal Z_2$, with real part $\operatorname{Re}(Q\bar Q^t)$. The Gram matrix is quaternionic of rank one.

**Where the slots come from.** The holomorphic flux map assigns each spacetime coordinate a triple $X^\mu+\omega_\mu\mathbf i$ with $\omega_\mu\in\mathbb C_{\mathbf k}$, i.e. an element of $\operatorname{span}(1,\mathbf i,\mathbf j)$. This gives a "3-slot" quaternionic structure. The quaternionic Hessian $H(u)=\sum e_\alpha\bar e_\beta\partial_{a\alpha}\partial_{b\beta}u$ is then the natural second-order operator to pair with $A$, because $A:H(u)=A^{\mathbb R}:D^2u$.

**The fourth slot** is added only in App. A, to obtain the quaternionic rigid motions used in the boundary $C^2$ estimate.

## 4.4 Why exactly $p_c$: the rank-one background

Let $\beta\approx\lambda$ be the single large eigenvalue of the background Hessian, and let $K$ be any further Hessian. By (3.7),
$$\det(\beta vv^*+K)=\det K+\beta\,\det(K|_{v^\perp}).$$
For $\beta$ huge, the equation $\det=\lambda\mu$ forces $\det(K|_{v^\perp})\approx\mu$, a Monge–Ampère equation in $d=n-1$ directions.

On a ball of radius $h$ where the density is about $m$, the comparison quadratic has curvature about $m^{1/d}$, so the depth is about $h^2m^{1/d}$. Hence
$$(-u)\mu\approx h^2m^{1+1/d},\qquad 2\Big(1+\frac1d\Big)=\frac{2n}{n-1}=p_c .$$
Since $\mu$ is comparable to $|\text{wave}|^2$, this is $|\text{wave}|^{p_c}$.

The upper bound uses a barrier $\propto rP_*$, where $P_*$ is the $d$-th root of a smoothed density. Young's inequality with exponent $1+1/d$ closes it. Without the background, the same equation would give $\mu^{1/n}$ and the exponent $2(n+1)/n<p_c$, which has no special meaning for the wave. The background's role is to remove exactly one degree from the determinant.

**One determinant controls the exponent.** No choice of $p$ is made anywhere; the degree of the determinant makes it. This is also why the argument cannot be pushed above $p_c$, where Prop 8.1 shows the estimate is false.

## 4.5 Why the determinant cancels and the boundary is harmless

**Setup.** Write $D_v=H_{\rm fib}(B_0+v)$ and define
$$E_z(v)=\frac1\lambda\int_{\Omega_z}v\int_0^1\det D_{\theta v}\,d\theta.$$
This is a quaternionic Aubin–Mabuchi energy.

**The determinant cancels.** Its parameter Laplacian is $\lambda^{-1}\int\det H(B_0+u)$, the *full* determinant including the $z$-directions. Two facts give this:
- the full cofactor $C(H(B_0+\theta u))$ is divergence-free (Lemma 3.3);
- $\partial_\theta\det=C:H(u)$.

The Laplacian of $X_0=\int(-u)\mu$ contains $-\lambda^{-1}\int\det H(B_0+u)$, via the Schur split (4.9). The two cancel.

**The boundary term vanishes.** The normal cofactor coefficient is $\nu^*C\nu=\det(H|_{\nu^\perp})$. On the wall, $H(B_0+\theta u)|_{\nu^\perp}=(B+\theta\rho H(r))|_{\nu^\perp}$. The fiber is designed so that the null line of $H(r)$ lies in $\nu^\perp\cap\ker B$. The compression is therefore singular and the term vanishes, whatever the sign of the indefinite full Hessian.

## 4.6 Why the potential is harmless

The energy's Laplacian is bounded by the wave flux through the wall, where $|W|\ge\sqrt3w_0/2$. The data live in a cap $|v|\le\gamma\ll w_0$, so the Gaussian $e^{-\frac12\kappa y|W-v|^2}$ gives $e^{-c_0s/h}$ there.

To bound the energy *below*, one inverts the Gaussian near the fiber center, which costs $e^{C\gamma^2s/h}$. Choosing $\gamma$ small makes the net factor $e^{-cs/h}$. The resulting potential $V(s)=Cs^{-2}e^{-cs/h}$ has $\int s\log(s/b)V\,ds=O(1)$, uniformly in $h$. This is exactly the hypothesis of Lemma 6.1.

## 4.7 Where the paper is least explained

1. **The holomorphic flux map $\Theta$** (2.7)–(2.9). It is presented as a definition; its role (turning Gram rows into $-\mathbf iQ,-\mathbf jQ$) is verified after the fact. A reader would benefit from the statement "Θ is the unique quadratic holomorphic map whose differential sends Gram rows to the desired slots".
2. **The fiber design (4.1)**: the ellipsoid in $(W,Q)$ and $Q=\Phi/a$. Its three roles are scattered across §4–§5:
   - tangency of the current to $\Phi$-levels;
   - the null line of $H(r)$;
   - the wall being Gaussian-small.
3. **App. A's fourth slot.** It is needed only for boundary estimates. A remark on whether 3-slot rigid motions suffice would help.

None of these is a gap; they are exposition debts. [05](05_DISTILLATION.md) repackages the argument so that each design choice appears as the hypothesis of an abstract theorem.

## 4.8 Sanity checks against known facts

- **Exponent.** The method yields exactly $p_c$, where the focusing example (Prop 8.1) shows the log is necessary. Above $p_c$ the paper interpolates with $L^\infty$, which is consistent with the sharp Sobolev order.
- **Planar case.** For $n=2$ the fiber problem is *linear* (App. D). This recovers sharp local smoothing at $p=4$ (Guth–Wang–Zhang) with a $\log^{1/2}$ loss. That is a known theorem, now reproved by an independent route, and a strong consistency check.
- **Variable coefficients.** Known negative results for curved metrics (Minicozzi–Sogge) and for variable-coefficient oscillatory integrals (Bourgain) would refute a variable-coefficient version. The method uses exact holomorphic extension of constant-coefficient waves and the polynomial flux map. Both fail for variable coefficients, so there is no conflict.
- **Logarithmic powers.** In $\mathbb R^2$ the Kakeya maximal power is $4/p_c=1\ge\frac12$ (Córdoba's sharp value), and the annular local smoothing power is $2/4=\frac12\ge\frac14$ (the sharp value). Both are consistent.
