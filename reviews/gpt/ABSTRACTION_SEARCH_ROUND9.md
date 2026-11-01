# Round 9: abstraction search, including the routes that failed

1 November 2026. The user requested a change of working method: search broadly, reduce the mechanisms to simpler objects, test speculative routes, and develop only a lead that survives. This is a research log, not a claim that the search was exhaustive.

**The lead developed:** [Kakeya 23](kakeya/23_GLOBAL_WAVE_CONE_RIGIDITY.md). A globally curved pencil at critical state dimension produces a sphere bundle whose total space is a sphere. The homotopy sequence forces exactly two spectral levels. Trace and variance then determine the entire dispersion, and its matrix algebra is a classical spin factor. The global hypotheses are essential.

## 1. Route map

| Starting abstraction | Probe | Outcome |
|---|---|---|
| Both proofs transport information through a simpler representation | Try to make their capacity barriers one universal information theorem | Too broad: already-local systems have no transport cost. No common worst-case class or useful new theorem emerged. |
| Biased energy is a decision made after random node deletion | Replace an average of optimized ratios by a ratio of averaged determinants | False. Exact $Q=3,D=1$ values are $503/1820$ versus $3/7$. Adaptive optimization matters. |
| Separate diminishing returns might be global concavity | Apply joint Jensen to all thinning probabilities at once | False, with an exact three-coordinate midpoint witness (§2). Separate concavity from ETH 21 survives. |
| Rational coordinates might encode much more transport per degree | Replace normalized vertex polynomials by degree-one rational functions | They can evade the degree-one polynomial obstruction, even without poles on the transport interval. But the common numerator reveals degree two (§3). No asymptotic algorithmic gain follows from this example. |
| Degree-one locality is an integrable height field | Turn edge equations into a potential with prescribed differences; inspect loop frustration | Useful reformulation and an explicit half-defect baseline (§4). The attempted sharper directional bound remains unproved. |
| Integrate amplitudes before squaring, as in interference | Replace divergent positive fibre mass by a compact contour period | The uncharged central Weyl period is zero: it erases the signal. A charged period survives, but no conserved quotient tensor has been constructed (§5). |
| Normalize redundant states by a preferred probability or “temperature” | Multiply by a Gaussian gauge cutoff | Positivity survives but the existing identity $\operatorname{div}(wB)=B\nabla w$ leaves a source. This is the previously known obstruction, not a new impossibility theorem. |
| Treat eigenvectors as states and directions as measurement settings | Ask whether every unit state acquires a unique setting when all local charts are globally regular | Productive: curvature plus critical dimension gives a global sphere bundle. Compactness is the missing ingredient. |
| A sphere divided into spheres should constrain dimensions | Use the homotopy sequence, then the opposite-direction spectrum | Productive: $n=m+1,\ N=2m$. No classification of arbitrary sphere fibrations is needed. |
| Two spectral outcomes should eliminate higher algebraic complexity | Recover mean and squared separation from traces | Productive: the pencil is a drifted anisotropic Clifford family; its Jordan algebra is a spin factor. This algebra is classical. The forcing argument is the addition to this repo. |

## 2. A precise failed shortcut: joint concavity

Here $\theta_i=4p_i(1-p_i)$ are the **thinning probabilities**, not the original biases. For the energy in ETH 21 with $Q=3,D=1,t_i=i,w_i=1$, take
$$
\theta^A=(3/4,1/4,3/4),\qquad
\theta^B=(3/4,3/4,1/4).
$$
Exact scalar-sector solves give
$$
E_A=\frac{26713}{80080},\qquad
E_B=\frac{8301}{34720},\qquad
E_{(A+B)/2}=\frac{247159}{865592}.
$$
Their joint Jensen gap is
$$
E_{(A+B)/2}-\frac{E_A+E_B}{2}
=-\frac{9374301}{11806674880}<0.
$$
Thus separate concavity cannot be promoted to joint concavity. No statement in ETH 21 needs changing: it explicitly asserted only separate concavity.

## 3. Rational coordinates can hide degree

At nodes $1,2$, assign according to vertex parity
$$
n_{\rm even}(T)=\frac{1+T/3}{1-T/3},\qquad
n_{\rm odd}(T)=1-3T.
$$
Both equal one at zero. Their values are $2,-2$ at the first node and $5,-5$ at the second, so every cube edge is local. The denominator is at least $1/3$ on $[0,2]$: poles on the transport interval are not necessary for this escape.

However,
$$
n_{\rm even}+n_{\rm odd}
=\frac{2-3T+T^2}{1-T/3}.
$$
Clearing the shared denominator gives an ordinary degree-two polynomial family, exactly enough to vanish at both nodes. The test refutes an unchanged-degree extension of the polynomial obstruction, not the effective-degree budget after clearing denominators. More general vertex-dependent denominators remain an open route requiring recovery and bit-complexity analysis.

## 4. Locality as a frustrated potential

For degree one, write $n_x(T)=1+a_xT$, and let $h_x=\operatorname{par}(x)a_x$. A good edge in direction $i$, oriented from an even vertex $x$ to an odd vertex $y$, requires
$$
h_x-h_y=-2/t_i.
$$
Hence the prescribed edge differences must sum to zero on every retained cycle. This recasts the search as deletion of inconsistent loops.

The elementary potential $h_x=2\sum_i x_i/t_i$ gives
$$
n_x(T)=1+2T\,\operatorname{par}(x)\sum_i x_i/t_i,\qquad
\nu_i(x)=1-\operatorname{par}(x)(1-2x_i).
$$
For $Q\ge2$, exactly half the edges in each direction fail. This is a rediscovered baseline, not a new sharpness theorem. A numerical mixed-integer search rejected a proposed improvement for $Q=4$; the $Q=5$ search reached its 20-second limit. Neither solver outcome is used as a mathematical proof.

## 5. Squaring after taking a period

At the central Weyl fibre, the relevant amplitude is
$$
\frac{\lambda+\lambda^{-1}}2.
$$
Its uncharged period against $d\lambda/\lambda$ vanishes: the Laurent constant coefficient is zero. Weighting by $\lambda^{-1}$ instead extracts $1/2$. The latter is a charged object; making it global and proving conservation, positivity and the required row estimates remains unresolved. The computation only rejects the naive uncharged period, not all oscillatory reductions.

## 6. Why the topology lead received the effort

The decisive shift was to stop asking how to integrate an infinite fibre and ask what a globally regular polarization family can be. Critical dimension makes evaluation from the unit eigenbundle a local diffeomorphism. Global compactness makes it surjective; convexity makes it injective. This identifies the total space with a sphere.

The next obstruction is tiny:
$$
0=\pi_m(S^k)\longrightarrow
\pi_{m-1}(S^{m-1})=\mathbb Z\longrightarrow
\pi_{m-1}(S^{k+m-1})=0
$$
cannot be exact when $k>m$. Opposite directions give $k\ge m$, so $k=m$. Two eigenspaces now exhaust the state space. Everything after that is linear algebra.

This is developed with the topology and Clifford ingredients explicitly sourced in [Kakeya 23](kakeya/23_GLOBAL_WAVE_CONE_RIGIDITY.md). The non-Clifford local example, the arrow pencil's south pole, and repeated Clifford representations test three distinct limits of the hypotheses.

## 7. Reproduction and open ends

~~~text
python reviews/gpt/experiments/round9_routes.py
python reviews/gpt/kakeya/scripts/verify_global_rigidity.py
~~~

The [route data](experiments/results/round9_routes.json) preserve exact counterexamples and constructions. The experiment script is separate from the verification runner and proves no asymptotic directional bound. The global theorem has a separate prose proof and [algebraic checks](kakeya/results/global_rigidity.json); a finite script does not verify its topology.

The productive next questions are local ones: how much of the spin-factor closure survives on a cap, whether a charged quotient can carry a conserved positive tensor, and whether loop frustration yields a sharper directional bound. None was settled in this turn.
