# Round 10: from fitted polynomials to minimal contradictions

1 November 2026. The user again requested abstraction and route searching before development. The lead that survived is [ETH 22](eth/22_CIRCUITS_AND_ERASURE.md): an affine system's optimized loss splits into its smallest contradictions after splitting weights without cancellation. This yields a sharp deletion/refitting theorem that does not require polynomials, Vandermonde matrices, squared loss, or full independence.

The elementary-vector ingredient is classical, as credited in the note. The contribution here is the explicit erasure theorem, its sharp independence threshold, and the separate interfaces with scalar energy and exact locality. Worldwide novelty is unestablished.

## 1. Routes and stopping decisions

| Abstraction | Test or obstruction | Decision |
|---|---|---|
| Degree might really mean number of fitted parameters | Replace powers of nodes by arbitrary regression features | Survived. Cauchy–Binet gives the same bound with no polynomial structure. |
| A contradiction should be the indivisible unit of fitting | Enumerate minimal inconsistent row sets and decompose the optimal dual residual | Stronger lead. Sign-compatible circuit splitting gives an exact decomposition of the full loss. |
| The square might be incidental | Repeat the splitting with the gradient of a convex loss | Survived. Every power $s>1$, and the stated differentiable convex row losses, work. Exact fourth-power fixtures supplement the general prose proof. |
| Independence might only be needed to keep a contradiction alive | Match joint retention probabilities on circuits, rather than on all masks | Survived. $c$-wise independence suffices when contradictions have size at most $c$. |
| Perhaps one fewer independence order suffices | Use a parity distribution that never retains the full contradiction | False, sharply: $(c-1)$-wise independent fair masks can make every sampled subsystem consistent. |
| Ambient rank might determine the loss exponent | Add unrelated contradictory pairs | Too coarse. Rank $d$ can coexist with circuit size two and exact factor $\theta$, rather than $\theta^d$. |
| The circuit hypergraph might determine the full energy law | Compare two rank-one systems with the same contradictory pairs | False: the normalized expected losses are $503/780$ and $5/8$. Geometry still matters. |
| A contradiction hypergraph alone might imply the same bound for arbitrary models | Minimize over six nonlinear residual states | False: size-two contradictions coexist with expected loss $5/16<\frac12(17/16)$. Affine structure is essential. |
| Frustration might be completely detected by square faces | Retain one inconsistent longer cube cycle and break all squares | False after deletion. The $Q=3$ circuit family has lengths 4, 6 and 8. |
| Repeated eigenvalues might always come from a common commuting gauge | Solve the real commutant exactly for the four critical Hopf pencils | False. Commutant dimensions are $1,2,4,1$ for multiplicities $1,2,4,8$. Octonionic multiplicity has no such eight-dimensional commuting symmetry. |
| Algebraic phase complexity might equal the number of necessary states | Factor the arrow pencil and span its top eigenvectors | False as a fixed-subspace compression rule. A quadratic phase factor uses top eigenvectors spanning three states. A weaker algebraic degree budget remains a candidate (§5). |

## 2. The shortest route through the successful idea

At an optimal fit, the weighted residual gradient $\lambda$ obeys
$$
A^{\mathsf T}\lambda=0.
$$
This is a balanced dependence among the equations. It can be split into support-minimal dependences, all with the same coordinate signs. Each component is a genuine contradiction, because its pairing with the original right side is positive.

Split the original row weights in the same proportions. The old fit remains optimal for every component, so **the full loss is exactly the sum of the component losses**. After random deletion, fitting the components independently can only reduce the loss. A component contributes only when its whole circuit survives.

For common retention probability $\theta$, a size-$j$ component survives with probability $\theta^j$, and inverse weighting divides its retained loss by $\theta$. Thus its factor is $\theta^{j-1}$. Summing gives
$$
\boxed{\mathbb E[\text{loss after deletion and refitting}]
\ge\theta^{c-1}[\text{original loss}],}
$$
provided the circuits have size at most $c$ and the required joint survival probabilities hold. Evaluation at the original fit gives the upper bound by the original loss.

This explains the old degree exponent: scalar degree-$D$ interpolation has contradictions of size exactly $D+1$. The exponent measures the size of an irreducible inconsistency, less the one power refunded by inverse weighting.

## 3. What survives, and what does not

The erasure theorem survives arbitrary feature matrices, nonnegative split weights, unequal retention probabilities, the stated convex losses, and limited independence. For squared loss, the determinant formula is a useful geometric refinement; fixed ridge penalties are added as equations that are never deleted.

The support problem uses the same minimal contradictions but a different optimization: delete enough equations to hit every contradiction. This exactly recasts the locality problem, including directional costs. It does not turn an energy bound into a support bound or produce a better asymptotic directional estimate.

Two boundaries were worth making explicit:

- Matching low-order retention moments does not determine the exact optimized expectation: the rank-one example has values $503/1820$ and $22/65$ under independent and pairwise independent masks.
- The loss theorem is not a theorem about an arbitrary family of states. The nonlinear six-state control has the right contradiction sizes and still violates the proposed bound. The cancellation-free affine dependence is doing real work.

At degree one, the circuit language is literal loop frustration. Exact enumeration of all deletion sets gives one inconsistent circuit at $Q=2$ and 24 at $Q=3$, with minimum largest directional counts one and two. Those finite values do not settle the sharp-direction target.

## 4. Rejected symmetry and compression shortcuts

For the critical Hopf pencils, every generator is a signed-permutation involution. The equations $XC_j=C_jX$ identify matrix entries up to signs. A connected component of those identifications contributes one free scalar unless a negative sign loop forces it to zero. This gives the exact commutant counts:

| Top multiplicity $m$ | Real state dimension $N$ | Dimension of real commutant |
|---:|---:|---:|
| 1 | 2 | 1 |
| 2 | 4 | 2 |
| 4 | 8 | 4 |
| 8 | 16 | 1 |

Hence a common commuting linear gauge cannot be the universal explanation of multiplicity or of the complex fibre's extra dimensions. This does not exclude nonlinear symmetries or other quotient constructions. It agrees with Kakeya 22's decision to compute the coarea directly without assuming a gauge group.

The real arrow pencil is
$$
C(x,y,z)=
\begin{pmatrix}z&x&y\\x&-z&0\\y&0&-z\end{pmatrix},\qquad
\det(tI-C)=(t+z)(t^2-x^2-y^2-z^2).
$$
Its curved phase has degree two. Yet top eigenvectors at the positive coordinate directions include
$$
(1,1,0),\quad(1,0,1),\quad(1,0,0),
$$
which span all three states. One cannot simply discard the linear spectral factor by restricting the original pencil to a fixed two-state subspace. At the south pole the top multiplicity changes, as already recorded in Kakeya 23.

## 5. A candidate retained for later checking

There may be a useful local **algebraic degree budget** behind the preceding failure. An irreducible characteristic factor of degree $h$ in the eigenvalue, occurring with generic multiplicity $m$, consumes at least $mh$ of the total state degree $N$. A genuinely curved phase cannot have degree one. At the limiting size $N=2m$, this suggests that a quadratic factor exhausts the spectrum; trace and variance would then force a wave cone locally, without the global sphere-bundle argument.

The relevant passage from a local constant-multiplicity branch to its homogeneous irreducible factor, and the treatment of other spectral factors, deserve a separate proof under the exact hypotheses. No new roadmap theorem is based on this candidate. The arrow example specifically blocks the stronger fixed-subspace compression claim, not the weaker factor-degree budget.

Another possible next step is a direction-sensitive circuit packing certificate for (22.16). Enumerating all circuits is not the goal; the useful outcome would be a short certificate that improves a proved directional bound. **T.ETH.7**, the local current/quotient problem **T.KAK.7**, and the broad cross-paper target **T.X.1** remain open.

## 6. Reproduction

~~~text
python reviews/gpt/eth/scripts/verify_circuit_erasure.py
python reviews/gpt/experiments/round10_routes.py
~~~

The first script is a deterministic verifier with [exact rational results](eth/results/circuit_erasure.json), including the false controls. The second retains the [commutant, arrow and geometry probes](experiments/results/round10_routes.json) separately from the verification runner. These are finite checks; the general assertions have their own proofs in ETH 22. No Lean build or improved headline exponent was attempted.
