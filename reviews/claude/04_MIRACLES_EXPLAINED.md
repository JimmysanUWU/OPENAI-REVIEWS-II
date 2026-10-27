# 4. Does the paper rely on miracles? What makes it work

## 4.1 Summary

**No step relies on an unexplained coincidence.** Every identity has a proof, and every proof uses only:
- linear algebra over $\mathbb Q$;
- Möbius inversion and Hermite interpolation;
- the binomial theorem and Hoeffding's inequality;
- standard derandomization;
- in Appendix D, 2-adic valuations.

Four places look like miracles at first reading. The table rates how well the paper explains each.

| Apparent miracle | Paper's explanation | Rating |
|---|---|---|
| A degree-$R$ deformation is "local modulo measurements" at $Q\gg R$ nodes | §2.7 Steps 2–3: eigenbasis of the swaps, binomial concentration | Good. The mechanism is shown, but not *why* the complete swap is the right centre, nor whether $R\approx\sqrt{Q\log}$ is optimal. Both are supplied below. |
| Exponentially many measurement *values* cost only $2^{O(n/r)}$ | §2.6: $\tau\log(2D+1)=O(1/r)$ | Adequate. The point is that there are few rows ($\tau n$) and only logarithmically many bits per row. |
| Taylor jets of a $2^n$-term sum from few substitutions | §2.5, Lemma 3.3 | Good. |
| Lifting $\mathbb F_2$ systems keeps linear degree | Lemma 5.2 and its $\ell=2$ example | Good. |
| Appendix D's 2-adic taper and modes | One overview paragraph; dense proofs | Weak. The reasons for tapering by $v_2(t)$ are left implicit. A reading is given in §6 below. |

## 4.2 The central mechanism, explained by characters

Throughout, $p=\frac{r-1}{2r}$ and $W\sim\mathrm{Bin}(Q,p)$.

The auxiliary array is $\mathcal B=([r]\times[r])^Q$. The group $G=(\mathbb Z/2)^Q$ acts on it by swapping the pair at chosen positions; $\Pi_s$ is the swap at position $s$, and $\Pi_S=\prod_{s\in S}\Pi_s$.

**Decomposition.** For one position, $\mathbb Q^{r\times r}$ splits into symmetric matrices (dimension $\binom{r+1}2$, swap eigenvalue $+1$) and antisymmetric matrices (dimension $\binom r2$, eigenvalue $-1$). Hence
$$
\mathbb Q^{\mathcal B}=\bigoplus_{\chi\in\{0,1\}^Q}V_\chi,\qquad \Pi_s|_{V_\chi}=(-1)^{\chi_s},\qquad \dim V_\chi=\tbinom{r+1}2^{Q-\lvert\chi\rvert}\tbinom r2^{\lvert\chi\rvert}=D\cdot\Pr[\text{a fixed }p\text{-biased string equals }\chi].
$$

**Locality.** At node $t$, coordinate $y$ lies in group $a_y(t)$, and the row for $u$ is assigned to group $d_u(t)$. The permutation $\Pi_{S^c}$ sends $u$ to $y$ with $a_y(s)=d_u(s)$ for every $s\notin S$.

So $\Pi_{S^c}$ is *local at every node outside $S$*. The complete swap $V_0=\Pi_{[Q]}$ is the only group element that is local at every node. The identity $\Pi_\emptyset$ is local at none, except where $a_u(t)=d_u(t)$.

**Deformation.** We need a matrix polynomial $P(T)$ with $P(0)=I$, such that $P(t)$ is local at $t$ modulo few measurements. Take a symmetric group-algebra element supported on $\{\Pi_{S^c}:\lvert S\rvert<R\}$:
$$
M_0=\sum_{k<R}\gamma_k\sum_{\lvert S\rvert=k}\Pi_{S^c},\qquad
P(T)=I-M_0+\sum_{k<R}\gamma_k\sum_{\lvert S\rvert=k}\Big(\prod_{s\in S}(1-T/s)\Big)\Pi_{S^c}.
$$
The three requirements are then checked as follows.
- $P(0)=I$.
- $\deg_TP<R$.
- At $T=t$, every surviving term has $t\notin S$ and is therefore local at $t$. The only nonlocal part left is the constant $I-M_0$.

**Rank.** On $V_\chi$ the element $M_0$ acts by the scalar $(-1)^{\lvert\chi\rvert}f(\lvert\chi\rvert)$, where $f=\sum_k\gamma_kK_k$ and $K_k$ is the Krawtchouk polynomial. As $\gamma$ ranges over $\mathbb Q^R$, $f$ ranges over all polynomials of degree $<R$. Hence
$$
\operatorname{rank}(I-M_0)=D\cdot\Pr\big[(-1)^Wf(W)\neq1\big].
$$
The measurements are the projections onto the isotypic pieces where $M_0\neq1$.

**The only analytic input.** Choose $f$ to interpolate $(-1)^j$ on the $R$ central integers. Then the rank fraction is the binomial tail outside a window of width $R$, namely $e^{-\Omega(R^2/Q)}$.

The paper's $M_0=V_0f(W_0)$ is exactly this element: $f(W_0)$ expands into $e_k(\Pi)$'s. The paper uses $(u,y)$-specific multipliers $\prod_{s\in E_{u,y}}$, whereas the formula above uses $S$-specific ones. Both are valid; the second is the closed form of [05, Theorem B](05_DISTILLATION.md), verified exactly in [`results/distilled_swap_formula.json`](results/distilled_swap_formula.json).

So the "miracle" is a single classical fact. **Parity is approximable in agreement, on the bulk of a binomial distribution, by polynomials of degree $\sqrt{Q\log(1/\varepsilon)}$.** This is the fact behind the non-rigidity of Walsh–Hadamard matrices. The swap action makes the $G$-isotypic pieces carry binomial weights, and the complete swap supplies the sign $(-1)^{\lvert\chi\rvert}$.

## 4.3 The three costs

Per variable, $\log_2(\text{time})/n$ is at most
$$
\underbrace{\tfrac{d}{n}\tfrac RQ\log_2\tfrac{eQ}R}_{\text{Taylor jets}}
+\underbrace{2\tau\log_2(2\lVert C\rVert_{\ell_1}+1)}_{\text{measurement values}}
+\underbrace{2/r}_{\text{one group}},
\qquad \varepsilon=\Pr[W\notin\mathcal I]\le\frac{\tau}{rb},\quad R\gtrsim\sqrt{8Q\log\tfrac{rb}{\tau}},\quad r^{2Q}=D\le n^{1/128}.
$$
For bounded-occurrence 3-CNF, $d=3n$. The paper's choice $r\approx\sqrt{\log n}/(\log\log n)^2$ and $Q\approx\log n/\log\log n$ makes the first and third terms equal, at $(\log\log n)^2/\sqrt{\log n}$. With the paper's constants, the Taylor term dominates at the break-even point (0.990 of the exponent; [`results/parameters.json`](results/parameters.json)).

## 4.4 Proposition (the interval choice is optimal among symmetric constructions)

*Let $f$ have degree $<R$. Every window of $R+1$ consecutive integers in $\{0,\dots,Q\}$ contains some $j$ with $(-1)^jf(j)\neq1$. Consequently every $\mathrm{Sym}(Q)$-invariant $M_0$ as in §4.2 has*
$$
\frac{\operatorname{rank}(I-M_0)}D\ \ge\ \max_A\ \min_{A\le j\le A+R}\Pr[W=j].
$$

*Proof.* Suppose $f(j)=(-1)^j$ for $j=A,\dots,A+R$. The $R$-th forward difference of $f$ then vanishes, because $\deg f<R$. But the $R$-th difference of $(-1)^j$ at $A$ is $(-1)^{A+R}2^R\neq0$, a contradiction. The rank formula of §4.2 finishes the proof. ∎

For $W\sim\mathrm{Bin}(Q,p)$, the lower bound is $\gtrsim Q^{-1/2}e^{-(R/2)^2/(2Qp(1-p))}$. The interval construction attains $\lesssim e^{-(R/2)^2/(2Qp(1-p))}$, so it is optimal up to a $\mathrm{poly}(Q)$ factor. In the tested cases the measured ratio of upper to lower bound is 2.1–8.6 ([`results/distilled_swap_formula.json`](results/distilled_swap_formula.json)).

Hence $R=\Theta(\sqrt{Q\log(1/\varepsilon)})$ is forced for symmetric constructions. Non-symmetric group-algebra constructions reduce to the *agreement degree of parity under the $p$-biased measure*. I expect the same $\Theta(\sqrt{Q\log(1/\varepsilon)})$ from known probabilistic-degree lower bounds for parity over $\mathbb R$, but I have not checked a reference.

## 4.5 Why the method stops at $\sqrt{\log n}$ (a heuristic barrier)

Assume the architecture of §3–§4:
- a universal array of size $D\ge r^{2Q}$ with $D\le n\cdot\mathrm{polylog}\,n$, since at least $r\varepsilon D$ measurement rows are needed;
- rank fraction $\varepsilon\lesssim1/r$;
- $R\gtrsim\sqrt{Q\log(1/\varepsilon)}$, from §4.4.

Then $Q\lesssim\log n/(2\log r)$, and
$$
\frac{\log_2\text{time}}n\ \gtrsim\ \max\Big(\frac1r,\ \frac{\log r\cdot\log\log n}{\sqrt{\log n}}\Big)\ \ge\ \Omega\Big(\frac{(\log\log n)^2}{\sqrt{\log n}}\Big)\quad\text{up to constants.}
$$

The paper's parameters are therefore optimal for its method up to constants. Appendix D removes the $\log(Q/R)$ factor of subset counting by counting 2-adic "modes" instead. It still meets $\sqrt{\log}$, because the number of usable 2-adic levels is $u\lesssim\log s$.

Going past $2^{n/\sqrt{\log n}}$ needs a genuinely different source of decoupling. This is a heuristic statement about the method, not a theorem about #SAT.

## 4.6 A reading of Appendix D's taper

Throughout, $\alpha=2^{a+1}$, $D_t=L/2^{v_2(t)}$ and $u$ is the number of retained 2-adic levels.

Nodes are $t\in S=\{-L,\dots,L-1\}$, a complete residue system mod $2^{a+1}$. Groups at node $t$ are read off bit $i=v_2(t)$ of the slot numbers. Adding $t$ flips bit $i$, so $W(t)$ (the shift by $t$) maps every port to a port in the same group. This is the 2-adic analogue of the complete swap.

Near $t$, the crossing terms of $W(t+\alpha\zeta)$ have coefficients divisible by $2^{b+1}$; the "mode" of weight $b+1$ appears in the geometric expansion (D.12). Two quantities depend on the node:
- **Precision demand.** Working mod $2^p$, only mode monomials of total weight $<p_t$ survive. Recovering $F(0)$ from the jets loses $\approx6hD_t$ bits, by (D.10).
- **Mode supply.** The number of modes of each weight is $\approx M/D_t$.

The count $\exp(O(\sqrt{Mp_t/D_t}))=\exp(O(M/\sqrt u))$ is therefore independent of the node. Tapering the multiplicity $q_t=h(u-v_2(t))^+$ is what keeps $p_t\asymp MD_t/u$ while guaranteeing $\sum q_t>\deg F$.

Two features do the work:
- Odd nodes ($D_t=L$) carry most of the interpolation.
- 2-adically deep nodes ($v_2(t)\ge u$) are dropped: their groups would be too coarse.

My exact checks confirm the local lemmas (Lemma D.4, (D.10), Lemma D.5, (D.12)). A full recount of the recursion remains open.
