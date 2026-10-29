# 11. Gapless corollaries: Max-k-SAT, colourings and matchings, circuit lower bounds, local maps

Claude, 28 October 2026. This file completes the ETH paper: it proves the corollaries of §1 and App. C from the counting theorems of [08](08_GAPLESS_PROOF.md)–[10](10_GAPLESS_TENSOR_CONTRACTION.md). The standard inputs are cited by name:
- the Tseitin encoding of circuits as 3-CNFs;
- Hoeffding's inequality;
- Žák's nondeterministic time hierarchy;
- Williams's witness-circuit argument;
- the Jahanjou–Miles–Viola local reduction;
- Jerrum–Valiant–Vazirani sampling.

**Verdict:** all five corollaries are correct. One dependency the paper does not flag: Cor. 1.4's bound $n\sqrt{\log n}$ needs App. D's fixed-density exponent. From the formally verified Prop. 6.4 alone, the same argument gives $c\,n\sqrt{\log n}/(\log\log n)^2$ (Remark 3.2).

**Which counting bound each corollary uses.**

| Result | Needs | Formal? |
|---|---|---|
| Cor. 1.2 (Max-$k$-SAT) | the decision version of Prop. 6.4 | yes |
| Cor. 1.3 (independent sets, colourings, matchings) | #$k$-SAT in $2^{o(n)}$ (Part I already suffices) | yes |
| Cor. C.1 ($E^{NP}\not\subseteq\mathrm{SIZE}(O(n))$) | 3-SAT in $2^{o(n)}$ | yes, plus Williams (cited) |
| Cor. 1.4 ($C_L(n)>cn\sqrt{\log n}$ i.o.) | 3-SAT at fixed density in $2^{O(n/\sqrt{\log n})}$ | no (App. D, proved in [10](10_GAPLESS_TENSOR_CONTRACTION.md)) |
| Cor. C.2 (local maps) | #$k$-SAT in $2^{O(\rho(n))}$ | yes |

Throughout, $\rho(n)=n(\log\log n)^2/\sqrt{\log n}$, and **Tseitin** means the following. A fan-in-two circuit with $n$ inputs and $S$ gates becomes a 3-CNF with one variable per input and per gate, at most four clauses per gate forcing the gate's value, and one unit clause forcing the output to $1$. It has $n+S$ variables and $O(S)$ clauses. Each input assignment accepted by the circuit extends to exactly one satisfying assignment, so the encoding preserves both satisfiability and counts.

---

## 1. Max-$k$-SAT (paper Cor. 1.2)

**Corollary 1.1.** For fixed $k\ge2$ and $\varepsilon>0$, a randomized algorithm finds an assignment satisfying at least $\mathrm{OPT}(F)-\varepsilon m$ clauses with probability $1-2^{-\Omega(n)}$. It uses $(n+m)^{O_{k,\varepsilon}(1)}2^{O_{k,\varepsilon}(\rho(n))}$ time and polynomial space.

*Proof.*
1. **Sampling.** Draw $s=\lceil8(n+1)/\varepsilon^2\rceil$ clause indices independently and uniformly. Let $\hat f(x)$ be the fraction of sampled clauses that $x$ satisfies, and $f(x)$ the true fraction. By Hoeffding, $\Pr[|\hat f(x)-f(x)|>\varepsilon/2]\le2e^{-s\varepsilon^2/2}\le2e^{-4(n+1)}$. A union bound over the $2^n$ assignments leaves failure probability $2^{-\Omega(n)}$.
2. **Exact maximization of $\hat f$.** For a threshold $b$, build a fan-in-two circuit that evaluates the $s$ sampled clauses ($O(ks)$ gates), adds the $s$ bits (an adder tree, $O(s)$ gates) and compares with $b$ ($O(\log s)$ gates). Its size is $O_{k,\varepsilon}(n)$.
   - Tseitin gives a 3-CNF with $N=O_{k,\varepsilon}(n)$ variables, satisfiable iff some $x$ has $s\hat f(x)\ge b$.
   - Binary search on $b$ ($O(\log s)$ decisions) finds $b^*=s\max\hat f$.
   - Self-reduction ($N$ more decisions) recovers an optimizer $\hat x$.
   - Each decision costs $2^{O(\rho(N))}=2^{O_{k,\varepsilon}(\rho(n))}$ by Prop. 6.4 ([09](09_GAPLESS_SHARP_EXPONENT.md)).
3. **Quality.** With $x^*$ optimal, $f(\hat x)\ge\hat f(\hat x)-\frac\varepsilon2\ge\hat f(x^*)-\frac\varepsilon2\ge f(x^*)-\varepsilon$.
4. **Sampling cost.** Rejection sampling of indices in $[m]$, truncated after $O(n+\log s)$ trials each, gives worst-case polynomial time and adds $2^{-\Omega(n)}$ to the failure probability. ∎

---

## 2. Independent sets, colourings, perfect matchings (paper Cor. 1.3)

**Corollary 2.1.** On a simple graph with $v$ vertices and $e$ edges, three counts are computable deterministically:
- the number of independent sets, in $2^{o(v)}$ time and $\mathrm{poly}(v,e)$ space;
- for fixed $q\ge2$, the number of proper $q$-colourings, in $2^{o(v)}$ time and $\mathrm{poly}_q(v,e)$ space;
- the number of perfect matchings, in $2^{o(v+e)}$ time and $\mathrm{poly}(v,e)$ space.

*Proof.* Each count is the number of satisfying assignments of a CNF built below; apply #$k$-SAT ([08](08_GAPLESS_PROOF.md), or Prop. 6.4).
- **Independent sets.** The 2-CNF $\bigwedge_{\{u,w\}\in E}(\neg x_u\vee\neg x_w)$ on $v$ variables. Its satisfying assignments are the indicator vectors of independent sets.
- **$q$-colourings.** Variables $x_{u,c}$, with the following clauses:
  - $(x_{u,1}\vee\dots\vee x_{u,q})$ and $(\neg x_{u,c}\vee\neg x_{u,c'})$ for $c<c'$: exactly one colour per vertex;
  - $(\neg x_{u,c}\vee\neg x_{w,c})$ per edge and colour.
  
  This is a $q$-CNF on $qv$ variables with $O_q(v+e)$ clauses, in bijection with proper colourings. The time is $2^{o(qv)}=2^{o(v)}$.
- **Perfect matchings.** One variable per edge, and a circuit of size $O(v+e)$ checking that each vertex meets exactly one selected edge. "Exactly one of $d$" costs $O(d)$ gates via prefix ORs. Tseitin gives $O(v+e)$ variables; the gate values are determined by the edge bits, so counts are preserved. ∎

---

## 3. Circuit lower bounds for $E^{NP}$ (paper Cors. C.1, 1.4)

**Corollary 3.1 (C.1).** $E^{NP}\not\subseteq\mathrm{SIZE}(O(n))$.

*Proof.* Williams proved: if 3-SAT on $n$ variables is decidable in $2^{o(n)}$ time, then $E^{NP}$ has no linear-size circuits. The paper cites this as [Wil13, Thm 6.1]; I cite it, not re-proved. The hypothesis holds by Prop. 6.4, which is also formally verified. ∎

**Corollary 3.2 (1.4).** Let $C_L(n)$ be the least number of non-input gates of a fan-in-two circuit, with arbitrary gates, computing $L$ on $\{0,1\}^n$. Some $L\in E^{NP}$ and $c>0$ have $C_L(n)>c\,n\sqrt{\log n}$ for infinitely many $n$.

*Proof.* The skeleton is Williams's witness-circuit argument, and the parameter bookkeeping is checked here.
1. **A hard language.** By Žák's hierarchy, fix a unary $H\in\mathrm{NTIME}(2^n)\setminus\mathrm{NTIME}(2^{n/2})$ on multitape machines.
2. **Succinct formulas.** The Jahanjou–Miles–Viola local reduction (with locality parameter $2$, as in the paper) maps $1^n$ to a 3-CNF $\varphi_n$ with $2^q$ clauses, $q=n+O(\log n)$. It comes with a decoder that maps a clause index to its three variable labels and signs. Each decoder output bit depends on at most $w=O(1)$ index bits; recovering each support and truth table by evaluating on low-weight strings gives a decoder circuit with $O(n)$ gates. Moreover $\varphi_n$ is satisfiable iff $1^n\in H$.
3. **The witness language.** Pad the labels to $dn$ bits. Let $W(i)$ be bit $i$ of the lexicographically first satisfying assignment of $\varphi_n$ when $|i|=dn$ and $\varphi_n$ is satisfiable, and $0$ otherwise.
   - Writing out $\varphi_n$ takes $2^{O(n)}$ time, and $2^{dn}$ NP-oracle queries recover the first assignment by prefixes.
   - So $W\in E^{NP}$.
4. **The fast decision step.** By Proposition 8.1 of [10](10_GAPLESS_TENSOR_CONTRACTION.md), 3-SAT with $V$ variables and at most $5V$ clauses is decided in time $V^b2^{AV/\sqrt{\log V}}$. A multitape simulation of the polynomial-space solver costs only a polynomial factor, absorbed in $b$.
5. **Contradiction.** Fix $\eta<1/(12A)$, put $s_n=\lfloor\eta n\lfloor\sqrt{\log n}\rfloor\rfloor$, and suppose $C_W(dn)\le s_n$ for all large $n$. Decide $H$ on $1^n$ nondeterministically:
   - Guess a circuit of $\le s_n$ gates on $dn$ inputs; it has a description of $n^{1+o(1)}$ bits.
   - Build the test circuit: decoder $\to$ three copies of the guess $\to$ "is the indexed clause violated?". It has $q$ inputs and $\le3s_n+O(n)$ gates.
   - Tseitin, with padding, gives $U_n=3s_n+Bn$ variables and $\le5U_n$ clauses.
   - Accept iff it is unsatisfiable. If $1^n\in H$, the true witness circuit makes it unsatisfiable. If $1^n\notin H$, every guess violates some clause.
   
   The time is $2^{AU_n/\sqrt{\log U_n}}\mathrm{poly}=2^{(3A\eta+o(1))n}\le2^{n/4+o(n)}$, since $U_n=(3\eta+o(1))n\sqrt{\log n}$ and $\log U_n\sim\log n$. So $H\in\mathrm{NTIME}(2^{n/2})$, with finitely many lengths hardwired. This contradicts step 1.
6. **Conclusion.** Hence $C_W(dn)>s_n\sim\eta n\sqrt{\log n}$ infinitely often. At input length $N=dn$ this is $C_W(N)>c\,N\sqrt{\log N}$ for any $c<\eta/d$, adjusting the constant for the base of the logarithm. ∎

**Remark 3.2 (what the formally verified bound gives).** Run step 4 with Prop. 6.4's time $2^{AV(\log\log V)^2/\sqrt{\log V}}$ instead. The same computation forces $s_n\asymp n\sqrt{\log n}/(\log\log n)^2$. So Cor. 1.4 as stated depends on App. D, which is prose-only and is proved in [10](10_GAPLESS_TENSOR_CONTRACTION.md). The weaker bound $C_L(n)>c\,n\sqrt{\log n}/(\log\log n)^2$ i.o. follows from formally verified counting plus the cited circuit-complexity theorems.

---

## 4. Inverting local maps (paper Cor. C.2)

**Corollary 4.1.** Fix $d$. Let $G:\{0,1\}^n\to\{0,1\}^m$ have output bits $G(x)_i=g_i(x|_{S_i})$ with $|S_i|\le d$. For a target $b$, $|G^{-1}(b)|$ is computed, and a preimage found when one exists, in $(n+m)^{O_d(1)}2^{O_d(\rho(n))}$ time and polynomial space.

*Proof.*
1. **Encoding.** For each $i$ and each local assignment $\alpha\in\{0,1\}^{S_i}$ with $g_i(\alpha)\ne b_i$, add the clause that forbids $\alpha$. An empty $S_i$ is checked directly. The result is a CNF of width $\le d$ with $\le2^dm$ clauses, whose satisfying assignments are exactly $G^{-1}(b)$.
2. **Counting.** Apply #$k$-SAT with $k=\max(2,d)$.
3. **Finding a preimage.** Variable-by-variable self-reduction uses $n$ further counts. Conditional counts also give exactly uniform sampling of preimages (Jerrum–Valiant–Vazirani). ∎

**Status.** With [08](08_GAPLESS_PROOF.md)–[11](11_GAPLESS_COROLLARIES.md), every claim in the ETH paper has a complete proof here, modulo the named classical theorems. The new consequence of finishing the paper is [12](12_FP_COUNTING.md), which counts over every prime field.
