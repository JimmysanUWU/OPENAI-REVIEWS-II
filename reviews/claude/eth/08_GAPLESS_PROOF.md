# 8. Gapless proof: #k-SAT in subexponential time, hence ETH is false

Claude, 27 October 2026. This file is a complete, self-contained proof of the ETH-refuting chain of *The Exponential Time Hypothesis Is False*. Every step is proved here, apart from three textbook facts named where they are used: Bertrand's postulate, the Chinese remainder theorem for comaximal ideals, and polynomial-overhead simulation of random-access machines by Turing machines.

The companion script [`scripts/verify_gapless_chain.py`](scripts/verify_gapless_chain.py) checks every lemma below in exact arithmetic ([`results/gapless_chain.json`](results/gapless_chain.json), ~30 s):
- Lemmas 2.1–2.3 and Proposition 2.4 at every node of the branching trees of 32 random formulas;
- Lemmas 3.1–3.2 on 40 formulas;
- the whole counting algorithm of §4–§7 on 16 instances, including cases with several copies ($h=2,3$), against brute force, with a negative control;
- the parameter inequalities of Lemma 7.2.

**What changed relative to the paper.** The mathematics is the paper's minimal chain (§2 overview, §3–4, §2.8, App. A), with four simplifications.
1. **Arithmetic modulo $p^\mu$** with a prime $p>Q$ replaces rational arithmetic. This removes App. B's bit-length bookkeeping.
2. **Closed-form swap rows** (Theorem 5.4) replace the rows (4.10).
3. **Each measurement row meets at most two groups per node** (Lemma 5.1(4)). This removes a factor $r$ from the measurement count.
4. **New parameters** $(r,Q,R)=(\lfloor\sqrt\lambda\rfloor,\lfloor\lambda/(4\nu)\rfloor,2\lceil\sqrt\lambda\rceil)$, with every inequality proved for all $\lambda\ge2^{17}+1$. They give the explicit exponent $13\,n\nu^2/\sqrt\lambda$. The bound beats $2^n$ from $\log_2n\approx9.3\cdot10^6$; the paper's parameters need about $7\cdot10^{10}$ ([`results/parameters.json`](results/parameters.json)).

**Scope.** The Main Theorem below has exponent $n(\log\log n)^{k+1}/(\log n)^{1/(2k)}$. That refutes ETH and #ETH. The paper's sharper exponent $n(\log\log n)^2/\sqrt{\log n}$ for every $k$ (§5–6) is **not** re-proved here. It is the statement of OpenAI's Lean proof, which I built and audited ([06](06_LEAN_STATEMENT_AUDIT.md) §6.A).

---

## 0. Statement

Throughout, $\log=\log_2$ and $\ln$ is the natural logarithm.

> **Main Theorem.** Fix $k\ge2$. There is a deterministic algorithm that, given a $k$-CNF $F$ with $n$ variables and $m$ clauses, outputs the number $\#F$ of satisfying assignments using
> $$2^{O_k\left(n(\log\log n)^{k+1}/(\log n)^{1/(2k)}\right)}\cdot\mathrm{poly}(n+m)$$
> bit operations and $\mathrm{poly}(n+m)$ bits of space.

> **Corollary (ETH is false).** 3-SAT on $n$ variables and $m$ clauses can be decided in time $2^{o(n)}\mathrm{poly}(m)$. Hence there is no $\delta>0$ such that 3-SAT requires time $2^{\delta n}$, which is the Exponential Time Hypothesis. The counting version (#ETH) fails as well.

*Proof of the corollary.* Take $k=3$ and output whether $\#F>0$. The exponent $n(\log\log n)^4/(\log n)^{1/6}$ is $o(n)$. ∎

**The chain.**
```
k-CNF, n variables
  │ §2 disjoint sparsification: ≤ (Cn+1)²·2^{σn} leaves, each with ≤ B occurrences per variable
  ▼
k-CNF leaves (counts add up exactly)
  │ §3 OR gates + Tovey cycles: ≤ 6kBn variables, width ≤ 3, each variable in 1..3 clauses
  ▼
bounded-occurrence 3-CNF
  │ §4 deformation counting  +  §5 swap interpolation  +  §6 transfer and schedules
  ▼
§7 Theorem 7.1: #F' in poly·2^{13 N ν²/√λ}  (λ = ⌈log N⌉, ν = ⌈log λ⌉)
```

## 1. Conventions

**Formulas.**
- A *literal* is $x_v$ or $\neg x_v$. A *clause* is a finite set of literals, read as their disjunction. A *CNF* is a finite set of clauses.
- A clause is *tautological* if it contains a literal and its negation. Its *width* is its size.
- $\#F$ counts the assignments of **all** $n$ declared variables that satisfy every clause. $\mathrm{Sat}(F)\subseteq\{0,1\}^n$ is the set of these assignments.
- A variable *occurs* in a clause if the clause contains one of its literals. Its *occurrence number* is the number of such clauses.

**Normalization.** $\mathrm{norm}(F)$ deletes tautological clauses and then every clause that strictly contains another clause. It does not change $\mathrm{Sat}(F)$, and its result is an *antichain*: no clause contains another. A normalized CNF of width $\le k$ has at most $\sum_{c\le k}(2n)^c$ clauses.

**Clause polynomials.** For a clause $C$ with variable set $I_C$, put
$$\Phi_C(Y)=1-\prod_{\ell\in C}\bigl(1-\ell(Y)\bigr)\in\mathbb Z[Y_v:v\in I_C],\qquad x_v(Y)=Y_v,\quad(\neg x_v)(Y)=1-Y_v .$$
Then $\deg\Phi_C=|C|$, and $\Phi_C(x)=[x\text{ satisfies }C]$ for $x\in\{0,1\}^{I_C}$.

**Cost.** Time is counted in bit operations of a random-access machine and space in bits. A multitape Turing machine simulates such a machine with polynomial overhead in time and space. This changes only the constants in bounds of the form $2^{O(\cdot)}\mathrm{poly}(\cdot)$.

**Notation.**
- $[n]=\{1,\dots,n\}$.
- $\binom xi=x(x-1)\cdots(x-i+1)/i!$ for any integer or indeterminate $x$; it is an integer whenever $x$ is.
- $C_q(n)=\sum_{j=0}^{\min(q,n)}\binom nj2^j$.
- For $S\subseteq[n]$ and $\xi\in\mathbb Q^S$, $\widetilde\xi\in\mathbb Q^n$ is the extension by zero, and $x|_S$ is the restriction.

---

## 2. Disjoint sparsification

Fix $k\ge2$ and an integer $M\ge2$. Put $\theta_j=M^j$ for $0\le j<k$. For a CNF $F$, a width $c$ and a literal set $H$, let
$$\deg_c(H)=\#\{C\in F:\ |C|=c,\ H\subseteq C\}.$$

**Branching rule.** Let $F$ be normalized and without the empty clause.
- A pair $(c,H)$ is **eligible** if $H$ is a nonempty proper subset of some width-$c$ clause and $\deg_c(H)\ge\theta_{c-|H|}$.
- If no pair is eligible, $F$ is a **leaf**.
- Otherwise choose the eligible pair with the least $c$, then the largest $|H|$, then the first in a fixed order.

Let $H\cup P_1,\dots,H\cup P_z$ be the width-$c$ clauses containing $H$, with $P_a\cap H=\emptyset$. Put $j=c-|H|$, so that $|P_a|=j$ and $z\ge\theta_j$. Let $G$ be the remaining clauses. The two children are
$$F_A=\mathrm{norm}\bigl(G\cup\{H\}\bigr),\qquad F_B=\mathrm{norm}\Bigl(G\cup\{P_1,\dots,P_z\}\cup\bigl\{\{\neg\ell\}:\ell\in H\bigr\}\Bigr).$$
A clause of a child is **new** if it is not a clause of the parent; it is *inserted* by that branch.

**Lemma 2.1 (one branch).** Let $F$ be normalized and without the empty clause, and let the rule choose $(c,H)$ with $j=c-|H|$.
1. $\mathrm{Sat}(F)=\mathrm{Sat}(F_A)\sqcup\mathrm{Sat}(F_B)$.
2. The new clauses are as follows:
   - for $F_A$, exactly $H$;
   - for $F_B$, the $z\ge\theta_j$ distinct sets $P_a$, together with those units $\{\neg\ell\}$ ($\ell\in H$) that are not in $F$.
   
   Every new clause is nonempty and has width $<c$.
3. Every clause of $F$ that is not a clause of a child strictly contains a clause that is new in that child.
4. Unit clauses are never removed, and the empty clause never appears.

*Proof.*
1. $F=G\wedge\bigwedge_a(H\vee P_a)$.
   - If some literal of $H$ is true, every $H\vee P_a$ holds, and $F$ is equivalent to $G\wedge H$.
   - If every literal of $H$ is false, $H\vee P_a$ is equivalent to $P_a$.

   The two cases are exclusive and exhaustive, and normalization does not change satisfying sets.
2. No clause of $F$ is a subset of $H$ or of any $P_a$: it would be a strict subset of the clause $H\cup P_a$ of the antichain $F$. Hence:
   - $H\notin F$, $P_a\notin F$, and no old clause is strictly inside $H$ or $P_a$.
   - $H\cup P_a$ is not tautological, so $P_a$ contains no literal of a variable of $H$. Hence no unit $\{\neg\ell\}$ with $\ell\in H$ lies inside $P_a$.
   - The $P_a$ are distinct (they are determined by the distinct clauses $H\cup P_a$) and have equal width, so none contains another.
   - A unit can strictly contain only the empty clause, which is absent.

   So every listed clause survives normalization. Widths: $|H|<c$; $|P_a|=j<c$; units have width $1<c$, since $c\ge2$.
3. A selected clause $H\cup P_a$ strictly contains $H$ (as $j\ge1$) and $P_a$ (as $|H|\ge1$). Any other clause of $F$ is removed only if it strictly contains a clause of the child. That clause is not in $F$, because $F$ is an antichain, so it is new.
4. A unit is never selected, since selected clauses have width $c\ge2$. A unit could only be removed for strictly containing the empty clause. By item 2 the empty clause is never inserted. ∎

**Insertion accounting.** Fix a root-to-node path. A clause is **marked** when a branch on the path inserts it. $T_j$ is the number of insertions of width-$j$ clauses along the path; it counts events, so a reinserted clause counts twice.

**Lemma 2.2.** At every node, for all $1\le i<j<k$ and every set $S$ of $i$ literals, fewer than $2\theta_{j-i}$ marked width-$j$ clauses of the node's formula contain $S$.

*Proof.* Initially nothing is marked, and deletions only lower the count. So consider a branch with selected $(c,H)$ that inserts width-$j$ clauses. By Lemma 2.1(2), $c>j$.

*Before the branch*, fewer than $\theta_{j-i}$ width-$j$ clauses contain $S$. Otherwise $S$ would be a nonempty proper subset ($1\le i<j$) of a width-$j$ clause with $\deg_j(S)\ge\theta_{j-i}$. Then $(j,S)$ would be eligible with $j<c$, contradicting the least-width priority.

*New width-$j$ clauses containing $S$:*
- Child $A$ inserts only $H$, so at most one new clause contains $S$. The total is then $<\theta_{j-i}+1\le2\theta_{j-i}$.
- In child $B$ the units have width $1<j$. If $|P_a|\ne j$ nothing of width $j$ is inserted, so assume $|P_a|=j$ and some $P_a\supseteq S$.
  - Then $S\cap H=\emptyset$, and $S\subseteq P_a$ iff $H\cup S\subseteq H\cup P_a$.
  - So the number of $a$ with $S\subseteq P_a$ equals $\deg_c(H\cup S)$.
  - $H\cup S$ is a nonempty proper subset of the width-$c$ clause $H\cup P_a$, larger than $H$, because $|H\cup S|=|H|+i<|H|+j=c$.
  - By the largest-$|H|$ priority it is not eligible: $\deg_c(H\cup S)<\theta_{c-|H\cup S|}=\theta_{j-i}$.
  
  The total is again $<2\theta_{j-i}$. ∎

**Lemma 2.3.** Along every path, $T_j\le4(3M)^{j-1}n$ for $1\le j<k$.

*Proof.* Write $T_j=L_j+D_j$, where:
- $L_j$ counts marked width-$j$ clauses still present;
- $D_j$ counts width-$j$ insertions that were later deleted.

*Width 1.* Units are never deleted (Lemma 2.1(4)), and a clause is inserted only if absent. So $T_1=L_1\le2n$.

*Width $j\ge2$.*
- **Present clauses.** Apply Lemma 2.2 with $i=1$ to each of the $2n$ literals. Every width-$j$ clause contains at least one literal, so $L_j<2n\cdot2\theta_{j-1}=4n\theta_{j-1}$.
- **Deleted clauses.** When a marked width-$j$ clause is deleted, it strictly contains a clause $N$ inserted by the same branch (Lemma 2.1(3)). Charge the deletion to the insertion of $N$; its width $i$ satisfies $1\le i<j$. Lemma 2.2, applied just before that branch with $S=N$, shows that one insertion of width $i$ receives fewer than $2\theta_{j-i}$ charges. Hence $D_j\le2\sum_{i<j}\theta_{j-i}T_i$.

So $T_j\le4n\theta_{j-1}+2\sum_{i=1}^{j-1}\theta_{j-i}T_i$. Dividing by $nM^{j-1}$ and inducting gives
$$\frac{T_j}{nM^{j-1}}\le4+2\sum_{i=1}^{j-1}4\cdot3^{i-1}=4\cdot3^{j-1}.\qquad ∎$$

**Proposition 2.4 (the tree).** Put
$$C=4\sum_{j=0}^{k-2}(3M)^j,\qquad\Gamma=2(3^{k-1}-1),$$
and assume $M\ge\Gamma$. Run the rule from $\mathrm{norm}(F)$.
1. Every path has length $\le Cn$, and at most $\Gamma n/M$ of its steps go to a $B$-child.
2. The tree has at most $(Cn+1)^2\,2^{(\Gamma n/M)\log(4CM)}$ nodes.
3. Every leaf is a normalized CNF of width $\le k$ in which each variable occurs in at most $B:=2\sum_{j=0}^{k-1}M^j\le4M^{k-1}$ clauses.
4. The leaves' satisfying sets partition $\mathrm{Sat}(F)$, so $\#F=\sum_{\text{leaves }L}\#L$.
5. A node is processed in $\mathrm{poly}_k(n)$ bit operations, and a depth-first traversal uses $\mathrm{poly}_k(n)$ space.

*Proof.*
1. Every branch inserts a clause (Lemma 2.1(2)), so the length is at most $\sum_{j<k}T_j\le Cn$ by Lemma 2.3. A $B$-step with $|P_a|=j$ inserts at least $M^j$ width-$j$ clauses, so there are at most $T_j/M^j\le4\cdot3^{j-1}n/M$ such steps. Summing over $j$ gives $\Gamma n/M$.
2. The rule is deterministic, so a node is determined by its depth $L\le Cn$ and the set of $B$-positions on its path. Let $a_*=\Gamma n/M$, so $a_*\le n\le Cn/4$.
   - For $1\le a\le a_*$, $\binom{Cn}a\le(eCn/a)^a\le(4Cn/a)^a\le(4Cn/a_*)^{a_*}$. The first step uses $a!\ge(a/e)^a$; the last holds because $a\mapsto a\ln(4Cn/a)$ increases on $(0,4Cn/e]$.
   - Hence the number of nodes is at most $(Cn+1)\sum_{a\le a_*}\binom{Cn}a\le(Cn+1)^2(4CM/\Gamma)^{\Gamma n/M}$. This is at most the stated bound since $\Gamma\ge1$.
3. At a leaf no pair is eligible. Fix a literal $\ell$ and a width $2\le c\le k$: then $\deg_c(\{\ell\})<\theta_{c-1}$. As a set, $F$ contains the unit $\{\ell\}$ at most once. So $\ell$ lies in at most $\sum_{j<k}M^j$ clauses, and a variable has two literals. Also $\sum_{j<k}M^j\le2M^{k-1}$ since $M\ge2$.
4. Apply Lemma 2.1(1) at every branch.
5. A normalized formula has $\mathrm{poly}_k(n)$ clauses, each with fewer than $2^k$ subsets, so counting $\deg_c$, selecting and normalizing are polynomial. Depth-first search stores one path of $\le Cn$ formulas. ∎

**Lemma 2.5 (choice of $M$).** Let $0<\sigma\le1$, and write $C_M$ for the $C$ of Proposition 2.4. Let $M$ be the least power of two with $M\ge\max(\Gamma,4)$ and $\Gamma\lceil\log(4C_MM)\rceil\le\sigma M$. Then:
- $M\le2\max(\Gamma,4,M_*)$, where $M_*=4A\log(2A)$ and $A=2\Gamma(k+2)/\sigma$; hence $M=O_k(\sigma^{-1}\log(2/\sigma))$;
- the tree has at most $(Cn+1)^22^{\sigma n}$ nodes, and $B=O_k\big((\sigma^{-1}\log(2/\sigma))^{k-1}\big)$.

*Proof.*
1. **A simpler sufficient condition.** For $M\ge4$:
   - $C_M\le4\cdot\tfrac{12}{11}(3M)^{k-2}\le5(3M)^{k-2}$;
   - hence $4C_MM\le20\cdot3^{k-2}M^{k-1}\le(3M)^{k+1}$;
   - hence $\lceil\log(4C_MM)\rceil\le(k+2)\log(3M)\le2(k+2)\log M$, since $\log(3M)\ge1$ and $3M\le M^2$.
   
   So the defining condition holds whenever $\log M\le M/A$.
2. **When that holds.** Note $A\ge32$. At $M=M_*$: $\log M_*\le\log(8A^2)=3+2\log A<4+4\log A=M_*/A$. Moreover $M/A-\log M$ increases for $M\ge A/\ln2$, and $M_*\ge4A\ge A/\ln2$. So every $M\ge\max(\Gamma,4,M_*)$ qualifies, and the least qualifying power of two is at most twice that maximum.
3. **Consequences.** The node bound follows from Proposition 2.4(2), and $B\le4M^{k-1}$. ∎

---

## 3. Reduction to at most three occurrences

**Lemma 3.1 (OR gate).** Let $C=\{\ell_1,\ell_2,\ell_3,\dots,\ell_w\}$ be a non-tautological clause with $w\ge4$, and let $z$ be a new variable. Replace $C$ by the four clauses
$$\{\neg\ell_1,z\},\quad\{\neg\ell_2,z\},\quad\{\ell_1,\ell_2,\neg z\},\quad\{z,\ell_3,\dots,\ell_w\}.$$
This preserves the number of satisfying assignments of the formula.

*Proof.* The first three clauses say $\ell_1\Rightarrow z$, $\ell_2\Rightarrow z$ and $z\Rightarrow\ell_1\vee\ell_2$, that is, $z=\ell_1\vee\ell_2$. So each old assignment has exactly one extension, and under it the last clause is equivalent to $C$. ∎

**Lemma 3.2 (Tovey cycle).** Let variable $x$ occur in clauses $D_1,\dots,D_p$ with $p\ge2$. Replace its occurrence in $D_i$ by a new variable $x_i$ with the same sign, and add the clauses $\{\neg x_i,x_{i+1}\}$ for $i\in[p]$, with $x_{p+1}=x_1$. This preserves the count.

*Proof.* The cycle $x_1\Rightarrow x_2\Rightarrow\dots\Rightarrow x_p\Rightarrow x_1$ forces all $x_i$ to be equal. This gives a bijection of satisfying assignments. ∎

**Proposition 3.3.** Let $L$ be a CNF of width $\le k$ on $n$ variables, with no empty clause, in which each variable occurs in at most $B$ clauses. In polynomial time one obtains:
- a CNF $L'$ on $N$ variables, with every clause of width $\le3$ and every variable occurring in $1$ to $3$ clauses;
- an integer $z\ge0$ with $\#L=2^z\,\#L'$;
- the bound $N\le6kBn$; for $k\le3$, $N\le Bn$.

*Proof.*
1. **Free variables.** Let $z$ be the number of variables occurring in no clause. Drop them; this divides the count by $2^z$.
2. **OR gates** (Lemma 3.1), applied until every clause has width $\le3$.
   - Each gate adds one variable and $6$ literal occurrences.
   - A clause of width $w$ needs $w-3$ gates.
   - The original occurrences number $O\le Bn$, and the number of clauses is at most $O$.
   
   So the total number of occurrences becomes at most $O+6(k-3)O\le6kBn$.
3. **Tovey cycles** (Lemma 3.2), applied to every variable occurring at least twice. Each copy then occurs exactly $3$ times, and a variable occurring once is left alone. The resulting $N$ is the total number of occurrences.

For $k\le3$ there are no gates, so $N\le O\le Bn$. ∎

---

## 4. Deformation counting

Fix a CNF $F$ on $[n]$ with clauses $C_f$ ($f\in\mathcal F$), scopes $I_f$, and $d:=\sum_f|C_f|$. Let $p$ be a prime and $\mu\ge1$. Put $A=\mathbb Z/p^\mu$ and let $\mathcal O\subset\mathbb Q$ be the rationals whose denominators are prime to $p$; reduction $\pi:\mathcal O\to A$ is a ring homomorphism.

**Definition 4.1 (measured partition, locality).** A *measured partition* of $[n]$ consists of:
- a partition $[n]=S_1\sqcup\dots\sqcup S_r$;
- integer matrices $C_g\in\mathbb Z^{c_g\times S_g}$ whose rows have $\ell^1$-norm $\le\beta$.

A row $\rho\in\mathbb Q^{1\times n}$ is **$g$-local** if, for every $h\ne g$, $\rho|_{S_h}$ lies in the $\mathbb Q$-row space of $C_h$. Equivalently, $\langle\rho,\widetilde\xi\rangle=0$ for every $h\ne g$ and every $\xi\in\mathbb Q^{S_h}$ with $C_h\xi=0$.

**Definition 4.2 (deformation).** A *deformation of degree $\le R'$ with $Q$ nodes* for $F$ consists of the following data.
- Rows $P_{f,v}(T)\in\mathcal O[T]^{1\times n}$ for $f\in\mathcal F$, $v\in I_f$, with $P_{f,v}(0)=e_v$ and $\deg_TP_{f,v}\le R'$.
- For each node $t\in[Q]$, a measured partition $(S_{t,g},C_{t,g})_{g\in[r]}$ and labels $\sigma_f(t)\in[r]$.

The requirement is that every $P_{f,v}(t)$ is $\sigma_f(t)$-local for the partition at $t$.

**Lemma 4.3 (Chinese remainder at the nodes).** Let $A$ be a commutative ring, $J\ge1$, and $t_1,\dots,t_Q\in A$ with every difference $t_a-t_b$ ($a\ne b$) a unit. Then the $A$-linear map
$$A[T]_{<QJ}\to\prod_{a=1}^QA[z]/(z^J),\qquad H\mapsto\bigl(H(t_a+z)\bmod z^J\bigr)_a$$
is bijective. Its matrix $V_{(a,j),u}=\binom uj t_a^{u-j}$ is therefore invertible.

If $A=\mathbb Z/p^\mu$, then $H(0)$ is computed from the jets with $O((QJ)^3)$ ring operations.

*Proof.*
- **Comaximal ideals.** $(T-t_a)-(T-t_b)$ is a unit, so the ideals $\big((T-t_a)^J\big)$ are pairwise comaximal. Their intersection is therefore the ideal generated by $N=\prod_a(T-t_a)^J$.
- **CRT.** The Chinese remainder theorem gives $A[T]/(N)\cong\prod_aA[T]/((T-t_a)^J)$.
- **Identifications.** Division by the monic $N$ of degree $QJ$ identifies $A[T]/(N)$ with $A[T]_{<QJ}$. The substitution $T=t_a+z$ identifies $A[T]/((T-t_a)^J)$ with $A[z]/(z^J)$. So the map is bijective.
- **Computation.** For $A=\mathbb Z/p^\mu$, $V$ is invertible, so it is invertible modulo $p$. Hence each column of the remaining submatrix has an entry nonzero modulo $p$, that is, a unit of $A$. Gaussian elimination pivoting on such entries solves $Vh=(\text{jets})$, and $H(0)=h_0$. ∎

**Lemma 4.4 (jets by Boolean substitution).** Let $\mathcal R$ be a commutative ring and $q\ge0$. For each $x\in\{0,1\}^n$ let $p_x\in\mathcal R[X_1,\dots,X_n]$ have all monomials of degree $\le q$. Then
$$\sum_{x\in\{0,1\}^n}p_x(x)=\sum_{E\subseteq[n],\ |E|\le q}\ \sum_{U\subseteq E}(-1)^{|E\setminus U|}\sum_{x:\ x_E=1}p_x(\mathbf 1_U).$$

*Proof.*
- **Multilinearize.** Replace every $X_v^a$ ($a\ge1$) by $X_v$; this does not raise degrees or change values on $\{0,1\}^n$. Write $p_x\equiv\sum_{|E|\le q}c_E(x)X^E$, where $X^E=\prod_{v\in E}X_v$.
- **Boolean values.** $x^E=[x_E=1]$, so $p_x(x)=\sum_Ec_E(x)[x_E=1]$.
- **Möbius inversion** on the Boolean lattice gives $c_E(x)=\sum_{U\subseteq E}(-1)^{|E\setminus U|}p_x(\mathbf 1_U)$. Indeed $p_x(\mathbf 1_U)=\sum_{E'\subseteq U}c_{E'}(x)$, and $\sum_{E'\subseteq U\subseteq E}(-1)^{|E\setminus U|}$ is $1$ if $E'=E$ and $0$ otherwise.

Substitute and exchange the sums. ∎

**Lemma 4.5 (separation by measurements).** Let:
- $\mathcal R$ be a commutative ring and $(S_g,C_g)_{g\in[r]}$ a measured partition with row norms $\le\beta$;
- $E\subseteq[n]$;
- $A_g:\mathbb Z^n\to\mathcal R$ ($g\in[r]$) satisfy $A_g(x+\widetilde\xi)=A_g(x)$ whenever $h\ne g$, $\xi\in\mathbb Z^{S_h}$ and $C_h\xi=0$.

For $y_h\in\mathbb Z^{c_h}$, define the fibre $\mathcal X_h(y_h)=\{x_h\in\{0,1\}^{S_h}:C_hx_h=y_h,\ x_h=1\text{ on }E\cap S_h\}$. Then
$$\sum_{x\in\{0,1\}^n,\ x_E=1}\ \prod_gA_g(x)=\sum_y\ \prod_g\ \sum_{x_g\in\mathcal X_g(y_g)}A_g\Bigl(\widetilde{x_g}+\sum_{h\ne g}\widetilde{\bar x_h}\Bigr).$$
Here $y=(y_h)_h$ runs over $\prod_h(\mathbb Z\cap[-\beta,\beta])^{c_h}$ with every fibre nonempty, and $\bar x_h\in\mathcal X_h(y_h)$ is any fixed element.

*Proof.*
- **Range of $y$.** Entries of $C_hx_h$ lie in $[-\beta,\beta]$ for Boolean $x_h$.
- **Bijection.** The map $x\mapsto\bigl((C_hx|_{S_h})_h,(x|_{S_g})_g\bigr)$ is a bijection from $\{x:x_E=1\}$ onto the pairs $(y,(x_g))$ with $x_g\in\mathcal X_g(y_g)$.
- **Invariance.** For such $x$ and $h\ne g$, $x|_{S_h}-\bar x_h\in\mathbb Z^{S_h}\cap\ker C_h$. Applying the invariance once for each $h\ne g$ gives $A_g(x)=A_g(\widetilde{x_g}+\sum_{h\ne g}\widetilde{\bar x_h})$.

So the summand depends on $x$ only through $y$ and $x_g$, factor by factor, and the sum factorizes. ∎

**Theorem 4.6 (deformation counting).** Let $F$ have a deformation of degree $\le R'$ with $Q<p$ nodes. Assume at every node:
- $r$ groups of size $\le s$;
- at most $c_*$ measurement rows in total, of norm $\le\beta$.

Put $q=\lfloor dR'/Q\rfloor$ and $J=q+1$. Given the deformation, $\#F\bmod p^\mu$ is computed with at most
$$Q\cdot C_q(n)\cdot(2\beta+1)^{c_*}\cdot2r\,2^{s}\cdot\mathrm{poly}(n,|\mathcal F|,J,\mu\log p)+\mathrm{poly}(Q,J,\mu\log p)$$
bit operations and polynomial space.

*Proof.*
1. **The global polynomial.** Let
   $$H(T)=\sum_{x\in\{0,1\}^n}\prod_f\Phi_{C_f}\bigl((\langle P_{f,v}(T),x\rangle)_{v\in I_f}\bigr)\in\mathcal O[T].$$
   Then $\deg H\le\sum_f|C_f|R'=dR'<QJ$, and $H(0)=\sum_x\prod_f\Phi_{C_f}(x|_{I_f})=\#F$. Let $\bar H=\pi(H)\in A[T]$.
2. **Jets at node $t$.** Let $\mathcal R=A[z]/(z^J)$ and $\Delta_{f,v}=P_{f,v}(t+z)-P_{f,v}(t)$, a row with every entry in $z\,\mathcal O[z]$. Define
   $$p_x(X)=\prod_f\Phi_{C_f}\bigl((\pi\langle P_{f,v}(t),x\rangle+\langle\pi\Delta_{f,v},X\rangle)_v\bigr)\in\mathcal R[X].$$
   - **Value.** $p_x(x)$ is the image of the $x$-summand of $H(t+z)$, so $\sum_xp_x(x)=\bar H(t+z)\bmod z^J$.
   - **Degree.** A monomial $X^\alpha$ of $p_x$ carries a coefficient divisible by $z^{|\alpha|}$, which vanishes in $\mathcal R$ when $|\alpha|\ge J$. So $p_x$ has degree $\le q$.
   
   Lemma 4.4 expresses the jet through the sums $\sum_{x:\,x_E=1}p_x(\mathbf 1_U)$.
3. **Separation.** For fixed $(E,U)$, $p_x(\mathbf 1_U)=\prod_gA_g(x)$, where
   $$A_g(x)=\prod_{f:\,\sigma_f(t)=g}\Phi_{C_f}\Bigl(\bigl(\pi\langle P_{f,v}(t),x\rangle+\langle\pi\Delta_{f,v},\mathbf 1_U\rangle\bigr)_v\Bigr).$$
   - **Invariance.** $P_{f,v}(t)$ is $g$-local for every factor in the product. So for $h\ne g$ and $\xi\in\mathbb Z^{S_h}\cap\ker C_h$ we get $\langle P_{f,v}(t),\widetilde\xi\rangle=0$ in $\mathbb Q$, hence in $A$, and $A_g$ is invariant.
   - **Evaluation.** Lemma 4.5 evaluates the sum. For each of the at most $(2\beta+1)^{c_*}$ tuples $y$, finding representatives and the $r$ group sums costs at most $2r2^s$ evaluations of some $A_g$, each a product of $\le|\mathcal F|$ factors in $\mathcal R$.
4. **Interpolation.** The differences of the nodes $1,\dots,Q$ are units of $A$ because $p>Q$. Lemma 4.3 recovers $\bar H(0)=\#F\bmod p^\mu$ from the $Q$ jets.
5. **Cost and space.** The loops over $t$, $(E,U)$ and $y$ give the time bound. Only the loop indices, one group sum, one partial product and the jet accumulators are stored. ∎

---

## 5. Swap interpolation

Fix integers $r\ge2$ and $Q\ge1$, and set up:
- the labels $\mathcal B=([r]\times[r])^Q$, with $D=|\mathcal B|=r^{2Q}$, writing $u(s)=(a_u(s),d_u(s))$;
- the involution $\theta_s$ of $\mathcal B$ that swaps the pair at position $s$, and the permutation matrix $\Pi_s$ with $(\Pi_s)_{u,y}=[y=\theta_su]$;
- for $S\subseteq[Q]$, the products $\theta_S$ and $\Pi_S$, and the complement $S^c=[Q]\setminus S$;
- $p_r=\tfrac{r-1}{2r}$ and $W\sim\mathrm{Bin}(Q,p_r)$.

**Lemma 5.1 (eigenbasis).** On $\mathbb Z^{[r]\times[r]}$ take the *local vectors*:
- $e_{\alpha\alpha}$ for $\alpha\in[r]$ (symmetric);
- $e_{\alpha\delta}+e_{\delta\alpha}$ for $\alpha<\delta$ (symmetric);
- $e_{\alpha\delta}-e_{\delta\alpha}$ for $\alpha<\delta$ (antisymmetric).

For local vectors $\psi_1,\dots,\psi_Q$ put $\psi(u)=\prod_s\psi_s(u(s))$, and let $j(\psi)$ be the number of antisymmetric factors. Then:
1. The $D$ vectors $\psi$ form an orthogonal basis of $\mathbb Q^{\mathcal B}$, with entries in $\{0,\pm1\}$ and $\lVert\psi\rVert_1\le2^Q$.
2. $\Pi_s\psi=\varepsilon_s\psi$, where $\varepsilon_s=-1$ if $\psi_s$ is antisymmetric and $\varepsilon_s=+1$ otherwise.
3. The number of $\psi$ with $j(\psi)=j$ is $N_j=\binom Qj\binom r2^j\binom{r+1}2^{Q-j}=D\Pr[W=j]$.
4. If $\psi_t$ is supported on $\{(\alpha,\alpha)\}$ or $\{(\alpha,\delta),(\delta,\alpha)\}$, then every $u\in\operatorname{supp}\psi$ has $a_u(t)\in\{\alpha,\delta\}$.

*Proof.*
1. The local vectors are pairwise orthogonal, $r^2$ in number, with entries in $\{0,\pm1\}$ and $\ell^1$-norms $\le2$. Tensor products inherit all of this.
2. $(\Pi_s\psi)(u)=\psi(\theta_su)$, and the swap acts on $\psi_s$ by $\varepsilon_s$.
3. There are $\binom r2$ antisymmetric and $r+\binom r2=\binom{r+1}2$ symmetric local vectors, and $\binom r2/r^2=p_r$.
4. This is immediate from the supports. ∎

**Krawtchouk polynomials.** $K_k(x)=\sum_{i=0}^k(-1)^i\binom xi\binom{Q-x}{k-i}$. It is a polynomial of degree $k$ with leading coefficient $(-2)^k/k!$, and it takes integer values at integers.

**Lemma 5.2.** If $j(\psi)=j$, then $\sum_{S\subseteq[Q],\,|S|=k}\Pi_{S^c}\psi=(-1)^jK_k(j)\psi$ for every $k\ge0$.

*Proof.* Let $X$ be the set of antisymmetric positions, $|X|=j$. By Lemma 5.1(2), $\Pi_{S^c}\psi=(-1)^{|X\setminus S|}\psi=(-1)^j(-1)^{|X\cap S|}\psi$. Counting the sets $S$ by $i=|X\cap S|$ gives $\sum_i(-1)^i\binom ji\binom{Q-j}{k-i}=K_k(j)$. For $k>Q$ both sides are $0$. ∎

**Lemma 5.3 (integrality).** Let $\mathcal I=\{j_0,\dots,j_0+R-1\}$ be $R\ge1$ consecutive integers. The unique polynomial $f$ of degree $<R$ with $f(j)=(-1)^j$ on $\mathcal I$ satisfies $f=\sum_{k<R}\gamma_kK_k$ with $\gamma_k\in\mathbb Z[\tfrac12]$.

*Proof.*
1. **$f$ is integer-valued.** $f(x)=(-1)^{j_0}\sum_{a<R}(-2)^a\binom{x-j_0}a$, since at $x=j_0+i$ the sum is $(1-2)^i$. So $f$ takes integer values at integers.
2. **Binomial basis.** An integer-valued polynomial has integer coordinates $\Delta^af(0)=\sum_i(-1)^{a-i}\binom aif(i)$ in the basis $\binom xa$ (Newton's forward-difference formula). So $f=\sum_a\beta_a\binom xa$ with $\beta\in\mathbb Z$, and $K_k=\sum_{a\le k}\kappa_{k,a}\binom xa$ with $\kappa\in\mathbb Z$ and $\kappa_{k,k}=(-2)^k$.
3. **Solve.** The system $\sum_k\gamma_k\kappa_{k,a}=\beta_a$ is triangular with diagonal $(-2)^k$. ∎

**Theorem 5.4 (swap interpolation, closed form).** Let $f,\gamma$ be as in Lemma 5.3, and $\varphi_S(T)=\prod_{s\in S}(1-T/s)$. Define
$$M_0=\sum_{k<R}\gamma_k\!\!\sum_{|S|=k}\!\Pi_{S^c},\qquad P(T)=I-M_0+\sum_{k<R}\gamma_k\!\!\sum_{|S|=k}\!\varphi_S(T)\,\Pi_{S^c}.$$
Let $C\in\mathbb Z^{c\times\mathcal B}$ have as rows the $\psi^{\mathsf T}$ with $j(\psi)\notin\mathcal I$, so $c=D\Pr[W\notin\mathcal I]$. Then:

- **(B1)** $P(0)=I$, $\deg_TP\le R-1$, and all coefficients lie in $\mathbb Z[\frac1{2\cdot Q!}]$.
- **(B2)** For $t\in[Q]$, the $(u,y)$ entry of $P(t)-(I-M_0)$ vanishes unless $a_y(t)=d_u(t)$.
- **(B3)** If $C\xi=0$ then $M_0\xi=\xi$.
- **(B4)** Let $t\in[Q]$, $g\in[r]$ and $d_u(t)\ne g$. If $\xi$ is supported on $\{y:a_y(t)=g\}$ and $C\xi=0$, then $P(t)_u\,\xi=0$.

*Proof.*
- **(B1)** $\varphi_S(0)=1$ and $\deg\varphi_S=|S|$. The denominators of $\varphi_S$ divide $Q!$, and $\gamma\in\mathbb Z[\frac12]$.
- **(B2)** If $t\in S$ then $\varphi_S(t)=0$. Otherwise row $u$ of $\Pi_{S^c}$ has its only $1$ at $y=\theta_{S^c}u$. This $y$ swaps position $t$, so $a_y(t)=d_u(t)$.
- **(B3)** By Lemma 5.2, $M_0\psi=(-1)^jf(j)\psi$ with $j=j(\psi)$, and $(-1)^jf(j)=1$ for $j\in\mathcal I$. By orthogonality, $\ker C$ is spanned by the $\psi$ with $j(\psi)\in\mathcal I$.
- **(B4)** By (B3), $(I-M_0)\xi=0$. By (B2), $P(t)_u-(I-M_0)_u$ is supported on $\{y:a_y(t)=d_u(t)\}$, which misses $\operatorname{supp}\xi$. ∎

**Lemma 5.5 (Hoeffding for the binomial).** For $a>0$, $\Pr[|W-Qp_r|\ge a]\le2e^{-2a^2/Q}$.

*Proof.*
1. **One coordinate.** For one Bernoulli($p$) coordinate $Z$, let $L(\lambda)=\ln\mathbb Ee^{\lambda(Z-p)}=-\lambda p+\ln(1-p+pe^\lambda)$. Then $L(0)=L'(0)=0$ and $L''=\vartheta(1-\vartheta)\le\frac14$, where $\vartheta=pe^\lambda/(1-p+pe^\lambda)$. So $L(\lambda)\le\lambda^2/8$.
2. **Markov.** By independence, $\Pr[W-Qp\ge a]\le e^{-\lambda a+Q\lambda^2/8}$. Take $\lambda=4a/Q$; the other tail is symmetric. ∎

**Corollary 5.6 (interval choice).** Let $R\ge4$ be even and $j_0=\lceil Qp_r\rceil-R/2$. Then $\mathcal I\supseteq\{j\in\mathbb Z:|j-Qp_r|\le R/4\}$, so
$$c/D=\Pr[W\notin\mathcal I]\le2e^{-R^2/(8Q)}.$$

*Proof.*
- **Lower end.** $j_0<Qp_r+1-R/2\le Qp_r-R/4$.
- **Upper end.** $j_0+R-1\ge Qp_r+R/2-1\ge Qp_r+R/4$.
- **Tail.** Apply Lemma 5.5 with $a=R/4$. ∎

---

## 6. Transfer to a formula, and balanced schedules

Let $F$ have clauses $C_f$ with scopes $I_f$. Set up:
- **incidences** $(f,v)$ with $v\in I_f$; there are $M$ of them;
- $b$, the largest occurrence number, and $w$, the largest width;
- **schedules** $a_v\in[r]^Q$ for each variable and $d_f\in[r]^Q$ for each clause;
- the **label** $u_{f,v}=(a_v(s),d_f(s))_{s\in[Q]}\in\mathcal B$, with $c_u$ incidences carrying label $u$ and $h=\max_uc_u$.

Number the incidences with label $u$ as $0,\dots,c_u-1$; this is the **copy** $\kappa(f,v)\in[h]$. Distinct incidences then hold distinct **positions** $(\kappa(f,v),u_{f,v})\in[h]\times\mathcal B$.

The **filling map** $\iota:\mathbb Q^n\to\mathbb Q^{[h]\times\mathcal B}$ is defined by $(\iota x)_{(\kappa,y)}=x_v$ if position $(\kappa,y)$ is held by an incidence of $v$, and $0$ otherwise. $\pi_\kappa$ projects onto copy $\kappa$.

**Lemma 6.1 (transfer).** Put
$$P_{f,v}(T)=P(T)_{u_{f,v}}\,\pi_{\kappa(f,v)}\,\iota,\qquad S_{t,g}=\{v:a_v(t)=g\},\qquad\sigma_f(t)=d_f(t),$$
and let $C_{t,g}$ consist of the nonzero rows among $(\psi^{\mathsf T}\pi_\kappa\iota)|_{S_{t,g}}$, for $\kappa\in[h]$ and rows $\psi^{\mathsf T}$ of $C$. Then:
- these data form a deformation of degree $\le R-1$ with $Q$ nodes;
- every measurement row has $\ell^1$-norm $\le2^Q$;
- each node has at most $2hc$ measurement rows in total.

*Proof.*
1. **Base point.** $P_{f,v}(0)=e_{u_{f,v}}^{\mathsf T}\pi_\kappa\iota=e_v^{\mathsf T}$.
2. **Norms.** $(\psi^{\mathsf T}\pi_\kappa\iota)_v=\sum_{y:\,(\kappa,y)\text{ held by }v}\psi(y)$, so the $\ell^1$-norm is at most $\lVert\psi\rVert_1\le2^Q$.
3. **Row count.** A position held by $v\in S_{t,g}$ has a label $y$ with $a_y(t)=a_v(t)=g$. By Lemma 5.1(4), $\psi^{\mathsf T}\pi_\kappa\iota$ vanishes on $S_{t,g}$ unless $g$ is one of the at most two indices of $\psi_t$. So each pair $(\kappa,\psi)$ contributes at most two rows per node.
4. **Locality.** Let $g\ne d_f(t)$ and $\xi\in\mathbb Q^{S_{t,g}}$ with $C_{t,g}\xi=0$. For each copy put $\eta_\kappa=\pi_\kappa\iota\widetilde\xi$.
   - **Support.** $\eta_\kappa$ is supported on positions held by variables in $S_{t,g}$, so on $\{y:a_y(t)=g\}$.
   - **Kernel.** For every row $\psi^{\mathsf T}$ of $C$, $\psi^{\mathsf T}\eta_\kappa=(\psi^{\mathsf T}\pi_\kappa\iota)|_{S_{t,g}}\,\xi=0$, so $C\eta_\kappa=0$.
   
   Since $d_{u_{f,v}}(t)=d_f(t)\ne g$, (B4) gives $P_{f,v}(t)\widetilde\xi=P(t)_{u_{f,v}}\eta_{\kappa(f,v)}=0$. ∎

**Lemma 6.2 (balanced schedules).** Assume $2(b+w)D^2\le M$ and $2Qr^2\le n$. Schedules with $h\le2M/D$ and $|S_{t,g}|\le2n/r$ for all $t,g$ can be computed deterministically with $\mathrm{poly}(n,M,D,Q,r)$ operations on rationals of polynomial bit length.

*Proof.*
1. **Random schedules.** Choose all coordinates $a_v(t),d_f(t)$ independently and uniformly in $[r]$. Put $n_{t,g}=|S_{t,g}|$ and
   $$V_1=\sum_{u\in\mathcal B}\bigl(c_u-M/D\bigr)^2,\qquad V_2=\sum_{t,g}\bigl(n_{t,g}-n/r\bigr)^2,\qquad\Psi=\frac{V_1}{M(b+w)}+\frac{V_2}{Qn}.$$
2. **Variance of $c_u$.** Write $c_u=\sum_\iota X_\iota$ with $X_{f,v}=[a_v=a(u)][d_f=d(u)]$, so $\mathbb EX_\iota=1/D$ and $\mathbb Ec_u=M/D$.
   - Two indicators sharing neither the variable nor the clause are independent.
   - Two sharing one of them have $\mathbb E[X_\iota X_{\iota'}]=r^{-3Q}=D^{-3/2}$.
   - There are at most $(b-1)M$ ordered pairs sharing a variable and $(w-1)M$ sharing a clause.
   
   Hence $\mathbb E(c_u-M/D)^2\le M/D+(b+w-2)MD^{-3/2}$, and $\mathbb EV_1\le M(b+w)$.
3. **Variance of $n_{t,g}$.** $n_{t,g}\sim\mathrm{Bin}(n,1/r)$, so $\mathbb EV_2=Qn(1-1/r)$. Hence $\mathbb E\Psi<2$.
4. **Conditional expectations.** Fix the coordinates one at a time, each to a value that does not increase $\mathbb E[\Psi\mid\text{fixed coordinates}]$; such a value exists because the conditional expectation is the average over the $r$ choices.
   - The conditional expectation is computable exactly: given fixed coordinates, the indicators above have explicit product probabilities, and pairs sharing nothing remain independent.
   - So $\mathbb E[V_1\mid\cdot]=\sum_u\bigl(\mathrm{Var}(c_u\mid\cdot)+(\mathbb E[c_u\mid\cdot]-M/D)^2\bigr)$ is a sum over $D$ labels of $O(M(b+w))$ products, and $V_2$ is similar.
5. **Conclusion.** At the end, $\Psi<2$. Hence:
   - $(c_u-M/D)^2<2M(b+w)$, so $c_u<M/D+\sqrt{2M(b+w)}\le2M/D$;
   - $(n_{t,g}-n/r)^2<2Qn$, so $n_{t,g}<n/r+\sqrt{2Qn}\le2n/r$. ∎

---

## 7. Counting with at most three occurrences

**Parameters.** For $n\ge2$ put
$$\lambda=\lceil\log n\rceil,\quad\nu=\lceil\log\lambda\rceil,\quad r=\lfloor\sqrt\lambda\rfloor,\quad Q=\Bigl\lfloor\frac\lambda{4\nu}\Bigr\rfloor,\quad R=2\lceil\sqrt\lambda\rceil,$$
with $p$ the least prime $>Q$ and $\mu=n+1$. Then $2^{\lambda-1}<n\le2^\lambda$ and $2^{\nu-1}<\lambda\le2^\nu$.

**Lemma 7.2 (parameter inequalities).** For $\lambda\ge2^{17}+1$ (so $\nu\ge18$):
- **(P0)** $3\le p\le2Q$ (Bertrand), so $p\nmid2\cdot Q!$.
- **(P1)** $r\ge2$.
- **(P2)** $R$ is even, and $6\le R\le Q/2$.
- **(P3)** $\log D\le\lambda/4$, so $12D^2<n$.
- **(P4)** $2Qr^2<n$.
- **(P5)** $R^2/(8Q)\ge2\nu$.
- **(P6)** $x:=3(R-1)/Q\le25\nu/\sqrt\lambda<e$.
- **(P7)** $24e^{-2\nu}(Q+2)\le1/\sqrt\lambda$.
- **(P8)** $2/r\le3/\sqrt\lambda$.

*Proof.*
- **(P1)** Clear.
- **(P2)** We have $R\le2\sqrt\lambda+2$ and $Q\ge\lambda/(4\nu)-1$, so it suffices that $g(\lambda):=\lambda/(8\nu)-2\sqrt\lambda-\frac52\ge0$.
  - On a block $2^{\nu-1}<\lambda\le2^\nu$, $g$ is increasing once $\sqrt\lambda\ge8\nu$, which holds for $\nu\ge18$.
  - At the left end $\lambda=2^{\nu-1}$, $g\ge0$ for $\nu\ge18$: the value is $910.2-724.1-2.5>0$ at $\nu=18$, and it increases with $\nu$.
- **(P3)** $r^2\le\lambda\le2^\nu$ and $4Q\nu\le\lambda$, so $\log D=2Q\log r\le Q\nu\le\lambda/4$. Then $12D^2\le12\cdot2^{\lambda/2}\le2^{\lambda-1}<n$ for $\lambda\ge10$.
- **(P4)** $2Qr^2\le\lambda^2/(2\nu)\le2^{\lambda-1}<n$.
- **(P5)** $R^2\ge4\lambda$ and $8Q\le2\lambda/\nu$.
- **(P6)** $x\le12\nu(2\sqrt\lambda+1)/(\lambda-4\nu)$, which is $\le25\nu/\sqrt\lambda$ iff $\lambda\ge12\sqrt\lambda+100\nu$. This holds for $\lambda\ge2^{17}$. Also $25\nu/\sqrt\lambda\le1.25<e$ there.
- **(P7)** $e^{-2\nu}\le\lambda^{-2/\ln2}$, so $24e^{-2\nu}(Q+2)\le24\lambda^{1-2.88}\le\lambda^{-1/2}$.
- **(P8)** $r\ge\sqrt\lambda-1\ge\frac23\sqrt\lambda$. ∎

The script checks the integer forms of (P1)–(P3), (P5), (P6) and (P8) exhaustively for $2^{17}<\lambda\le2^{21}$ and at both ends of every block up to $\nu=400$. It checks (P7) and the exponent bounds below on a grid up to $2^{120}$. Just below $2^{17}$ only (P2) fails, so the threshold is where the stated argument needs it.

**Theorem 7.1.** There is a deterministic algorithm that, given a CNF $F$ on $n$ variables in which every clause has width $\le3$ and every variable occurs in $1$ to $3$ clauses, outputs $\#F$ using $\mathrm{poly}(n)$ space and the following time:
- for $n>2^{2^{17}}$: $\mathrm{poly}(n)\cdot2^{13\,n\nu^2/\sqrt\lambda}\le\mathrm{poly}(n)\cdot2^{13n(\log\log n+2)^2/\sqrt{\log n}}$ bit operations;
- for smaller $n$: $2^n\mathrm{poly}(n)$, by exhaustive search.

*Proof.* Here $M=d=\sum_f|C_f|\in[n,3n]$ and $b,w\le3$.
1. **Algorithm.**
   - Compute the schedules of Lemma 6.2, whose hypotheses are (P3) and (P4).
   - Compute $\gamma$ by Lemma 5.3, then $P(T)$ and $C$ by Theorem 5.4 with the interval of Corollary 5.6 (by (P2), $R\ge4$ is even).
   - Form the deformation of Lemma 6.1, with entries reduced modulo $p^\mu$; this is legitimate by (B1) and (P0).
   - Run Theorem 4.6.
   
   The tables have $D^2\le(2n)^{1/2}$ entries, each a sum of at most $2^Q\le D$ terms, so preprocessing is polynomial.
2. **Correctness.** Theorem 4.6 returns $\#F\bmod p^\mu$. Since $\#F\le2^n<3^{n+1}\le p^\mu$, this is $\#F$.
3. **Time.** The parameters of Theorem 4.6 are:
   - $R'=R-1$;
   - $q\le3n(R-1)/Q=xn$;
   - $s\le2n/r$ (Lemma 6.2);
   - $\beta=2^Q$;
   - $c_*\le2hc\le2\cdot\frac{2M}D\cdot2De^{-R^2/(8Q)}\le24n\,e^{-2\nu}$, by Lemma 6.1, Lemma 6.2, Corollary 5.6 and (P5).

   Lemma A.1 gives $\log C_q(n)\le n\,x\log(e^2/x)$, using $x<e$ from (P6). Hence
   $$\frac{\log(\mathrm{time})}n\le\underbrace{x\log\frac{e^2}x}_{\text{Taylor jets}}+\underbrace{24e^{-2\nu}(Q+2)}_{\text{measurements}}+\underbrace{\frac2r}_{\text{groups}}+O\Bigl(\frac{\log n}n\Bigr).$$
   - **Taylor jets.** $x\log(e^2/x)$ increases on $(0,e)$. With (P6), it is at most $\frac{25\nu}{\sqrt\lambda}\log\frac{e^2\sqrt\lambda}{25\nu}\le\frac{25\nu}{\sqrt\lambda}\cdot\frac\nu2$.
   - **Measurements and groups.** By (P7) and (P8) these contribute $\le4/\sqrt\lambda$.
   
   The total is $\le(12.5\nu^2+4)/\sqrt\lambda\le13\nu^2/\sqrt\lambda$. Finally, $\nu\le\log\log n+2$ and $\lambda\ge\log n$.
4. **Space.** Every stored object (schedules, tables, matrices, loop indices, accumulators) has polynomial size. ∎

**Remark (where the bound becomes nontrivial).** The exact bound $x\log(e^2/x)+24e^{-2\nu}(Q+2)+2/r$ is below $1$ for every $\lambda>9{,}297{,}119$.
- It is checked exhaustively for $\lambda\le2\cdot10^7$, which covers the block $\nu=24$.
- For $\nu\ge25$, on each block the majorant $\hat x\log(e^2/\hat x)+2/(\sqrt\lambda-1)$, with $\hat x=12\nu(2\sqrt\lambda+1)/(\lambda-4\nu)\ge x$, decreases. At the block's start it is at most $0.83$, and the measurement term is below $10^{-14}$.

The paper's parameters reach this point only near $\lambda\approx7\cdot10^{10}$. Both thresholds are artifacts of crude constants (see [07](07_TOOLKIT.md) E4).

---

## 8. Proof of the Main Theorem

Let $F$ be a $k$-CNF on $n$ variables and $m$ clauses. For $n$ below a constant $n_0(k)$, use exhaustive search. Otherwise:
1. **Normalize** in $\mathrm{poly}(n+m)$ time. If the empty clause is present, output $0$. Now $m\le2(2n)^k$.
2. **Sparsify.** Put $\sigma=1/\lceil(\log n)^{1/(2k)}\rceil$, take $M$ from Lemma 2.5, and traverse the tree of Proposition 2.4 depth first. It has at most $(Cn+1)^22^{\sigma n}$ nodes, $C=\mathrm{poly}_k(\log n)$.
3. **Count each leaf $L$.** Proposition 3.3 gives $L'$ on $N\le6kBn$ variables with $\#L=2^z\#L'$. Count $L'$ by Theorem 7.1, or by exhaustive search if $N\le2^{2^{17}}$, at cost $O(1)$.
4. **Add** the leaf counts (Proposition 2.4(4)) in an accumulator of $n+1$ bits.

**Time.** Let $g(x)=13x(\log\log x+2)^2/\sqrt{\log x}$, which is increasing for $x\ge4$. Each leaf costs $O(1)+\mathrm{poly}(n)2^{g(6kBn)}$.

By Lemma 2.5:
- $M=O_k\big((\log n)^{1/(2k)}\log\log n\big)$;
- $B\le4M^{k-1}=O_k\big((\log n)^{(k-1)/(2k)}(\log\log n)^{k-1}\big)$.

For large $n$, $6kB\le n$, so $\log(6kBn)\le2\log n$ and
$$g(6kBn)\le78k\,Bn\,\frac{(\log\log n+3)^2}{\sqrt{\log n}}=O_k\Bigl(n\,\frac{(\log\log n)^{k+1}}{(\log n)^{1/(2k)}}\Bigr).$$
Also $\sigma n\le n/(\log n)^{1/(2k)}$. Multiplying by the node count gives the time bound.

**Space.** Depth-first search stores one path, one leaf computation and the accumulator, all polynomial (Proposition 2.4(5) and Theorem 7.1). ∎

---

## 9. What this proof covers, and what it does not

- **Covered, with every step proved:** subexponential exact counting for $k$-CNF; the refutation of ETH and #ETH; polynomial space.
- **Not re-proved: the paper's §5–6.** That route is lifted parity tests, Tanner codes on $\varepsilon$-biased Cayley graphs, and cross-equation 2-adic truncation. It removes the sparsification loss and gives exponent $n(\log\log n)^2/\sqrt{\log n}$ for every $k$. This is the theorem OpenAI formalized; its Lean proof builds and is audited in [06](06_LEAN_STATEMENT_AUDIT.md). The new ingredient, the cross-equation truncation, is proved for every prime in [07](07_TOOLKIT.md) E2.
- **Not re-proved: App. D**, which removes one $\log\log n$ factor (verified numerically in [`results/appendix_d.json`](results/appendix_d.json); ~80% in [01](01_VERDICT.md)).
- **Why the exponent here is weaker.** Tovey's reduction multiplies the number of variables by $B=\mathrm{polylog}(n)$, and the Taylor cost is linear in the total degree. Removing that factor needs the §5–6 degree compression.

## 10. Cross-references

| This file | Paper | Closest Lean counterpart (`OAI/Algorithms/SubexponentialSAT/`) | Check |
|---|---|---|---|
| Lemma 2.1 | Lemma A.1 | `SunflowerBranch` (children, antichain, partition) | Part S |
| Lemma 2.2 | Lemma A.2 | `InsertionInvariant` | Part S |
| Lemma 2.3 | Lemma A.3 | `RemovalCharging` | Part S |
| Prop. 2.4, Lemma 2.5 | App. A.3 | `QuantitativeBranchTree`, `SparsificationConstants` | Parts S, A |
| Lemmas 3.1–3.2 | §2.8 | none (the formal route avoids Tovey) | Part R |
| Lemma 4.3 | §3, B.3 | `JetCRT`, `SubexpCRTAtZero` | Part K |
| Lemma 4.4 | Lemma 3.3 | `Counting`, `Moments` (different coefficient extraction) | Part K |
| Lemma 4.5 | Lemma 3.2 | `Separation`, `SubexpSeparation`, `SubexpAlgorithm` | Part K |
| Lemma 5.1 | Lemma 4.1 | `LocalSpectrum`, `Subexp/TensorSpectrum`, `BinomialMultiplicity`, `SubexpSwapMultiplicity` | Part K |
| Lemmas 5.2–5.3, Thm 5.4 | Lemmas 4.1–4.3, (4.10) | `SpectralRank`, `RowInterpolation`, `Annihilation` | Part K; [`verify_distilled_swap_formula.py`](scripts/verify_distilled_swap_formula.py) |
| Cor. 5.6 | Lemma 4.8 | `Subexp/TailInterval` (a second-moment tail with larger $R$) | Part K |
| Lemma 6.1 | §4.5 | `SlotAssignment`, `Layout`, `SubexpAnnihilation` | Part K |
| Lemma 6.2 | Lemma 4.5 | `ConditionalExpectation`, `Variance`, `TypeDistribution`, `SubexpTypeVariance` | Part K |
| Lemma 7.2, Thm 7.1 | Lemma 4.8, (2.19), App. B | `SubexpSchedule`, `SubexpParameters`; arithmetic in `ActualLiftRing` | Parts K, P |

The Lean route differs in its details: multicolour types, Galois-ring arithmetic, packing, and §5–6 in place of Tovey. The table maps mathematical content, not proofs line by line.

---

## Appendix A. An elementary estimate

**Lemma A.1.** For $1\le q\le n$, $C_q(n)\le(e^2n/q)^q$. Hence, for integers $0\le q\le xn$ with $0<x\le e$, $\log C_q(n)\le n\,x\log(e^2/x)$.

*Proof.* For $0<\vartheta\le1$,
$$C_q(n)\le\vartheta^{-q}\sum_j\binom nj(2\vartheta)^j=\vartheta^{-q}(1+2\vartheta)^n.$$
Take $\vartheta=q/n$ and use $(1+2q/n)^n\le e^{2q}$. For the second claim:
- the function $q\mapsto q\log(e^2n/q)$ increases for $q<en$, and $C_0(n)=1$;
- if $q>n$ then $x>1$ and $C_q(n)=3^n$, while $x\log(e^2/x)\ge\log e^2>\log3$ for $1\le x\le e$. ∎
