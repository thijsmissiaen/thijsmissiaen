---
name: second-order-debranges-lehmer-lifting
description: Scientific knowledge skill for Paper XIII of the Transverse Wronskian series on Riemann's Xi-function. Contains the proven analytical identities, 60-decimal mpmath and 550-decimal Lehmer deformation tables, and 35-decimal mpmath.iv interval enclosures for the second-order Hermite-Biehler de Branges structure function E_{alpha,beta}(z) = Xi(z) - i*alpha*Xi'(z) - beta*Xi''(z) (Lagarias 2006; Ki 2008; Suzuki 2011), the exact cross-term-free reproducing kernel splitting K_{alpha,beta}(t-id, t-id) == (alpha/pi)*W_tilde_{0,1}(d,t) + (alpha*beta/pi)*W_tilde_{1,2}(d,t), the d = 0 lifting of the Lehmer near-coalescence degeneracy at order d^0, the radial 4D-Laplace PDE and universal dipole recursion c_{2m-1}^{(1)}(t) == (1/4)*(2m(2m+1)*c_{2m+1}(t) + c_{2m-1}''(t)), and the universal order-M Hurwitz chain coalescence law pi*K_alpha(t_0, t_0) = +alpha_{m-1}*alpha_m*(Xi^{(m)}(t_0))^2 > 0.
---

# Second-Order Hermite–Biehler De Branges Kernels, the Derivative-Wronskian Splitting $K _ {\alpha,\beta} \equiv \frac{\alpha}{\pi}\widetilde{W} _ {0,1} + \frac{\alpha\beta}{\pi}\widetilde{W} _ {1,2}$, and $d = 0$ Lifting of the Lehmer Near-Coalescence Degeneracy (Paper XIII)

## 1. Scope and Mathematical Architecture

This skill encapsulates the analytical theorems, exact reproducing-kernel decompositions, directed-rounding `mpmath.iv` interval-arithmetic enclosures, and GPU/vectorized continuum verification procedures of **Paper XIII** in the series on the Transverse Wronskian of Riemann's $\Xi$-Function (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Paper XIII resolves how the second-order (and general order $M$) canonical Hermite–Biehler de Branges structure functions studied by Lagarias (2006), Ki (2008), and Suzuki (2011) interact with the parity-reduced transverse Wronskian hierarchy and lift the critical-line Lehmer near-coalescence degeneracy already at order $d^0$:

1. **The Exact Second-Order De Branges Reproducing Kernel Identity and Cauchy–Pick Gram Splitting (Section 2, Theorem 2.1 & Proposition 2.2)**: For the second-order Hermite–Biehler structure function $E _ {\alpha,\beta}(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z) - \beta\Xi''(z)$ ($\alpha > 0$, $\beta > 0$) and its Schwarz reflection $E _ {\alpha,\beta}^{\sharp}(z) := \Xi(z) + \mathrm{i}\alpha\Xi'(z) - \beta\Xi''(z)$, evaluating the numerator $\lvert E _ {\alpha,\beta}(w)\rvert^2 - \lvert E _ {\alpha,\beta}^{\sharp}(w)\rvert^2$ ($w = t - \mathrm{i}d$) annihilates not only the three quadratic squares $\lvert\Xi(w)\rvert^2$, $\alpha^2\lvert\Xi'(w)\rvert^2$, and $\beta^2\lvert\Xi''(w)\rvert^2$, but also the potentially sign-indefinite $(0,2)$ cross-term $-2\beta\mathrm{Re}(\overline{\Xi(w)}\Xi''(w)) \equiv 0$. Consequently, the diagonal reproducing kernel splits without cross-term contamination into the sum of the primary parity-reduced Wronskian $\widetilde{W} _ {0,1}(d,t)$ of $\Xi$ and the parity-reduced Lucas–Jensen derivative Wronskian $\widetilde{W} _ {1,2}(d,t)$ of $\Xi'$:

$$K _ {\alpha,\beta}(t - \mathrm{i}d, t - \mathrm{i}d) := \frac{\lvert E _ {\alpha,\beta}(t - \mathrm{i}d)\rvert^2 - \lvert E _ {\alpha,\beta}^{\sharp}(t - \mathrm{i}d)\rvert^2}{4\pi d} \equiv \frac{\alpha}{\pi}\widetilde{W} _ {0,1}(d,t) + \frac{\alpha\beta}{\pi}\widetilde{W} _ {1,2}(d,t).$$

Moreover, any $M \times M$ off-diagonal Cauchy–Pick Gram matrix splits identically as $\mathbf{G} _ {\alpha,\beta}(z_1,\dots,z_M) \equiv \alpha\mathbf{G}^{(0,1)} + \alpha\beta\mathbf{G}^{(1,2)}$ with zero non-adjacent jump terms.

2. **Lifting the $d = 0$ Lehmer Near-Coalescence Degeneracy at Order $d^0$ and Global Square Interlocking (Section 3, Theorem 3.1)**: In the first-order case ($\beta = 0$, Paper IV), a hypothetical double critical-line zero ($\Xi(t_0) = \Xi'(t_0) = 0$, $\Xi''(t_0) \neq 0$) forces $K _ \alpha(t_0, t_0) = 0$ on $d = 0$, so strict positivity emerges only at order $d^2$ via $c_3(t_0) = \frac{1}{2}(\Xi''(t_0))^2 > 0$; similarly, at the midpoint $t _ {\mathrm{mid}}$ of a closely spaced Lehmer pair $t _ \pm = t _ {\mathrm{mid}} \pm \delta$ ($\delta \to 0$), $\pi K _ \alpha(t _ {\mathrm{mid}}, t _ {\mathrm{mid}}) = \alpha c_1(t _ {\mathrm{mid}}) = \mathcal{O}(\delta^2) \to 0$. For the second-order kernel ($\beta > 0$), taking $d \to 0^+$ yields $\pi K _ {\alpha,\beta}(t,t) \equiv \alpha c_1(t) + \alpha\beta c_1^{(1)}(t)$, where $c_1^{(1)}(t) = (\Xi''(t))^2 - \Xi'(t)\Xi'''(t)$. At any critical-line extremum $t _ {\mathrm{mid}}$ where $\Xi'(t _ {\mathrm{mid}}) = 0$, the product $-\Xi'(t _ {\mathrm{mid}})\Xi'''(t _ {\mathrm{mid}}) \equiv 0$ vanishes identically, collapsing $c_1^{(1)}(t _ {\mathrm{mid}})$ to the pure square $(\Xi''(t _ {\mathrm{mid}}))^2 > 0$ and keeping $K _ {\alpha,\beta}(t _ {\mathrm{mid}}, t _ {\mathrm{mid}})$ strictly positive on $d = 0$ itself (already at order $d^0$):

$$\pi K _ {\alpha,\beta}(t _ {\mathrm{mid}}, t _ {\mathrm{mid}}) \equiv -\alpha\Xi(t _ {\mathrm{mid}})\Xi''(t _ {\mathrm{mid}}) + \alpha\beta\left(\Xi''(t _ {\mathrm{mid}})\right)^2 \xrightarrow{\delta \to 0} +\alpha\beta\left(\Xi''(t_0)\right)^2 > 0.$$

3. **The Radial 4D-Laplace–Wigner Identity, Universal Dipole Recursion, and Order $M$ Hurwitz Chains (Section 4, Theorems 4.1 & 4.2)**: Because $(\partial_d^2 + \partial_t^2)\lvert\Xi\rvert^2 = 4\lvert\Xi'\rvert^2$, the derivative Wronskian $\widetilde{W} _ {1,2}(d,t)$ is coupled to $\widetilde{W} _ {0,1}(d,t)$ by the exact radial four-dimensional Laplace operator $\widetilde{W} _ {1,2}(d,t) \equiv \frac{1}{4}(\partial_d^2 + \frac{2}{d}\partial_d + \partial_t^2)\widetilde{W} _ {0,1}(d,t)$, which induces the universal dipole recursion $c _ {2m-1}^{(1)}(t) \equiv \frac{1}{4}(2m(2m+1)c _ {2m+1}(t) + c _ {2m-1}''(t))$ for all $m \ge 1$ (in particular, $c_1^{(1)}(t) \equiv \frac{3}{2}c_3(t) + \frac{1}{4}c_1''(t)$ for $m = 1$). In Wigner phase space, $I^{(1)}(y,t) = \frac{1}{4}(y^2 + \partial_t^2)I(y,t)$, which enhances the Riemann–Siegel diagonal to $\sim \frac{1}{80}(\ln\frac{t}{2\pi})^5$ and applies quartic $\mathcal{O}(N(t)^{-5}) = \mathcal{O}(t^{-5/2})$ self-damping to the slowest boundary pair $(N-1, N)$. Finally, for any real Hurwitz polynomial $P_M(X) = \sum _ {r=0}^M \alpha_r X^r$ ($\alpha_0 = 1$, $\alpha_r > 0$), $E _ {\boldsymbol{\alpha}}(z) := P_M(-\mathrm{i}\partial_z)\Xi(z)$ satisfies the universal $d = 0$ coalescence law $\pi K _ {\boldsymbol{\alpha}}(t_0, t_0) = +\alpha _ {m-1}\alpha_m(\Xi^{(m)}(t_0))^2 + \mathcal{O}(\alpha _ {>m}) > 0$ for every multiplicity $1 \le m \le M$.

## 2. Canonical Identities and Theorems

### 2.1 Second-Order De Branges Splitting and Cauchy–Pick Gram Decomposition (Theorem 2.1 & Proposition 2.2)

For $w = t - \mathrm{i}d$ ($d \neq 0$) and real parameters $\alpha, \beta \in \mathbb{R}$, decompose $E _ {\alpha,\beta}(z) = A _ \beta(z) - \mathrm{i}\alpha B(z)$ and $E _ {\alpha,\beta}^{\sharp}(z) = A _ \beta(z) + \mathrm{i}\alpha B(z)$ into their real-entire even-order and odd-order parts $A _ \beta(z) := \Xi(z) - \beta\Xi''(z)$ and $B(z) := \Xi'(z)$. For any pair of points $z_1, z_2 \in \mathbb{C}$ with $z_1 \neq \bar{z} _ 2$, Hermitian polarization cancels $A _ \beta(z_1)\overline{A _ \beta(z_2)} + \alpha^2 B(z_1)\overline{B(z_2)}$ identically, yielding the two-point reproducing kernel identity

$$K _ {\alpha,\beta}(z_1, z_2) := \frac{E _ {\alpha,\beta}(z_1)\overline{E _ {\alpha,\beta}(z_2)} - E _ {\alpha,\beta}^{\sharp}(z_1)\overline{E _ {\alpha,\beta}^{\sharp}(z_2)}}{2\pi\mathrm{i}(z_1 - \bar{z} _ 2)} \equiv \frac{\alpha}{\pi}\frac{\Xi(z_1)\overline{\Xi'(z_2)} - \Xi'(z_1)\overline{\Xi(z_2)}}{z_1 - \bar{z} _ 2} + \frac{\alpha\beta}{\pi}\frac{\Xi'(z_1)\overline{\Xi''(z_2)} - \Xi''(z_1)\overline{\Xi'(z_2)}}{z_1 - \bar{z} _ 2}.$$

On the diagonal $z_1 = z_2 = t - \mathrm{i}d$, where $z_1 - \bar{z} _ 2 = -2\mathrm{i}d$, this reduces to

$$K _ {\alpha,\beta}(t - \mathrm{i}d, t - \mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W} _ {0,1}(d,t) + \frac{\alpha\beta}{\pi}\widetilde{W} _ {1,2}(d,t), \qquad \widetilde{W} _ {k,k+1}(d,t) := \frac{\mathrm{Im}\left(\overline{\Xi^{(k)}(t - \mathrm{i}d)}\Xi^{(k+1)}(t - \mathrm{i}d)\right)}{d} = \frac{1}{2d}\frac{\partial}{\partial d}\left\lvert\Xi^{(k)}(t - \mathrm{i}d)\right\rvert^2.$$

### 2.2 Critical-Line $d = 0$ Lehmer Coalescence Lifting and Square Interlocking (Theorem 3.1)

Expanding $\widetilde{W} _ {1,2}(d,t) = \sum _ {m=1}^{\infty} c _ {2m-1}^{(1)}(t)d^{2m-2}$ around $d = 0$ yields the critical-line identity

$$\pi K _ {\alpha,\beta}(t,t) \equiv \alpha\left(\left(\Xi'(t)\right)^2 - \Xi(t)\Xi''(t)\right) + \alpha\beta\left(\left(\Xi''(t)\right)^2 - \Xi'(t)\Xi'''(t)\right) = \alpha c_1(t) + \alpha\beta c_1^{(1)}(t).$$

Across the critical line $d = 0$, the three families of special points—simple zeros $\gamma_k$ ($\Xi(\gamma_k) = 0$), local extrema $\tau_k$ ($\Xi'(\tau_k) = 0$), and inflection points $\nu_k$ ($\Xi''(\nu_k) = 0$)—exhibit exact topological square interlocking:

$$\pi K _ {\alpha,\beta}(\gamma_k, \gamma_k) = \alpha\left(\Xi'(\gamma_k)\right)^2 + \alpha\beta c_1^{(1)}(\gamma_k), \qquad \pi K _ {\alpha,\beta}(\tau_k, \tau_k) = -\alpha\Xi(\tau_k)\Xi''(\tau_k) + \alpha\beta\left(\Xi''(\tau_k)\right)^2, \qquad \pi K _ {\alpha,\beta}(\nu_k, \nu_k) = \alpha\left(\Xi'(\nu_k)\right)^2 - \alpha\beta\Xi'(\nu_k)\Xi'''(\nu_k).$$

### 2.3 Radial 4D-Laplace PDE, Dipole Recursion, and Order $M$ Hurwitz Chains (Theorems 4.1 & 4.2)

The primary and derivative parity-reduced Wronskians satisfy the global elliptic PDE and Taylor coefficient recursion

$$\widetilde{W} _ {1,2}(d,t) \equiv \frac{1}{4}\left(\frac{\partial^2\widetilde{W} _ {0,1}}{\partial d^2} + \frac{2}{d}\frac{\partial\widetilde{W} _ {0,1}}{\partial d} + \frac{\partial^2\widetilde{W} _ {0,1}}{\partial t^2}\right), \qquad c _ {2m-1}^{(1)}(t) \equiv \frac{1}{4}\left(2m(2m+1)c _ {2m+1}(t) + c _ {2m-1}''(t)\right) \qquad (\forall m \ge 1).$$

For a general degree $M$ Hurwitz polynomial $P_M(X) = \sum _ {r=0}^M \alpha_r X^r = P _ {\mathrm{even}}(X^2) + X P _ {\mathrm{odd}}(X^2)$ ($\alpha_0 = 1$, $\alpha_r > 0$) and $E _ {\boldsymbol{\alpha}}(z) := P_M(-\mathrm{i}\partial_z)\Xi(z) = F _ {\mathrm{even}}(z) - \mathrm{i}F _ {\mathrm{odd}}(z)$, the critical-line kernel equals the 1D Wronskian $\pi K _ {\boldsymbol{\alpha}}(t,t) \equiv F _ {\mathrm{even}}'(t)F _ {\mathrm{odd}}(t) - F _ {\mathrm{even}}(t)F _ {\mathrm{odd}}'(t)$, which at any multiplicity $m$ zero $t_0 \in \mathbb{R}$ ($1 \le m \le M$) collapses to the strictly positive pure square

$$\pi K _ {\boldsymbol{\alpha}}(t_0, t_0) = +\alpha _ {m-1}\alpha_m\left(\Xi^{(m)}(t_0)\right)^2 + \mathcal{O}(\alpha _ {>m}) > 0.$$

## 3. Reproducible Verification Workflow

Run the standalone verification and figure-generation suite from the repository root:

```bash
python scripts/verify_paper13_second_order_debranges.py
```
