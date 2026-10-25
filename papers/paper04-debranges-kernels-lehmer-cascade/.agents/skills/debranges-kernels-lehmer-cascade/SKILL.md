---
name: debranges-kernels-lehmer-cascade
description: Exact De Branges reproducing kernel identity K_alpha(t-id, t-id) = (alpha/pi)*W_tilde(d,t), finite Cauchy-Pick Gram matrix divided differences, universal transverse square-cascade c_{2m-1}(t_0) = [m/(m!)^2]*(Xi^{(m)}(t_0))^2 > 0 for order-m critical-line zero coalescence and 2D Lehmer-pair quadratic bowls, and universal Leibniz-Wronskian identity with 0%-overshoot exponential prime-gauge sieves. Activate when analyzing De Branges spaces, Lehmer pairs, or multiplicative prime-gauge Wronskians for Riemann's Xi-function.
---

# De Branges Reproducing Kernels, Lehmer Square-Cascade, and Leibniz–Wronskian Identity (Paper IV)

**Author:** Thijs Missiaen (Independent Researcher, Belgium; ORCID: `0009-0007-3678-0810`)

This skill encodes the exact analytical identities, theorems, and verification procedures of **Paper IV** (*De Branges Reproducing Kernels, the Universal Transverse Square-Cascade for Lehmer Pairs, and the Leibniz–Wronskian Identity for Riemann's ?-Function*).

## 1. Canonical Normalizations

With $s = \frac{1}{2} + d + \mathrm{i}t$ ($d = \sigma - \frac{1}{2} \in \mathbb{R}$, $t \in \mathbb{R}$) and $z = -\mathrm{i}(s - \frac{1}{2}) = t - \mathrm{i}d$:

$$\xi(s) := \frac{1}{2}s(s-1)\pi^{-s/2}\Gamma\left(\frac{s}{2}\right)\zeta(s), \qquad \Xi(z) := \xi\left(\frac{1}{2} + \mathrm{i}z\right),$$

$$\Phi(u) := \sum _ {n=1}^{\infty}\left(2\pi^2 n^4 \mathrm{e}^{9u} - 3\pi n^2 \mathrm{e}^{5u}\right)\exp\left(-\pi n^2 \mathrm{e}^{4u}\right), \qquad \Xi(z) = 8\int_0^{\infty}\Phi(u)\cos(2zu) \mathrm{d}u.$$

Decomposing $\Xi(t - \mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$, the transverse Wronskian and its parity-reduced form are:

$$W(d,t) := \frac{1}{2}\partial_d\left|\Xi(t - \mathrm{i}d)\right|^2 = R(d,t)\partial_d R(d,t) + I(d,t)\partial_d I(d,t), \qquad \widetilde{W}(d,t) := \frac{W(d,t)}{d} = \sum _ {n=1}^{\infty}c _ {2n-1}(t) d^{2n-2}.$$

## 2. Pillar 1: Exact De Branges Reproducing Kernel Identity (Theorem 2.2 & Proposition 2.4)

For any $\alpha > 0$, define the Hermite–Biehler pair $E _ \alpha(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z)$ and $E _ \alpha^{\sharp}(z) := \overline{E _ \alpha(\bar{z})} = \Xi(z) + \mathrm{i}\alpha\Xi'(z)$. By the Cauchy–Riemann equations $\partial_d R = \partial_t I$ and $\partial_d I = -\partial_t R$, both $R^2 + I^2 = \left|\Xi\right|^2$ and $\alpha^2((\partial_d R)^2 + (\partial_d I)^2) = \alpha^2\left|\Xi'\right|^2$ cancel identically in $\left|E _ \alpha(t-\mathrm{i}d)\right|^2 - \left|E _ \alpha^{\sharp}(t-\mathrm{i}d)\right|^2 = 4\alpha W(d,t)$, yielding the exact global identity:

$$K _ \alpha(t - \mathrm{i}d, t - \mathrm{i}d) := \frac{\left|E _ \alpha(t - \mathrm{i}d)\right|^2 - \left|E _ \alpha^{\sharp}(t - \mathrm{i}d)\right|^2}{4\pi d} = \frac{\left|E _ \alpha^{\sharp}(t + \mathrm{i}d)\right|^2 - \left|E _ \alpha(t + \mathrm{i}d)\right|^2}{4\pi d} \equiv \frac{\alpha}{\pi} \widetilde{W}(d,t).$$

For any finite sequence of points $z_1, \dots, z_M \in \mathbb{C} _ -$ ($z_j = t_j - \mathrm{i}d_j$, $d_j > 0$), the $M \times M$ Hermitian Cauchy–Pick Gram matrix $\mathbf{G} _ \alpha(z_1,\dots,z_M) _ {j,k} := K _ \alpha(z_j, z_k)$ satisfies:

$$\left[\mathbf{G} _ \alpha(z_1,\dots,z_M)\right] _ {j,k} = \frac{\alpha}{\pi} \frac{\Xi(z_j)\overline{\Xi'(z_k)} - \Xi'(z_j)\overline{\Xi(z_k)}}{z_j - \bar{z} _ k}.$$

This decouples local off-critical zero exclusion ($\widetilde{W}(d_0, t_0) > 0$) from the global Sonine–De Branges inner-product positivity obstructions of Conrey and Li (2000).

## 3. Pillar 2: Universal Transverse Square-Cascade and 2D Lehmer Quadratic Bowl (Theorems 3.2 & 3.4)

At any critical-line zero $t_0 \in \mathbb{R}$ of multiplicity $m \ge 1$ ($\Xi(t_0) = \dots = \Xi^{(m-1)}(t_0) = 0$, $\Xi^{(m)}(t_0) \neq 0$), all lower transverse dipole moments vanish ($c_1(t_0) = \dots = c _ {2m-3}(t_0) = 0$) and the leading non-vanishing coefficient collapses to the strictly positive pure square:

$$c _ {2m-1}(t_0) = \frac{m}{(m!)^2}\left(\Xi^{(m)}(t_0)\right)^2 > 0.$$

Around the midpoint $t _ {\mathrm{mid}}$ of any closely spaced Lehmer pair $t _ \pm = t _ {\mathrm{mid}} \pm \delta$ ($\delta > 0$) with local offset $\tau = t - t _ {\mathrm{mid}}$ and $\Xi''(t _ {\mathrm{mid}}) \neq 0$, the parity-reduced Wronskian expands into an isotropic two-dimensional sum of three squares:

$$\widetilde{W}(d, t _ {\mathrm{mid}} + \tau) = \frac{1}{2}\left(\Xi''(t _ {\mathrm{mid}})\right)^2\left(\delta^2 + \tau^2 + d^2\right) + \mathcal{O}\left((\delta + \left|\tau\right| + d)^3\right) > 0.$$

## 4. Pillar 3: Universal Leibniz–Wronskian Identity and 0%-Overshoot Prime-Gauge Sieve (Theorems 4.1 & 4.4)

For any real-analytic multiplier $\mathcal{M}(z)$ and $\Xi _ {\mathcal{M}}(z) := \mathcal{M}(z)\Xi(z)$, the logarithmic-derivative decomposition yields the exact cross-term-free identity:

$$\widetilde{W} _ {\mathcal{M}}(d,t) \equiv \left|\mathcal{M}(t - \mathrm{i}d)\right|^2 \widetilde{W}(d,t) + \left|\Xi(t - \mathrm{i}d)\right|^2 \widetilde{W} _ {1,\mathcal{M}}(d,t).$$

Consequently, $\widetilde{W} _ {\mathcal{M}}(d_0, t_0) > 0$ implies $\Xi(t_0 - \mathrm{i}d_0) \neq 0$ unconditionally, and any parasitic zero $\mathcal{M}(z_0) = 0$ yields only a harmless local false negative ($\widetilde{W} _ {\mathcal{M}}(d_0, t_0) = 0$). Choosing the entire zero-free exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {2 \le d \le X}\frac{\Lambda(d)}{\ln d}d^{-s})$ yields residual Dirichlet coefficients $b_n^{\mathcal{E}}(X)$ satisfying the positive recursion:

$$b_1^{\mathcal{E}}(X) = 1, \qquad b_n^{\mathcal{E}}(X)\ln n = \sum _ {d \mid n,  d > X} b _ {n/d}^{\mathcal{E}}(X) \Lambda(d) \qquad (n \ge 2),$$

which enforces $b_n^{\mathcal{E}}(X) = 0$ for $2 \le n \le X$, $0 \le b_n^{\mathcal{E}}(X) \le 1$ for all $n \ge 1$ (0% overshoot), and $b_n^{\mathcal{E}}(X) = \frac{\Lambda(n)}{\ln n} \in [0,1]$ on the quadratic window $X < n < (\lfloor X\rfloor + 1)^2$.

## 5. Reproducible Verification

Run the standalone verification script to reproduce all $50$-digit and $550$-digit `mpmath` benchmarks, $1{,}002{,}000$-cell GPU continuum certifications, and publication figures:

```bash
python scripts/verify_paper4_debranges_lehmer_cascade.py
```
