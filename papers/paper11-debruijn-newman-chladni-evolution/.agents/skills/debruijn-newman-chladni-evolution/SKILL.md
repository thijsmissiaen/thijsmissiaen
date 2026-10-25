---
name: debruijn-newman-chladni-evolution
description: Scientific knowledge skill for Paper XI of the Transverse Wronskian series on Riemann's Xi-function. Contains the proven analytical identities and 60-decimal mpmath / 35-decimal mpmath.iv certifications for the De Bruijn-Newman initial zero-velocity identity d(gamma_k)/d(tau)|_{tau=0} = (1/4)*kappa(gamma_k), the dual extremum wall law d(tau_k)/d(tau)|_{tau=0} = (3/4)*kappa_Gamma(tau_k), the closed-form evaluation of the principal-value Calogero-Moser / Dyson interaction sum, the cross-term-free Gaussian Wigner factorization G_{y,tau}(x) = exp(tau*y^2/2)*exp(tau*x^2/2)*G_y(x) with pitchfork saddle-coalescence shift d(y_c)/d(tau)|_{tau=0} = eta_c / d_y L_{y_c}'(eta_c) > 0, and the Universal Transverse Dipole Evolution PDE Hierarchy d_tau c_{2m-1}(t;tau) + (1/8)*d_t^2 c_{2m-1}(t;tau) = (m(2m+1)/4)*c_{2m+1}(t;tau) governing Lehmer-pair backward coalescence thresholds.
---

# De Bruijn–Newman Heat-Flow Velocity $\left.\frac{\mathrm{d}\gamma_k}{\mathrm{d}\tau}\right| _ {\tau=0} \equiv \frac{1}{4}\kappa(\gamma_k)$, Cross-Term-Free Gaussian Wigner Factorization, and the Transverse Dipole Evolution PDE Hierarchy (Paper XI)

## 1. Scope and Mathematical Architecture

This skill encapsulates the exact analytical identities, partial differential hierarchies, directed-rounding `mpmath.iv` interval-arithmetic enclosures, and GPU continuum verification procedures of **Paper XI** in the series on the Transverse Wronskian of Riemann's $\Xi$-Function (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Paper XI establishes three exact structural bridges between the classical De Bruijn–Newman heat-flow theory (De Bruijn 1950, Newman 1976, Csordas–Norfolk–Varga 1988, Csordas–Smith–Varga 1994, Rodgers–Tao 2020, Polymath 15 2019) and the transverse Chladni, Wigner, and dipole-cascade geometry developed in Papers I–IV and X:

1. **Backwards Heat-Flow PDE, Parabolic Chladni Curvature $\kappa(\gamma_k)$, and Closed-Form Calogero–Moser / Dyson Sum (Section 2)**: For the heat-deformed Jacobi–Riemann kernel $\Phi _ \tau(u) := \Phi(u)\mathrm{e}^{\tau u^2}$ and $\Xi _ \tau(z) := 8\int_0^{\infty} \Phi(u)\mathrm{e}^{\tau u^2}\cos(2zu) \mathrm{d}u$ (where $\tau = 4t _ {\mathrm{DBN}}$ and $\tau _ {\mathrm{DBN}} = 4\Lambda$), double-exponential decay $\Phi(u) \ll \mathrm{e}^{9\lvert u\rvert}\exp(-\pi\mathrm{e}^{4\lvert u\rvert})$ justifies differentiating under the integral sign to all orders, yielding the backwards heat-flow equation $\partial _ \tau\Xi _ \tau(z) \equiv -\frac{1}{4}\partial_z^2\Xi _ \tau(z)$. Because both the 2D harmonic Laplace equation $\partial_d^2 R = -\partial_t^2 R$ on the critical line $d = 0$ and the 1D backwards heat equation $\partial _ \tau\Xi _ \tau = -\frac{1}{4}\partial_t^2\Xi _ \tau$ at $\tau = 0$ are governed by the same operator $-\partial_t^2$, the initial De Bruijn–Newman velocity of every simple zero $\gamma_k(\tau)$ is identically equal to one-fourth of the static parabolic Chladni curvature $\kappa(\gamma_k) := t _ {\mathcal{R} _ k}''(0) = \Xi''(\gamma_k)/\Xi'(\gamma_k)$ of the real zero curve $R(d, t _ {\mathcal{R} _ k}(d)) = 0$ entering the critical strip $d > 0$, while every critical-line extremum $\tau_k(\tau)$ ($\partial_t\Xi _ \tau(\tau_k(\tau)) = 0$) drifts with initial velocity $\frac{3}{4}\kappa _ \Gamma(\tau_k)$, where $\kappa _ \Gamma(\tau_k) = \frac{\Xi'''(\tau_k)}{3\Xi''(\tau_k)}$ is the curvature of the imaginary Chladni chamber wall $\Gamma_k : I(d, \tau _ {\Gamma_k}(d)) = 0$. Coupling this identity with the Hadamard partial-fraction expansion evaluates the conditionally convergent Calogero–Moser / Dyson interaction sum $\mathrm{P.V.}\sum _ {\gamma_j \neq \gamma_k} (\gamma_k - \gamma_j)^{-1}$ in exact closed form.
2. **Cross-Term-Free Gaussian Wigner Factorization and Pitchfork Saddle-Coalescence Shift (Section 3)**: In mirror phase-space coordinates $x = u + v$ and $y = u - v$ ($u = \frac{x+y}{2}, v = \frac{x-y}{2}$), the parallelogram law $u^2 + v^2 \equiv \frac{1}{2}(x^2 + y^2)$ contains zero $xy$ cross-term. Consequently, the deformed mirror kernel $G _ {y,\tau}(x) := \Phi _ \tau(\frac{x+y}{2})\Phi _ \tau(\frac{x-y}{2})$ factors into an outer hyperbolic shift weight $\mathrm{e}^{\frac{\tau}{2}y^2}$ and an inner holomorphic contour factor $\mathrm{e}^{\frac{\tau}{2}x^2}G_y(x)$. Because $\exp(-2\pi\cosh(2y)\mathrm{e}^{2L}\cos 2\eta)$ dominates $\mathrm{e}^{\frac{\tau}{2}L^2}$ for all $\tau \in \mathbb{R}$, the holomorphic strip $\lvert\mathrm{Im}(x)\rvert < \frac{\pi}{4}$ and Hermitian mirror symmetry $G _ {y,\tau}(-\xi + \mathrm{i}\eta) = \overline{G _ {y,\tau}(\xi + \mathrm{i}\eta)}$ remain invariant for all $\tau \in \mathbb{R}$, while the pitchfork saddle-point coalescence threshold $y_c(\omega; \tau)$ shifts with strictly positive initial velocity $\left.\frac{\partial y_c(\omega;\tau)}{\partial\tau}\right| _ {\tau=0} = \frac{\eta_c(\omega;0)}{\partial_y L _ {y_c}'(\eta_c(\omega;0))} > 0$.
3. **The Universal Transverse Dipole Evolution PDE Hierarchy and Lehmer Coalescence Threshold (Section 4)**: Differentiating the parity-reduced Wronskian $\widetilde{W} _ \tau(d,t) = \sum _ {m=1}^{\infty} c _ {2m-1}(t;\tau)d^{2m-2}$ under the backwards heat flow yields a fourth-order radial wave-heat PDE and proves that every transverse dipole moment $c _ {2m-1}(t;\tau)$ obeys the universal ladder PDE $\partial _ \tau c _ {2m-1}(t;\tau) + \frac{1}{8}\partial_t^2 c _ {2m-1}(t;\tau) \equiv \frac{m(2m+1)}{4}c _ {2m+1}(t;\tau)$ for all $m \ge 1$. At any $m$-fold coalescence point $t_0$ ($m \ge 2$), this ladder reduces via the universal square-cascade to $\partial _ \tau c _ {2m-3}(t_0; \tau _ {\mathrm{coal}}) \equiv \frac{m}{4}c _ {2m-1}(t_0; \tau _ {\mathrm{coal}}) = \frac{m^2}{4(m!)^2}(\partial_t^m\Xi _ {\tau _ {\mathrm{coal}}}(t_0))^2 > 0$, and across any closely spaced Lehmer pair $t _ {\pm} = t _ {\mathrm{mid}} \pm \delta$, it recovers the Csordas–Norfolk–Varga backward coalescence threshold in terms of the ratio of the first two transverse dipole invariants: $\tau _ {\mathrm{coal}}^{\mathrm{Lehmer}}(t _ {\mathrm{mid}}) = -2\frac{c_1(t _ {\mathrm{mid}};0)}{c_3(t _ {\mathrm{mid}};0)}(1 + \mathcal{O}(\delta^2)) = -2\delta^2(1 + \mathcal{O}(\delta^2)) < 0$.

## 2. Canonical Identities and Theorems

### 2.1 Backwards Heat-Flow PDE, Chladni Curvature Duality, and Closed-Form Dyson Sum (Theorem 2.1 & Corollary 2.2)

For $\Phi _ \tau(u) := \Phi(u)\mathrm{e}^{\tau u^2}$ and $\Xi _ \tau(z) := 8\int_0^{\infty} \Phi(u)\mathrm{e}^{\tau u^2}\cos(2zu) \mathrm{d}u$, differentiation under the integral sign yields

$$\frac{\partial \Xi _ \tau}{\partial \tau}(z) \equiv -\frac{1}{4}\frac{\partial^2 \Xi _ \tau}{\partial z^2}(z) \qquad (\forall (z,\tau) \in \mathbb{C} \times \mathbb{R}).$$

At any simple real zero $\gamma_k(\tau)$ of $\Xi _ \tau$ and any non-degenerate critical-line extremum $\tau_k(\tau)$ ($\partial_t\Xi _ \tau(\tau_k(\tau)) = 0$), implicit differentiation at $\tau = 0$ yields the exact Chladni curvature dualities

$$\left.\frac{\mathrm{d}\gamma_k}{\mathrm{d}\tau}\right| _ {\tau=0} = \frac{\Xi''(\gamma_k)}{4\Xi'(\gamma_k)} \equiv \frac{1}{4}\kappa(\gamma_k) \equiv \frac{1}{2}c _ {\mathrm{par}}(\gamma_k) \equiv -\frac{1}{4}\mathrm{Im}\left(2g_1\left(\frac{1}{2}+\mathrm{i}\gamma_k\right) + \frac{\zeta''(\frac{1}{2}+\mathrm{i}\gamma_k)}{\zeta'(\frac{1}{2}+\mathrm{i}\gamma_k)}\right),$$

$$\left.\frac{\mathrm{d}\tau_k}{\mathrm{d}\tau}\right| _ {\tau=0} = \frac{\Xi'''(\tau_k)}{4\Xi''(\tau_k)} \equiv \frac{3}{4}\kappa _ \Gamma(\tau_k), \qquad \kappa _ \Gamma(\tau_k) := \tau _ {\Gamma_k}''(0) = \frac{\Xi'''(\tau_k)}{3\Xi''(\tau_k)},$$

while the principal-value Calogero–Moser / Dyson interaction sum over all non-zero roots $\gamma _ {\pm m} = \pm\gamma_m$ ($m \ge 1$) excluding $+\gamma_k$ is evaluated in closed form by

$$\mathrm{P.V.}\sum _ {\gamma_j \neq \gamma_k} \frac{1}{\gamma_k - \gamma_j} := \frac{1}{2\gamma_k} + \sum _ {\substack{m=1 \\ m \neq k}}^{\infty} \frac{2\gamma_k}{\gamma_k^2 - \gamma_m^2} \equiv \frac{\Xi''(\gamma_k)}{2\Xi'(\gamma_k)} \equiv \frac{1}{2}\kappa(\gamma_k) \equiv c _ {\mathrm{par}}(\gamma_k).$$

### 2.2 Cross-Term-Free Gaussian Wigner Factorization and Pitchfork Shift (Theorem 3.1)

In mirror coordinates $x = u + v, y = u - v$, the identity $u^2 + v^2 \equiv \frac{1}{2}(x^2 + y^2)$ yields the exact factorization and Wigner correlation PDE ($\omega = 2t$):

$$G _ {y,\tau}(x) := \Phi _ \tau\left(\frac{x+y}{2}\right)\Phi _ \tau\left(\frac{x-y}{2}\right) \equiv \mathrm{e}^{\frac{\tau}{2}y^2}\mathrm{e}^{\frac{\tau}{2}x^2}G_y(x), \qquad \partial _ \tau I _ \tau(y,\omega) \equiv \frac{1}{2}\left(y^2 - \partial _ \omega^2\right)I _ \tau(y,\omega) \equiv \frac{1}{8}\left(4y^2 - \partial_t^2\right)I _ \tau(y,2t).$$

On the imaginary axis $x = \mathrm{i}\eta$ ($\lvert\eta\rvert < \frac{\pi}{4}$) with $L_y(\eta) := \ln G_y(\mathrm{i}\eta) = 2\ln\lvert\Phi(\frac{y+\mathrm{i}\eta}{2})\rvert$, the deformed pitchfork coalescence locus $(y_c(\omega;\tau), \eta_c(\omega;\tau))$ satisfies $L _ {y_c}'(\eta_c) - \tau\eta_c = \omega$ and $L _ {y_c}''(\eta_c) = \tau$, with exact initial velocities at $\tau = 0$:

$$\left.\frac{\partial y_c(\omega;\tau)}{\partial\tau}\right| _ {\tau=0} = \frac{\eta_c(\omega;0)}{\partial_y L _ {y_c}'\left(\eta_c(\omega;0)\right)} > 0, \qquad \left.\frac{\partial \eta_c(\omega;\tau)}{\partial\tau}\right| _ {\tau=0} = \frac{1 - \left(\partial_y L _ {y_c}''(\eta_c)\right)\partial _ \tau y_c(\omega;0)}{L _ {y_c}'''\left(\eta_c(\omega;0)\right)} < 0.$$

### 2.3 Universal Transverse Dipole Evolution PDE Hierarchy and Lehmer Coalescence (Theorem 4.1 & Corollary 4.2)

For all $(d,t,\tau) \in \mathbb{R}^3$, the parity-reduced Wronskian $\widetilde{W} _ \tau(d,t) = \sum _ {m=1}^{\infty} c _ {2m-1}(t;\tau)d^{2m-2}$ and derivative Wronskian $\widetilde{W} _ {1,2;\tau}(d,t) := \frac{1}{2d}\partial_d\lvert\Xi _ \tau'(t-\mathrm{i}d)\rvert^2$ satisfy the exact fourth-order PDEs

$$\frac{\partial \widetilde{W} _ \tau}{\partial\tau} + \frac{1}{8}\left(\frac{\partial^2\widetilde{W} _ \tau}{\partial t^2} - \frac{\partial^2\widetilde{W} _ \tau}{\partial d^2} - \frac{2}{d}\frac{\partial\widetilde{W} _ \tau}{\partial d}\right) \equiv 0, \qquad \frac{\partial \widetilde{W} _ \tau}{\partial\tau} + \frac{1}{4}\frac{\partial^2\widetilde{W} _ \tau}{\partial t^2} \equiv \frac{1}{8}\left(\frac{\partial^2\widetilde{W} _ \tau}{\partial t^2} + \frac{\partial^2\widetilde{W} _ \tau}{\partial d^2} + \frac{2}{d}\frac{\partial\widetilde{W} _ \tau}{\partial d}\right) \equiv \frac{1}{2}\widetilde{W} _ {1,2;\tau}(d,t),$$

and every transverse dipole invariant $c _ {2m-1}(t;\tau)$ obeys the **Universal Transverse Dipole Evolution PDE Hierarchy**

$$\frac{\partial c _ {2m-1}}{\partial\tau}(t;\tau) + \frac{1}{8}\frac{\partial^2 c _ {2m-1}}{\partial t^2}(t;\tau) \equiv \frac{m(2m+1)}{4}c _ {2m+1}(t;\tau) \qquad (\forall m \ge 1).$$

At any $m$-fold coalescence point $t_0$ ($m \ge 2$), $\partial _ \tau c _ {2m-3}(t_0; \tau _ {\mathrm{coal}}) \equiv \frac{m}{4}c _ {2m-1}(t_0; \tau _ {\mathrm{coal}}) = \frac{m^2}{4(m!)^2}(\partial_t^m\Xi _ {\tau _ {\mathrm{coal}}}(t_0))^2 > 0$, while across any closely spaced Lehmer pair $t _ {\pm} = t _ {\mathrm{mid}} \pm \delta$, $\left.\frac{\partial c_1}{\partial\tau}\right| _ {(t _ {\mathrm{mid}},0)} \equiv \frac{1}{4}(\Xi''(t _ {\mathrm{mid}}))^2 + \frac{1}{4}\Xi(t _ {\mathrm{mid}})\Xi^{(4)}(t _ {\mathrm{mid}}) > 0$, yielding the backward coalescence threshold $\tau _ {\mathrm{coal}}^{\mathrm{Lehmer}}(t _ {\mathrm{mid}}) = -2\frac{c_1(t _ {\mathrm{mid}};0)}{c_3(t _ {\mathrm{mid}};0)}(1 + \mathcal{O}(\delta^2)) = -2\delta^2(1 + \mathcal{O}(\delta^2)) < 0$.

## 3. Reproducible Verification Workflow

Run the standalone verification and figure-generation suite from the repository root:

```bash
python scripts/verify_paper11_debruijn_newman_chladni.py
```
