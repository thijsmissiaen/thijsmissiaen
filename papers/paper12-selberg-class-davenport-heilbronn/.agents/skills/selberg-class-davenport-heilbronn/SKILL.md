---
name: selberg-class-davenport-heilbronn
description: Scientific knowledge skill for Paper XII of the Transverse Wronskian series on Riemann's Xi-function. Contains the proven analytical identities and 60-decimal mpmath / 35-decimal mpmath.iv certifications for Universal Mother-Wavelet Dilation and Symmetric Branch Damping in the Selberg Class, the Phase-Independent Confocal Ellipse Identity, the Multiplicativity Defect +-(1+kappa^2) and Breakdown of Prime-Gauge Purity at N=4,6 in the Davenport-Heilbronn counterexample, and the Four-Fold Diagnostic Detection of the first off-line zero (d_*, gamma_*).
---

# Universal Mother-Wavelet Dilation and the Phase-Independent Confocal Ellipse Identity in the Selberg Class versus Breakdown of Prime-Gauge Purity and Chladni Chamber Topology in the Davenport–Heilbronn Counterexample (Paper XII)

## 1. Scope and Mathematical Architecture

This skill encapsulates the analytical theorems, exact mode-by-mode algebraic decompositions, directed-rounding `mpmath.iv` interval-arithmetic enclosures, and GPU continuum verification procedures of **Paper XII** in the series on the Transverse Wronskian of Riemann's Ξ-Function (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Paper XII resolves the classical **Davenport–Heilbronn litmus test** (Titchmarsh 1986, Section 10.25; Davenport & Heilbronn 1936; Bombieri & Ghosh 2011; Balanzario & Sánchez-Ortiz 2007) by establishing a sharp structural dichotomy between universal harmonic laws that hold for every self-dual Dirichlet series with a Gamma functional equation and the Euler-product closure mechanisms that fail for the non-multiplicative Davenport–Heilbronn counterexample $\mathcal{D} _ {\mathrm{DH}}(s)$:

1. **Universal Mother-Wavelet Dilation and Symmetric Branch Damping in the Selberg Class (Section 2)**: For any self-dual Dirichlet series $L(s,f) = \sum _ {n=1}^{\infty} a_f(n)n^{-s}$ ($a_f(n) \in \mathbb{R}$, $a_f(1)=1$) with completed function $\xi_f(s) = Q_f^s\gamma_f(s)L(s,f) = \xi_f(1-s)$ of degree $d_f > 0$ (including Riemann's $\zeta(s)$, odd primitive Dirichlet $L(s,\chi_5)$, Ramanujan's modular cusp form $L(s,\Delta)$ of degree $d _ {\Delta} = 2$, and $\mathcal{D} _ {\mathrm{DH}}(s){}$), the logarithmic kernel $\Phi_f(u) = \Phi_f(-u)$ decomposes into dilated, amplitude-weighted copies of a single class mother wavelet $\varphi_1^{(f)}(w)$ via $\phi_m^{(f)}(u) \equiv a_f(m)m^{-1/2}\varphi_1^{(f)}(u + \frac{1}{2}\ln m)$, while strict trigamma convexity enforces the Symmetric Branch-Damping Inequality $\rho _ {\Gamma,f}(d,t)\lvert a_f(n)\rvert n^{-1/2+d} < \mathrm{e}^{-2d\omega _ {n,Q_f}(t)}\lvert a_f(n)\rvert n^{-1/2-d} \le \lvert a_f(n)\rvert n^{-1/2-d}$ across all active modes $1 \le n \le N_f(t) = \lfloor\exp(\vartheta_f'(t))\rfloor$.
2. **The Universal Phase-Independent Confocal Ellipse Identity (Section 3)**: Regardless of multiplicativity or an Euler product, every combined forward and reflected Riemann–Siegel mode $T_n^{(f)}(d,t) = 2a_n n^{-1/2}(\cosh(d\omega _ {n,Q})\cos\phi_n + \mathrm{i}\sinh(d\omega _ {n,Q})\sin\phi_n)$ satisfies the exact phase-independent increment law $\lvert T_n^{(f)}(d,t)\rvert^2 - \lvert T_n^{(f)}(0,t)\rvert^2 \equiv 4\frac{\lvert a_n\rvert^2}{n}\sinh^2(d\omega _ {n,Q}(t)) \ge 0$ for all $\phi_n(t) \in \mathbb{R}$, tracing a family of confocal ellipses with $d$-invariant focal distance $c_n^{(f)} \equiv 2\lvert a_n\rvert n^{-1/2}$, symplectic area $2\pi\frac{\lvert a_n\rvert^2}{n}\sinh(2d\lvert\omega _ {n,Q}\rvert)$, and diagonal ratio bounded by $\lvert a_n\rvert^2 n^{-(1+2d)}$.
3. **Multiplicativity Defect $\pm(1+\kappa^2)$ and Breakdown of Prime-Gauge Purity at $N = 4, 6$ in $\mathcal{D} _ {\mathrm{DH}}(s)$ (Section 4)**: Whereas any Euler-product $L$-function satisfies Quadratic Prime-Power Purity $b_N^{\mathcal{E},(f)}(X) = \frac{\Lambda_f(N)}{\ln N}$ on $(X, X^2]$ (annihilating all non-prime-power composites), the Davenport–Heilbronn function $\mathcal{D} _ {\mathrm{DH}}(s)$ ($\kappa = \frac{\sqrt{10-2\sqrt{5}}-2}{\sqrt{5}-1} \approx 0.2840790438$) suffers an exact algebraic multiplicativity defect $\delta _ {\mathrm{DH}}(m,n) = a _ {\mathrm{DH}}(mn) - a _ {\mathrm{DH}}(m)a _ {\mathrm{DH}}(n) = \mp(1+\kappa^2) = \mp 1.0807009031\dots \neq 0$ at $N = 4 = 2^2$, $N = 6 = 2\cdot 3$, $N = 9 = 3^2$, and $N = 14 = 2\cdot 7$, leaving unsieved composite residuals $b_4^{\mathcal{E}}(X) = -(1+\kappa^2)$ and $b_6^{\mathcal{E}}(X) = +(1+\kappa^2)$. At the first off-line zero $(d _ {\star}, \gamma _ {\star}) = (0.3085171825, 85.6993484854)$ ($N_5(\gamma _ {\star}) = 8$), all active primes $p \in \lbrace 2,3,5,7\rbrace$ combined contribute only 38.82% of $\lvert T_1^{(\mathrm{DH})}\rvert$ (incapable of creating a zero), whereas the non-multiplicative composites $n \in \lbrace 4,6,8\rbrace$ (enhanced by $1/\kappa^2 = 12.3914\times$) supply the remaining 54.03% destructive interference.
4. **Four-Fold Diagnostic Detection of the Off-Line Zero $(d _ {\star}, \gamma _ {\star})$ (Section 5)**: At $(d _ {\star}, \gamma _ {\star})$, all four diagnostic instruments register the topological and spectral breakdown: (i) saddle recombination creates a false critical point $z_c = \tau _ {\mathrm{false}} = 85.7472008265$ where $\Xi _ {\mathrm{DH}}(z_c) < 0$ and $\Xi _ {\mathrm{DH}}''(z_c) < 0$ share the same negative sign, reversing Cauchy–Riemann wall monotonicity $\frac{\mathrm{d}\lvert R _ {\mathrm{DH}}\rvert}{\mathrm{d}s}\rvert _ {\Gamma _ {\mathrm{false}}} = -\lvert\nabla I _ {\mathrm{DH}}\rvert < 0$ until $\Gamma _ {\mathrm{false}}$ collides with the detached real nodal parabola $R _ {\mathrm{DH}} = 0$ ($d _ {\mathrm{tip}} = 0.2980980120$); (ii) $\widetilde{W} _ {\mathrm{DH}}(d,\gamma _ {\star}) < 0$ on $[0, d _ {\star})$ with negative Transverse Integral Law; (iii) the Hadamard flux jumps by $\Phi _ {\mathrm{Had}}(d; d _ {\star}) = \pi\Theta(d - d _ {\star})$; and (iv) the de Branges inner function crosses $\lvert\Theta _ {\alpha}(\gamma _ {\star} + \mathrm{i}d _ {\star})\rvert = 1$ with Sylvester inertia jump $\nu _ - : 0 \to 1$.

## 2. Canonical Identities and Theorems

### 2.1 Universal Mother-Wavelet Dilation and Symmetric Branch Damping (Theorem 2.1)

For any self-dual completed $L$-function $\xi_f(s) = \mathcal{G} _ f(s)L(s,f) = \xi_f(1-s)$ with Archimedean factor $\mathcal{G} _ f(s) = Q_f^s\gamma_f(s)$ of degree $d_f > 0$, define the class mother wavelet $\varphi_1^{(f)}(w)$ by inverse bilateral Laplace–Mellin transformation $\mathcal{G} _ f(s) = 4\int _ {-\infty}^{\infty} \varphi_1^{(f)}(w)\mathrm{e}^{2(s-1/2)w}\mathrm{d}w$. Explicitly, for Riemann's $\zeta(s)$ ($d _ {\zeta}=1$), the odd modulo $q$ class ($d_q=1$, including $L(s,\chi_5)$ and $\mathcal{D} _ {\mathrm{DH}}(s){}$), and Ramanujan's modular cusp form $L(s,\Delta)$ ($d _ {\Delta}=2$), the mother wavelets are

$$\varphi_1^{(\zeta)}(w) = \left(2\pi^2\mathrm{e}^{9w} - 3\pi\mathrm{e}^{5w}\right)\mathrm{e}^{-\pi\mathrm{e}^{4w}}, \qquad \varphi_1^{(q)}(w) = \mathrm{e}^{3w}\mathrm{e}^{-\frac{\pi}{q}\mathrm{e}^{4w}}, \qquad \varphi_1^{(\Delta)}(w) = \frac{1}{2}\mathrm{e}^{12w}\mathrm{e}^{-2\pi\mathrm{e}^{2w}}.$$

Then the full kernel $\Phi_f(u) = \Phi_f(-u)$ in $\Xi_f(z) = 4\int _ {-\infty}^{\infty} \Phi_f(u)\mathrm{e}^{2\mathrm{i}zu}\mathrm{d}u$ decomposes into dilated shells

$$\Phi_f(u) = \sum _ {m=1}^{\infty} \phi_m^{(f)}(u), \qquad \phi_m^{(f)}(u) \equiv \frac{a_f(m)}{\sqrt{m}}\varphi_1^{(f)}\left(u + \frac{1}{2}\ln m\right) \qquad (\forall m \ge 1,   u \in \mathbb{C}),$$

and for every active mode $1 \le n \le N_f(t) = \lfloor\exp(\vartheta_f'(t))\rfloor$ ($\omega _ {n,Q_f}(t) := \vartheta_f'(t) - \ln n \ge 0$), strict trigamma curvature enforces

$$\rho _ {\Gamma,f}(d,t)\lvert a_f(n)\rvert n^{-1/2+d} < \mathrm{e}^{-2d\omega _ {n,Q_f}(t)}\lvert a_f(n)\rvert n^{-1/2-d} \le \lvert a_f(n)\rvert n^{-1/2-d}.$$

### 2.2 The Universal Phase-Independent Confocal Ellipse Identity (Theorem 3.1)

For every combined forward and reflected Riemann–Siegel mode $T_n^{(f)}(d,t) := 2a_n n^{-1/2}(\cosh(d\omega _ {n,Q}(t))\cos\phi_n(t) + \mathrm{i}\sinh(d\omega _ {n,Q}(t))\sin\phi_n(t))$ with $a_n \in \mathbb{R}$ and $\phi_n(t) \in \mathbb{R}$, the identity $\cosh^2(d\omega _ {n,Q}) \equiv 1 + \sinh^2(d\omega _ {n,Q})$ combined with $\cos^2\phi_n + \sin^2\phi_n \equiv 1$ yields the exact phase-independent increment law

$$\left\lvert T_n^{(f)}(d,t)\right\rvert^2 - \left\lvert T_n^{(f)}(0,t)\right\rvert^2 \equiv 4\frac{\lvert a_n\rvert^2}{n}\sinh^2\left(d\omega _ {n,Q}(t)\right) \ge 0 \qquad (\forall \phi_n(t) \in \mathbb{R}),$$

with $d$-invariant focal distance $c_n^{(f)} = \sqrt{A_n(d,t)^2 - B_n(d,t)^2} \equiv 2\lvert a_n\rvert n^{-1/2}$, symplectic area $\mathrm{Area}(\mathcal{E} _ n^{(f)}(d,t)) \equiv 2\pi\frac{\lvert a_n\rvert^2}{n}\sinh(2d\lvert\omega _ {n,Q}(t)\rvert)$, and relative increment ratio bounded by $\lvert a_n\rvert^2 n^{-(1+2d)}$.

### 2.3 Multiplicativity Defect $\pm(1+\kappa^2)$ and Prime-Gauge Breakdown (Theorem 4.1)

For the Davenport–Heilbronn function $\mathcal{D} _ {\mathrm{DH}}(s) = \frac{1-\mathrm{i}\kappa}{2}L(s,\chi_5) + \frac{1+\mathrm{i}\kappa}{2}L(s,\overline{\chi_5})$ with $\kappa = \frac{\sqrt{10-2\sqrt{5}}-2}{\sqrt{5}-1} \approx 0.2840790438$, the 5-periodic coefficients $a _ {\mathrm{DH}}(1)=1, a _ {\mathrm{DH}}(2)=\kappa, a _ {\mathrm{DH}}(3)=-\kappa, a _ {\mathrm{DH}}(4)=-1, a _ {\mathrm{DH}}(5)=0$ exhibit the exact multiplicativity defect

$$\delta _ {\mathrm{DH}}(2,2) = -(1 + \kappa^2), \qquad \delta _ {\mathrm{DH}}(2,3) = +(1 + \kappa^2), \qquad \delta _ {\mathrm{DH}}(3,3) = -(1 + \kappa^2), \qquad 1 + \kappa^2 = 1.08070090314928289127\dots \neq 0,$$

leaving unsieved composite residuals $b_4^{\mathcal{E}}(X) = -(1+\kappa^2)$ and $b_6^{\mathcal{E}}(X) = +(1+\kappa^2)$ under any prime-generated gauge factor $E_X^{\mathcal{E}}(s)$ ($X \ge 5$).

### 2.4 Four-Fold Diagnostic Detection of $(d _ {\star}, \gamma _ {\star})$ (Theorem 5.1)

At the first off-line zero $(d _ {\star}, \gamma _ {\star}) = (0.3085171824566374, 85.69934848537759)$, $\Xi _ {\mathrm{DH}}(z)$ exhibits: (i) a critical-line saddle bifurcation at $z_c = \tau _ {\mathrm{false}} = 85.7472008265$ where $\Xi _ {\mathrm{DH}}(z_c)\Xi _ {\mathrm{DH}}''(z_c) > 0$, reversing wall monotonicity $\frac{\mathrm{d}\lvert R _ {\mathrm{DH}}\rvert}{\mathrm{d}s}\rvert _ {\Gamma _ {\mathrm{false}}} = -\lvert\nabla I _ {\mathrm{DH}}\rvert < 0$ until intersection with the detached parabola $R _ {\mathrm{DH}} = 0$ ($d _ {\mathrm{tip}} = 0.2980980120$); (ii) $\widetilde{W} _ {\mathrm{DH}}(d,\gamma _ {\star}) < 0$ on $[0, d _ {\star})$ with $\int_0^{d _ {\star}} u\widetilde{W} _ {\mathrm{DH}}(u,\gamma _ {\star})\mathrm{d}u = -\frac{1}{2}\lvert\Xi _ {\mathrm{DH}}(\gamma _ {\star})\rvert^2 < 0$; (iii) the quantized Hadamard flux jump $\Phi _ {\mathrm{Had}}(d; d _ {\star}) = \int _ {-\infty}^{\infty} K(d, d _ {\star}, u)\mathrm{d}u = \pi\Theta(d - d _ {\star})$ via exact arcsine cancellation $\mp\arcsin(d/d _ {\star})$; and (iv) the de Branges inner-function crossing $\lvert\Theta _ {\alpha}(\gamma _ {\star} + \mathrm{i}d _ {\star})\rvert = 1$ with Sylvester negative inertia jump $\nu _ - : 0 \to 1$.

## 3. Reproducible Verification Workflow

Run the standalone verification and figure-generation suite from the repository root:

```bash
python scripts/verify_paper12_selberg_davenport_heilbronn.py
```
