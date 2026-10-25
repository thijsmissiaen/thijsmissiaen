# The Chladni Portal-Frame Synthesis: Codimension-Three Saddle Overdetermination, Green–Dirichlet–Neumann Chamber Equipartition, and Euler Dilogarithm $\mathrm{Li} _ 2(2\theta)$ Horizon Damping for Riemann's $\Xi$ Function

**Author:** Thijs Missiaen ([thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com))  
**Affiliation:** Independent Researcher, Belgium  
**ORCID:** [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)  
**Date:** October 25, 2026  
**Manuscript:** [`paper/paper14_portal_dilogarithm_synthesis.pdf`](paper/paper14_portal_dilogarithm_synthesis.pdf) (22 pages, AMS-LaTeX / Tectonic `pdfTeX`)

> [!IMPORTANT]
> **Resolved Classical Bottleneck / Open Problem Context (Three-Scale Portal-Frame Dichotomy, Chamber Flux Equipartition, and Semiprime Horizon Damping)**
>
> Analyzing the zeros of Riemann's completed zeta function $\Xi(t - \mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$ inside a harmonic Chladni nodal chamber $\Omega_k(d_0)$ links three structural questions across macro-, meso-, and micro-scales: how off-critical zeros can arise along a Chladni wall $I(d,t)=0$, how boundary flux is partitioned between $R$ and $I$ across curved chamber walls, and how the Riemann–Siegel cutoff $pq \le N _ \star(t) = \sqrt{t/(2\pi)}$ modifies the quadratic prime-gauge expansion at the cubic window $\theta = 1/6$.
>
> **What this paper proves:**
> 1. **Macro-Scale Dichotomy and Codimension-Three Overdetermination (Theorems 2.1–2.2):** Along any Chladni wall $I=0$, $\mathrm{d}R/\mathrm{d}s = \pm \lvert \Xi'\rvert$ and $\frac{1}{2}\partial_d^2[R^2] _ {d=0} = -\Xi(\tau _ \star)\Xi''(\tau _ \star) = c_1(\tau _ \star)$, forcing any off-critical zero into either **Case I** ($c_1(\tau _ \star) \le 0$ on $d=0$, satisfied by 100% of the first five off-critical zeros of Davenport–Heilbronn $\Xi _ {\mathrm{DH}}$) or **Case II** (an interior saddle satisfying the codimension-three system $I = \partial_d I = \partial_t I = 0$).
> 2. **Meso-Scale 50%–50% Green–Dirichlet–Neumann Equipartition and Wronskian–Slope Reciprocity (Theorems 3.1–3.2):** Homogeneous Dirichlet ($I=0$) and Neumann ($\partial_n R = 0$) wall conditions enforce $\int R \partial_d R \mathrm{d}t \equiv \int I \partial_d I \mathrm{d}t \equiv \iint _ {\Omega_k} \lvert \Xi'\rvert^2 \mathrm{d}d \mathrm{d}t = \frac{1}{2}\mathcal{W} _ {\Omega_k}(d_0) > 0$, together with $\tau_k'(d) \equiv R \partial_t R/(d \widetilde{W})$.
> 3. **Micro-Scale Euler Dilogarithm $\mathrm{Li} _ 2(2\theta)$ Capture Laws and Hyperbola Duality (Theorems 4.1–4.3):** Quadratic Wronskian horizon damping ($p=2$) raises the forward-triangle capture fraction at $\theta = 1/6$ from $\eta^{(0)}(1/6) = 0.243955 < 1/2$ to $\eta^{(2)}(1/6) = 0.676891 > 2/3 > 1/2$ ($\Delta = -0.353781 < 0$, $\theta _ {*,50\text{ pct}}^{(2)} = 0.190550$), while the reflected corner obeys fourth-order pinching $(\omega_p\omega_q)^2 \le \frac{1}{16}(\ln m _ {\mathrm{dual}})^4$.

## Abstract

We unify the macro-, meso-, and micro-scale Cauchy–Riemann geometry of Riemann's completed zeta function $\Xi(t-\mathrm{i}d)=R(d,t)+\mathrm{i}I(d,t)$ ($d\ge 0$) around a harmonic Chladni nodal chamber $\Omega_k(d_0)$. At the macro-scale, we prove that any off-critical zero on a transverse Chladni wall $I(d,t)=0$ is governed by an exact dichotomy between a critical-line Laguerre defect $c_1(\tau _ \star)=-\Xi(\tau _ \star)\Xi''(\tau _ \star)\le 0$ and a codimension-three interior saddle overdetermination $\lbrace I=\partial_d I=\partial_t I=0\rbrace$ ($3$ real equations in $2$ unknowns), with 100% of the first five off-critical zeros of the Davenport–Heilbronn counterexample $\Xi _ {\mathrm{DH}}(z)$ forced into the critical-line defect case. At the meso-scale, Green's first identity enforces an exact 50%–50% Dirichlet–Neumann boundary-flux equipartition

$$\int _ {\tau_k(d_0)}^{\tau _ {k+1}(d_0)} R \partial_d R \mathrm{d}t \equiv \int _ {\tau_k(d_0)}^{\tau _ {k+1}(d_0)} I \partial_d I \mathrm{d}t \equiv \iint _ {\Omega_k(d_0)} \lvert \Xi'(t-\mathrm{i}d)\rvert^2 \mathrm{d}d \mathrm{d}t = \frac{1}{2}\mathcal{W} _ {\Omega_k}(d_0) > 0$$

coupled with the Wronskian–slope reciprocity $\tau_k'(d)\equiv R \partial_t R/(d \widetilde{W})$. At the micro-scale, Euler dilogarithm $\mathrm{Li} _ 2(2\theta)$ horizon damping raises the quadratic Wronskian forward-triangle capture fraction at $\theta=1/6$ from $\eta^{(0)}(1/6)=0.243955$ to $\eta^{(2)}(1/6)=0.676891>2/3>1/2$ ($\Delta(\eta^{(2)}(1/6))=-0.353781<0$) with fourth-order Dirichlet hyperbola apex pinching $(\omega_p\omega_q)^2\le \frac{1}{16}(\ln m _ {\mathrm{dual}})^4$.

---

## Key Analytical Theorems & Certified Constants

| Scale | Theorem | Core Analytical Identity / Result | Certified Benchmark Value (65-Decimal `mpmath` / `iv`) |
| :--- | :--- | :--- | :--- |
| **Macro (Sec. 2)** | **Theorem 2.1** | Harmonic curvature conservation and codimension-3 dichotomy | `0.5 * d_d^2[R^2]_{d=0} = -0.5 * d_t^2[Xi^2]_{t=tau*} = c_1(tau*)` |
| **Macro (Sec. 2)** | **Theorem 2.2** | 100% Case-I forcing of `Xi_DH(z)` off-line zeros | `S_m^2 * c_1^DH(tau_false,m) in [-873.26, -34.21] < 0` (`m = 1..5`) |
| **Meso (Sec. 3)** | **Theorem 3.1** | Boxed 50%–50% Green–Dirichlet–Neumann equipartition | `R_Omega(d_0) == I_Omega(d_0) == E_Omega(d_0) = 0.5 * W_Omega(d_0) > 0` (`< 1e-57` rel. err.) |
| **Meso (Sec. 3)** | **Theorem 3.2** | Wronskian–slope reciprocity and fold-tip `W_tilde = 0` | `(d_tip, t_tip)_DH = (0.2980980120, 85.6254944428)`, `abs(W_tilde_DH) = 2.47e-144` |
| **Micro (Sec. 4)** | **Theorem 4.1** | Attenuated spring floor for `eta > 1/2` (`Delta = 1 - 2*eta < 0`) | `abs(1 + U + 0.5*eta*U^2)^2 >= (1 - 1/(2*eta))^2 = +0.0682917228 > 0` at `eta^(2)(1/6)` |
| **Micro (Sec. 4)** | **Theorem 4.2** | Exact closed-form `Li_2(2*theta)` capture fractions `eta^(p)(1/6)` | `eta^(0) = 0.24395510`, `eta^(1) = 0.49559874`, `eta^(2) = 0.67689062`, `eta^(3) = 0.79809062` |
| **Micro (Sec. 4)** | **Theorem 4.2(iii)** | Quadratic Wronskian 50% stability horizon `theta_{*,50%}^(2)` | `theta_{*,50%}^(2) = 0.190549846717 > theta_{W-pure}^inf = 0.183944870826 > 1/(2e) > 1/6` |
| **Micro (Sec. 4)** | **Theorem 4.3** | Hyperbola reflection duality and 4th-order apex pinching | `omega_p + omega_q == ln(m_dual(pq))`, `(omega_p * omega_q)^2 <= (1/16) * (ln m_dual)^4 -> 0` |

---

## Reproducible Verification & Compilation

### 1. Run the 65-Decimal `mpmath` & Interval-Arithmetic Verification Suite
```bash
python scripts/verify_paper14_portal_dilogarithm_synthesis.py
```
This script executes all analytical assertions for Sections 2, 3, and 4 at $65$-decimal working precision (`mp.dps = 65`) with $35$-decimal rigorous interval certification (`mp.iv.dps = 35`), regenerates all four publication figure files in `paper/figures/`, and exits with code `0`.

### 2. Compile the AMS-LaTeX Manuscript
```bash
tectonic paper/paper14_portal_dilogarithm_synthesis.tex
```

## License & Citation

Released under the MIT License. Please cite this work using the metadata in [`CITATION.cff`](CITATION.cff).
