# Harmonic Nodal Chambers, Parabolic Zero-Curve Curvature, and the Transverse Integral Law for Riemann's $\Xi$-Function

**Author:** Thijs Missiaen (Independent Researcher, Belgium; ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810))  
**Repository & Version:** [`https://github.com/thijsmissiaen/paper3-harmonic-nodal-chambers`](https://github.com/thijsmissiaen/paper3-harmonic-nodal-chambers) (v1.0.0)  
**Manuscript:** [`paper/paper3_harmonic_nodal_chambers.pdf`](paper/paper3_harmonic_nodal_chambers.pdf) (16 pages, AMS-LaTeX / pdfTeX)  
**MSC 2020:** Primary 11M06, 30C15; Secondary 30D15, 31A05, 65G30  
**License:** MIT

## Abstract

Writing $\Xi(t - \mathrm{i}d) = \xi(\frac{1}{2} + d + \mathrm{i}t) = R(d,t) + \mathrm{i}I(d,t)$ as a harmonic potential ($\Delta R = \Delta I = 0$), we prove that the imaginary nodal walls $\Gamma_k := \lbrace (d,t) : I(d,t) = 0,  d > 0\rbrace$ departing orthogonally from $\Xi'(\tau_k) = 0$ partition $d > 0$ into oscillation-free harmonic nodal chambers $\Omega_k$ along whose boundaries $\lvert R\rvert$ grows strictly with arc length $s$ via $\frac{\mathrm{d}}{\mathrm{d}s}\lvert R(d(s),t(s))\rvert = \lvert\nabla I\rvert = \lvert\Xi'(t(s)-\mathrm{i}d(s))\rvert > 0$. Implicit Cauchy—Riemann differentiation shows that each real zero curve $R(d,t) = 0$ departs orthogonally from a simple zero $\gamma_k$ with parabolic curvature $\kappa(\gamma_k) = t''(0) = \Xi''(\gamma_k)/\Xi'(\gamma_k)$, while the transverse integral law $\lvert\Xi(t-\mathrm{i}d)\rvert^2 = \lvert\Xi(t)\rvert^2 + 2\int_0^d u \widetilde{W}(u,t) \mathrm{d}u$ transforms $\mathrm{sinc} _ {\mathrm{hyp}}(2u y)$ into the pure square $\mathrm{sinc} _ {\mathrm{hyp}}^2(d y) > 0$ via the hyperbolic low-pass filter $\mathcal{F} _ {\mathrm{int}}(2d y) = \frac{\tanh(d y)}{d y} \in (0,1)$ and forces a Rolle turning point $d _ {\mathrm{peak}} \in (0,d _ {\star})$ at any hypothetical off-critical zero. These identities are verified across $\gamma_1, \dots, \gamma _ {15}$, $\Omega_1, \dots, \Omega _ {25}$, and $t \in [14, 1000]$.

## Main Analytical Results

1. **Harmonic Nodal Chambers and Deductive Euler—Stirling Right Boundary (Propositions 2.1 and 2.2):** On the critical line $d = 0$, $I(0,t) \equiv 0$ (Dirichlet nodal wall), $\partial_d R(0,t) \equiv 0$ (Neumann mirror wall), and $\int _ {t_1}^{t_2} R(0,t) \partial_d R(0,t) \mathrm{d}t \equiv 0$. Transverse nodal walls $\Gamma_k \subset \lbrace I(d,t) = 0,  d > 0\rbrace$ depart orthogonally from $d = 0$ precisely at the extrema $\Xi'(\tau_k) = 0$ with curvature $\kappa _ \Gamma(\tau_k) = t _ {\Gamma_k}''(0) = \frac{\Xi'''(\tau_k)}{3\Xi''(\tau_k)}$. By the harmonic maximum principle ($\Delta I = 0$), no nodal branch can form a closed loop in $d \ge 0$, and each chamber $\Omega_k$ is closed on the right at the deductive Euler—Stirling boundary $d = d _ {\mathrm{Euler}}(t) = \frac{1}{2} + [\frac{1}{2}\ln\frac{t}{2\pi} + \gamma _ {\mathrm{E}} + \mathcal{O}((\ln t)^{-1})]^{-1} \to \frac{1}{2}^+$.
2. **Cauchy—Riemann Arc-Length Monotonicity and Parabolic Zero-Curve Curvature (Theorems 3.1 and 3.3):** Along any smooth wall $\Gamma_k$ parameterized by arc length $s$ moving into $d > 0$ from $(0, \tau_k)$, Cauchy—Riemann orthogonality ($\nabla R \perp \nabla I$) and $c_1(\tau_k) = -\Xi(\tau_k)\Xi''(\tau_k) > 0$ yield strict monotonicity:

$$\frac{\mathrm{d}}{\mathrm{d}s}\lvert R(d(s), t(s))\rvert = +\lvert\nabla I(d(s), t(s))\rvert = \lvert\Xi'(t(s) - \mathrm{i}d(s))\rvert > 0.$$

Around each simple real zero $\gamma_k$, the real zero curve $\mathcal{R} _ k := \lbrace (d, t _ {\mathcal{R} _ k}(d)) : R(d,t) = 0\rbrace$ departs orthogonally from $d = 0$ with exact parabolic curvature $\kappa(\gamma_k) := t _ {\mathcal{R} _ k}''(0)$ and parabolic Taylor coefficient $c _ {\mathrm{par}}(\gamma_k) = \frac{1}{2}\kappa(\gamma_k)$:

$$t _ {\mathcal{R} _ k}(d) - \gamma_k = c _ {\mathrm{par}}(\gamma_k) d^2 + c _ {\mathrm{quart}}(\gamma_k) d^4 + \mathcal{O}(d^6), \qquad \kappa(\gamma_k) = \frac{\Xi''(\gamma_k)}{\Xi'(\gamma_k)}, \qquad c _ {\mathrm{par}}(\gamma_k) = \frac{\Xi''(\gamma_k)}{2 \Xi'(\gamma_k)}.$$

3. **Transverse Integral Law, Rolle Reversal, and Hyperbolic Low-Pass Filter (Theorems 4.1 and 4.2):** Integrating $\partial_u\lvert\Xi(t-\mathrm{i}u)\rvert^2 = 2u \widetilde{W}(u,t)$ over $u \in [0,d]$ yields the Transverse Integral Law and complementary $90^\circ$ phase-locked pure squares:

$$\lvert\Xi(t - \mathrm{i}d)\rvert^2 = \lvert\Xi(t)\rvert^2 + 2\int_0^d u \widetilde{W}(u,t) \mathrm{d}u = \Big[\Xi(t)^2 + d^2\big(\Xi'(t)\big)^2\Big] - d^2\Xi(t)\Xi''(t) + 2\int_0^d u\big(\widetilde{W}(u,t) - \widetilde{W}(0,t)\big) \mathrm{d}u,$$

while the normalized integral operator $\frac{2}{d^2}\int_0^d u (\cdot) \mathrm{d}u$ transforms $\mathrm{sinc} _ {\mathrm{hyp}}(2u y)$ into the pure square $\mathrm{sinc} _ {\mathrm{hyp}}^2(d y) > 0$ via the strictly contractive hyperbolic low-pass filter

$$\mathcal{F} _ {\mathrm{int}}(2d y) := \frac{\mathrm{sinc} _ {\mathrm{hyp}}^2(d y)}{\mathrm{sinc} _ {\mathrm{hyp}}(2d y)} = \frac{\tanh(d y)}{d y} \in (0, 1) \qquad (\forall  d y > 0).$$

## Summary of Certified Quantitative Benchmarks

| Verification Tier | Domain / Scope | Key Certified Invariant / Metric | Verified Result |
| :--- | :--- | :--- | :--- |
| **Tier 1A: Harmonic Boundary** | $d = 0$, $t \in [14, 90]$ | Max residual of $\Delta R = \Delta I = 0$ and $I(0,t) = \partial_d R(0,t) = 0$ | $1.77 \times 10^{-50}$ (`mp.dps = 50`) |
| **Tier 1B: Parabolic Curvature** | $\gamma_1, \dots, \gamma _ {15}$ | First zero curvature $\kappa(\gamma_1) = \Xi''(\gamma_1)/\Xi'(\gamma_1)$ and $c _ {\mathrm{par}}(\gamma_1) = \frac{1}{2}\kappa(\gamma_1)$ | $\kappa(\gamma_1) = -1.1592611550$, $c _ {\mathrm{par}}(\gamma_1) = -0.5796305775$ |
| **Tier 2: Chamber Tracking** | $\Omega_1, \dots, \Omega _ {25}$, $d \in (0, 1.65]$ | Wall growth $R(d _ {\mathrm{Euler}})/\Xi(\tau_k)$ and interior sign purity ($103{,}750$ points) | Growth $1.81\times$ to $6.20\times$; 100.0000% sign-pure |
| **Tier 3A: Integral Law** | $t \in [\gamma_1, 1000]$, $d \in [0.25, 0.50]$ | Max relative residual of boxed Transverse Integral Law | $8.75 \times 10^{-51}$ (`mp.dps = 50`) |
| **Tier 3B: Rolle Reversal** | $(d _ {\star}, t _ {\star}) = (0.25, 14.20)$ | Mandatory Rolle turning point $d _ {\mathrm{peak}} \in (0, d _ {\star})$ and integral residual | $d _ {\mathrm{peak}} = 0.1342212372$; residual $5.22 \times 10^{-54}$ |
| **Tier 3C: Continuum Scan** | $t \in [14, 1000]$, $60{,}000$ points | Mean amplification $\widetilde{\mathcal{E}}(d,t)/\widetilde{W}(d,t)$ ($d = 0.05$ to $0.35$) and $\min \mathcal{R} _ {\mathrm{sq}}$ | $77.19\times$ ($d=0.05$), $2.11\times$ ($d=0.35$); $\min \mathcal{R} _ {\mathrm{sq}} = 0.6544 > 0$ |

## Repository Structure

```text
paper3-harmonic-nodal-chambers/
├── paper/
│   ├── paper3_harmonic_nodal_chambers.tex
│   ├── paper3_harmonic_nodal_chambers.pdf
│   └── figures/
│       ├── fig1_harmonic_chladni_chambers_and_curvature.pdf
│       ├── fig1_harmonic_chladni_chambers_and_curvature.png
│       ├── fig2_transverse_integral_law_and_lowpass.pdf
│       └── fig2_transverse_integral_law_and_lowpass.png
├── scripts/
│   └── verify_paper3_harmonic_nodal_chambers.py
├── .agents/
│   └── skills/
│       └── harmonic-nodal-chambers-xi/
│           └── SKILL.md
├── AGENT_SKILL_MANIFEST.json
├── CITATION.cff
├── .zenodo.json
├── LICENSE
└── README.md
```

## Reproducing All Verifications and Figures

Requires Python 3.10+ with `mpmath`, `numpy`, `scipy`, `matplotlib`, and `torch` (optional CUDA acceleration). Run from the repository root:

```bash
python scripts/verify_paper3_harmonic_nodal_chambers.py
```
