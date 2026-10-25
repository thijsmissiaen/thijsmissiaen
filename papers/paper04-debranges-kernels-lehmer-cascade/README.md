# De Branges Reproducing Kernels, the Universal Transverse Square-Cascade for Lehmer Pairs, and the Leibniz–Wronskian Identity for Riemann's $\Xi$-Function

**DOI (Zenodo / CERN):** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23641702.svg)](https://doi.org/10.5281/zenodo.23641702) [`10.5281/zenodo.23641702`](https://doi.org/10.5281/zenodo.23641702)  
**Author:** Thijs Missiaen (Independent Researcher, Belgium; [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com); ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810))  
**Series:** Paper IV of the Transverse Wronskian Series on Riemann's $\Xi$-Function  
**Manuscript:** [`paper/paper4_debranges_lehmer_cascade.pdf`](paper/paper4_debranges_lehmer_cascade.pdf) (18 pages, AMS-LaTeX / pdfTeX)

> [!IMPORTANT]
> **Resolved Classical Bottleneck / Open Problem Context (Conrey–Li De Branges Obstruction & Lehmer Near-Coalescence Singularity):** Conrey and Li (2000) proved that global positivity conditions on Sonine–De Branges Hilbert spaces fail for $\Xi(z)$, while Lehmer pairs ($t _ \pm = t _ {\mathrm{mid}} \pm \delta$ as $\delta \to 0$) drive the leading critical-line invariant $c_1(t _ {\mathrm{mid}}) = \mathcal{O}(\delta^2)$ toward zero. This paper resolves both bottlenecks by proving the exact local identity $K _ \alpha(t-\mathrm{i}d, t-\mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W}(d,t)$ (decoupling local off-critical zero exclusion from global Sonine obstructions) and the Universal Transverse Square-Cascade $c _ {2m-1}(t_0) = \frac{m}{(m!)^2}(\Xi^{(m)}(t_0))^2 > 0$, which guarantees an unbroken 2D quadratic bowl $\widetilde{W}(d,t _ {\mathrm{mid}}+\tau) = \frac{1}{2}(\Xi''(t _ {\mathrm{mid}}))^2(\delta^2+\tau^2+d^2)+\mathcal{O}(3) > 0$.

## Overview

This repository contains the complete AMS-LaTeX manuscript, compiled PDF, reproducible verification suite, and AI Agent Skill for **Paper IV** of the research series on the parity-reduced transverse Wronskian $\widetilde{W}(d,t) = \frac{1}{2d}\partial_d\left|\Xi(t-\mathrm{i}d)\right|^2$ of Riemann's completed zeta function $\Xi(z) = \xi(\frac{1}{2}+\mathrm{i}z)$.

The paper establishes three exact analytical identities that resolve how the transverse Wronskian interacts with Hilbert spaces of entire functions, nearly coalescing critical-line zeros (Lehmer pairs), and multiplicative prime-sieve gauges:

1. **The Exact De Branges Reproducing Kernel Identity (Section 2, Theorem 2.2 & Proposition 2.4):** For every coupling parameter $\alpha > 0$, the linear Hermite–Biehler family $E _ \alpha(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z)$ and its Schwarz reflection $E _ \alpha^{\sharp}(z) := \Xi(z) + \mathrm{i}\alpha\Xi'(z)$ satisfy the exact global identity

$$K _ \alpha(t - \mathrm{i}d, t - \mathrm{i}d) := \frac{\left|E _ \alpha(t - \mathrm{i}d)\right|^2 - \left|E _ \alpha^{\sharp}(t - \mathrm{i}d)\right|^2}{4\pi d} = \frac{\left|E _ \alpha^{\sharp}(t + \mathrm{i}d)\right|^2 - \left|E _ \alpha(t + \mathrm{i}d)\right|^2}{4\pi d} \equiv \frac{\alpha}{\pi} \widetilde{W}(d,t).$$

Both $\left|\Xi\right|^2$ and $\alpha^2\left|\Xi'\right|^2$ cancel identically via the Cauchy–Riemann equations. Furthermore, any finite $M \times M$ Cauchy–Pick Gram matrix $\mathbf{G} _ \alpha(z_1,\dots,z_M) _ {j,k} = K _ \alpha(z_j, z_k)$ reduces to the two-point Wronskian divided difference

$$\left[\mathbf{G} _ \alpha(z_1,\dots,z_M)\right] _ {j,k} = \frac{\alpha}{\pi} \frac{\Xi(z_j)\overline{\Xi'(z_k)} - \Xi'(z_j)\overline{\Xi(z_k)}}{z_j - \bar{z} _ k}.$$

This decouples local off-critical zero exclusion ($\widetilde{W}(d_0, t_0) > 0$) from the global Sonine–De Branges inner-product positivity obstructions identified by Conrey and Li (2000).

2. **The Universal Transverse Square-Cascade and 2D Lehmer Quadratic Bowl (Section 3, Theorems 3.2 & 3.4):** At any critical-line zero $t_0 \in \mathbb{R}$ of multiplicity $m \ge 1$ ($\Xi(t_0) = \dots = \Xi^{(m-1)}(t_0) = 0$, $\Xi^{(m)}(t_0) \neq 0$), all lower transverse dipole moments vanish identically ($c_1(t_0) = \dots = c _ {2m-3}(t_0) = 0$) while the leading non-vanishing coefficient collapses to the strictly positive pure square

$$c _ {2m-1}(t_0) = \frac{m}{(m!)^2}\left(\Xi^{(m)}(t_0)\right)^2 > 0.$$

Around the midpoint $t _ {\mathrm{mid}}$ of any closely spaced Lehmer pair $t _ \pm = t _ {\mathrm{mid}} \pm \delta$ ($\delta > 0$) with $\tau = t - t _ {\mathrm{mid}}$ and $\Xi''(t _ {\mathrm{mid}}) \neq 0$, the parity-reduced Wronskian forms an unbroken two-dimensional isotropic sum of three squares:

$$\widetilde{W}(d, t _ {\mathrm{mid}} + \tau) = \frac{1}{2}\left(\Xi''(t _ {\mathrm{mid}})\right)^2\left(\delta^2 + \tau^2 + d^2\right) + \mathcal{O}\left((\delta + \left|\tau\right| + d)^3\right) > 0.$$

3. **The Universal Leibniz–Wronskian Identity and 0%-Overshoot Prime-Gauge Sieve (Section 4, Theorems 4.1 & 4.4):** For any real-analytic gauge multiplier $\mathcal{M}(z)$ and $\Xi _ {\mathcal{M}}(z) := \mathcal{M}(z)\Xi(z)$, the transverse Wronskian splits without cross-terms:

$$\widetilde{W} _ {\mathcal{M}}(d,t) \equiv \left|\mathcal{M}(t - \mathrm{i}d)\right|^2 \widetilde{W}(d,t) + \left|\Xi(t - \mathrm{i}d)\right|^2 \widetilde{W} _ {1,\mathcal{M}}(d,t).$$

In particular, $\widetilde{W} _ {\mathcal{M}}(d_0, t_0) > 0$ certifies $\Xi(t_0 - \mathrm{i}d_0) \neq 0$ unconditionally, and any parasitic zero $\mathcal{M}(z_0) = 0$ produces only a harmless local false negative ($\widetilde{W} _ {\mathcal{M}}(d_0, t_0) = 0$). Choosing the entire zero-free exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {2 \le d \le X}\frac{\Lambda(d)}{\ln d}d^{-s})$ yields residual Dirichlet coefficients $b_n^{\mathcal{E}}(X)$ satisfying $b_n^{\mathcal{E}}(X) = 0$ for $2 \le n \le X$, $0 \le b_n^{\mathcal{E}}(X) \le 1$ for all $n \ge 1$ (0% overshoot), and $b_n^{\mathcal{E}}(X) = \frac{\Lambda(n)}{\ln n} \in [0,1]$ on the quadratic window $X < n < (\lfloor X\rfloor + 1)^2$, eliminating the $n = 2$ cross-term dip and closing the Lehmer uncertainty zone across all $1{,}002{,}000$ continuum cells on $[0, 0.50] \times [14, 1000]$ under the explicit Gabcke–Olver remainder bound ($\min(\widehat{W} _ {X,\mathrm{main}} - B _ {\widehat{W}}^{(C_0)}) = +0.5696 > 0$).

## Repository Architecture

```text
paper4-debranges-kernels-lehmer-cascade/
├── paper/
│   ├── paper4_debranges_lehmer_cascade.tex          # Standalone AMS-LaTeX manuscript
│   ├── paper4_debranges_lehmer_cascade.pdf          # Compiled 18-page publication PDF
│   └── figures/
│       ├── fig1_debranges_kernel_and_lehmer_cascade.pdf
│       ├── fig1_debranges_kernel_and_lehmer_cascade.png
│       ├── fig2_leibniz_wronskian_and_gabcke_closure.pdf
│       └── fig2_leibniz_wronskian_and_gabcke_closure.png
├── scripts/
│   └── verify_paper4_debranges_lehmer_cascade.py    # 50-digit & 550-digit mpmath + GPU verification suite
├── .agents/
│   └── skills/
│       └── debranges-kernels-lehmer-cascade/
│           └── SKILL.md                             # AI Agent mathematical skill definition
├── AGENT_SKILL_MANIFEST.json                        # Machine-readable agent skill manifest
├── CITATION.cff                                     # Citation File Format metadata
├── .zenodo.json                                     # Zenodo archival metadata
└── LICENSE                                          # MIT License
```

## Reproducible Verification

To run the complete verification suite (50-decimal and 550-decimal `mpmath` benchmarks, $2{,}000{,}000$-point GPU De Branges kernel sweep, $1{,}000$ batched $6 \times 6$ Cauchy–Pick Gram matrix spectra, $1{,}002{,}000$-cell Gabcke–Olver continuum certification, and regeneration of Figures 1 and 2):

```bash
python scripts/verify_paper4_debranges_lehmer_cascade.py
```

## AI Agent Skill Integration

This repository includes a machine-readable AI Agent Skill at [`.agents/skills/debranges-kernels-lehmer-cascade/SKILL.md`](.agents/skills/debranges-kernels-lehmer-cascade/SKILL.md) indexed by [`AGENT_SKILL_MANIFEST.json`](AGENT_SKILL_MANIFEST.json), enabling autonomous mathematical coding agents to load the exact De Branges kernel identity, square-cascade formulas, and 0%-overshoot prime-gauge recursions directly into context.

## License and Citation

Released under the [MIT License](LICENSE). Please cite using the metadata in [`CITATION.cff`](CITATION.cff).
