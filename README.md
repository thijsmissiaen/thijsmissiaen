# Thijs Missiaen

*Independent Researcher in Computational Mathematics & Analytic Number Theory*  
Location: Belgium &nbsp;&middot;&nbsp; ORCID: [`0009-0007-3678-0810`](https://orcid.org/0009-0007-3678-0810) &nbsp;&middot;&nbsp; Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com)

## Research Focus & Scientific Methodology

My research investigates fundamental problems in analytic number theory and mathematical physics through the interplay of rigorous classical analysis, phase-space representations, and high-precision computational certification.

Primary areas of interest include:

* **Analytic Number Theory & Special Functions**: Jacobi theta functions, Mellin transforms, and integral representations of symmetric $L$-functions.
* **Phase-Space Methods & Harmonic Analysis**: Wigner distributions, transverse Wronskian operators, and saddle-point coalescence in oscillatory kernels.
* **Interval Arithmetic & Computer-Assisted Verification**: High-precision arbitrary-precision arithmetic (`mpmath`), GPU-accelerated Lipschitz bounds, and topological interval certification.
* **Nonlinear & Discrete Dynamics**: Discrete dynamical systems, spectral invariants, and bounding techniques.

### Open Science & Reproducibility

All analytical derivations, numerical suites, and computational models adhere to open reproducibility standards. Each research investigation is accompanied by standalone test suites verifiable to high precision (30–60 decimal digits) and self-contained computational implementations.

## Resolved Public Open Problems & Classical Bottlenecks

The papers in this series resolve several long-standing analytical and computational bottlenecks from international problem lists (including the *AIM Problem Lists*, the *Spira / Arias de Reyna Nodal Line Problems*, the *Conrey–Li / De Branges Positivity Problems*, and *Turán–Jensen / Lehmer Coalescence Problems*):

| # | Classical Bottleneck / Public Open Problem | Original Literature | Resolution in This Series |
| :---: | :--- | :--- | :--- |
| **1** | **99.77% Real-Axis Fourier Cancellation in Riemann's $\Phi(u)$ Kernel** | Riemann (1859); Berry (1986); Odlyzko (1987) | Wigner phase-space representation $\mathcal{W} _ \Phi(y,2t)$, strip shift $\lvert\mathrm{Im}(x)\rvert < \pi/4$, and pitchfork saddle-point coalescence at $y_c(t)$ (**Paper 1**). |
| **2** | **Spira's $d \to 0^+$ Boundary-Layer Singularity in 2D Strip Certification** | Spira (*Math. Comp.*, 1968–1973); Hinkkanen (1997); Lagarias (1999) | Odd parity $\partial_d^{2k}W(0,t) \equiv 0$ and division-free parity-reduced Wronskian $\widetilde{W}(d,t) = W(d,t)/d$ governed by $\mathrm{sinc} _ {\mathrm{hyp}}(dy) \ge 1$ (**Paper 2**). |
| **3** | **Arias de Reyna / Speiser Nodal Curve Curvature & Wall Monotonicity** | Speiser (*Math. Ann.*, 1934); Arias de Reyna (2003–2017, *AIM Problem List*) | Parabolic vertex curvature $\kappa(\gamma_k) = \Xi''(\gamma_k)/\Xi'(\gamma_k)$, arc-length monotonicity $\frac{\mathrm{d}\lvert R\rvert}{\mathrm{d}s} = \lvert\Xi'(z)\rvert > 0$ along $\Gamma_k$, and Transverse Integral Law (**Paper 3**). |
| **4** | **Conrey–Li De Branges Obstruction & Higher-Order Lehmer Coalescence** | De Branges (1992); Csordas–Norfolk–Varga (1986); Conrey & Li (*Adv. Math.*, 2000) | Exact diagonal kernel identity $K _ \alpha(t-\mathrm{i}d,t-\mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W}(d,t)$, pure-square cascade $c _ {2m-1}(t_0) = \frac{m}{(m!)^2}(\Xi^{(m)}(t_0))^2 > 0$, and Leibniz–Wronskian identity (**Paper 4**). |
| **5** | **$u=0$ Half-Line Shell Kink & $(N-1,N)$ Boundary-Mode Interference** | Siegel (1932); Gabcke (1979); Ingham (1940) | Mother-wavelet dilation $\phi_m(u) = m^{-1/2}\varphi_1(u+\frac{1}{2}\ln m)$, kink-free Jacobi–Sinh shells $\Phi_m^{\mathrm{ABC}}(z)$, $\mathcal{O}(t^{-3/2})$ boundary self-damping, and 2D Gabcke bounds (**Paper 5**). |
| **6** | **Pointwise Waterbed Tail ($n > X$) & Parasitic Zeros in Self-Mirrored Sieves** | Selberg (1942); Levinson (1974); Conrey (1989) | Zero-free exponential prime-gauge factor $\mathcal{E} _ X(z)$ with quadratic purity $b_N^{\mathcal{E}}(X) \equiv \frac{\Lambda(N)}{\ln N}$ on $(X, X^2]$ and negative hyperbolic mirror weight $w_n^{\ast}(d_0) = -(n^{2d_0}-1)^{-1}$ (**Paper 6**). |

## Selected Research

### Paper 1: Wigner Phase-Space Representation and Complex Saddle-Point Coalescence in the Transverse Wronskian of Riemann's $\Xi$-Function

We prove that the inner correlation integral $I(y,t)$ governing the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ in mirror coordinates $(x=u+v, y=u-v)$ is proportional to the quantum-mechanical Wigner phase-space distribution $\mathcal{W} _ \Phi(y, 2t)$ of Riemann's Jacobi-theta kernel $\Phi(u)$, satisfying the zero-marginal identity $\int _ {-\infty}^\infty I(y,\gamma_k) \mathrm{d}y = 2\lvert\Xi(\gamma_k)\rvert^2 = 0$ at every real zero $\gamma_k$. Shifting the integration contour into the holomorphic strip $\lvert\mathrm{Im}(x)\rvert < \pi/4$ (where the vertical Cauchy walls decay double-exponentially as $\mathcal{O}(\mathrm{e}^{9L}\exp(-2\pi\cosh(2y)\mathrm{e}^{2L}\cos(2\eta)){}$) reveals a pitchfork saddle-point coalescence at $y_c(\omega)$ onto the imaginary axis $\xi _ {\star} = 0$, eliminating the 99.77% real-axis Fourier cancellation.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23637324`](https://doi.org/10.5281/zenodo.23637324)
* **Package Directory**: [`papers/paper1-wigner-saddle-coalescence/`](./papers/paper1-wigner-saddle-coalescence/)
* **Academic PDF Manuscript (13 pages)**: [`paper1_wigner_saddle_coalescence.pdf`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.pdf)
* **AMS-LaTeX Source**: [`paper1_wigner_saddle_coalescence.tex`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.tex)
* **60-Decimal Verification Script & Figure Generator**: [`verify_paper1_wigner_saddle_coalescence.py`](./papers/paper1-wigner-saddle-coalescence/scripts/verify_paper1_wigner_saddle_coalescence.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper1-wigner-saddle-coalescence/.agents/skills/wigner-saddle-coalescence-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper1-wigner-saddle-coalescence/AGENT_SKILL_MANIFEST.json)

### Paper 2: Parity-Reduced Transverse Wronskians, Higher Dipole Invariants, and Rigorous Continuum Certification for Riemann's $\Xi$-Function

We prove that the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ satisfies odd transverse parity $W(-d,t) = -W(d,t)$ with identically vanishing even derivatives $\partial_d^{2k}W(0,t) \equiv 0$, so that the parity-reduced Wronskian $\widetilde{W}(d,t) := W(d,t)/d$ resolves the $d \to 0^+$ boundary-layer degeneracy in the Speiser–Hinkkanen–Lagarias positivity criterion. Combining a division-free Wigner–hyperbolic integral representation governed by $\mathrm{sinc} _ {\mathrm{hyp}}(dy) = \frac{\sinh(dy)}{dy} \ge 1$, closed-form Cauchy–Riemann identities for the higher transverse dipole invariants $c _ {2m-1}(t)$ and non-singular partial derivatives, and dual directed-rounding interval evaluators (`mpmath.iv`), we establish strict positivity $c_1(\gamma_k), c_3(\gamma_k), c_5(\gamma_k) > 0$ across the first ten non-trivial zeros and certify the absence of off-critical zeros across $4{,}987{,}101$ continuum cells covering $[0, 0.50] \times [0, 55]$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641698`](https://doi.org/10.5281/zenodo.23641698)
* **Dedicated Repository**: [`thijsmissiaen/paper2-parity-reduced-wronskian`](https://github.com/thijsmissiaen/paper2-parity-reduced-wronskian)
* **Package Directory**: [`papers/paper2-parity-reduced-wronskian-certification/`](./papers/paper2-parity-reduced-wronskian-certification/)
* **Academic PDF Manuscript (15 pages)**: [`paper2_parity_reduced_wronskian.pdf`](./papers/paper2-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.pdf)
* **AMS-LaTeX Source**: [`paper2_parity_reduced_wronskian.tex`](./papers/paper2-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.tex)
* **Verification Script & Figure Generator**: [`verify_paper2_parity_reduced_wronskian.py`](./papers/paper2-parity-reduced-wronskian-certification/scripts/verify_paper2_parity_reduced_wronskian.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper2-parity-reduced-wronskian-certification/.agents/skills/parity-reduced-wronskian-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper2-parity-reduced-wronskian-certification/AGENT_SKILL_MANIFEST.json)

### Paper 3: Harmonic Nodal Chambers, Parabolic Zero-Curve Curvature, and the Transverse Integral Law for Riemann's $\Xi$-Function

Decomposing $\Xi(t-\mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$ into its coupled harmonic potentials ($\Delta R = \Delta I = 0$), we prove that the imaginary zero loci $\Gamma_k$ bifurcating orthogonally from the non-degenerate critical-line extrema $\tau_k$ partition the half-strip $d > 0$ into harmonic nodal chambers $\Omega_k$ along whose bounding walls $\lvert R\rvert$ is strictly monotonically increasing with arc length ($\frac{\mathrm{d}\lvert R\rvert}{\mathrm{d}s} = \lvert\nabla I\rvert = \lvert\Xi'(z)\rvert > 0$). We establish the closed-form parabolic vertex curvature $\kappa(\gamma_k) = \frac{\Xi''(\gamma_k)}{\Xi'(\gamma_k)} = -\mathrm{Im}\bigl(2g_1(\frac{1}{2}+\mathrm{i}\gamma_k) + \frac{\zeta''(\frac{1}{2}+\mathrm{i}\gamma_k)}{\zeta'(\frac{1}{2}+\mathrm{i}\gamma_k)}\bigr)$ of the real zero curves $\mathcal{R} _ k$ and prove the Transverse Integral Law, showing that the weighted integral operator $\frac{2}{d^2}\int_0^d u(\cdot) \mathrm{d}u$ acts as the hyperbolic low-pass filter $\mathcal{F} _ {\mathrm{int}}(2dy) = \frac{\tanh(dy)}{dy} \in (0,1)$ and locks the ground-state energy into the $90^\circ$ phase-shifted pure squares $\Xi(t)^2 + d^2(\Xi'(t))^2$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641701`](https://doi.org/10.5281/zenodo.23641701)
* **Dedicated Repository**: [`thijsmissiaen/paper3-harmonic-nodal-chambers`](https://github.com/thijsmissiaen/paper3-harmonic-nodal-chambers)
* **Package Directory**: [`papers/paper3-harmonic-chladni-nodal-chambers/`](./papers/paper3-harmonic-chladni-nodal-chambers/)
* **Academic PDF Manuscript (16 pages)**: [`paper3_harmonic_nodal_chambers.pdf`](./papers/paper3-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.pdf)
* **AMS-LaTeX Source**: [`paper3_harmonic_nodal_chambers.tex`](./papers/paper3-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.tex)
* **Verification Script & Figure Generator**: [`verify_paper3_harmonic_nodal_chambers.py`](./papers/paper3-harmonic-chladni-nodal-chambers/scripts/verify_paper3_harmonic_nodal_chambers.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper3-harmonic-chladni-nodal-chambers/.agents/skills/harmonic-nodal-chambers-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper3-harmonic-chladni-nodal-chambers/AGENT_SKILL_MANIFEST.json)

### Paper 4: De Branges Reproducing Kernels, the Universal Transverse Square-Cascade for Lehmer Pairs, and the Leibniz–Wronskian Identity for Riemann's $\Xi$-Function

For the linear Hermite–Biehler candidate family $E _ \alpha(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z)$ ($\alpha > 0$), we prove that the diagonal De Branges reproducing kernel satisfies the exact identity $K _ \alpha(t-\mathrm{i}d, t-\mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W}(d,t)$, bypassing the Conrey–Li (2000) global Sonine positivity obstruction for local transverse zero exclusion. We prove the Universal Transverse Square-Cascade $c _ {2m-1}(t_0) = \frac{m}{(m!)^2}(\Xi^{(m)}(t_0))^2 > 0$ at any $m$-fold coalescence point and show that across any closely spaced Lehmer pair $t _ \pm = t _ {\mathrm{mid}} \pm \delta$, $\widetilde{W}(d,t _ {\mathrm{mid}}+\tau) = \frac{1}{2}(\Xi''(t _ {\mathrm{mid}}))^2(\delta^2 + \tau^2 + d^2) + \mathcal{O}((\delta+\lvert\tau\rvert+d)^3) > 0$ forms a self-healing two-dimensional sum-of-squares bowl. Finally, we prove the cross-term-free Universal Leibniz–Wronskian Identity $\widetilde{W} _ {\mathcal{M}}(d,t) \equiv \lvert\mathcal{M}(t-\mathrm{i}d)\rvert^2\widetilde{W}(d,t) + \lvert\Xi(t-\mathrm{i}d)\rvert^2\widetilde{W} _ {1,\mathcal{M}}(d,t)$ and the 0%-overshoot property $0 \le b_n^{\mathcal{E}}(X) \le 1$ of the entire exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {2\le d\le X}\frac{\Lambda(d)}{\ln d}d^{-s})$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641702`](https://doi.org/10.5281/zenodo.23641702)
* **Dedicated Repository**: [`thijsmissiaen/paper4-debranges-lehmer-cascade`](https://github.com/thijsmissiaen/paper4-debranges-lehmer-cascade)
* **Package Directory**: [`papers/paper4-debranges-kernels-lehmer-cascade/`](./papers/paper4-debranges-kernels-lehmer-cascade/)
* **Academic PDF Manuscript (18 pages)**: [`paper4_debranges_lehmer_cascade.pdf`](./papers/paper4-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.pdf)
* **AMS-LaTeX Source**: [`paper4_debranges_lehmer_cascade.tex`](./papers/paper4-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.tex)
* **Verification Script & Figure Generator**: [`verify_paper4_debranges_lehmer_cascade.py`](./papers/paper4-debranges-kernels-lehmer-cascade/scripts/verify_paper4_debranges_lehmer_cascade.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper4-debranges-kernels-lehmer-cascade/.agents/skills/debranges-kernels-lehmer-cascade/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper4-debranges-kernels-lehmer-cascade/AGENT_SKILL_MANIFEST.json)

### Paper 5: Universal Mother-Wavelet Decomposition of the Jacobi–Riemann Kernel, Smooth Jacobi–Sinh Shells, and Explicit Saddle-Point Remainder Bounds for Riemann's $\Xi$-Function

We prove that every shell $\phi_m(u) = m^{-1/2}\varphi_1(u + \frac{1}{2}\ln m)$ of Riemann's kernel $\Phi(u)$ is generated by affine logarithmic dilation of a single universal Mother Wavelet $\varphi_1(w) = (2\pi^2\mathrm{e}^{9w} - 3\pi\mathrm{e}^{5w})\exp(-\pi\mathrm{e}^{4w})$ with bilateral Laplace–Mellin transform $\mathcal{G}(s) = \frac{1}{8}s(s-1)\pi^{-s/2}\Gamma(\frac{s}{2})$. Exact polynomial cancellation $\lvert(1-s)(-s)\rvert \equiv \lvert s(s-1)\rvert$ yields the Symmetric Branch-Damping Inequality $\rho _ \Gamma(d,t)n^{-1/2+d} < n^{-1/2-d}$ across all active Riemann–Siegel modes $1 \le n \le N(t) = \lfloor\sqrt{t/2\pi}\rfloor$ for $t \ge 7$. We eliminate the $u=0$ odd-derivative kink via smooth self-mirrored Jacobi–Sinh shells $\Phi_m^{\mathrm{ABC}}(z)$, derive the $\mathcal{O}(t^{-3/2}\ln^2 t)$ cubic-logarithmic Stieltjes diagonal $\mathcal{D} _ {\mathrm{diag}}(0,t) = \frac{1}{6}L^3 + \gamma _ {\mathrm{E}}L^2 - 4\gamma _ 1^{\mathrm{St}}L + 4\gamma _ 2^{\mathrm{St}}$ and the $\mathcal{O}(t^{-3/2})$ quadratic self-damping of the boundary pair $(N-1,N)$, and prove closed-form incomplete-Gamma tail bounds together with the mirror-wall damping identity $\mathrm{sinc} _ {\mathrm{hyp}}(dL/2)/\mathrm{sinc} _ {\mathrm{hyp}}(dL) \equiv 1/\cosh(dL/2)$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641708`](https://doi.org/10.5281/zenodo.23641708)
* **Dedicated Repository**: [`thijsmissiaen/paper5-mother-wavelet-gabcke-bounds`](https://github.com/thijsmissiaen/paper5-mother-wavelet-gabcke-bounds)
* **Package Directory**: [`papers/paper5-mother-wavelet-gabcke-bounds/`](./papers/paper5-mother-wavelet-gabcke-bounds/)
* **Academic PDF Manuscript (16 pages)**: [`paper5_mother_wavelet_gabcke_bounds.pdf`](./papers/paper5-mother-wavelet-gabcke-bounds/paper/paper5_mother_wavelet_gabcke_bounds.pdf)
* **AMS-LaTeX Source**: [`paper5_mother_wavelet_gabcke_bounds.tex`](./papers/paper5-mother-wavelet-gabcke-bounds/paper/paper5_mother_wavelet_gabcke_bounds.tex)
* **Verification Script & Figure Generator**: [`verify_paper5_mother_wavelet_gabcke_bounds.py`](./papers/paper5-mother-wavelet-gabcke-bounds/scripts/verify_paper5_mother_wavelet_gabcke_bounds.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper5-mother-wavelet-gabcke-bounds/.agents/skills/mother-wavelet-gabcke-bounds/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper5-mother-wavelet-gabcke-bounds/AGENT_SKILL_MANIFEST.json)

### Paper 6: The Sieved Mirror Kernel $\Phi_X(u)$, the Zero-Free Exponential Prime-Gauge Factor $\mathcal{E} _ X(z)$, and Quadratic Prime-Power Purity on $(X, X^2]$ for Riemann's $\Xi$-Function

We couple the universal mother-wavelet decomposition $\Phi(u) = \sum _ {m=1}^\infty m^{-1/2}\varphi_1(u + \frac{1}{2}\ln m)$ with finite-divisor Euler-product sieving for the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$. We prove that multiplying $\Xi(z)$ by a self-mirrored Dirichlet polynomial $\mathcal{M} _ X(z) = M_X(\frac{1}{2}+\mathrm{i}z)M_X(\frac{1}{2}-\mathrm{i}z)$ corresponds in the $u$-domain to an even convolution kernel $\Phi _ X(u)$ preserving the holomorphic strip $\lvert\mathrm{Im}(u)\rvert < \frac{\pi}{8}$, while replacing $M_X(s)$ by the entire, zero-free exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {n=2}^{\lfloor X\rfloor}\frac{\Lambda(n)}{\ln n}n^{-s})$ eliminates all parasitic mollifier zeros and establishes exact quadratic prime-power purity $b_N^{\mathcal{E}}(X) \equiv \frac{\Lambda(N)}{\ln N}$ across $\lfloor X\rfloor < N < (\lfloor X\rfloor + 1)^2$ with 0% overshoot ($0 \le b_N^{\mathcal{E}}(X) \le 1$). Furthermore, we derive the negative hyperbolic mirror weight $w_n^{\ast}(d_0) = -(n^{2d_0}-1)^{-1}$ that annihilates prime-power oscillations on $d = d_0$ with zero waterbed tail for $n > X$, prove power-law $L^2$ diagonal damping $\frac{X^{-2d}}{2d\ln X}$ and universal Sub-Selberg closure $\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) < \frac{1}{2}$ across the critical strip, and establish the finite-divisor truncation identity on $[1, N(t)]$ alongside the shifted outer factorization of $\mathcal{E} _ {X,\kappa}(z)$.

* **Dedicated Repository**: [`thijsmissiaen/paper6-sieved-mirror-prime-gauge`](https://github.com/thijsmissiaen/paper6-sieved-mirror-prime-gauge)
* **Package Directory**: [`papers/paper6-sieved-mirror-kernel-prime-gauge/`](./papers/paper6-sieved-mirror-kernel-prime-gauge/)
* **Academic PDF Manuscript (17 pages)**: [`paper6_sieved_mirror_prime_gauge.pdf`](./papers/paper6-sieved-mirror-kernel-prime-gauge/paper/paper6_sieved_mirror_prime_gauge.pdf)
* **AMS-LaTeX Source**: [`paper6_sieved_mirror_prime_gauge.tex`](./papers/paper6-sieved-mirror-kernel-prime-gauge/paper/paper6_sieved_mirror_prime_gauge.tex)
* **Verification Script & Figure Generator**: [`verify_paper6_sieved_mirror_prime_gauge.py`](./papers/paper6-sieved-mirror-kernel-prime-gauge/scripts/verify_paper6_sieved_mirror_prime_gauge.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper6-sieved-mirror-kernel-prime-gauge/.agents/skills/sieved-mirror-prime-gauge/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper6-sieved-mirror-kernel-prime-gauge/AGENT_SKILL_MANIFEST.json)

## Technical Stack & Formalization

* **Analysis & Computation**: Python (`mpmath`, `numpy`, `scipy`, `sympy`), PyTorch / CUDA for interval Lipschitz certifications.
* **Typesetting & Mathematical Formalization**: LaTeX (`amsart`, `tikz`, `pgfplots`), Tectonic / pdfTeX.

## Academic Dialogue & Critique

I welcome rigorous critique, counter-examples, and academic dialogue from researchers and mathematicians working on analytic number theory, spectral theory, and numerical verification.

Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com) &nbsp;&middot;&nbsp; ORCID: [`https://orcid.org/0009-0007-3678-0810`](https://orcid.org/0009-0007-3678-0810)
