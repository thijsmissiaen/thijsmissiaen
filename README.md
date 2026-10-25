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
| **7** | **de la Vallée Poussin $\sigma = 1$ Pole Barrier & Trudgian Variational Coupling** | de la Vallée Poussin (1899); Guinand (1948); Trudgian (2014) | Mirror-paired Hadamard kernel $K(d, d _ \rho, u)$ with universal conservation law $\int _ {\mathbb{R}} K \mathrm{d}u \equiv \pi$, minorant $K _ {\min}(d,u)$ on $d \ge \frac{1}{2}$, and $d$-independent stationary trajectory $H _ {\mathrm{opt}}(kt) \equiv 2\pi B_S(kt)/\ln(kt/2\pi)$ (**Paper 7**). |
| **8** | **$\mathcal{O}(N(t)^2)$ Cross-Term Explosion, $1+\sqrt{2}$ Root Penalty & $L^1$ Barrier $\frac{1}{2}+d < 1$** | Montgomery & Vaughan (1974); Preissmann (1984); Conrey (1989) | Bilinear rank-2 factorization $w _ {mn}^{\mathrm{int}} \equiv \frac{2}{d^2}(\cosh(d\omega_m)\cosh(d\omega_n)-1)$, pure-square completion $+\frac{1}{d^2}\lVert\Psi _ {\ge 2}\rVert^2 \ge 0$, and $1+2d > 1$ $L^2$ exponent jump (**Paper 8**). |
| **9** | **Phase-Dependent Riemann–Siegel Branch Interference & $\mathcal{O}(t^{1/2})$ Prime-Sum Ceiling** | Siegel (1932); Rosser & Schoenfeld (1962); Dusart (2018) | Phase-independent confocal ellipse identity $\lvert T_n(d,t)\rvert^2 - \lvert T_n(0,t)\rvert^2 \equiv 4n^{-1}\sinh^2(d\omega_n)$, three-part $u$-domain bridge ($\theta_c < \frac{1}{4}$), and explicit Stieltjes bounds on $[10^8, \infty)$ (**Paper 9**). |
| **10** | **Bohr–Kronecker Torus Anti-Phase Paradox & Off-Critical Hadamard Flux Defect** | Bohr & Jessen (1930); Voronin (1975); Krein (1949) | Kac–Weyl recurrence barrier $\exp(C\sqrt{t}/\ln t) \gg t$, semiprime sum-of-squares floor $\ge \frac{1}{4}$ on $V=0$, quantized flux step $\mathcal{F}(d,d _ \star) = \pi\Theta(d-d _ \star)$ via $\mp\arcsin(d/d _ \star)$, and 1D Chladni wall identity (**Paper 10**). |

## Selected Research

### Paper 1: Wigner Phase-Space Representation and Complex Saddle-Point Coalescence in the Transverse Wronskian of Riemann's $\Xi$-Function

We prove that the inner correlation integral $I(y,t)$ governing the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ in mirror coordinates $(x=u+v, y=u-v)$ is proportional to the quantum-mechanical Wigner phase-space distribution $\mathcal{W} _ \Phi(y, 2t)$ of Riemann's Jacobi-theta kernel $\Phi(u)$, satisfying the zero-marginal identity $\int _ {-\infty}^\infty I(y,\gamma_k) \mathrm{d}y = 2\lvert\Xi(\gamma_k)\rvert^2 = 0$ at every real zero $\gamma_k$. Shifting the integration contour into the holomorphic strip $\lvert\mathrm{Im}(x)\rvert < \pi/4$ (where the vertical Cauchy walls decay double-exponentially as $\mathcal{O}(\mathrm{e}^{9L}\exp(-2\pi\cosh(2y)\mathrm{e}^{2L}\cos(2\eta)){}$) reveals a pitchfork saddle-point coalescence at $y_c(\omega)$ onto the imaginary axis $\xi _ {\star} = 0$, eliminating the 99.77% real-axis Fourier cancellation.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23637324`](https://doi.org/10.5281/zenodo.23637324)
* **Package Directory**: [`papers/paper01-wigner-saddle-coalescence/`](./papers/paper01-wigner-saddle-coalescence/)
* **Academic PDF Manuscript (13 pages)**: [`paper1_wigner_saddle_coalescence.pdf`](./papers/paper01-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.pdf)
* **AMS-LaTeX Source**: [`paper1_wigner_saddle_coalescence.tex`](./papers/paper01-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.tex)
* **60-Decimal Verification Script & Figure Generator**: [`verify_paper1_wigner_saddle_coalescence.py`](./papers/paper01-wigner-saddle-coalescence/scripts/verify_paper1_wigner_saddle_coalescence.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper01-wigner-saddle-coalescence/.agents/skills/wigner-saddle-coalescence-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper01-wigner-saddle-coalescence/AGENT_SKILL_MANIFEST.json)

### Paper 2: Parity-Reduced Transverse Wronskians, Higher Dipole Invariants, and Rigorous Continuum Certification for Riemann's $\Xi$-Function

We prove that the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ satisfies odd transverse parity $W(-d,t) = -W(d,t)$ with identically vanishing even derivatives $\partial_d^{2k}W(0,t) \equiv 0$, so that the parity-reduced Wronskian $\widetilde{W}(d,t) := W(d,t)/d$ resolves the $d \to 0^+$ boundary-layer degeneracy in the Speiser–Hinkkanen–Lagarias positivity criterion. Combining a division-free Wigner–hyperbolic integral representation governed by $\mathrm{sinc} _ {\mathrm{hyp}}(dy) = \frac{\sinh(dy)}{dy} \ge 1$, closed-form Cauchy–Riemann identities for the higher transverse dipole invariants $c _ {2m-1}(t)$ and non-singular partial derivatives, and dual directed-rounding interval evaluators (`mpmath.iv`), we establish strict positivity $c_1(\gamma_k), c_3(\gamma_k), c_5(\gamma_k) > 0$ across the first ten non-trivial zeros and certify the absence of off-critical zeros across $4{,}987{,}101$ continuum cells covering $[0, 0.50] \times [0, 55]$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641698`](https://doi.org/10.5281/zenodo.23641698)
* **Dedicated Repository**: [`thijsmissiaen/paper2-parity-reduced-wronskian`](https://github.com/thijsmissiaen/paper2-parity-reduced-wronskian)
* **Package Directory**: [`papers/paper02-parity-reduced-wronskian-certification/`](./papers/paper02-parity-reduced-wronskian-certification/)
* **Academic PDF Manuscript (15 pages)**: [`paper2_parity_reduced_wronskian.pdf`](./papers/paper02-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.pdf)
* **AMS-LaTeX Source**: [`paper2_parity_reduced_wronskian.tex`](./papers/paper02-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.tex)
* **Verification Script & Figure Generator**: [`verify_paper2_parity_reduced_wronskian.py`](./papers/paper02-parity-reduced-wronskian-certification/scripts/verify_paper2_parity_reduced_wronskian.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper02-parity-reduced-wronskian-certification/.agents/skills/parity-reduced-wronskian-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper02-parity-reduced-wronskian-certification/AGENT_SKILL_MANIFEST.json)

### Paper 3: Harmonic Nodal Chambers, Parabolic Zero-Curve Curvature, and the Transverse Integral Law for Riemann's $\Xi$-Function

Decomposing $\Xi(t-\mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$ into its coupled harmonic potentials ($\Delta R = \Delta I = 0$), we prove that the imaginary zero loci $\Gamma_k$ bifurcating orthogonally from the non-degenerate critical-line extrema $\tau_k$ partition the half-strip $d > 0$ into harmonic nodal chambers $\Omega_k$ along whose bounding walls $\lvert R\rvert$ is strictly monotonically increasing with arc length ($\frac{\mathrm{d}\lvert R\rvert}{\mathrm{d}s} = \lvert\nabla I\rvert = \lvert\Xi'(z)\rvert > 0$). We establish the closed-form parabolic vertex curvature $\kappa(\gamma_k) = \frac{\Xi''(\gamma_k)}{\Xi'(\gamma_k)} = -\mathrm{Im}\bigl(2g_1(\frac{1}{2}+\mathrm{i}\gamma_k) + \frac{\zeta''(\frac{1}{2}+\mathrm{i}\gamma_k)}{\zeta'(\frac{1}{2}+\mathrm{i}\gamma_k)}\bigr)$ of the real zero curves $\mathcal{R} _ k$ and prove the Transverse Integral Law, showing that the weighted integral operator $\frac{2}{d^2}\int_0^d u(\cdot) \mathrm{d}u$ acts as the hyperbolic low-pass filter $\mathcal{F} _ {\mathrm{int}}(2dy) = \frac{\tanh(dy)}{dy} \in (0,1)$ and locks the ground-state energy into the $90^\circ$ phase-shifted pure squares $\Xi(t)^2 + d^2(\Xi'(t))^2$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641701`](https://doi.org/10.5281/zenodo.23641701)
* **Dedicated Repository**: [`thijsmissiaen/paper3-harmonic-nodal-chambers`](https://github.com/thijsmissiaen/paper3-harmonic-nodal-chambers)
* **Package Directory**: [`papers/paper03-harmonic-chladni-nodal-chambers/`](./papers/paper03-harmonic-chladni-nodal-chambers/)
* **Academic PDF Manuscript (16 pages)**: [`paper3_harmonic_nodal_chambers.pdf`](./papers/paper03-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.pdf)
* **AMS-LaTeX Source**: [`paper3_harmonic_nodal_chambers.tex`](./papers/paper03-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.tex)
* **Verification Script & Figure Generator**: [`verify_paper3_harmonic_nodal_chambers.py`](./papers/paper03-harmonic-chladni-nodal-chambers/scripts/verify_paper3_harmonic_nodal_chambers.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper03-harmonic-chladni-nodal-chambers/.agents/skills/harmonic-nodal-chambers-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper03-harmonic-chladni-nodal-chambers/AGENT_SKILL_MANIFEST.json)

### Paper 4: De Branges Reproducing Kernels, the Universal Transverse Square-Cascade for Lehmer Pairs, and the Leibniz–Wronskian Identity for Riemann's $\Xi$-Function

For the linear Hermite–Biehler candidate family $E _ \alpha(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z)$ ($\alpha > 0$), we prove that the diagonal De Branges reproducing kernel satisfies the exact identity $K _ \alpha(t-\mathrm{i}d, t-\mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W}(d,t)$, bypassing the Conrey–Li (2000) global Sonine positivity obstruction for local transverse zero exclusion. We prove the Universal Transverse Square-Cascade $c _ {2m-1}(t_0) = \frac{m}{(m!)^2}(\Xi^{(m)}(t_0))^2 > 0$ at any $m$-fold coalescence point and show that across any closely spaced Lehmer pair $t _ \pm = t _ {\mathrm{mid}} \pm \delta$, $\widetilde{W}(d,t _ {\mathrm{mid}}+\tau) = \frac{1}{2}(\Xi''(t _ {\mathrm{mid}}))^2(\delta^2 + \tau^2 + d^2) + \mathcal{O}((\delta+\lvert\tau\rvert+d)^3) > 0$ forms a self-healing two-dimensional sum-of-squares bowl. Finally, we prove the cross-term-free Universal Leibniz–Wronskian Identity $\widetilde{W} _ {\mathcal{M}}(d,t) \equiv \lvert\mathcal{M}(t-\mathrm{i}d)\rvert^2\widetilde{W}(d,t) + \lvert\Xi(t-\mathrm{i}d)\rvert^2\widetilde{W} _ {1,\mathcal{M}}(d,t)$ and the 0%-overshoot property $0 \le b_n^{\mathcal{E}}(X) \le 1$ of the entire exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {2\le d\le X}\frac{\Lambda(d)}{\ln d}d^{-s})$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641702`](https://doi.org/10.5281/zenodo.23641702)
* **Dedicated Repository**: [`thijsmissiaen/paper4-debranges-lehmer-cascade`](https://github.com/thijsmissiaen/paper4-debranges-lehmer-cascade)
* **Package Directory**: [`papers/paper04-debranges-kernels-lehmer-cascade/`](./papers/paper04-debranges-kernels-lehmer-cascade/)
* **Academic PDF Manuscript (18 pages)**: [`paper4_debranges_lehmer_cascade.pdf`](./papers/paper04-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.pdf)
* **AMS-LaTeX Source**: [`paper4_debranges_lehmer_cascade.tex`](./papers/paper04-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.tex)
* **Verification Script & Figure Generator**: [`verify_paper4_debranges_lehmer_cascade.py`](./papers/paper04-debranges-kernels-lehmer-cascade/scripts/verify_paper4_debranges_lehmer_cascade.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper04-debranges-kernels-lehmer-cascade/.agents/skills/debranges-kernels-lehmer-cascade/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper04-debranges-kernels-lehmer-cascade/AGENT_SKILL_MANIFEST.json)

### Paper 5: Universal Mother-Wavelet Decomposition of the Jacobi–Riemann Kernel, Smooth Jacobi–Sinh Shells, and Explicit Saddle-Point Remainder Bounds for Riemann's $\Xi$-Function

We prove that every shell $\phi_m(u) = m^{-1/2}\varphi_1(u + \frac{1}{2}\ln m)$ of Riemann's kernel $\Phi(u)$ is generated by affine logarithmic dilation of a single universal Mother Wavelet $\varphi_1(w) = (2\pi^2\mathrm{e}^{9w} - 3\pi\mathrm{e}^{5w})\exp(-\pi\mathrm{e}^{4w})$ with bilateral Laplace–Mellin transform $\mathcal{G}(s) = \frac{1}{8}s(s-1)\pi^{-s/2}\Gamma(\frac{s}{2})$. Exact polynomial cancellation $\lvert(1-s)(-s)\rvert \equiv \lvert s(s-1)\rvert$ yields the Symmetric Branch-Damping Inequality $\rho _ \Gamma(d,t)n^{-1/2+d} < n^{-1/2-d}$ across all active Riemann–Siegel modes $1 \le n \le N(t) = \lfloor\sqrt{t/2\pi}\rfloor$ for $t \ge 7$. We eliminate the $u=0$ odd-derivative kink via smooth self-mirrored Jacobi–Sinh shells $\Phi_m^{\mathrm{ABC}}(z)$, derive the $\mathcal{O}(t^{-3/2}\ln^2 t)$ cubic-logarithmic Stieltjes diagonal $\mathcal{D} _ {\mathrm{diag}}(0,t) = \frac{1}{6}L^3 + \gamma _ {\mathrm{E}}L^2 - 4\gamma _ 1^{\mathrm{St}}L + 4\gamma _ 2^{\mathrm{St}}$ and the $\mathcal{O}(t^{-3/2})$ quadratic self-damping of the boundary pair $(N-1,N)$, and prove closed-form incomplete-Gamma tail bounds together with the mirror-wall damping identity $\mathrm{sinc} _ {\mathrm{hyp}}(dL/2)/\mathrm{sinc} _ {\mathrm{hyp}}(dL) \equiv 1/\cosh(dL/2)$.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23641708`](https://doi.org/10.5281/zenodo.23641708)
* **Dedicated Repository**: [`thijsmissiaen/paper5-mother-wavelet-gabcke-bounds`](https://github.com/thijsmissiaen/paper5-mother-wavelet-gabcke-bounds)
* **Package Directory**: [`papers/paper05-mother-wavelet-gabcke-bounds/`](./papers/paper05-mother-wavelet-gabcke-bounds/)
* **Academic PDF Manuscript (16 pages)**: [`paper5_mother_wavelet_gabcke_bounds.pdf`](./papers/paper05-mother-wavelet-gabcke-bounds/paper/paper5_mother_wavelet_gabcke_bounds.pdf)
* **AMS-LaTeX Source**: [`paper5_mother_wavelet_gabcke_bounds.tex`](./papers/paper05-mother-wavelet-gabcke-bounds/paper/paper5_mother_wavelet_gabcke_bounds.tex)
* **Verification Script & Figure Generator**: [`verify_paper5_mother_wavelet_gabcke_bounds.py`](./papers/paper05-mother-wavelet-gabcke-bounds/scripts/verify_paper5_mother_wavelet_gabcke_bounds.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper05-mother-wavelet-gabcke-bounds/.agents/skills/mother-wavelet-gabcke-bounds/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper05-mother-wavelet-gabcke-bounds/AGENT_SKILL_MANIFEST.json)

### Paper 6: The Sieved Mirror Kernel $\Phi_X(u)$, the Zero-Free Exponential Prime-Gauge Factor $\mathcal{E} _ X(z)$, and Quadratic Prime-Power Purity on $(X, X^2]$ for Riemann's $\Xi$-Function

We couple the universal mother-wavelet decomposition $\Phi(u) = \sum _ {m=1}^\infty m^{-1/2}\varphi_1(u + \frac{1}{2}\ln m)$ with finite-divisor Euler-product sieving for the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$. We prove that multiplying $\Xi(z)$ by a self-mirrored Dirichlet polynomial $\mathcal{M} _ X(z) = M_X(\frac{1}{2}+\mathrm{i}z)M_X(\frac{1}{2}-\mathrm{i}z)$ corresponds in the $u$-domain to an even convolution kernel $\Phi _ X(u)$ preserving the holomorphic strip $\lvert\mathrm{Im}(u)\rvert < \frac{\pi}{8}$, while replacing $M_X(s)$ by the entire, zero-free exponential prime-gauge factor $E_X(s) = \exp(-\sum _ {n=2}^{\lfloor X\rfloor}\frac{\Lambda(n)}{\ln n}n^{-s})$ eliminates all parasitic mollifier zeros and establishes exact quadratic prime-power purity $b_N^{\mathcal{E}}(X) \equiv \frac{\Lambda(N)}{\ln N}$ across $\lfloor X\rfloor < N < (\lfloor X\rfloor + 1)^2$ with 0% overshoot ($0 \le b_N^{\mathcal{E}}(X) \le 1$). Furthermore, we derive the negative hyperbolic mirror weight $w_n^{\ast}(d_0) = -(n^{2d_0}-1)^{-1}$ that annihilates prime-power oscillations on $d = d_0$ with zero waterbed tail for $n > X$, prove power-law $L^2$ diagonal damping $\frac{X^{-2d}}{2d\ln X}$ and universal Sub-Selberg closure $\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) < \frac{1}{2}$ across the critical strip, and establish the finite-divisor truncation identity on $[1, N(t)]$ alongside the shifted outer factorization of $\mathcal{E} _ {X,\kappa}(z)$.

* **Dedicated Repository**: [`thijsmissiaen/paper6-sieved-mirror-prime-gauge`](https://github.com/thijsmissiaen/paper6-sieved-mirror-prime-gauge)
* **Package Directory**: [`papers/paper06-sieved-mirror-kernel-prime-gauge/`](./papers/paper06-sieved-mirror-kernel-prime-gauge/)
* **Academic PDF Manuscript (17 pages)**: [`paper6_sieved_mirror_prime_gauge.pdf`](./papers/paper06-sieved-mirror-kernel-prime-gauge/paper/paper6_sieved_mirror_prime_gauge.pdf)
* **AMS-LaTeX Source**: [`paper6_sieved_mirror_prime_gauge.tex`](./papers/paper06-sieved-mirror-kernel-prime-gauge/paper/paper6_sieved_mirror_prime_gauge.tex)
* **Verification Script & Figure Generator**: [`verify_paper6_sieved_mirror_prime_gauge.py`](./papers/paper06-sieved-mirror-kernel-prime-gauge/scripts/verify_paper6_sieved_mirror_prime_gauge.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper06-sieved-mirror-kernel-prime-gauge/.agents/skills/sieved-mirror-prime-gauge/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper06-sieved-mirror-kernel-prime-gauge/AGENT_SKILL_MANIFEST.json)

### Paper 7: The Mirror-Paired Hadamard Square-Denominator Identity, the Universal Conservation Law $\int _ {\mathbb{R}} K \mathrm{d}u \equiv \pi$, and Stationary Trudgian Variation in the Critical Strip

We establish an exact dual representation linking the normalized transverse Wronskian $\widehat{W}(d,t) = W(d,t)/\lvert\Xi(t-\mathrm{i}d)\rvert^2$ of Riemann's $\Xi$-function to Fejér–Lukács pure-square trigonometric multiplets $\mathcal{S} _ {\pm}(d,t)$—resolving the classical de la Vallée Poussin $\sigma = 1$ ($d = \frac{1}{2}$) pole barrier—on the Euler-product side and to mirror-paired Hadamard zero sums on the spectral side. Pairing each nontrivial zero $\rho = \frac{1}{2} + d _ \rho + \mathrm{i}\gamma _ \rho$ with its functional-equation mirror $1 - \bar{\rho} = \frac{1}{2} - d _ \rho + \mathrm{i}\gamma _ \rho$ yields the square-denominator kernel $K(d, d _ \rho, u) = d(d^2 - d _ \rho^2 + u^2)/[((d^2 - d _ \rho^2 + u^2)^2 + 4u^2 d _ \rho^2)]$, which satisfies the universal Guinand–Weil / Hadamard mirror-pair conservation law $\int _ {-\infty}^{\infty} K(d, d _ \rho, u) \mathrm{d}u \equiv \pi$ for all $d > \lvert d _ \rho\rvert$ and admits on $d \ge \frac{1}{2}$ the hypothesis-free minorant $K _ {\min}(d,u) := K(d, \frac{1}{2}, u)$ throughout $\lvert u\rvert \ge H _ {\min}(d) = \sqrt{(d^2 - 1/4)/3}$. Integrating by parts against the Riemann–von Mangoldt decomposition $\mathrm{d}N(\gamma) = \frac{1}{2\pi}\ln(\gamma/2\pi) \mathrm{d}\gamma + \mathrm{d}S(\gamma)$ with Trudgian's explicit bound $\lvert S(T)\rvert \le B_S(T)$ reveals that $\partial_H K _ {\min}(d,H) < 0$ factors out identically from stationarity, pinning the optimal window half-width to the $d$-independent trajectory $H _ {\mathrm{opt}}(kt) \equiv 2\pi B_S(kt)/\ln(kt/2\pi) \to 2\pi(0.112) = 0.70371675\dots$.

* **Dedicated Repository**: [`thijsmissiaen/paper7-mirror-hadamard-variation`](https://github.com/thijsmissiaen/paper7-mirror-hadamard-variation)
* **Package Directory**: [`papers/paper07-mirror-hadamard-variational-law/`](./papers/paper07-mirror-hadamard-variational-law/)
* **Academic PDF Manuscript (18 pages)**: [`paper7_mirror_hadamard_variation.pdf`](./papers/paper07-mirror-hadamard-variational-law/paper/paper7_mirror_hadamard_variation.pdf)
* **AMS-LaTeX Source**: [`paper7_mirror_hadamard_variation.tex`](./papers/paper07-mirror-hadamard-variational-law/paper/paper7_mirror_hadamard_variation.tex)
* **Verification Script & Figure Generator**: [`verify_paper7_mirror_hadamard_variation.py`](./papers/paper07-mirror-hadamard-variational-law/scripts/verify_paper7_mirror_hadamard_variation.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper07-mirror-hadamard-variational-law/.agents/skills/mirror-hadamard-variational-law/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper07-mirror-hadamard-variational-law/AGENT_SKILL_MANIFEST.json)

### Paper 8: Bilinear Rank-2 Factorization of Integrated Wronskian Kernels, Pure-Square Completion, and the $1+2d > 1$ Montgomery–Vaughan Exponent Jump in the Critical Strip

We resolve four classical bottlenecks governing the $\binom{N(t)}{2}$ off-diagonal cross-pairs of the order-$r$ Mobius-sieved Riemann–Siegel main sum up to $N(t) = \lfloor\sqrt{t/2\pi}\rfloor$. Integrating the parity-reduced Wronskian transversely over $u \in [0,d]$ factors the coupled hyperbolic mode weight into the exact bilinear rank-2 identity $w _ {mn}^{\mathrm{int}}(d,t) \equiv \frac{2}{d^2}(\cosh(d\omega_m)\cosh(d\omega_n)-1)$, collapsing $\mathcal{O}(N(t)^2)$ double sums into $\mathcal{O}(N(t))$ prefix sums and revealing a sharp Euler–Maclaurin phase transition $\propto N(t)^{1-4d}(\ln N(t))^{2(r-1)}$ at $d = \frac{1}{4}$ ($\sigma = \frac{3}{4}$). Furthermore, the $-1$ lower-limit boundary terms cancel pairwise against the $d=0$ critical-line ground state, completing the $(m,n \ge 2)$ sector into the non-negative Euclidean square $+\frac{1}{d^2}\lVert\Psi _ {\ge 2}\rVert^2 \ge 0$ (eliminating the $1+\sqrt{2}$ root-isolation penalty), while passing from $L^1$ to $L^2$ via the Montgomery–Vaughan–Preissmann inequality shifts the exponent from $\frac{1}{2}+d < 1$ to $2(\frac{1}{2}+d) = 1+2d > 1$.

* **Dedicated Repository**: [`thijsmissiaen/paper8-bilinear-rank2-wronskian`](https://github.com/thijsmissiaen/paper8-bilinear-rank2-wronskian)
* **Package Directory**: [`papers/paper08-bilinear-rank2-montgomery-vaughan/`](./papers/paper08-bilinear-rank2-montgomery-vaughan/)
* **Academic PDF Manuscript (18 pages)**: [`paper8_bilinear_rank2_wronskian.pdf`](./papers/paper08-bilinear-rank2-montgomery-vaughan/paper/paper8_bilinear_rank2_wronskian.pdf)
* **AMS-LaTeX Source**: [`paper8_bilinear_rank2_wronskian.tex`](./papers/paper08-bilinear-rank2-montgomery-vaughan/paper/paper8_bilinear_rank2_wronskian.tex)
* **Verification Script & Figure Generator**: [`verify_paper8_bilinear_rank2_wronskian.py`](./papers/paper08-bilinear-rank2-montgomery-vaughan/scripts/verify_paper8_bilinear_rank2_wronskian.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper08-bilinear-rank2-montgomery-vaughan/.agents/skills/bilinear-rank2-wronskian/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper08-bilinear-rank2-montgomery-vaughan/AGENT_SKILL_MANIFEST.json)

### Paper 9: The Phase-Independent Confocal Ellipse Identity $\lvert T_n(d,t)\rvert^2 - \lvert T_n(0,t)\rvert^2 \equiv 4n^{-1}\sinh^2(d\omega_n)$, the Three-Part $u$-Domain Bridge Identity, and Explicit Rosser–Schoenfeld Stieltjes Bounds

We establish the Three-Part $u$-Domain Bridge Identity $\frac{1}{2}\Xi(t-\mathrm{i}d)E_X(s) \equiv \mathcal{I} _ +^{\mathcal{E} _ X}(s) + E_X(s)F _ +(1-s) + \mathcal{B} _ {\mathrm{bridge}}^{\mathcal{E} _ X}(s)$ across shifted Jacobi mirror points $w = u + \frac{1}{2}\ln k = 0$, proving that each half-line shell integral $\mathcal{J} _ n(s)$ admits a closed incomplete-Gamma form decaying super-exponentially as $\mathcal{O}(\pi n^2\mathrm{e}^{-\pi n^2})$ and that the pure Wronskian excess $\cosh(d\omega_n)-1 = 2\sinh^2(\frac{1}{2}d\omega_n)$ induces $(\ln N(t))^{-3}$ horizon dilution ($\theta_c < \frac{1}{4}$). Using the exact algebraic factorization $\sinh^2(\frac{1}{2}d\omega_n)/\sinh^2(\frac{1}{2}d\omega_1) \equiv n^{-d}((1-(n/N _ \star)^d)/(1-N _ \star^{-d}))^2 \le n^{-d}$, we derive closed-form Rosser–Schoenfeld Stieltjes upper bounds across $t \in [10^8, \infty)$ without summing primes, and prove the Phase-Independent Confocal Ellipse Identity $\lvert T_n(d,t)\rvert^2 - \lvert T_n(0,t)\rvert^2 \equiv 4n^{-1}\sinh^2(d\omega_n(t)) > 0$ ($\forall\phi_n \in \mathbb{R}$) with $d$-invariant focal distance $c_n \equiv \pm 2n^{-1/2}$.

* **Dedicated Repository**: [`thijsmissiaen/paper9-confocal-ellipse-synthesis`](https://github.com/thijsmissiaen/paper9-confocal-ellipse-synthesis)
* **Package Directory**: [`papers/paper09-confocal-ellipse-bridge-synthesis/`](./papers/paper09-confocal-ellipse-bridge-synthesis/)
* **Academic PDF Manuscript (16 pages)**: [`paper9_confocal_ellipse_synthesis.pdf`](./papers/paper09-confocal-ellipse-bridge-synthesis/paper/paper9_confocal_ellipse_synthesis.pdf)
* **AMS-LaTeX Source**: [`paper9_confocal_ellipse_synthesis.tex`](./papers/paper09-confocal-ellipse-bridge-synthesis/paper/paper9_confocal_ellipse_synthesis.tex)
* **Verification Script & Figure Generator**: [`verify_paper9_confocal_ellipse_synthesis.py`](./papers/paper09-confocal-ellipse-bridge-synthesis/scripts/verify_paper9_confocal_ellipse_synthesis.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper09-confocal-ellipse-bridge-synthesis/.agents/skills/confocal-ellipse-bridge-synthesis/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper09-confocal-ellipse-bridge-synthesis/AGENT_SKILL_MANIFEST.json)

### Paper 10: Beyond the Bohr–Kronecker Torus Barrier: Semiprime Phase-Locking, Quantized Hadamard–Krein Flux Defects $\pi\Theta(d - d _ {\star})$, and the Exact Chladni–Wronskian Wall Identity for Riemann's $\Xi$-Function

We establish three structural closure mechanisms for the parity-reduced transverse Wronskian $\widetilde{W}(d,t) = \frac{1}{2d}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ beyond the classical Bohr–Kronecker and Voronin torus approximation barriers. First, we prove that the ergodic (Kac–Weyl) mean Bohr–Kronecker recurrence horizon for simultaneous $\varepsilon$-alignment of the $\pi(N(t)) - \pi(X)$ active prime phases grows super-polynomially as $\exp(C\frac{\sqrt{t}}{\ln t}(1+o(1))) \gg t$, while multiplicative phase-locking $\phi _ {pq}(t) = \phi_p(t) + \phi_q(t)$ in $\zeta(s)E_X(s) = \exp(\mathcal{P} _ {>X}(s))$ enforces the exact sum-of-squares identity $\lvert 1 + \mathcal{P} _ {>X} + \frac{1}{2}\mathcal{P} _ {>X}^2\rvert^2 = \frac{1}{4}((1+U)^2 + 1 - V^2)^2 + V^2(1+U)^2 \ge \frac{1}{4}$ on the real anti-phase axis $V = 0$. Second, we prove that the mirror-paired Hadamard–Poisson kernel of any off-critical zero pair $z _ {\pm} = \gamma _ {\star} \mp \mathrm{i}d _ {\star}$ integrates over $u \in (-\infty,\infty)$ to the quantized step function $\pi\Theta(d - d _ {\star})$ via exact arcsine lobe cancellation $\mp\arcsin(d/d _ {\star})$. Third, along every imaginary Chladni nodal wall $\Gamma_k : I(d, \tau_k(d)) \equiv 0$, implicit Cauchy–Riemann differentiation yields the exact 1D identity $\widetilde{W}(d, \tau_k(d)) \equiv R(d, \tau_k(d))\partial_t I(d, \tau_k(d))/d$ without Dirichlet triangle inequalities.

* **Dedicated Repository**: [`thijsmissiaen/paper10-spectral-resonance-chladni-closure`](https://github.com/thijsmissiaen/paper10-spectral-resonance-chladni-closure)
* **Package Directory**: [`papers/paper10-spectral-resonance-chladni-closure/`](./papers/paper10-spectral-resonance-chladni-closure/)
* **Academic PDF Manuscript (13 pages)**: [`paper10_spectral_resonance_chladni_closure.pdf`](./papers/paper10-spectral-resonance-chladni-closure/paper/paper10_spectral_resonance_chladni_closure.pdf)
* **AMS-LaTeX Source**: [`paper10_spectral_resonance_chladni_closure.tex`](./papers/paper10-spectral-resonance-chladni-closure/paper/paper10_spectral_resonance_chladni_closure.tex)
* **Verification Script & Figure Generator**: [`verify_paper10_spectral_resonance_chladni_closure.py`](./papers/paper10-spectral-resonance-chladni-closure/scripts/verify_paper10_spectral_resonance_chladni_closure.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper10-spectral-resonance-chladni-closure/.agents/skills/spectral-resonance-chladni-closure/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper10-spectral-resonance-chladni-closure/AGENT_SKILL_MANIFEST.json)

## Technical Stack & Formalization

* **Analysis & Computation**: Python (`mpmath`, `numpy`, `scipy`, `sympy`), PyTorch / CUDA for interval Lipschitz certifications.
* **Typesetting & Mathematical Formalization**: LaTeX (`amsart`, `tikz`, `pgfplots`), Tectonic / pdfTeX.

## Academic Dialogue & Critique

I welcome rigorous critique, counter-examples, and academic dialogue from researchers and mathematicians working on analytic number theory, spectral theory, and numerical verification.

Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com) &nbsp;&middot;&nbsp; ORCID: [`https://orcid.org/0009-0007-3678-0810`](https://orcid.org/0009-0007-3678-0810)
