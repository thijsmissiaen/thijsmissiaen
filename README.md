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

## Selected Research

### Paper 1: Wigner Phase-Space Representation and Complex Saddle-Point Coalescence in the Transverse Wronskian of Riemann's $\Xi$-Function

We prove that the inner correlation integral $I(y,t)$ governing the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ in mirror coordinates $(x=u+v, y=u-v)$ is proportional to the quantum-mechanical Wigner phase-space distribution $\mathcal{W}_\Phi(y, 2t)$ of Riemann's Jacobi-theta kernel $\Phi(u)$, satisfying the zero-marginal identity $\int_{-\infty}^\infty I(y,\gamma_k)\,\mathrm{d}y = 2\lvert\Xi(\gamma_k)\rvert^2 = 0$ at every real zero $\gamma_k$. Shifting the integration contour into the holomorphic strip $\lvert\mathrm{Im}(x)\rvert < \pi/4$ (where the vertical Cauchy walls decay double-exponentially as $\mathcal{O}(\mathrm{e}^{9L}\exp(-2\pi\cosh(2y)\mathrm{e}^{2L}\cos(2\eta)))$) reveals a pitchfork saddle-point coalescence at $y_c(\omega)$ onto the imaginary axis $\xi_* = 0$, eliminating the $99.77\%$ real-axis Fourier cancellation.

* **DOI (Zenodo / CERN)**: [`10.5281/zenodo.23637324`](https://doi.org/10.5281/zenodo.23637324)
* **Package Directory**: [`papers/paper1-wigner-saddle-coalescence/`](./papers/paper1-wigner-saddle-coalescence/)
* **Academic PDF Manuscript (13 pages)**: [`paper1_wigner_saddle_coalescence.pdf`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.pdf)
* **AMS-LaTeX Source**: [`paper1_wigner_saddle_coalescence.tex`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.tex)
* **60-Decimal Verification Script & Figure Generator**: [`verify_paper1_wigner_saddle_coalescence.py`](./papers/paper1-wigner-saddle-coalescence/scripts/verify_paper1_wigner_saddle_coalescence.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper1-wigner-saddle-coalescence/.agents/skills/wigner-saddle-coalescence-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper1-wigner-saddle-coalescence/AGENT_SKILL_MANIFEST.json)

### Paper 2: Parity-Reduced Transverse Wronskians, Higher Dipole Invariants, and Rigorous Continuum Certification for Riemann's $\Xi$-Function

We prove that the transverse Wronskian $W(d,t) = \frac{1}{2}\partial_d\lvert\Xi(t-\mathrm{i}d)\rvert^2$ satisfies odd transverse parity $W(-d,t) = -W(d,t)$ with identically vanishing even derivatives $\partial_d^{2k}W(0,t) \equiv 0$, so that the parity-reduced Wronskian $\widetilde{W}(d,t) := W(d,t)/d$ resolves the $d \to 0^+$ boundary-layer degeneracy in the Speiser–Hinkkanen–Lagarias positivity criterion. Combining a division-free Wigner–hyperbolic integral representation governed by $\mathrm{sinc}_{\mathrm{hyp}}(dy) = \frac{\sinh(dy)}{dy} \ge 1$, closed-form Cauchy–Riemann identities for the higher transverse dipole invariants $c_{2m-1}(t)$ and non-singular partial derivatives, and dual directed-rounding interval evaluators (`mpmath.iv`), we establish strict positivity $c_1(\gamma_k), c_3(\gamma_k), c_5(\gamma_k) > 0$ across the first ten non-trivial zeros and certify the absence of off-critical zeros across $4{,}987{,}101$ continuum cells covering $[0, 0.50] \times [0, 55]$.

* **Dedicated Repository**: [`thijsmissiaen/paper2-parity-reduced-wronskian`](https://github.com/thijsmissiaen/paper2-parity-reduced-wronskian)
* **Package Directory**: [`papers/paper2-parity-reduced-wronskian-certification/`](./papers/paper2-parity-reduced-wronskian-certification/)
* **Academic PDF Manuscript (15 pages)**: [`paper2_parity_reduced_wronskian.pdf`](./papers/paper2-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.pdf)
* **AMS-LaTeX Source**: [`paper2_parity_reduced_wronskian.tex`](./papers/paper2-parity-reduced-wronskian-certification/paper/paper2_parity_reduced_wronskian.tex)
* **Verification Script & Figure Generator**: [`verify_paper2_parity_reduced_wronskian.py`](./papers/paper2-parity-reduced-wronskian-certification/scripts/verify_paper2_parity_reduced_wronskian.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper2-parity-reduced-wronskian-certification/.agents/skills/parity-reduced-wronskian-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper2-parity-reduced-wronskian-certification/AGENT_SKILL_MANIFEST.json)

### Paper 3: Harmonic Nodal Chambers, Parabolic Zero-Curve Curvature, and the Transverse Integral Law for Riemann's $\Xi$-Function

Decomposing $\Xi(t-\mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$ into its coupled harmonic potentials ($\Delta R = \Delta I = 0$), we prove that the imaginary zero loci $\Gamma_k$ bifurcating orthogonally from the non-degenerate critical-line extrema $\tau_k$ partition the half-strip $d > 0$ into harmonic nodal chambers $\Omega_k$ along whose bounding walls $\lvert R\rvert$ is strictly monotonically increasing with arc length ($\frac{\mathrm{d}\lvert R\rvert}{\mathrm{d}s} = \lvert\nabla I\rvert = \lvert\Xi'(z)\rvert > 0$). We establish the closed-form parabolic vertex curvature $\kappa(\gamma_k) = \frac{\Xi''(\gamma_k)}{\Xi'(\gamma_k)} = -\mathrm{Im}\bigl(2g_1(\frac{1}{2}+\mathrm{i}\gamma_k) + \frac{\zeta''(\frac{1}{2}+\mathrm{i}\gamma_k)}{\zeta'(\frac{1}{2}+\mathrm{i}\gamma_k)}\bigr)$ of the real zero curves $\mathcal{R}_k$ and prove the Transverse Integral Law, showing that the weighted integral operator $\frac{2}{d^2}\int_0^d u(\cdot)\,\mathrm{d}u$ acts as the hyperbolic low-pass filter $\mathcal{F}_{\mathrm{int}}(2dy) = \frac{\tanh(dy)}{dy} \in (0,1)$ and locks the ground-state energy into the $90^\circ$ phase-shifted pure squares $\Xi(t)^2 + d^2(\Xi'(t))^2$.

* **Dedicated Repository**: [`thijsmissiaen/paper3-harmonic-nodal-chambers`](https://github.com/thijsmissiaen/paper3-harmonic-nodal-chambers)
* **Package Directory**: [`papers/paper3-harmonic-chladni-nodal-chambers/`](./papers/paper3-harmonic-chladni-nodal-chambers/)
* **Academic PDF Manuscript (16 pages)**: [`paper3_harmonic_nodal_chambers.pdf`](./papers/paper3-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.pdf)
* **AMS-LaTeX Source**: [`paper3_harmonic_nodal_chambers.tex`](./papers/paper3-harmonic-chladni-nodal-chambers/paper/paper3_harmonic_nodal_chambers.tex)
* **Verification Script & Figure Generator**: [`verify_paper3_harmonic_nodal_chambers.py`](./papers/paper3-harmonic-chladni-nodal-chambers/scripts/verify_paper3_harmonic_nodal_chambers.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper3-harmonic-chladni-nodal-chambers/.agents/skills/harmonic-nodal-chambers-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper3-harmonic-chladni-nodal-chambers/AGENT_SKILL_MANIFEST.json)

### Paper 4: De Branges Reproducing Kernels, the Universal Transverse Square-Cascade for Lehmer Pairs, and the Leibniz–Wronskian Identity for Riemann's $\Xi$-Function

For the linear Hermite–Biehler candidate family $E_\alpha(z) := \Xi(z) - \mathrm{i}\alpha\Xi'(z)$ ($\alpha > 0$), we prove that the diagonal De Branges reproducing kernel satisfies the exact identity $K_\alpha(t-\mathrm{i}d, t-\mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W}(d,t)$, bypassing the Conrey–Li (2000) global Sonine positivity obstruction for local transverse zero exclusion. We prove the Universal Transverse Square-Cascade $c_{2m-1}(t_0) = \frac{m}{(m!)^2}(\Xi^{(m)}(t_0))^2 > 0$ at any $m$-fold coalescence point and show that across any closely spaced Lehmer pair $t_\pm = t_{\mathrm{mid}} \pm \delta$, $\widetilde{W}(d,t_{\mathrm{mid}}+\tau) = \frac{1}{2}(\Xi''(t_{\mathrm{mid}}))^2(\delta^2 + \tau^2 + d^2) + \mathcal{O}((\delta+\lvert\tau\rvert+d)^3) > 0$ forms a self-healing two-dimensional sum-of-squares bowl. Finally, we prove the cross-term-free Universal Leibniz–Wronskian Identity $\widetilde{W}_{\mathcal{M}}(d,t) \equiv \lvert\mathcal{M}(t-\mathrm{i}d)\rvert^2\widetilde{W}(d,t) + \lvert\Xi(t-\mathrm{i}d)\rvert^2\widetilde{W}_{1,\mathcal{M}}(d,t)$ and the $0\%$-overshoot property $0 \le b_n^{\mathcal{E}}(X) \le 1$ of the entire exponential prime-gauge factor $E_X(s) = \exp(-\sum_{2\le d\le X}\frac{\Lambda(d)}{\ln d}d^{-s})$.

* **Dedicated Repository**: [`thijsmissiaen/paper4-debranges-lehmer-cascade`](https://github.com/thijsmissiaen/paper4-debranges-lehmer-cascade)
* **Package Directory**: [`papers/paper4-debranges-kernels-lehmer-cascade/`](./papers/paper4-debranges-kernels-lehmer-cascade/)
* **Academic PDF Manuscript (18 pages)**: [`paper4_debranges_lehmer_cascade.pdf`](./papers/paper4-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.pdf)
* **AMS-LaTeX Source**: [`paper4_debranges_lehmer_cascade.tex`](./papers/paper4-debranges-kernels-lehmer-cascade/paper/paper4_debranges_lehmer_cascade.tex)
* **Verification Script & Figure Generator**: [`verify_paper4_debranges_lehmer_cascade.py`](./papers/paper4-debranges-kernels-lehmer-cascade/scripts/verify_paper4_debranges_lehmer_cascade.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper4-debranges-kernels-lehmer-cascade/.agents/skills/debranges-kernels-lehmer-cascade/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper4-debranges-kernels-lehmer-cascade/AGENT_SKILL_MANIFEST.json)

## Technical Stack & Formalization

* **Analysis & Computation**: Python (`mpmath`, `numpy`, `scipy`, `sympy`), PyTorch / CUDA for interval Lipschitz certifications.
* **Typesetting & Mathematical Formalization**: LaTeX (`amsart`, `tikz`, `pgfplots`), Tectonic / pdfTeX.

## Academic Dialogue & Critique

I welcome rigorous critique, counter-examples, and academic dialogue from researchers and mathematicians working on analytic number theory, spectral theory, and numerical verification.

Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com) &nbsp;&middot;&nbsp; ORCID: [`https://orcid.org/0009-0007-3678-0810`](https://orcid.org/0009-0007-3678-0810)
