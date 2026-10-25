# Thijs Missiaen

*Independent Researcher in Computational Mathematics & Analytic Number Theory*  
Location: Belgium &nbsp;&middot;&nbsp; Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com)

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

* **Package Directory**: [`papers/paper1-wigner-saddle-coalescence/`](./papers/paper1-wigner-saddle-coalescence/)
* **Academic PDF Manuscript (13 pages)**: [`paper1_wigner_saddle_coalescence.pdf`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.pdf)
* **AMS-LaTeX Source**: [`paper1_wigner_saddle_coalescence.tex`](./papers/paper1-wigner-saddle-coalescence/paper/paper1_wigner_saddle_coalescence.tex)
* **60-Decimal Verification Script & Figure Generator**: [`verify_paper1_wigner_saddle_coalescence.py`](./papers/paper1-wigner-saddle-coalescence/scripts/verify_paper1_wigner_saddle_coalescence.py)
* **AI Agent Skill & Manifest**: [`SKILL.md`](./papers/paper1-wigner-saddle-coalescence/.agents/skills/wigner-saddle-coalescence-xi/SKILL.md) &middot; [`AGENT_SKILL_MANIFEST.json`](./papers/paper1-wigner-saddle-coalescence/AGENT_SKILL_MANIFEST.json)

## Technical Stack & Formalization

* **Analysis & Computation**: Python (`mpmath`, `numpy`, `scipy`, `sympy`), PyTorch / CUDA for interval Lipschitz certifications.
* **Typesetting & Mathematical Formalization**: LaTeX (`amsart`, `tikz`, `pgfplots`), Tectonic / pdfTeX.

## Academic Dialogue & Critique

I welcome rigorous critique, counter-examples, and academic dialogue from researchers and mathematicians working on analytic number theory, spectral theory, and numerical verification.

Contact: [thijs.missiaen@gmail.com](mailto:thijs.missiaen@gmail.com)
