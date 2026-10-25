---
name: chladni-portal-dilogarithm-synthesis
description: Scientific knowledge skill for Flagship Paper XIV containing the exact analytical identities and 65-decimal certified benchmarks for the Chladni Portal-Frame Synthesis: Macro-Scale Codimension-Three Saddle Overdetermination and Case-I Davenport-Heilbronn Forcing, Meso-Scale 50%-50% Green-Dirichlet-Neumann Chamber Equipartition and Wronskian-Slope Reciprocity, and Micro-Scale Euler Dilogarithm Li2(2*theta) Horizon Damping and Dirichlet Hyperbola Reflection Duality.
---

# The Chladni Portal-Frame & Dilogarithm Horizon Synthesis Skill (Paper XIV)

Activate this skill when working with **Paper XIV** ("The Chladni Portal-Frame Synthesis: Codimension-Three Saddle Overdetermination, Green–Dirichlet–Neumann Chamber Equipartition, and Euler Dilogarithm $\mathrm{Li} _ 2(2\theta)$ Horizon Damping for Riemann's $\Xi$ Function").

## 1. Normalization & Coordinate Dictionary

- **Spectral coordinates:** $s = \sigma + \mathrm{i}t = \frac{1}{2} + d + \mathrm{i}t$, rotated coordinate $z = -\mathrm{i}(s - \frac{1}{2}) = t - \mathrm{i}d$ ($d \ge 0$).
- **Harmonic decomposition:** $\Xi(t - \mathrm{i}d) = R(d,t) + \mathrm{i}I(d,t)$, where $R(-d,t) = R(d,t)$ (even) and $I(-d,t) = -I(d,t)$ (odd, $I(0,t) \equiv 0$).
- **Critical-line Laguerre invariant:** $c_1(t) := (\Xi'(t))^2 - \Xi(t)\Xi''(t)$.
- **Reduced transverse Wronskian:**

$$\widetilde{W}(d,t) := \frac{R(d,t) \partial_t I(d,t) - I(d,t) \partial_t R(d,t)}{d} = \frac{R(d,t) \partial_d R(d,t) + I(d,t) \partial_d I(d,t)}{d}, \qquad \lim _ {d\to 0^+} \widetilde{W}(d,t) = c_1(t).$$

- **Full transverse Wronskian:** $W(d,t) = 2(R \partial_t I - I \partial_t R) \equiv 2d \widetilde{W}(d,t)$.

---

## 2. Macro-Scale Identities (Section 2)

1. **Tangential Cauchy–Riemann Monotonicity Law (Theorem 2.1(a)):** Along any smooth Chladni nodal arc $\Gamma \subset \lbrace I(d,t)=0\rbrace$ parameterized by arc length $s$,

$$\frac{\mathrm{d}}{\mathrm{d}s}R\bigl(d(s),t(s)\bigr) = \pm \bigl\lvert \nabla I\bigl(d(s),t(s)\bigr)\bigr\rvert = \pm \bigl\lvert \Xi'\bigl(t(s)-\mathrm{i}d(s)\bigr)\bigr\rvert, \qquad \frac{\partial R}{\partial n}(\Gamma) = 0.$$

2. **Harmonic Curvature Conservation Law (Theorem 2.1(c)):** At any simple critical-line extremum $\tau _ {\star}$ ($\Xi'(\tau _ {\star})=0$, $\Xi(\tau _ {\star})\ne 0$),

$$\frac{1}{2}\frac{\mathrm{d}^2}{\mathrm{d}d^2}\Bigl[\bigl\lvert R\bigl(d,\tau _ {\star}(d)\bigr)\bigr\rvert^2\Bigr] _ {d=0} = -\frac{1}{2}\frac{\mathrm{d}^2}{\mathrm{d}t^2}\Bigl[\Xi(t)^2\Bigr] _ {t=\tau _ {\star}} = -\Xi(\tau _ {\star})\Xi''(\tau _ {\star}) = c_1(\tau _ {\star}),$$

   with Taylor expansion $\lvert R(d,\tau _ {\star}(d))\rvert^2 = \Xi(\tau _ {\star})^2 + c_1(\tau _ {\star})d^2 + O(d^4)$.
3. **Wall-Zero Dichotomy & Codimension-Three Saddle Overdetermination (Theorem 2.1(d) & Theorem 2.2):** An off-critical zero $z _ {\star} = \gamma _ {\star} - \mathrm{i}d _ {\star}$ ($d _ {\star} > 0$) on a Chladni wall $\Gamma _ {\star}$ anchored at $(0,\tau _ {\star})$ requires either:
   - **Case I (1D Critical-Line Laguerre Defect):** $c_1(\tau _ {\star}) = -\Xi(\tau _ {\star})\Xi''(\tau _ {\star}) \le 0$ on $d=0$. In the Davenport–Heilbronn counterexample $\Xi _ {\mathrm{DH}}(z)$, **100% of the first 5 off-critical zeros** in $t \in [10,260]$ are forced into Case I ($S_m^2 c_1^{\mathrm{DH}}(\tau _ {\mathrm{false},m}) \in [-873.26, -34.21] < 0$).
   - **Case II (Interior Codimension-3 Saddle Bifurcation):** $c_1(\tau _ {\star}) > 0$, requiring an interior complex critical point $z_c = t_c - \mathrm{i}d_c$ ($d_c > 0$) on the wall satisfying the **3-equation system in 2 real unknowns** $I(d_c,t_c) = \partial_d I(d_c,t_c) = \partial_t I(d_c,t_c) = 0$.

---

## 3. Meso-Scale Identities (Section 3)

1. **Exact 50%–50% Green–Dirichlet–Neumann Chamber Equipartition (Theorem 3.1):** On every Chladni chamber $\Omega_k(d_0) = \lbrace (d,t) : 0 < d < d_0, \tau_k(d) < t < \tau _ {k+1}(d)\rbrace$, the curved side walls $\Gamma_k,\Gamma _ {k+1}$ act simultaneously as homogeneous Dirichlet boundaries ($I=0$) for $I$ and homogeneous Neumann boundaries ($\partial_n R = 0$) for $R$, enforcing the exact identity:

$$\int _ {\tau_k(d_0)}^{\tau _ {k+1}(d_0)} R(d_0,t) \partial_d R(d_0,t) \mathrm{d}t \equiv \int _ {\tau_k(d_0)}^{\tau _ {k+1}(d_0)} I(d_0,t) \partial_d I(d_0,t) \mathrm{d}t \equiv \iint _ {\Omega_k(d_0)} \bigl\lvert \Xi'(t - \mathrm{i}d)\bigr\rvert^2 \mathrm{d}d \mathrm{d}t = \frac{1}{2}\mathcal{W} _ {\Omega_k}(d_0) > 0.$$

2. **Wronskian–Slope Reciprocity & Wall-Amplitude Evolution (Theorem 3.2):** Along $\Gamma_k = \lbrace (d,\tau_k(d))\rbrace$,

$$d \widetilde{W}\bigl(d,\tau_k(d)\bigr) \equiv R\bigl(d,\tau_k(d)\bigr) \partial_t I\bigl(d,\tau_k(d)\bigr) \equiv R\bigl(d,\tau_k(d)\bigr) \partial_d R\bigl(d,\tau_k(d)\bigr), \qquad \tau_k'(d) \equiv \frac{R\bigl(d,\tau_k(d)\bigr) \partial_t R\bigl(d,\tau_k(d)\bigr)}{d \widetilde{W}\bigl(d,\tau_k(d)\bigr)}, \qquad \frac{\mathrm{d}}{\mathrm{d}d}\Bigl[\frac{1}{2}R\bigl(d,\tau_k(d)\bigr)^2\Bigr] \equiv d \widetilde{W} \bigl(1 + \tau_k'(d)^2\bigr).$$

   Along any real nodal contour $R(d,t_R(d))=0$, $t_R'(d) \equiv I \partial_t I/(d \widetilde{W})$. At the detached real parabolic fold tip of $\Xi _ {\mathrm{DH}}(z)$ at $(d _ {\mathrm{tip}},t _ {\mathrm{tip}}) = (0.2980980120, 85.6254944428)$, $\widetilde{W} _ {\mathrm{DH}}(d _ {\mathrm{tip}},t _ {\mathrm{tip}}) = 0$ and $\lvert t_R'(d _ {\mathrm{tip}})\rvert = \infty$.

## 4. Micro-Scale Identities (Section 4)

1. **Attenuated Quadratic Semiprime Spring (Theorem 4.1):** For $\mathcal{Q} _ \eta(U) = 1 + U + \frac{\eta}{2}U^2$ with capture fraction $\eta \in (1/2, 1]$, the discriminant is $\Delta(\eta) = 1 - 2\eta < 0$ and the anti-phase floor is $\lvert \mathcal{Q} _ \eta(U)\rvert^2 \ge (1 - \frac{1}{2\eta})^2 > 0$.
2. **Exact Closed-Form Euler Dilogarithm $\mathrm{Li} _ 2(2\theta)$ Capture Laws (Theorem 4.2):** With $\tau = 2\theta \in (0,1/2)$, $I_p(\theta) = \int _ \theta^{1/2} (1-2\alpha)^p \alpha^{-1} \mathrm{d}\alpha$, and $\mathcal{M} _ {\mathrm{sq}}^{(p)}(\theta) = (I_p(\theta))^2$, the forward-triangle integrals $\mathcal{M} _ {\mathrm{tri}}^{(p)}(\theta)$ for $p \in \lbrace 0,1,2,3\rbrace$ are:
   - $\mathcal{M} _ {\mathrm{tri}}^{(0)}(\theta) = \frac{\pi^2}{6} - 2 \mathrm{Li} _ 2(\tau) + \ln^2\tau - 2(\ln\tau)\ln(1-\tau)$
   - $\mathcal{M} _ {\mathrm{tri}}^{(1)}(\theta) = \mathcal{M} _ {\mathrm{tri}}^{(0)}(\theta) - 2(1-\tau)\ln(\frac{1-\tau}{\tau}) + (1-2\tau)$
   - $\mathcal{M} _ {\mathrm{tri}}^{(2)}(\theta) = \mathcal{M} _ {\mathrm{tri}}^{(0)}(\theta) - (1-\tau)(3-\tau)\ln(\frac{1-\tau}{\tau}) + (\frac{17}{6} - \frac{13}{6}\tau + \frac{1}{3}\tau^2)(1-2\tau)$
   - $\mathcal{M} _ {\mathrm{tri}}^{(3)}(\theta) = \mathcal{M} _ {\mathrm{tri}}^{(0)}(\theta) - \frac{1}{3}(1-\tau)(11 - 7\tau + 2\tau^2)\ln(\frac{1-\tau}{\tau}) + (\frac{1801}{360} - \frac{1199}{180}\tau + \frac{218}{45}\tau^2 - \frac{89}{45}\tau^3 + \frac{49}{90}\tau^4 - \frac{1}{9}\tau^5)(1-2\tau)$
3. **Canonical Values at $\theta = 1/6$ ($\tau = 1/3$):**
   - $\eta^{(0)}(1/6) = 0.2439550996 < 1/2$ ($\Delta = +0.512090 > 0$, under-damped)
   - $\eta^{(1)}(1/6) = 0.4955987356 < 1/2$ ($\Delta = +0.008803 > 0$, near-critical)
   - $\eta^{(2)}(1/6) = 0.6768906167 > 2/3 > 1/2$ ($\Delta = -0.353781 < 0$, spring floor $+0.0682917228 > 0$, $\theta _ {*,50\text{ pct}}^{(2)} = 0.1905498467 > \theta _ {\mathrm{W\text{-}pure}}^\infty = 0.1839448708 > \frac{1}{2\mathrm{e}} > \frac{1}{6}$)
   - $\eta^{(3)}(1/6) = 0.7980906218 > 3/4$ ($\Delta = -0.596181 < 0$, spring floor $+0.1395067855 > 0$, $\theta _ {*,50\text{ pct}}^{(3)} = 0.2039297718 > 1/5$)
4. **Dirichlet Hyperbola Reflection Duality & 4th-Order Apex Pinching (Theorem 4.3):** For $m _ {\mathrm{dual}}(pq) = N _ \star(t)^2/(pq)$, $\omega_p(t) + \omega_q(t) \equiv \ln m _ {\mathrm{dual}}(pq)$ and $(\omega_p\omega_q)^2 \le \frac{1}{16}(\ln m _ {\mathrm{dual}}(pq))^4 \to 0$ as $pq \to N _ \star(t)^2$.

## 5. Verification Command

```bash
python scripts/verify_paper14_portal_dilogarithm_synthesis.py
```
