---
name: mirror-hadamard-variational-law
description: Scientific knowledge skill for Paper VII of the Transverse Wronskian series on Riemann's Xi-function. Contains the proven analytical identities for pure-square trigonometric Wronskian multiplets S_pm(d,t) and the Fejer-Lukacs Toeplitz Rayleigh maximum Lambda_max(K) = 2 cos(pi/(K+2)), the mirror-paired Hadamard square-denominator identity K(d, d_rho, u) = d(d^2 - d_rho^2 + u^2)/[((d^2 - d_rho^2 + u^2)^2 + 4u^2 d_rho^2)], monotone critical-strip penetration under finite prime-power sieving, the universal conservation law int_R K(d, d_rho, u) du = pi, the exact quartic difference factorization proving the hypothesis-free minorant K_min(d,u) = K(d, 1/2, u) on |u| >= H_min(d), and the d-independent stationary variational law H_opt(kt) = 2*pi*B_S(kt)/ln(kt/2pi) -> 2*pi(0.112) = 0.70371675... driven by Trudgian's explicit S(T) bound.
---

# The Mirror-Paired Hadamard Square-Denominator Identity, the Universal Conservation Law $\int _ {\mathbb{R}} K \mathrm{d}u \equiv \pi$, and Stationary Trudgian Variation in the Critical Strip (Paper VII)

## 1. Scope and Mathematical Architecture

This skill encapsulates the analytical theorems, exact Toeplitz spectral formulas, partial-fraction difference factorizations, directed-rounding `mpmath.iv` interval-arithmetic certifications, and multi-harmonic variational procedures of **Paper VII** in the series on the Transverse Wronskian of Riemann's $\Xi$-Function (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Paper VII bridges the arithmetic Dirichlet-series representation of the normalized transverse Wronskian $\widehat{W}(d,t) = W(d,t)/\lvert\Xi(t-\mathrm{i} d)\rvert^2 = \mathrm{Re}(\xi'(\frac{1}{2}+d+\mathrm{i} t)/\xi(\frac{1}{2}+d+\mathrm{i} t))$ with its spectral Hadamard partial-fraction expansion across four pillars:

1. **Pure-Square Trigonometric Multiplets and Fejér–Lukács Pole Suppression (Section 2)**: Coupling $\widehat{W}(d,kt)$ ($1 \le k \le K$) with a non-negative trigonometric polynomial $T_K(\theta) = \lvert\sum _ {j=0}^K a_j \mathrm{e}^{\mathrm{i} j\theta}\rvert^2 = c_0 + \sum _ {k=1}^K c_k\cos(k\theta) \ge 0$ eliminates all phase reversals in the Euler-product prime series on $d > \frac{1}{2}$. Over all amplitude vectors $\mathbf{a} \in \mathbb{R}^{K+1} \setminus \lbrace 0\rbrace$, the fundamental ratio $c_1/c_0$ is bounded above by the Toeplitz Rayleigh maximum $\Lambda _ {\max}(K) = 2\cos(\frac{\pi}{K+2})$, attained uniquely by the Fejér–Lukács sine eigenvector $a_j^{\mathrm{FL}} = \sin(\frac{(j+1)\pi}{K+2})$, which achieves the closed-form $\sigma = 1$ pole-suppression ratio $(\sum _ {k=1}^K c_k^{\mathrm{FL}})/c_0^{\mathrm{FL}} = \frac{2}{K+2}\cot^2(\frac{\pi}{2(K+2)}) - 1$.
2. **The Mirror-Paired Hadamard Square-Denominator Identity and Sieved Strip Penetration (Section 3)**: Pairing each nontrivial zero $\rho = \frac{1}{2} + d _ \rho + \mathrm{i}\gamma _ \rho$ with its functional-equation mirror $1 - \bar{\rho} = \frac{1}{2} - d _ \rho + \mathrm{i}\gamma _ \rho$ at the same ordinate $\gamma _ \rho$ yields the mean kernel $K(d, d _ \rho, u) = d(d^2 - d _ \rho^2 + u^2)/D(d, d _ \rho, u)$, whose denominator $D(d, d _ \rho, u) \equiv (d^2 - d _ \rho^2 + u^2)^2 + 4u^2 d _ \rho^2 \equiv (d^2 + d _ \rho^2 + u^2)^2 - 4d^2 d _ \rho^2$ is a sum of two squares and strictly positive for all $u \neq 0$. Subtracting the finite pure-square prime sum $\Delta B_X(d,t)$ at sieve cutoff $X \ge 1$ contracts the certified critical-strip threshold $d _ {\mathrm{strip},X}^{\mathrm{mirror}}(t) < \frac{1}{2}$ monotonically in $X$, reaching $0.000000$ at $t \in \lbrace \gamma_1, \gamma_2, \gamma_3, 60, 100\rbrace$ for $X = 1000$.
3. **The Universal Conservation Law $\int _ {\mathbb{R}} K \mathrm{d}u \equiv \pi$ and Hypothesis-Free Minorant $K _ {\min}(d,u)$ (Section 4)**: For every $d > \lvert d _ \rho\rvert$, the vertical integral of $K(d, d _ \rho, u)$ over $\mathbb{R}$ equals $\pi$ identically, independently of $(d, d _ \rho)$. Moreover, placing $K(d, d _ \rho, u) - K(d, \frac{1}{2}, u)$ over its positive common denominator yields an exact quartic factorization in $u$, proving that for all $d \ge \frac{1}{2}$ the strip-boundary kernel $K _ {\min}(d,u) := K(d, \frac{1}{2}, u) = \frac{d(d^2 - 1/4 + u^2)}{(d^2 - 1/4 + u^2)^2 + u^2}$ is a strictly positive, hypothesis-free minorant $0 < K _ {\min}(d,u) \le K(d, d _ \rho, u)$ for all $\lvert d _ \rho\rvert < \frac{1}{2}$ and all $\lvert u\rvert \ge H _ {\min}(d) := \sqrt{(d^2 - 1/4)/3}$ ($u \neq 0$).
4. **The $d$-Independent Stationary Variational Law $H _ {\mathrm{opt}}(kt) \equiv \frac{2\pi B_S(kt)}{\ln(kt/2\pi)}$ (Section 5)**: Under Stieltjes integration by parts of $\mathrm{d}N(\gamma) = \frac{1}{2\pi}\ln\frac{\gamma}{2\pi} \mathrm{d}\gamma + \mathrm{d}S(\gamma)$ with Trudgian's explicit bound $\lvert S(T)\rvert \le B_S(T)$, the continuous plateau functional $\Phi _ {\mathrm{plat}}(d,H) = \pi\alpha _ {\mathrm{tail}}(d,H) + 2H K _ {\min}(d,H)$ exhibits exact boundary-jump cancellation $\partial_H \Phi _ {\mathrm{plat}}(d,H) \equiv 2H \partial_H K _ {\min}(d,H) < 0$. Consequently, $2 \partial_H K _ {\min}(d,H) \neq 0$ factors out of $\partial_H \mathcal{B} _ {\mathrm{bg}}^{\mathrm{plat}}(d,kt;H) = 0$, pinning the optimal window half-width to the $d$-independent stationary trajectory $H _ {\mathrm{opt}}(kt) \equiv \frac{2\pi B_S(kt)}{\ln(kt/2\pi)} \to 2\pi(0.112) = 0.7037167544\dots$ and annihilating the Stieltjes boundary penalty at the optimum.

## 2. Canonical Identities and Theorems

### 2.1 Pure-Square Trigonometric Multiplets and Fejér–Lukács Toeplitz Bounds (Lemma 2.1 & Theorem 2.2)

For $s = \frac{1}{2} + d + \mathrm{i} t$ ($d > 0$), $\widehat{W}(d,t) = \mathcal{A}(d,t) - \mathcal{Z}(d,t)$, where $\mathcal{A}(d,t) = \frac{1}{2}\mathrm{Re} \psi(\frac{1}{4}+\frac{d}{2}+\mathrm{i}\frac{t}{2}) - \frac{1}{2}\ln\pi + \frac{d+1/2}{(d+1/2)^2+t^2} + \frac{d-1/2}{(d-1/2)^2+t^2}$ and $\mathcal{Z}(d,t) = -\mathrm{Re}(\zeta'(s)/\zeta(s))$, with majorant $M(d) = -\zeta'(\frac{1}{2}+d)/\zeta(\frac{1}{2}+d)$ for $d > \frac{1}{2}$. For $T_K(\theta) = \lvert\sum _ {j=0}^K a_j \mathrm{e}^{\mathrm{i} j\theta}\rvert^2 = c_0 + \sum _ {k=1}^K c_k\cos(k\theta) \ge 0$:

$$\mathcal{S} _ {\pm}(d,t) := c_0 M(d) \pm \sum _ {k=1}^K (\pm 1)^{k-1} c_k \mathcal{Z}(d,kt) \equiv \sum _ {n=2}^{\infty} \frac{\Lambda(n)}{n^{1/2+d}}\left\lvert\sum _ {j=0}^K (\pm 1)^j a_j n^{\mathrm{i} j t}\right\rvert^2 \ge 0,$$

$$\sum _ {k=1}^K c_k \widehat{W}(d,kt) + \mathcal{S} _ +(d,t) \equiv \sum _ {k=1}^K c_k \mathcal{A}(d,kt) + c_0 M(d).$$

The Fejér–Lukács sine eigenvector $a_j^{\mathrm{FL}} = \sin(\frac{(j+1)\pi}{K+2})$ ($0 \le j \le K$) uniquely maximizes $c_1/c_0$:

$$\max _ {\mathbf{a} \neq 0} \frac{c_1(\mathbf{a})}{c_0(\mathbf{a})} = \Lambda _ {\max}(K) = 2\cos\left(\frac{\pi}{K+2}\right), \qquad \frac{1}{c_0^{\mathrm{FL}}}\sum _ {k=1}^K c_k^{\mathrm{FL}} = \frac{2}{K+2} \cot^2\left(\frac{\pi}{2(K+2)}\right) - 1.$$

### 2.2 Mirror-Paired Hadamard Square-Denominator Identity and Sieved Decomposition (Theorems 3.1 & 3.3)

For any zero pair $\lbrace \frac{1}{2} + d _ \rho + \mathrm{i}\gamma _ \rho, \frac{1}{2} - d _ \rho + \mathrm{i}\gamma _ \rho\rbrace$ with $u = t - \gamma _ \rho$:

$$K(d, d _ \rho, u) = \frac{d (d^2 - d _ \rho^2 + u^2)}{D(d, d _ \rho, u)}, \qquad D(d, d _ \rho, u) \equiv \left(d^2 - d _ \rho^2 + u^2\right)^2 + 4u^2 d _ \rho^2 \equiv \left(d^2 + d _ \rho^2 + u^2\right)^2 - 4d^2 d _ \rho^2 > 0 \quad (u \neq 0).$$

Splitting $\mathcal{S} _ +(d,t) = \Delta B_X(d,t) + \mathcal{S} _ {+,>X}(d,t)$ at cutoff $X \ge 1$ yields the effective budget $U _ {T,X}(d,t) = \sum _ {k=1}^K c_k\mathcal{A}(d,kt) + c_0 M(d) - \Delta B_X(d,t)$ and the monotone critical-strip bound

$$\lvert d _ {\rho_0}\rvert \le d _ {\mathrm{strip},X}^{\mathrm{mirror}}(t) := \inf _ {d > 1/2} \sqrt{\max\left(0,  d^2 - \frac{2c_1 d}{U _ {T,X}(d,t)}\right)} \le d _ {\mathrm{strip},1}^{\mathrm{mirror}}(t) < \frac{1}{2}.$$

### 2.3 Universal $\pi$-Conservation Law and Hypothesis-Free Minorant $K _ {\min}(d,u)$ (Theorems 4.1 & 4.2)

For all $d > \lvert d _ \rho\rvert$ and $H \ge 0$:

$$\int _ {-H}^H K(d, d _ \rho, u) \mathrm{d}u = \arctan\left(\frac{H}{d - d _ \rho}\right) + \arctan\left(\frac{H}{d + d _ \rho}\right), \qquad \int _ {-\infty}^{\infty} K(d, d _ \rho, u) \mathrm{d}u \equiv \pi.$$

With $W := d^2 + u^2$ and $K _ {\min}(d,u) := K(d, \frac{1}{2}, u) = \frac{d(d^2 - 1/4 + u^2)}{(d^2 - 1/4 + u^2)^2 + u^2}$, the exact difference factorization holds:

$$K(d, d _ \rho, u) - K _ {\min}(d,u) \equiv \frac{d \left(\frac{1}{4} - d _ \rho^2\right)\left[3u^4 + (2d^2 + d _ \rho^2 + \frac{1}{4})u^2 - (d^2 - \frac{1}{4})(d^2 - d _ \rho^2)\right]}{\left((d^2 - d _ \rho^2 + u^2)^2 + 4u^2 d _ \rho^2\right)\left((d^2 - \frac{1}{4} + u^2)^2 + u^2\right)} \ge 0$$

for all $d \ge \frac{1}{2}$, all $d _ \rho \in (-\frac{1}{2}, \frac{1}{2})$, and all $\lvert u\rvert \ge H _ {\min}(d) := \sqrt{(d^2 - 1/4)/3}$ ($u \neq 0$).

### 2.4 Stationary Variational Law $H _ {\mathrm{opt}}(kt)$ and Multi-Harmonic Strip Contraction (Theorems 5.1 & 5.2)

With $\alpha _ {\mathrm{tail}}(d,H) = 1 - \frac{1}{\pi}[\arctan\frac{H}{d-1/2} + \arctan\frac{H}{d+1/2}]$ and $\Phi _ {\mathrm{plat}}(d,H) = \pi\alpha _ {\mathrm{tail}}(d,H) + 2H K _ {\min}(d,H)$:

$$\frac{\partial\Phi _ {\mathrm{plat}}}{\partial H}(d,H) \equiv 2H \frac{\partial K _ {\min}}{\partial H}(d,H) < 0, \qquad \frac{\partial\mathcal{B} _ {\mathrm{bg}}^{\mathrm{plat}}}{\partial H}(d,\tau;H) \equiv 2 \frac{\partial K _ {\min}}{\partial H}(d,H)\left[\frac{H}{2\pi}\ln\left(\frac{\tau}{2\pi}\right) - B_S(\tau)\right],$$

yielding the unique $d$-independent stationary window and penalty-free optimum:

$$H _ {\mathrm{opt}}(\tau) \equiv \frac{2\pi B_S(\tau)}{\ln(\tau/2\pi)} \underset{\tau\to\infty}{\longrightarrow} 2\pi A_S = 0.70371675440411368541\dots, \qquad \mathcal{B} _ {\mathrm{bg}}^{\mathrm{opt}}(d,\tau) \equiv \frac{1}{2} \alpha _ {\mathrm{tail}}\bigl(d, H _ {\mathrm{opt}}(\tau)\bigr)\ln\left(\frac{\tau}{2\pi}\right).$$

## 3. Reproducible Verification Workflow

Run the standalone verification and figure-generation suite from the repository root:

```bash
python scripts/verify_paper7_mirror_hadamard_variation.py
```
