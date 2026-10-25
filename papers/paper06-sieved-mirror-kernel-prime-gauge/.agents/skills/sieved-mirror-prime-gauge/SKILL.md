---
name: sieved-mirror-prime-gauge
description: Scientific knowledge skill for Paper VI of the Transverse Wronskian series on Riemann's Xi-function. Contains the proven analytical identities for the sieved mirror kernel Phi_X(u) in |Im(u)| < pi/8, the 3x3 Schur-complement optimal mollifier filter P_opt(x), the zero-free exponential prime-gauge factor E_X(z) with bijective zero preservation and quadratic prime-power purity on (X, X^2], the exact negative hyperbolic mirror weight w_n*(d_0) = -(n^{2d_0}-1)^{-1} without waterbed tail, power-law L^2 diagonal damping X^{-2d}/(2d ln X) with universal Sub-Selberg closure theta_E^{sq}(d,t) < 1/2, the exact finite-divisor identity on [1, N(t)] with composite-gap Gaussian tail suppression, and the shifted self-mirrored prime-gauge factor E_{X,kappa}(z).
---

# The Sieved Mirror Kernel $\Phi_X(u)$, the Zero-Free Exponential Prime-Gauge Factor $\mathcal{E} _ X(z)$, and Quadratic Prime-Power Purity on $(X, X^2]$ (Paper VI)

## 1. Scope and Mathematical Architecture

This skill encapsulates the analytical theorems, exact arithmetic recursions, directed-rounding `mpmath.iv` interval-arithmetic certifications, and GPU continuum verification procedures of **Paper VI** in the series on the Transverse Wronskian of Riemann's $\Xi$-Function (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Paper VI couples the universal mother-wavelet decomposition $\Phi(u) = \sum _ {m=1}^\infty m^{-1/2}\varphi_1(u + \frac{1}{2}\ln m)$ established in Paper V with finite-divisor Euler-product sieving across four analytical pillars:

1. **The Sieved Mirror Kernel $\Phi_X(u)$ and $3\times 3$ Schur-Complement Optimization $P _ {\mathrm{opt}}(x)$ (Section 2)**: Multiplying $\Xi(z)$ by the holomorphic self-mirrored Dirichlet polynomial $\mathcal{M} _ X(z) := M_X(\frac{1}{2}+\mathrm{i}z)M_X(\frac{1}{2}-\mathrm{i}z)$ preserves all zeros of $\Xi(z)$, maintains exact critical-line mirror symmetry, and corresponds in the logarithmic $u$-domain to a finite discrete convolution kernel $\Phi_X(u)$ that is holomorphic throughout the full strip $\lvert\mathrm{Im}(u)\rvert < \frac{\pi}{8}$. Within the cubic Conrey–Selberg filter family $P(x) = x + a_2 x(1-x) + a_3 x(1-x)^2$, the unique $3\times 3$ Schur-complement minimizer $(a_2^{\ast}, a_3^{\ast}) = (+0.81746499, +2.13817017)$ reduces quadratic cross-term energy by a factor of $31.4\times$ relative to Conrey's cubic filter $P_3(x)$.
2. **The Zero-Free Exponential Prime-Gauge Factor $\mathcal{E} _ X(z)$, Quadratic Prime-Power Purity on $(X, X^2]$, and Negative Hyperbolic Mirror Annihilation $w_n^{\ast}(d_0)$ (Section 3)**: Replacing the polynomial mollifier $M_X(s)$ by the exponential prime-gauge factor $E_X(s) := \exp(-\sum _ {n=2}^{\lfloor X\rfloor} \frac{\Lambda(n)}{\ln n}n^{-s})$ and $\mathcal{E} _ X(z) := E_X(\frac{1}{2}+\mathrm{i}z)E_X(\frac{1}{2}-\mathrm{i}z)$ eliminates all parasitic mollifier zeros and poles in $\mathbb{C}$ ($\mathcal{Z}(\Xi _ {\mathcal{E} _ X}) = \mathcal{Z}(\Xi)$ bijectively). In $\zeta(s)E_X(s) = \sum _ {N=1}^\infty b_N^{\mathcal{E}}(X)N^{-s}$, we prove $b_1^{\mathcal{E}}(X) = 1$, $b_N^{\mathcal{E}}(X) = 0$ for $2 \le N \le \lfloor X\rfloor$, universal bounds $0 \le b_N^{\mathcal{E}}(X) \le 1$, and **Quadratic Prime-Power Purity** $b_N^{\mathcal{E}}(X) \equiv \frac{\Lambda(N)}{\ln N}$ across the entire quadratic window $\lfloor X\rfloor < N < (\lfloor X\rfloor + 1)^2$. Furthermore, in the weighted self-mirrored gauge $\mathcal{E} _ {X,\mathbf{w}}(z)$, the **negative hyperbolic mirror weight** $w_n^{\ast}(d_0) := -(n^{2d_0}-1)^{-1} < 0$ annihilates the effective prime-power amplitude $C_n(d_0, w_n^{\ast}(d_0)) \equiv 0$ on $d = d_0$ while introducing zero terms for $n > X$ (0% waterbed tail).
3. **Power-Law $L^2$ Damping $\frac{X^{-2d}}{2d\ln X}$ and Universal Sub-Selberg Closure $\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) < \frac{1}{2}$ (Section 4)**: Whenever $X \in [\lfloor\sqrt{N(t)}\rfloor, N(t)]$ ($\theta \in [\frac{1}{4}, \frac{1}{2}]$), the entire Riemann–Siegel main sum $1 \le n \le N(t) = \lfloor\sqrt{t/2\pi}\rfloor$ lies inside the quadratic purity window $(X, (\lfloor X\rfloor+1)^2)$, replacing logarithmic $(\ln X)^{-2r}$ decay by the power-law decay $\mathcal{E} _ {\mathrm{diag}}^{\mathcal{E},\mathrm{fwd}}(d,t;X) < \frac{X^{-2d}}{2d\ln X} \to 0$ and guaranteeing $X _ {\mathrm{krit}}^{\mathcal{E},\mathrm{sq}}(d,t) < N(t)$ ($\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) < \frac{1}{2}$) across all $d \in (0, \frac{1}{2}]$ and $t \ge 8\pi$ ($N(t) \ge 2$, with $\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) = 0 < \frac{1}{2}$ on $2\pi < t < 8\pi$).
4. **Exact Finite-Divisor Identity on $[1, N(t)]$, Composite-Gap Gaussian Tail Barrier, and Shifted Self-Mirrored Gauge $\mathcal{E} _ {X,\kappa}(z)$ (Section 5)**: The finite truncation $M _ {X,N(t)}^{\mathcal{E}}(s) = \sum _ {k=1}^{N(t)} c_k^{\mathcal{E}}(X)k^{-s}$ generates the exact same coefficients $b_n^{\mathcal{E}}(X)$ on $1 \le n \le N(t)$ as the infinite series $E_X(s)$, while composite tail gaps deepen the Gaussian saddle-point barrier by $10^{-64{,}798}$ at $t = 10^8$. Finally, the shifted gauge $\mathcal{E} _ {X,\kappa}(z) := E_X(\frac{1}{2}+\kappa+\mathrm{i}z)E_X(\frac{1}{2}+\kappa-\mathrm{i}z)$ satisfies an exact outer factorization and eliminates conjugate $n^{+d}$ growth whenever $\kappa \ge d$.

## 2. Canonical Identities and Theorems

### 2.1 Sieved Mirror Kernel $\Phi_X(u)$ and Mother-Wavelet Shell Annihilation (Lemma 2.1, Theorem 2.2 & Proposition 2.3)

For a real filter $P : [0,1] \to \mathbb{R}$ with $P(0) = 0$ and $P(1) = 1$, define $M_X(s) := \sum _ {k=1}^{\lfloor X\rfloor} c_k k^{-s}$ with $c_k := \mu(k)P(\frac{\ln(X/k)}{\ln X})$ and the self-mirrored multiplier:

$$\mathcal{M} _ X(z) := M_X\left(\frac{1}{2} + \mathrm{i}z\right)M_X\left(\frac{1}{2} - \mathrm{i}z\right) = \sum _ {j,k \le X} \frac{c_j c_k}{\sqrt{jk}}\cos\left(z\ln\frac{j}{k}\right), \qquad \Xi_X(z) := \Xi(z) \mathcal{M} _ X(z).$$

In the canonical Mellin–Jacobi normalization $\Xi(z) = 8\int_0^\infty \Phi(u)\cos(2zu) \mathrm{d}u$ (equivalently $\Xi(z) = 2\int_0^\infty \widetilde{\Phi}(\tilde{u})\cos(\tilde{u}z) \mathrm{d}\tilde{u}$ with $\widetilde{\Phi}(\tilde{u}) = 2\Phi(\tilde{u}/2){}$), $\Xi_X(z) = 8\int_0^\infty \Phi_X(u)\cos(2zu) \mathrm{d}u$, where the **Sieved Mirror Kernel**:

$$\Phi_X(u) := \sum _ {j,k \le X} \frac{c_j c_k}{2\sqrt{jk}}\left[\Phi\left(u + \frac{1}{2}\ln\frac{j}{k}\right) + \Phi\left(u - \frac{1}{2}\ln\frac{j}{k}\right)\right]$$

is even ($\Phi_X(-u) = \Phi_X(u){}$) and holomorphic throughout $\lvert\mathrm{Im}(u)\rvert < \frac{\pi}{8}$. Under the forward shift operator, the Mother Wavelet $\varphi_1(w)$ satisfies:

$$\mathcal{T} _ X[\Phi] (u) := \sum _ {k \le X} \frac{c_k}{\sqrt{k}} \Phi\left(u + \frac{1}{2}\ln k\right) = \sum _ {N=1}^{\infty} \frac{b_N}{\sqrt{N}} \varphi_1\left(u + \frac{1}{2}\ln N\right), \qquad b_N := \sum _ {k \mid N,   k \le X} c_k, \qquad b_1 \equiv 1.$$

For the cubic family $P(x) = x + a_2 x(1-x) + a_3 x(1-x)^2$ ($\mathbf{a} = (1, a_2, a_3)^T$), the $2\times 2$ Schur-complement minimizer of $\mathcal{E} _ {\mathrm{quad}}(\mathbf{a}) = \mathbf{a}^T\mathbf{Q}\mathbf{a}$ on $t \in [14, 300]$ ($d = 0.20, X = 50$) is:

$$a_2^{\ast} = +0.81746499, \qquad a_3^{\ast} = +2.13817017, \qquad \mathcal{E} _ {\mathrm{quad}}(\mathbf{a}^{\ast}) = 2.4843 \times 10^{-4}    (31.4\times \text{ below } P_3).$$

### 2.2 Zero-Free Exponential Prime-Gauge Factor $\mathcal{E} _ X(z)$ and Quadratic Purity (Theorems 3.1 & 3.2)

Define the entire, zero-free exponential prime-gauge factor and its self-mirrored counterpart:

$$E_X(s) := \exp\left(-\sum _ {n=2}^{\lfloor X\rfloor} \frac{\Lambda(n)}{\ln n} n^{-s}\right) = \sum _ {k=1}^{\infty} c_k^{\mathcal{E}}(X) k^{-s}, \qquad \mathcal{E} _ X(z) := E_X\left(\frac{1}{2} + \mathrm{i}z\right)E_X\left(\frac{1}{2} - \mathrm{i}z\right).$$

In the sieved product $\zeta(s)E_X(s) = \exp(\sum _ {n > \lfloor X\rfloor} \frac{\Lambda(n)}{\ln n}n^{-s}) = \sum _ {N=1}^\infty b_N^{\mathcal{E}}(X)N^{-s}$:
1. $c_k^{\mathcal{E}}(X) \equiv \mu(k)$ for $1 \le k \le \lfloor X\rfloor$, $b_1^{\mathcal{E}}(X) = 1$, and $b_N^{\mathcal{E}}(X) = 0$ for $2 \le N \le \lfloor X\rfloor$.
2. $b _ {MN}^{\mathcal{E}}(X) = b_M^{\mathcal{E}}(X)b_N^{\mathcal{E}}(X)$ for $\gcd(M,N) = 1$, and $0 \le b_N^{\mathcal{E}}(X_2) \le b_N^{\mathcal{E}}(X_1) \le 1$ for $X_2 \ge X_1 \ge 1$.
3. **Quadratic Prime-Power Purity**: For all $\lfloor X\rfloor < N < (\lfloor X\rfloor + 1)^2$:

$$b_N^{\mathcal{E}}(X) \equiv \frac{\Lambda(N)}{\ln N} = \begin{cases} 1 & \text{if } N = p \text{ is prime}, \\ 1/k & \text{if } N = p^k \text{ is a prime power } (k \ge 2), \\ 0 & \text{if } N \text{ is composite with } \ge 2 \text{ distinct prime factors}. \end{cases}$$

In the weighted gauge $\mathcal{E} _ {X,\mathbf{w}}(z)$ with $E _ {X,\mathbf{w}}(s) := \exp(-\sum _ {n=2}^{\lfloor X\rfloor} w_n\frac{\Lambda(n)}{\ln n}n^{-s})$, the effective prime-power amplitude in $\frac{W _ {\mathcal{E} _ {X,\mathbf{w}}}(d,t)}{\lvert\Xi _ {\mathcal{E} _ {X,\mathbf{w}}}(t-\mathrm{i}d)\rvert^2}$ is $C_n(d, w_n) = (1 - w_n)n^{-(1/2+d)} + w_n n^{-(1/2-d)} = n^{-(1/2+d)} + 2w_n n^{-1/2}\sinh(d\ln n)$. Choosing the **negative hyperbolic mirror weight**:

$$w_n^{\ast}(d_0) := -\frac{n^{-d_0}}{n^{+d_0} - n^{-d_0}} = -\frac{1}{n^{2d_0} - 1} < 0 \qquad (\forall  d_0 > 0,   2 \le n \le \lfloor X\rfloor)$$

yields $C_n(d_0, w_n^{\ast}(d_0)) \equiv 0$ identically on $d = d_0$ with zero waterbed tail for $n > X$.

### 2.3 Power-Law $L^2$ Damping and Universal Sub-Selberg Closure (Theorem 4.1)

For any $X \in [\lfloor\sqrt{N(t)}\rfloor, N(t)]$ ($\theta = \frac{\ln X}{\ln(t/2\pi)} \in [\frac{1}{4}, \frac{1}{2}]$), the entire Riemann–Siegel main sum $1 \le n \le N(t) = \lfloor\sqrt{t/2\pi}\rfloor$ lies inside $(X, (\lfloor X\rfloor+1)^2)$, so the forward $L^2$ diagonal decays as a power law:

$$\mathcal{E} _ {\mathrm{diag}}^{\mathcal{E},\mathrm{fwd}}(d,t;X) = \sum _ {\lfloor X\rfloor < n \le N(t)} \left(\frac{\Lambda(n)}{\ln n}\right)^2 n^{-(1+2d)} \sim \frac{X^{-2d} - N(t)^{-2d}}{2d \ln X} < \frac{X^{-2d}}{2d \ln X} \longrightarrow 0,$$

and monotonicity together with $Y^{\mathcal{E}}(d,t;N(t)) = 0$ ensures universal Sub-Selberg closure $X _ {\mathrm{krit}}^{\mathcal{E},\mathrm{sq}}(d,t) < N(t) \iff \theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) < \frac{1}{2}$ for all $d \in (0, \frac{1}{2}]$ and all $t \ge 8\pi$ ($N(t) \ge 2$, with $\theta _ {\mathcal{E}}^{\mathrm{sq}}(d,t) = 0 < \frac{1}{2}$ on $2\pi < t < 8\pi$).

### 2.4 Finite-Divisor Identity, Composite-Gap Tail, and Shifted Gauge $\mathcal{E} _ {X,\kappa}(z)$ (Theorems 5.1 & 5.2)

Because only divisors $1 \le k \le n \le N(t)$ contribute to $b_n^{\mathcal{E}}(X) = \sum _ {k\mid n} c_k^{\mathcal{E}}(X)$, the finite truncation $M _ {X,N(t)}^{\mathcal{E}}(s) := \sum _ {k=1}^{N(t)} c_k^{\mathcal{E}}(X)k^{-s}$ satisfies $b_n^{\mathcal{E},\mathrm{trunc}}(X;N(t)) \equiv b_n^{\mathcal{E}}(X)$ for all $1 \le n \le N(t)$. For $n \ge N(t)+1$, $0 \le b_n^{\mathcal{E}}(X) \le 1$ prevents any tail inflation, and composite gaps $N(t)+1 \le n < q _ {\mathrm{next}}(X,t)$ suppress the Gaussian tail barrier by $\exp(-\pi(q _ {\mathrm{next}}^2 - (N(t)+1)^2)\cos 4\eta^{\ast})$ ($10^{-64{,}798}$ at $t = 10^8$, where $3990 \dots 4000$ vanish identically up to $p _ {\mathrm{next}} = 4001$).

Finally, for any shift $\kappa \ge 0$, the shifted self-mirrored gauge $\mathcal{E} _ {X,\kappa}(z) := E_X(\frac{1}{2}+\kappa+\mathrm{i}z)E_X(\frac{1}{2}+\kappa-\mathrm{i}z)$ satisfies the exact outer factorization:

$$\left\lvert\Xi _ {\mathcal{E} _ {X,\kappa}}(t - \mathrm{i}d)\right\rvert^2 \equiv \left\lvert E_X\left(\frac{1}{2} + \kappa - d - \mathrm{i}t\right)\right\rvert^2 \cdot \left\lvert\Xi(t - \mathrm{i}d) E_X\left(\frac{1}{2} + \kappa + d + \mathrm{i}t\right)\right\rvert^2,$$

suppressing conjugate $n^{+d}$ growth whenever $\kappa \ge d$, while the integrated Gabcke–Olver remainder obeys $\frac{2}{d^2}\int_0^d u \mathrm{sinc} _ {\mathrm{hyp}}(uL/2) \mathrm{d}u \equiv \mathrm{sinc} _ {\mathrm{hyp}}^2(dL/4)$ and $\mathrm{sinc} _ {\mathrm{hyp}}^2(dL/4)/\mathrm{sinc} _ {\mathrm{hyp}}^2(dL/2) \equiv \cosh^{-2}(dL/4)$.

## 3. Reproducible Verification Workflow

Run the standalone verification and figure-generation suite from the repository root:

```bash
python scripts/verify_paper6_sieved_mirror_prime_gauge.py
```
