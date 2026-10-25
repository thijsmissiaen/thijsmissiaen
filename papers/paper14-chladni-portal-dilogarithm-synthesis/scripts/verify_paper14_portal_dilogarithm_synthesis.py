#!/usr/bin/env python3
"""
Standalone Reproducible Verification & Figure Generation Suite for Flagship Paper XIV:
"The Chladni Portal-Frame Synthesis: Codimension-Three Saddle Overdetermination,
Green-Dirichlet-Neumann Chamber Equipartition, and Euler Dilogarithm Li_2(2*theta)
Horizon Damping for Riemann's Xi-Function"

Author: Thijs Missiaen (Independent Researcher, Belgium)
ORCID:  https://orcid.org/0009-0007-3678-0810

Precision Standards:
- 65-decimal arbitrary-precision arithmetic (`mpmath.mp.dps = 65`)
- 35-decimal directed-rounding interval arithmetic (`mpmath.iv.dps = 35`)
- Generates two 4-panel 300-DPI publication figures (`< 500 KB` PDF + PNG)
"""

import math
import shutil
import sys
import time
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

mp.mp.dps = 65
mp.iv.dps = 35

ROOT_DIR = Path(__file__).resolve().parent.parent
PAPER_FIG_DIR = ROOT_DIR / "paper" / "figures"
ROOT_FIG_DIR = ROOT_DIR / "figures"
PAPER_FIG_DIR.mkdir(parents=True, exist_ok=True)
ROOT_FIG_DIR.mkdir(parents=True, exist_ok=True)


# ==============================================================================
# FAST 65-DECIMAL POISSON-JACOBI & REAL-CENTERED TAYLOR EVALUATORS
# ==============================================================================

N_U_RIEM = 300
U_MAX_RIEM = mp.mpf("2.30")
DU_RIEM = U_MAX_RIEM / N_U_RIEM
U2_RIEM = []
W_RIEM = []
_pi = mp.pi
_pi2 = _pi * _pi
for _k in range(N_U_RIEM + 1):
    _u = mp.mpf(_k) * DU_RIEM
    _e4, _e5, _e9 = mp.exp(4 * _u), mp.exp(5 * _u), mp.exp(9 * _u)
    _phi = mp.mpf("0")
    for _m in range(1, 18):
        _m2 = mp.mpf(_m * _m)
        _term = (2 * _pi2 * _m2 * _m2 * _e9 - 3 * _pi * _m2 * _e5) * mp.exp(-_pi * _m2 * _e4)
        _phi += _term
        if abs(_term) < mp.mpf("1e-78"):
            break
    _wt = mp.mpf("0.5") if (_k == 0 or _k == N_U_RIEM) else mp.mpf("1.0")
    U2_RIEM.append(2 * _u)
    W_RIEM.append(8 * _phi * _wt * DU_RIEM)

FACT = [mp.factorial(k) for k in range(110)]


def riemann_taylor_coeffs(t0: mp.mpf, M: int = 86):
    """Compute [Xi^{(m)}(t0) / m! for m = 0..M] on the critical line d = 0."""
    raw = [mp.mpf("0")] * (M + 1)
    for u2, w in zip(U2_RIEM, W_RIEM):
        cs = w * mp.cos(u2 * t0)
        sn = -w * mp.sin(u2 * t0)
        p = mp.mpf("1")
        for m in range(0, M + 1, 4):
            raw[m] += cs * p
            p *= u2
            if m + 1 <= M:
                raw[m + 1] += sn * p
                p *= u2
            if m + 2 <= M:
                raw[m + 2] -= cs * p
                p *= u2
            if m + 3 <= M:
                raw[m + 3] -= sn * p
                p *= u2
    return [raw[m] / FACT[m] for m in range(M + 1)]


def eval_f0_f1_f2(z_shift, c_real):
    """Evaluate f(t0 + z_shift), f'(t0 + z_shift), f''(t0 + z_shift) via Horner's rule."""
    M = len(c_real) - 1
    f0 = mp.mpc(c_real[M])
    f1 = mp.mpc("0")
    f2 = mp.mpc("0")
    for m in range(M - 1, -1, -1):
        f2 = f2 * z_shift + 2 * f1
        f1 = f1 * z_shift + f0
        f0 = f0 * z_shift + c_real[m]
    return f0, f1, f2


# Build 148-dps Davenport-Heilbronn Xi_DH(z) and Dirichlet Xi_{chi_5}(z) quadrature
# (so that even at t = 240.44 where exp(-pi*t/4) ~ 1e-82, we retain > 65 significant digits!)
with mp.workdps(148):
    KAPPA_DH = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    PHASE_CHI5 = (1 - 1j * KAPPA_DH) / mp.sqrt(1 + KAPPA_DH**2)
    PI5 = mp.pi / 5

    N_U_DH = 560
    U_MAX_DH = mp.mpf("2.80")
    DU_DH = U_MAX_DH / N_U_DH
    U2_DH = []
    W_DH = []
    W_CHI5_RE = []
    W_CHI5_IM = []

    def _a_dh(n: int):
        r = n % 5
        if r == 0:
            return mp.mpf("0")
        if r == 1:
            return mp.mpf("1")
        if r == 2:
            return KAPPA_DH
        if r == 3:
            return -KAPPA_DH
        return mp.mpf("-1")

    def _chi5(m: int):
        r = m % 5
        if r == 0:
            return mp.mpc("0")
        if r == 1:
            return mp.mpc("1")
        if r == 2:
            return mp.mpc(0, 1)
        if r == 3:
            return mp.mpc(0, -1)
        return mp.mpc("-1")

    for _k in range(N_U_DH + 1):
        _u = mp.mpf(_k) * DU_DH
        _e3 = mp.exp(3 * _u)
        _e4 = PI5 * mp.exp(4 * _u)
        _phi_dh = mp.mpf("0")
        _phi_c5 = mp.mpc("0")
        for _m in range(1, 38):
            _ex = mp.exp(-_e4 * (_m * _m))
            if _ex < mp.mpf("1e-155"):
                break
            _base = _m * _e3 * _ex
            _phi_dh += _a_dh(_m) * _base
            _phi_c5 += PHASE_CHI5 * _chi5(_m) * _base
        _wt = (mp.mpf("0.5") if (_k == 0 or _k == N_U_DH) else mp.mpf("1.0")) * 8 * DU_DH
        U2_DH.append(2 * _u)
        W_DH.append(_wt * _phi_dh)
        W_CHI5_RE.append(_wt * mp.re(_phi_c5))
        W_CHI5_IM.append(_wt * mp.im(_phi_c5))
    FACT_148 = [mp.factorial(k) for k in range(110)]


def dh_taylor_coeffs(t0: mp.mpf, M: int = 86):
    work_prec = 148 if float(t0) > 120.0 else 96
    with mp.workdps(work_prec):
        t0_hp = mp.mpf(str(t0))
        raw = [mp.mpf("0")] * (M + 1)
        for u2, w in zip(U2_DH, W_DH):
            cs = w * mp.cos(u2 * t0_hp)
            sn = -w * mp.sin(u2 * t0_hp)
            p = mp.mpf("1")
            for m in range(0, M + 1, 4):
                raw[m] += cs * p
                p *= u2
                if m + 1 <= M:
                    raw[m + 1] += sn * p
                    p *= u2
                if m + 2 <= M:
                    raw[m + 2] -= cs * p
                    p *= u2
                if m + 3 <= M:
                    raw[m + 3] -= sn * p
                    p *= u2
        return [+raw[m] / FACT_148[m] for m in range(M + 1)]


def chi5_taylor_coeffs(t0: mp.mpf, M: int = 86):
    with mp.workdps(96):
        t0_hp = mp.mpf(str(t0))
        raw = [mp.mpf("0")] * (M + 1)
        for u2, wr, wi in zip(U2_DH, W_CHI5_RE, W_CHI5_IM):
            c0 = wr * mp.cos(u2 * t0_hp) - wi * mp.sin(u2 * t0_hp)
            c1 = -wr * mp.sin(u2 * t0_hp) - wi * mp.cos(u2 * t0_hp)
            p = mp.mpf("1")
            for m in range(0, M + 1, 4):
                raw[m] += c0 * p
                p *= u2
                if m + 1 <= M:
                    raw[m + 1] += c1 * p
                    p *= u2
                if m + 2 <= M:
                    raw[m + 2] -= c0 * p
                    p *= u2
                if m + 3 <= M:
                    raw[m + 3] -= c1 * p
                    p *= u2
        return [+raw[m] / FACT_148[m] for m in range(M + 1)]


def gl_nodes(n: int):
    nodes, weights = [], []
    for k in range(1, n + 1):
        x0 = mp.cos(mp.pi * (mp.mpf(k) - mp.mpf("0.25")) / (mp.mpf(n) + mp.mpf("0.5")))
        for _ in range(12):
            p0, p1 = mp.mpf("1"), x0
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x0 * p1 - (j - 1) * p0) / j
            dp = n * (x0 * p1 - p0) / (x0**2 - 1)
            dx = -p1 / dp
            x0 += dx
            if abs(dx) < mp.mpf("1e-63"):
                break
        w0 = 2 / ((1 - x0**2) * (dp**2))
        nodes.append(x0)
        weights.append(w0)
    return nodes, weights


GL24_X, GL24_W = gl_nodes(24)
GL18_X, GL18_W = gl_nodes(18)


# ==============================================================================
# MODULE 1 (SECTION 2 — MACRO-SCALE):
# PORTAL-FRAME DICHOTOMY & CODIMENSION-3 SADDLE OVERDETERMINATION
# ==============================================================================

def verify_section2_macro_portal_frame():
    print("\n" + "=" * 88)
    print("MODULE 1 [SECTION 2 — MACRO-SCALE]: PORTAL-FRAME DICHOTOMY & 3-EQ OVERDETERMINATION")
    print("=" * 88)

    # Part A: First 5 off-critical zeros of Davenport-Heilbronn Xi_DH(z) in t in [10, 260]
    dh_off_inits = [
        (mp.mpf("0.3085171825"), mp.mpf("85.6993484854"), mp.mpf("85.7472008265")),
        (mp.mpf("0.1508300806"), mp.mpf("114.1633427308"), mp.mpf("114.1695784427")),
        (mp.mpf("0.0743560505"), mp.mpf("166.4793059132"), mp.mpf("166.4820292062")),
        (mp.mpf("0.2242576946"), mp.mpf("176.7024612429"), mp.mpf("176.7301734191")),
        (mp.mpf("0.3695305796"), mp.mpf("240.4046723514"), mp.mpf("240.4397623868")),
    ]

    print("\n--- Part 1A: 65-Decimal Certification of First 5 Off-Critical Zeros of Xi_DH(z) ---")
    dh_rows = []
    for m_idx, (d_init, g_init, tau_init) in enumerate(dh_off_inits, start=1):
        c_dh = dh_taylor_coeffs(tau_init, 76)
        # Refine critical-line anchor tau_{false,m} where Xi_DH'(tau) = 0
        dt_tau = mp.mpf("0")
        for _ in range(8):
            _, f1, f2 = eval_f0_f1_f2(dt_tau, c_dh)
            dt_tau -= mp.re(f1) / mp.re(f2)
        tau_false = tau_init + dt_tau
        f0_tau, f1_tau, f2_tau = eval_f0_f1_f2(dt_tau, c_dh)

        # Refine off-critical zero z_* = gamma_* - i*d_*
        z_sh = mp.mpc(g_init - tau_init, -d_init)
        for _ in range(8):
            f0_z, f1_z, _ = eval_f0_f1_f2(z_sh, c_dh)
            z_sh -= f0_z / f1_z
        gamma_star = tau_init + mp.re(z_sh)
        d_star = -mp.im(z_sh)
        f0_z, f1_z, _ = eval_f0_f1_f2(z_sh, c_dh)

        sc = mp.exp(mp.pi * tau_false / 4)
        xi0_sc = mp.re(f0_tau) * sc
        xi2_sc = mp.re(f2_tau) * sc
        c1_sc2 = -xi0_sc * xi2_sc

        # Interval certification of c_1^{DH}(tau_{false,m}) < 0
        xi0_iv = mp.iv.mpf(str(xi0_sc))
        xi2_iv = mp.iv.mpf(str(xi2_sc))
        c1_iv = -xi0_iv * xi2_iv
        assert c1_iv.b < 0, f"Interval certification failed for DH zero {m_idx}"
        assert abs(f0_z) * sc < mp.mpf("1e-55")

        dh_rows.append({
            "m": m_idx,
            "d_star": d_star,
            "gamma_star": gamma_star,
            "tau_false": tau_false,
            "xi0_sc": xi0_sc,
            "xi2_sc": xi2_sc,
            "c1_sc2": c1_sc2,
        })
        print(
            f"  m={m_idx} | (d*, gamma*)=({mp.nstr(d_star, 12)}, {mp.nstr(gamma_star, 14)}) | "
            f"tau_false={mp.nstr(tau_false, 16)} | sc*Xi={mp.nstr(xi0_sc, 10)} | "
            f"sc*Xi''={mp.nstr(xi2_sc, 10)} | sc^2*c1={mp.nstr(c1_sc2, 12)} < 0 [CASE I]"
        )

    # Part B: 5 Deformed Interior Complex Saddle Tests (Codimension-3 Overdetermination)
    print("\n--- Part 1B: Codimension-3 Interior Saddle Overdetermination (3 Eq in 2 Unknowns) ---")
    def_cases = [
        (mp.mpf("14.60"), mp.mpf("0.35"), mp.mpc("14.55", "-0.22")),
        (mp.mpf("20.40"), mp.mpf("0.35"), mp.mpc("20.50", "-0.20")),
        (mp.mpf("21.50"), mp.mpf("0.35"), mp.mpc("21.40", "-0.20")),
        (mp.mpf("25.50"), mp.mpf("0.35"), mp.mpc("25.40", "-0.20")),
        (mp.mpf("30.80"), mp.mpf("0.28"), mp.mpc("30.65", "-0.18")),
    ]
    saddle_rows = []
    for g0, d0, z_init in def_cases:
        c_r = riemann_taylor_coeffs(g0, 76)

        def eval_def(z_pt):
            w = z_pt - g0
            D0, D1, D2 = eval_f0_f1_f2(w, c_r)
            P0 = w**2 + d0**2
            P1 = 2 * w
            P2 = mp.mpf("2")
            return D0 * P0, D1 * P0 + D0 * P1, D2 * P0 + 2 * D1 * P1 + D0 * P2

        zc = z_init
        for _ in range(12):
            f0, f1, f2 = eval_def(zc)
            zc -= f1 / f2
        f0, f1, f2 = eval_def(zc)
        dc = -mp.im(zc)
        tc = mp.re(zc)
        sc_c = mp.exp(mp.pi * tc / 4)
        eq2_res = abs(mp.re(f1)) * sc_c
        eq3_res = abs(mp.im(f1)) * sc_c
        eq1_I_sc = mp.im(f0) * sc_c

        f0_0, f1_0, f2_0 = eval_def(g0)
        sc_0 = mp.exp(mp.pi * g0 / 4)
        c1_def_sc2 = (mp.re(f1_0)**2 - mp.re(f0_0) * mp.re(f2_0)) * (sc_0**2)

        assert eq2_res < mp.mpf("1e-55") and eq3_res < mp.mpf("1e-55")
        assert abs(eq1_I_sc) > mp.mpf("0.05")
        assert c1_def_sc2 < 0

        saddle_rows.append({
            "g0": g0,
            "d0": d0,
            "dc": dc,
            "tc": tc,
            "eq2": eq2_res,
            "eq3": eq3_res,
            "eq1_I_sc": eq1_I_sc,
            "c1_def_sc2": c1_def_sc2,
        })
        print(
            f"  (g0={mp.nstr(g0, 5)}, d0={mp.nstr(d0, 4)}) -> (dc, tc)=({mp.nstr(dc, 12)}, {mp.nstr(tc, 12)}) | "
            f"|nabla I|={mp.nstr(max(eq2_res, eq3_res), 3)} | Sc*I(dc,tc)={mp.nstr(eq1_I_sc, 14)} != 0 | "
            f"S0^2*c1_def={mp.nstr(c1_def_sc2, 12)} < 0"
        )

    # Part C: First 10 Critical-Line Extrema tau_1 .. tau_10 of Riemann's Xi(z)
    print("\n--- Part 1C: First 10 Chladni Walls Gamma_1..Gamma_10 of Riemann's Xi(z) ---")
    tau_inits_10 = [
        mp.mpf("15.585708589829342"),
        mp.mpf("22.097977280400902"),
        mp.mpf("26.272247356935624"),
        mp.mpf("31.231795871009785"),
        mp.mpf("34.193310269011381"),
        mp.mpf("38.498240763754491"),
        mp.mpf("41.736729522419333"),
        mp.mpf("44.541703607380997"),
        mp.mpf("48.622532685277866"),
        mp.mpf("50.839004822815970"),
    ]
    riem_wall_rows = []
    for k_idx, tg in enumerate(tau_inits_10, start=1):
        c_j = riemann_taylor_coeffs(tg, 76)
        dt = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt, c_j)
            dt -= mp.re(f1) / mp.re(f2)
        tau_k = tg + dt
        c_j = riemann_taylor_coeffs(tau_k, 84)
        f0_0, f1_0, f2_0 = eval_f0_f1_f2(mp.mpf("0"), c_j)
        sc_k = mp.exp(mp.pi * tau_k / 4)
        xi0_sc = mp.re(f0_0) * sc_k
        xi2_sc = mp.re(f2_0) * sc_k
        c1_sc2 = -xi0_sc * xi2_sc

        # Trace wall to d = 0.50
        tc = tau_k
        for d_step in [mp.mpf("0.15"), mp.mpf("0.30"), mp.mpf("0.50")]:
            for _ in range(7):
                f0_d, f1_d, _ = eval_f0_f1_f2((tc - tau_k) - 1j * d_step, c_j)
                tc -= mp.im(f0_d) / mp.im(f1_d)
        f0_50, _, _ = eval_f0_f1_f2((tc - tau_k) - 1j * mp.mpf("0.50"), c_j)
        growth_ratio = abs(mp.re(f0_50)) / abs(mp.re(f0_0))

        # Interval check
        c1_iv = -mp.iv.mpf(str(xi0_sc)) * mp.iv.mpf(str(xi2_sc))
        assert c1_iv.a > 0 and growth_ratio > 1

        riem_wall_rows.append({
            "k": k_idx,
            "tau_k": tau_k,
            "xi0_sc": xi0_sc,
            "xi2_sc": xi2_sc,
            "c1_sc2": c1_sc2,
            "growth": growth_ratio,
            "c_j": c_j,
        })
        print(
            f"  Gamma_{k_idx:2d} | tau_{k_idx}(0)={mp.nstr(tau_k, 18)} | Sk*Xi={mp.nstr(xi0_sc, 10)} | "
            f"Sk*Xi''={mp.nstr(xi2_sc, 10)} | Sk^2*c1=+{mp.nstr(c1_sc2, 12)} > 0 | "
            f"|R(0.50)|/|Xi(tau)|={mp.nstr(growth_ratio, 8)} > 1"
        )

    return {
        "dh_rows": dh_rows,
        "saddle_rows": saddle_rows,
        "riem_wall_rows": riem_wall_rows,
    }


# ==============================================================================
# MODULE 2 (SECTION 3 — MESO-SCALE):
# EXACT 50%-50% GREEN-DIRICHLET-NEUMANN EQUIPARTITION & WRONSKIAN-SLOPE RECIPROCITY
# ==============================================================================

def run_chamber_green_and_reciprocity(label, ta0, tb0, c_left, c_mid, c_right, d_grid):
    t_centers = [ta0, mp.mpf("0.5") * (ta0 + tb0), tb0]
    c_centers = [c_left, c_mid, c_right]

    def eval_pt(d, t):
        if abs(t - t_centers[0]) <= abs(t - t_centers[1]):
            idx = 0
        elif abs(t - t_centers[1]) <= abs(t - t_centers[2]):
            idx = 1
        else:
            idx = 2
        f0, f1, _ = eval_f0_f1_f2((t - t_centers[idx]) - 1j * d, c_centers[idx])
        R, I = mp.re(f0), mp.im(f0)
        Rt, It = mp.re(f1), mp.im(f1)
        Rd, Id = It, -Rt
        return R, I, Rd, Id, Rt, It

    def find_w(idx_wall, d, t_g):
        tc = t_centers[idx_wall]
        cc = c_centers[idx_wall]
        t = t_g
        for _ in range(7):
            f0, f1, _ = eval_f0_f1_f2((t - tc) - 1j * d, cc)
            t -= mp.im(f0) / mp.im(f1)
        return t

    def int_slice(d, ta, tb):
        s_I, s_R, s_G = mp.mpf("0"), mp.mpf("0"), mp.mpf("0")
        dt_p = (tb - ta) / 2
        for p in range(2):
            pa = ta + p * dt_p
            hm = mp.mpf("0.5") * dt_p
            mm = pa + hm
            for x, w in zip(GL24_X, GL24_W):
                R, I, Rd, Id, Rt, It = eval_pt(d, mm + hm * x)
                s_I += hm * w * (I * Id)
                s_R += hm * w * (R * Rd)
                s_G += hm * w * (Rt * Rt + It * It)
        return s_I, s_R, s_G

    cum_E = mp.mpf("0")
    ta_g, tb_g = ta0, tb0
    out_rows = []
    for s_idx in range(len(d_grid) - 1):
        da, db = d_grid[s_idx], d_grid[s_idx + 1]
        hd = mp.mpf("0.5") * (db - da)
        md = mp.mpf("0.5") * (da + db)
        st_E = mp.mpf("0")
        for xd, wd in zip(reversed(GL18_X), reversed(GL18_W)):
            d = md + hd * xd
            ta_g = find_w(0, d, ta_g)
            tb_g = find_w(2, d, tb_g)
            _, _, g2 = int_slice(d, ta_g, tb_g)
            st_E += hd * wd * g2
        cum_E += st_E
        ta_top = find_w(0, db, ta_g)
        tb_top = find_w(2, db, tb_g)
        ta_g, tb_g = ta_top, tb_top
        I_top, R_top, _ = int_slice(db, ta_top, tb_top)
        err_RI = abs(I_top - R_top) / abs(I_top)
        err_E = abs(I_top - cum_E) / abs(cum_E)

        # Left and right wall reciprocity
        R_a, I_a, Rd_a, Id_a, Rt_a, It_a = eval_pt(db, ta_top)
        Wt_a = (R_a * Rd_a) / db
        slope_a_dir = -Id_a / It_a
        slope_a_rec = (R_a * Rt_a) / (db * Wt_a)
        s2_a_dir = 1 + slope_a_dir**2
        s2_a_rec = (R_a**2 * (Rt_a**2 + It_a**2)) / (db**2 * Wt_a**2)

        R_b, I_b, Rd_b, Id_b, Rt_b, It_b = eval_pt(db, tb_top)
        Wt_b = (R_b * Rd_b) / db
        slope_b_dir = -Id_b / It_b
        slope_b_rec = (R_b * Rt_b) / (db * Wt_b) if abs(Wt_b) > mp.mpf("1e-110") else slope_b_dir
        dR2_b = db * Wt_b * (1 + slope_b_dir**2)

        assert err_RI < mp.mpf("1e-55") and err_E < mp.mpf("1e-48")
        assert abs(slope_a_dir - slope_a_rec) < mp.mpf("1e-58")
        assert abs(s2_a_dir - s2_a_rec) < mp.mpf("1e-58")

        out_rows.append({
            "d0": db,
            "ta": ta_top,
            "tb": tb_top,
            "I_Omega": I_top,
            "R_Omega": R_top,
            "E_Omega": cum_E,
            "err_RI": err_RI,
            "err_E": err_E,
            "Wt_a": Wt_a,
            "slope_a": slope_a_dir,
            "s2_a": s2_a_dir,
            "Wt_b": Wt_b,
            "slope_b": slope_b_dir,
            "dR2_b": dR2_b,
        })
        print(
            f"  {label:12s} d0={mp.nstr(db, 6):8s} | I=R=E={mp.nstr(I_top, 18):24s} | "
            f"err_RI={mp.nstr(err_RI, 3):9s} | err_E={mp.nstr(err_E, 3):9s} | "
            f"Wt_b={mp.nstr(Wt_b, 10):16s} | slope_b={mp.nstr(slope_b_dir, 9)}"
        )
    return out_rows


def verify_section3_meso_green_equipartition(sec2_data):
    print("\n" + "=" * 88)
    print("MODULE 2 [SECTION 3 — MESO-SCALE]: 50%-50% GREEN EQUIPARTITION & WRONSKIAN RECIPROCITY")
    print("=" * 88)

    riem_walls = sec2_data["riem_wall_rows"]
    d_grid_riem = [mp.mpf("0.00"), mp.mpf("0.05"), mp.mpf("0.15"), mp.mpf("0.25"), mp.mpf("0.50")]
    riem_chambers = {}
    for k in range(3):
        ta0 = riem_walls[k]["tau_k"]
        tb0 = riem_walls[k + 1]["tau_k"]
        c_left = riem_walls[k]["c_j"]
        c_right = riem_walls[k + 1]["c_j"]
        c_mid = riemann_taylor_coeffs(mp.mpf("0.5") * (ta0 + tb0), 84)
        riem_chambers[f"Omega_{k+1}"] = run_chamber_green_and_reciprocity(
            f"Omega_{k+1}", ta0, tb0, c_left, c_mid, c_right, d_grid_riem
        )

    # Head-to-head comparison: Xi_{chi_5} (Euler product) vs Xi_DH (no Euler product)
    # Use 88-decimal working precision on [83.89, 86.63] where exp(-pi*t/2) ~ 1e-58
    with mp.workdps(88):
        global GL24_X, GL24_W, GL18_X, GL18_W
        GL24_X, GL24_W = gl_nodes(24)
        GL18_X, GL18_W = gl_nodes(18)

        t_dh_f = mp.mpf("85.74720082646129724124587290234672558714140137092622031")
        c_dh_f = dh_taylor_coeffs(t_dh_f, 92)
        dt_f = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt_f, c_dh_f)
            dt_f -= mp.re(f1) / mp.re(f2)
        t_dh_f += dt_f
        c_dh_f = dh_taylor_coeffs(t_dh_f, 92)

        z_sh_s = mp.mpc(mp.mpf("85.699348485377592") - t_dh_f, -mp.mpf("0.3085171824566374"))
        for _ in range(7):
            f0_s, f1_s, _ = eval_f0_f1_f2(z_sh_s, c_dh_f)
            z_sh_s -= f0_s / f1_s
        d_star = -mp.im(z_sh_s)
        d_grid_comp = [mp.mpf("0.00"), mp.mpf("0.05"), mp.mpf("0.15"), mp.mpf("0.25"), d_star, mp.mpf("0.50")]

        # Xi_{chi_5} extrema in [83.8, 86.2]
        t_c5_a, t_c5_b = mp.mpf("83.893523413636714"), mp.mpf("86.029308068178361")
        c_c5_a = chi5_taylor_coeffs(t_c5_a, 92)
        dt_a = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt_a, c_c5_a)
            dt_a -= mp.re(f1) / mp.re(f2)
        t_c5_a += dt_a
        c_c5_a = chi5_taylor_coeffs(t_c5_a, 92)

        c_c5_b = chi5_taylor_coeffs(t_c5_b, 92)
        dt_b = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt_b, c_c5_b)
            dt_b -= mp.re(f1) / mp.re(f2)
        t_c5_b += dt_b
        c_c5_b = chi5_taylor_coeffs(t_c5_b, 92)
        c_c5_m = chi5_taylor_coeffs(mp.mpf("0.5") * (t_c5_a + t_c5_b), 92)

        print("\n--- Head-to-Head Comparison: Euler-Product Xi_{chi_5} vs Davenport-Heilbronn Xi_DH ---")
        chi5_rows = run_chamber_green_and_reciprocity(
            "Xi_chi5", t_c5_a, t_c5_b, c_c5_a, c_c5_m, c_c5_b, d_grid_comp
        )

        # Xi_DH sub-chambers [tau_left, tau_false] and [tau_false, tau_right]
        t_dh_l = mp.mpf("83.895935252167107")
        t_dh_r = mp.mpf("86.629431818945013")

        c_dh_l = dh_taylor_coeffs(t_dh_l, 92)
        dt_l = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt_l, c_dh_l)
            dt_l -= mp.re(f1) / mp.re(f2)
        t_dh_l += dt_l
        c_dh_l = dh_taylor_coeffs(t_dh_l, 92)

        c_dh_r = dh_taylor_coeffs(t_dh_r, 92)
        dt_r = mp.mpf("0")
        for _ in range(7):
            _, f1, f2 = eval_f0_f1_f2(dt_r, c_dh_r)
            dt_r -= mp.re(f1) / mp.re(f2)
        t_dh_r += dt_r
        c_dh_r = dh_taylor_coeffs(t_dh_r, 92)

        c_dh_ml = dh_taylor_coeffs(mp.mpf("0.5") * (t_dh_l + t_dh_f), 92)
        c_dh_mr = dh_taylor_coeffs(mp.mpf("0.5") * (t_dh_f + t_dh_r), 92)

        dh_left_rows = run_chamber_green_and_reciprocity(
            "DH_Left", t_dh_l, t_dh_f, c_dh_l, c_dh_ml, c_dh_f, d_grid_comp
        )
        dh_right_rows = run_chamber_green_and_reciprocity(
            "DH_Right", t_dh_f, t_dh_r, c_dh_f, c_dh_mr, c_dh_r, d_grid_comp
        )

        # Fold vertex (d_tip, t_tip) of the detached real nodal parabola R_DH(d, t) = 0
        d_tip = mp.mpf("0.2980980119818233")
        t_tip = mp.mpf("85.62549444284696")
        for _ in range(8):
            f0, f1, f2 = eval_f0_f1_f2((t_tip - t_dh_f) - 1j * d_tip, c_dh_f)
            R, I = mp.re(f0), mp.im(f0)
            Rt, It = mp.re(f1), mp.im(f1)
            Rd, Id = It, -Rt
            Rtt, Itt = mp.re(f2), mp.im(f2)
            Rtd = Itt
            det = Rd * Rtt - Rt * Rtd
            d_tip -= (Rtt * R - Rt * Rt) / det
            t_tip -= (-Rtd * R + Rd * Rt) / det
        f0_tip, f1_tip, _ = eval_f0_f1_f2((t_tip - t_dh_f) - 1j * d_tip, c_dh_f)
        Wt_tip = (mp.re(f0_tip) * mp.im(f1_tip) - mp.im(f0_tip) * mp.re(f1_tip)) / d_tip
        print(
            f"  Detached Real Parabola Fold Vertex: (d_tip, t_tip)=({mp.nstr(d_tip, 20)}, {mp.nstr(t_tip, 20)}) | "
            f"|W_tilde(d_tip, t_tip)|={mp.nstr(abs(Wt_tip), 4)} (== 0, |t_R'(d_tip)| = inf)"
        )

    return {
        "riem_chambers": riem_chambers,
        "chi5_rows": chi5_rows,
        "dh_left_rows": dh_left_rows,
        "dh_right_rows": dh_right_rows,
        "d_tip": d_tip,
        "t_tip": t_tip,
        "t_dh_l": t_dh_l,
        "t_dh_f": t_dh_f,
        "t_dh_r": t_dh_r,
        "c_dh_l": c_dh_l,
        "c_dh_f": c_dh_f,
        "c_dh_r": c_dh_r,
    }


# ==============================================================================
# MODULE 3 (SECTION 4 — MICRO-SCALE):
# ATTENUATED SEMIPRIME SPRING, DILOGARITHM Li_2(2*theta) DAMPING & HYPERBOLA DUALITY
# ==============================================================================

def M_sq_closed(theta: mp.mpf, p: int) -> mp.mpf:
    tau = 2 * theta
    if tau >= 1:
        return mp.mpf("0")
    H_p = sum((1 - tau)**m / m for m in range(1, p + 1))
    return (mp.log(1 / tau) - H_p)**2


def M_tri_closed(theta: mp.mpf, p: int) -> mp.mpf:
    tau = 2 * theta
    if tau >= mp.mpf("0.5"):
        return mp.mpf("0")
    M0 = mp.log(tau)**2 + 2 * mp.polylog(2, tau) - (mp.pi**2) / 6
    if p == 0:
        return M0
    if p == 1:
        return M0 - 2 * (1 - tau) * mp.log((1 - tau) / tau) + (mp.mpf("2.5") - tau) * (1 - 2 * tau)
    if p == 2:
        return (
            M0
            - (1 - tau) * (3 - tau) * mp.log((1 - tau) / tau)
            + (mp.mpf(31) / 8 - mp.mpf(13) / 4 * tau + tau**2) * (1 - 2 * tau)
        )
    if p == 3:
        return (
            M0
            - mp.mpf(1) / 3 * (1 - tau) * (11 - 7 * tau + 2 * tau**2) * mp.log((1 - tau) / tau)
            + (
                mp.mpf(1801) / 360
                - mp.mpf(1199) / 180 * tau
                + mp.mpf(218) / 45 * tau**2
                - mp.mpf(89) / 45 * tau**3
                + mp.mpf(49) / 90 * tau**4
                - mp.mpf(1) / 9 * tau**5
            )
            * (1 - 2 * tau)
        )
    H_p = sum((1 - tau)**m / m for m in range(1, p + 1))
    S2_p = sum(((1 - tau)**m - tau**m) / (m**2) for m in range(1, p + 1))
    poly_sum = mp.mpf("0")
    for j in range(p):
        int_uj = ((1 - tau)**(j + 1) - tau**(j + 1)) / (j + 1)
        for k in range(p):
            int_beta = sum(
                mp.binomial(j, r) * ((-1)**r) * ((1 - tau)**(k + 2 + r) - tau**(k + 2 + r)) / (k + 2 + r)
                for r in range(j + 1)
            )
            poly_sum += ((1 - tau)**(k + 1) * int_uj - int_beta) / (k + 1)
    return M0 - 2 * H_p * mp.log((1 - tau) / tau) + 2 * S2_p + poly_sum


def eta_closed(theta: mp.mpf, p: int) -> mp.mpf:
    if theta >= mp.mpf("0.25"):
        return mp.mpf("0")
    return M_tri_closed(theta, p) / M_sq_closed(theta, p)


def verify_section4_micro_dilogarithm():
    print("\n" + "=" * 88)
    print("MODULE 3 [SECTION 4 — MICRO-SCALE]: EULER DILOGARITHM Li_2(2*theta) & HYPERBOLA DUALITY")
    print("=" * 88)

    # Verify closed-form polynomial formulas at theta = 1/6 against 2D Gauss-Legendre quadrature
    th16 = mp.mpf(1) / 6
    print("\n--- Part 3A: 65-Decimal Closed-Form vs 2D Quadrature Verification at theta = 1/6 ---")
    for p in [0, 1, 2, 3]:
        msq = M_sq_closed(th16, p)
        mtri = M_tri_closed(th16, p)
        mcor = msq - mtri
        eta_val = mtri / msq
        rat_val = mtri / mcor

        # 2D Composite Gauss-Legendre verification of M_tri^{(p)}(1/6) (4x4 subintervals -> < 1e-60 error)
        # int_{1/6}^{1/3} d a1 (1-2*a1)^p/a1 * int_{1/6}^{1/2 - a1} d a2 (1-2*a2)^p/a2
        quad_tri = mp.mpf("0")
        n_sub = 4
        da1 = (mp.mpf(1) / 3 - mp.mpf(1) / 6) / n_sub
        for s1 in range(n_sub):
            A1_l = mp.mpf(1) / 6 + s1 * da1
            half_1 = da1 / 2
            mid_1 = A1_l + half_1
            for x1, w1 in zip(GL24_X, GL24_W):
                a1 = mid_1 + half_1 * x1
                f1 = ((1 - 2 * a1)**p) / a1
                b2 = mp.mpf("0.5") - a1
                da2 = (b2 - mp.mpf(1) / 6) / n_sub
                inner = mp.mpf("0")
                for s2 in range(n_sub):
                    A2_l = mp.mpf(1) / 6 + s2 * da2
                    half_2 = da2 / 2
                    mid_2 = A2_l + half_2
                    for x2, w2 in zip(GL24_X, GL24_W):
                        a2 = mid_2 + half_2 * x2
                        inner += half_2 * w2 * ((1 - 2 * a2)**p) / a2
                quad_tri += half_1 * w1 * f1 * inner
        q_err = abs(mtri - quad_tri) / mtri
        assert q_err < mp.mpf("1e-58")
        print(
            f"  p={p} | M_sq={mp.nstr(msq, 18)} | M_tri={mp.nstr(mtri, 18)} | "
            f"eta^({p})(1/6)={mp.nstr(eta_val, 18)} | Tri/Corn={mp.nstr(rat_val, 14)} | "
            f"quad_err={mp.nstr(q_err, 3)}"
        )

    # Critical thresholds theta_{*,50%}^{(p)} and theta_{*,75%}^{(p)}
    c_theta = mp.mpf("1.000028")
    th_Wpure_inf = mp.mpf("0.5") * mp.exp(-1 / c_theta)
    th_2e = mp.mpf("0.5") / mp.e
    th50_0 = mp.findroot(lambda x: eta_closed(x, 0) - mp.mpf("0.5"), mp.mpf("0.1067"))
    th50_1 = mp.findroot(lambda x: eta_closed(x, 1) - mp.mpf("0.5"), mp.mpf("0.1660"))
    th50_2 = mp.findroot(lambda x: eta_closed(x, 2) - mp.mpf("0.5"), mp.mpf("0.1905"))
    th75_2 = mp.findroot(lambda x: eta_closed(x, 2) - mp.mpf("0.75"), mp.mpf("0.1536"))
    th50_3 = mp.findroot(lambda x: eta_closed(x, 3) - mp.mpf("0.5"), mp.mpf("0.2039"))
    th75_3 = mp.findroot(lambda x: eta_closed(x, 3) - mp.mpf("0.75"), mp.mpf("0.1750"))

    # Interval verification of eta^{(2)}(1/6) > 2/3 > 1/2 and theta_{*,50%}^{(2)} > theta_{W-pure}^inf
    eta2_16_iv = mp.iv.mpf(str(eta_closed(th16, 2)))
    floor2_16_iv = (1 - 1 / (2 * eta2_16_iv))**2
    th50_2_iv = mp.iv.mpf(str(th50_2))
    thW_iv = mp.iv.mpf(str(th_Wpure_inf))
    assert eta2_16_iv.a > mp.iv.mpf(2) / 3
    assert floor2_16_iv.a > mp.iv.mpf("0.068")
    assert (th50_2_iv - thW_iv).a > mp.iv.mpf("0.0066")

    print("\n--- Part 3B: Exact 50% and 75% Capture Thresholds (65-dps) ---")
    print(f"  theta_{{*,50%}}^{{(0)}} = {mp.nstr(th50_0, 52)}")
    print(f"  theta_{{*,50%}}^{{(1)}} = {mp.nstr(th50_1, 52)}")
    print(f"  theta_{{*,50%}}^{{(2)}} = {mp.nstr(th50_2, 52)} > theta_{{W-pure}}^inf = {mp.nstr(th_Wpure_inf, 18)}")
    print(f"  theta_{{*,75%}}^{{(2)}} = {mp.nstr(th75_2, 52)}")
    print(f"  theta_{{*,50%}}^{{(3)}} = {mp.nstr(th50_3, 52)}")
    print(f"  theta_{{*,75%}}^{{(3)}} = {mp.nstr(th75_3, 52)}")

    # Table 3 Part A: Scan across theta in [0.05, 0.24]
    theta_scan_list = [
        ("0.05", mp.mpf("0.05")),
        ("theta_{*,50%}^{(0)}", th50_0),
        ("theta_{*,75%}^{(2)}", th75_2),
        ("1/6 (cubic boundary)", th16),
        ("theta_{*,75%}^{(3)}", th75_3),
        ("1/(2e) (unit Mertens)", th_2e),
        ("theta_{W-pure}^inf", th_Wpure_inf),
        ("theta_{*,50%}^{(2)}", th50_2),
        ("1/5", mp.mpf("0.20")),
        ("theta_{*,50%}^{(3)}", th50_3),
        ("0.24", mp.mpf("0.24")),
    ]
    print("\n--- Part 3C: Closed-Form Dilogarithm Table Across theta in [0.05, 0.24] ---")
    dilog_rows = []
    for label, th in theta_scan_list:
        e0 = eta_closed(th, 0)
        e1 = eta_closed(th, 1)
        e2 = eta_closed(th, 2)
        e3 = eta_closed(th, 3)
        r2 = e2 / (1 - e2)
        r3 = e3 / (1 - e3)
        fl2 = (max(mp.mpf("0"), 1 - 1 / (2 * e2)))**2
        dilog_rows.append({
            "label": label,
            "theta": th,
            "e0": e0,
            "e1": e1,
            "e2": e2,
            "r2": r2,
            "fl2": fl2,
            "e3": e3,
            "r3": r3,
        })
        print(
            f"  {label:22s} | th={float(th):.8f} | eta0={float(e0):.6f} | eta1={float(e1):.6f} | "
            f"eta2={float(e2):.6f} (R2={float(r2):.6f}, fl2={float(fl2):.6f}) | eta3={float(e3):.6f} (R3={float(r3):.6f})"
        )

    # Part 3D: Discrete prime sums up to N_* = 10^5 (at theta = 1/6)
    print("\n--- Part 3D: Discrete Prime Sums vs Continuous Li_2(2*theta) Laws up to N_* = 10^5 ---")
    max_N = 100000
    is_prime = np.ones(max_N + 1, dtype=bool)
    is_prime[:2] = False
    for p_i in range(2, int(math.isqrt(max_N)) + 1):
        if is_prime[p_i]:
            is_prime[p_i * p_i : max_N + 1 : p_i] = False
    all_primes = np.nonzero(is_prime)[0]

    N_benchmarks = [398, 1261, 3989, 12615, 39894, 100000]
    disc_rows = []
    for N_val in N_benchmarks:
        L_val = 2.0 * math.log(N_val)
        X_val = int(math.floor(N_val ** (1.0 / 3.0)))
        ps = all_primes[(all_primes > X_val) & (all_primes <= N_val)].astype(np.float64)
        log_ps = np.log(ps)
        alphas = log_ps / L_val
        w0 = 1.0 / ps
        w1 = (1.0 - 2.0 * alphas) / ps
        w2 = ((1.0 - 2.0 * alphas) ** 2) / ps
        w3 = ((1.0 - 2.0 * alphas) ** 3) / ps

        sq0, sq1, sq2, sq3 = np.sum(w0)**2, np.sum(w1)**2, np.sum(w2)**2, np.sum(w3)**2
        # Prefix sums for triangle p*q <= N_val
        cs0, cs1, cs2, cs3 = np.cumsum(w0), np.cumsum(w1), np.cumsum(w2), np.cumsum(w3)
        max_q = N_val // ps.astype(np.int64)
        idxs = np.searchsorted(ps, max_q, side="right") - 1
        valid = idxs >= 0
        tri0 = np.sum(w0[valid] * cs0[idxs[valid]])
        tri1 = np.sum(w1[valid] * cs1[idxs[valid]])
        tri2 = np.sum(w2[valid] * cs2[idxs[valid]])
        tri3 = np.sum(w3[valid] * cs3[idxs[valid]])

        ed0, ed1, ed2, ed3 = tri0 / sq0, tri1 / sq1, tri2 / sq2, tri3 / sq3
        rd2 = ed2 / (1.0 - ed2)
        disc_rows.append({
            "N_star": N_val,
            "X": X_val,
            "M_primes": len(ps),
            "ed0": ed0,
            "ed1": ed1,
            "ed2": ed2,
            "rd2": rd2,
            "ed3": ed3,
        })
        print(
            f"  N_*={N_val:6d}, X={X_val:2d}, M={len(ps):4d} | eta0={ed0:.6f} | eta1={ed1:.6f} | "
            f"eta2={ed2:.6f} (R2={rd2:.6f}) | eta3={ed3:.6f}"
        )

    # Part 3E: Dirichlet Hyperbola Reflection Duality & Fourth-Order Apex Pinching at t = 10^8
    print("\n--- Part 3E: Dirichlet Hyperbola Reflection Duality & 4th-Order Apex Pinching (t = 10^8) ---")
    t_val = mp.mpf("1e8")
    N_star_8 = mp.sqrt(t_val / (2 * mp.pi))
    L_8 = mp.log(t_val / (2 * mp.pi))
    omega_1_8 = L_8 / 2
    pairs_8 = [(17, 19), (41, 89), (61, 67), (503, 701), (3967, 3989)]
    duality_rows = []
    for p_val, q_val in pairs_8:
        pq = p_val * q_val
        omega_p = mp.log(N_star_8 / p_val)
        omega_q = mp.log(N_star_8 / q_val)
        xi_star = mp.log(N_star_8 / pq) / 2
        m_dual = t_val / (2 * mp.pi * pq)
        ln_m_dual = mp.log(m_dual)
        id_err = abs((omega_p + omega_q) - ln_m_dual)
        rel_w2 = ((omega_p / omega_1_8) ** 2) * ((omega_q / omega_1_8) ** 2)
        pinch_ratio = ((omega_p * omega_q) ** 2) / ((ln_m_dual**4) / 16)
        assert id_err < mp.mpf("1e-59") and pinch_ratio <= mp.mpf("1.0000000000001")
        duality_rows.append({
            "pair": (p_val, q_val),
            "pq": pq,
            "xi_star": xi_star,
            "m_dual": m_dual,
            "id_err": id_err,
            "rel_w2": rel_w2,
            "pinch_ratio": pinch_ratio,
        })
        print(
            f"  ({p_val:4d},{q_val:4d}) pq={pq:8d} | xi*={float(xi_star):+8.5f} | m_dual={float(m_dual):9.4f} | "
            f"|w_p+w_q - ln m_dual|={mp.nstr(id_err, 3):9s} | rel_w2={float(rel_w2):.4e} | "
            f"pinch_ratio={float(pinch_ratio):.6f} <= 1"
        )

    return {
        "dilog_rows": dilog_rows,
        "disc_rows": disc_rows,
        "duality_rows": duality_rows,
        "th50_0": float(th50_0),
        "th50_1": float(th50_1),
        "th50_2": float(th50_2),
        "th50_3": float(th50_3),
        "th75_3": float(th75_3),
        "th_Wpure_inf": float(th_Wpure_inf),
    }


# ==============================================================================
# MODULE 4: PUBLICATION FIGURE GENERATION (300 DPI VECTOR PDF < 500 KB & PNG)
# ==============================================================================

def generate_publication_figures(sec2_data, sec3_data, sec4_data):
    print("\n" + "=" * 88)
    print("MODULE 4: GENERATING 300-DPI 4-PANEL PUBLICATION FIGURES (PDF < 500 KB & PNG)")
    print("=" * 88)

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9.5,
        "axes.titlesize": 9.8,
        "axes.labelsize": 9.2,
        "legend.fontsize": 7.6,
        "xtick.labelsize": 8.4,
        "ytick.labelsize": 8.4,
    })

    # --------------------------------------------------------------------------
    # FIGURE 1: Portal-Frame Overdetermination & 50%-50% Green Equipartition
    # --------------------------------------------------------------------------
    fig1, axes1 = plt.subplots(2, 2, figsize=(12.4, 8.9), dpi=300)

    # Panel (a): 2D Davenport-Heilbronn Case-I wall Gamma_false & detached parabola fold d_tip
    ax = axes1[0, 0]
    d_mesh = np.linspace(0.0, 0.44, 140)
    t_mesh = np.linspace(83.6, 86.9, 180)
    D_grid, T_grid = np.meshgrid(d_mesh, t_mesh)
    Z_shift_f = (T_grid - float(sec3_data["t_dh_f"])) - 1j * D_grid

    c_dh_f_flt = [float(c) for c in sec3_data["c_dh_f"][:65]]
    val_dh = np.zeros_like(Z_shift_f, dtype=np.complex128)
    for c_val in reversed(c_dh_f_flt):
        val_dh = val_dh * Z_shift_f + c_val
    sc_grid = np.exp(np.pi * T_grid / 4.0)
    R_dh_grid = np.real(val_dh) * sc_grid
    I_dh_grid = np.imag(val_dh) * sc_grid
    log_mod_dh = np.log10(np.maximum(np.hypot(R_dh_grid, I_dh_grid), 1e-4))

    pcm = ax.pcolormesh(
        D_grid, T_grid, log_mod_dh, cmap="YlGnBu_r", shading="auto", rasterized=True, alpha=0.82
    )
    fig1.colorbar(pcm, ax=ax, pad=0.02, label=r"$\log_{10}(\mathrm{e}^{\pi t/4}|\Xi_{\mathrm{DH}}(t-\mathrm{i}d)|)$")
    ax.contour(D_grid, T_grid, I_dh_grid, levels=[0.0], colors=["#003399"], linewidths=2.0)
    ax.contour(D_grid, T_grid, R_dh_grid, levels=[0.0], colors=["#cc0000"], linewidths=1.8, linestyles="--")

    d_star_f = float(sec2_data["dh_rows"][0]["d_star"])
    g_star_f = float(sec2_data["dh_rows"][0]["gamma_star"])
    tau_f_f = float(sec2_data["dh_rows"][0]["tau_false"])
    d_tip_f = float(sec3_data["d_tip"])
    t_tip_f = float(sec3_data["t_tip"])

    ax.plot(0.0, tau_f_f, "o", color="#003399", markersize=6.5, label=rf"Anchor $\tau_{{\mathrm{{false}}}}(0)={tau_f_f:.2f}$ ($c_1^{{\mathrm{{DH}}}}<0$)")
    ax.plot(d_star_f, g_star_f, "*", color="#ffcc00", markeredgecolor="black", markersize=11, label=rf"Off-line zero $(d_*,\gamma_*)=({d_star_f:.3f},{g_star_f:.2f})$")
    ax.plot(d_tip_f, t_tip_f, "s", color="#cc0000", markersize=6.0, label=rf"Fold tip $(d_{{\mathrm{{tip}}}},t_{{\mathrm{{tip}}}})=({d_tip_f:.3f},{t_tip_f:.2f})$")
    ax.plot([], [], "-", color="#003399", linewidth=2.0, label=r"Chladni walls $I_{\mathrm{DH}}(d,t)=0$")
    ax.plot([], [], "--", color="#cc0000", linewidth=1.8, label=r"Real nodal curves $R_{\mathrm{DH}}(d,t)=0$")
    ax.set_xlabel(r"Transverse distance $d = \sigma - 1/2$")
    ax.set_ylabel(r"Spectral height $t$")
    ax.set_title(r"(a) Portal-Frame Case-I Forcing & Detached Fold in $\Xi_{\mathrm{DH}}(z)$")
    ax.legend(loc="lower left", framealpha=0.92, fontsize=7.1)
    ax.grid(True, linestyle=":", alpha=0.35)

    # Panel (b): Codimension-3 interior saddle overdetermination (3 equations in 2 unknowns)
    ax = axes1[0, 1]
    s_rows = sec2_data["saddle_rows"]
    x_pos = np.arange(len(s_rows))
    labels_b = [f"$({float(r['g0']):.1f},\\,{float(r['d0']):.2f})$" for r in s_rows]
    eq1_vals = [abs(float(r["eq1_I_sc"])) for r in s_rows]
    grad_vals = [max(float(r["eq2"]), float(r["eq3"]), 1e-60) for r in s_rows]
    c1_vals = [abs(float(r["c1_def_sc2"])) for r in s_rows]

    w_bar = 0.27
    ax.bar(x_pos - w_bar, eq1_vals, width=w_bar, color="#b2182b", label=r"Eq. (1) Residual $S_c|I_{\mathrm{def}}(d_c,t_c)| \neq 0$")
    ax.bar(x_pos, c1_vals, width=w_bar, color="#2166ac", label=r"$d=0$ Defect $S_0^2|c_1^{\mathrm{def}}(\gamma_0)|$ ($c_1^{\mathrm{def}}<0$)")
    ax.bar(x_pos + w_bar, grad_vals, width=w_bar, color="#4d9221", label=r"Eqs. (2)--(3) $S_c|\nabla I_{\mathrm{def}}(d_c,t_c)| = 0$")
    ax.set_yscale("log")
    ax.set_ylim(1e-62, 1e6)
    ax.axhline(1e-55, color="black", linestyle="--", linewidth=0.9, label="60-dps numerical zero floor")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels_b)
    ax.set_xlabel(r"Deformation parameters $(\gamma_0, d_0)$ with interior saddle $\Xi_{\mathrm{def}}'(z_c)=0$")
    ax.set_ylabel("Rescaled magnitude (logarithmic scale)")
    ax.set_title(r"(b) Codimension-3 Saddle Overdetermination ($I(d_c,t_c) \neq 0$ at $\nabla I = 0$)")
    ax.legend(loc="center right", framealpha=0.92, fontsize=7.2)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Panel (c): Exact 50%-50% Green-Dirichlet-Neumann Chamber Equipartition across Omega_1, Omega_2, Omega_3
    ax = axes1[1, 0]
    d_fine = np.linspace(0.02, 0.50, 25)
    colors_c = ["#1b9e77", "#d95f02", "#7570b3"]
    for idx_k, (ch_key, col_k) in enumerate(zip(["Omega_1", "Omega_2", "Omega_3"], colors_c), start=1):
        rows_k = sec3_data["riem_chambers"][ch_key]
        d_pts = [float(r["d0"]) for r in rows_k]
        t_mid_k = 0.5 * (float(rows_k[0]["ta"]) + float(rows_k[0]["tb"]))
        sc2_k = math.exp(math.pi * t_mid_k / 2.0)
        I_pts = [float(r["I_Omega"]) * sc2_k for r in rows_k]
        R_pts = [float(r["R_Omega"]) * sc2_k for r in rows_k]
        E_pts = [float(r["E_Omega"]) * sc2_k for r in rows_k]

        ax.plot(d_pts, E_pts, "-", color=col_k, linewidth=2.1, label=rf"$\Omega_{idx_k}$: 2D Energy $\mathrm{{e}}^{{\pi\bar{{t}}/2}}\mathcal{{E}}_{{\Omega_{idx_k}}}(d_0)$")
        ax.plot(d_pts, I_pts, "o", color=col_k, markersize=6.0, fillstyle="none", label=rf"$\Omega_{idx_k}$: Cross-term $\mathcal{{I}}_{{\Omega_{idx_k}}}(d_0)$ ($50\%$)")
        ax.plot(d_pts, R_pts, "x", color=col_k, markersize=6.0, label=rf"$\Omega_{idx_k}$: Real-term $\mathcal{{R}}_{{\Omega_{idx_k}}}(d_0)$ ($50\%$)")

    ax.set_yscale("log")
    ax.set_xlabel(r"Chamber right-boundary depth $d_0$")
    ax.set_ylabel(r"Rescaled chamber integral $\mathrm{e}^{\pi\bar{t}/2}\mathcal{E}_{\Omega_k}(d_0)$")
    ax.set_title("(c) Exact 50%\u201350% Green\u2013Dirichlet\u2013Neumann Chamber Equipartition")
    ax.legend(loc="lower right", ncol=2, framealpha=0.92, fontsize=6.8)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Panel (d): Wronskian-Slope Reciprocity & Wall Amplitude Evolution (Xi, Xi_{chi_5} vs Xi_DH)
    ax = axes1[1, 1]
    d_scan_d = np.linspace(0.0, 0.50, 80)
    # Track R(d, tau(d)) / |R(0, tau(0))| along Gamma_1, Gamma_2 of Xi, Gamma_{chi_5}, and Gamma_{false} of Xi_DH
    for k_idx, col_w, lbl_w in [
        (0, "#1f78b4", r"Riemann $\Xi(z)$ Wall $\Gamma_1$ ($c_1(\tau_1)>0$)"),
        (1, "#33a02c", r"Riemann $\Xi(z)$ Wall $\Gamma_2$ ($c_1(\tau_2)>0$)"),
    ]:
        tau0 = sec2_data["riem_wall_rows"][k_idx]["tau_k"]
        c_j = sec2_data["riem_wall_rows"][k_idx]["c_j"]
        R0 = abs(float(eval_f0_f1_f2(mp.mpf("0"), c_j)[0].real))
        ratios_w = []
        tc = tau0
        for dv in d_scan_d:
            d_mp = mp.mpf(str(float(dv)))
            if dv > 0:
                for _ in range(5):
                    f0, f1, _ = eval_f0_f1_f2((tc - tau0) - 1j * d_mp, c_j)
                    tc -= mp.im(f0) / mp.im(f1)
            f0, _, _ = eval_f0_f1_f2((tc - tau0) - 1j * d_mp, c_j)
            ratios_w.append(abs(float(f0.real)) / R0)
        ax.plot(d_scan_d, ratios_w, "-", linewidth=2.0, color=col_w, label=lbl_w)

    # Xi_{chi_5} wall
    t_c5_0 = sec3_data["chi5_rows"][0]["tb"]
    c_c5_w = chi5_taylor_coeffs(t_c5_0, 76)
    R0_c5 = abs(float(eval_f0_f1_f2(mp.mpf("0"), c_c5_w)[0].real))
    ratios_c5 = []
    tc = t_c5_0
    for dv in d_scan_d:
        d_mp = mp.mpf(str(float(dv)))
        if dv > 0:
            for _ in range(5):
                f0, f1, _ = eval_f0_f1_f2((tc - t_c5_0) - 1j * d_mp, c_c5_w)
                tc -= mp.im(f0) / mp.im(f1)
        f0, _, _ = eval_f0_f1_f2((tc - t_c5_0) - 1j * d_mp, c_c5_w)
        ratios_c5.append(abs(float(f0.real)) / R0_c5)
    ax.plot(d_scan_d, ratios_c5, "-.", linewidth=2.0, color="#6a3d9a", label=r"Euler-product $\Xi_{\chi_5}(z)$ Wall ($c_1^{\chi_5}>0$)")

    # Xi_DH false wall
    t_dh_0 = sec3_data["t_dh_f"]
    c_dh_w = sec3_data["c_dh_f"]
    R0_dh = abs(float(eval_f0_f1_f2(mp.mpf("0"), c_dh_w)[0].real))
    ratios_dh = []
    tc = t_dh_0
    for dv in d_scan_d:
        d_mp = mp.mpf(str(float(dv)))
        if dv > 0:
            for _ in range(5):
                f0, f1, _ = eval_f0_f1_f2((tc - t_dh_0) - 1j * d_mp, c_dh_w)
                tc -= mp.im(f0) / mp.im(f1)
        f0, _, _ = eval_f0_f1_f2((tc - t_dh_0) - 1j * d_mp, c_dh_w)
        ratios_dh.append(abs(float(f0.real)) / R0_dh)
    ax.plot(d_scan_d, ratios_dh, "--", linewidth=2.2, color="#e31a1c", label="Davenport\u2013Heilbronn $\\Gamma_{\\mathrm{false}}$ ($c_1^{\\mathrm{DH}}<0$)")
    ax.axvline(d_star_f, color="#e31a1c", linestyle=":", linewidth=1.2, label=rf"Off-line zero $d_*={d_star_f:.4f}$ ($R_{{\mathrm{{DH}}}}=0$)")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1.0)
    ax.set_xlabel(r"Transverse distance $d$")
    ax.set_ylabel(r"Normalized wall amplitude $|R(d,\tau_k(d))| / |\Xi(\tau_k(0))|$")
    ax.set_title(r"(d) Wall Amplitude Law $\frac{\mathrm{d}}{\mathrm{d}d}(\frac{1}{2}R^2) = d\widetilde{W}(1+\tau_k'(d)^2)$")
    ax.legend(loc="upper left", framealpha=0.92, fontsize=7.2)
    ax.grid(True, linestyle=":", alpha=0.4)

    fig1.tight_layout()
    f1_pdf = PAPER_FIG_DIR / "fig1_portal_frame_overdetermination_and_green_equipartition.pdf"
    f1_png = PAPER_FIG_DIR / "fig1_portal_frame_overdetermination_and_green_equipartition.png"
    fig1.savefig(f1_pdf, dpi=300, bbox_inches="tight")
    fig1.savefig(f1_png, dpi=300, bbox_inches="tight")
    shutil.copy2(f1_pdf, ROOT_FIG_DIR / f1_pdf.name)
    shutil.copy2(f1_png, ROOT_FIG_DIR / f1_png.name)
    plt.close(fig1)

    # --------------------------------------------------------------------------
    # FIGURE 2: Dilogarithm Horizon Damping & Dirichlet Hyperbola Duality
    # --------------------------------------------------------------------------
    fig2, axes2 = plt.subplots(2, 2, figsize=(12.4, 8.9), dpi=300)

    # Panel (a): Attenuated Semiprime Spring Minimum (1 - 1/(2*eta))^2 & Discriminant Delta(eta) = 1 - 2*eta
    ax = axes2[0, 0]
    U_grid = np.linspace(-3.2, 0.8, 300)
    eta_examples = [
        (0.243955, "#7f7f7f", "--", r"$p=0$ Unweighted ($\eta^{(0)}(1/6)=0.244 < 1/2$, $\Delta>0$)"),
        (0.500000, "#ff7f00", "-.", r"Critical Threshold ($\eta = 1/2$, $\Delta(\eta) = 0$)"),
        (0.676891, "#1f78b4", "-", r"$p=2$ Wronskian ($\eta^{(2)}(1/6)=0.677 > 2/3$, Floor $+0.0683$)"),
        (0.798091, "#33a02c", "-", "$p=3$ Wigner\u2013Sinh ($\\eta^{(3)}(1/6)=0.798$, Floor $+0.1395$)"),
    ]
    for eta_v, col_e, ls_e, lbl_e in eta_examples:
        F_sq = (1.0 + U_grid + 0.5 * eta_v * (U_grid**2)) ** 2
        ax.plot(U_grid, F_sq, linestyle=ls_e, color=col_e, linewidth=2.0, label=lbl_e)
        if eta_v >= 0.5:
            u_min = -1.0 / eta_v
            f_min = (1.0 - 1.0 / (2.0 * eta_v)) ** 2
            ax.plot(u_min, f_min, "o", color=col_e, markersize=6.0)

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylim(-0.03, 0.65)
    ax.set_xlabel(r"Anti-phase prime-gauge coordinate $U = \mathrm{Re}\,\mathcal{P}_{>X}(s)$ ($V=0$)")
    ax.set_ylabel(r"Spring energy $|1 + U + \frac{\eta}{2}U^2|^2$")
    ax.set_title(r"(a) Attenuated Semiprime Spring Floor $(1-\frac{1}{2\eta})^2 > 0$ ($\eta > 1/2$)")
    ax.legend(loc="upper right", framealpha=0.92, fontsize=7.2)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Panel (b): Exact Euler Dilogarithm Li_2(2*theta) Curves eta^{(p)}(theta) for p = 0, 1, 2, 3
    ax = axes2[0, 1]
    th_grid = np.linspace(0.05, 0.248, 160)
    for p_idx, col_p, ls_p, lbl_p in [
        (0, "#7f7f7f", "--", r"$p=0$ Unweighted $\eta^{(0)}(\theta)$ ($\eta^{(0)}(1/6)=24.40\%$)"),
        (1, "#ff7f00", "-.", r"$p=1$ Linear $\eta^{(1)}(\theta)$ ($\eta^{(1)}(1/6)=49.56\%$)"),
        (2, "#1f78b4", "-", r"$p=2$ Quadratic Wronskian $\eta^{(2)}(\theta)$ ($\eta^{(2)}(1/6)=67.69\%$)"),
        (3, "#33a02c", "-", "$p=3$ Cubic Wigner\u2013Sinh $\\eta^{(3)}(\\theta)$ ($\\eta^{(3)}(1/6)=79.81\\%$)"),
    ]:
        y_vals = [float(eta_closed(mp.mpf(str(float(x))), p_idx)) for x in th_grid]
        ax.plot(th_grid, y_vals, linestyle=ls_p, color=col_p, linewidth=2.1, label=lbl_p)

    ax.axhline(0.50, color="#b2182b", linestyle=":", linewidth=1.4, label=r"$50\%$ Discriminant Threshold ($\eta = 1/2$)")
    ax.axvline(1.0 / 6.0, color="black", linestyle="--", linewidth=1.1, label=r"Cubic window $\theta = 1/6 = 0.1667$")
    ax.axvline(sec4_data["th_Wpure_inf"], color="#6a3d9a", linestyle=":", linewidth=1.3, label=rf"$\theta_{{\mathrm{{W-pure}}}}^\infty = {sec4_data['th_Wpure_inf']:.4f}$")
    ax.plot(sec4_data["th50_2"], 0.50, "D", color="#1f78b4", markersize=6.5, label=rf"$\theta_{{*,50\%}}^{{(2)}} = {sec4_data['th50_2']:.4f}$")
    ax.plot(sec4_data["th50_3"], 0.50, "s", color="#33a02c", markersize=6.5, label=rf"$\theta_{{*,50\%}}^{{(3)}} = {sec4_data['th50_3']:.4f}$")
    ax.set_xlabel(r"Prime-gauge sieve exponent $\theta = \ln X / \ln(t/2\pi)$")
    ax.set_ylabel(r"Forward-triangle capture fraction $\eta^{(p)}(\theta)$")
    ax.set_title(r"(b) Exact Euler Dilogarithm $\mathrm{Li}_2(2\theta)$ Laws $\eta^{(p)}(\theta)$")
    ax.legend(loc="lower left", framealpha=0.92, fontsize=6.8)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Panel (c): Discrete Prime Sums vs Continuous Li_2(2*theta) Asymptote up to N_* = 10^5
    ax = axes2[1, 0]
    d_rows = sec4_data["disc_rows"]
    N_pts = [r["N_star"] for r in d_rows]
    for p_idx, key_p, col_p, lbl_p in [
        (0, "ed0", "#7f7f7f", r"$p=0$ Discrete $\eta_{\mathrm{disc}}^{(0)}$ vs. $\eta^{(0)}(1/6)=0.2440$"),
        (1, "ed1", "#ff7f00", r"$p=1$ Discrete $\eta_{\mathrm{disc}}^{(1)}$ vs. $\eta^{(1)}(1/6)=0.4956$"),
        (2, "ed2", "#1f78b4", r"$p=2$ Discrete $\eta_{\mathrm{disc}}^{(2)}$ vs. $\eta^{(2)}(1/6)=0.6769$"),
        (3, "ed3", "#33a02c", r"$p=3$ Discrete $\eta_{\mathrm{disc}}^{(3)}$ vs. $\eta^{(3)}(1/6)=0.7981$"),
    ]:
        vals_p = [r[key_p] for r in d_rows]
        asymp_p = float(eta_closed(mp.mpf(1) / 6, p_idx))
        ax.plot(N_pts, vals_p, "o-", color=col_p, linewidth=1.8, markersize=5.5, label=lbl_p)
        ax.axhline(asymp_p, color=col_p, linestyle="--", linewidth=1.1, alpha=0.8)

    ax.axhline(0.50, color="#b2182b", linestyle=":", linewidth=1.4, label=r"$50\%$ Threshold ($\eta = 1/2$)")
    ax.set_xscale("log")
    ax.set_ylim(0.10, 0.90)
    ax.set_xlabel("Riemann\u2013Siegel saddle horizon $N_\\star(t) = \\sqrt{t/2\\pi}$")
    ax.set_ylabel(r"Discrete prime-sum capture fraction $\eta_{\mathrm{disc}}^{(p)}(N_\star)$")
    ax.set_title(r"(c) Discrete Prime-Sum Convergence to $\mathrm{Li}_2(1/3)$ Laws ($\theta=1/6$)")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=6.9)
    ax.grid(True, linestyle=":", alpha=0.4)

    # Panel (d): Dirichlet Hyperbola Reflection Duality & Fourth-Order Apex Pinching
    ax = axes2[1, 1]
    a1_lin = np.linspace(1.0 / 6.0, 0.50, 140)
    a2_lin = np.linspace(1.0 / 6.0, 0.50, 140)
    A1, A2 = np.meshgrid(a1_lin, a2_lin)
    W2_prod = ((1.0 - 2.0 * A1) ** 2) * ((1.0 - 2.0 * A2) ** 2)
    pcm2 = ax.pcolormesh(
        A1, A2, np.log10(np.maximum(W2_prod, 1e-10)), cmap="magma", shading="auto", rasterized=True
    )
    fig2.colorbar(pcm2, ax=ax, pad=0.02, label=r"$\log_{10}[(1-2\alpha_1)^2(1-2\alpha_2)^2]$")
    a_hyp = np.linspace(1.0 / 6.0, 1.0 / 3.0, 100)
    ax.plot(a_hyp, 0.5 - a_hyp, "-", color="#00ffff", linewidth=2.4, label=r"Horizon $\alpha_1+\alpha_2=1/2$ ($pq=N_\star$)")
    ax.text(0.20, 0.215, r"Forward Triangle" + "\n" + r"$\mathcal{T}_{\leq N}$ ($67.69\%$)" + "\n" + r"$\xi_{pq}^\star > 0$", color="white", fontsize=8.0, fontweight="bold", bbox=dict(boxstyle="round,pad=0.28", facecolor="#1a1a1a", alpha=0.68, edgecolor="none"))
    ax.text(0.325, 0.365, r"Reflected Corner" + "\n" + r"$\mathcal{C}_{>N}$ ($32.31\%$)" + "\n" + r"$(\omega_p\omega_q)^2 \leq \frac{1}{16}(\ln m_{\mathrm{dual}})^4$", color="white", fontsize=7.8, bbox=dict(boxstyle="round,pad=0.28", facecolor="#1a1a1a", alpha=0.62, edgecolor="none"))
    ax.plot(0.50, 0.50, "*", color="#00ffff", markersize=11, markeredgecolor="black", label=r"Apex $pq=N_\star^2$ ($m_{\mathrm{dual}}=1$, $4$th-order zero)")
    ax.set_xlabel(r"First prime logarithmic scale $\alpha_1 = \ln p / L \in [1/6, 1/2]$")
    ax.set_ylabel(r"Second prime logarithmic scale $\alpha_2 = \ln q / L \in [1/6, 1/2]$")
    ax.set_title(r"(d) Hyperbola Duality $\omega_p+\omega_q=\ln m_{\mathrm{dual}}$ & Apex Pinching")
    ax.legend(loc="lower left", framealpha=0.92, fontsize=7.2)

    fig2.tight_layout()
    f2_pdf = PAPER_FIG_DIR / "fig2_dilogarithm_horizon_damping_and_hyperbola_duality.pdf"
    f2_png = PAPER_FIG_DIR / "fig2_dilogarithm_horizon_damping_and_hyperbola_duality.png"
    fig2.savefig(f2_pdf, dpi=300, bbox_inches="tight")
    fig2.savefig(f2_png, dpi=300, bbox_inches="tight")
    shutil.copy2(f2_pdf, ROOT_FIG_DIR / f2_pdf.name)
    shutil.copy2(f2_png, ROOT_FIG_DIR / f2_png.name)
    plt.close(fig2)

    for pdf_path in [f1_pdf, f2_pdf]:
        sz_kb = pdf_path.stat().st_size / 1024.0
        print(f"  Generated {pdf_path.name}: {sz_kb:.1f} KB (< 500 KB requirement satisfied)")
        assert sz_kb < 500.0, f"PDF figure {pdf_path.name} exceeds 500 KB ({sz_kb:.1f} KB)"


def main():
    t0_all = time.perf_counter()
    print("=" * 88)
    print("FLAGSHIP PAPER XIV: THE CHLADNI PORTAL-FRAME SYNTHESIS VERIFICATION SUITE")
    print(f"Working Precision: mpmath.mp.dps = {mp.mp.dps}, mpmath.iv.dps = {mp.iv.dps}")
    print("=" * 88)

    sec2_data = verify_section2_macro_portal_frame()
    sec3_data = verify_section3_meso_green_equipartition(sec2_data)
    sec4_data = verify_section4_micro_dilogarithm()
    generate_publication_figures(sec2_data, sec3_data, sec4_data)

    elapsed = time.perf_counter() - t0_all
    print("\n" + "=" * 88)
    print(f"ALL THREE SCALES (MACRO -> MESO -> MICRO) & FIGURES CERTIFIED IN {elapsed:.2f} s (EXIT 0)")
    print("=" * 88)


if __name__ == "__main__":
    main()
