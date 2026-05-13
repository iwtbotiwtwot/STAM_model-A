#!/usr/bin/env python3
"""
G110_axial_QNM_higher_modes.py

Extension of G108 / G109 to higher modes:
   ell = 2, 3, 4  (fundamentals)
   n   = 0, 1     (fundamental + first overtone)

Setup
=====
Same V_exact = V_geom + (sqrt f)'' / sqrt f as G108, with the analytic
(sqrt f)'' computed by chain rule.  Same convergence-tested grid choice:
    rs_min = -500, N = 11001    (G109 row 2; eps_min = 4.3e-4)

For each ell:
    1. Time-domain evolve V_GR(ell) and V_exact(ell) with an early-time
       narrow pulse that excites overtones.
    2. Calibrate against Leaver gold for (ell, n=0) and (ell, n=1).
    3. Use a TWO-MODE damped-sinusoid fit on an early-time window to
       extract both n=0 and n=1 simultaneously.
    4. Report omega_GR and omega_exact for each (ell, n), and their shifts.

Leaver gold reference (Schwarzschild axial, M = 1)
==================================================
    (ell, n)  omega
    (2, 0)    0.373672 - 0.088962 i
    (2, 1)    0.346711 - 0.273915 i
    (3, 0)    0.599443 - 0.092703 i
    (3, 1)    0.582644 - 0.281298 i
    (4, 0)    0.809178 - 0.094164 i
    (4, 1)    0.796632 - 0.284334 i

Outputs
=======
    results/G110_higher_modes_summary.md
    plots/G110_higher_modes.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# ===========================================================================
# Constants and Leaver gold table
# ===========================================================================

M = 1.0
CALIBRATION_THRESHOLD = 0.02      # 2% per mode (n=1 is intrinsically harder)

LEAVER_GOLD = {
    (2, 0): complex(0.373672, -0.088962),
    (2, 1): complex(0.346711, -0.273915),
    (3, 0): complex(0.599443, -0.092703),
    (3, 1): complex(0.582644, -0.281298),
    (4, 0): complex(0.809178, -0.094164),
    (4, 1): complex(0.796632, -0.284334),
}


# ===========================================================================
# Metric, f(Sigma), V_geom, V_exact   (copied from G108/G109)
# ===========================================================================

def h_fn(r): return 1.0 - 2.0 * M / r
def hp_fn(r): return 2.0 * M / r**2
def y_of_r(r): return 6.0 * M / r - 2.0
def F_quintic(y): return 1.0 - 5.0 * y**4 + 4.0 * y**5
def F_first(y): return -20.0 * y**3 + 20.0 * y**4

def k_STAM(r):
    y = y_of_r(r)
    return h_fn(r) * F_quintic(y) if y > 0 else h_fn(r)

def kp_STAM(r):
    y = y_of_r(r)
    if y > 0:
        return hp_fn(r) * F_quintic(y) + h_fn(r) * F_first(y) * (-6.0 * M / r**2)
    return hp_fn(r)


def dlnf_dA_numeric(A):
    if A <= 2.0 / 3.0:
        return 0.0
    y = 3.0 * A - 2.0
    F = 1.0 - 5.0 * y**4 + 4.0 * y**5
    return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)

def d2lnf_dA2_numeric(A, hstep=1e-6):
    return (dlnf_dA_numeric(A + hstep) - dlnf_dA_numeric(A - hstep)) / (2.0 * hstep)


def build_f_grid(r_grid_STAM):
    ln_f = np.zeros_like(r_grid_STAM)
    inside = r_grid_STAM < 3.0 * M
    if not np.any(inside):
        return np.ones_like(r_grid_STAM)
    r_inside = r_grid_STAM[inside]
    r_descending = np.sort(r_inside)[::-1]
    def rhs(rv, yv):
        A = 2.0 * M / rv
        return [dlnf_dA_numeric(A) * (-2.0 * M / rv**2)]
    sol = solve_ivp(rhs, [3.0 * M - 1e-12, r_descending[-1]], [0.0],
                    t_eval=r_descending, method="RK45",
                    rtol=1e-12, atol=1e-14, max_step=0.005)
    ln_f_inside = np.empty_like(r_inside)
    order = np.argsort(r_inside)
    ln_f_inside[order] = sol.y[0][::-1]
    ln_f[inside] = ln_f_inside
    return np.exp(ln_f)


def V_geom_grid(r_grid, ell):
    out = np.empty_like(r_grid)
    for i, rv in enumerate(r_grid):
        h_val = h_fn(rv)
        k_val = k_STAM(rv)
        out[i] = h_val * (
            ell * (ell + 1) / rv**2
            - 2.0 * (1.0 - k_val) / rv**2
            - kp_STAM(rv) / (2.0 * rv)
            - k_val * hp_fn(rv) / (2.0 * h_val * rv)
        )
    return out


def analytic_correction(r_grid_STAM):
    out = np.zeros_like(r_grid_STAM)
    for i, rv in enumerate(r_grid_STAM):
        if rv >= 3.0 * M:
            out[i] = 0.0
            continue
        A = 2.0 * M / rv
        dA_dr = -2.0 * M / rv**2
        d2A_dr2 = 4.0 * M / rv**3
        u = dlnf_dA_numeric(A)
        up = d2lnf_dA2_numeric(A)
        dlnf_dr = u * dA_dr
        d_dlnf_dr = up * (dA_dr**2) + u * d2A_dr2
        hk = h_fn(rv) * k_STAM(rv)
        if hk <= 0:
            out[i] = 0.0
            continue
        sqrt_hk = np.sqrt(hk)
        d_hk = hp_fn(rv) * k_STAM(rv) + h_fn(rv) * kp_STAM(rv)
        P = sqrt_hk * dlnf_dr
        d_sqrt_hk_dr = 0.5 * d_hk / sqrt_hk
        dP_dr = d_sqrt_hk_dr * dlnf_dr + sqrt_hk * d_dlnf_dr
        P_star = sqrt_hk * dP_dr
        out[i] = 0.5 * P_star + 0.25 * P**2
    return out


def r_of_rstar_GR(rstar):
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset
    def f(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target
    if rstar > 0:
        return brentq(f, 3.0 * M - 1e-10, 1e8)
    if rstar == 0:
        return 3.0 * M
    if rstar < -40:
        arg = (rstar + 3.0 * M + 2.0 * M * np.log(0.5) - 2.0 * M) / (2.0 * M)
        r_guess = 2.0 * M + 2.0 * M * np.exp(arg)
        r_lo = max(2.0 * M + 1e-300, r_guess * 1e-3)
        r_hi = min(3.0 * M - 1e-12, max(r_guess * 1e3, 2.0 * M + 1e-6))
        try:
            return brentq(f, r_lo, r_hi)
        except ValueError:
            return r_guess
    return brentq(f, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


def build_r_of_rstar_STAM(rstar_grid):
    def rhs(rs, y):
        rv = y[0]
        return [np.sqrt(h_fn(rv) * k_STAM(rv))]
    rs_pos = rstar_grid[rstar_grid > 0]
    rs_neg = rstar_grid[rstar_grid < 0]
    r_arr = np.zeros_like(rstar_grid)
    if len(rs_pos) > 0:
        order = np.argsort(rs_pos)
        x = rs_pos[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=1.0)
        vals = np.zeros_like(rs_pos)
        vals[order] = sol.y[0]
        r_arr[rstar_grid > 0] = vals
    if len(rs_neg) > 0:
        order = np.argsort(rs_neg)[::-1]
        x = rs_neg[order]
        sol = solve_ivp(rhs, [0, x[-1]], [3.0 * M], t_eval=x,
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.02)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar_grid < 0] = vals
    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


def evolve_TD(V_grid, dr_star, T_final, dt, psi0, i_obs):
    psi_old = psi0.copy()
    psi = psi0.copy()
    N_steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star
    t_arr = np.zeros(N_steps)
    sig_arr = np.zeros(N_steps)
    for step in range(N_steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2
        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1])
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])
        psi_old = psi
        psi = psi_new
        t_arr[step] = (step + 1) * dt
        sig_arr[step] = psi[i_obs]
    return t_arr, sig_arr


# ===========================================================================
# Single-mode and two-mode fitting
# ===========================================================================

def one_mode(t, A, om_r, om_i, ph, off):
    tau = t - t[0]
    return A * np.exp(-om_i * tau) * np.cos(om_r * tau + ph) + off


def two_mode(t, A0, om0r, om0i, ph0, A1, om1r, om1i, ph1, off):
    tau = t - t[0]
    return (A0 * np.exp(-om0i * tau) * np.cos(om0r * tau + ph0)
          + A1 * np.exp(-om1i * tau) * np.cos(om1r * tau + ph1) + off)


def fit_one_mode(t, sig, t0, t1, om_guess):
    mask = (t >= t0) & (t <= t1)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 10:
        return None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, om_guess.real, abs(om_guess.imag), 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 5.0, 5.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(one_mode, tf, sf, p0=p0, bounds=bounds, maxfev=40000)
        return complex(popt[1], -popt[2])
    except Exception:
        return None


def fit_two_modes(t, sig, t0, t1, om0_guess, om1_guess):
    mask = (t >= t0) & (t <= t1)
    tf, sf = t[mask], sig[mask]
    if len(tf) < 20:
        return None, None
    A_g = 0.5 * (np.max(sf) - np.min(sf))
    p0 = [A_g, om0_guess.real, abs(om0_guess.imag), 0.0,
          A_g, om1_guess.real, abs(om1_guess.imag), 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 5.0, 5.0, np.pi, np.inf, 5.0, 5.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(two_mode, tf, sf, p0=p0, bounds=bounds, maxfev=80000)
        om0 = complex(popt[1], -popt[2])
        om1 = complex(popt[5], -popt[6])
        # By Im(omega) magnitude: smaller |Im| is n=0, larger |Im| is n=1
        if abs(om0.imag) > abs(om1.imag):
            om0, om1 = om1, om0
        return om0, om1
    except Exception:
        return None, None


# ===========================================================================
# Main pipeline
# ===========================================================================

def run_for_ell(ell, rstar, r_grid_GR, r_grid_STAM, T_final=350.0):
    """Run one ell value through TD on V_GR and V_exact.  Return (t, sig_GR, sig_exact)."""
    dr_star = rstar[1] - rstar[0]

    V_GR_grid = h_fn(r_grid_GR) * (ell * (ell + 1) / r_grid_GR**2
                                    - 6.0 * M / r_grid_GR**3)
    V_g = V_geom_grid(r_grid_STAM, ell)
    corr = analytic_correction(r_grid_STAM)
    V_exact = V_g + corr

    # Pulse choice: narrow + close to PS to excite overtones.
    # For higher ell (sharper barrier), the WKB peak shifts; broader pulse OK.
    pulse_center = 15.0
    pulse_sigma = 2.0 if ell <= 3 else 1.5
    psi0 = np.exp(-(rstar - pulse_center)**2 / (2.0 * pulse_sigma**2))
    i_obs = int(np.argmin(np.abs(rstar - 40.0)))
    dt = 0.5 * dr_star

    t_arr_GR, sig_GR = evolve_TD(V_GR_grid, dr_star, T_final, dt, psi0, i_obs)
    t_arr_e, sig_e = evolve_TD(V_exact, dr_star, T_final, dt, psi0, i_obs)
    # Use same t_arr (they're identical structurally)
    return t_arr_GR, sig_GR, sig_e


def main():
    print("=" * 88)
    print("G110  --  higher-mode QNM under V_exact analytic action correction")
    print("=" * 88)
    print()

    # Grid: G109 row 2 (rs_min = -500, N = 11001).  Converged for n=0; for n=1
    # we may need to revisit if the higher decay rate sampling is inadequate.
    rs_min, rs_max, N = -500.0, 300.0, 11001
    rstar = np.linspace(rs_min, rs_max, N)
    dr_star = rstar[1] - rstar[0]
    print(f"r* grid: N = {N}, range [{rs_min}, {rs_max}], dr* = {dr_star:.4e}")
    print()

    print("Building r_GR(r*), r_STAM(r*), f(Sigma) ...", flush=True)
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rstar])
    r_grid_STAM = build_r_of_rstar_STAM(rstar)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"  STAM eps_min = {eps_min:.4e}")
    print()

    # Loop over ell, run TD, extract n=0 and n=1 from each
    results = []
    for ell in [2, 3, 4]:
        print("-" * 88)
        print(f"ell = {ell}")
        print("-" * 88, flush=True)
        t_arr, sig_GR, sig_exact = run_for_ell(ell, rstar, r_grid_GR, r_grid_STAM)

        # Single-mode fit for n=0 calibration (late-time window)
        om_GR_n0 = fit_one_mode(t_arr, sig_GR, 80.0, 250.0, LEAVER_GOLD[(ell, 0)])
        if om_GR_n0 is None:
            print("  one-mode fit (n=0) failed")
            continue
        err_GR_n0 = abs(om_GR_n0 - LEAVER_GOLD[(ell, 0)]) / abs(LEAVER_GOLD[(ell, 0)])
        print(f"  n=0 GR calib: omega = {om_GR_n0},  Leaver = {LEAVER_GOLD[(ell, 0)]}"
              f",  err = {err_GR_n0*100:.3f}%")

        # Two-mode fit on EARLY-TIME window for n=0 + n=1
        # n=1 decays with Im(om) ~ 0.27-0.28, so by t=30 amplitude is exp(-30*0.28) ~ 2e-4
        # We want a window where n=1 still has measurable amplitude relative to n=0.
        # Use early window t in [10, 50]:
        om_GR_n0_tm, om_GR_n1 = fit_two_modes(
            t_arr, sig_GR, 10.0, 60.0,
            LEAVER_GOLD[(ell, 0)], LEAVER_GOLD[(ell, 1)])
        if om_GR_n1 is None:
            print(f"  two-mode fit (GR) failed for ell={ell}")
            err_GR_n1 = None
        else:
            err_GR_n1 = abs(om_GR_n1 - LEAVER_GOLD[(ell, 1)]) / abs(LEAVER_GOLD[(ell, 1)])
            print(f"  n=1 GR calib: omega = {om_GR_n1},  Leaver = {LEAVER_GOLD[(ell, 1)]}"
                  f",  err = {err_GR_n1*100:.3f}%")

        # V_exact extractions
        om_exact_n0 = fit_one_mode(t_arr, sig_exact, 80.0, 250.0, LEAVER_GOLD[(ell, 0)])
        if om_exact_n0 is None:
            print("  one-mode fit (n=0, V_exact) failed")
            continue
        shift_n0 = abs(om_exact_n0 - om_GR_n0) / abs(om_GR_n0)
        print(f"  n=0 V_exact: omega = {om_exact_n0},  shift = {shift_n0*100:.4f}%")

        om_exact_n0_tm, om_exact_n1 = fit_two_modes(
            t_arr, sig_exact, 10.0, 60.0,
            LEAVER_GOLD[(ell, 0)], LEAVER_GOLD[(ell, 1)])
        if om_exact_n1 is None:
            shift_n1 = None
            print("  two-mode fit (V_exact) failed")
        else:
            # Use the two-mode n=0 from V_exact paired with two-mode n=0 from GR
            # for the n=1 comparison (consistent extraction method)
            if om_GR_n1 is not None:
                shift_n1 = abs(om_exact_n1 - om_GR_n1) / abs(om_GR_n1)
                print(f"  n=1 V_exact: omega = {om_exact_n1},  shift = {shift_n1*100:.4f}%")
            else:
                shift_n1 = None
                print(f"  n=1 V_exact: omega = {om_exact_n1}  (no GR n=1 comparison)")

        results.append({
            "ell": ell,
            "om_GR_n0": om_GR_n0, "om_exact_n0": om_exact_n0,
            "shift_n0": shift_n0,
            "om_GR_n0_tm": om_GR_n0_tm,
            "om_GR_n1": om_GR_n1, "om_exact_n1": om_exact_n1,
            "shift_n1": shift_n1,
            "err_GR_n0": err_GR_n0, "err_GR_n1": err_GR_n1,
            "t": t_arr, "sig_GR": sig_GR, "sig_exact": sig_exact,
        })
        print()

    # ----- Summary -----
    print("=" * 88)
    print("Summary table")
    print("=" * 88)
    print()
    print(f"  {'ell':>4}  {'n':>2}  {'calib err':>10}  "
          f"{'omega_GR':>26}  {'omega_exact':>26}  {'shift %':>10}")
    print("-" * 90)
    for r in results:
        ell = r["ell"]
        print(f"  {ell:>4}  {'0':>2}  {r['err_GR_n0']*100:>9.3f}%  "
              f"{str(r['om_GR_n0']):>26}  {str(r['om_exact_n0']):>26}  "
              f"{r['shift_n0']*100:>9.4f}%")
        if r['om_GR_n1'] is not None and r['shift_n1'] is not None:
            print(f"  {ell:>4}  {'1':>2}  "
                  f"{(r['err_GR_n1'] or 0)*100:>9.3f}%  "
                  f"{str(r['om_GR_n1']):>26}  {str(r['om_exact_n1']):>26}  "
                  f"{r['shift_n1']*100:>9.4f}%")
        elif r['om_GR_n1'] is not None:
            print(f"  {ell:>4}  {'1':>2}  "
                  f"{(r['err_GR_n1'] or 0)*100:>9.3f}%  "
                  f"{str(r['om_GR_n1']):>26}  {'no shift':>26}  ?")
        else:
            print(f"  {ell:>4}  {'1':>2}  n=1 extraction FAILED")
    print()

    # ----- Plot -----
    n_results = len(results)
    if n_results > 0:
        fig, axes = plt.subplots(2, n_results, figsize=(5*n_results, 8))
        if n_results == 1:
            axes = axes[:, np.newaxis]
        for j, r in enumerate(results):
            ax = axes[0, j]
            ax.plot(r["t"], np.abs(r["sig_GR"]) + 1e-30, "tab:blue", linewidth=0.8, label="GR")
            ax.plot(r["t"], np.abs(r["sig_exact"]) + 1e-30, "tab:green", linewidth=0.8, label="V_exact")
            ax.set_yscale("log")
            ax.set_xlabel("t / M")
            ax.set_ylabel("|psi|")
            ax.set_title(f"ell = {r['ell']}")
            ax.set_xlim(0, 250)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            # Highlight fit windows
            ax.axvspan(10, 60, alpha=0.1, color="orange", label="n=1 fit window")
            ax.axvspan(80, 250, alpha=0.08, color="blue", label="n=0 fit window")

            ax = axes[1, j]
            # Bar plot of omegas
            n_labels = []
            n0_re, n0_im, e0_re, e0_im = [], [], [], []
            n1_re, n1_im, e1_re, e1_im = [], [], [], []
            if r["om_GR_n0"] is not None:
                n_labels.append("n=0")
                n0_re.append(LEAVER_GOLD[(r["ell"], 0)].real)
                n0_im.append(abs(LEAVER_GOLD[(r["ell"], 0)].imag))
                e0_re.append(r["om_exact_n0"].real)
                e0_im.append(abs(r["om_exact_n0"].imag))
            if r["om_GR_n1"] is not None and r["om_exact_n1"] is not None:
                n_labels.append("n=1")
                n0_re.append(LEAVER_GOLD[(r["ell"], 1)].real)
                n0_im.append(abs(LEAVER_GOLD[(r["ell"], 1)].imag))
                e0_re.append(r["om_exact_n1"].real)
                e0_im.append(abs(r["om_exact_n1"].imag))
            x = np.arange(len(n_labels))
            width = 0.35
            ax.bar(x - width/2, n0_re, width, label="GR Re(omega)", color="tab:blue")
            ax.bar(x + width/2, e0_re, width, label="V_exact Re(omega)", color="tab:green")
            ax.set_xticks(x)
            ax.set_xticklabels(n_labels)
            ax.set_ylabel("Re(omega)")
            ax.set_title(f"Re(omega) shifts, ell={r['ell']}")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = PLOTS / "G110_higher_modes.png"
        plt.savefig(out_png, dpi=160)
        plt.close()
        print(f"Plot saved: {out_png}")

    # ----- Markdown -----
    md = ["# G110 - Higher-mode axial QNM under V_exact analytic action correction\n"]
    md.append("\nExtends G108/G109 to ell = 3, 4 fundamentals and n = 1 first overtone.\n")
    md.append(f"\n- Grid: rs_min = {rs_min}, N = {N}, dr* = {dr_star:.4e}, eps_min = {eps_min:.4e}\n")
    md.append("- Fit method: late-time single-mode for n=0; early-time two-mode for n=1.\n")
    md.append("\n## Results\n")
    md.append("| ell | n | omega_GR_fit | Leaver | calib % | omega_exact | shift |\n")
    md.append("|---:|---:|---|---|---:|---|---:|\n")
    for r in results:
        ell = r["ell"]
        md.append(f"| {ell} | 0 | {r['om_GR_n0']} | {LEAVER_GOLD[(ell, 0)]} | "
                  f"{r['err_GR_n0']*100:.3f}% | {r['om_exact_n0']} | "
                  f"{r['shift_n0']*100:.4f}% |\n")
        if r['om_GR_n1'] is not None and r['shift_n1'] is not None:
            md.append(f"| {ell} | 1 | {r['om_GR_n1']} | {LEAVER_GOLD[(ell, 1)]} | "
                      f"{(r['err_GR_n1'] or 0)*100:.3f}% | {r['om_exact_n1']} | "
                      f"{r['shift_n1']*100:.4f}% |\n")
        else:
            md.append(f"| {ell} | 1 | --- | {LEAVER_GOLD[(ell, 1)]} | n=1 extraction failed | --- | --- |\n")
    md.append("\n## Files\n")
    md.append("- [scripts/G110_axial_QNM_higher_modes.py](../scripts/G110_axial_QNM_higher_modes.py)\n")
    md.append("- [plots/G110_higher_modes.png](../plots/G110_higher_modes.png)\n")
    out_md = RESULTS / "G110_higher_modes_summary.md"
    out_md.write_text("".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
