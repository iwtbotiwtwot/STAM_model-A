#!/usr/bin/env python3
"""
G91_time_domain_rstar.py

Time-domain QNM extraction in the proper tortoise coordinate, with
Sommerfeld absorbing boundary conditions.

Method  (the fixes from G88's 56% calibration failure)
======================================================
1. Build r*_GR(r) analytically:  r*_GR = r + 2M ln(r/2M - 1)
2. Build r*_STAM(r) numerically by ODE integration:  dr/dr* = sqrt(h k)
3. Use a UNIFORM r* grid (the wave equation is simplest in r*).
4. Inner boundary FAR DOWN the r* throat -- for STAM, this means very
   negative r* (the STAM throat is power-law-deep, not log-deep).
5. Sommerfeld BCs:
       Right edge:  (∂_t + ∂_r*) ψ = 0   (outgoing)
       Left edge:   (∂_t − ∂_r*) ψ = 0   (ingoing toward horizon)
6. Calibrate Schwarzschild to < 1% vs Leaver gold BEFORE running STAM.

Wave equation in r* (clean, separable):
    ∂²ψ/∂t² − ∂²ψ/∂r*² + V(r(r*)) ψ = 0
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


M = 1.0


# ---------------------------------------------------------------------------
# Metric functions
# ---------------------------------------------------------------------------

def h_fn(r):
    return 1.0 - 2.0 * M / r


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def y_of_r(r):
    return 6.0 * M / r - 2.0


def k_GR(r):
    return h_fn(r)


def k_STAM(r):
    y = y_of_r(r)
    h_val = h_fn(r)
    if y > 0:
        return h_val * F_quintic(y)
    return h_val


def V_RW(r, ell, k_fn):
    h_val = h_fn(r)
    k_val = k_fn(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


# ---------------------------------------------------------------------------
# Tortoise coordinate construction  (anchor: r*(3M) = 0)
# ---------------------------------------------------------------------------

def r_of_rstar_GR(rstar):
    """Invert r* = r + 2M ln(r/2M − 1) − [3M + 2M ln(1/2)] for r."""
    from scipy.optimize import brentq
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def f(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    # f(r) is monotonic in r > 2M, asymptotic to -∞ at r=2M+ and +∞ at r=∞
    if rstar > 0:
        return brentq(f, 3.0 * M - 1e-10, 1e6)
    else:
        # rstar < 0:  r ∈ (2M, 3M)
        return brentq(f, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


def build_r_of_rstar_STAM(rstar_grid):
    """Build r(r*) for the STAM metric via ODE integration anchored at
       r*(3M) = 0.
    """
    def rhs(rs, y):
        r_val = y[0]
        return [np.sqrt(h_fn(r_val) * k_STAM(r_val))]

    rs_pos = rstar_grid[rstar_grid > 0]
    rs_neg = rstar_grid[rstar_grid < 0]
    rs_zero = rstar_grid[rstar_grid == 0]

    r_arr = np.zeros_like(rstar_grid)

    if len(rs_pos) > 0:
        # Sort ascending and integrate outward
        sort_pos = np.argsort(rs_pos)
        rs_pos_sorted = rs_pos[sort_pos]
        sol = solve_ivp(rhs, [0, rs_pos_sorted[-1]], [3.0 * M],
                        t_eval=rs_pos_sorted, method='RK45',
                        rtol=1e-12, atol=1e-14, max_step=1.0)
        # Map back to original positions
        r_pos = np.zeros_like(rs_pos)
        r_pos[sort_pos] = sol.y[0]
        r_arr[rstar_grid > 0] = r_pos

    if len(rs_neg) > 0:
        # Sort DESCENDING (closest to 0 first) and integrate inward
        sort_neg = np.argsort(rs_neg)[::-1]
        rs_neg_sorted = rs_neg[sort_neg]
        sol = solve_ivp(rhs, [0, rs_neg_sorted[-1]], [3.0 * M],
                        t_eval=rs_neg_sorted, method='RK45',
                        rtol=1e-12, atol=1e-14, max_step=0.05)
        r_neg = np.zeros_like(rs_neg)
        r_neg[sort_neg] = sol.y[0]
        r_arr[rstar_grid < 0] = r_neg

    if len(rs_zero) > 0:
        r_arr[rstar_grid == 0] = 3.0 * M

    return r_arr


# ---------------------------------------------------------------------------
# Time-domain evolution with Sommerfeld BCs
# ---------------------------------------------------------------------------

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, psidot0=None,
              i_obs=None, save_snapshots=False):
    """Wave eq:  ∂²ψ/∂t² = ∂²ψ/∂r*² − V ψ
       Sommerfeld BCs:
          left  (r* → −∞):  (∂_t − ∂_r*) ψ = 0
          right (r* → +∞):  (∂_t + ∂_r*) ψ = 0
    """
    N = len(V_grid)
    if psidot0 is None:
        psidot0 = np.zeros_like(psi0)

    # For leapfrog, need ψ at two consecutive times
    psi_old = psi0 - dt * psidot0  # ψ(−dt)
    psi = psi0.copy()                # ψ(0)

    N_steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star

    t_arr = np.zeros(N_steps)
    sig_arr = np.zeros(N_steps)

    snapshots = []
    snapshot_times = []

    for step in range(N_steps):
        # Compute ∂²ψ/∂r*² via central differences (interior)
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2

        # Wave eq update (interior)
        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1])

        # Sommerfeld BCs (first-order outflow)
        # Right edge (outgoing in +r*):
        #   (∂_t + ∂_r*) ψ = 0  =>  ψ[N-1]_{n+1} = ψ[N-1]_n - (dt/dr*)(ψ[N-1]_n - ψ[N-2]_n)
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        # Left edge (outgoing in -r* = ingoing toward horizon):
        #   (∂_t − ∂_r*) ψ = 0  =>  ψ[0]_{n+1} = ψ[0]_n - (dt/dr*)(ψ[1]_n - ψ[0]_n)
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])

        psi_old = psi
        psi = psi_new

        t_arr[step] = (step + 1) * dt
        if i_obs is not None:
            sig_arr[step] = psi[i_obs]
        if save_snapshots and (step % (N_steps // 10) == 0):
            snapshots.append(psi.copy())
            snapshot_times.append(t_arr[step])

    return t_arr, sig_arr, snapshots, snapshot_times


# ---------------------------------------------------------------------------
# QNM extraction
# ---------------------------------------------------------------------------

def damped_sinusoid(t, A, omega_r, omega_i, phi, offset):
    return A * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + offset


def extract_QNM(t, signal, t_fit_start, t_fit_end,
                init_omega_r=0.37, init_omega_i=0.09):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    t_fit = t[mask]
    s_fit = signal[mask]

    A_guess = (np.max(s_fit) - np.min(s_fit)) / 2
    p0 = [A_guess, init_omega_r, init_omega_i, 0.0, 0.0]
    bounds = ([-np.inf, 0, 0, -np.pi, -np.inf],
              [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, pcov = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0,
                                bounds=bounds, maxfev=30000)
        A_fit, omega_r_fit, omega_i_fit, phi_fit, off_fit = popt
        omega = omega_r_fit - 1j * omega_i_fit
        return omega, popt
    except Exception as ex:
        print(f"  curve_fit failed: {ex}", flush=True)
        return None, None


# ---------------------------------------------------------------------------
# G91a:  Schwarzschild calibration
# ---------------------------------------------------------------------------

print("=" * 80)
print("G91a: Schwarzschild calibration in r*_GR coordinate")
print("=" * 80, flush=True)
print()

# Build r*_GR grid: span well beyond the potential peak.  Values much below
# -60 invert to r - 2M smaller than double precision can represent.
rs_min_GR = -60.0
rs_max_GR =  400.0
N_GR = 6001
rs_grid_GR = np.linspace(rs_min_GR, rs_max_GR, N_GR)
dr_star_GR = rs_grid_GR[1] - rs_grid_GR[0]
print(f"Grid: N = {N_GR}, r* ∈ [{rs_min_GR}, {rs_max_GR}], dr* = {dr_star_GR:.4e}", flush=True)

# Build r(r*) for GR
print("Inverting r*_GR(r) to build r(r*) ...", flush=True)
r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rs_grid_GR])
print(f"  r range: [{r_grid_GR.min():.6f}, {r_grid_GR.max():.6f}]", flush=True)

# Build V on grid
ell = 2
V_grid_GR = np.array([V_RW(r, ell, k_GR) for r in r_grid_GR])
print(f"  V max at r* = {rs_grid_GR[np.argmax(V_grid_GR)]:.4f}, "
      f"r = {r_grid_GR[np.argmax(V_grid_GR)]:.4f}", flush=True)

# Initial Gaussian pulse at r* = 30 (well outside potential peak at r*=0)
rs_pulse = 30.0
sigma_pulse = 3.0
psi0 = np.exp(-(rs_grid_GR - rs_pulse)**2 / (2 * sigma_pulse**2))

# Observation at r* = 50
i_obs_GR = int(np.argmin(np.abs(rs_grid_GR - 50.0)))
print(f"  Initial pulse at r* = {rs_pulse}, σ = {sigma_pulse}", flush=True)
print(f"  Observation at r* = {rs_grid_GR[i_obs_GR]:.4f}", flush=True)

dt = 0.5 * dr_star_GR  # CFL safety factor
T_final = 350.0
N_steps_est = int(T_final / dt)
print(f"  T_final = {T_final}, dt = {dt:.4e}, N_steps = {N_steps_est}", flush=True)
print()

print("Evolving Schwarzschild ...", flush=True)
t_GR, sig_GR, _, _ = evolve_TD(V_grid_GR, dr_star_GR, T_final, dt,
                                psi0, i_obs=i_obs_GR)
print(f"  Done.  max|ψ| = {np.max(np.abs(sig_GR)):.4e}", flush=True)
print()


# Fit late-time QNM
omega_gold = complex(0.373672, -0.088962)
t_fit_start = 100.0
t_fit_end = 250.0
print(f"Fitting QNM in t ∈ [{t_fit_start}, {t_fit_end}] M", flush=True)
omega_GR, popt_GR = extract_QNM(t_GR, sig_GR, t_fit_start, t_fit_end)
print(f"  ω_GR     = {omega_GR}")
print(f"  ω_gold   = {omega_gold}")
if omega_GR is not None:
    err_GR = abs(omega_GR - omega_gold) / abs(omega_gold)
    print(f"  |Δω/ω|   = {err_GR * 100:.4f} %", flush=True)
    if err_GR < 0.01:
        print("  ✓ PASS: calibration < 1%", flush=True)
        calib_ok = True
    elif err_GR < 0.05:
        print("  ⚠ ACCEPTABLE: 1-5%, can refine grid for sub-percent", flush=True)
        calib_ok = True
    else:
        print("  ✗ FAIL: > 5% calibration error. Refine before running STAM.", flush=True)
        calib_ok = False
print()


# ---------------------------------------------------------------------------
# G91b:  STAM run  (only if calibration passes)
# ---------------------------------------------------------------------------

if 'calib_ok' in dir() and calib_ok:
    print("=" * 80)
    print("G91b: STAM run in r*_STAM coordinate")
    print("=" * 80, flush=True)
    print()

    # Build r*_STAM grid.  STAM throat is power-law-deep:
    # at r* = -200, r ≈ 2M + 0.001M ≈ 2.001 M (still well inside shell)
    rs_min_STAM = -500.0   # deeper than Schwarzschild needed
    rs_max_STAM = 400.0
    N_STAM = 10001
    rs_grid_STAM = np.linspace(rs_min_STAM, rs_max_STAM, N_STAM)
    dr_star_STAM = rs_grid_STAM[1] - rs_grid_STAM[0]
    print(f"Grid: N = {N_STAM}, r*_STAM ∈ [{rs_min_STAM}, {rs_max_STAM}], "
          f"dr* = {dr_star_STAM:.4e}", flush=True)

    print("Building r(r*)_STAM via ODE integration ...", flush=True)
    r_grid_STAM = build_r_of_rstar_STAM(rs_grid_STAM)
    print(f"  r range: [{r_grid_STAM.min():.6f}, {r_grid_STAM.max():.6f}]", flush=True)
    eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
    print(f"  Innermost ε = (r-2M)/2M = {eps_min:.4e}", flush=True)

    V_grid_STAM = np.array([V_RW(r, ell, k_STAM) for r in r_grid_STAM])
    print(f"  V max at r* = {rs_grid_STAM[np.argmax(V_grid_STAM)]:.4f}, "
          f"r = {r_grid_STAM[np.argmax(V_grid_STAM)]:.4f}", flush=True)

    # Same pulse and observation r* position
    psi0_STAM = np.exp(-(rs_grid_STAM - rs_pulse)**2 / (2 * sigma_pulse**2))
    i_obs_STAM = int(np.argmin(np.abs(rs_grid_STAM - 50.0)))
    dt_STAM = 0.5 * dr_star_STAM

    print(f"  Observation at r* = {rs_grid_STAM[i_obs_STAM]:.4f}", flush=True)
    print(f"  dt = {dt_STAM:.4e}, N_steps = {int(T_final / dt_STAM)}", flush=True)
    print()

    print("Evolving STAM ...", flush=True)
    t_STAM, sig_STAM, _, _ = evolve_TD(V_grid_STAM, dr_star_STAM, T_final, dt_STAM,
                                        psi0_STAM, i_obs=i_obs_STAM)
    print(f"  Done.  max|ψ| = {np.max(np.abs(sig_STAM)):.4e}", flush=True)
    print()

    omega_STAM, popt_STAM = extract_QNM(t_STAM, sig_STAM, t_fit_start, t_fit_end)
    print(f"  ω_STAM   = {omega_STAM}", flush=True)
    if omega_STAM is not None and omega_GR is not None:
        dw_over_w = abs(omega_STAM - omega_GR) / abs(omega_GR)
        print(f"  |Δω/ω|   = {dw_over_w * 100:.4f} %  (vs Schwarzschild)", flush=True)
        re_shift = (omega_STAM.real - omega_GR.real) / omega_GR.real
        im_shift = (omega_STAM.imag - omega_GR.imag) / omega_GR.imag
        print(f"  Re shift: {re_shift * 100:+.4f} %", flush=True)
        print(f"  Im shift: {im_shift * 100:+.4f} %", flush=True)
    print()

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(rs_grid_GR, V_grid_GR, 'tab:blue', linewidth=2, label='V_GR(r*_GR)')
    ax.plot(rs_grid_STAM, V_grid_STAM, 'tab:orange', linewidth=2, label='V_STAM(r*_STAM)')
    ax.set_xlabel('r* / M')
    ax.set_ylabel('V_RW^proxy (ℓ=2)')
    ax.set_title('Effective potential vs r*')
    ax.set_xlim(-50, 100)
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(t_GR, np.abs(sig_GR), 'tab:blue', linewidth=1, label='Schwarzschild')
    ax.plot(t_STAM, np.abs(sig_STAM), 'tab:orange', linewidth=1, label='STAM')
    ax.axvspan(t_fit_start, t_fit_end, alpha=0.15, color='gray', label='fit window')
    ax.set_xlabel('t / M')
    ax.set_ylabel('|ψ(t, r_obs)|')
    ax.set_title('Time-domain signal (log)')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    mask_GR = (t_GR >= t_fit_start) & (t_GR <= t_fit_end)
    mask_STAM = (t_STAM >= t_fit_start) & (t_STAM <= t_fit_end)
    ax.plot(t_GR[mask_GR], sig_GR[mask_GR], 'tab:blue', linewidth=1, label='Schw')
    if popt_GR is not None:
        ax.plot(t_GR[mask_GR], damped_sinusoid(t_GR[mask_GR], *popt_GR),
                'k--', alpha=0.7, label='Schw fit')
    ax.plot(t_STAM[mask_STAM], sig_STAM[mask_STAM], 'tab:orange', linewidth=1, label='STAM')
    if popt_STAM is not None:
        ax.plot(t_STAM[mask_STAM], damped_sinusoid(t_STAM[mask_STAM], *popt_STAM),
                'r--', alpha=0.7, label='STAM fit')
    ax.set_xlabel('t / M')
    ax.set_ylabel('ψ(t, r_obs)')
    ax.set_title('Fit window with damped-sinusoid overlay')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.axis('off')
    txt = [
        f"G91 results — time-domain in r* coordinate",
        "",
        f"Calibration (Schwarzschild ℓ=2, n=0):",
        f"  ω_TD     = {omega_GR}",
        f"  ω_Leaver = {omega_gold}",
        f"  |Δω/ω|   = {err_GR*100:.4f} %",
        "",
        f"STAM ℓ=2, n=0:",
        f"  ω_TD     = {omega_STAM}",
        "",
        f"Framework shift vs Schwarzschild:",
        f"  |Δω/ω|   = {dw_over_w*100:.4f} %",
        "",
        f"Grid: dr*_GR = {dr_star_GR:.3e}, dr*_STAM = {dr_star_STAM:.3e}",
        f"Sommerfeld BCs, fit window t ∈ [{t_fit_start}, {t_fit_end}] M",
    ]
    ax.text(0.05, 0.95, '\n'.join(txt), transform=ax.transAxes,
            fontsize=10, family='monospace', verticalalignment='top')

    plt.tight_layout()
    out_png = PLOTS / "G91_time_domain_rstar.png"
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"Plot saved: {out_png}", flush=True)

    # Summary
    md = []
    md.append("# G91 — Time-domain QNM in r* coordinate (proper Sommerfeld)\n")
    md.append("**Date: 2026-05-13.**\n")
    md.append("## Schwarzschild calibration (ℓ = 2, n = 0)\n")
    md.append(f"- ω_TD     = {omega_GR.real:.6f} − {abs(omega_GR.imag):.6f} i\n")
    md.append(f"- ω_Leaver = {omega_gold.real:.6f} − {abs(omega_gold.imag):.6f} i\n")
    md.append(f"- **|Δω/ω| = {err_GR*100:.4f} %**\n")
    md.append("\n")
    md.append("## STAM result (ℓ = 2, n = 0)\n")
    md.append(f"- ω_STAM   = {omega_STAM.real:.6f} {'-' if omega_STAM.imag < 0 else '+'} "
              f"{abs(omega_STAM.imag):.6f} i\n")
    md.append(f"- **|Δω/ω| (vs Schwarzschild) = {dw_over_w*100:.4f} %**\n")
    md.append(f"- Re shift: {re_shift*100:+.4f} %\n")
    md.append(f"- Im shift: {im_shift*100:+.4f} %\n")
    md.append("\n")
    md.append("## Method\n")
    md.append("- Wave eq in r*: ∂²ψ/∂t² − ∂²ψ/∂r*² + V ψ = 0.\n")
    md.append(f"- GR grid:   N = {N_GR}, r* ∈ [{rs_min_GR}, {rs_max_GR}].\n")
    md.append(f"- STAM grid: N = {N_STAM}, r* ∈ [{rs_min_STAM}, {rs_max_STAM}] "
              f"(deeper to span the power-law STAM throat).\n")
    md.append("- Sommerfeld first-order outflow BCs at both edges.\n")
    md.append(f"- Initial Gaussian pulse at r* = {rs_pulse}, σ = {sigma_pulse}.\n")
    md.append(f"- Fit window: t ∈ [{t_fit_start}, {t_fit_end}] M.\n")
    md.append("\n")
    md.append("## Files\n")
    md.append("- [scripts/G91_time_domain_rstar.py](../scripts/G91_time_domain_rstar.py)\n")
    md.append("- [plots/G91_time_domain_rstar.png](../plots/G91_time_domain_rstar.png)\n")

    out_md = RESULTS / "G91_time_domain_rstar_summary.md"
    out_md.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}", flush=True)
else:
    print("Schwarzschild calibration failed; not running STAM.", flush=True)
    print("Refine grid resolution or boundary placement before retrying.", flush=True)
