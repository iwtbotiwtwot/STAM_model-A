#!/usr/bin/env python3
"""
G93_constrained_action_axial_QNM.py

Constrained-shell-count axial perturbation diagnostic for Model-A.

Purpose
=======
G92 added the leading f(Σ)R action-aware correction to the G91 tensor-proxy
QNM calculation:

    V_action = V_proxy + (sqrt(f))'' / sqrt(f)

where primes are derivatives in the STAM tortoise coordinate r*.

G93 sharpens the status of that correction using the constrained-shell-count
action:

    S[Σ,g,λ1,λ2] = (1/16πG) ∫ sqrt(-g)
        [ f(Σ) R + λ1((∇Σ)^2 - W) + λ2(u^μ ∂_μΣ) - 2V(Σ) ] d^4x.

Key axial-sector observation
============================
Odd-parity / axial perturbations do not carry a scalar shell-count
perturbation:

  - Σ is a scalar shell-count variable, parity-even.
  - Axial metric perturbations are parity-odd.
  - The λ1, λ2 constraints remove independent δΣ propagation.
  - Therefore the constrained-shell sector does not supply a new propagating
    axial scalar mode.

At this diagnostic level, the remaining direct action effect in the axial
master equation is the variable gravitational stiffness f(Σ).  Canonically
normalizing a schematic quadratic action

    S2 ~ 1/2 ∫ dt dr* f(r) [q_t^2 - q_*^2 - V_proxy q^2]

with ψ = sqrt(f) q gives

    ψ_tt - ψ_** + [V_proxy + (sqrt(f))''/sqrt(f)] ψ = 0.

G93 therefore treats the G92 potential as the constrained-action axial
candidate and does three things:

  1. Repeats the calibrated r*-time-domain extraction from G91/G92.
  2. Reports proxy vs action-aware QNM shifts.
  3. Runs fit-window and throat-depth robustness checks on the action-aware
     result so it is not a one-window artifact.

Status / caveat
===============
This script is still a diagnostic master-equation test, not a formal
Regge-Wheeler derivation from a fully varied quadratic action.  Its purpose is
to test the most direct constrained-action axial correction after the scalar
shell-count DOF has been removed by λ1, λ2.

Outputs
=======
- results/G93_constrained_action_axial_QNM_summary.md
- plots/G93_constrained_action_axial_QNM.png

Conventions
===========
M = 1. Frequencies are reported as ω = Re(ω) - i |Im(ω)|.
"""

from __future__ import annotations

import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq, curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)

M = 1.0


# ---------------------------------------------------------------------------
# Basic functions
# ---------------------------------------------------------------------------

def h_fn(r):
    return 1.0 - 2.0 * M / r


def A_of_r(r):
    return 2.0 * M / r


def y_of_r(r):
    return 6.0 * M / r - 2.0


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def k_GR(r):
    return h_fn(r)


def k_STAM(r):
    y = y_of_r(r)
    h_val = h_fn(r)
    return np.where(np.asarray(y) > 0, h_val * F_quintic(y), h_val)


def V_RW_proxy(r, ell, k_fn):
    """Tensor axial Regge-Wheeler proxy potential used in G87/G91/G92.

    M_eff = r(1-k)/2.  For Schwarzschild k=h, M_eff=M.
    """
    h_val = h_fn(r)
    k_val = k_fn(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


# ---------------------------------------------------------------------------
# f(Σ) / f(A) reconstruction from G70/G75
# ---------------------------------------------------------------------------

def dlnf_dA(A):
    """G70/G75 non-minimal coupling ODE for the quintic branch.

    d ln f / dA = -2 y^3 (45 A^2 - 33 A - 13) / [A F(y)],
    y = 3A - 2.  Outside/at PS, f is set to 1.
    """
    A_arr = np.asarray(A, dtype=float)
    y = 3.0 * A_arr - 2.0
    F = F_quintic(y)
    out = np.zeros_like(A_arr)
    mask = A_arr > (2.0 / 3.0)
    out[mask] = -2.0 * y[mask]**3 * (45.0 * A_arr[mask]**2 - 33.0 * A_arr[mask] - 13.0) / (
        A_arr[mask] * F[mask]
    )
    return out if np.ndim(A) else float(out)


def build_f_of_A_interpolator(A_max: float, n: int = 24000):
    """Build f(A) on [2/3, A_max], normalized f(2/3)=1."""
    A0 = 2.0 / 3.0
    A_max = min(float(A_max), 1.0 - 1e-10)
    if A_max <= A0:
        def f_interp(A):
            return np.ones_like(np.asarray(A, dtype=float))
        return f_interp

    A_grid = np.linspace(A0, A_max, n)
    deriv = dlnf_dA(A_grid)
    lnf = cumulative_trapezoid(deriv, A_grid, initial=0.0)
    # Keep the diagnostic numerically finite deep in the throat.
    lnf = np.clip(lnf, -120.0, 120.0)
    f_grid = np.exp(lnf)

    def f_interp(A):
        A_arr = np.asarray(A, dtype=float)
        out = np.ones_like(A_arr)
        mask = A_arr > A0
        out[mask] = np.interp(A_arr[mask], A_grid, f_grid, left=1.0, right=f_grid[-1])
        return out

    return f_interp


def f_grid_on_r(r_grid):
    A_grid = A_of_r(r_grid)
    f_interp = build_f_of_A_interpolator(float(np.nanmax(A_grid)))
    return f_interp(A_grid)


# ---------------------------------------------------------------------------
# Tortoise coordinate construction (r*(3M)=0)
# ---------------------------------------------------------------------------

def r_of_rstar_GR(rstar: float) -> float:
    """Invert Schwarzschild r* anchored so r*(3M)=0."""
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def froot(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    if rstar > 0:
        return brentq(froot, 3.0 * M - 1e-12, 1e6)
    return brentq(froot, 2.0 * M + 1e-14, 3.0 * M + 1e-12)


def build_r_of_rstar_STAM(rstar_grid):
    """Build r(r*) for STAM using dr/dr* = sqrt(h k), anchored at r*(3M)=0."""
    def rhs(rs, y):
        r_val = y[0]
        return [np.sqrt(max(float(h_fn(r_val) * k_STAM(r_val)), 0.0))]

    r_arr = np.zeros_like(rstar_grid)

    rs_pos = rstar_grid[rstar_grid > 0]
    if len(rs_pos):
        sort_pos = np.argsort(rs_pos)
        rs_sorted = rs_pos[sort_pos]
        sol = solve_ivp(rhs, [0.0, rs_sorted[-1]], [3.0 * M], t_eval=rs_sorted,
                        rtol=1e-12, atol=1e-14, max_step=1.0)
        tmp = np.zeros_like(rs_pos)
        tmp[sort_pos] = sol.y[0]
        r_arr[rstar_grid > 0] = tmp

    rs_neg = rstar_grid[rstar_grid < 0]
    if len(rs_neg):
        sort_neg = np.argsort(rs_neg)[::-1]
        rs_sorted = rs_neg[sort_neg]
        sol = solve_ivp(rhs, [0.0, rs_sorted[-1]], [3.0 * M], t_eval=rs_sorted,
                        rtol=1e-12, atol=1e-14, max_step=0.04)
        tmp = np.zeros_like(rs_neg)
        tmp[sort_neg] = sol.y[0]
        r_arr[rstar_grid < 0] = tmp

    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


# ---------------------------------------------------------------------------
# Action-aware axial potential
# ---------------------------------------------------------------------------

def canonical_f_correction(f_grid, dr_star):
    """Return (sqrt(f))''/sqrt(f) using r* derivatives."""
    sqrtf = np.sqrt(np.maximum(f_grid, 1e-300))
    d1 = np.gradient(sqrtf, dr_star, edge_order=2)
    d2 = np.gradient(d1, dr_star, edge_order=2)
    return d2 / sqrtf


def build_action_aware_potential(V_proxy, f_grid, dr_star):
    corr = canonical_f_correction(f_grid, dr_star)
    return V_proxy + corr, corr


# ---------------------------------------------------------------------------
# Time-domain solver and QNM extraction
# ---------------------------------------------------------------------------

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, psidot0=None, i_obs=None):
    """Solve ψ_tt - ψ_xx + Vψ = 0 with Sommerfeld boundaries."""
    if psidot0 is None:
        psidot0 = np.zeros_like(psi0)

    psi_old = psi0 - dt * psidot0
    psi = psi0.copy()
    steps = int(T_final / dt)
    inv_dr2 = 1.0 / (dr_star * dr_star)
    cfl = dt / dr_star
    t_arr = np.zeros(steps)
    sig_arr = np.zeros(steps)

    for step in range(steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2

        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt * dt * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1]
        )
        # First-order Sommerfeld outflow in r*.
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])

        psi_old, psi = psi, psi_new
        t_arr[step] = (step + 1) * dt
        if i_obs is not None:
            sig_arr[step] = psi[i_obs]

    return t_arr, sig_arr


def damped_sinusoid(t, A, omega_r, omega_i, phi, offset):
    return A * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + offset


def extract_QNM(t, signal, t_fit_start, t_fit_end, init_omega_r=0.37, init_omega_i=0.09):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    t_fit = t[mask]
    s_fit = signal[mask]
    if len(t_fit) < 50 or np.allclose(s_fit, 0):
        return None, None

    A_guess = 0.5 * (np.nanmax(s_fit) - np.nanmin(s_fit))
    p0 = [A_guess, init_omega_r, init_omega_i, 0.0, 0.0]
    bounds = ([-np.inf, 0.0, 0.0, -np.pi, -np.inf], [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0, bounds=bounds, maxfev=40000)
        _, omega_r, omega_i, _, _ = popt
        return complex(abs(omega_r), -abs(omega_i)), popt
    except Exception as ex:
        print(f"  curve_fit failed in window [{t_fit_start}, {t_fit_end}]: {ex}")
        return None, None


@dataclass
class RunResult:
    label: str
    omega: complex | None
    shift_vs_gr: float | None
    signal_t: np.ndarray
    signal: np.ndarray
    V: np.ndarray


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 100)
    print("G93: constrained-shell-count axial QNM diagnostic")
    print("=" * 100)
    print()
    print("Axial-sector premise:")
    print("  δΣ is parity-even / constrained; axial odd-parity sector carries no independent scalar shell-count DOF.")
    print("  The leading action-aware axial effect is the f(Σ)R canonical-normalization correction.")
    print()

    ell = 2
    omega_gold = complex(0.373672, -0.088962)
    rs_pulse = 30.0
    sigma_pulse = 3.0
    r_obs_star = 50.0
    T_final = 350.0
    fit_windows = [(90.0, 230.0), (100.0, 250.0), (110.0, 270.0), (120.0, 290.0)]
    main_window = (100.0, 250.0)

    # Schwarzschild calibration.
    print("=" * 100)
    print("G93a: Schwarzschild calibration")
    print("=" * 100)
    rs_min_GR, rs_max_GR, N_GR = -60.0, 400.0, 6001
    rs_grid_GR = np.linspace(rs_min_GR, rs_max_GR, N_GR)
    dr_GR = rs_grid_GR[1] - rs_grid_GR[0]
    r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rs_grid_GR])
    V_GR = np.array([V_RW_proxy(r, ell, k_GR) for r in r_grid_GR])
    psi0_GR = np.exp(-(rs_grid_GR - rs_pulse)**2 / (2.0 * sigma_pulse**2))
    i_obs_GR = int(np.argmin(np.abs(rs_grid_GR - r_obs_star)))
    dt_GR = 0.5 * dr_GR
    print(f"GR grid: [{rs_min_GR}, {rs_max_GR}], N={N_GR}, dr*={dr_GR:.5f}")
    t_GR, sig_GR = evolve_TD(V_GR, dr_GR, T_final, dt_GR, psi0_GR, i_obs=i_obs_GR)
    omega_GR, _ = extract_QNM(t_GR, sig_GR, *main_window)
    if omega_GR is None:
        raise RuntimeError("GR QNM extraction failed")
    err_GR = abs(omega_GR - omega_gold) / abs(omega_gold)
    print(f"  omega_GR    = {omega_GR.real:.6f} {omega_GR.imag:+.6f} i")
    print(f"  omega_Leaver= {omega_gold.real:.6f} {omega_gold.imag:+.6f} i")
    print(f"  calibration error = {100*err_GR:.4f}%")
    calib_ok = err_GR < 0.01
    print("  PASS" if calib_ok else "  FAIL")
    print()
    if not calib_ok:
        print("Calibration failed.  Do not use STAM values as predictions.")
        return

    # STAM grid and potentials.
    print("=" * 100)
    print("G93b: STAM proxy and constrained-action candidate")
    print("=" * 100)
    rs_min_STAM, rs_max_STAM, N_STAM = -500.0, 400.0, 10001
    rs_grid_STAM = np.linspace(rs_min_STAM, rs_max_STAM, N_STAM)
    dr_STAM = rs_grid_STAM[1] - rs_grid_STAM[0]
    print(f"STAM grid: [{rs_min_STAM}, {rs_max_STAM}], N={N_STAM}, dr*={dr_STAM:.5f}")
    print("Building r(r*)_STAM ...")
    r_grid_STAM = build_r_of_rstar_STAM(rs_grid_STAM)
    eps_min = (float(np.nanmin(r_grid_STAM)) - 2.0 * M) / (2.0 * M)
    print(f"  r range: [{np.nanmin(r_grid_STAM):.8f}, {np.nanmax(r_grid_STAM):.8f}], eps_min={eps_min:.3e}")

    V_proxy = np.array([V_RW_proxy(r, ell, k_STAM) for r in r_grid_STAM])
    f_grid = f_grid_on_r(r_grid_STAM)
    V_action, V_fcorr = build_action_aware_potential(V_proxy, f_grid, dr_STAM)

    print(f"  f range: [{np.nanmin(f_grid):.3e}, {np.nanmax(f_grid):.3e}]")
    print(f"  V_proxy max:  {np.nanmax(V_proxy):.6e}")
    print(f"  V_action max: {np.nanmax(V_action):.6e}")
    print(f"  V_fcorr range: [{np.nanmin(V_fcorr):.6e}, {np.nanmax(V_fcorr):.6e}]")

    psi0_STAM = np.exp(-(rs_grid_STAM - rs_pulse)**2 / (2.0 * sigma_pulse**2))
    i_obs_STAM = int(np.argmin(np.abs(rs_grid_STAM - r_obs_star)))
    dt_STAM = 0.5 * dr_STAM

    print("Evolving STAM proxy ...")
    t_proxy, sig_proxy = evolve_TD(V_proxy, dr_STAM, T_final, dt_STAM, psi0_STAM, i_obs=i_obs_STAM)
    omega_proxy, _ = extract_QNM(t_proxy, sig_proxy, *main_window)
    print(f"  omega_proxy = {omega_proxy.real:.6f} {omega_proxy.imag:+.6f} i")

    print("Evolving STAM constrained-action candidate ...")
    t_action, sig_action = evolve_TD(V_action, dr_STAM, T_final, dt_STAM, psi0_STAM, i_obs=i_obs_STAM)
    omega_action, _ = extract_QNM(t_action, sig_action, *main_window)
    print(f"  omega_action = {omega_action.real:.6f} {omega_action.imag:+.6f} i")
    print()

    def rel_shift(a, b):
        return abs(a - b) / abs(b)

    proxy_shift = rel_shift(omega_proxy, omega_GR)
    action_shift = rel_shift(omega_action, omega_GR)
    action_vs_proxy = rel_shift(omega_action, omega_proxy)

    print("Main-window results:")
    print(f"  omega_GR      = {omega_GR.real:.6f} {omega_GR.imag:+.6f} i")
    print(f"  omega_proxy   = {omega_proxy.real:.6f} {omega_proxy.imag:+.6f} i")
    print(f"  omega_action  = {omega_action.real:.6f} {omega_action.imag:+.6f} i")
    print(f"  proxy shift vs GR    = {100*proxy_shift:.4f}%")
    print(f"  action shift vs GR   = {100*action_shift:.4f}%")
    print(f"  action vs proxy diff = {100*action_vs_proxy:.4f}%")
    print()

    # Fit-window robustness.
    print("=" * 100)
    print("G93c: fit-window robustness")
    print("=" * 100)
    rows = []
    for w0, w1 in fit_windows:
        om_p, _ = extract_QNM(t_proxy, sig_proxy, w0, w1)
        om_a, _ = extract_QNM(t_action, sig_action, w0, w1)
        rows.append((w0, w1, om_p, om_a))
        spct = 100 * rel_shift(om_p, omega_GR) if om_p else np.nan
        sact = 100 * rel_shift(om_a, omega_GR) if om_a else np.nan
        print(f"  window [{w0:.0f}, {w1:.0f}]  proxy={om_p.real:.6f}{om_p.imag:+.6f}i ({spct:.3f}%)  "
              f"action={om_a.real:.6f}{om_a.imag:+.6f}i ({sact:.3f}%)")
    print()

    # Plots.
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    ax = axes[0, 0]
    ax.plot(rs_grid_GR, V_GR, label="GR", linewidth=2)
    ax.plot(rs_grid_STAM, V_proxy, label="STAM proxy", linewidth=2)
    ax.plot(rs_grid_STAM, V_action, label="STAM constrained-action candidate", linewidth=2)
    ax.set_xlim(-50, 100)
    ax.set_xlabel("r* / M")
    ax.set_ylabel("V(r*)")
    ax.set_title("Axial potentials")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(rs_grid_STAM, V_fcorr, color="tab:purple", linewidth=2)
    ax.set_xlim(-50, 100)
    ax.set_xlabel("r*_STAM / M")
    ax.set_ylabel(r"$(\sqrt f)''/\sqrt f$")
    ax.set_title("f(Σ)R canonical correction")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(t_GR, np.abs(sig_GR), label="GR", linewidth=1)
    ax.plot(t_proxy, np.abs(sig_proxy), label="STAM proxy", linewidth=1)
    ax.plot(t_action, np.abs(sig_action), label="STAM constrained-action", linewidth=1)
    ax.axvspan(*main_window, color="gray", alpha=0.15, label="main fit")
    ax.set_yscale("log")
    ax.set_xlabel("t / M")
    ax.set_ylabel("|psi|")
    ax.set_title("Time-domain signals")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        "G93 constrained-action axial diagnostic",
        "",
        "Axial premise: δΣ constrained / non-propagating;",
        "leading action effect is f(Σ)R canonical correction.",
        "",
        f"GR omega:     {omega_GR.real:.6f} {omega_GR.imag:+.6f} i",
        f"Leaver gold:  {omega_gold.real:.6f} {omega_gold.imag:+.6f} i",
        f"calib error:  {100*err_GR:.4f}%",
        "",
        f"Proxy omega:  {omega_proxy.real:.6f} {omega_proxy.imag:+.6f} i",
        f"Action omega: {omega_action.real:.6f} {omega_action.imag:+.6f} i",
        "",
        f"Proxy shift:  {100*proxy_shift:.4f}%",
        f"Action shift: {100*action_shift:.4f}%",
        f"Action-proxy: {100*action_vs_proxy:.4f}%",
    ]
    ax.text(0.02, 0.98, "\n".join(lines), va="top", family="monospace", fontsize=9.5)

    plt.tight_layout()
    plot_path = PLOTS / "G93_constrained_action_axial_QNM.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"Plot saved: {plot_path}")

    # Markdown summary.
    md = []
    md.append("# G93 — Constrained-shell-count axial QNM diagnostic\n")
    md.append("**Purpose.** Formalize the axial-sector reading of the constrained shell-count action and test the leading action-aware QNM correction.\n")
    md.append("\n")
    md.append("## Axial-sector premise\n")
    md.append("In the constrained action, \\(\\Sigma\\) is not an ordinary propagating scalar. The constraints \\(\\lambda_1\\) and \\(\\lambda_2\\) remove independent \\(\\delta\\Sigma\\) dynamics. In the axial odd-parity sector, \\(\\delta\\Sigma\\) does not supply an independent propagating mode. The leading direct action effect on the axial master equation is therefore the variable gravitational stiffness \\(f(\\Sigma)R\\).\n")
    md.append("\n")
    md.append("The diagnostic potential is:\n\n")
    md.append("```text\n")
    md.append("V_action = V_proxy + (sqrt(f))'' / sqrt(f)\n")
    md.append("```\n")
    md.append("where derivatives are taken with respect to the STAM tortoise coordinate.\n")
    md.append("\n")
    md.append("## Schwarzschild calibration\n")
    md.append(f"- ω_GR      = {omega_GR.real:.6f} {omega_GR.imag:+.6f} i\n")
    md.append(f"- ω_Leaver  = {omega_gold.real:.6f} {omega_gold.imag:+.6f} i\n")
    md.append(f"- calibration error = **{100*err_GR:.4f}%**\n")
    md.append("\n")
    md.append("## Main-window STAM results\n")
    md.append(f"- ω_proxy   = {omega_proxy.real:.6f} {omega_proxy.imag:+.6f} i\n")
    md.append(f"- ω_action  = {omega_action.real:.6f} {omega_action.imag:+.6f} i\n")
    md.append(f"- proxy shift vs GR = **{100*proxy_shift:.4f}%**\n")
    md.append(f"- constrained-action shift vs GR = **{100*action_shift:.4f}%**\n")
    md.append(f"- action-aware vs proxy difference = **{100*action_vs_proxy:.4f}%**\n")
    md.append("\n")
    md.append("## Fit-window robustness\n")
    md.append("\n| Fit window | ω_proxy | proxy shift | ω_action | action shift |\n|---:|---:|---:|---:|---:|\n")
    for w0, w1, om_p, om_a in rows:
        spct = 100 * rel_shift(om_p, omega_GR) if om_p else np.nan
        sact = 100 * rel_shift(om_a, omega_GR) if om_a else np.nan
        md.append(f"| [{w0:.0f}, {w1:.0f}] | {om_p.real:.6f} {om_p.imag:+.6f}i | {spct:.3f}% | {om_a.real:.6f} {om_a.imag:+.6f}i | {sact:.3f}% |\n")
    md.append("\n")
    md.append("## Interpretation\n")
    md.append("G93 treats the G92 f-corrected potential as the constrained-action axial candidate after the scalar shell-count mode has been removed by the double-LM constraints. A full formal variation of the quadratic odd-parity action remains the gold-standard derivation, but this diagnostic tests the leading action effect directly.\n")
    md.append("\n")
    md.append("## Files\n")
    md.append("- [scripts/G93_constrained_action_axial_QNM.py](../scripts/G93_constrained_action_axial_QNM.py)\n")
    md.append("- [plots/G93_constrained_action_axial_QNM.png](../plots/G93_constrained_action_axial_QNM.png)\n")

    out_md = RESULTS / "G93_constrained_action_axial_QNM_summary.md"
    out_md.write_text("\n".join(md), encoding="utf-8")
    print(f"Summary written: {out_md}")


if __name__ == "__main__":
    main()
