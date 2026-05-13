#!/usr/bin/env python3
"""
G94_QNM_robustness_sweep.py

Robustness sweep for the calibrated low-ell QNM estimate from G91-G93.

Purpose
=======
G91 produced the first calibrated low-ell time-domain QNM estimate using the
proper STAM tortoise coordinate.  G92/G93 added the leading action-aware axial
correction from the constrained shell-count action:

    V_action = V_proxy + (sqrt(f))'' / sqrt(f)

where derivatives are taken with respect to r*_STAM.

G94 asks whether the action-aware percent-level shift is stable under changes
in numerical/extraction choices:

  1. fit-window choice
  2. STAM throat depth / inner r* boundary
  3. outer r* boundary
  4. grid resolution
  5. pulse center and width
  6. observer location
  7. f-normalization

Important status
================
This is still a diagnostic axial-potential calculation, not the final formal
odd-parity perturbation theorem from the fully varied constrained-shell-count
action.  It is meant to determine whether the G93 action-aware QNM shift is a
stable numerical result or a one-window artifact.

Outputs
=======
- results/G94_QNM_robustness_sweep_summary.md
- results/G94_QNM_robustness_sweep.csv
- plots/G94_QNM_robustness_sweep.png

Conventions
===========
M = 1. Frequencies are reported as omega = Re(omega) - i |Im(omega)|.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

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
ELL = 2
OMEGA_LEAVER = complex(0.373672, -0.088962)


# ---------------------------------------------------------------------------
# Metric / potential functions
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
    """Tensor axial Regge-Wheeler proxy potential used in G87/G91/G92/G93.

    M_eff = r(1-k)/2.  For Schwarzschild k=h, M_eff=M.
    """
    h_val = h_fn(r)
    k_val = k_fn(r)
    M_eff = 0.5 * r * (1.0 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6.0 * M_eff / r**3)


# ---------------------------------------------------------------------------
# f(Sigma) / f(A) diagnostic reconstruction
# ---------------------------------------------------------------------------

def dlnf_dA(A):
    """G70/G75 non-minimal coupling ODE for the quintic branch.

    d ln f / dA = -2 y^3 (45 A^2 - 33 A - 13) / [A F(y)].
    Outside/at PS, f is set to 1.
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


def build_f_of_A_interpolator(A_max: float, n: int = 24000, f_scale: float = 1.0):
    """Build f(A) on [2/3, A_max], normalized f(2/3)=f_scale.

    The canonical correction (sqrt(f))''/sqrt(f) should be invariant under
    constant f_scale.  G94 explicitly checks this.
    """
    A0 = 2.0 / 3.0
    A_max = min(float(A_max), 1.0 - 1e-10)
    if A_max <= A0:
        def f_interp(A):
            return f_scale * np.ones_like(np.asarray(A, dtype=float))
        return f_interp

    A_grid = np.linspace(A0, A_max, n)
    deriv = dlnf_dA(A_grid)
    lnf = cumulative_trapezoid(deriv, A_grid, initial=0.0)
    lnf = np.clip(lnf, -120.0, 120.0)
    f_grid = f_scale * np.exp(lnf)

    def f_interp(A):
        A_arr = np.asarray(A, dtype=float)
        out = f_scale * np.ones_like(A_arr)
        mask = A_arr > A0
        out[mask] = np.interp(A_arr[mask], A_grid, f_grid, left=f_scale, right=f_grid[-1])
        return out

    return f_interp


def f_grid_on_r(r_grid, f_scale: float = 1.0):
    A_grid = A_of_r(r_grid)
    f_interp = build_f_of_A_interpolator(float(np.nanmax(A_grid)), f_scale=f_scale)
    return f_interp(A_grid)


# ---------------------------------------------------------------------------
# Tortoise coordinate construction, anchored r*(3M)=0
# ---------------------------------------------------------------------------

def r_of_rstar_GR(rstar: float) -> float:
    """Invert Schwarzschild r*, anchored so r*(3M)=0."""
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def froot(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    if rstar > 0:
        return brentq(froot, 3.0 * M - 1e-12, 1e6)
    return brentq(froot, 2.0 * M + 1e-14, 3.0 * M + 1e-12)


def build_r_of_rstar_STAM(rstar_grid, max_step_in: float = 0.04):
    """Build r(r*) for STAM using dr/dr* = sqrt(h k)."""
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
                        rtol=1e-12, atol=1e-14, max_step=max_step_in)
        tmp = np.zeros_like(rs_neg)
        tmp[sort_neg] = sol.y[0]
        r_arr[rstar_grid < 0] = tmp

    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


# ---------------------------------------------------------------------------
# Action-aware potential
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
# Time-domain solver and extraction
# ---------------------------------------------------------------------------

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, psidot0=None, i_obs=None):
    """Solve psi_tt - psi_xx + V psi = 0 with Sommerfeld boundaries."""
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
        return None

    A_guess = 0.5 * (np.nanmax(s_fit) - np.nanmin(s_fit))
    p0 = [A_guess, init_omega_r, init_omega_i, 0.0, 0.0]
    bounds = ([-np.inf, 0.0, 0.0, -np.pi, -np.inf], [np.inf, 2.0, 2.0, np.pi, np.inf])
    try:
        popt, _ = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0, bounds=bounds, maxfev=50000)
        _, omega_r, omega_i, _, _ = popt
        return complex(abs(omega_r), -abs(omega_i))
    except Exception as ex:
        print(f"    curve_fit failed in window [{t_fit_start}, {t_fit_end}]: {ex}")
        return None


def rel_shift(a: complex, b: complex) -> float:
    return abs(a - b) / abs(b)


# ---------------------------------------------------------------------------
# Config / variants
# ---------------------------------------------------------------------------

@dataclass
class Config:
    name: str
    rs_min_stam: float = -500.0
    rs_max_stam: float = 400.0
    n_stam: int = 10001
    rs_min_gr: float = -60.0
    rs_max_gr: float = 400.0
    n_gr: int = 6001
    rs_pulse: float = 30.0
    sigma_pulse: float = 3.0
    r_obs_star: float = 50.0
    T_final: float = 350.0
    fit_start: float = 100.0
    fit_end: float = 250.0
    f_scale: float = 1.0


@dataclass
class Row:
    name: str
    group: str
    omega_gr_re: float
    omega_gr_im: float
    gr_calib_pct: float
    omega_proxy_re: float
    omega_proxy_im: float
    proxy_shift_pct: float
    omega_action_re: float
    omega_action_im: float
    action_shift_pct: float
    action_vs_proxy_pct: float
    status: str


def omega_to_str(z: complex | None) -> str:
    if z is None:
        return "None"
    return f"{z.real:.6f} {z.imag:+.6f}i"


# ---------------------------------------------------------------------------
# Simulation wrappers
# ---------------------------------------------------------------------------

def run_GR(cfg: Config):
    rs_grid = np.linspace(cfg.rs_min_gr, cfg.rs_max_gr, cfg.n_gr)
    dr = rs_grid[1] - rs_grid[0]
    r_grid = np.array([r_of_rstar_GR(rs) for rs in rs_grid])
    V = np.array([V_RW_proxy(r, ELL, k_GR) for r in r_grid])
    psi0 = np.exp(-(rs_grid - cfg.rs_pulse)**2 / (2.0 * cfg.sigma_pulse**2))
    i_obs = int(np.argmin(np.abs(rs_grid - cfg.r_obs_star)))
    dt = 0.5 * dr
    t, sig = evolve_TD(V, dr, cfg.T_final, dt, psi0, i_obs=i_obs)
    omega = extract_QNM(t, sig, cfg.fit_start, cfg.fit_end)
    return omega


def run_STAM(cfg: Config):
    rs_grid = np.linspace(cfg.rs_min_stam, cfg.rs_max_stam, cfg.n_stam)
    dr = rs_grid[1] - rs_grid[0]
    r_grid = build_r_of_rstar_STAM(rs_grid)
    V_proxy = np.array([V_RW_proxy(r, ELL, k_STAM) for r in r_grid])
    f_grid = f_grid_on_r(r_grid, f_scale=cfg.f_scale)
    V_action, _ = build_action_aware_potential(V_proxy, f_grid, dr)

    psi0 = np.exp(-(rs_grid - cfg.rs_pulse)**2 / (2.0 * cfg.sigma_pulse**2))
    i_obs = int(np.argmin(np.abs(rs_grid - cfg.r_obs_star)))
    dt = 0.5 * dr

    t_p, sig_p = evolve_TD(V_proxy, dr, cfg.T_final, dt, psi0, i_obs=i_obs)
    omega_proxy = extract_QNM(t_p, sig_p, cfg.fit_start, cfg.fit_end)

    t_a, sig_a = evolve_TD(V_action, dr, cfg.T_final, dt, psi0, i_obs=i_obs)
    omega_action = extract_QNM(t_a, sig_a, cfg.fit_start, cfg.fit_end)
    return omega_proxy, omega_action


def run_variant(cfg: Config, group: str) -> Row:
    print(f"\n--- Running {group}: {cfg.name} ---", flush=True)
    try:
        omega_gr = run_GR(cfg)
        if omega_gr is None:
            raise RuntimeError("GR extraction failed")
        calib = 100.0 * rel_shift(omega_gr, OMEGA_LEAVER)
        print(f"  GR     {omega_to_str(omega_gr)}  calib={calib:.3f}%", flush=True)

        omega_proxy, omega_action = run_STAM(cfg)
        if omega_proxy is None or omega_action is None:
            raise RuntimeError("STAM extraction failed")

        proxy_shift = 100.0 * rel_shift(omega_proxy, omega_gr)
        action_shift = 100.0 * rel_shift(omega_action, omega_gr)
        action_vs_proxy = 100.0 * rel_shift(omega_action, omega_proxy)

        print(f"  proxy  {omega_to_str(omega_proxy)}  shift={proxy_shift:.3f}%", flush=True)
        print(f"  action {omega_to_str(omega_action)}  shift={action_shift:.3f}%", flush=True)
        status = "PASS" if calib < 1.0 else "CALIB_FAIL"
        return Row(cfg.name, group,
                   omega_gr.real, omega_gr.imag, calib,
                   omega_proxy.real, omega_proxy.imag, proxy_shift,
                   omega_action.real, omega_action.imag, action_shift,
                   action_vs_proxy, status)
    except Exception as ex:
        print(f"  FAILED: {ex}", flush=True)
        return Row(cfg.name, group, np.nan, np.nan, np.nan,
                   np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                   np.nan, f"FAILED: {ex}")


def make_variants(full: bool) -> list[tuple[str, Config]]:
    base = Config(name="baseline")
    variants: list[tuple[str, Config]] = [("baseline", base)]

    # Fit-window sweep: same evolution parameters; changed extraction window.
    for w0, w1 in [(90, 230), (100, 250), (110, 270), (120, 290)]:
        variants.append(("fit_window", Config(name=f"fit_{w0}_{w1}", fit_start=float(w0), fit_end=float(w1))))

    # Throat depth / inner r* boundary. Keep dr* close to baseline by adjusting N.
    for rs_min in ([-300.0, -500.0, -700.0] if not full else [-250.0, -400.0, -500.0, -700.0, -900.0]):
        n = int(round((400.0 - rs_min) / 0.09)) + 1
        variants.append(("throat_depth", Config(name=f"rsmin_{rs_min:g}", rs_min_stam=rs_min, n_stam=n)))

    # Outer boundary.
    for rs_max in ([300.0, 400.0, 500.0] if not full else [250.0, 300.0, 400.0, 500.0, 650.0]):
        n = int(round((rs_max - (-500.0)) / 0.09)) + 1
        variants.append(("outer_boundary", Config(name=f"rsmax_{rs_max:g}", rs_max_stam=rs_max, n_stam=n,
                                                   rs_max_gr=rs_max, n_gr=int(round((rs_max - (-60.0)) / 0.0767)) + 1)))

    # Grid resolution.
    for n_stam in ([8001, 10001, 12001] if not full else [7001, 9001, 10001, 12001, 15001]):
        variants.append(("grid_resolution", Config(name=f"Nstam_{n_stam}", n_stam=n_stam)))

    # Pulse center/width.
    pulse_cases = [(20.0, 3.0), (30.0, 3.0), (40.0, 3.0), (30.0, 2.0), (30.0, 4.0)]
    if full:
        pulse_cases += [(20.0, 2.0), (20.0, 4.0), (40.0, 2.0), (40.0, 4.0)]
    for center, width in pulse_cases:
        variants.append(("pulse", Config(name=f"pulse_{center:g}_sig_{width:g}", rs_pulse=center, sigma_pulse=width)))

    # Observer location.
    for obs in ([40.0, 50.0, 70.0] if not full else [30.0, 40.0, 50.0, 70.0, 90.0]):
        variants.append(("observer", Config(name=f"obs_{obs:g}", r_obs_star=obs)))

    # f-normalization: should leave V_action invariant.
    for scale in [0.1, 1.0, 10.0]:
        variants.append(("f_normalization", Config(name=f"fscale_{scale:g}", f_scale=scale)))

    # Remove exact duplicate baseline-like variants by name/group while preserving order.
    seen = set()
    out = []
    for group, cfg in variants:
        key = (group, cfg.name)
        if key not in seen:
            seen.add(key)
            out.append((group, cfg))
    return out


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(rows: list[Row], path: Path):
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_summary(rows: list[Row], path: Path, full: bool):
    pass_rows = [r for r in rows if r.status == "PASS"]
    action_shifts = np.array([r.action_shift_pct for r in pass_rows if np.isfinite(r.action_shift_pct)])
    proxy_shifts = np.array([r.proxy_shift_pct for r in pass_rows if np.isfinite(r.proxy_shift_pct)])
    calib = np.array([r.gr_calib_pct for r in pass_rows if np.isfinite(r.gr_calib_pct)])

    md = []
    md.append("# G94 — QNM robustness sweep\n")
    md.append("**Purpose.** Test whether the G93 action-aware low-ell QNM shift is robust under fit-window, grid, throat-depth, outer-boundary, pulse, observer, and f-normalization choices.\n")
    md.append("\n")
    md.append(f"Mode: {'full' if full else 'quick'} sweep.\n")
    md.append("\n")
    md.append("## Summary statistics for calibrated PASS rows\n")
    if len(pass_rows):
        md.append(f"- PASS rows: {len(pass_rows)} / {len(rows)}\n")
        md.append(f"- Schwarzschild calibration error: mean {np.nanmean(calib):.3f}%, range {np.nanmin(calib):.3f}%–{np.nanmax(calib):.3f}%\n")
        md.append(f"- Proxy shift: mean {np.nanmean(proxy_shifts):.3f}%, range {np.nanmin(proxy_shifts):.3f}%–{np.nanmax(proxy_shifts):.3f}%\n")
        md.append(f"- Action-aware shift: mean {np.nanmean(action_shifts):.3f}%, range {np.nanmin(action_shifts):.3f}%–{np.nanmax(action_shifts):.3f}%\n")
        md.append(f"- Action-aware shift std: {np.nanstd(action_shifts):.3f}%\n")
    else:
        md.append("- No calibrated PASS rows.\n")
    md.append("\n")
    md.append("## Results table\n\n")
    md.append("| Group | Variant | Calib % | ω_GR | Proxy shift % | Action shift % | Status |\n")
    md.append("|---|---|---:|---:|---:|---:|---|\n")
    for r in rows:
        omgr = f"{r.omega_gr_re:.6f}{r.omega_gr_im:+.6f}i" if np.isfinite(r.omega_gr_re) else "—"
        calib_s = f"{r.gr_calib_pct:.3f}" if np.isfinite(r.gr_calib_pct) else "—"
        proxy_s = f"{r.proxy_shift_pct:.3f}" if np.isfinite(r.proxy_shift_pct) else "—"
        action_s = f"{r.action_shift_pct:.3f}" if np.isfinite(r.action_shift_pct) else "—"
        md.append(f"| {r.group} | {r.name} | {calib_s} | {omgr} | {proxy_s} | {action_s} | {r.status} |\n")
    md.append("\n")
    md.append("## Interpretation guide\n")
    md.append("A robust action-aware QNM estimate should keep Schwarzschild calibration below 1% while leaving the action-aware shift in a narrow band under numerical/extraction variations. f-normalization should not change the action-aware result because (sqrt(f))''/sqrt(f) is invariant under constant rescaling of f.\n")
    md.append("\n")
    md.append("## Files\n")
    md.append("- `scripts/G94_QNM_robustness_sweep.py`\n")
    md.append("- `results/G94_QNM_robustness_sweep.csv`\n")
    md.append("- `plots/G94_QNM_robustness_sweep.png`\n")
    path.write_text("".join(md), encoding="utf-8")


def make_plot(rows: list[Row], path: Path):
    pass_rows = [r for r in rows if r.status == "PASS" and np.isfinite(r.action_shift_pct)]
    if not pass_rows:
        return
    groups = sorted(set(r.group for r in pass_rows))
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=False)

    # Top: action shifts by variant.
    labels = [f"{r.group}\n{r.name}" for r in pass_rows]
    x = np.arange(len(pass_rows))
    action = [r.action_shift_pct for r in pass_rows]
    proxy = [r.proxy_shift_pct for r in pass_rows]
    calib = [r.gr_calib_pct for r in pass_rows]

    ax = axes[0]
    ax.plot(x, action, marker='o', label='action-aware shift vs GR')
    ax.plot(x, proxy, marker='s', label='proxy shift vs GR', alpha=0.75)
    ax.axhline(np.nanmean(action), color='tab:blue', linestyle='--', alpha=0.5, label='mean action shift')
    ax.set_ylabel('QNM shift (%)')
    ax.set_title('G94 robustness sweep — proxy vs action-aware shifts')
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1]
    ax.plot(x, calib, marker='o', color='tab:red', label='Schwarzschild calibration error')
    ax.axhline(1.0, color='black', linestyle='--', label='1% gate')
    ax.set_ylabel('Calibration error (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=70, ha='right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="G94 QNM robustness sweep")
    parser.add_argument("--full", action="store_true", help="Run larger sweep. Default is quicker but still multi-variant.")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of variants for testing.")
    args = parser.parse_args()

    print("=" * 100)
    print("G94: QNM robustness sweep")
    print("=" * 100)
    print()
    print("Baseline target from G93: proxy shift ~1.2%, action-aware shift ~5.35%.")
    print("G94 tests whether that action-aware shift survives numerical/extraction changes.")

    variants = make_variants(full=args.full)
    if args.limit is not None:
        variants = variants[:args.limit]
    print(f"Variants to run: {len(variants)}")

    rows: list[Row] = []
    for group, cfg in variants:
        rows.append(run_variant(cfg, group))

    csv_path = RESULTS / "G94_QNM_robustness_sweep.csv"
    md_path = RESULTS / "G94_QNM_robustness_sweep_summary.md"
    plot_path = PLOTS / "G94_QNM_robustness_sweep.png"
    write_csv(rows, csv_path)
    write_summary(rows, md_path, full=args.full)
    make_plot(rows, plot_path)

    print("\n" + "=" * 100)
    print("G94 summary")
    print("=" * 100)
    pass_rows = [r for r in rows if r.status == "PASS" and np.isfinite(r.action_shift_pct)]
    if pass_rows:
        shifts = np.array([r.action_shift_pct for r in pass_rows])
        print(f"PASS rows: {len(pass_rows)} / {len(rows)}")
        print(f"Action-aware shift range: {np.nanmin(shifts):.3f}%–{np.nanmax(shifts):.3f}%")
        print(f"Action-aware shift mean ± std: {np.nanmean(shifts):.3f}% ± {np.nanstd(shifts):.3f}%")
    else:
        print("No PASS rows.")
    print(f"CSV:     {csv_path}")
    print(f"Summary: {md_path}")
    print(f"Plot:    {plot_path}")


if __name__ == "__main__":
    main()
