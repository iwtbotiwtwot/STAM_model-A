#!/usr/bin/env python3
"""
G88_leaver_proxy_QNM_l2.py

Phase 2 of G87: time-domain wave evolution to extract the fundamental
ℓ = 2, n = 0 QNM frequency on the framework's tensor RW proxy potential,
benchmarked against Schwarzschild.

Method (Option B from G87's two-phase plan):
============================================
Evolve the wave equation in r-coordinate:

    ∂²ψ/∂t² = (h k) ∂²ψ/∂r² + (1/2)(d(h k)/dr) ∂ψ/∂r - V_RW^proxy(r) ψ

(This is the r-coordinate form of  ∂²ψ/∂t² - ∂²ψ/∂r*² + V ψ = 0
with r* the tortoise coordinate defined by dr/dr* = sqrt(h k).)

Initial condition: Gaussian pulse centered at moderate r, at rest.
Boundary condition: absorbing (waves leave the domain).
Observation point: a fixed r between the photon sphere and the outer edge.
After the initial transient, the signal ringdowns at the dominant QNM
frequency.

Frequency extraction: fit late-time signal ψ(t) at the observation point
to A · exp(-ω_i t) · cos(ω_r t + φ).  Cross-check via Prony / matrix-pencil
if curve_fit struggles.

Validation:
- For Schwarzschild ℓ = 2, n = 0, the gold-standard Leaver value is
    ω_22^Schw = 0.373672 - 0.088962 i  (M = 1)
  The time-domain run should reproduce this to ~3-5% accuracy with
  modest grid (sufficient to distinguish framework vs GR if Δω/ω > 5%).

- For framework (STAM), the same method gives ω_22^STAM.

The Δω/ω at ℓ = 2, n = 0 from this run is the deliverable.  If > 10%,
framework has a near-term LIGO ringdown prediction.  If small (< 1%),
G86 / G87 Phase 1 were WKB noise at low ℓ.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Metric and potential
# ---------------------------------------------------------------------------

def h_fn(r, M=1.0):
    return 1 - 2 * M / r


def F_quintic(y):
    return 1 - 5 * y**4 + 4 * y**5


def k_Schw(r, M=1.0):
    return h_fn(r, M)


def k_STAM(r, M=1.0):
    A = 2 * M / r
    y = 3 * A - 2
    F = F_quintic(y)
    h_val = h_fn(r, M)
    # Inside PS (y > 0): k = h * F.  Outside PS (y <= 0): k = h.
    return np.where(y > 0, h_val * F, h_val)


def V_RW_proxy(r, ell, k_fn, M=1.0):
    """Effective tensor axial Regge-Wheeler proxy potential."""
    h_val = h_fn(r, M)
    k_val = k_fn(r, M)
    M_eff = (r / 2) * (1 - k_val)
    return h_val * (ell * (ell + 1) / r**2 - 6 * M_eff / r**3)


# ---------------------------------------------------------------------------
# Time-domain evolution
# ---------------------------------------------------------------------------

def evolve(r_grid, V, hk, hkp, T_final, dt, psi0, psidot0=None, i_obs=None):
    """Leapfrog time evolution of ψ(t, r).

    Wave eq: ∂²ψ/∂t² = hk ∂²ψ/∂r² + (1/2) hk' ∂ψ/∂r - V ψ.

    Returns: time array, signal at i_obs (ψ vs t).
    """
    dr = r_grid[1] - r_grid[0]
    N_steps = int(T_final / dt)
    psi = psi0.copy()
    if psidot0 is None:
        psidot0 = np.zeros_like(psi0)
    psi_old = psi - dt * psidot0  # one step backward for leapfrog

    # Precompute coefficients
    inv_dr2 = 1.0 / (dr * dr)
    inv_2dr = 1.0 / (2.0 * dr)

    t_arr = np.zeros(N_steps)
    sig = np.zeros(N_steps)

    for step in range(N_steps):
        # Spatial derivatives (central differences, interior)
        psi_rr = np.zeros_like(psi)
        psi_r = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) * inv_dr2
        psi_r[1:-1] = (psi[2:] - psi[:-2]) * inv_2dr

        # Wave equation
        psi_tt = hk * psi_rr + 0.5 * hkp * psi_r - V * psi

        # Leapfrog update
        psi_new = 2 * psi - psi_old + dt * dt * psi_tt

        # Absorbing boundary conditions (Sommerfeld-like, copy outgoing)
        psi_new[0] = psi[1]   # near horizon: outflow
        psi_new[-1] = psi[-2]  # far field: outflow

        psi_old = psi
        psi = psi_new

        t_arr[step] = (step + 1) * dt
        if i_obs is not None:
            sig[step] = psi[i_obs]

    return t_arr, sig


# ---------------------------------------------------------------------------
# QNM extraction via late-time damped sinusoid fit
# ---------------------------------------------------------------------------

def damped_sinusoid(t, A, omega_r, omega_i, phi, offset):
    return A * np.exp(-omega_i * (t - t[0])) * np.cos(omega_r * (t - t[0]) + phi) + offset


def extract_QNM(t, signal, t_fit_start, t_fit_end,
                init_omega_r=0.4, init_omega_i=0.1):
    """Fit late-time signal to damped sinusoid.  Returns omega = omega_r - i omega_i."""
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    t_fit = t[mask]
    s_fit = signal[mask]

    # Initial guess
    A_guess = (np.max(s_fit) - np.min(s_fit)) / 2
    p0 = [A_guess, init_omega_r, init_omega_i, 0.0, 0.0]

    try:
        popt, pcov = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0, maxfev=20000)
        A_fit, omega_r_fit, omega_i_fit, phi_fit, off_fit = popt
        omega_r_fit = abs(omega_r_fit)  # absorb sign into phi
        omega_i_fit = abs(omega_i_fit)  # damping is positive
        omega_qnm = omega_r_fit - 1j * omega_i_fit  # convention Im < 0 for damped
        # Standard deviations
        perr = np.sqrt(np.diag(pcov))
        return omega_qnm, popt, perr
    except Exception as ex:
        print(f"  curve_fit failed: {ex}")
        return None, None, None


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

print("=" * 80)
print("G88: Time-domain QNM extraction for ℓ = 2, n = 0")
print("=" * 80)
print()
print("Method: leapfrog evolution of ∂²ψ/∂t² = hk ψ_rr + (1/2)hk' ψ_r - V ψ")
print("Comparison: Schwarzschild vs STAM proxy potential")
print()

# Grid setup
M = 1.0
r_min, r_max = 2.05, 60.0
N_grid = 3000
r_grid = np.linspace(r_min, r_max, N_grid)
dr = r_grid[1] - r_grid[0]
ell = 2
print(f"Grid built.", flush=True)

# Precompute coefficients for both
hk_S = h_fn(r_grid) * k_Schw(r_grid)
hk_T = h_fn(r_grid) * k_STAM(r_grid)
hkp_S = np.gradient(hk_S, dr)
hkp_T = np.gradient(hk_T, dr)
V_S = V_RW_proxy(r_grid, ell, k_Schw, M)
V_T = V_RW_proxy(r_grid, ell, k_STAM, M)

# Initial condition: Gaussian centered at r = 7M (well outside PS)
r_pulse = 7.0
sigma_pulse = 1.0
psi0 = np.exp(-(r_grid - r_pulse)**2 / (2 * sigma_pulse**2))

# Observation: read off at r = 20M
r_obs = 20.0
i_obs = int(np.argmin(np.abs(r_grid - r_obs)))
print(f"Grid: N = {N_grid}, r ∈ [{r_min}, {r_max}] M, dr = {dr:.4e}")
print(f"Initial pulse: Gaussian at r = {r_pulse} M, σ = {sigma_pulse}")
print(f"Observation: r = {r_grid[i_obs]:.4f} M (index {i_obs})")
print()

# Time step (CFL): dt < dr / max|sqrt(hk)|.  max(hk) ~ 1, so dt < dr.
dt = 0.4 * dr
T_final = 200.0
N_steps_est = int(T_final / dt)
print(f"Time evolution: T_final = {T_final} M, dt = {dt:.4e}, N_steps = {N_steps_est}", flush=True)
print(flush=True)

print("Evolving Schwarzschild ...", flush=True)
t_S, sig_S = evolve(r_grid, V_S, hk_S, hkp_S, T_final, dt, psi0, i_obs=i_obs)
print(f"  Done.  max |ψ| = {np.max(np.abs(sig_S)):.4e}", flush=True)

print("Evolving STAM ...", flush=True)
t_T, sig_T = evolve(r_grid, V_T, hk_T, hkp_T, T_final, dt, psi0, i_obs=i_obs)
print(f"  Done.  max |ψ| = {np.max(np.abs(sig_T)):.4e}", flush=True)
print(flush=True)


# Extract QNM from late-time fit
# Skip initial transient (waves crossing the potential), then fit
# Fit window: avoid early transient (t < 40 M) and very late-time tail (t > 150 M)
t_fit_start = 40.0
t_fit_end = 130.0

print("=" * 80)
print(f"QNM extraction (fit window: t ∈ [{t_fit_start}, {t_fit_end}] M)")
print("=" * 80)
print()

omega_S, popt_S, perr_S = extract_QNM(t_S, sig_S, t_fit_start, t_fit_end,
                                       init_omega_r=0.37, init_omega_i=0.09)
omega_T, popt_T, perr_T = extract_QNM(t_T, sig_T, t_fit_start, t_fit_end,
                                       init_omega_r=0.37, init_omega_i=0.09)

# Gold standard for Schwarzschild ℓ=2, n=0 (Leaver):
omega_gold = 0.373672 - 0.088962j

print(f"Schwarzschild  ℓ = 2, n = 0:")
print(f"  ω (time-domain) = {omega_S}")
print(f"  ω (Leaver gold) = {omega_gold}")
if omega_S is not None:
    rel_err_S = abs(omega_S - omega_gold) / abs(omega_gold)
    print(f"  Method calibration error |Δω/ω| = {rel_err_S * 100:.3f} %")
    if perr_S is not None:
        print(f"  Fit uncertainties: Δω_r ~ {perr_S[1]:.4e}, Δω_i ~ {perr_S[2]:.4e}")
print()

print(f"STAM           ℓ = 2, n = 0:")
print(f"  ω (time-domain) = {omega_T}")
if omega_T is not None:
    if perr_T is not None:
        print(f"  Fit uncertainties: Δω_r ~ {perr_T[1]:.4e}, Δω_i ~ {perr_T[2]:.4e}")
print()

if omega_S is not None and omega_T is not None:
    dw_over_w = abs(omega_T - omega_S) / abs(omega_S)
    print(f"Framework shift: |ω_STAM - ω_Schw| / |ω_Schw| = {dw_over_w * 100:.2f} %")
    print()


# ---------------------------------------------------------------------------
# Diagnostic: Prony / multi-mode fit as cross-check
# ---------------------------------------------------------------------------

def matrix_pencil_qnm(t, signal, n_modes=2):
    """Simple matrix pencil method for extracting damped oscillation modes."""
    M = len(signal)
    L = M // 2
    # Build Hankel matrices
    H0 = np.array([signal[i:i+L] for i in range(M - L)])
    H1 = np.array([signal[i+1:i+L+1] for i in range(M - L)])
    # Solve generalized eigenvalue problem H1 = H0 Z
    # Use SVD-based pseudo-inverse
    U, s, Vt = np.linalg.svd(H0, full_matrices=False)
    # Truncate to n_modes
    s_inv = np.zeros_like(s)
    s_inv[:n_modes] = 1.0 / s[:n_modes]
    H0_inv = Vt.T @ np.diag(s_inv) @ U.T
    Z = H0_inv @ H1
    eigvals = np.linalg.eigvals(Z)
    # Convert: signal mode Z_k = exp(s_k dt), so s_k = log(Z_k)/dt
    dt_local = t[1] - t[0]
    poles = np.log(eigvals) / dt_local
    # Sort by damping (smaller |Im(s)|)
    poles = sorted(poles, key=lambda p: abs(p.real))
    # In our convention: s_k = -ω_i + i ω_r (since signal ~ exp(-ω_i t) cos(ω_r t))
    # So ω = ω_r - i ω_i = Im(s_k) - i (-Re(s_k)) = Im(s) + i Re(s).
    # ie omega = ω_r - i ω_i;  with s_k = -ω_i + i ω_r, so Im(s_k) = ω_r, Re(s_k) = -ω_i
    # ω = ω_r - i ω_i = Im(s_k) - i(-Re(s_k)) = Im(s_k) + i Re(s_k).  Hmm let me think again.
    # signal ~ exp(s t) = exp((-ω_i + i ω_r) t) = exp(-ω_i t)(cos(ω_r t) + i sin(ω_r t))
    # Real part: exp(-ω_i t) cos(ω_r t)
    # So s = -ω_i + i ω_r, and we want ω = ω_r - i ω_i (Im < 0 for damped).
    # ω = Im(s) - i (-Re(s)) ... = Im(s) + i Re(s).  Re(s) < 0 for damped, so Im(ω) < 0. ✓
    omega_modes = [pole.imag + 1j * pole.real for pole in poles]
    return omega_modes


print("=" * 80)
print("Cross-check: matrix-pencil mode extraction (2-mode fit on fit window)")
print("=" * 80)
print()

mask = (t_S >= t_fit_start) & (t_S <= t_fit_end)
t_fit_arr = t_S[mask]
sig_S_fit = sig_S[mask]
sig_T_fit = sig_T[mask]

try:
    modes_S = matrix_pencil_qnm(t_fit_arr, sig_S_fit, n_modes=4)
    modes_T = matrix_pencil_qnm(t_fit_arr, sig_T_fit, n_modes=4)
    print("Schwarzschild dominant modes (sorted by damping):")
    for m in modes_S[:4]:
        if 0.1 < m.real < 2 and -1 < m.imag < 0:
            print(f"  ω = {m.real:.5f} {'+' if m.imag>=0 else '-'} {abs(m.imag):.5f} i")
    print()
    print("STAM dominant modes (sorted by damping):")
    for m in modes_T[:4]:
        if 0.1 < m.real < 2 and -1 < m.imag < 0:
            print(f"  ω = {m.real:.5f} {'+' if m.imag>=0 else '-'} {abs(m.imag):.5f} i")
    print()
except Exception as ex:
    print(f"  Matrix pencil failed: {ex}")
    print()


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Full time series, log scale
ax = axes[0, 0]
ax.plot(t_S, np.abs(sig_S), 'tab:blue', linewidth=1, alpha=0.8, label='Schwarzschild')
ax.plot(t_T, np.abs(sig_T), 'tab:orange', linewidth=1, alpha=0.8, label='STAM')
ax.set_xlabel('t / M')
ax.set_ylabel(r'$|\psi(t, r_{\rm obs})|$')
ax.set_title('Time-domain signal at observation point (log scale)')
ax.set_yscale('log')
ax.axvline(t_fit_start, color='red', linestyle=':', alpha=0.5, label='fit window')
ax.axvline(t_fit_end, color='red', linestyle=':', alpha=0.5)
ax.legend()
ax.grid(True, alpha=0.3)

# (2) Fit window with damped-sinusoid fit overlay
ax = axes[0, 1]
mask = (t_S >= t_fit_start) & (t_S <= t_fit_end)
ax.plot(t_S[mask], sig_S[mask], 'tab:blue', linewidth=1, label='Schw signal')
if popt_S is not None:
    fit_S = damped_sinusoid(t_S[mask], *popt_S)
    ax.plot(t_S[mask], fit_S, 'k--', linewidth=1, alpha=0.7, label='Schw fit')
ax.plot(t_T[mask], sig_T[mask], 'tab:orange', linewidth=1, label='STAM signal')
if popt_T is not None:
    fit_T = damped_sinusoid(t_T[mask], *popt_T)
    ax.plot(t_T[mask], fit_T, 'r--', linewidth=1, alpha=0.7, label='STAM fit')
ax.set_xlabel('t / M')
ax.set_ylabel(r'$\psi(t, r_{\rm obs})$')
ax.set_title('Fit window with damped-sinusoid overlay')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (3) Potential comparison
ax = axes[1, 0]
ax.plot(r_grid, V_S, 'tab:blue', linewidth=2, label='Schw V_RW^proxy')
ax.plot(r_grid, V_T, 'tab:orange', linewidth=2, label='STAM V_RW^proxy')
ax.axvline(3, color='black', linestyle=':', alpha=0.5, label='PS r=3M')
ax.axvline(2, color='red', linestyle=':', alpha=0.5, label='horizon r=2M')
ax.set_xlabel('r / M')
ax.set_ylabel('V_RW^proxy (ℓ=2)')
ax.set_title('Tensor RW proxy potential, ℓ = 2')
ax.set_xlim(2, 10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# (4) Summary text
ax = axes[1, 1]
ax.axis('off')
txt_lines = [
    "G88 results: time-domain QNM extraction",
    "",
    f"Method: leapfrog evolution, fit damped sinusoid in t ∈ [{t_fit_start}, {t_fit_end}] M",
    f"Grid: N = {N_grid}, dr = {dr:.4e}, dt = {dt:.4e}",
    "",
    f"Schwarzschild ℓ = 2, n = 0:",
    f"  Time-domain:  ω = {omega_S}",
    f"  Leaver gold:  ω = {omega_gold}",
]
if omega_S is not None:
    txt_lines.append(f"  Method calibration: |Δω/ω| = {abs(omega_S - omega_gold)/abs(omega_gold)*100:.2f} %")
txt_lines += ["", f"STAM ℓ = 2, n = 0:"]
if omega_T is not None:
    txt_lines.append(f"  Time-domain:  ω = {omega_T}")
txt_lines.append("")
if omega_S is not None and omega_T is not None:
    txt_lines.append(f"Framework shift: |Δω/ω| = {abs(omega_T - omega_S)/abs(omega_S)*100:.2f} %")
ax.text(0.05, 0.95, '\n'.join(txt_lines),
        transform=ax.transAxes, fontsize=10, family='monospace',
        verticalalignment='top')

plt.tight_layout()
out_png = PLOTS / "G88_time_domain_QNM.png"
plt.savefig(out_png, dpi=200)
plt.close()
print(f"Plot saved: {out_png}")
print()


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

print("=" * 80)
print("VERDICT")
print("=" * 80)
print()

# Method calibration
if omega_S is not None:
    rel_err = abs(omega_S - omega_gold) / abs(omega_gold)
    print(f"Method calibration (Schwarzschild vs Leaver gold):")
    print(f"  |Δω/ω|_calibration = {rel_err * 100:.2f} %")
    if rel_err < 0.05:
        print("  -> Method is reliable to ~5% for this grid.")
    else:
        print("  -> Larger calibration error; framework prediction has corresponding uncertainty.")
    print()

if omega_S is not None and omega_T is not None:
    dw = abs(omega_T - omega_S) / abs(omega_S)
    print(f"Framework prediction:")
    print(f"  ω^Schw  = {omega_S}")
    print(f"  ω^STAM  = {omega_T}")
    print(f"  |Δω/ω|  = {dw * 100:.2f} %")
    print()
    if dw > 0.10:
        print("  =>  > 10%:  FRAMEWORK HAS A NEAR-TERM LIGO-FACING PREDICTION.")
        print("              The ℓ = 2, n = 0 ringdown frequency on the STAM proxy")
        print("              potential differs from Schwarzschild at the level of")
        print("              current observational precision.")
    elif dw > 0.01:
        print("  =>  1-10%:  Framework predicts a percent-level deviation at ℓ = 2.")
        print("              Within reach of high-precision LIGO/LISA ringdown analysis.")
    else:
        print("  =>  < 1%:   Framework matches Schwarzschild ringdown to high precision.")
        print("              G86 / G87 large WKB shifts were WKB low-ℓ noise.")
    print()

print("Caveat:  this uses the TENSOR PROXY potential, not the rigorous axial")
print("perturbation potential for the constrained-Σ action.  The proxy captures")
print("how QNM responds to the framework's k(r) modification via the standard")
print("M_eff formulation.  A rigorous derivation (G89+) would include f(Σ)R-")
print("coupling corrections and could shift the numerical answer; the structural")
print("conclusion (eikonal matches, deviation at higher ℓ derivatives) is robust.")
print()


# ---------------------------------------------------------------------------
# Summary file
# ---------------------------------------------------------------------------

md = []
md.append("# G88 — Time-domain QNM extraction for ℓ = 2, n = 0\n")
md.append("**Date: 2026-05-13.**  Phase 2 of the QNM track: directly evolve the wave "
          "equation on the framework's tensor RW proxy potential and extract the "
          "ringdown frequency.\n")
md.append("## Method\n")
md.append("Leapfrog time-domain evolution of\n")
md.append("```\n"
          "∂²ψ/∂t² = h(r) k(r) ∂²ψ/∂r² + (1/2)(d(hk)/dr) ∂ψ/∂r - V_RW^proxy(r) ψ\n"
          "```\n")
md.append("with absorbing boundary conditions, Gaussian initial pulse at r = 7M, "
          "and signal extraction at r = 20M. Late-time signal fit to "
          "ψ(t) = A exp(-ω_i t) cos(ω_r t + φ).\n")
md.append("Grid: N = {}, r ∈ [{}, {}] M.  Time: T_final = {} M, dt = {:.4e}.\n".format(
    N_grid, r_min, r_max, T_final, dt))
md.append("\n")
md.append("## Results\n")
md.append("| Background | ω (time-domain) | Leaver gold (Schw) | Method error |\n|---|---|---|---:|\n")
if omega_S is not None:
    rel_err = abs(omega_S - omega_gold) / abs(omega_gold)
    md.append(f"| Schwarzschild | {omega_S.real:.5f} {'−' if omega_S.imag < 0 else '+'} "
              f"{abs(omega_S.imag):.5f} i | {omega_gold.real:.5f} − {abs(omega_gold.imag):.5f} i | "
              f"{rel_err*100:.2f} % |\n")
if omega_T is not None:
    md.append(f"| **STAM proxy** | {omega_T.real:.5f} {'−' if omega_T.imag < 0 else '+'} "
              f"{abs(omega_T.imag):.5f} i | — | — |\n")
md.append("\n")
md.append("## Framework prediction\n")
if omega_S is not None and omega_T is not None:
    dw_pct = abs(omega_T - omega_S) / abs(omega_S) * 100
    md.append(f"**|Δω/ω| at ℓ = 2, n = 0:  {dw_pct:.2f} %**\n")
    md.append("\n")
    if dw_pct > 10:
        md.append("**Verdict:** > 10% → framework has a near-term LIGO-facing prediction.\n")
    elif dw_pct > 1:
        md.append("**Verdict:** 1-10% → percent-level deviation; within reach of "
                  "high-precision LIGO/LISA ringdown analysis.\n")
    else:
        md.append("**Verdict:** < 1% → framework matches Schwarzschild ringdown to "
                  "high precision; G86/G87 WKB shifts were low-ℓ WKB noise.\n")
md.append("\n")
md.append("## Caveat: proxy vs rigorous derivation\n")
md.append("This uses the tensor RW PROXY potential V_RW^proxy(r) = h(r)[ℓ(ℓ+1)/r² - "
          "6 M_eff(r)/r³] with M_eff = (r/2)(1-k(r)). The rigorous axial perturbation "
          "potential for the constrained-Σ action (linearizing f(Σ)R + λ_1((∇Σ)² - W) "
          "+ λ_2(u^μ ∂_μΣ) - 2V(Σ) around the background) is a separate computation "
          "(G89+). The proxy captures how QNM respond to the framework's k(r) "
          "modification through the standard M_eff substitution; the rigorous "
          "derivation could shift the numerical answer but the structural "
          "conclusions (eikonal match, deviation at higher derivatives of V) are "
          "robust.\n")
md.append("\n")
md.append("## Files\n")
md.append("- [scripts/G88_leaver_proxy_QNM_l2.py](../scripts/G88_leaver_proxy_QNM_l2.py)\n")
md.append("- [plots/G88_time_domain_QNM.png](../plots/G88_time_domain_QNM.png)\n")

out_md = RESULTS / "G88_time_domain_QNM_summary.md"
out_md.write_text("\n".join(md), encoding="utf-8")
print(f"Summary written: {out_md}")
