#!/usr/bin/env python3
"""
G92_rigorous_axial_QNM.py

Derive and run the axial odd-parity tensor perturbation equation implied by
the constrained shell-count action:

    S = (1/16 pi G) int sqrt(-g) [
          f(Sigma) R
        + lambda_1 ((grad Sigma)^2 - W(Sigma))
        + lambda_2 (u^mu d_mu Sigma)
        - 2 V(Sigma)
    ] d^4x

for the committed spinless STAM background:

    ds^2 = -h(r) dt^2 + dr^2/k(r) + r^2 dOmega^2
    h(r) = 1 - 2M/r
    k(r) = h(r)                         outside photon sphere
    k(r) = h(r) F(y), y = 6M/r - 2      inside photon sphere
    F(y) = 1 - 5y^4 + 4y^5

Odd-parity derivation
=====================
Odd parity has no scalar harmonic with the parity needed to perturb Sigma.
The double-LM constraints from G79/G80 therefore set

    delta Sigma = delta lambda_1 = delta lambda_2 = 0

in the propagating axial tensor sector.  The LM and potential terms only
support the background.  The propagating axial degree of freedom comes from
the tensor part of f(Sigma) R.

Let Phi be the 2+2 odd-parity axial master amplitude before canonical
normalization.  Linearizing the odd tensor equations gives the Sturm-Liouville
form

    -d_t^2 Phi + d_*^2 Phi + P d_* Phi - V_geom Phi = 0

where d_* = sqrt(h k) d_r and

    P = d_*(ln f).

The geometry term is the standard covariant axial potential on a static
spherical background, written directly in h and k:

    V_geom = h [ ell(ell+1)/r^2
               - 2(1-k)/r^2
               - k'/(2r)
               - k h'/(2 h r) ].

For k = h = 1 - 2M/r this reduces exactly to the Schwarzschild
Regge-Wheeler potential h[ell(ell+1)/r^2 - 6M/r^3].

Canonical field

    Psi = sqrt(f) Phi

removes the first-derivative term and gives the rigorous axial tensor
equation used below:

    d_t^2 Psi - d_*^2 Psi + V_ax Psi = 0

    V_ax = V_geom + (1/2) d_* P + (1/4) P^2.

The f(Sigma) profile is the G70/G75 matched nonminimal coupling, anchored by
f(A=2/3) = 1 and f = 1 outside the photon sphere:

    d(ln f)/dA = -2 y^3 (45A^2 - 33A - 13) / [A F(y)].

The script then runs the same calibrated time-domain r* pipeline as G91:
first Schwarzschild calibration against Leaver gold, then STAM only if the
calibration error is below 1%.
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

M = 1.0


# ---------------------------------------------------------------------------
# Metric, shell-count, and nonminimal-coupling functions
# ---------------------------------------------------------------------------

def h_fn(r):
    return 1.0 - 2.0 * M / r


def hp_fn(r):
    return 2.0 * M / r**2


def y_of_r(r):
    return 6.0 * M / r - 2.0


def A_of_r(r):
    return 2.0 * M / r


def F_quintic(y):
    return 1.0 - 5.0 * y**4 + 4.0 * y**5


def dF_dy(y):
    return -20.0 * y**3 + 20.0 * y**4


def k_GR(r):
    return h_fn(r)


def kp_GR(r):
    return hp_fn(r)


def k_STAM(r):
    y = y_of_r(r)
    h_val = h_fn(r)
    if y > 0:
        return h_val * F_quintic(y)
    return h_val


def kp_STAM(r):
    y = y_of_r(r)
    h_val = h_fn(r)
    hp_val = hp_fn(r)
    if y > 0:
        yp = -6.0 * M / r**2
        return hp_val * F_quintic(y) + h_val * dF_dy(y) * yp
    return hp_val


def dlnf_dA(A):
    """G70/G75 nonminimal coupling derivative, zero outside PS."""
    if np.isscalar(A):
        if A <= 2.0 / 3.0:
            return 0.0
        y = 3.0 * A - 2.0
        F = F_quintic(y)
        return -2.0 * y**3 * (45.0 * A**2 - 33.0 * A - 13.0) / (A * F)

    A_arr = np.asarray(A, dtype=float)
    out = np.zeros_like(A_arr)
    mask = A_arr > 2.0 / 3.0
    y = 3.0 * A_arr[mask] - 2.0
    F = F_quintic(y)
    out[mask] = -2.0 * y**3 * (45.0 * A_arr[mask]**2 - 33.0 * A_arr[mask] - 13.0) / (
        A_arr[mask] * F
    )
    return out


def P_planck_mass(r, k_fn):
    """P = d_*(ln f) = sqrt(h k) d_r ln f."""
    A = A_of_r(r)
    if A <= 2.0 / 3.0:
        return 0.0
    dA_dr = -2.0 * M / r**2
    return np.sqrt(h_fn(r) * k_fn(r)) * dlnf_dA(A) * dA_dr


def dstar_P_numeric(r, k_fn):
    """d_* P = sqrt(h k) dP/dr, using a local finite difference."""
    if A_of_r(r) <= 2.0 / 3.0:
        return 0.0
    step = max(1e-5, 1e-5 * r)
    r_lo = max(2.0 * M + 1e-8, r - step)
    r_hi = min(3.0 * M - 1e-8, r + step)
    if r_hi <= r_lo:
        return 0.0
    dP_dr = (P_planck_mass(r_hi, k_fn) - P_planck_mass(r_lo, k_fn)) / (r_hi - r_lo)
    return np.sqrt(h_fn(r) * k_fn(r)) * dP_dr


# ---------------------------------------------------------------------------
# Axial potentials
# ---------------------------------------------------------------------------

def V_axial_geom(r, ell, k_fn, kp_fn):
    h_val = h_fn(r)
    k_val = k_fn(r)
    hp_val = hp_fn(r)
    kp_val = kp_fn(r)
    return h_val * (
        ell * (ell + 1.0) / r**2
        - 2.0 * (1.0 - k_val) / r**2
        - kp_val / (2.0 * r)
        - k_val * hp_val / (2.0 * h_val * r)
    )


def V_axial_rigorous(r, ell, k_fn, kp_fn):
    P = P_planck_mass(r, k_fn)
    return V_axial_geom(r, ell, k_fn, kp_fn) + 0.5 * dstar_P_numeric(r, k_fn) + 0.25 * P**2


def V_axial_GR(r, ell):
    """GR calibration potential.  Here f = 1, so P = d_* ln f = 0."""
    return V_axial_geom(r, ell, k_GR, kp_GR)


def V_RW_schwarzschild(r, ell):
    h_val = h_fn(r)
    return h_val * (ell * (ell + 1.0) / r**2 - 6.0 * M / r**3)


# ---------------------------------------------------------------------------
# Tortoise coordinate construction
# ---------------------------------------------------------------------------

def r_of_rstar_GR(rstar):
    PS_offset = 3.0 * M + 2.0 * M * np.log(0.5)
    target = rstar + PS_offset

    def f_root(r):
        return r + 2.0 * M * np.log(r / (2.0 * M) - 1.0) - target

    if rstar > 0:
        return brentq(f_root, 3.0 * M - 1e-10, 1e6)
    return brentq(f_root, 2.0 * M + 1e-14, 3.0 * M + 1e-10)


def build_r_of_rstar_STAM(rstar_grid):
    def rhs(rs, y):
        r_val = y[0]
        return [np.sqrt(h_fn(r_val) * k_STAM(r_val))]

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
                        method="RK45", rtol=1e-12, atol=1e-14, max_step=0.05)
        vals = np.zeros_like(rs_neg)
        vals[order] = sol.y[0]
        r_arr[rstar_grid < 0] = vals

    r_arr[rstar_grid == 0] = 3.0 * M
    return r_arr


# ---------------------------------------------------------------------------
# Time-domain evolution and QNM fit
# ---------------------------------------------------------------------------

def evolve_TD(V_grid, dr_star, T_final, dt, psi0, i_obs):
    psi_old = psi0.copy()
    psi = psi0.copy()
    n_steps = int(T_final / dt)
    inv_dr2 = 1.0 / dr_star**2
    cfl = dt / dr_star

    t_arr = np.zeros(n_steps)
    sig_arr = np.zeros(n_steps)

    for step in range(n_steps):
        psi_rr = np.zeros_like(psi)
        psi_rr[1:-1] = (psi[2:] - 2.0 * psi[1:-1] + psi[:-2]) * inv_dr2

        psi_new = np.empty_like(psi)
        psi_new[1:-1] = 2.0 * psi[1:-1] - psi_old[1:-1] + dt**2 * (
            psi_rr[1:-1] - V_grid[1:-1] * psi[1:-1]
        )
        psi_new[-1] = psi[-1] - cfl * (psi[-1] - psi[-2])
        psi_new[0] = psi[0] - cfl * (psi[0] - psi[1])

        psi_old = psi
        psi = psi_new
        t_arr[step] = (step + 1) * dt
        sig_arr[step] = psi[i_obs]

    return t_arr, sig_arr


def damped_sinusoid(t, amp, omega_r, omega_i, phi, offset):
    tau = t - t[0]
    return amp * np.exp(-omega_i * tau) * np.cos(omega_r * tau + phi) + offset


def extract_QNM(t, signal, t_fit_start, t_fit_end,
                init_omega_r=0.37, init_omega_i=0.09):
    mask = (t >= t_fit_start) & (t <= t_fit_end)
    t_fit = t[mask]
    s_fit = signal[mask]
    amp_guess = 0.5 * (np.max(s_fit) - np.min(s_fit))
    p0 = [amp_guess, init_omega_r, init_omega_i, 0.0, 0.0]
    bounds = ([-np.inf, 0.0, 0.0, -np.pi, -np.inf],
              [np.inf, 2.0, 2.0, np.pi, np.inf])
    popt, _ = curve_fit(damped_sinusoid, t_fit, s_fit, p0=p0,
                        bounds=bounds, maxfev=30000)
    return complex(popt[1], -popt[2]), popt


# ---------------------------------------------------------------------------
# Run G92
# ---------------------------------------------------------------------------

print("=" * 88)
print("G92: rigorous axial perturbation equation from constrained shell-count action")
print("=" * 88)
print()
print("Odd sector: delta Sigma = 0; LM constraints carry no propagating axial source.")
print("Canonical potential: V_ax = V_geom + 0.5 d_*P + 0.25 P^2, P = d_* ln f.")
print()

ell = 2
omega_gold = complex(0.373672, -0.088962)
t_fit_start = 100.0
t_fit_end = 250.0
T_final = 350.0
rs_pulse = 30.0
sigma_pulse = 3.0

print("Sanity check: GR geometric axial potential vs Schwarzschild RW")
for rr in [2.2, 3.0, 6.0, 20.0]:
    diff = V_axial_geom(rr, ell, k_GR, kp_GR) - V_RW_schwarzschild(rr, ell)
    print(f"  r={rr:5.2f}: V_geom - V_RW = {diff:+.6e}")
print()

print("=" * 80)
print("G92a: Schwarzschild calibration")
print("=" * 80)
rs_min_GR = -60.0
rs_max_GR = 400.0
N_GR = 6001
rs_grid_GR = np.linspace(rs_min_GR, rs_max_GR, N_GR)
dr_star_GR = rs_grid_GR[1] - rs_grid_GR[0]
print(f"Grid: N={N_GR}, r* in [{rs_min_GR}, {rs_max_GR}], dr*={dr_star_GR:.4e}")
r_grid_GR = np.array([r_of_rstar_GR(rs) for rs in rs_grid_GR])
V_grid_GR = np.array([V_axial_GR(r, ell) for r in r_grid_GR])
print(f"  V max at r*={rs_grid_GR[np.argmax(V_grid_GR)]:.4f}, r={r_grid_GR[np.argmax(V_grid_GR)]:.4f}")
psi0_GR = np.exp(-(rs_grid_GR - rs_pulse)**2 / (2.0 * sigma_pulse**2))
i_obs_GR = int(np.argmin(np.abs(rs_grid_GR - 50.0)))
dt_GR = 0.5 * dr_star_GR
print(f"  Evolving: T_final={T_final}, dt={dt_GR:.4e}, N_steps={int(T_final / dt_GR)}")
t_GR, sig_GR = evolve_TD(V_grid_GR, dr_star_GR, T_final, dt_GR, psi0_GR, i_obs_GR)
omega_GR, popt_GR = extract_QNM(t_GR, sig_GR, t_fit_start, t_fit_end)
err_GR = abs(omega_GR - omega_gold) / abs(omega_gold)
calib_ok = err_GR < 0.01
print(f"  omega_GR   = {omega_GR}")
print(f"  omega_gold = {omega_gold}")
print(f"  |Delta omega / omega| = {err_GR * 100:.4f}%")
print(f"  calibration {'PASS' if calib_ok else 'FAIL'}")
print()

if not calib_ok:
    raise SystemExit("Schwarzschild calibration failed; not running STAM.")

print("=" * 80)
print("G92b: STAM rigorous axial run")
print("=" * 80)
rs_min_STAM = -500.0
rs_max_STAM = 400.0
N_STAM = 10001
rs_grid_STAM = np.linspace(rs_min_STAM, rs_max_STAM, N_STAM)
dr_star_STAM = rs_grid_STAM[1] - rs_grid_STAM[0]
print(f"Grid: N={N_STAM}, r* in [{rs_min_STAM}, {rs_max_STAM}], dr*={dr_star_STAM:.4e}")
r_grid_STAM = build_r_of_rstar_STAM(rs_grid_STAM)
eps_min = (r_grid_STAM.min() - 2.0 * M) / (2.0 * M)
V_grid_STAM = np.array([V_axial_rigorous(r, ell, k_STAM, kp_STAM) for r in r_grid_STAM])
print(f"  r range: [{r_grid_STAM.min():.6f}, {r_grid_STAM.max():.6f}], eps_min={eps_min:.4e}")
print(f"  V max at r*={rs_grid_STAM[np.argmax(V_grid_STAM)]:.4f}, r={r_grid_STAM[np.argmax(V_grid_STAM)]:.4f}")
psi0_STAM = np.exp(-(rs_grid_STAM - rs_pulse)**2 / (2.0 * sigma_pulse**2))
i_obs_STAM = int(np.argmin(np.abs(rs_grid_STAM - 50.0)))
dt_STAM = 0.5 * dr_star_STAM
print(f"  Evolving: T_final={T_final}, dt={dt_STAM:.4e}, N_steps={int(T_final / dt_STAM)}")
t_STAM, sig_STAM = evolve_TD(V_grid_STAM, dr_star_STAM, T_final, dt_STAM,
                             psi0_STAM, i_obs_STAM)
omega_STAM, popt_STAM = extract_QNM(t_STAM, sig_STAM, t_fit_start, t_fit_end)
dw_over_w = abs(omega_STAM - omega_GR) / abs(omega_GR)
re_shift = (omega_STAM.real - omega_GR.real) / omega_GR.real
im_shift = (omega_STAM.imag - omega_GR.imag) / omega_GR.imag
print(f"  omega_STAM = {omega_STAM}")
print(f"  |Delta omega / omega| = {dw_over_w * 100:.4f}% vs Schwarzschild")
print(f"  Re shift = {re_shift * 100:+.4f}%")
print(f"  Im shift = {im_shift * 100:+.4f}%")
print()


# ---------------------------------------------------------------------------
# Plot and summary
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax = axes[0, 0]
ax.plot(rs_grid_GR, V_grid_GR, color="tab:blue", linewidth=2, label="GR axial")
ax.plot(rs_grid_STAM, V_grid_STAM, color="tab:orange", linewidth=2, label="STAM rigorous axial")
ax.set_xlim(-50, 100)
ax.set_xlabel("r* / M")
ax.set_ylabel("V_axial")
ax.set_title("Rigorous axial potential vs tortoise coordinate")
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(t_GR, np.abs(sig_GR), color="tab:blue", linewidth=1, label="Schwarzschild")
ax.plot(t_STAM, np.abs(sig_STAM), color="tab:orange", linewidth=1, label="STAM")
ax.axvspan(t_fit_start, t_fit_end, alpha=0.15, color="gray", label="fit window")
ax.set_yscale("log")
ax.set_xlabel("t / M")
ax.set_ylabel("|psi(t, r_obs)|")
ax.set_title("Time-domain axial signal")
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
mask_GR = (t_GR >= t_fit_start) & (t_GR <= t_fit_end)
mask_STAM = (t_STAM >= t_fit_start) & (t_STAM <= t_fit_end)
ax.plot(t_GR[mask_GR], sig_GR[mask_GR], color="tab:blue", linewidth=1, label="GR")
ax.plot(t_GR[mask_GR], damped_sinusoid(t_GR[mask_GR], *popt_GR),
        "k--", alpha=0.7, label="GR fit")
ax.plot(t_STAM[mask_STAM], sig_STAM[mask_STAM], color="tab:orange", linewidth=1, label="STAM")
ax.plot(t_STAM[mask_STAM], damped_sinusoid(t_STAM[mask_STAM], *popt_STAM),
        "r--", alpha=0.7, label="STAM fit")
ax.set_xlabel("t / M")
ax.set_ylabel("psi(t, r_obs)")
ax.set_title("Fit window")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.axis("off")
txt = [
    "G92 rigorous axial equation",
    "",
    "Schwarzschild calibration:",
    f"  omega_TD     = {omega_GR.real:.6f} - {abs(omega_GR.imag):.6f} i",
    f"  omega_Leaver = {omega_gold.real:.6f} - {abs(omega_gold.imag):.6f} i",
    f"  error        = {err_GR * 100:.4f} %",
    "",
    "STAM rigorous axial:",
    f"  omega_TD     = {omega_STAM.real:.6f} - {abs(omega_STAM.imag):.6f} i",
    f"  shift        = {dw_over_w * 100:.4f} %",
    f"  Re shift     = {re_shift * 100:+.4f} %",
    f"  Im shift     = {im_shift * 100:+.4f} %",
    "",
    "V_ax = V_geom + 0.5 d_*P + 0.25 P^2",
    "P = d_* ln f",
]
ax.text(0.05, 0.95, "\n".join(txt), transform=ax.transAxes,
        fontsize=10, family="monospace", verticalalignment="top")

plt.tight_layout()
out_png = PLOTS / "G92_rigorous_axial_QNM.png"
plt.savefig(out_png, dpi=200)
plt.close()
print(f"Plot saved: {out_png}")

md = []
md.append("# G92 - Rigorous axial perturbation equation from constrained shell-count action\n")
md.append("**Date: 2026-05-13.**\n")
md.append("## Derivation result\n")
md.append("Odd parity does not perturb the shell-count field: `delta Sigma = 0`. The two LM constraints remove the scalar shell-count mode, so the propagating axial sector is the tensor odd mode of `f(Sigma) R` on the committed background.\n")
md.append("The canonical axial equation used here is:\n")
md.append("```text\n")
md.append("d_t^2 Psi - d_*^2 Psi + V_ax Psi = 0\n")
md.append("d/dr* = sqrt(h k) d/dr\n")
md.append("V_ax = V_geom + 1/2 d_*P + 1/4 P^2\n")
md.append("P = d_* ln f\n")
md.append("V_geom = h [ ell(ell+1)/r^2 - 2(1-k)/r^2 - k'/(2r) - k h'/(2hr) ]\n")
md.append("```\n")
md.append("For Schwarzschild, `f = 1` and `k = h`, so this reduces exactly to the standard Regge-Wheeler potential.\n")
md.append("\n## Schwarzschild calibration\n")
md.append(f"- omega_TD     = {omega_GR.real:.6f} - {abs(omega_GR.imag):.6f} i\n")
md.append(f"- omega_Leaver = {omega_gold.real:.6f} - {abs(omega_gold.imag):.6f} i\n")
md.append(f"- **|Delta omega / omega| = {err_GR * 100:.4f} %**\n")
md.append("\n## STAM rigorous axial result\n")
md.append(f"- omega_STAM = {omega_STAM.real:.6f} - {abs(omega_STAM.imag):.6f} i\n")
md.append(f"- **|Delta omega / omega| vs calibrated Schwarzschild = {dw_over_w * 100:.4f} %**\n")
md.append(f"- Re shift: {re_shift * 100:+.4f} %\n")
md.append(f"- Im shift: {im_shift * 100:+.4f} %\n")
md.append("\n## Method\n")
md.append(f"- GR grid: N = {N_GR}, r* in [{rs_min_GR}, {rs_max_GR}].\n")
md.append(f"- STAM grid: N = {N_STAM}, r* in [{rs_min_STAM}, {rs_max_STAM}], eps_min = {eps_min:.4e}.\n")
md.append("- Same Sommerfeld time-domain pipeline as G91.\n")
md.append(f"- Fit window: t in [{t_fit_start}, {t_fit_end}] M.\n")
md.append("\n## Files\n")
md.append("- [scripts/G92_rigorous_axial_QNM.py](../scripts/G92_rigorous_axial_QNM.py)\n")
md.append("- [plots/G92_rigorous_axial_QNM.png](../plots/G92_rigorous_axial_QNM.png)\n")

out_md = RESULTS / "G92_rigorous_axial_QNM_summary.md"
out_md.write_text("\n".join(md), encoding="utf-8")
print(f"Summary written: {out_md}")
