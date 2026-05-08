"""
Script 42 - "Einstein was right, we calculated wrong" - the math

Sean's concept (2026-05-07):
  Einstein's matter-only cosmology is correct. The supernova candle data
  is real. We just compute distances WRONG when we assume LCDM. The same
  m(z) data should yield Einstein-de Sitter distances if we propagate
  light correctly through STAM's A field.

Setup:
  d_L_observed(z)  = what SN data measures
  d_L_FRW_EdS(z)   = Einstein matter-only FRW distance:
                     (1+z) * (2c/H0) * [1 - 1/sqrt(1+z)]
  F(z)             = d_L_observed(z) / d_L_FRW_EdS(z)
                   = the multiplicative factor STAM needs to explain

If STAM's photon-A propagation in matter-only cosmology produces
exactly F(z), the framework is self-consistent. We don't need dark
energy in cosmology; we need an A-field optical-effect formula.

This script:
  1. Compute d_L_FRW_EdS(z) and d_L_LCDM(z) (what data fits).
  2. Compute F(z) = d_L_LCDM/d_L_FRW_EdS.
  3. Try to fit F(z) to candidate physical forms:
        ln F(z) = alpha * D_C(z)/D_H              (Beer's law, optical depth ~ path)
        F(z) = (1+z)^p                            (power law in (1+z))
        F(z) = exp[alpha * (1 - 1/sqrt(1+z))]     (alpha-times-EdS-distance form)
  4. Identify the cleanest functional form -- this is the F(z) that STAM's
     photon-A physics has to produce.
  5. State what physical coupling strength is required.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

C_KMS = 299792.458
H0_KMS_MPC = 67.4
D_H_Mpc = C_KMS / H0_KMS_MPC                      # Hubble distance ~4448 Mpc
D_H_Mly = D_H_Mpc * 3.262                         # in Mly

OMEGA_M = 0.315
OMEGA_R = 9.2e-5

def cumtrapz0(y, x):
    dx = np.diff(x); midy = 0.5*(y[1:]+y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy*dx)])

def d_L_LCDM(z):
    z_grid = np.linspace(0.0, z.max()*1.05, 10000)
    a = 1.0/(1.0+z_grid)
    h2 = OMEGA_M*a**-3 + OMEGA_R*a**-4 + (1.0-OMEGA_M-OMEGA_R)
    DC = D_H_Mpc * cumtrapz0(1.0/np.sqrt(h2), z_grid)
    return (1.0 + z) * np.interp(z, z_grid, DC)

def d_L_EdS(z):
    """Einstein-de Sitter (flat, Omega_m = 1, no DE). Closed-form."""
    return (2.0 * D_H_Mpc) * (1.0 + z) * (1.0 - 1.0/np.sqrt(1.0 + z))

def DC_EdS(z):
    return (2.0 * D_H_Mpc) * (1.0 - 1.0/np.sqrt(1.0 + z))

# ==================================================================
# Compute F(z) = d_L_LCDM(z) / d_L_EdS(z)
# This is the multiplicative factor STAM needs to derive.
# ==================================================================
z = np.linspace(0.001, 3.0, 600)
dL_obs   = d_L_LCDM(z)        # what data fits
dL_EdS   = d_L_EdS(z)         # Einstein matter-only
F_needed = dL_obs / dL_EdS

# Distance modulus difference (for context):
Delta_mu = 5.0 * np.log10(F_needed)

print("=" * 72)
print("Required STAM photon-A correction factor F(z)")
print("=" * 72)
print(f"{'z':>6} {'d_L_EdS [Mpc]':>14} {'d_L_LCDM [Mpc]':>15} {'F(z)':>8} {'Delta_mu [mag]':>15}")
for zi in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5]:
    i = np.argmin(np.abs(z - zi))
    print(f"{z[i]:>6.2f} {dL_EdS[i]:>14.1f} {dL_obs[i]:>15.1f} "
          f"{F_needed[i]:>8.4f} {Delta_mu[i]:>+15.4f}")

print()

# ==================================================================
# Fit candidate functional forms
# ==================================================================
print("Functional fits to F(z) (excluding z=0 to avoid log(0)):")

# Form A: ln F(z) = alpha * D_C_EdS(z) / D_H
# i.e., Beer's law optical depth proportional to comoving distance
mask = z > 0.05
DC_over_DH = DC_EdS(z) / D_H_Mpc
ln_F = np.log(F_needed)

# Linear fit (no intercept, since F(0) = 1 → ln F(0) = 0):
alpha_A = np.sum(ln_F[mask] * DC_over_DH[mask]) / np.sum(DC_over_DH[mask]**2)
ln_F_A = alpha_A * DC_over_DH
rms_A = np.sqrt(np.mean((ln_F[mask] - ln_F_A[mask])**2))

# Form B: F(z) = (1+z)^p  ->  ln F = p * ln(1+z)
ln1pz = np.log(1.0 + z)
p_B = np.sum(ln_F[mask] * ln1pz[mask]) / np.sum(ln1pz[mask]**2)
ln_F_B = p_B * ln1pz
rms_B = np.sqrt(np.mean((ln_F[mask] - ln_F_B[mask])**2))

# Form C: F(z) = exp[alpha * (1 - 1/sqrt(1+z))]
# Note: 1 - 1/sqrt(1+z) = D_C_EdS(z)/(2 D_H), so this is just form A scaled by 2
# but let's report it explicitly because it's a clean closed form.
g_z = 1.0 - 1.0/np.sqrt(1.0 + z)
alpha_C = np.sum(ln_F[mask] * g_z[mask]) / np.sum(g_z[mask]**2)
ln_F_C = alpha_C * g_z
rms_C = np.sqrt(np.mean((ln_F[mask] - ln_F_C[mask])**2))

# Form D: F(z) = exp[alpha * z]  (linear-in-z optical depth)
alpha_D = np.sum(ln_F[mask] * z[mask]) / np.sum(z[mask]**2)
ln_F_D = alpha_D * z
rms_D = np.sqrt(np.mean((ln_F[mask] - ln_F_D[mask])**2))

print(f"  (A) ln F = alpha * D_C_EdS/D_H        :  alpha = {alpha_A:.4f}, RMS residual = {rms_A:.5f}")
print(f"  (B) F = (1+z)^p                       :  p     = {p_B:.4f}, RMS residual = {rms_B:.5f}")
print(f"  (C) F = exp[alpha * (1 - 1/sqrt(1+z))]:  alpha = {alpha_C:.4f}, RMS residual = {rms_C:.5f}")
print(f"  (D) F = exp[alpha * z]                :  alpha = {alpha_D:.4f}, RMS residual = {rms_D:.5f}")

best = min([
    ('A: Beer\'s law in D_C',        rms_A, alpha_A),
    ('B: power law (1+z)^p',         rms_B, p_B),
    ('C: alpha * EdS distance',      rms_C, alpha_C),
    ('D: linear-in-z optical depth', rms_D, alpha_D),
], key=lambda x: x[1])
print(f"\nBest single-parameter fit: {best[0]}  with parameter = {best[2]:.4f}")
print()

# ==================================================================
# What does the best-fit alpha tell us PHYSICALLY?
#
# Form C: ln F = alpha * (1 - 1/sqrt(1+z))
# The factor (1 - 1/sqrt(1+z)) = D_C_EdS(z)/(2 D_H) is HALF the EdS comoving
# distance in Hubble units.
#
# Equivalent Beer's law: optical depth tau(z) = alpha * (1 - 1/sqrt(1+z))
# where alpha is a dimensionless coupling strength.
#
# As z -> infty, tau -> alpha (asymptotic optical depth across all of cosmic history).
# Coupling per unit Hubble distance (form A): alpha_A
# ==================================================================
print("=" * 72)
print("Physical interpretation")
print("=" * 72)
print(f"Best fit (form C):  ln F(z) = {alpha_C:.4f} * [1 - 1/sqrt(1+z)]")
print()
print(f"  Asymptotic optical depth as z -> infty: tau_inf = {alpha_C:.4f}")
print(f"  Equivalently, the cumulative photon-A attenuation summed over")
print(f"  all of cosmic history is e^{alpha_C:.4f} = {np.exp(alpha_C):.4f}")
print()
print(f"Form A interpretation: Beer's law with optical depth per Hubble distance:")
print(f"  alpha_per_D_H = {alpha_A:.4f}")
print(f"  Mean free path = D_H / alpha_per_D_H = {1.0/alpha_A:.3f} * D_H")
print(f"                = {1.0/alpha_A * D_H_Mly/1000:.1f} Gly")

# ==================================================================
# Compare to A_0 = 0.0265 calibrated bridge term:
#   If "naive" coupling is alpha_naive = A_0 * (some O(1) factor),
#   how does alpha_naive compare to alpha_required?
# ==================================================================
A0 = 0.026514
# Naive Shapiro-like coupling: alpha_naive = A_0 (constant fractional dimming
# per Hubble distance, factor 2.17 in magnitudes)
# In our form C: alpha_C should equal "effective A_LoS"
print(f"\nComparison to A_0 = {A0:.4f} (bridge-calibrated cosmic mean):")
print(f"  Required alpha (form C)     = {alpha_C:.4f}")
print(f"  Naive A_0 prediction        = {A0:.4f}")
print(f"  Required / naive ratio      = {alpha_C/A0:.1f}")
print(f"  -> Effective integrated A along line of sight needs to be")
print(f"     ~{alpha_C/A0:.0f}x larger than the bridge-calibrated cosmic mean A_0.")
print(f"     This is the 'structure-amplified A_LoS' hypothesis quantified.")

# ==================================================================
# Plot
# ==================================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# (a) d_L curves
ax = axes[0,0]
ax.plot(z, dL_EdS, lw=2.2, color='C2', label='d_L Einstein matter-only (EdS)')
ax.plot(z, dL_obs, lw=2.2, color='k', label='d_L LCDM (data fits this)')
ax.set_xlabel('redshift z'); ax.set_ylabel('d_L (Mpc)')
ax.set_title('Einstein-EdS vs LCDM luminosity distance')
ax.set_xscale('log'); ax.legend(); ax.grid(alpha=0.3)

# (b) F(z) = ratio
ax = axes[0,1]
ax.plot(z, F_needed, lw=2.5, color='C0',
        label='F(z) = d_L_LCDM / d_L_EdS (needed)')
ax.plot(z, np.exp(ln_F_A), lw=1.5, ls='--', color='C1',
        label=f'(A) Beer\'s law: alpha={alpha_A:.3f} per D_H, RMS={rms_A:.4f}')
ax.plot(z, np.exp(ln_F_B), lw=1.5, ls=':', color='C3',
        label=f'(B) (1+z)^p: p={p_B:.3f}, RMS={rms_B:.4f}')
ax.plot(z, np.exp(ln_F_C), lw=1.5, ls='-.', color='C4',
        label=f'(C) alpha=({alpha_C:.3f})*(1-1/sqrt(1+z)), RMS={rms_C:.4f}')
ax.set_xlabel('redshift z'); ax.set_ylabel('F(z)')
ax.set_title('STAM correction factor needed: F(z) = d_L_obs / d_L_EdS')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# (c) ln F(z) — the optical-depth view
ax = axes[1,0]
ax.plot(z, ln_F, lw=2.5, color='C0', label='ln F(z) needed')
ax.plot(z, ln_F_A, lw=1.5, ls='--', color='C1',
        label=f'(A) tau = {alpha_A:.3f} D_C/D_H')
ax.plot(z, ln_F_C, lw=1.5, ls='-.', color='C4',
        label=f'(C) tau = {alpha_C:.3f}(1 - 1/sqrt(1+z))')
ax.set_xlabel('redshift z'); ax.set_ylabel('ln F (= optical depth)')
ax.set_title('"Optical depth" interpretation of STAM correction')
ax.legend(); ax.grid(alpha=0.3)

# (d) residuals from each fit
ax = axes[1,1]
ax.plot(z[mask], (ln_F - ln_F_A)[mask], lw=1.5, color='C1', label='(A) residual')
ax.plot(z[mask], (ln_F - ln_F_B)[mask], lw=1.5, color='C3', label='(B) residual')
ax.plot(z[mask], (ln_F - ln_F_C)[mask], lw=1.5, color='C4', label='(C) residual')
ax.plot(z[mask], (ln_F - ln_F_D)[mask], lw=1.5, color='C5', label='(D) residual')
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xlabel('redshift z'); ax.set_ylabel('ln F residual')
ax.set_title('Fit residuals (smaller = better functional form)')
ax.legend(); ax.grid(alpha=0.3)

plt.suptitle(
    f"Einstein matter-only + STAM photon-A: required correction F(z) = {F_needed[np.argmin(np.abs(z-1))]:.3f} at z=1",
    fontsize=13
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/script_42')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir/'einstein_was_right_math.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir/'einstein_was_right_math.pdf', bbox_inches='tight')

print(f"\nSaved: {outdir.resolve()}")
print("=" * 72)
