"""Compute w_eff(z) for V_3 under numerical KG dynamics, comparing to the
tracking proxy. This pins down where Model-A's actual matter-to-field
transition peak lands (Sean's z=0.35 intuition vs the LCDM-style z=0.295
density crossover).
"""
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
ALPHA_OVER_BETA = (A0 / (1.0 - A0))**2
BETA_TILDE_NUM = 0.426474   # from G8c numerical-KG closure


def Y_pot(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return alpha_tilde / A + beta_tilde / (1.0 - A)


def Yp(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return -alpha_tilde / A**2 + beta_tilde / (1.0 - A)**2


def Ypp(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return 2.0 * alpha_tilde / A**3 + 2.0 * beta_tilde / (1.0 - A)**3


# Load the numerical KG trajectory from G8c output
traj = pd.read_csv('reports/G8c/v3_kg_trajectory.csv')
lna = traj['lna'].values
A = traj['A'].values
omega = traj['omega'].values
h = traj['h'].values
a = np.exp(lna)
z = 1.0 / a - 1.0

# rho_DE = Y_pot(A) + omega^2 / 6 (potential + kinetic of the field)
rho_DE_num = Y_pot(A, BETA_TILDE_NUM) + omega**2 / 6.0

# w_eff = -1 - (1/3) d ln rho_DE / d ln a
ln_a = lna
ln_rho_DE = np.log(np.abs(rho_DE_num))   # |abs| in case of zero crossings
w_eff_num = -1.0 - (1.0 / 3.0) * np.gradient(ln_rho_DE, ln_a)

# === Tracking proxy for comparison ===
def A_eq_at(a_val, beta_tilde):
    rho_m = OMEGA_M * a_val**-3
    f = lambda Av: Yp(Av, beta_tilde) - rho_m
    return brentq(f, A0 + 1e-12, 1.0 - 1e-9)


# Tracking with original tracking-calibrated beta
BETA_TILDE_TRK = 0.647255
a_trk = np.linspace(np.exp(lna[0]), 1.0, 600)
A_eq_trk = np.array([A_eq_at(av, BETA_TILDE_TRK) for av in a_trk])
rho_m_trk = OMEGA_M * a_trk**-3
Y_val_trk = Y_pot(A_eq_trk, BETA_TILDE_TRK)
Ypp_val_trk = Ypp(A_eq_trk, BETA_TILDE_TRK)
K_trk = (3.0 / 2.0) * rho_m_trk**2 / Ypp_val_trk**2
h2_trk = (rho_m_trk + OMEGA_R * a_trk**-4 + Y_val_trk) / (1.0 - K_trk)
omega_trk = -3.0 * np.sqrt(h2_trk) * Yp(A_eq_trk, BETA_TILDE_TRK) / Ypp_val_trk
rho_DE_trk = Y_val_trk + omega_trk**2 / 6.0
ln_a_trk = np.log(a_trk)
ln_rho_DE_trk = np.log(np.abs(rho_DE_trk))
w_eff_trk = -1.0 - (1.0 / 3.0) * np.gradient(ln_rho_DE_trk, ln_a_trk)
z_trk = 1.0 / a_trk - 1.0

# Find peak of w_eff in numerical-KG case for the z range of interest
mask = (z >= 0.05) & (z <= 1.0)
peak_idx = np.argmax(w_eff_num[mask])
z_peak_num = z[mask][peak_idx]
w_peak_num = w_eff_num[mask][peak_idx]

# Find peak in tracking
mask_trk = (z_trk >= 0.05) & (z_trk <= 1.0)
peak_idx_trk = np.argmax(w_eff_trk[mask_trk])
z_peak_trk = z_trk[mask_trk][peak_idx_trk]
w_peak_trk = w_eff_trk[mask_trk][peak_idx_trk]

print("V_3 w_eff(z) — numerical KG vs tracking proxy")
print("=" * 72)
print(f"Numerical-KG cosmology (G8c, beta_tilde = {BETA_TILDE_NUM:.4f}):")
print(f"  Peak in w_eff(z) over z=[0.05,1.0]: z = {z_peak_num:.4f}, w = {w_peak_num:+.4f}")
print(f"  w_eff at z=0:    {w_eff_num[np.argmin(np.abs(z))]:+.4f}")
print(f"  w_eff at z=0.30: {np.interp(0.30, z[::-1], w_eff_num[::-1]):+.4f}")
print(f"  w_eff at z=0.35: {np.interp(0.35, z[::-1], w_eff_num[::-1]):+.4f}")
print(f"  w_eff at z=0.50: {np.interp(0.50, z[::-1], w_eff_num[::-1]):+.4f}")
print(f"  w_eff at z=1.00: {np.interp(1.00, z[::-1], w_eff_num[::-1]):+.4f}")
print()
print(f"Tracking proxy cosmology (beta_tilde = {BETA_TILDE_TRK:.4f}):")
print(f"  Peak in w_eff(z) over z=[0.05,1.0]: z = {z_peak_trk:.4f}, w = {w_peak_trk:+.4f}")

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(z_trk, w_eff_trk, lw=2, color='tab:purple', alpha=0.6, ls='--',
        label=f'V_3 tracking proxy (peak at z={z_peak_trk:.3f}, w={w_peak_trk:+.3f})')
ax.plot(z, w_eff_num, lw=2.5, color='C0',
        label=f'V_3 numerical KG (peak at z={z_peak_num:.3f}, w={w_peak_num:+.3f})')
ax.axhline(-1, ls=':', color='k', alpha=0.5, label='w = -1 reference')
ax.axvline(0.295, ls=':', color='gray', alpha=0.6, label='LCDM matter-DE crossover (z=0.295)')
ax.axvspan(0.30, 0.35, alpha=0.15, color='red', label='Sean\'s 0.30-0.35 region')
ax.set_xlim(0, 3)
ax.set_xlabel('redshift z')
ax.set_ylabel('w_eff(z)')
ax.set_title('Model-A V_3 effective dark-energy equation of state\nnumerical KG vs tracking proxy')
ax.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3)

outdir = Path('reports/G8c')
outdir.mkdir(parents=True, exist_ok=True)
plt.tight_layout()
plt.savefig(outdir / 'v3_w_eff_num_vs_tracking.png', dpi=140, bbox_inches='tight')
print(f"\nSaved: {outdir / 'v3_w_eff_num_vs_tracking.png'}")
