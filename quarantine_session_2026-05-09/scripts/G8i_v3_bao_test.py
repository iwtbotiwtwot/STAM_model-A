"""
G8i - Model-A V_3 numerical KG vs DESI DR1 BAO

Independent observational test of the V_3 numerical KG cosmology at
intermediate redshift (z = 0.5 - 2.3) using BAO standard-ruler data.
SN distance fits (G8b/G8c) constrain the cosmology at low-to-mid z;
BAO constrains it at intermediate z via a totally different physical
ruler (sound horizon at the drag epoch). A model that fits both
regimes is on much firmer footing than one that fits only SN.

The previously-quarantined script 45 ran this test against the SU
polynomial cosmology (rejected at 22 sigma). The current Model-A
cosmology (V_3 numerical KG, G8c) is a completely different functional
form and has not been tested against BAO.

DESI DR1 BAO observables: D_M(z) / r_d (transverse comoving distance
divided by sound horizon at drag epoch).

For each framework:
  1. Compute D_M(z) under that cosmology at the BAO redshifts.
  2. Fit a single nuisance r_d (sound horizon) such that
     D_M_pred / r_d matches D_M_obs / r_d in chi^2.
  3. Report chi^2 / dof.

Per the two-tables framing, both Model-A V_3 and LCDM are tested
independently against the same data.
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ===================================================================
# Constants and cosmologies
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
C_KMS = 299792.458
H0_KMS_MPC = 73.04
H0_LCDM_KMS_MPC = 67.4

print("=" * 72)
print("G8i - Model-A V_3 numerical KG vs DESI DR1 BAO")
print("=" * 72)

traj = pd.read_csv('reports/G8c/v3_kg_trajectory.csv')
idx = np.argsort(traj['lna'].values)
lna_sorted = traj['lna'].values[idx]
h2_sorted = traj['h2'].values[idx]


def h2_v3_num(a):
    a = np.atleast_1d(np.asarray(a, dtype=float))
    lna = np.log(a)
    lna_clip = np.clip(lna, lna_sorted[0], lna_sorted[-1])
    return np.interp(lna_clip, lna_sorted, h2_sorted)


def h2_LCDM(a):
    Om = 0.315
    Or = 9.2e-5
    OL = 1.0 - Om - Or
    return Om * a**-3 + Or * a**-4 + OL


def cumtrapz0(y, x):
    dx = np.diff(x)
    midy = 0.5 * (y[1:] + y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy * dx)])


def D_M(z_obs, h2_func, H0_used, n_grid=12000):
    """Comoving (transverse) distance D_M(z) in Mpc, flat universe."""
    z_obs = np.atleast_1d(z_obs)
    z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    return np.interp(z_obs, z_grid, DC_grid)


# ===================================================================
# DESI DR1 BAO data (z, D_M/r_d, sigma) — same as script 45
# ===================================================================
# These are the DESI DR1 D_M/r_d measurements.
BAO = [
    (0.510, 13.62, 0.25),   # LRG1
    (0.706, 16.85, 0.32),   # LRG2
    (0.930, 21.71, 0.28),   # LRG3+ELG1
    (1.317, 27.79, 0.69),   # ELG2
    (2.330, 39.71, 0.94),   # Lyman-alpha
]
zB = np.array([b[0] for b in BAO])
DM_obs = np.array([b[1] for b in BAO])
sigma_obs = np.array([b[2] for b in BAO])

print(f"\nDESI DR1 BAO data points: {len(BAO)}")
for z, dm, s in BAO:
    print(f"  z = {z:.3f}  D_M/r_d = {dm:>5.2f} +/- {s:.2f}")


# ===================================================================
# Fit single r_d nuisance for each cosmology
# ===================================================================
def fit_rd_chi2(zB, DM_obs, sigma, h2_func, H0_used):
    """Compute D_M_pred(z) and find r_d that minimizes chi^2.

    Model: D_M_obs = D_M_pred / r_d. With diagonal sigma,
    1/r_d_best = sum(w * DM_obs * DM_pred) / sum(w * DM_pred^2)
    chi^2 = sum(w * (DM_pred / r_d - DM_obs)^2)
    """
    DM_pred = D_M(zB, h2_func, H0_used)
    w = 1.0 / sigma**2
    inv_rd = float(np.sum(w * DM_obs * DM_pred) / np.sum(w * DM_pred**2))
    rd = 1.0 / inv_rd
    resid = DM_pred / rd - DM_obs
    chi2 = float(np.sum(w * resid**2))
    dof = len(zB) - 1
    return {'rd_Mpc': rd, 'chi2': chi2, 'dof': dof,
            'chi2_per_dof': chi2 / dof, 'DM_pred': DM_pred,
            'resid': resid}


print("\n" + "=" * 72)
print("Fits (single nuisance r_d each, no cosmology fitting)")
print("=" * 72)

res_V3 = fit_rd_chi2(zB, DM_obs, sigma_obs, h2_v3_num, H0_KMS_MPC)
res_LCDM = fit_rd_chi2(zB, DM_obs, sigma_obs, h2_LCDM, H0_LCDM_KMS_MPC)

print(f"\nModel-A V_3 numerical KG (H_0 = {H0_KMS_MPC}):")
print(f"  best-fit r_d   = {res_V3['rd_Mpc']:.2f} Mpc")
print(f"  chi^2          = {res_V3['chi2']:.3f}")
print(f"  chi^2 / dof    = {res_V3['chi2_per_dof']:.4f} ({res_V3['dof']} dof)")
for i, z in enumerate(zB):
    pred = res_V3['DM_pred'][i] / res_V3['rd_Mpc']
    print(f"  z = {z:.3f}  obs = {DM_obs[i]:.2f}  pred = {pred:>5.2f}  "
          f"resid = {res_V3['resid'][i]:>+5.3f} ({res_V3['resid'][i] / sigma_obs[i]:>+.2f} sigma)")

print(f"\nLCDM (H_0 = {H0_LCDM_KMS_MPC}):")
print(f"  best-fit r_d   = {res_LCDM['rd_Mpc']:.2f} Mpc")
print(f"  chi^2          = {res_LCDM['chi2']:.3f}")
print(f"  chi^2 / dof    = {res_LCDM['chi2_per_dof']:.4f} ({res_LCDM['dof']} dof)")
for i, z in enumerate(zB):
    pred = res_LCDM['DM_pred'][i] / res_LCDM['rd_Mpc']
    print(f"  z = {z:.3f}  obs = {DM_obs[i]:.2f}  pred = {pred:>5.2f}  "
          f"resid = {res_LCDM['resid'][i]:>+5.3f} ({res_LCDM['resid'][i] / sigma_obs[i]:>+.2f} sigma)")

# Reference: Planck CMB-derived sound horizon r_d ~ 147 Mpc
# Local-distance-ladder-consistent r_d ~ 137 Mpc
# This isn't a constraint we impose, but worth noting
print()
print("Reference: Planck-CMB-inferred r_d ~ 147 Mpc; local-anchored ~ 137 Mpc.")
print("Each framework's fitted r_d is its own best alignment to the data.")

# ===================================================================
# Two-tables comparison
# ===================================================================
print("\n" + "=" * 72)
print("Two-tables comparison (independent fits to DESI DR1 BAO)")
print("=" * 72)
print()
print("Model-A V_3 numerical KG table:")
print(f"  r_d (best fit)    = {res_V3['rd_Mpc']:.2f} Mpc")
print(f"  chi^2 / dof       = {res_V3['chi2_per_dof']:.4f}")
print(f"  largest residual  = {max(np.abs(res_V3['resid'] / sigma_obs)):.2f} sigma "
      f"(at z = {zB[np.argmax(np.abs(res_V3['resid'] / sigma_obs))]:.3f})")
print()
print("LCDM table:")
print(f"  r_d (best fit)    = {res_LCDM['rd_Mpc']:.2f} Mpc")
print(f"  chi^2 / dof       = {res_LCDM['chi2_per_dof']:.4f}")
print(f"  largest residual  = {max(np.abs(res_LCDM['resid'] / sigma_obs)):.2f} sigma "
      f"(at z = {zB[np.argmax(np.abs(res_LCDM['resid'] / sigma_obs))]:.3f})")

# ===================================================================
# Plot
# ===================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# (0): D_M(z)/r_d data + both predictions
ax = axes[0]
z_smooth = np.linspace(0.05, 2.6, 400)
DM_V3_smooth = D_M(z_smooth, h2_v3_num, H0_KMS_MPC) / res_V3['rd_Mpc']
DM_LCDM_smooth = D_M(z_smooth, h2_LCDM, H0_LCDM_KMS_MPC) / res_LCDM['rd_Mpc']

ax.errorbar(zB, DM_obs, yerr=sigma_obs, fmt='ko', ms=8, capsize=4,
            label='DESI DR1 BAO', zorder=5)
ax.plot(z_smooth, DM_V3_smooth, lw=2, color='C0',
        label=f'Model-A V_3 num KG (r_d = {res_V3["rd_Mpc"]:.1f} Mpc, '
              f'chi²/dof = {res_V3["chi2_per_dof"]:.3f})')
ax.plot(z_smooth, DM_LCDM_smooth, lw=2, ls='--', color='C3',
        label=f'LCDM (r_d = {res_LCDM["rd_Mpc"]:.1f} Mpc, '
              f'chi²/dof = {res_LCDM["chi2_per_dof"]:.3f})')
ax.set_xlabel('redshift z')
ax.set_ylabel('D_M(z) / r_d')
ax.set_title('DESI DR1 BAO: data + framework predictions')
ax.legend(loc='upper left', fontsize=9)
ax.grid(alpha=0.3)

# (1): residuals
ax = axes[1]
ax.errorbar(zB, res_V3['resid'], yerr=sigma_obs, fmt='s', color='C0',
            ms=8, capsize=4, label='V_3 numerical KG residuals')
ax.errorbar(zB, res_LCDM['resid'], yerr=sigma_obs, fmt='o', color='C3',
            ms=8, capsize=4, label='LCDM residuals')
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xlabel('redshift z')
ax.set_ylabel('D_M_pred/r_d - D_M_obs')
ax.set_title('BAO residuals (each framework, independent r_d fit)')
ax.legend(loc='best', fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle(
    "Model-A V_3 numerical KG vs DESI DR1 BAO — apples-to-apples r_d fit",
    fontsize=12,
)
plt.tight_layout(rect=[0, 0, 1, 0.94])

outdir = Path('reports/G8i')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_vs_lcdm_bao.png', dpi=140, bbox_inches='tight')

# Save summary
summary = pd.DataFrame([
    {'framework': 'Model-A V_3 numerical KG', 'r_d_Mpc': res_V3['rd_Mpc'],
     'chi2': res_V3['chi2'], 'dof': res_V3['dof'],
     'chi2_per_dof': res_V3['chi2_per_dof']},
    {'framework': 'LCDM', 'r_d_Mpc': res_LCDM['rd_Mpc'],
     'chi2': res_LCDM['chi2'], 'dof': res_LCDM['dof'],
     'chi2_per_dof': res_LCDM['chi2_per_dof']},
])
summary.to_csv(outdir / 'v3_vs_lcdm_bao_summary.csv', index=False)

# Per-point residuals
detail = pd.DataFrame({
    'z': zB,
    'DM_over_rd_obs': DM_obs,
    'sigma': sigma_obs,
    'V3_DM_pred': res_V3['DM_pred'],
    'V3_DM_over_rd_pred': res_V3['DM_pred'] / res_V3['rd_Mpc'],
    'V3_resid': res_V3['resid'],
    'V3_resid_sigma': res_V3['resid'] / sigma_obs,
    'LCDM_DM_pred': res_LCDM['DM_pred'],
    'LCDM_DM_over_rd_pred': res_LCDM['DM_pred'] / res_LCDM['rd_Mpc'],
    'LCDM_resid': res_LCDM['resid'],
    'LCDM_resid_sigma': res_LCDM['resid'] / sigma_obs,
})
detail.to_csv(outdir / 'v3_vs_lcdm_bao_per_point.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 72)
