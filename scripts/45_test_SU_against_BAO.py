"""
Script 45 - Test STAM SU cosmology against DESI DR1 BAO

Sean's point: STAM predicts a flatter d_L curve than LCDM. SN data fits
LCDM well because LCDM was tuned to fit it. The honest discriminator is
INDEPENDENT data. BAO is the cleanest cross-check: standard ruler with
fully different physics (sound horizon at drag epoch projected onto
transverse comoving distance).

Setup:
  STAM cosmology:  H(z) = H_0 * (1+z)^(2-A_local)  with A_local = 3/10, H_0 = 73.05
  LCDM:            H(z) = H_0 * sqrt(Omega_m (1+z)^3 + Omega_Lambda)
  EdS:             H(z) = H_0 * (1+z)^(3/2)

For each: compute D_M(z), then fit single nuisance r_d (sound horizon)
to DESI BAO observables D_M/r_d. Compare chi^2.

This isolates SHAPE. r_d is a single overall scale that gets absorbed.
The remaining chi^2 measures whether each cosmology has the right
RELATIVE D_M(z) shape across the BAO redshifts.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

C_KMS = 299792.458
H0_SH0ES = 73.05
H0_PLANCK = 67.4
A_LOCAL = 3.0/10.0

# ==================================================================
# DESI DR1 BAO data (5 unique redshifts, D_M/r_d observable)
# Pulled from project's results/bao_geometric_split file
# ==================================================================
BAO_data = [
    ('LRG_0p510',     0.510, 13.62, 0.25),
    ('LRG_0p706',     0.706, 16.85, 0.32),
    ('LRG_ELG_0p930', 0.930, 21.71, 0.28),
    ('ELG_1p317',     1.317, 27.79, 0.69),
    ('Lya_2p330',     2.330, 39.71, 0.94),
]
df_bao = pd.DataFrame(BAO_data, columns=['sample', 'z', 'DM_over_rd', 'sigma'])
print("=" * 72)
print("DESI DR1 BAO data:")
print("=" * 72)
print(df_bao.to_string(index=False))
print()

# ==================================================================
# Cosmology functions
# ==================================================================
def D_M_STAM(z, H0, A):
    """STAM: H(z) = H_0 (1+z)^(2-A), so
       D_M = (c/H_0)/(1-A) * [1 - (1+z)^(-(1-A))]"""
    return (C_KMS/H0)/(1.0 - A) * (1.0 - (1.0 + z)**(-(1.0 - A)))

def D_M_EdS(z, H0):
    return (2.0*C_KMS/H0) * (1.0 - 1.0/np.sqrt(1.0 + z))

def D_M_LCDM(z_arr, H0, Om=0.315):
    """LCDM with flat universe."""
    if np.isscalar(z_arr):
        z_arr = np.array([z_arr])
    OL = 1.0 - Om - 9.2e-5
    z_max = z_arr.max()*1.05
    z_grid = np.linspace(0, z_max, 8000)
    h2 = Om*(1+z_grid)**3 + 9.2e-5*(1+z_grid)**4 + OL
    inv_h = 1.0/np.sqrt(h2)
    DC_grid = (C_KMS/H0) * np.concatenate([[0.0],
              np.cumsum(0.5*(inv_h[1:]+inv_h[:-1])*np.diff(z_grid))])
    return np.interp(z_arr, z_grid, DC_grid)

# ==================================================================
# Predictions for D_M(z) at BAO redshifts
# ==================================================================
z_BAO = df_bao['z'].values
DM_obs = df_bao['DM_over_rd'].values
sigma  = df_bao['sigma'].values

DM_STAM_h73 = D_M_STAM(z_BAO, H0_SH0ES, A_LOCAL)
DM_STAM_h67 = D_M_STAM(z_BAO, H0_PLANCK, A_LOCAL)
DM_EdS_h73  = D_M_EdS(z_BAO, H0_SH0ES)
DM_EdS_h67  = D_M_EdS(z_BAO, H0_PLANCK)
DM_LCDM_h73 = D_M_LCDM(z_BAO, H0_SH0ES)
DM_LCDM_h67 = D_M_LCDM(z_BAO, H0_PLANCK)

# ==================================================================
# Best-fit r_d (sound horizon) per model: marginalize as nuisance
# Observable: DM_obs = D_M(z) / r_d
# So r_d_best = sum(w_i * D_M_pred_i / DM_obs_i) / sum(w_i / DM_obs_i^2)
# Wait, the right least-squares for r_d:
#   chi^2 = sum [(DM_pred/r_d - DM_obs)/sigma]^2
#   d/d(1/r_d) = 0  -->  1/r_d = sum(w_i*DM_obs_i*DM_pred_i)/sum(w_i*DM_pred_i^2)
# ==================================================================
def fit_rd(DM_pred, DM_obs, sigma):
    w = 1.0/sigma**2
    inv_rd = np.sum(w * DM_obs * DM_pred) / np.sum(w * DM_pred**2)
    rd = 1.0/inv_rd
    pred_obs = DM_pred / rd
    chi2 = np.sum(w * (pred_obs - DM_obs)**2)
    return rd, chi2, pred_obs

models = {}
for name, DM in [
    ('STAM (H0=73.05)', DM_STAM_h73),
    ('STAM (H0=67.4) ', DM_STAM_h67),
    ('EdS  (H0=73.05)', DM_EdS_h73),
    ('EdS  (H0=67.4) ', DM_EdS_h67),
    ('LCDM (H0=73.05)', DM_LCDM_h73),
    ('LCDM (H0=67.4) ', DM_LCDM_h67),
]:
    rd, chi2, pred_obs = fit_rd(DM, DM_obs, sigma)
    models[name] = (rd, chi2, pred_obs, DM)

print("=" * 72)
print("BAO fit: DM_obs = D_M_pred(z) / r_d, marginalize r_d as nuisance")
print("=" * 72)
print(f"{'Model':<22} {'r_d (Mpc)':>10} {'chi^2':>10} {'chi^2/dof':>11}")
for name, (rd, chi2, _, _) in models.items():
    dof = len(z_BAO) - 1   # 1 nuisance: r_d
    print(f"{name:<22} {rd:>10.2f} {chi2:>10.4f} {chi2/dof:>11.4f}")
print()
print("Notes:")
print("  - Standard r_d (Planck LCDM) ~ 147 Mpc")
print("  - dof = 5 - 1 = 4 after r_d marginalization")
print()

# ==================================================================
# Detailed per-point comparison for STAM and LCDM
# ==================================================================
rd_STAM, chi2_STAM, pred_STAM, _ = models['STAM (H0=73.05)']
rd_LCDM, chi2_LCDM, pred_LCDM, _ = models['LCDM (H0=67.4) ']

print("=" * 72)
print("Per-point residuals (sigma units):")
print("=" * 72)
print(f"{'z':>7} {'D_M/r_d obs':>12} {'sigma':>7} "
      f"{'STAM pred':>10} {'STAM/sig':>10} {'LCDM pred':>10} {'LCDM/sig':>10}")
for i, zi in enumerate(z_BAO):
    sP = (pred_STAM[i] - DM_obs[i])/sigma[i]
    sL = (pred_LCDM[i] - DM_obs[i])/sigma[i]
    print(f"{zi:>7.3f} {DM_obs[i]:>12.3f} {sigma[i]:>7.3f} "
          f"{pred_STAM[i]:>10.3f} {sP:>+10.2f} {pred_LCDM[i]:>10.3f} {sL:>+10.2f}")
print()

# ==================================================================
# Plot
# ==================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# (a) data + predictions
ax = axes[0]
ax.errorbar(z_BAO, DM_obs, yerr=sigma, fmt='o', color='k', ms=8, label='DESI DR1 BAO')
zsmooth = np.linspace(0.4, 2.5, 200)
ax.plot(zsmooth, D_M_STAM(zsmooth, H0_SH0ES, A_LOCAL)/rd_STAM,
        lw=2, color='C1', label=f"STAM A=3/10, H0=73 (chi^2/dof = {chi2_STAM/4:.2f})")
ax.plot(zsmooth, D_M_LCDM(zsmooth, H0_PLANCK)/rd_LCDM,
        lw=2, color='k', ls='--', label=f"LCDM Planck (chi^2/dof = {chi2_LCDM/4:.2f})")
rd_EdS, chi2_EdS_h73, _, _ = models['EdS  (H0=73.05)']
ax.plot(zsmooth, D_M_EdS(zsmooth, H0_SH0ES)/rd_EdS,
        lw=1.5, color='C2', ls=':', label=f"EdS H0=73 (chi^2/dof = {chi2_EdS_h73/4:.2f})")
ax.set_xlabel('z')
ax.set_ylabel('D_M / r_d')
ax.set_title('DESI DR1 BAO: D_M/r_d vs z')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (b) residuals (sigma units)
ax = axes[1]
ax.bar(np.arange(len(z_BAO))-0.2, [(pred_STAM[i]-DM_obs[i])/sigma[i] for i in range(len(z_BAO))],
       width=0.4, label='STAM', color='C1')
ax.bar(np.arange(len(z_BAO))+0.2, [(pred_LCDM[i]-DM_obs[i])/sigma[i] for i in range(len(z_BAO))],
       width=0.4, label='LCDM', color='k')
ax.axhline(0, color='gray', alpha=0.5)
ax.axhline(1, ls=':', color='gray', alpha=0.4)
ax.axhline(-1, ls=':', color='gray', alpha=0.4)
ax.set_xticks(range(len(z_BAO)))
ax.set_xticklabels([f'{z:.2f}' for z in z_BAO])
ax.set_xlabel('z (BAO redshift bin)')
ax.set_ylabel('residual (sigma)')
ax.set_title('Per-point BAO residuals')
ax.legend()
ax.grid(alpha=0.3, axis='y')

plt.suptitle("DESI DR1 BAO: STAM-SU vs LCDM (fully INDEPENDENT of SN data)",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])

outdir = Path('reports/script_45')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir/'BAO_test.png', dpi=130, bbox_inches='tight')

# ==================================================================
# Verdict
# ==================================================================
print("=" * 72)
print("VERDICT")
print("=" * 72)
ranked = sorted(models.items(), key=lambda x: x[1][1])
print("Models ranked by BAO chi^2:")
for name, (rd, chi2, _, _) in ranked:
    print(f"  {name}: r_d={rd:.1f} Mpc, chi^2/dof={chi2/4:.3f}")
print()
print("Standard Planck-LCDM r_d ~ 147 Mpc. Models giving very different r_d")
print("are inconsistent with BBN+CMB sound-horizon physics.")
print()
print(f"Saved: {outdir.resolve()}")
print("=" * 72)
