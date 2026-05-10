"""
G8h - Catalog self-consistency check for Model-A V_3 and LCDM

Two-tables analog of script 39. For each framework independently, fit
the absolute-magnitude offset ΔM per catalog (Pantheon+, Union3, DES)
and report the cross-catalog spread. A framework whose μ(z) shape is
right at all z gets the same ΔM in every catalog up to noise; a
framework whose shape is wrong gets catalog-dependent ΔM because each
catalog's z-distribution probes the model differently. So the cross-
catalog ΔM spread is a model-self-consistency diagnostic that tests
each framework on its own terms.

What the original script 39 did: framed STAM as predicting the catalog
ΔM tension that LCDM produces. Quarantined under the two-tables
framing — Model-A is researched on its own terms; LCDM is a comparison
framework, not the reference Model-A explains.

What this script does (reframed):
  1. Independently fits ΔM per catalog under Model-A V_3 numerical KG.
  2. Independently fits ΔM per catalog under LCDM.
  3. Reports the cross-catalog ΔM spread for each framework.
  4. Both are tested for self-consistency on the same data.
  5. No claim that one framework "explains" the other's spread.

Cosmology: V_3 numerical KG (G8c).
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits

# ===================================================================
# Cosmology
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
C_KMS = 299792.458
H0_KMS_MPC = 73.04
H0_LCDM_KMS_MPC = 67.4

print("=" * 72)
print("G8h - Catalog self-consistency: Model-A V_3 numerical KG vs LCDM")
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


def mu_predict(z_obs, h2_func, H0_used, n_grid=8000):
    z_obs = np.atleast_1d(z_obs)
    z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0


def fit_DM_diag(z, mu, sigma, h2_func, H0_used):
    mu_model = mu_predict(z, h2_func, H0_used)
    w = 1.0 / sigma**2
    DM = float(np.sum(w * (mu - mu_model)) / np.sum(w))
    chi2 = float(np.sum(w * (mu - mu_model - DM)**2))
    dof = len(z) - 1
    return DM, chi2, dof


def fit_DM_full_cov(z, mu, cov, h2_func, H0_used):
    mu_model = mu_predict(z, h2_func, H0_used)
    delta = mu - mu_model
    ones = np.ones_like(delta)
    L = np.linalg.cholesky(cov)
    Cinv_delta = np.linalg.solve(L.T, np.linalg.solve(L, delta))
    Cinv_ones = np.linalg.solve(L.T, np.linalg.solve(L, ones))
    a = float(ones @ Cinv_ones)
    b = float(ones @ Cinv_delta)
    chi2_unmarg = float(delta @ Cinv_delta)
    DM = b / a
    chi2 = chi2_unmarg - b * b / a
    dof = len(z) - 1
    return DM, chi2, dof


def fit_DM_inv_cov(z, mu, inv_cov, h2_func, H0_used):
    mu_model = mu_predict(z, h2_func, H0_used)
    delta = mu - mu_model
    ones = np.ones_like(delta)
    Cinv_delta = inv_cov @ delta
    Cinv_ones = inv_cov @ ones
    a = float(ones @ Cinv_ones)
    b = float(ones @ Cinv_delta)
    chi2_unmarg = float(delta @ Cinv_delta)
    DM = b / a
    chi2 = chi2_unmarg - b * b / a
    dof = len(z) - 1
    return DM, chi2, dof


# ===================================================================
# Load all three catalogs with proper covariances
# ===================================================================
DATA = Path('data')

print("\nLoading Pantheon+ (STAT+SYS covariance) ...")
pantheon = pd.read_csv(DATA / 'pantheon.csv')
mP = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
mask_idx_P = np.where(mP.values)[0]
zP = pantheon.loc[mP, 'zCMB'].values
muP = pantheon.loc[mP, 'MU_SH0ES'].values
with open(DATA / 'Pantheon+SH0ES_STAT+SYS.cov', 'r') as f:
    Nfull = int(f.readline().strip())
    flat = np.fromfile(f, sep='\n', dtype=np.float64, count=Nfull * Nfull)
C_full = flat.reshape(Nfull, Nfull)
C_full = 0.5 * (C_full + C_full.T)
covP = C_full[np.ix_(mask_idx_P, mask_idx_P)]
print(f"  N = {len(zP)}, median z = {float(np.median(zP)):.4f}")

print("Loading Union3 (FITS inverse covariance) ...")
with fits.open(DATA / 'mu_mat_union3_cosmo2_mu.fits') as hdul:
    fits_data = hdul[0].data
zU = fits_data[0, 1:].astype(float)
muU = fits_data[1:, 0].astype(float)
invCovU = fits_data[1:, 1:].astype(float)
print(f"  N = {len(zU)}, median z = {float(np.median(zU)):.4f}")

print("Loading DES (diagonal MUERR) ...")
des = pd.read_csv(DATA / 'des.csv')
mD = (des['zCMB'] > 0.01) & (des['MUERR'] > 0) & np.isfinite(des['MU'])
zD = des.loc[mD, 'zCMB'].values
muD = des.loc[mD, 'MU'].values
sigmaD = des.loc[mD, 'MUERR'].values
print(f"  N = {len(zD)}, median z = {float(np.median(zD)):.4f}")

# ===================================================================
# Fit per-catalog ΔM under each framework, independently
# ===================================================================
print("\n" + "=" * 72)
print("Per-catalog ΔM fits — each framework on its own table")
print("=" * 72)

# Pantheon+: full STAT+SYS covariance
DM_V3_P, chi2_V3_P, dof_P = fit_DM_full_cov(zP, muP, covP, h2_v3_num, H0_KMS_MPC)
DM_L_P,  chi2_L_P,  _    = fit_DM_full_cov(zP, muP, covP, h2_LCDM, H0_LCDM_KMS_MPC)

# Union3: published inverse-covariance Mahalanobis
DM_V3_U, chi2_V3_U, dof_U = fit_DM_inv_cov(zU, muU, invCovU, h2_v3_num, H0_KMS_MPC)
DM_L_U,  chi2_L_U,  _    = fit_DM_inv_cov(zU, muU, invCovU, h2_LCDM, H0_LCDM_KMS_MPC)

# DES: diagonal MUERR
DM_V3_D, chi2_V3_D, dof_D = fit_DM_diag(zD, muD, sigmaD, h2_v3_num, H0_KMS_MPC)
DM_L_D,  chi2_L_D,  _    = fit_DM_diag(zD, muD, sigmaD, h2_LCDM, H0_LCDM_KMS_MPC)

print()
print("Model-A V_3 numerical KG table:")
print(f"  Pantheon+  ΔM = {DM_V3_P:+.4f}  chi²/dof = {chi2_V3_P / dof_P:.4f}")
print(f"  Union3     ΔM = {DM_V3_U:+.4f}  chi²/dof = {chi2_V3_U / dof_U:.4f}")
print(f"  DES        ΔM = {DM_V3_D:+.4f}  chi²/dof = {chi2_V3_D / dof_D:.4f}")
print()
print("LCDM table:")
print(f"  Pantheon+  ΔM = {DM_L_P:+.4f}  chi²/dof = {chi2_L_P / dof_P:.4f}")
print(f"  Union3     ΔM = {DM_L_U:+.4f}  chi²/dof = {chi2_L_U / dof_U:.4f}")
print(f"  DES        ΔM = {DM_L_D:+.4f}  chi²/dof = {chi2_L_D / dof_D:.4f}")

# ===================================================================
# Cross-catalog ΔM spread per framework
# ===================================================================
def spread_stats(label, dms):
    mean_dm = float(np.mean(dms))
    range_dm = float(np.max(dms) - np.min(dms))
    std_dm = float(np.std(dms, ddof=1))
    return {
        'framework': label,
        'mean_DM': mean_dm,
        'range_DM': range_dm,
        'std_DM': std_dm,
    }


spread_V3 = spread_stats('Model-A V_3 numerical KG',
                          [DM_V3_P, DM_V3_U, DM_V3_D])
spread_LCDM = spread_stats('LCDM',
                            [DM_L_P, DM_L_U, DM_L_D])

print("\n" + "=" * 72)
print("Cross-catalog ΔM spread (each framework's own self-consistency)")
print("=" * 72)
print(f"  {'framework':<28}  {'mean ΔM':>10}  {'range ΔM':>10}  {'std ΔM':>10}")
for s in [spread_V3, spread_LCDM]:
    print(f"  {s['framework']:<28}  {s['mean_DM']:>+10.4f}  "
          f"{s['range_DM']:>10.4f}  {s['std_DM']:>10.4f}")

# Per-pair tensions per framework
print("\nPer-pair ΔM differences (framework's own per-catalog tensions):")
print(f"  {'pair':<28}  {'V_3':>10}  {'LCDM':>10}")
for a, b, label, da_V3, db_V3, da_L, db_L in [
    ('Pantheon+', 'Union3', 'Pantheon+ - Union3', DM_V3_P, DM_V3_U, DM_L_P, DM_L_U),
    ('Pantheon+', 'DES', 'Pantheon+ - DES', DM_V3_P, DM_V3_D, DM_L_P, DM_L_D),
    ('Union3', 'DES', 'Union3 - DES', DM_V3_U, DM_V3_D, DM_L_U, DM_L_D),
]:
    print(f"  {label:<28}  {da_V3 - db_V3:>+10.4f}  {da_L - db_L:>+10.4f}")

# ===================================================================
# Plot
# ===================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

cats = ['Pantheon+', 'Union3', 'DES']
DM_V3 = [DM_V3_P, DM_V3_U, DM_V3_D]
DM_L = [DM_L_P, DM_L_U, DM_L_D]
chi_V3 = [chi2_V3_P / dof_P, chi2_V3_U / dof_U, chi2_V3_D / dof_D]
chi_L = [chi2_L_P / dof_P, chi2_L_U / dof_U, chi2_L_D / dof_D]
xs = np.arange(len(cats))

# Panel 0: per-catalog ΔM
ax = axes[0]
w = 0.35
ax.bar(xs - w / 2, DM_V3, w, label='Model-A V_3 numerical KG', color='C0')
ax.bar(xs + w / 2, DM_L, w, label='LCDM', color='k')
ax.axhline(spread_V3['mean_DM'], ls=':', color='C0', alpha=0.6,
           label=f"V_3 mean = {spread_V3['mean_DM']:+.3f}")
ax.axhline(spread_LCDM['mean_DM'], ls=':', color='k', alpha=0.6,
           label=f"LCDM mean = {spread_LCDM['mean_DM']:+.3f}")
ax.set_xticks(xs)
ax.set_xticklabels(cats)
ax.set_ylabel('best-fit ΔM (mag)')
ax.set_title('Per-catalog ΔM under each framework (independent fits)')
ax.legend(loc='best', fontsize=9)
ax.grid(alpha=0.3, axis='y')
for i, (v, l) in enumerate(zip(DM_V3, DM_L)):
    offset = 0.005 * np.sign(v) if v != 0 else 0.005
    ax.text(i - w / 2, v + offset, f'{v:+.3f}', ha='center', fontsize=8)
    offset_l = 0.005 * np.sign(l) if l != 0 else 0.005
    ax.text(i + w / 2, l + offset_l, f'{l:+.3f}', ha='center', fontsize=8)

# Panel 1: per-catalog chi²/dof
ax = axes[1]
ax.bar(xs - w / 2, chi_V3, w, label='Model-A V_3 numerical KG', color='C0')
ax.bar(xs + w / 2, chi_L, w, label='LCDM', color='k')
ax.axhline(1.0, ls=':', color='red', alpha=0.7, label='ideal = 1')
ax.set_xticks(xs)
ax.set_xticklabels(cats)
ax.set_ylabel('chi²/dof')
ax.set_title('Per-catalog goodness of fit (independent fits)')
ax.legend(loc='best', fontsize=9)
ax.grid(alpha=0.3, axis='y')
for i, (v, l) in enumerate(zip(chi_V3, chi_L)):
    ax.text(i - w / 2, v + 0.05, f'{v:.3f}', ha='center', fontsize=8)
    ax.text(i + w / 2, l + 0.05, f'{l:.3f}', ha='center', fontsize=8)

plt.suptitle(
    "Each framework's own per-catalog self-consistency — "
    "no cross-framework explanatory claims",
    fontsize=12,
)
plt.tight_layout(rect=[0, 0, 1, 0.94])

outdir = Path('reports/G8h')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_lcdm_catalog_self_consistency.png', dpi=140,
            bbox_inches='tight')

# Save summary
summary = pd.DataFrame([
    {'framework': 'Model-A V_3', 'catalog': 'Pantheon+',
     'DeltaM': DM_V3_P, 'chi2_per_dof': chi2_V3_P / dof_P, 'N': len(zP)},
    {'framework': 'Model-A V_3', 'catalog': 'Union3',
     'DeltaM': DM_V3_U, 'chi2_per_dof': chi2_V3_U / dof_U, 'N': len(zU)},
    {'framework': 'Model-A V_3', 'catalog': 'DES',
     'DeltaM': DM_V3_D, 'chi2_per_dof': chi2_V3_D / dof_D, 'N': len(zD)},
    {'framework': 'LCDM', 'catalog': 'Pantheon+',
     'DeltaM': DM_L_P, 'chi2_per_dof': chi2_L_P / dof_P, 'N': len(zP)},
    {'framework': 'LCDM', 'catalog': 'Union3',
     'DeltaM': DM_L_U, 'chi2_per_dof': chi2_L_U / dof_U, 'N': len(zU)},
    {'framework': 'LCDM', 'catalog': 'DES',
     'DeltaM': DM_L_D, 'chi2_per_dof': chi2_L_D / dof_D, 'N': len(zD)},
])
summary.to_csv(outdir / 'v3_lcdm_catalog_self_consistency.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 72)
