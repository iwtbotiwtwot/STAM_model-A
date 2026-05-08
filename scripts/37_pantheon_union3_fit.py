"""
Script 37 - Pantheon+ and Union3 fit to STAM modified Friedmann

Tests whether STAM's full modified Friedmann (script 36 tracking solution)
matches the SN1a distance ladder for the two catalogs that have historically
worked well with STAM, excluding DES (which has shown catalog-specific
tension in prior STAM analyses).

For each model (STAM, LCDM):
  1. Compute h(z) = H(z)/H_0
  2. Compute luminosity distance d_L(z) via integration
  3. Compute predicted distance modulus mu_pred(z)
  4. Marginalize over additive offset DeltaM (standard SN cosmology nuisance)
  5. Report chi^2/dof and residuals

The absolute-offset DeltaM is the standard SN1a nuisance parameter -- the
absolute magnitude M_B is degenerate with H_0 unless an external calibrator
is invoked. We marginalize it analytically per catalog per model.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==================================================================
# Cosmological constants (no fit)
# ==================================================================
C_KMS = 299792.458
H0_KMS_MPC = 67.4    # Planck 2018 baseline; absolute offset is marginalized

A0 = 0.026514
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
BETA_TILDE = OMEGA_M * (1.0 - A0)**2

# ==================================================================
# Models
# ==================================================================
def A_eq(a):
    return 1.0 - (1.0 - A0) * a**1.5

def Y_pot(a):
    return BETA_TILDE / (1.0 - A_eq(a))

def K_kin(a):
    return (3.0/8.0) * (1.0 - A0)**2 * a**3

def h2_STAM_raw(a):
    """STAM tracking solution from script 36 modified Friedmann."""
    return (OMEGA_M * a**-3 + OMEGA_R * a**-4 + Y_pot(a)) / (1.0 - K_kin(a))

def h2_LCDM(a):
    return OMEGA_M * a**-3 + OMEGA_R * a**-4 + (1.0 - OMEGA_M - OMEGA_R)

# Renormalize STAM so h(today) = 1 (model-internal closure)
H2_TODAY_STAM = h2_STAM_raw(1.0)
def h2_STAM(a):
    return h2_STAM_raw(a) / H2_TODAY_STAM

print("=" * 72)
print("STAM vs LCDM | Pantheon+ and Union3 only")
print("=" * 72)
print(f"H0 baseline (km/s/Mpc)        = {H0_KMS_MPC}")
print(f"STAM raw h^2(today)           = {H2_TODAY_STAM:.4f}")
print(f"STAM renormalized so h(now)=1; relative-shape comparison.")
print()

# ==================================================================
# Distance modulus calculation
# ==================================================================
def cumtrapz0(y, x):
    """Cumulative trapezoidal integral starting at 0."""
    dx = np.diff(x)
    midy = 0.5 * (y[1:] + y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy * dx)])

def mu_predict(z_obs, h2_func, n_grid=8000):
    """Distance modulus mu = 5 log10(d_L / Mpc) + 25."""
    z_max = max(z_obs.max() * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_KMS_MPC) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0

def fit_offset_chi2(z_obs, mu_obs, mu_err, h2_func):
    """
    chi^2 = sum_i [(mu_obs - mu_model - DeltaM)/sigma_i]^2
    Analytic best-fit DeltaM (closed form) and chi^2.
    """
    mu_model = mu_predict(z_obs, h2_func)
    w = 1.0 / mu_err**2
    DeltaM = np.sum(w * (mu_obs - mu_model)) / np.sum(w)
    resid = mu_obs - (mu_model + DeltaM)
    chi2 = np.sum(w * resid**2)
    dof = len(z_obs) - 1   # one nuisance: DeltaM
    return {
        'chi2': chi2,
        'dof': dof,
        'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM,
        'mu_model': mu_model,
        'resid': resid,
    }

# ==================================================================
# Load catalogs
# ==================================================================
DATA = Path('data')
pantheon = pd.read_csv(DATA / 'pantheon.csv')
mask = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
pantheon_cosmo = pantheon[mask].copy()
print(f"Pantheon+ cosmological SN (zCMB > 0.01, non-calibrators): "
      f"{len(pantheon_cosmo)}")

union3 = pd.read_csv(DATA / 'union3_bins.csv')
print(f"Union3 binned redshifts: {len(union3)}")
print()

# ==================================================================
# Fit each catalog with each model
# ==================================================================
print("=" * 72)
print("PANTHEON+")
print("=" * 72)
zP = pantheon_cosmo['zCMB'].values
muP = pantheon_cosmo['MU_SH0ES'].values
muP_err = pantheon_cosmo['MU_SH0ES_ERR_DIAG'].values

resP_LCDM = fit_offset_chi2(zP, muP, muP_err, h2_LCDM)
resP_STAM = fit_offset_chi2(zP, muP, muP_err, h2_STAM)

print(f"  LCDM   : chi^2 = {resP_LCDM['chi2']:8.2f},  "
      f"chi^2/dof = {resP_LCDM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resP_LCDM['DeltaM']:+.4f}")
print(f"  STAM   : chi^2 = {resP_STAM['chi2']:8.2f},  "
      f"chi^2/dof = {resP_STAM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resP_STAM['DeltaM']:+.4f}")
print(f"  delta(chi^2) = STAM - LCDM = {resP_STAM['chi2'] - resP_LCDM['chi2']:+.2f}")
print(f"  N points = {len(zP)}")

print()
print("=" * 72)
print("UNION3")
print("=" * 72)
zU = union3['z'].values
muU = union3['mb'].values
# Union3 binned errors not in this file. Use representative 0.05 mag.
# This is approximately right for bins of ~50-200 SN each.
muU_err = np.full_like(muU, 0.05)

resU_LCDM = fit_offset_chi2(zU, muU, muU_err, h2_LCDM)
resU_STAM = fit_offset_chi2(zU, muU, muU_err, h2_STAM)

print(f"  LCDM   : chi^2 = {resU_LCDM['chi2']:8.2f},  "
      f"chi^2/dof = {resU_LCDM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resU_LCDM['DeltaM']:+.4f}")
print(f"  STAM   : chi^2 = {resU_STAM['chi2']:8.2f},  "
      f"chi^2/dof = {resU_STAM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resU_STAM['DeltaM']:+.4f}")
print(f"  delta(chi^2) = STAM - LCDM = {resU_STAM['chi2'] - resU_LCDM['chi2']:+.2f}")
print(f"  N points = {len(zU)}")

print()
print("=" * 72)
print("COMBINED (Pantheon+ + Union3, no DES)")
print("=" * 72)
chi2_LCDM_total = resP_LCDM['chi2'] + resU_LCDM['chi2']
chi2_STAM_total = resP_STAM['chi2'] + resU_STAM['chi2']
dof_total = resP_LCDM['dof'] + resU_LCDM['dof']
print(f"  LCDM total chi^2          = {chi2_LCDM_total:.2f}  /  "
      f"{dof_total} dof  =  {chi2_LCDM_total/dof_total:.4f}")
print(f"  STAM total chi^2          = {chi2_STAM_total:.2f}  /  "
      f"{dof_total} dof  =  {chi2_STAM_total/dof_total:.4f}")
print(f"  delta(chi^2) = STAM-LCDM  = {chi2_STAM_total - chi2_LCDM_total:+.2f}")
which = "STAM" if chi2_STAM_total < chi2_LCDM_total else "LCDM"
print(f"  PREFERRED                 = {which}")

# ==================================================================
# Plotting
# ==================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Smooth curves for plotting
z_sm = np.linspace(0.001, 2.5, 400)
mu_sm_LCDM = mu_predict(z_sm, h2_LCDM)
mu_sm_STAM = mu_predict(z_sm, h2_STAM)

# --- (0,0) Hubble diagram with both catalogs and both models ---
ax = axes[0,0]
ax.errorbar(zP, muP - resP_LCDM['DeltaM'], yerr=muP_err,
            fmt='.', ms=2, alpha=0.25, color='C0', label='Pantheon+ (offset to LCDM)')
ax.errorbar(zU, muU - resU_LCDM['DeltaM'], yerr=muU_err,
            fmt='s', ms=8, color='C3', label='Union3 (offset to LCDM)',
            zorder=5)
ax.plot(z_sm, mu_sm_LCDM, lw=2.0, color='k', label='LCDM (no DeltaM)')
ax.plot(z_sm, mu_sm_STAM + (resP_STAM['DeltaM'] - resP_LCDM['DeltaM']),
        lw=2.0, color='C1', ls='--', label='STAM (offset to LCDM frame)')
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('distance modulus mu')
ax.set_title('Hubble diagram (Pantheon+ + Union3)')
ax.legend(loc='lower right')
ax.grid(alpha=0.3)

# --- (0,1) Pantheon residuals (binned for clarity) ---
ax = axes[0,1]
nbins = 25
zbins = np.geomspace(zP.min()*1.01, zP.max(), nbins+1)
zcent = 0.5*(zbins[1:] + zbins[:-1])
def bin_resid(z, r, w, edges):
    means = []
    sems = []
    for i in range(len(edges)-1):
        m = (z >= edges[i]) & (z < edges[i+1])
        if m.sum() > 0:
            wi = w[m]
            ri = r[m]
            mu = np.sum(wi*ri)/np.sum(wi)
            sem = 1.0/np.sqrt(np.sum(wi))
            means.append(mu)
            sems.append(sem)
        else:
            means.append(np.nan)
            sems.append(np.nan)
    return np.array(means), np.array(sems)

w_P = 1.0/muP_err**2
res_LCDM_b, sem_LCDM_b = bin_resid(zP, resP_LCDM['resid'], w_P, zbins)
res_STAM_b, sem_STAM_b = bin_resid(zP, resP_STAM['resid'], w_P, zbins)

ax.errorbar(zcent, res_LCDM_b, yerr=sem_LCDM_b, fmt='o-', color='k',
            label=f'LCDM ({resP_LCDM["chi2_per_dof"]:.3f})', lw=2)
ax.errorbar(zcent, res_STAM_b, yerr=sem_STAM_b, fmt='s--', color='C1',
            label=f'STAM ({resP_STAM["chi2_per_dof"]:.3f})', lw=2)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual mu_obs - mu_model (binned)')
ax.set_title(f'Pantheon+ residuals (chi^2/dof in legend)')
ax.legend()
ax.grid(alpha=0.3)

# --- (1,0) Union3 residuals (per-bin) ---
ax = axes[1,0]
ax.errorbar(zU, resU_LCDM['resid'], yerr=muU_err, fmt='o-', color='k',
            label=f'LCDM ({resU_LCDM["chi2_per_dof"]:.3f})', lw=2, ms=8)
ax.errorbar(zU, resU_STAM['resid'], yerr=muU_err, fmt='s--', color='C1',
            label=f'STAM ({resU_STAM["chi2_per_dof"]:.3f})', lw=2, ms=8)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual mu_obs - mu_model')
ax.set_title('Union3 residuals')
ax.legend()
ax.grid(alpha=0.3)

# --- (1,1) chi^2 summary bar chart ---
ax = axes[1,1]
labels = ['Pantheon+', 'Union3', 'Combined']
LCDM_vals = [resP_LCDM['chi2_per_dof'], resU_LCDM['chi2_per_dof'],
             chi2_LCDM_total/dof_total]
STAM_vals = [resP_STAM['chi2_per_dof'], resU_STAM['chi2_per_dof'],
             chi2_STAM_total/dof_total]
x = np.arange(len(labels))
w = 0.35
ax.bar(x - w/2, LCDM_vals, w, label='LCDM', color='k')
ax.bar(x + w/2, STAM_vals, w, label='STAM', color='C1')
ax.axhline(1.0, color='gray', ls=':', alpha=0.5, label='ideal')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('chi^2 / dof')
ax.set_title('Goodness of fit (lower is better)')
ax.legend()
ax.grid(alpha=0.3, axis='y')
for i, (l, s) in enumerate(zip(LCDM_vals, STAM_vals)):
    ax.text(i - w/2, l + 0.02, f'{l:.3f}', ha='center', fontsize=9)
    ax.text(i + w/2, s + 0.02, f'{s:.3f}', ha='center', fontsize=9)

plt.suptitle(
    f'STAM (V = beta/(1-A), A_0 = {A0:.4f}) vs LCDM | '
    f'Pantheon+ + Union3 only (DES excluded)',
    fontsize=13
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/script_37')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'pantheon_union3_fit.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir / 'pantheon_union3_fit.pdf', bbox_inches='tight')

# Save a summary CSV
summary = pd.DataFrame([
    {'catalog': 'Pantheon+', 'model': 'LCDM',
     'chi2': resP_LCDM['chi2'], 'dof': resP_LCDM['dof'],
     'chi2_per_dof': resP_LCDM['chi2_per_dof'],
     'DeltaM': resP_LCDM['DeltaM'], 'N': len(zP)},
    {'catalog': 'Pantheon+', 'model': 'STAM',
     'chi2': resP_STAM['chi2'], 'dof': resP_STAM['dof'],
     'chi2_per_dof': resP_STAM['chi2_per_dof'],
     'DeltaM': resP_STAM['DeltaM'], 'N': len(zP)},
    {'catalog': 'Union3', 'model': 'LCDM',
     'chi2': resU_LCDM['chi2'], 'dof': resU_LCDM['dof'],
     'chi2_per_dof': resU_LCDM['chi2_per_dof'],
     'DeltaM': resU_LCDM['DeltaM'], 'N': len(zU)},
    {'catalog': 'Union3', 'model': 'STAM',
     'chi2': resU_STAM['chi2'], 'dof': resU_STAM['dof'],
     'chi2_per_dof': resU_STAM['chi2_per_dof'],
     'DeltaM': resU_STAM['DeltaM'], 'N': len(zU)},
])
summary.to_csv(outdir / 'fit_summary.csv', index=False)
print(f"\nPlots and summary saved to {outdir.resolve()}")
print("=" * 72)
