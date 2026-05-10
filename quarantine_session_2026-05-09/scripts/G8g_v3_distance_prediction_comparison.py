"""
G8g - Model-A V_3 vs LCDM: distance prediction comparison

Two-tables analog of script 38. Computes and plots μ(z) under two
independent frameworks — Model-A V_3 numerical KG (G8c) and LCDM —
without framing either as the reference. Reports the functional shape
of the prediction difference Delta_mu(z) = μ_V3(z) - μ_LCDM(z) as a model
comparison, not as Model-A's "derived adjustment to bridge to LCDM."

What the original script 38 did: characterized "the derived adjustment
needed to bring STAM's curve onto LCDM's curve." Quarantined under the
two-tables framing — Model-A is researched on its own terms; LCDM is a
comparison framework, not the reference Model-A must justify itself
against.

What this script does (reframed):
  1. Computes μ(z) under Model-A V_3 numerical KG.
  2. Computes μ(z) under LCDM.
  3. Plots both predictions side-by-side.
  4. Computes Delta_mu(z) = μ_V3(z) - μ_LCDM(z) as the difference between
     two independent predictions.
  5. Fits simple functional forms to Delta_mu(z) (constant, z-linear,
     log(1+z), quadratic) — characterizes the shape difference.
  6. Plots Pantheon+ residuals under both frameworks for visual
     side-by-side comparison.
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ===================================================================
# Cosmology setup
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
C_KMS = 299792.458
H0_KMS_MPC = 73.04
H0_LCDM_KMS_MPC = 67.4

print("=" * 72)
print("G8g - V_3 numerical KG vs LCDM distance prediction comparison")
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


def mu_predict(z_obs, h2_func, H0_used, n_grid=8000, z_max=None):
    z_obs = np.atleast_1d(z_obs)
    if z_max is None:
        z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0


# ===================================================================
# Compute both predictions on a smooth grid
# ===================================================================
z_grid = np.linspace(0.001, 2.5, 600)
mu_V3 = mu_predict(z_grid, h2_v3_num, H0_KMS_MPC, z_max=3.0)
mu_LCDM = mu_predict(z_grid, h2_LCDM, H0_LCDM_KMS_MPC, z_max=3.0)
Delta_mu = mu_V3 - mu_LCDM

print(f"\nPrediction-difference Delta_mu(z) = mu_V3(z) - mu_LCDM(z):")
for z in [0.05, 0.1, 0.3, 0.5, 1.0, 2.0]:
    val = float(np.interp(z, z_grid, Delta_mu))
    print(f"  z = {z:.2f}: Delta_mu = {val:+.4f} mag")

# ===================================================================
# Functional form fits to Delta_mu(z)
# ===================================================================
print("\n" + "=" * 72)
print("Functional fits to Delta_mu(z) — characterizing the shape difference")
print("=" * 72)


def fit_polyz(z, dmu, deg):
    coefs = np.polyfit(z, dmu, deg)
    pred = np.polyval(coefs, z)
    rms = float(np.sqrt(np.mean((dmu - pred)**2)))
    return coefs, pred, rms


def fit_log1pz(z, dmu):
    X = np.log1p(z)
    coefs = np.polyfit(X, dmu, 1)
    pred = np.polyval(coefs, X)
    rms = float(np.sqrt(np.mean((dmu - pred)**2)))
    return coefs, pred, rms


c0, p0, rms0 = fit_polyz(z_grid, Delta_mu, 0)
c1, p1, rms1 = fit_polyz(z_grid, Delta_mu, 1)
c2, p2, rms2 = fit_polyz(z_grid, Delta_mu, 2)
clog, plog, rmslog = fit_log1pz(z_grid, Delta_mu)

print(f"  constant         : Delta_mu = {c0[0]:+.4f}")
print(f"     RMS residual:  {rms0:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p0)):.4f})")
print(f"  linear in z      : Delta_mu = {c1[1]:+.4f} {c1[0]:+.4f} z")
print(f"     RMS residual:  {rms1:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p1)):.4f})")
print(f"  quadratic in z   : Delta_mu = {c2[2]:+.4f} {c2[1]:+.4f} z {c2[0]:+.4f} z²")
print(f"     RMS residual:  {rms2:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p2)):.4f})")
print(f"  linear ln(1+z)   : Delta_mu = {clog[1]:+.4f} {clog[0]:+.4f} ln(1+z)")
print(f"     RMS residual:  {rmslog:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - plog)):.4f})")

best = min(
    [('constant', rms0), ('linear_z', rms1),
     ('quadratic', rms2), ('linear_ln1pz', rmslog)],
    key=lambda x: x[1])
print(f"\n  Best simple fit  : {best[0]}  (RMS = {best[1]:.4f} mag)")

# Fractional luminosity-distance difference
fractional_dL = 10.0**(Delta_mu / 5.0) - 1.0

# ===================================================================
# Pantheon+ residuals under both frameworks (for visual reference)
# ===================================================================
print("\n" + "=" * 72)
print("Pantheon+ residuals — visual side-by-side")
print("=" * 72)
pantheon = pd.read_csv('data/pantheon.csv')
mP = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
zP = pantheon.loc[mP, 'zCMB'].values
muP_obs = pantheon.loc[mP, 'MU_SH0ES'].values
sigmaP = pantheon.loc[mP, 'MU_SH0ES_ERR_DIAG'].values

mu_V3_at_P = mu_predict(zP, h2_v3_num, H0_KMS_MPC)
mu_LCDM_at_P = mu_predict(zP, h2_LCDM, H0_LCDM_KMS_MPC)
w_P = 1.0 / sigmaP**2
DM_V3 = float(np.sum(w_P * (muP_obs - mu_V3_at_P)) / np.sum(w_P))
DM_LCDM = float(np.sum(w_P * (muP_obs - mu_LCDM_at_P)) / np.sum(w_P))
resid_V3 = muP_obs - mu_V3_at_P - DM_V3
resid_LCDM = muP_obs - mu_LCDM_at_P - DM_LCDM
print(f"  Pantheon+ V_3   : RMS resid = {np.sqrt(np.mean(resid_V3**2)):.4f}")
print(f"  Pantheon+ LCDM  : RMS resid = {np.sqrt(np.mean(resid_LCDM**2)):.4f}")

# ===================================================================
# Plot
# ===================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (0,0): both μ(z) predictions, raw (no offset)
ax = axes[0, 0]
ax.plot(z_grid, mu_V3, lw=2, color='C0', label='Model-A V_3 numerical KG (no offset)')
ax.plot(z_grid, mu_LCDM, lw=2, ls='--', color='k', label='LCDM (no offset)')
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('μ (mag)')
ax.set_title('Distance modulus — independent predictions, raw')
ax.legend(loc='lower right')
ax.grid(alpha=0.3)

# (0,1): Delta_mu(z) = V_3 - LCDM with functional-form fits overlaid
ax = axes[0, 1]
ax.plot(z_grid, Delta_mu, lw=3, color='C3',
        label='Delta_mu = μ_V3 - μ_LCDM (prediction difference)')
ax.plot(z_grid, p1, lw=1.2, ls='--', color='C2',
        label=f'linear: {c1[1]:+.3f} {c1[0]:+.3f} z (RMS={rms1:.4f})')
ax.plot(z_grid, plog, lw=1.2, ls=':', color='C4',
        label=f'ln(1+z): {clog[1]:+.3f} {clog[0]:+.3f} ln(1+z) (RMS={rmslog:.4f})')
ax.axhline(0, color='gray', alpha=0.4)
ax.set_xlabel('redshift z')
ax.set_ylabel('Delta_mu (mag)')
ax.set_title('Shape of the prediction difference between two frameworks')
ax.legend(loc='best', fontsize=9)
ax.grid(alpha=0.3)

# (1,0): Pantheon+ binned residuals under both
ax = axes[1, 0]
nbins = 25
zbins = np.geomspace(zP.min() * 1.01, zP.max(), nbins + 1)
zcent = 0.5 * (zbins[1:] + zbins[:-1])
res_V3_b = []
res_LCDM_b = []
sem_b = []
for i in range(nbins):
    m = (zP >= zbins[i]) & (zP < zbins[i + 1])
    if m.sum() > 0:
        wi = 1.0 / sigmaP[m]**2
        res_V3_b.append(np.sum(wi * resid_V3[m]) / np.sum(wi))
        res_LCDM_b.append(np.sum(wi * resid_LCDM[m]) / np.sum(wi))
        sem_b.append(1.0 / np.sqrt(np.sum(wi)))
    else:
        res_V3_b.append(np.nan); res_LCDM_b.append(np.nan); sem_b.append(np.nan)
res_V3_b = np.array(res_V3_b); res_LCDM_b = np.array(res_LCDM_b); sem_b = np.array(sem_b)

ax.errorbar(zcent, res_V3_b, yerr=sem_b, fmt='s-', color='C0',
            label='V_3 numerical KG residuals', lw=2)
ax.errorbar(zcent, res_LCDM_b, yerr=sem_b, fmt='o--', color='k',
            label='LCDM residuals', lw=2)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('mean residual (mag, per-catalog ΔM marginalized)')
ax.set_title('Pantheon+ binned residuals — both frameworks')
ax.legend()
ax.grid(alpha=0.3)

# (1,1): fractional luminosity-distance difference
ax = axes[1, 1]
ax.plot(z_grid, fractional_dL * 100, lw=2.5, color='C5')
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xlabel('redshift z')
ax.set_ylabel('(d_L_V3 - d_L_LCDM) / d_L_LCDM  [%]')
ax.set_title('Fractional luminosity-distance difference')
ax.grid(alpha=0.3)

plt.suptitle(
    "Two-tables comparison — Model-A V_3 numerical KG and LCDM as "
    "independent distance predictions",
    fontsize=12,
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/G8g')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_lcdm_distance_comparison.png', dpi=140,
            bbox_inches='tight')

# Save Delta_mu(z) table
df = pd.DataFrame({
    'z': z_grid,
    'mu_V3': mu_V3,
    'mu_LCDM': mu_LCDM,
    'Delta_mu_V3_minus_LCDM': Delta_mu,
    'fractional_dL_pct': fractional_dL * 100,
})
df.to_csv(outdir / 'v3_lcdm_distance_difference.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 72)
