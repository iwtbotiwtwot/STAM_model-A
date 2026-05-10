"""
G8d - Continuity scatter test (Sean's prediction)

Sean's prediction (2026-05-09):
  "If continuity occurs at a certain point, measurements beyond will be
  consistent, measurements below will not, this would result in a scatter
  of sorts in the data at the threshold of continuity."

In Model-A V_3 numerical KG, the cosmic energy budget transitions from
matter-dominated to field-coasting near z ~ 0.30-0.35. Above the
transition: matter-dominated regime, cosmology behaves consistently in
one way. Below: field-coasting regime, consistently in another.
Through the transition, the field is rolling fastest, so individual
sight lines sample slightly different effective expansions, producing
enhanced apparent-distance scatter.

Test:
  1. Compute Pantheon+ residuals (mu_obs - mu_pred - DeltaM_marg) under
     Model-A V_3 numerical KG and under LCDM.
  2. Bin in z. For each bin compute variance of *standardized* residuals
     (r_i / sigma_i). Under a perfect-fit null hypothesis with correct
     reported errors, this variance should be ~1 in every bin.
  3. Look for an excess at z ~ 0.30-0.35.

If the excess is present and concentrated at the transition redshift,
Sean's continuity-scatter prediction is supported. If it's flat, no
detection (which doesn't falsify Model-A, just means the predicted
excess is below catalog noise).

Comparison with LCDM is informative: LCDM has no transition mechanism
(constant Lambda), so any z-localized excess scatter that LCDM also
shows is intrinsic catalog structure, not a Model-A-distinctive signal.

Author: STAM Model-A research line, 2026-05-09.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq
from astropy.io import fits

# ===================================================================
# Load V_3 numerical KG cosmology from G8c trajectory
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
ALPHA_OVER_BETA = (A0 / (1.0 - A0))**2
BETA_TILDE_NUM = 0.426474   # G8c numerical-KG calibration

C_KMS = 299792.458
H0_KMS_MPC = 73.04
H0_LCDM_KMS_MPC = 67.4

print("=" * 72)
print("G8d - Continuity scatter test (Sean's prediction)")
print("=" * 72)

traj = pd.read_csv('reports/G8c/v3_kg_trajectory.csv')
lna_traj = traj['lna'].values
h2_traj = traj['h2'].values

idx_sort = np.argsort(lna_traj)
lna_sorted = lna_traj[idx_sort]
h2_sorted = h2_traj[idx_sort]


def h2_v3_num(a):
    a = np.atleast_1d(np.asarray(a, dtype=float))
    lna = np.log(a)
    lna_clipped = np.clip(lna, lna_sorted[0], lna_sorted[-1])
    return np.interp(lna_clipped, lna_sorted, h2_sorted)


def h2_LCDM(a):
    Om = 0.315
    Or = 9.2e-5
    OL = 1.0 - Om - Or
    return Om * a**-3 + Or * a**-4 + OL


# ===================================================================
# Distance modulus
# ===================================================================
def cumtrapz0(y, x):
    dx = np.diff(x)
    midy = 0.5 * (y[1:] + y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy * dx)])


def mu_predict(z_obs, h2_func, H0_used, n_grid=8000):
    z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0


# ===================================================================
# Load Pantheon+ data and STAT+SYS covariance for proper sigma
# ===================================================================
print("Loading Pantheon+ data and STAT+SYS covariance...")
pantheon = pd.read_csv('data/pantheon.csv')
mask = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
mask_idx = np.where(mask.values)[0]
zP = pantheon.loc[mask, 'zCMB'].values
muP = pantheon.loc[mask, 'MU_SH0ES'].values

with open('data/Pantheon+SH0ES_STAT+SYS.cov', 'r') as f:
    N = int(f.readline().strip())
    flat = np.fromfile(f, sep='\n', dtype=np.float64, count=N * N)
C_full = flat.reshape(N, N)
C_full = 0.5 * (C_full + C_full.T)
covP = C_full[np.ix_(mask_idx, mask_idx)]
sigma_P = np.sqrt(np.diag(covP))   # diagonal sigma per SN

print(f"  Pantheon+: N = {len(zP)} SNe")
print(f"  z range:   [{zP.min():.4f}, {zP.max():.4f}]")
print(f"  sigma per SN: median = {np.median(sigma_P):.4f}")


# ===================================================================
# Compute residuals under each cosmology with marginalized DeltaM
# ===================================================================
def fit_DeltaM_full_cov(z_obs, mu_obs, cov, h2_func, H0_used):
    mu_model = mu_predict(z_obs, h2_func, H0_used)
    delta = mu_obs - mu_model
    ones = np.ones_like(delta)
    L = np.linalg.cholesky(cov)
    Cinv_delta = np.linalg.solve(L.T, np.linalg.solve(L, delta))
    Cinv_ones = np.linalg.solve(L.T, np.linalg.solve(L, ones))
    a = float(ones @ Cinv_ones)
    b = float(ones @ Cinv_delta)
    DeltaM = b / a
    resid = delta - DeltaM
    return resid, DeltaM, mu_model


print("\nComputing residuals under V_3 numerical KG and LCDM...")
resid_V3, dM_V3, mu_V3 = fit_DeltaM_full_cov(zP, muP, covP, h2_v3_num, H0_KMS_MPC)
resid_L, dM_L, mu_L = fit_DeltaM_full_cov(zP, muP, covP, h2_LCDM, H0_LCDM_KMS_MPC)
print(f"  V_3 numerical KG : DeltaM = {dM_V3:+.4f}, RMS resid = {np.sqrt(np.mean(resid_V3**2)):.4f} mag")
print(f"  LCDM             : DeltaM = {dM_L:+.4f}, RMS resid = {np.sqrt(np.mean(resid_L**2)):.4f} mag")

# Standardized residuals (r_i / sigma_i)
std_V3 = resid_V3 / sigma_P
std_L = resid_L / sigma_P

# ===================================================================
# Bin in z and compute variance of standardized residuals
# ===================================================================
# Use logarithmic z bins to balance low-z density and high-z reach.
# Test bin counts: ~30 bins gives reasonable resolution at the
# transition while preserving enough SNe per bin for variance estimate.
n_bins = 30
zbins = np.geomspace(zP.min() * 1.001, zP.max() * 1.001, n_bins + 1)
zcent = np.sqrt(zbins[1:] * zbins[:-1])  # geometric mean


def bin_variance(z, std_resid, edges):
    """For each bin, compute variance of standardized residuals.

    Variance = sum(r_i^2) / N where r_i = (resid_i / sigma_i).
    Under perfect-fit null with correct sigma, this should be ~1.
    Returns: bin_center, variance, N_per_bin, chi2_per_dof_per_bin (=variance)
    """
    var = np.full(len(edges) - 1, np.nan)
    Ns = np.zeros(len(edges) - 1, dtype=int)
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1])
        N = m.sum()
        Ns[i] = N
        if N > 5:   # need at least 5 to have a meaningful variance
            var[i] = np.mean(std_resid[m]**2)
    return var, Ns


var_V3, Ns = bin_variance(zP, std_V3, zbins)
var_L, _ = bin_variance(zP, std_L, zbins)

# Expected variance under null: 1.0 in every bin.
# Statistical fluctuation of variance with N points: Chi^2_N / N
# has variance 2/N, so 1-sigma error on the bin variance is sqrt(2/N).
sigma_var = np.where(Ns > 5, np.sqrt(2.0 / np.where(Ns > 0, Ns, 1)),
                     np.nan)

# ===================================================================
# Compute excess significance per bin
# ===================================================================
print("\n" + "=" * 72)
print("Per-bin variance of standardized residuals (r_i / sigma_i)^2")
print("Under perfect-fit null: variance should be ~1.0; >1 indicates excess scatter.")
print("=" * 72)
print(f"{'z_lo':>8} {'z_hi':>8} {'N':>5} "
      f"{'var V_3':>9} {'var LCDM':>9} {'1sig':>8} "
      f"{'V_3 excess sigma':>18}")
for i in range(len(zbins) - 1):
    if Ns[i] >= 5:
        excess_sigma_V3 = (var_V3[i] - 1.0) / sigma_var[i]
        marker = " ***" if (zbins[i] >= 0.20 and zbins[i + 1] <= 0.50) else ""
        print(f"{zbins[i]:>8.4f} {zbins[i+1]:>8.4f} {Ns[i]:>5d} "
              f"{var_V3[i]:>9.4f} {var_L[i]:>9.4f} {sigma_var[i]:>8.4f} "
              f"{excess_sigma_V3:>+15.2f} sigma{marker}")

# ===================================================================
# Highlight the transition region z = 0.20 - 0.50
# ===================================================================
print("\n" + "=" * 72)
print("Aggregate variance in transition region z = 0.20 - 0.50 vs others:")
print("=" * 72)
mask_trans = (zP >= 0.20) & (zP <= 0.50)
mask_lo = zP < 0.20
mask_hi = zP > 0.50

for label, m in [('Below transition (z<0.20)', mask_lo),
                 ('Transition       (0.20<=z<=0.50)', mask_trans),
                 ('Above transition (z>0.50)', mask_hi)]:
    N = m.sum()
    if N > 0:
        var_V3_reg = float(np.mean(std_V3[m]**2))
        var_L_reg = float(np.mean(std_L[m]**2))
        sigma_reg = math.sqrt(2.0 / N)
        print(f"  {label:<40}  N={N:>5}  V_3 var={var_V3_reg:.4f}  "
              f"LCDM var={var_L_reg:.4f}  1sig={sigma_reg:.4f}")

# ===================================================================
# Plot
# ===================================================================
fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

# Top panel: standardized residuals scatter (per SN)
ax = axes[0]
ax.scatter(zP, std_V3, s=4, alpha=0.4, color='C0', label='V_3 numerical KG')
ax.scatter(zP, std_L, s=4, alpha=0.4, color='gray', label='LCDM')
ax.axhline(0, color='k', alpha=0.5)
ax.axhline(1, color='red', ls=':', alpha=0.4)
ax.axhline(-1, color='red', ls=':', alpha=0.4)
ax.axvspan(0.20, 0.50, alpha=0.12, color='red',
           label='transition region')
ax.set_xscale('log')
ax.set_ylim(-6, 6)
ax.set_ylabel('standardized residual (mu_obs - mu_pred - DM) / sigma')
ax.set_title('Pantheon+ standardized residuals per SN')
ax.legend(loc='upper left', fontsize=9)
ax.grid(alpha=0.3)

# Bottom panel: bin-wise variance of standardized residuals
ax = axes[1]
ax.errorbar(zcent, var_V3, yerr=sigma_var, fmt='s-', color='C0',
            label='V_3 numerical KG  (variance per bin)', lw=2, ms=6)
ax.errorbar(zcent, var_L, yerr=sigma_var, fmt='o--', color='gray',
            label='LCDM  (variance per bin)', lw=2, ms=6, alpha=0.7)
ax.axhline(1.0, color='red', ls=':', alpha=0.7,
           label='null expectation (variance = 1)')
ax.axvspan(0.20, 0.50, alpha=0.12, color='red',
           label='transition region (Sean\'s prediction)')
ax.axvline(0.295, ls=':', color='gray', alpha=0.6,
           label='matter-DE crossover z=0.295')
ax.set_xscale('log')
ax.set_xlabel('redshift z (bin center, geometric mean)')
ax.set_ylabel('variance of standardized residuals')
ax.set_title('Per-bin variance of (mu_obs - mu_pred) / sigma — '
             'continuity-scatter prediction')
ax.legend(loc='upper left', fontsize=9)
ax.grid(alpha=0.3)
ax.set_ylim(0, max(1.0, np.nanmax(var_V3) * 1.1, np.nanmax(var_L) * 1.1))

plt.tight_layout()

outdir = Path('reports/G8d')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'continuity_scatter_test.png', dpi=140, bbox_inches='tight')
plt.savefig(outdir / 'continuity_scatter_test.pdf', bbox_inches='tight')

# Save bin-wise CSV
df = pd.DataFrame({
    'z_lo': zbins[:-1],
    'z_hi': zbins[1:],
    'z_center': zcent,
    'N': Ns,
    'var_V3_num_KG': var_V3,
    'var_LCDM': var_L,
    '1sigma_var': sigma_var,
})
df.to_csv(outdir / 'continuity_scatter_per_bin.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 72)
