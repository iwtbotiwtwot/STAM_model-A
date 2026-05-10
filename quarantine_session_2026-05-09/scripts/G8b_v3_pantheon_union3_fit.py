"""
G8b - Model-A V_3 fit to Pantheon+ and Union3 supernovae

V_3 retrofit of script 37. The cosmology is the V_3 modified Friedmann
in the tracking-ansatz form (same phenomenological choice as V_1 in
script 37), so the V_1-vs-V_3 SN-distance comparison is on equal footing.

Inputs (no fitting at the cosmology level):
    A_0 = 1/(12 pi)                       G7 commitment
    alpha/beta = [A_0/(1-A_0)]^2          V_3 structural minimum at A_0
    Omega_m = 0.315, Omega_r = 9.2e-5     Model-A internal content
    beta_tilde calibrated by h^2(today) = 1 closure (G8a)

Per supernova catalog the only nuisance parameter is the absolute-magnitude
offset DeltaM, marginalized analytically as standard SN cosmology practice.

The other framework (LCDM) is computed and reported on a separate row for
side-by-side independent fit comparison. Per the two-tables framing,
Model-A is on its own table; LCDM is on its table; both are independent
predictions tested against the same data, neither is the reference for
the other.

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
# V_3 cosmology (matches G8a)
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
ALPHA_OVER_BETA = (A0 / (1.0 - A0))**2
H0_KMS_MPC = 73.04                    # Model-A's H_0
H0_LCDM_KMS_MPC = 67.4                # other framework's H_0; absolute-M is marginalized so this is plot label only
C_KMS = 299792.458


def Y_pot(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return alpha_tilde / A + beta_tilde / (1.0 - A)


def Yp(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return -alpha_tilde / A**2 + beta_tilde / (1.0 - A)**2


def Ypp(A, beta_tilde):
    alpha_tilde = beta_tilde * ALPHA_OVER_BETA
    return 2.0 * alpha_tilde / A**3 + 2.0 * beta_tilde / (1.0 - A)**3


def A_eq_at(a, beta_tilde):
    rho_m = OMEGA_M * a**-3
    f = lambda A: Yp(A, beta_tilde) - rho_m
    return brentq(f, A0 + 1e-12, 1.0 - 1e-9)


def h2_v3_tracking(a, beta_tilde):
    A_eq = A_eq_at(a, beta_tilde)
    rho_m = OMEGA_M * a**-3
    rho_r = OMEGA_R * a**-4
    Y_val = Y_pot(A_eq, beta_tilde)
    Ypp_val = Ypp(A_eq, beta_tilde)
    K = (3.0 / 2.0) * rho_m**2 / Ypp_val**2
    return (rho_m + rho_r + Y_val) / (1.0 - K)


# Calibrate beta_tilde by h^2(today) = 1
def find_beta_tilde():
    def f(b):
        return h2_v3_tracking(1.0, b) - 1.0
    return brentq(f, 0.1, 5.0, xtol=1e-10)


BETA_TILDE = find_beta_tilde()

# Vectorize for arrays
def h2_v3_array(a_arr, beta_tilde=BETA_TILDE):
    return np.array([h2_v3_tracking(a, beta_tilde) for a in np.atleast_1d(a_arr)])


def h2_LCDM(a):
    """Other framework, LCDM at H_0 = 67.4, Omega_m = 0.315, Omega_r = 9.2e-5."""
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
    """Distance modulus mu = 5 log10(d_L / Mpc) + 25.

    h2_func(a) is dimensionless; H0_used sets the scale.
    """
    z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0


def fit_offset_chi2(z_obs, mu_obs, mu_err, h2_func, H0_used):
    """Diagonal chi^2 with analytic best-fit DeltaM (standard SN nuisance).

    Used for Pantheon+ where we have per-SN diagonal errors.
    """
    mu_model = mu_predict(z_obs, h2_func, H0_used)
    w = 1.0 / mu_err**2
    DeltaM = np.sum(w * (mu_obs - mu_model)) / np.sum(w)
    resid = mu_obs - (mu_model + DeltaM)
    chi2 = np.sum(w * resid**2)
    dof = len(z_obs) - 1
    return {
        'chi2': chi2,
        'dof': dof,
        'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM,
        'mu_model': mu_model,
        'resid': resid,
    }


def fit_offset_chi2_mahalanobis(z_obs, mu_obs, inv_cov, h2_func, H0_used):
    """Full-covariance Mahalanobis chi^2 with analytic offset marginalization.

    chi^2(DM) = (Delta - 1*DM)^T C^{-1} (Delta - 1*DM)
    where Delta = mu_obs - mu_pred and 1 is the all-ones vector.
    Analytic minimum:
      DM_best = (1^T C^{-1} Delta) / (1^T C^{-1} 1)
      chi^2_min = Delta^T C^{-1} Delta - (1^T C^{-1} Delta)^2 / (1^T C^{-1} 1)
    Used for Union3 where the published inverse-covariance has off-diagonal
    bin-bin correlations that diagonal-only forms drop.
    """
    mu_model = mu_predict(z_obs, h2_func, H0_used)
    delta = mu_obs - mu_model
    ones = np.ones_like(delta)
    Cinv_delta = inv_cov @ delta
    Cinv_ones = inv_cov @ ones
    a = float(ones @ Cinv_ones)
    b = float(ones @ Cinv_delta)
    chi2_unmarg = float(delta @ Cinv_delta)
    DeltaM = b / a
    chi2 = chi2_unmarg - b * b / a
    resid = delta - DeltaM
    dof = len(z_obs) - 1
    return {
        'chi2': chi2,
        'dof': dof,
        'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM,
        'mu_model': mu_model,
        'resid': resid,
    }


def load_union3_fits(path='data/mu_mat_union3_cosmo2_mu.fits'):
    """Load Union3 binned distances + inverse covariance from the official FITS.

    Returns: z (22,), mu (22,), inv_cov (22, 22).
    The FITS layout: 23x23 array. data[0, 1:] = z values; data[1:, 0] = mu values;
    data[1:, 1:] = inverse covariance matrix.
    """
    with fits.open(path) as hdul:
        data = hdul[0].data
    z = np.asarray(data[0, 1:], dtype=float)
    mu = np.asarray(data[1:, 0], dtype=float)
    inv_cov = np.asarray(data[1:, 1:], dtype=float)
    return z, mu, inv_cov


def load_pantheon_full_cov(cov_path='data/Pantheon+SH0ES_STAT+SYS.cov'):
    """Load the published Pantheon+ stat+sys covariance matrix.

    Format: first line is N=1701, then N*N floats one per line in row-major
    order. Row ordering matches pantheon.csv exactly.
    Returns: C (1701, 1701) symmetrized to clean text-format float asymmetry.
    """
    with open(cov_path, 'r') as f:
        N = int(f.readline().strip())
        flat = np.fromfile(f, sep='\n', dtype=np.float64, count=N * N)
    C = flat.reshape(N, N)
    return 0.5 * (C + C.T)


def fit_offset_chi2_mahalanobis_full(z_obs, mu_obs, cov, h2_func, H0_used):
    """Mahalanobis chi^2 with analytic offset marginalization, using full
    covariance (not inverse). Solves C x = delta via Cholesky for stability.
    """
    mu_model = mu_predict(z_obs, h2_func, H0_used)
    delta = mu_obs - mu_model
    ones = np.ones_like(delta)
    L = np.linalg.cholesky(cov)
    Cinv_delta = np.linalg.solve(L.T, np.linalg.solve(L, delta))
    Cinv_ones = np.linalg.solve(L.T, np.linalg.solve(L, ones))
    a = float(ones @ Cinv_ones)
    b = float(ones @ Cinv_delta)
    chi2_unmarg = float(delta @ Cinv_delta)
    DeltaM = b / a
    chi2 = chi2_unmarg - b * b / a
    resid = delta - DeltaM
    dof = len(z_obs) - 1
    return {
        'chi2': chi2,
        'dof': dof,
        'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM,
        'mu_model': mu_model,
        'resid': resid,
    }


# ===================================================================
# Header
# ===================================================================
print("=" * 72)
print("Model-A V_3 vs LCDM | Pantheon+ + Union3 SN distance fit (G8b)")
print("=" * 72)
print(f"V_3:     A_0 = 1/(12 pi) = {A0:.8f}")
print(f"         alpha/beta     = {ALPHA_OVER_BETA:.6e}")
print(f"         beta_tilde     = {BETA_TILDE:.6f}  (calibrated by Model-A internal closure)")
print(f"         h^2(today, V_3 tracking) = {h2_v3_tracking(1.0, BETA_TILDE):.6f}")
print(f"         H_0 (Model-A, plot label) = {H0_KMS_MPC} km/s/Mpc")
print(f"LCDM:    Omega_m = 0.315, Omega_L = 0.685, H_0 = {H0_LCDM_KMS_MPC} km/s/Mpc")
print(f"         (absolute magnitude DeltaM marginalized analytically per catalog)")
print()

# ===================================================================
# Load catalogs
# ===================================================================
DATA = Path('data')
pantheon = pd.read_csv(DATA / 'pantheon.csv')
mask = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
mask_idx = np.where(mask.values)[0]
pantheon_cosmo = pantheon[mask].copy()
zP = pantheon_cosmo['zCMB'].values
muP = pantheon_cosmo['MU_SH0ES'].values
muPe = pantheon_cosmo['MU_SH0ES_ERR_DIAG'].values

# Load published STAT+SYS covariance and submatrix to cosmological mask
print(f"Pantheon+ cosmological SN (zCMB > 0.01, non-calibrators): {len(zP)}")
print("Loading Pantheon+ STAT+SYS covariance...")
C_full = load_pantheon_full_cov(DATA / 'Pantheon+SH0ES_STAT+SYS.cov')
print(f"  loaded {C_full.shape}, submatrixing to {len(zP)} cosmological rows")
covP = C_full[np.ix_(mask_idx, mask_idx)]
sigma_diag_P = np.sqrt(np.diag(covP))
print(f"  per-SN sigma (sqrt diag of cov): "
      f"min={sigma_diag_P.min():.4f}, median={float(np.median(sigma_diag_P)):.4f}, "
      f"max={sigma_diag_P.max():.4f}")

zU, muU, invCovU = load_union3_fits(DATA / 'mu_mat_union3_cosmo2_mu.fits')
sigma_diag_U = np.sqrt(np.diag(np.linalg.inv(invCovU)))
print(f"Union3 binned redshifts: {len(zU)}  (proper inverse covariance from "
      f"data/mu_mat_union3_cosmo2_mu.fits)")
print(f"  per-bin sigma (sqrt diag of cov): "
      f"min={sigma_diag_U.min():.3f}, median={float(np.median(sigma_diag_U)):.3f}, "
      f"max={sigma_diag_U.max():.3f}")
# Also load the legacy per-bin CSV (no errors) for plotting only
muUe = sigma_diag_U
print()

# ===================================================================
# Fit each catalog with each model independently
# ===================================================================
print("=" * 72)
print("PANTHEON+ (Mahalanobis chi^2 with published STAT+SYS covariance)")
print("=" * 72)
resP_V3 = fit_offset_chi2_mahalanobis_full(zP, muP, covP, h2_v3_array, H0_KMS_MPC)
resP_LCDM = fit_offset_chi2_mahalanobis_full(zP, muP, covP, h2_LCDM, H0_LCDM_KMS_MPC)

print(f"  Model-A V_3  : chi^2 = {resP_V3['chi2']:8.2f},  "
      f"chi^2/dof = {resP_V3['chi2_per_dof']:.4f},  "
      f"DeltaM = {resP_V3['DeltaM']:+.4f}")
print(f"  LCDM         : chi^2 = {resP_LCDM['chi2']:8.2f},  "
      f"chi^2/dof = {resP_LCDM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resP_LCDM['DeltaM']:+.4f}")
print(f"  delta(chi^2) = V_3 - LCDM = {resP_V3['chi2'] - resP_LCDM['chi2']:+.2f}")
print(f"  N points = {len(zP)}")
print()

print("=" * 72)
print("UNION3 (Mahalanobis chi^2 with proper bin-bin covariance)")
print("=" * 72)
resU_V3 = fit_offset_chi2_mahalanobis(zU, muU, invCovU, h2_v3_array, H0_KMS_MPC)
resU_LCDM = fit_offset_chi2_mahalanobis(zU, muU, invCovU, h2_LCDM, H0_LCDM_KMS_MPC)

print(f"  Model-A V_3  : chi^2 = {resU_V3['chi2']:8.2f},  "
      f"chi^2/dof = {resU_V3['chi2_per_dof']:.4f},  "
      f"DeltaM = {resU_V3['DeltaM']:+.4f}")
print(f"  LCDM         : chi^2 = {resU_LCDM['chi2']:8.2f},  "
      f"chi^2/dof = {resU_LCDM['chi2_per_dof']:.4f},  "
      f"DeltaM = {resU_LCDM['DeltaM']:+.4f}")
print(f"  delta(chi^2) = V_3 - LCDM = {resU_V3['chi2'] - resU_LCDM['chi2']:+.2f}")
print(f"  N points = {len(zU)}")
print()

print("=" * 72)
print("COMBINED (Pantheon+ + Union3)")
print("=" * 72)
chi2_V3_total = resP_V3['chi2'] + resU_V3['chi2']
chi2_LCDM_total = resP_LCDM['chi2'] + resU_LCDM['chi2']
dof_total = resP_V3['dof'] + resU_V3['dof']
print(f"  Model-A V_3 total chi^2  = {chi2_V3_total:.2f}  /  "
      f"{dof_total} dof  =  {chi2_V3_total / dof_total:.4f}")
print(f"  LCDM        total chi^2  = {chi2_LCDM_total:.2f}  /  "
      f"{dof_total} dof  =  {chi2_LCDM_total / dof_total:.4f}")
print(f"  delta(chi^2) = V_3 - LCDM = {chi2_V3_total - chi2_LCDM_total:+.2f}")

print()
print("Two-tables comparison (independent fits to the same data):")
print()
print("  Model-A V_3 table:")
print(f"    Pantheon+ chi^2/dof = {resP_V3['chi2_per_dof']:.4f}")
print(f"    Union3    chi^2/dof = {resU_V3['chi2_per_dof']:.4f}")
print(f"    Combined  chi^2/dof = {chi2_V3_total / dof_total:.4f}")
print()
print("  LCDM table:")
print(f"    Pantheon+ chi^2/dof = {resP_LCDM['chi2_per_dof']:.4f}")
print(f"    Union3    chi^2/dof = {resU_LCDM['chi2_per_dof']:.4f}")
print(f"    Combined  chi^2/dof = {chi2_LCDM_total / dof_total:.4f}")

# ===================================================================
# V_1 reference comparison: re-run script 37's V_1 calibration here
# for direct V_1-vs-V_3 chi^2 comparison
# ===================================================================
print()
print("=" * 72)
print("Reference: V_1 = beta/(1-A) cosmology (matches script 37)")
print("=" * 72)
A0_V1 = 0.026514                              # script 37 used the empirical A_0
BETA_TILDE_V1 = OMEGA_M * (1.0 - A0_V1)**2     # V_1 calibration: V'(A_0) = kappa rho_m_today


def A_eq_V1(a):
    return 1.0 - (1.0 - A0_V1) * a**1.5


def Y_pot_V1(a):
    return BETA_TILDE_V1 / (1.0 - A_eq_V1(a))


def K_kin_V1(a):
    return (3.0 / 8.0) * (1.0 - A0_V1)**2 * a**3


def h2_V1_raw(a):
    return (OMEGA_M * a**-3 + OMEGA_R * a**-4 + Y_pot_V1(a)) / (1.0 - K_kin_V1(a))


h2_V1_today = h2_V1_raw(1.0)


def h2_V1(a):
    return h2_V1_raw(a) / h2_V1_today


resP_V1 = fit_offset_chi2_mahalanobis_full(zP, muP, covP, h2_V1, H0_LCDM_KMS_MPC)
resU_V1 = fit_offset_chi2_mahalanobis(zU, muU, invCovU, h2_V1, H0_LCDM_KMS_MPC)
chi2_V1_total = resP_V1['chi2'] + resU_V1['chi2']
print(f"  V_1 Pantheon+ chi^2/dof = {resP_V1['chi2_per_dof']:.4f}")
print(f"  V_1 Union3    chi^2/dof = {resU_V1['chi2_per_dof']:.4f}")
print(f"  V_1 combined  chi^2/dof = {chi2_V1_total / dof_total:.4f}")
print()
print("Direct V_1 vs V_3 comparison (Model-A internal):")
print(f"  V_1 - V_3 combined chi^2 = {chi2_V1_total - chi2_V3_total:+.2f}")
which = "V_3" if chi2_V3_total < chi2_V1_total else "V_1"
print(f"  preferred (Model-A internal) = {which}")

# ===================================================================
# Plotting
# ===================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

z_sm = np.linspace(0.001, 2.5, 400)
mu_sm_V3 = mu_predict(z_sm, h2_v3_array, H0_KMS_MPC)
mu_sm_LCDM = mu_predict(z_sm, h2_LCDM, H0_LCDM_KMS_MPC)
mu_sm_V1 = mu_predict(z_sm, h2_V1, H0_LCDM_KMS_MPC)

# (0,0) Hubble diagram
ax = axes[0, 0]
ax.errorbar(zP, muP - resP_V3['DeltaM'], yerr=muPe,
            fmt='.', ms=2, alpha=0.25, color='C0', label='Pantheon+ (offset to V_3)')
ax.errorbar(zU, muU - resU_V3['DeltaM'], yerr=muUe,
            fmt='s', ms=8, color='C3', label='Union3 (offset to V_3)', zorder=5)
ax.plot(z_sm, mu_sm_V3, lw=2.0, color='C0', label='Model-A V_3 (no offset)')
ax.plot(z_sm, mu_sm_LCDM + (resP_V3['DeltaM'] - resP_LCDM['DeltaM']),
        lw=2.0, color='k', ls='--', label='LCDM (offset to V_3 frame)')
ax.plot(z_sm, mu_sm_V1 + (resP_V3['DeltaM'] - resP_V1['DeltaM']),
        lw=1.5, color='C2', ls=':', label='V_1 reference (offset to V_3 frame)')
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('distance modulus mu')
ax.set_title('Hubble diagram (Pantheon+ + Union3)')
ax.legend(loc='lower right', fontsize=9)
ax.grid(alpha=0.3)

# (0,1) Pantheon residuals binned
ax = axes[0, 1]
nbins = 25
zbins = np.geomspace(zP.min() * 1.01, zP.max(), nbins + 1)
zcent = 0.5 * (zbins[1:] + zbins[:-1])


def bin_resid(z, r, w, edges):
    means = []
    sems = []
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1])
        if m.sum() > 0:
            wi = w[m]
            ri = r[m]
            mu = np.sum(wi * ri) / np.sum(wi)
            sem = 1.0 / np.sqrt(np.sum(wi))
            means.append(mu)
            sems.append(sem)
        else:
            means.append(np.nan)
            sems.append(np.nan)
    return np.array(means), np.array(sems)


w_P = 1.0 / muPe**2
res_V3_b, sem_V3_b = bin_resid(zP, resP_V3['resid'], w_P, zbins)
res_LCDM_b, sem_LCDM_b = bin_resid(zP, resP_LCDM['resid'], w_P, zbins)

ax.errorbar(zcent, res_V3_b, yerr=sem_V3_b, fmt='s-', color='C0',
            label=f'Model-A V_3 ({resP_V3["chi2_per_dof"]:.3f})', lw=2)
ax.errorbar(zcent, res_LCDM_b, yerr=sem_LCDM_b, fmt='o--', color='k',
            label=f'LCDM ({resP_LCDM["chi2_per_dof"]:.3f})', lw=2)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual mu_obs - mu_model (binned)')
ax.set_title(f'Pantheon+ residuals (chi^2/dof in legend)')
ax.legend()
ax.grid(alpha=0.3)

# (1,0) Union3 residuals
ax = axes[1, 0]
ax.errorbar(zU, resU_V3['resid'], yerr=muUe, fmt='s-', color='C0',
            label=f'Model-A V_3 ({resU_V3["chi2_per_dof"]:.3f})', lw=2, ms=8)
ax.errorbar(zU, resU_LCDM['resid'], yerr=muUe, fmt='o--', color='k',
            label=f'LCDM ({resU_LCDM["chi2_per_dof"]:.3f})', lw=2, ms=8)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual mu_obs - mu_model')
ax.set_title('Union3 residuals')
ax.legend()
ax.grid(alpha=0.3)

# (1,1) chi^2 bar chart
ax = axes[1, 1]
labels = ['Pantheon+', 'Union3', 'Combined']
V3_vals = [resP_V3['chi2_per_dof'], resU_V3['chi2_per_dof'],
           chi2_V3_total / dof_total]
LCDM_vals = [resP_LCDM['chi2_per_dof'], resU_LCDM['chi2_per_dof'],
             chi2_LCDM_total / dof_total]
V1_vals = [resP_V1['chi2_per_dof'], resU_V1['chi2_per_dof'],
           chi2_V1_total / dof_total]
x = np.arange(len(labels))
w = 0.27
ax.bar(x - w, V3_vals, w, label='Model-A V_3', color='C0')
ax.bar(x, LCDM_vals, w, label='LCDM', color='k')
ax.bar(x + w, V1_vals, w, label='V_1 (reference)', color='C2', alpha=0.6)
ax.axhline(1.0, color='gray', ls=':', alpha=0.5, label='ideal')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('chi^2 / dof')
ax.set_title('Goodness of fit (lower is better)')
ax.legend()
ax.grid(alpha=0.3, axis='y')
for i, (v3v, Lv, V1v) in enumerate(zip(V3_vals, LCDM_vals, V1_vals)):
    ax.text(i - w, v3v + 0.005, f'{v3v:.3f}', ha='center', fontsize=8)
    ax.text(i, Lv + 0.005, f'{Lv:.3f}', ha='center', fontsize=8)
    ax.text(i + w, V1v + 0.005, f'{V1v:.3f}', ha='center', fontsize=8)

plt.suptitle(
    f'Model-A V_3 (β̃ = {BETA_TILDE:.4f}, A_0 = 1/(12π)) vs LCDM | '
    f'Pantheon+ + Union3 SN distance fit',
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/G8b')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_pantheon_union3_fit.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir / 'v3_pantheon_union3_fit.pdf', bbox_inches='tight')

# Summary CSV
summary = pd.DataFrame([
    {'catalog': 'Pantheon+', 'model': 'Model-A V_3',
     'chi2': resP_V3['chi2'], 'dof': resP_V3['dof'],
     'chi2_per_dof': resP_V3['chi2_per_dof'],
     'DeltaM': resP_V3['DeltaM'], 'N': len(zP)},
    {'catalog': 'Pantheon+', 'model': 'LCDM',
     'chi2': resP_LCDM['chi2'], 'dof': resP_LCDM['dof'],
     'chi2_per_dof': resP_LCDM['chi2_per_dof'],
     'DeltaM': resP_LCDM['DeltaM'], 'N': len(zP)},
    {'catalog': 'Pantheon+', 'model': 'V_1 (reference)',
     'chi2': resP_V1['chi2'], 'dof': resP_V1['dof'],
     'chi2_per_dof': resP_V1['chi2_per_dof'],
     'DeltaM': resP_V1['DeltaM'], 'N': len(zP)},
    {'catalog': 'Union3', 'model': 'Model-A V_3',
     'chi2': resU_V3['chi2'], 'dof': resU_V3['dof'],
     'chi2_per_dof': resU_V3['chi2_per_dof'],
     'DeltaM': resU_V3['DeltaM'], 'N': len(zU)},
    {'catalog': 'Union3', 'model': 'LCDM',
     'chi2': resU_LCDM['chi2'], 'dof': resU_LCDM['dof'],
     'chi2_per_dof': resU_LCDM['chi2_per_dof'],
     'DeltaM': resU_LCDM['DeltaM'], 'N': len(zU)},
    {'catalog': 'Union3', 'model': 'V_1 (reference)',
     'chi2': resU_V1['chi2'], 'dof': resU_V1['dof'],
     'chi2_per_dof': resU_V1['chi2_per_dof'],
     'DeltaM': resU_V1['DeltaM'], 'N': len(zU)},
])
summary.to_csv(outdir / 'v3_fit_summary.csv', index=False)

print(f"\nPlots and summary saved to {outdir.resolve()}")
print("=" * 72)
