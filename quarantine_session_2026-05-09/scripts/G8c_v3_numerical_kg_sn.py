"""
G8c - Model-A V_3 with full numerical Klein-Gordon dynamics, SN distance test

Replaces the tracking-ansatz proxy used in G8a/G8b with the actual V_3
cosmology from full KG integration. The tracking ansatz is the slow-roll
attractor for V_1 but is NOT the dynamics V_3 honors — V_3 has a true
minimum at A_0 so the field oscillates around it under Hubble damping.
G8a found tracking-vs-numerical disagreement of factor ~6 in A and ~18
in omega today. This script uses the numerical KG cosmology directly.

Setup:
- A_0 = 1/(12 pi)                          (G7 commitment)
- alpha/beta = [A_0 / (1-A_0)]^2           (V_3 minimum at A_0)
- Omega_m = 0.315, Omega_r = 9.2e-5        (Model-A internal content)
- beta_tilde calibrated such that the NUMERICAL KG integration produces
  h^2(today) = 1, with tracking initial conditions at a = 1e-3.
- h(a) for SN distance integration comes from interpolating the
  numerical KG trajectory, not from tracking.

Apples-to-apples covariance for both catalogs:
- Pantheon+: published STAT+SYS covariance, submatrixed to 1578 cosmological SNe
- Union3:    published inverse-covariance from FITS

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
# Constants
# ===================================================================
A0 = 1.0 / (12.0 * math.pi)
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
ALPHA_OVER_BETA = (A0 / (1.0 - A0))**2
H0_KMS_MPC = 73.04
H0_LCDM_KMS_MPC = 67.4
C_KMS = 299792.458

# Integration parameters
LNA_INIT = math.log(1e-3)
TARGET_DLNA = 0.003
MAX_STEPS = 500_000


# ===================================================================
# V_3 potential
# ===================================================================
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


# ===================================================================
# Klein-Gordon + Friedmann RK4 integration
# ===================================================================
def integrate_kg(beta_tilde, lna_init=LNA_INIT, target_dlna=TARGET_DLNA,
                 max_steps=MAX_STEPS):
    """Integrate the coupled KG + Friedmann system in tau = H_0 t.

    State: y = [A, omega, lna] with omega = dA/dtau.
    Equations:
        dA/dtau     = omega
        domega/dtau = -3 h omega - 3 (Y'(A) - Omega_m a^-3)
        dlna/dtau   = h
        h^2         = Omega_m a^-3 + Omega_r a^-4 + omega^2/6 + Y(A)

    Initial condition: A in tracking slow-roll equilibrium at lna_init,
    omega from the tracking velocity formula. This is the natural matter-
    dominated-era starting point.

    Returns: (lna_arr, A_arr, omega_arr, h_arr) sampled along the trajectory.
    Final entry interpolated to lna=0 (today).
    """
    a_init = math.exp(lna_init)
    A_init = A_eq_at(a_init, beta_tilde)
    rho_m_init = OMEGA_M * a_init**-3
    rho_r_init = OMEGA_R * a_init**-4
    Y_init = Y_pot(A_init, beta_tilde)
    # h_init from Friedmann (potential negligible at very early times)
    h_init = math.sqrt(rho_m_init + rho_r_init + Y_init)
    # tracking velocity
    omega_init = -3.0 * h_init * Yp(A_init, beta_tilde) / Ypp(A_init, beta_tilde)

    def rhs(y):
        A, omega, lna = y
        a_ = math.exp(lna)
        Y_v = Y_pot(A, beta_tilde)
        Yp_v = Yp(A, beta_tilde)
        rho_m = OMEGA_M * a_**-3
        rho_r = OMEGA_R * a_**-4
        h2_ = rho_m + rho_r + omega * omega / 6.0 + Y_v
        h_ = math.sqrt(max(h2_, 1e-30))
        dA = omega
        domega = -3.0 * h_ * omega - 3.0 * (Yp_v - rho_m)
        dlna = h_
        return np.array([dA, domega, dlna])

    y = np.array([A_init, omega_init, lna_init])
    traj_lna = [lna_init]
    traj_A = [A_init]
    traj_omega = [omega_init]
    traj_h = [h_init]

    for step in range(max_steps):
        h_now = math.sqrt(
            OMEGA_M * math.exp(-3 * y[2]) + OMEGA_R * math.exp(-4 * y[2])
            + y[1]**2 / 6.0 + Y_pot(y[0], beta_tilde)
        )
        dtau = target_dlna / max(h_now, 1e-10)
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dtau * k1)
        k3 = rhs(y + 0.5 * dtau * k2)
        k4 = rhs(y + dtau * k3)
        y_new = y + (dtau / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if y_new[2] >= 0.0:
            frac = -y[2] / (y_new[2] - y[2])
            y_today = y + frac * (y_new - y)
            h_today = math.sqrt(
                OMEGA_M + OMEGA_R + y_today[1]**2 / 6.0
                + Y_pot(y_today[0], beta_tilde)
            )
            traj_lna.append(0.0)
            traj_A.append(y_today[0])
            traj_omega.append(y_today[1])
            traj_h.append(h_today)
            return (np.array(traj_lna), np.array(traj_A),
                    np.array(traj_omega), np.array(traj_h))
        y = y_new
        traj_lna.append(y[2])
        traj_A.append(y[0])
        traj_omega.append(y[1])
        h_step = math.sqrt(
            OMEGA_M * math.exp(-3 * y[2]) + OMEGA_R * math.exp(-4 * y[2])
            + y[1]**2 / 6.0 + Y_pot(y[0], beta_tilde)
        )
        traj_h.append(h_step)
    raise RuntimeError(f"KG integration did not reach today in {max_steps} steps")


def h2_today_residual(beta_tilde):
    _, _, _, h_traj = integrate_kg(beta_tilde)
    return h_traj[-1]**2 - 1.0


# ===================================================================
# Calibrate beta_tilde such that NUMERICAL h(today) = 1
# ===================================================================
print("=" * 72)
print("Model-A V_3 numerical Klein-Gordon SN distance test (G8c)")
print("=" * 72)
print(f"A_0 = 1/(12 pi)  = {A0:.8f}")
print(f"alpha/beta       = {ALPHA_OVER_BETA:.6e}")
print(f"Omega_m          = {OMEGA_M}")
print()
print("Calibrating beta_tilde so the NUMERICAL KG gives h^2(today) = 1...")
print("(tracking-ansatz value was 0.6473 but tracking does not honor V_3)")

# Probe to find a bracket
betas_probe = [0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90]
residuals_probe = []
for b in betas_probe:
    try:
        r = h2_today_residual(b)
        residuals_probe.append(r)
        print(f"  beta_tilde = {b:.4f}  h^2(today) = {r + 1:.6f}  "
              f"residual = {r:+.6f}")
    except Exception as e:
        residuals_probe.append(float('nan'))
        print(f"  beta_tilde = {b:.4f}  INTEGRATION FAILED: {e}")

residuals_probe = np.array(residuals_probe)
# Find sign change
sign_changes = []
for i in range(len(residuals_probe) - 1):
    a_, b_ = residuals_probe[i], residuals_probe[i + 1]
    if np.isfinite(a_) and np.isfinite(b_) and (a_ * b_ < 0):
        sign_changes.append(i)

if not sign_changes:
    raise RuntimeError("No bracket found for numerical KG h^2(today)=1 closure.")

i = sign_changes[0]
beta_low, beta_high = betas_probe[i], betas_probe[i + 1]
print(f"\nBracketed in beta_tilde in [{beta_low}, {beta_high}]; refining via brentq...")
BETA_TILDE_NUM = brentq(h2_today_residual, beta_low, beta_high, xtol=1e-6)
print(f"NUMERICAL-KG calibrated beta_tilde = {BETA_TILDE_NUM:.6f}")

# Final integration trajectory
lna_arr, A_arr, omega_arr, h_arr = integrate_kg(BETA_TILDE_NUM)
a_arr = np.exp(lna_arr)
h2_arr = h_arr**2
print(f"\nFinal trajectory:")
print(f"  steps integrated: {len(lna_arr)}")
print(f"  h^2(today)        = {h2_arr[-1]:.8f}")
print(f"  A(today)          = {A_arr[-1]:.6f}")
print(f"  omega(today)      = {omega_arr[-1]:+.6f}")
print(f"  Y_pot(today)      = {Y_pot(A_arr[-1], BETA_TILDE_NUM):.6f}")
print(f"  Omega_kin(today)  = {omega_arr[-1]**2 / 6.0:.6f}")
print(f"  Omega_DE_total    = {Y_pot(A_arr[-1], BETA_TILDE_NUM) + omega_arr[-1]**2/6.0:.6f}")

# ===================================================================
# h^2(a) interpolator for distance integration
# ===================================================================
# Sort by lna ascending to ensure monotonic for interpolation
idx_sort = np.argsort(lna_arr)
lna_sorted = lna_arr[idx_sort]
h2_sorted = h2_arr[idx_sort]


def h2_v3_num(a):
    """Numerical KG h^2(a). Uses log-linear interpolation in lna."""
    a = np.atleast_1d(a)
    lna = np.log(a)
    # clip to integration range
    lna_clipped = np.clip(lna, lna_sorted[0], lna_sorted[-1])
    return np.interp(lna_clipped, lna_sorted, h2_sorted)


def h2_LCDM(a):
    Om = 0.315
    Or = 9.2e-5
    OL = 1.0 - Om - Or
    return Om * a**-3 + Or * a**-4 + OL


# ===================================================================
# Distance modulus + Mahalanobis chi^2
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


def fit_offset_chi2_mahalanobis_full_cov(z_obs, mu_obs, cov, h2_func, H0_used):
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
        'chi2': chi2, 'dof': dof, 'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM, 'mu_model': mu_model, 'resid': resid,
    }


def fit_offset_chi2_mahalanobis_inv_cov(z_obs, mu_obs, inv_cov, h2_func, H0_used):
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
        'chi2': chi2, 'dof': dof, 'chi2_per_dof': chi2 / dof,
        'DeltaM': DeltaM, 'mu_model': mu_model, 'resid': resid,
    }


# ===================================================================
# Load data with proper covariance for both catalogs
# ===================================================================
DATA = Path('data')

# Pantheon+
print()
print("Loading Pantheon+ data and STAT+SYS covariance...")
pantheon = pd.read_csv(DATA / 'pantheon.csv')
mask = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
mask_idx = np.where(mask.values)[0]
zP = pantheon.loc[mask, 'zCMB'].values
muP = pantheon.loc[mask, 'MU_SH0ES'].values

with open(DATA / 'Pantheon+SH0ES_STAT+SYS.cov', 'r') as f:
    N = int(f.readline().strip())
    flat = np.fromfile(f, sep='\n', dtype=np.float64, count=N * N)
C_full = flat.reshape(N, N)
C_full = 0.5 * (C_full + C_full.T)
covP = C_full[np.ix_(mask_idx, mask_idx)]
print(f"  Pantheon+: {len(zP)} cosmological SNe, cov = {covP.shape}")

# Union3
print("Loading Union3 FITS...")
with fits.open(DATA / 'mu_mat_union3_cosmo2_mu.fits') as hdul:
    fits_data = hdul[0].data
zU = fits_data[0, 1:].astype(float)
muU = fits_data[1:, 0].astype(float)
invCovU = fits_data[1:, 1:].astype(float)
print(f"  Union3:    {len(zU)} bins, inv_cov = {invCovU.shape}")

# ===================================================================
# Compute chi^2 for V_3 numerical, V_3 tracking (reference), and LCDM
# ===================================================================
print()
print("=" * 72)
print("PANTHEON+ (Mahalanobis with STAT+SYS covariance)")
print("=" * 72)
resP_V3num = fit_offset_chi2_mahalanobis_full_cov(zP, muP, covP, h2_v3_num,
                                                   H0_KMS_MPC)
resP_LCDM = fit_offset_chi2_mahalanobis_full_cov(zP, muP, covP, h2_LCDM,
                                                  H0_LCDM_KMS_MPC)
print(f"  Model-A V_3 (numerical KG) : chi^2 = {resP_V3num['chi2']:8.2f}, "
      f"chi^2/dof = {resP_V3num['chi2_per_dof']:.4f}, "
      f"DeltaM = {resP_V3num['DeltaM']:+.4f}")
print(f"  LCDM                       : chi^2 = {resP_LCDM['chi2']:8.2f}, "
      f"chi^2/dof = {resP_LCDM['chi2_per_dof']:.4f}, "
      f"DeltaM = {resP_LCDM['DeltaM']:+.4f}")
print(f"  delta(chi^2) V_3 - LCDM    = {resP_V3num['chi2'] - resP_LCDM['chi2']:+.2f}")

print()
print("=" * 72)
print("UNION3 (Mahalanobis with bin-bin covariance)")
print("=" * 72)
resU_V3num = fit_offset_chi2_mahalanobis_inv_cov(zU, muU, invCovU, h2_v3_num,
                                                  H0_KMS_MPC)
resU_LCDM = fit_offset_chi2_mahalanobis_inv_cov(zU, muU, invCovU, h2_LCDM,
                                                 H0_LCDM_KMS_MPC)
print(f"  Model-A V_3 (numerical KG) : chi^2 = {resU_V3num['chi2']:8.2f}, "
      f"chi^2/dof = {resU_V3num['chi2_per_dof']:.4f}, "
      f"DeltaM = {resU_V3num['DeltaM']:+.4f}")
print(f"  LCDM                       : chi^2 = {resU_LCDM['chi2']:8.2f}, "
      f"chi^2/dof = {resU_LCDM['chi2_per_dof']:.4f}, "
      f"DeltaM = {resU_LCDM['DeltaM']:+.4f}")
print(f"  delta(chi^2) V_3 - LCDM    = {resU_V3num['chi2'] - resU_LCDM['chi2']:+.2f}")

print()
print("=" * 72)
print("COMBINED")
print("=" * 72)
chi2_V3num_total = resP_V3num['chi2'] + resU_V3num['chi2']
chi2_LCDM_total = resP_LCDM['chi2'] + resU_LCDM['chi2']
dof_total = resP_V3num['dof'] + resU_V3num['dof']
print(f"  Model-A V_3 (numerical KG) : chi^2 = {chi2_V3num_total:8.2f}  / "
      f"{dof_total} dof  =  {chi2_V3num_total / dof_total:.4f}")
print(f"  LCDM                       : chi^2 = {chi2_LCDM_total:8.2f}  / "
      f"{dof_total} dof  =  {chi2_LCDM_total / dof_total:.4f}")
print(f"  delta(chi^2) V_3 - LCDM    = {chi2_V3num_total - chi2_LCDM_total:+.2f}")
print()
print("Two-tables comparison (apples-to-apples Mahalanobis, both catalogs):")
print()
print("  Model-A V_3 (numerical KG):")
print(f"    Pantheon+ chi^2/dof = {resP_V3num['chi2_per_dof']:.4f}")
print(f"    Union3    chi^2/dof = {resU_V3num['chi2_per_dof']:.4f}")
print(f"    Combined  chi^2/dof = {chi2_V3num_total / dof_total:.4f}")
print()
print("  LCDM:")
print(f"    Pantheon+ chi^2/dof = {resP_LCDM['chi2_per_dof']:.4f}")
print(f"    Union3    chi^2/dof = {resU_LCDM['chi2_per_dof']:.4f}")
print(f"    Combined  chi^2/dof = {chi2_LCDM_total / dof_total:.4f}")

# ===================================================================
# Plotting
# ===================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 0,0: numerical KG trajectory in A vs lna
ax = axes[0, 0]
ax.plot(lna_arr, A_arr, lw=2, color='C0', label='V_3 numerical KG A(a)')
ax.axhline(A0, ls=':', color='gray', alpha=0.7,
           label=f'A_0 = 1/(12π) = {A0:.5f}')
A_eq_traj = np.array([A_eq_at(a, BETA_TILDE_NUM) for a in a_arr])
ax.plot(lna_arr, A_eq_traj, lw=1.5, ls='--', color='C1', alpha=0.6,
        label='tracking A_eq (proxy, for reference)')
ax.set_xlabel('log10(a) implied by lna')
ax.set_ylabel('A')
ax.set_title('V_3 field trajectory: numerical KG vs tracking proxy')
ax.legend()
ax.grid(alpha=0.3)

# Panel 0,1: h(z) for V_3 num and LCDM
ax = axes[0, 1]
z_grid = np.linspace(0, 5, 600)
a_grid = 1.0 / (1.0 + z_grid)
ax.plot(z_grid, np.sqrt(h2_v3_num(a_grid)) * H0_KMS_MPC, lw=2, color='C0',
        label=f'V_3 numerical KG (H_0={H0_KMS_MPC})')
ax.plot(z_grid, np.sqrt(h2_LCDM(a_grid)) * H0_LCDM_KMS_MPC, lw=2, ls='--',
        color='k', label=f'LCDM (H_0={H0_LCDM_KMS_MPC})')
ax.set_xlabel('redshift z')
ax.set_ylabel('H(z) [km/s/Mpc]')
ax.set_title('Hubble parameter — apples-to-apples cosmology comparison')
ax.legend()
ax.grid(alpha=0.3)

# Panel 1,0: Pantheon+ binned residuals
ax = axes[1, 0]
nbins = 25
zbins = np.geomspace(zP.min() * 1.01, zP.max(), nbins + 1)
zcent = 0.5 * (zbins[1:] + zbins[:-1])
sigma_diag_P = np.sqrt(np.diag(covP))
w_P = 1.0 / sigma_diag_P**2

def bin_resid(z, r, w, edges):
    means, sems = [], []
    for i in range(len(edges) - 1):
        m = (z >= edges[i]) & (z < edges[i + 1])
        if m.sum() > 0:
            wi, ri = w[m], r[m]
            means.append(np.sum(wi * ri) / np.sum(wi))
            sems.append(1.0 / np.sqrt(np.sum(wi)))
        else:
            means.append(np.nan)
            sems.append(np.nan)
    return np.array(means), np.array(sems)

res_V3_b, sem_V3_b = bin_resid(zP, resP_V3num['resid'], w_P, zbins)
res_LCDM_b, sem_LCDM_b = bin_resid(zP, resP_LCDM['resid'], w_P, zbins)
ax.errorbar(zcent, res_V3_b, yerr=sem_V3_b, fmt='s-', color='C0',
            label=f'V_3 numerical KG ({resP_V3num["chi2_per_dof"]:.3f})', lw=2)
ax.errorbar(zcent, res_LCDM_b, yerr=sem_LCDM_b, fmt='o--', color='k',
            label=f'LCDM ({resP_LCDM["chi2_per_dof"]:.3f})', lw=2)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual (mag)')
ax.set_title('Pantheon+ binned residuals (Mahalanobis-weighted bins)')
ax.legend()
ax.grid(alpha=0.3)

# Panel 1,1: Union3 residuals per bin
ax = axes[1, 1]
sigma_U = np.sqrt(np.diag(np.linalg.inv(invCovU)))
ax.errorbar(zU, resU_V3num['resid'], yerr=sigma_U, fmt='s-', color='C0',
            label=f'V_3 numerical KG ({resU_V3num["chi2_per_dof"]:.3f})',
            lw=2, ms=8)
ax.errorbar(zU, resU_LCDM['resid'], yerr=sigma_U, fmt='o--', color='k',
            label=f'LCDM ({resU_LCDM["chi2_per_dof"]:.3f})',
            lw=2, ms=8)
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('redshift z')
ax.set_ylabel('residual (mag)')
ax.set_title('Union3 residuals per bin')
ax.legend()
ax.grid(alpha=0.3)

plt.suptitle(
    f'Model-A V_3 numerical KG (β̃_num = {BETA_TILDE_NUM:.4f}) vs LCDM | '
    f'apples-to-apples Mahalanobis χ² for both catalogs',
    fontsize=12
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/G8c')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'v3_numerical_kg_sn_fit.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir / 'v3_numerical_kg_sn_fit.pdf', bbox_inches='tight')

# Summary CSV
summary = pd.DataFrame([
    {'catalog': 'Pantheon+', 'model': 'V_3 numerical KG',
     'chi2': resP_V3num['chi2'], 'dof': resP_V3num['dof'],
     'chi2_per_dof': resP_V3num['chi2_per_dof'],
     'DeltaM': resP_V3num['DeltaM'], 'N': len(zP)},
    {'catalog': 'Pantheon+', 'model': 'LCDM',
     'chi2': resP_LCDM['chi2'], 'dof': resP_LCDM['dof'],
     'chi2_per_dof': resP_LCDM['chi2_per_dof'],
     'DeltaM': resP_LCDM['DeltaM'], 'N': len(zP)},
    {'catalog': 'Union3', 'model': 'V_3 numerical KG',
     'chi2': resU_V3num['chi2'], 'dof': resU_V3num['dof'],
     'chi2_per_dof': resU_V3num['chi2_per_dof'],
     'DeltaM': resU_V3num['DeltaM'], 'N': len(zU)},
    {'catalog': 'Union3', 'model': 'LCDM',
     'chi2': resU_LCDM['chi2'], 'dof': resU_LCDM['dof'],
     'chi2_per_dof': resU_LCDM['chi2_per_dof'],
     'DeltaM': resU_LCDM['DeltaM'], 'N': len(zU)},
])
summary.to_csv(outdir / 'v3_numerical_kg_fit_summary.csv', index=False)

# Save trajectory CSV
traj_df = pd.DataFrame({
    'lna': lna_arr, 'a': a_arr, 'A': A_arr, 'omega': omega_arr,
    'h': h_arr, 'h2': h2_arr,
})
traj_df.to_csv(outdir / 'v3_kg_trajectory.csv', index=False)

print(f"\nPlots, summary, and trajectory saved to {outdir.resolve()}")
print("=" * 72)
