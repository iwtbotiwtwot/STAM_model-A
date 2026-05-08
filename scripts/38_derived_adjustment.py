"""
Script 38 - The "derived adjustment" between STAM and LCDM

Sean's framing (2026-05-07):
  "STAM's honest prediction is that supernovae lie on a flatter curve, and
   unadjusted distances will be different. But if we can get to LCDM through
   a derived adjustment, the derived adjustment suggests the source of the
   difference."

Approach:
  1. Compute mu_STAM(z) and mu_LCDM(z) at the same H_0 and Omega_m -- no
     DeltaM marginalization, no nuisance parameter.
  2. Form Delta_mu(z) = mu_LCDM(z) - mu_STAM(z). This is the "derived
     adjustment" needed to bring STAM's curve onto LCDM's curve.
  3. Characterize Delta_mu(z) shape: is it constant (calibration only),
     linear in z (like the historical bridge term b*z), quadratic, log,
     or structured around some z?
  4. Compare against catalog residuals (Pantheon+, Union3) to see whether
     the derived adjustment matches what data prefer.

If Delta_mu(z) has a clean functional form, that form IS STAM's prediction
of what mechanism the LCDM dark-energy fit is implicitly absorbing.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==================================================================
# Setup (same as script 37)
# ==================================================================
C_KMS = 299792.458
H0_KMS_MPC = 67.4

A0 = 0.026514
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
BETA_TILDE = OMEGA_M * (1.0 - A0)**2

def A_eq(a):
    return 1.0 - (1.0 - A0) * a**1.5
def Y_pot(a):
    return BETA_TILDE / (1.0 - A_eq(a))
def K_kin(a):
    return (3.0/8.0) * (1.0 - A0)**2 * a**3
def h2_STAM_raw(a):
    return (OMEGA_M * a**-3 + OMEGA_R * a**-4 + Y_pot(a)) / (1.0 - K_kin(a))
H2_TODAY_STAM = h2_STAM_raw(1.0)
def h2_STAM(a):
    return h2_STAM_raw(a) / H2_TODAY_STAM
def h2_LCDM(a):
    return OMEGA_M * a**-3 + OMEGA_R * a**-4 + (1.0 - OMEGA_M - OMEGA_R)

def cumtrapz0(y, x):
    dx = np.diff(x)
    midy = 0.5 * (y[1:] + y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy * dx)])

def mu_predict(z_obs, h2_func, n_grid=8000, z_max=None):
    if z_max is None:
        z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_KMS_MPC) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0

# ==================================================================
# Compute Delta_mu(z) on a smooth grid, characterize shape
# ==================================================================
z_grid = np.linspace(1e-3, 2.5, 600)
mu_LCDM = mu_predict(z_grid, h2_LCDM, z_max=3.0)
mu_STAM = mu_predict(z_grid, h2_STAM, z_max=3.0)
Delta_mu = mu_LCDM - mu_STAM

# Fit several candidate functional forms to Delta_mu(z):
# (1) constant         A
# (2) linear in z      A + B z
# (3) linear in log(1+z) A + B ln(1+z)
# (4) quadratic in z   A + B z + C z^2
# (5) saturating       A + B z / (1 + z/z_t)
def fit_polyz(z, dmu, deg):
    coefs = np.polyfit(z, dmu, deg)
    pred = np.polyval(coefs, z)
    rms = np.sqrt(np.mean((dmu - pred)**2))
    return coefs, pred, rms

def fit_log1pz(z, dmu):
    X = np.log1p(z)
    A, B = np.polyfit(X, dmu, 1)[::-1]   # slope, intercept reversed
    coefs = np.polyfit(X, dmu, 1)
    pred = np.polyval(coefs, X)
    rms = np.sqrt(np.mean((dmu - pred)**2))
    return coefs, pred, rms

c0, p0, rms0 = fit_polyz(z_grid, Delta_mu, 0)
c1, p1, rms1 = fit_polyz(z_grid, Delta_mu, 1)
c2, p2, rms2 = fit_polyz(z_grid, Delta_mu, 2)
clog, plog, rmslog = fit_log1pz(z_grid, Delta_mu)

print("=" * 72)
print("Functional fits to Delta_mu(z) = mu_LCDM(z) - mu_STAM(z)")
print("=" * 72)
print(f"  Constant:           Delta_mu = {c0[0]:+.4f}")
print(f"     RMS residual:    {rms0:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p0)):.4f})")
print(f"  Linear in z:        Delta_mu = {c1[1]:+.4f} + {c1[0]:+.4f}*z")
print(f"     RMS residual:    {rms1:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p1)):.4f})")
print(f"  Quadratic in z:     Delta_mu = {c2[2]:+.4f} + {c2[1]:+.4f}*z + {c2[0]:+.4f}*z^2")
print(f"     RMS residual:    {rms2:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - p2)):.4f})")
print(f"  Linear in ln(1+z):  Delta_mu = {clog[1]:+.4f} + {clog[0]:+.4f}*ln(1+z)")
print(f"     RMS residual:    {rmslog:.4f} mag  (max abs = {np.max(np.abs(Delta_mu - plog)):.4f})")

# Find which fit is best
fits = [('constant', rms0, p0),
        ('linear_z', rms1, p1),
        ('quadratic', rms2, p2),
        ('linear_log1pz', rmslog, plog)]
best = min(fits, key=lambda x: x[1])
print(f"\n  Best simple fit:    {best[0]}  (RMS = {best[1]:.4f} mag)")

# ==================================================================
# What does Delta_mu(z) physically correspond to?
# Convert to luminosity-distance ratio:
# d_L_LCDM/d_L_STAM = 10^(Delta_mu/5)
# Or as a fractional difference at each z:
# (d_L_LCDM - d_L_STAM)/d_L_STAM = 10^(Delta_mu/5) - 1
# ==================================================================
fractional_dL = 10.0**(Delta_mu/5.0) - 1.0
# in Mpc:
DC_LCDM_grid = mu_predict(z_grid, h2_LCDM, z_max=3.0) - 25.0
dL_LCDM_grid = 10.0**(DC_LCDM_grid/5.0)
dL_STAM_grid = 10.0**((mu_STAM - 25.0)/5.0)
diff_dL_Mpc = dL_LCDM_grid - dL_STAM_grid

# z=0 limit of the difference
print(f"\n  Delta_mu(z->0)              = {Delta_mu[0]:+.4f} mag")
print(f"  Delta_mu(z=0.5)             = {np.interp(0.5, z_grid, Delta_mu):+.4f} mag")
print(f"  Delta_mu(z=1.0)             = {np.interp(1.0, z_grid, Delta_mu):+.4f} mag")
print(f"  Delta_mu(z=2.0)             = {np.interp(2.0, z_grid, Delta_mu):+.4f} mag")

# ==================================================================
# Compare to catalog residuals (no DeltaM marginalization)
# ==================================================================
DATA = Path('data')
pantheon = pd.read_csv(DATA / 'pantheon.csv')
mask = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
pantheon_cosmo = pantheon[mask].copy()
zP = pantheon_cosmo['zCMB'].values
muP = pantheon_cosmo['MU_SH0ES'].values
muP_err = pantheon_cosmo['MU_SH0ES_ERR_DIAG'].values

union3 = pd.read_csv(DATA / 'union3_bins.csv')
zU = union3['z'].values
muU = union3['mb'].values
muU_err = np.full_like(muU, 0.05)

mu_LCDM_P = mu_predict(zP, h2_LCDM, z_max=zP.max()*1.05)
mu_STAM_P = mu_predict(zP, h2_STAM, z_max=zP.max()*1.05)
mu_LCDM_U = mu_predict(zU, h2_LCDM)
mu_STAM_U = mu_predict(zU, h2_STAM)

# raw residuals (no DeltaM)
raw_res_LCDM_P = muP - mu_LCDM_P
raw_res_STAM_P = muP - mu_STAM_P
raw_res_LCDM_U = muU - mu_LCDM_U
raw_res_STAM_U = muU - mu_STAM_U

# Fit a constant DeltaM (zero-point) to each model+catalog so we can see
# the SHAPE of the leftover residual (after just the calibration is removed):
def best_const(res, err):
    w = 1.0/err**2
    return np.sum(w*res)/np.sum(w)

dM_LCDM_P = best_const(raw_res_LCDM_P, muP_err)
dM_STAM_P = best_const(raw_res_STAM_P, muP_err)
dM_LCDM_U = best_const(raw_res_LCDM_U, muU_err)
dM_STAM_U = best_const(raw_res_STAM_U, muU_err)

# residual after constant offset removed -- this shows shape only
shape_LCDM_P = raw_res_LCDM_P - dM_LCDM_P
shape_STAM_P = raw_res_STAM_P - dM_STAM_P
shape_LCDM_U = raw_res_LCDM_U - dM_LCDM_U
shape_STAM_U = raw_res_STAM_U - dM_STAM_U

# Bin Pantheon residuals for visual clarity
def bin_w(z, r, w, edges):
    out_z, out_r, out_e = [], [], []
    for i in range(len(edges)-1):
        m = (z >= edges[i]) & (z < edges[i+1])
        if m.sum() > 2:
            wi = w[m]; ri = r[m]
            mean_z = np.sum(wi*z[m])/np.sum(wi)
            mean_r = np.sum(wi*ri)/np.sum(wi)
            sem = 1.0/np.sqrt(np.sum(wi))
            out_z.append(mean_z)
            out_r.append(mean_r)
            out_e.append(sem)
    return np.array(out_z), np.array(out_r), np.array(out_e)

zbins = np.geomspace(zP.min()*1.01, zP.max(), 26)
zb_P, sLb_P, sLe_P = bin_w(zP, shape_LCDM_P, 1.0/muP_err**2, zbins)
_,    sSb_P, sSe_P = bin_w(zP, shape_STAM_P, 1.0/muP_err**2, zbins)

# ==================================================================
# Plotting
# ==================================================================
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.28)

# Panel A: mu(z) absolute, both models
ax = fig.add_subplot(gs[0, 0])
ax.plot(z_grid, mu_LCDM, lw=2.0, color='k', label='LCDM (no offset)')
ax.plot(z_grid, mu_STAM, lw=2.0, ls='--', color='C1', label='STAM (no offset)')
ax.set_xlabel('redshift z')
ax.set_ylabel('mu (mag)')
ax.set_title('Predicted distance modulus (raw, same H_0 = 67.4)')
ax.legend()
ax.grid(alpha=0.3)
ax.set_xscale('log')

# Panel B: Delta_mu(z) and fits
ax = fig.add_subplot(gs[0, 1])
ax.plot(z_grid, Delta_mu, lw=3.0, color='C0', label='Delta_mu = LCDM - STAM')
ax.plot(z_grid, p1, lw=1.5, ls='--', color='C3',
        label=f'linear: {c1[1]:+.3f} + {c1[0]:+.3f} z')
ax.plot(z_grid, p2, lw=1.5, ls=':', color='C2',
        label=f'quadratic (RMS={rms2:.4f})')
ax.plot(z_grid, plog, lw=1.5, ls='-.', color='C4',
        label=f'linear in ln(1+z): {clog[1]:+.3f} + {clog[0]:+.3f}*ln(1+z)')
ax.axhline(0, color='gray', alpha=0.4)
ax.set_xlabel('redshift z')
ax.set_ylabel('Delta_mu (mag)')
ax.set_title('THE DERIVED ADJUSTMENT  (LCDM minus STAM)')
ax.legend(loc='best', fontsize=9)
ax.grid(alpha=0.3)

# Panel C: Pantheon shape residuals (offset removed) vs Delta_mu prediction
ax = fig.add_subplot(gs[1, 0])
ax.errorbar(zb_P, sLb_P, yerr=sLe_P, fmt='o-', color='k',
            label='Pantheon+ residual vs LCDM (offset removed)', lw=1.5)
ax.errorbar(zb_P, sSb_P, yerr=sSe_P, fmt='s--', color='C1',
            label='Pantheon+ residual vs STAM (offset removed)', lw=1.5)
ax.axhline(0, color='gray', alpha=0.4)
ax.set_xlabel('redshift z')
ax.set_ylabel('mu residual (mag)')
ax.set_title('Pantheon+ binned residuals')
ax.set_xscale('log')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel D: Union3 shape residuals
ax = fig.add_subplot(gs[1, 1])
ax.errorbar(zU, shape_LCDM_U, yerr=muU_err, fmt='o-', color='k',
            label='Union3 residual vs LCDM (offset removed)', lw=1.5, ms=8)
ax.errorbar(zU, shape_STAM_U, yerr=muU_err, fmt='s--', color='C1',
            label='Union3 residual vs STAM (offset removed)', lw=1.5, ms=8)
ax.axhline(0, color='gray', alpha=0.4)
ax.set_xlabel('redshift z')
ax.set_ylabel('mu residual (mag)')
ax.set_title('Union3 residuals')
ax.set_xscale('log')
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel E: fractional d_L difference (more physical)
ax = fig.add_subplot(gs[2, 0])
ax.plot(z_grid, fractional_dL * 100, lw=2.5, color='C0')
ax.axhline(0, color='gray', alpha=0.4)
ax.set_xlabel('redshift z')
ax.set_ylabel('(d_L_LCDM - d_L_STAM) / d_L_STAM  [%]')
ax.set_title('Luminosity distance fractional adjustment')
ax.grid(alpha=0.3)

# Panel F: w_eff context
ax = fig.add_subplot(gs[2, 1])
# recompute w_eff briefly for context
def rho_DE(a_):
    return Y_pot(a_) + K_kin(a_) * h2_STAM_raw(a_)
a_arr = np.logspace(-2, -0.001, 300)
z_arr = 1/a_arr - 1
w_eff = -1.0 - (1.0/3.0) * np.gradient(np.log(rho_DE(a_arr)), np.log(a_arr))
ax.plot(z_arr, w_eff, lw=2.5, color='C2', label='STAM w_eff(z)')
ax.axhline(-1, ls=':', color='k', alpha=0.6, label='LCDM w = -1')
ax.set_xlabel('redshift z')
ax.set_ylabel('w_eff(z)')
ax.set_title('STAM effective dark-energy EoS')
ax.legend()
ax.grid(alpha=0.3)
ax.set_xlim(0, 2)
ax.set_ylim(-1.4, 0.0)

plt.suptitle(
    "STAM 'derived adjustment' to bridge to LCDM | "
    f"max |Delta_mu| in 0<z<2.5 = {np.max(np.abs(Delta_mu)):.3f} mag",
    fontsize=13
)

outdir = Path('reports/script_38')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'derived_adjustment.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir / 'derived_adjustment.pdf', bbox_inches='tight')

# Save adjustment table
adj_df = pd.DataFrame({
    'z': z_grid,
    'mu_LCDM': mu_LCDM,
    'mu_STAM': mu_STAM,
    'Delta_mu': Delta_mu,
    'fractional_dL_pct': fractional_dL * 100,
})
adj_df.to_csv(outdir / 'derived_adjustment.csv', index=False)
print(f"\nSaved: {outdir.resolve()}")
print("=" * 72)
