"""
Script 39b - Cross-catalog tension test with DES top-6 outliers removed

Same as 39, but removes the 6 SN in DES-Y5 with the largest absolute
residual against the best-fit LCDM. This tests whether the residual
DES-vs-Pantheon ~110 mmag tension is being driven by a small number
of extreme outliers (DES has known catastrophic outliers).
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ==================================================================
# Cosmological setup (same as scripts 36-39)
# ==================================================================
C_KMS = 299792.458
H0_KMS_MPC = 67.4
A0 = 0.026514
OMEGA_M = 0.315
OMEGA_R = 9.2e-5
BETA_TILDE = OMEGA_M * (1.0 - A0)**2

def A_eq(a):     return 1.0 - (1.0 - A0) * a**1.5
def Y_pot(a):    return BETA_TILDE / (1.0 - A_eq(a))
def K_kin(a):    return (3.0/8.0) * (1.0 - A0)**2 * a**3
def h2_STAM_raw(a):
    return (OMEGA_M*a**-3 + OMEGA_R*a**-4 + Y_pot(a))/(1.0 - K_kin(a))
H2_TODAY_STAM = h2_STAM_raw(1.0)
def h2_STAM(a): return h2_STAM_raw(a) / H2_TODAY_STAM
def h2_LCDM(a):
    return OMEGA_M*a**-3 + OMEGA_R*a**-4 + (1.0 - OMEGA_M - OMEGA_R)

def cumtrapz0(y, x):
    dx = np.diff(x); midy = 0.5*(y[1:]+y[:-1])
    return np.concatenate([[0.0], np.cumsum(midy*dx)])

def mu_predict(z_obs, h2_func, n_grid=8000, z_max=None):
    if z_max is None:
        z_max = max(np.max(z_obs)*1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0/(1.0+z_grid)
    inv_h = 1.0/np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS/H0_KMS_MPC)*cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0+z_obs)*DC_obs
    return 5.0*np.log10(dL_obs)+25.0

def fit_offset_chi2(z_obs, mu_obs, mu_err, h2_func):
    mu_model = mu_predict(z_obs, h2_func)
    w = 1.0/mu_err**2
    DeltaM = np.sum(w*(mu_obs-mu_model))/np.sum(w)
    resid = mu_obs - (mu_model + DeltaM)
    chi2 = np.sum(w*resid**2)
    return {'chi2':chi2, 'dof':len(z_obs)-1, 'chi2_per_dof':chi2/(len(z_obs)-1),
            'DeltaM':DeltaM, 'mu_model':mu_model, 'resid':resid,
            'z':z_obs, 'mu':mu_obs, 'mu_err':mu_err}

# ==================================================================
# Load catalogs
# ==================================================================
DATA = Path('data')
pantheon = pd.read_csv(DATA/'pantheon.csv')
mP = (pantheon['IS_CALIBRATOR']==0) & (pantheon['zCMB']>0.01)
zP = pantheon[mP]['zCMB'].values
muP = pantheon[mP]['MU_SH0ES'].values
muPe = pantheon[mP]['MU_SH0ES_ERR_DIAG'].values

union3 = pd.read_csv(DATA/'union3_bins.csv')
zU = union3['z'].values
muU = union3['mb'].values
muUe = np.full_like(muU, 0.05)

des = pd.read_csv(DATA/'des.csv')
mD = (des['zCMB']>0.01)
des_cosmo = des[mD].copy().reset_index(drop=True)

# ==================================================================
# Identify DES outliers from LCDM fit
# ==================================================================
zD_full = des_cosmo['zCMB'].values
muD_full = des_cosmo['MU'].values
muDe_full = des_cosmo['MUERR'].values

fit_DES_full = fit_offset_chi2(zD_full, muD_full, muDe_full, h2_LCDM)
abs_resid = np.abs(fit_DES_full['resid'])
order = np.argsort(abs_resid)[::-1]    # largest first
N_OUTLIERS = 6
outlier_idx = order[:N_OUTLIERS]
keep_idx = order[N_OUTLIERS:]

print("=" * 72)
print(f"Top {N_OUTLIERS} DES outliers identified by |residual| from LCDM fit:")
print("=" * 72)
print(f"{'CID':>10} {'z':>8} {'MU':>10} {'MU_pred':>10} {'resid':>10} {'|resid|':>10}")
for k, i in enumerate(outlier_idx):
    cid = des_cosmo['CID'].iloc[i]
    z = zD_full[i]
    mu = muD_full[i]
    mu_pred = fit_DES_full['mu_model'][i] + fit_DES_full['DeltaM']
    res = fit_DES_full['resid'][i]
    print(f"{cid:>10} {z:>8.4f} {mu:>10.3f} {mu_pred:>10.3f} {res:>+10.3f} {abs_resid[i]:>10.3f}")

# Also show how they compare to literature outliers we tracked earlier:
prior_top6 = [1299503, 1289982, 1291149, 1294436, 1900801, 1319870]
print("\nFor reference, top 6 from prior analysis (results/supernova_discrepancy):")
for cid in prior_top6:
    print(f"  CID {cid}")

print()
print(f"DES-Y5 size: full = {len(zD_full)}, after removing top-{N_OUTLIERS} = {len(keep_idx)}")
print()

# Build trimmed DES arrays
zD = zD_full[keep_idx]
muD = muD_full[keep_idx]
muDe = muDe_full[keep_idx]

# ==================================================================
# Independent fits (with DES trimmed)
# ==================================================================
fits_LCDM = {
    'Pantheon+': fit_offset_chi2(zP, muP, muPe, h2_LCDM),
    'Union3'   : fit_offset_chi2(zU, muU, muUe, h2_LCDM),
    'DES-Y5'   : fit_offset_chi2(zD, muD, muDe, h2_LCDM),
}
fits_STAM = {
    'Pantheon+': fit_offset_chi2(zP, muP, muPe, h2_STAM),
    'Union3'   : fit_offset_chi2(zU, muU, muUe, h2_STAM),
    'DES-Y5'   : fit_offset_chi2(zD, muD, muDe, h2_STAM),
}

print("Independent best-fit DeltaM (offset) per catalog (DES trimmed):")
print(f"{'catalog':<12} {'LCDM DeltaM':>14} {'STAM DeltaM':>14} "
      f"{'LCDM chi2/dof':>14} {'STAM chi2/dof':>14}")
for cat in ('Pantheon+', 'Union3', 'DES-Y5'):
    rL = fits_LCDM[cat]; rS = fits_STAM[cat]
    print(f"{cat:<12} {rL['DeltaM']:>+14.4f} {rS['DeltaM']:>+14.4f} "
          f"{rL['chi2_per_dof']:>14.4f} {rS['chi2_per_dof']:>14.4f}")

# ==================================================================
# Cross-catalog tensions
# ==================================================================
print()
print("OBSERVED inter-catalog DeltaM tension (under LCDM fit, DES trimmed):")
pairs = [('Pantheon+','Union3'), ('Pantheon+','DES-Y5'), ('Union3','DES-Y5')]
obs_tensions = {}
for a,b in pairs:
    diff = fits_LCDM[a]['DeltaM'] - fits_LCDM[b]['DeltaM']
    obs_tensions[(a,b)] = diff
    print(f"  DeltaM({a}) - DeltaM({b}) = {diff:+.4f} mag")

def Delta_mu_LCDM_minus_STAM(z):
    return mu_predict(z, h2_LCDM) - mu_predict(z, h2_STAM)

def weighted_mean_Delta_mu(z, mu_err):
    w = 1.0/mu_err**2
    return np.sum(w*Delta_mu_LCDM_minus_STAM(z))/np.sum(w)

mean_dmu = {
    'Pantheon+': weighted_mean_Delta_mu(zP, muPe),
    'Union3'   : weighted_mean_Delta_mu(zU, muUe),
    'DES-Y5'   : weighted_mean_Delta_mu(zD, muDe),
}

print()
print("STAM-predicted inter-catalog DeltaM tension:")
pred_tensions = {}
for a,b in pairs:
    pred = -(mean_dmu[a] - mean_dmu[b])
    pred_tensions[(a,b)] = pred
    print(f"  pred[DeltaM({a}) - DeltaM({b})] = {pred:+.4f} mag")

print()
print("=" * 72)
print(f"HEAD-TO-HEAD (DES-Y5 minus top-{N_OUTLIERS} outliers):")
print("=" * 72)
print(f"{'pair':<28} {'observed':>10} {'predicted':>10} {'ratio':>8} {'sign':>8}")
for a,b in pairs:
    obs = obs_tensions[(a,b)]
    pred = pred_tensions[(a,b)]
    ratio = obs/pred if abs(pred)>1e-6 else float('nan')
    sign = "yes" if obs*pred>0 else "NO"
    print(f"  {a:<10} - {b:<10}    {obs:>+10.4f} {pred:>+10.4f} {ratio:>8.3f} {sign:>8}")

# ==================================================================
# Plotting (same layout as 39 for direct comparison)
# ==================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# A: Per-catalog DeltaM
ax = axes[0,0]
cats = ['Pantheon+','Union3','DES-Y5']
xs = np.arange(len(cats))
LCDM_dM = [fits_LCDM[c]['DeltaM'] for c in cats]
STAM_dM = [fits_STAM[c]['DeltaM'] for c in cats]
ax.bar(xs-0.2, LCDM_dM, 0.4, label='LCDM', color='k')
ax.bar(xs+0.2, STAM_dM, 0.4, label='STAM', color='C1')
ax.axhline(np.mean(LCDM_dM), ls=':', color='gray', alpha=0.5,
           label=f'mean LCDM DeltaM = {np.mean(LCDM_dM):.3f}')
ax.set_xticks(xs); ax.set_xticklabels(cats)
ax.set_ylabel('best-fit DeltaM (mag)')
ax.set_title(f'Per-catalog DeltaM (DES top-{N_OUTLIERS} removed)')
ax.legend(loc='best')
ax.grid(alpha=0.3, axis='y')

# B: cross-cat tension comparison
ax = axes[0,1]
labels = [f'{a}\n vs {b}' for a,b in pairs]
obs_vals = [obs_tensions[p] for p in pairs]
pred_vals = [pred_tensions[p] for p in pairs]
xs2 = np.arange(len(pairs))
ax.bar(xs2-0.2, obs_vals, 0.4, label='Observed (LCDM-fit)', color='C0')
ax.bar(xs2+0.2, pred_vals, 0.4, label='STAM-predicted', color='C2')
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xticks(xs2); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel('DeltaM tension (mag)')
ax.set_title('Cross-catalog tension (observed vs STAM)')
ax.legend(loc='best')
ax.grid(alpha=0.3, axis='y')

# C: residual histogram for DES (full vs trimmed)
ax = axes[1,0]
bins = np.linspace(-3.5, 3.5, 80)
ax.hist(fit_DES_full['resid'], bins=bins, alpha=0.45, label='DES full',
        color='gray')
fit_DES_trim = fits_LCDM['DES-Y5']
ax.hist(fit_DES_trim['resid'], bins=bins, alpha=0.7,
        label=f'DES (top-{N_OUTLIERS} removed)', color='C2')
ax.axvline(0, color='k', alpha=0.4)
ax.set_yscale('log')
ax.set_xlabel('LCDM-fit residual (mag)')
ax.set_ylabel('count')
ax.set_title('DES residual distribution')
ax.legend()
ax.grid(alpha=0.3)

# D: Delta_mu(z) with each catalog mean (after trim)
ax = axes[1,1]
zg = np.linspace(0.001, 2.0, 400)
dmug = Delta_mu_LCDM_minus_STAM(zg)
ax.plot(zg, dmug, lw=2.5, color='C0', label='Delta_mu(z) = LCDM - STAM')
medians = {'Pantheon+': np.median(zP), 'Union3': np.median(zU), 'DES-Y5': np.median(zD)}
for cat, color in zip(cats, ['C0','C3','C2']):
    ax.scatter([medians[cat]], [mean_dmu[cat]], s=120, color=color, zorder=10,
               label=f'{cat}: <Delta_mu> = {mean_dmu[cat]:.3f}')
ax.axhline(0, color='gray', alpha=0.5)
ax.set_xlabel('redshift z')
ax.set_ylabel('Delta_mu (mag)')
ax.set_title('STAM Delta_mu(z) and per-catalog means')
ax.legend()
ax.grid(alpha=0.3)

plt.suptitle(
    f"Cross-catalog tension after removing {N_OUTLIERS} DES outliers",
    fontsize=13
)
plt.tight_layout(rect=[0, 0, 1, 0.96])

outdir = Path('reports/script_39b')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir/'catalog_tension_DES_trimmed.png', dpi=130, bbox_inches='tight')
plt.savefig(outdir/'catalog_tension_DES_trimmed.pdf', bbox_inches='tight')

# Save numerical summary
summary = pd.DataFrame([
    {'pair': f'{a} vs {b}',
     'observed_mag': obs_tensions[(a,b)],
     'predicted_mag': pred_tensions[(a,b)],
     'ratio_obs_pred': obs_tensions[(a,b)]/pred_tensions[(a,b)]
                       if abs(pred_tensions[(a,b)])>1e-6 else float('nan'),
     'sign_match': obs_tensions[(a,b)]*pred_tensions[(a,b)]>0}
    for a,b in pairs
])
summary.to_csv(outdir/'cross_catalog_DES_trimmed_summary.csv', index=False)

# Save list of removed outliers
out_df = pd.DataFrame({
    'CID': [des_cosmo['CID'].iloc[i] for i in outlier_idx],
    'z': zD_full[outlier_idx],
    'MU': muD_full[outlier_idx],
    'MUERR': muDe_full[outlier_idx],
    'residual_LCDM': fit_DES_full['resid'][outlier_idx],
    'abs_residual_LCDM': abs_resid[outlier_idx],
})
out_df.to_csv(outdir/'removed_outliers.csv', index=False)

print(f"\nSaved: {outdir.resolve()}")
print("=" * 72)
