"""
G8e - Continuity-region scatter test (Pantheon+ and DES, broad bins)

Looks for excess scatter in per-SN residuals across coarse z bins, with
a narrow special bin at z = 0.28 - 0.32 covering the suspected matter-A
transition region. If the framework's continuity prediction holds, the
transition bin should show enhanced rms residual relative to neighboring
control bins, beyond what reported per-SN errors predict.

Bin scheme (mostly delta z = 0.1; transition narrowed to 0.28-0.32):
    0.0-0.1, 0.1-0.2, 0.2-0.28,
    0.28-0.32 (TRANSITION, narrow),
    0.32-0.40, 0.40-0.50,
    0.5-0.6, 0.6-0.7, 0.7-0.8, 0.8-0.9, 0.9-1.0, 1.0-1.5, 1.5-2.5

Cosmology: Model-A V_3 numerical KG (G8c).
Per-catalog ΔM marginalized analytically with diagonal weights.

Author: STAM Model-A research line, 2026-05-09.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

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
print("G8e - Continuity-region scatter test")
print("Pantheon+ and DES, V_3 numerical KG, broad z bins with 0.30-0.50 isolated")
print("=" * 72)

# Load V_3 numerical KG trajectory from G8c
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
    z_max = max(np.max(z_obs) * 1.05, 0.01)
    z_grid = np.linspace(0.0, z_max, n_grid)
    a_grid = 1.0 / (1.0 + z_grid)
    inv_h = 1.0 / np.sqrt(h2_func(a_grid))
    DC_grid = (C_KMS / H0_used) * cumtrapz0(inv_h, z_grid)
    DC_obs = np.interp(z_obs, z_grid, DC_grid)
    dL_obs = (1.0 + z_obs) * DC_obs
    return 5.0 * np.log10(dL_obs) + 25.0


def fit_DeltaM_diag(z_obs, mu_obs, sigma, h2_func, H0_used):
    mu_model = mu_predict(z_obs, h2_func, H0_used)
    w = 1.0 / sigma**2
    DM = float(np.sum(w * (mu_obs - mu_model)) / np.sum(w))
    resid = mu_obs - mu_model - DM
    return resid, DM


# ===================================================================
# Pantheon+ — STAT+SYS diagonal sigma per SN
# ===================================================================
print("\nLoading Pantheon+ ...")
pantheon = pd.read_csv('data/pantheon.csv')
mP = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
mask_idx_P = np.where(mP.values)[0]
zP = pantheon.loc[mP, 'zCMB'].values
muP = pantheon.loc[mP, 'MU_SH0ES'].values
with open('data/Pantheon+SH0ES_STAT+SYS.cov', 'r') as f:
    N_full = int(f.readline().strip())
    flat = np.fromfile(f, sep='\n', dtype=np.float64, count=N_full * N_full)
C_full = flat.reshape(N_full, N_full)
C_full = 0.5 * (C_full + C_full.T)
covP = C_full[np.ix_(mask_idx_P, mask_idx_P)]
sigmaP = np.sqrt(np.diag(covP))
print(f"  N = {len(zP)}, sigma_diag median = {np.median(sigmaP):.4f}")

# ===================================================================
# DES — diagonal MUERR per SN
# ===================================================================
print("\nLoading DES ...")
des = pd.read_csv('data/des.csv')
des_mask = (des['zCMB'] > 0.01) & (des['MUERR'] > 0) & np.isfinite(des['MU'])
zD = des.loc[des_mask, 'zCMB'].values
muD = des.loc[des_mask, 'MU'].values
sigmaD = des.loc[des_mask, 'MUERR'].values
print(f"  N = {len(zD)}, z range [{zD.min():.4f}, {zD.max():.4f}], "
      f"sigma_diag median = {np.median(sigmaD):.4f}")


# ===================================================================
# Compute per-SN residuals under V_3 numerical KG
# ===================================================================
resid_P_V3, DM_P_V3 = fit_DeltaM_diag(zP, muP, sigmaP, h2_v3_num, H0_KMS_MPC)
resid_D_V3, DM_D_V3 = fit_DeltaM_diag(zD, muD, sigmaD, h2_v3_num, H0_KMS_MPC)
resid_P_L, DM_P_L = fit_DeltaM_diag(zP, muP, sigmaP, h2_LCDM, H0_LCDM_KMS_MPC)
resid_D_L, DM_D_L = fit_DeltaM_diag(zD, muD, sigmaD, h2_LCDM, H0_LCDM_KMS_MPC)

print(f"\nResiduals (per-catalog DeltaM marginalized):")
print(f"  Pantheon+ V_3  : DM = {DM_P_V3:+.4f}, RMS = {np.sqrt(np.mean(resid_P_V3**2)):.4f}")
print(f"  DES       V_3  : DM = {DM_D_V3:+.4f}, RMS = {np.sqrt(np.mean(resid_D_V3**2)):.4f}")
print(f"  Pantheon+ LCDM : DM = {DM_P_L:+.4f}, RMS = {np.sqrt(np.mean(resid_P_L**2)):.4f}")
print(f"  DES       LCDM : DM = {DM_D_L:+.4f}, RMS = {np.sqrt(np.mean(resid_D_L**2)):.4f}")


# ===================================================================
# Bin scheme — broad bins, transition isolated
# ===================================================================
bin_edges = [0.00, 0.10, 0.20, 0.28, 0.32, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.50, 2.50]
bin_labels = []
for i in range(len(bin_edges) - 1):
    lo, hi = bin_edges[i], bin_edges[i + 1]
    if lo == 0.28 and hi == 0.32:
        bin_labels.append(f"{lo:.2f}-{hi:.2f} *TRANSITION*")
    else:
        bin_labels.append(f"{lo:.2f}-{hi:.2f}")


def bin_stats(z, resid, sigma, edges):
    rows = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = (z >= lo) & (z < hi)
        N = int(m.sum())
        if N >= 2:
            r = resid[m]
            s = sigma[m]
            rms = float(np.sqrt(np.mean(r**2)))
            expected_rms = float(np.sqrt(np.mean(s**2)))
            chi2_per_sn = float(np.mean((r / s)**2))
            sigma_chi = math.sqrt(2.0 / N)  # 1-sigma error on chi2_per_sn
            excess_sigma = (chi2_per_sn - 1.0) / sigma_chi if sigma_chi > 0 else np.nan
            mean_r = float(np.mean(r))
        else:
            rms = expected_rms = chi2_per_sn = sigma_chi = excess_sigma = mean_r = np.nan
        rows.append({
            'z_lo': lo, 'z_hi': hi, 'label': bin_labels[i], 'N': N,
            'mean_resid': mean_r, 'rms_resid': rms,
            'expected_rms': expected_rms,
            'ratio_rms_to_expected': (rms / expected_rms) if expected_rms > 0 else np.nan,
            'chi2_per_sn': chi2_per_sn, 'sigma_chi2': sigma_chi,
            'excess_sigma': excess_sigma,
        })
    return pd.DataFrame(rows)


def print_table(df, header):
    print("\n" + "=" * 72)
    print(header)
    print("=" * 72)
    print(f"{'bin':<22} {'N':>5} {'mean_r':>9} {'rms_r':>9} {'<sig>':>9} "
          f"{'rms/<sig>':>11} {'chi2/SN':>9} {'excess':>10}")
    for _, row in df.iterrows():
        if row['N'] >= 2:
            print(f"{row['label']:<22} {int(row['N']):>5d} "
                  f"{row['mean_resid']:>+9.4f} {row['rms_resid']:>9.4f} "
                  f"{row['expected_rms']:>9.4f} "
                  f"{row['ratio_rms_to_expected']:>11.4f} "
                  f"{row['chi2_per_sn']:>9.4f} "
                  f"{row['excess_sigma']:>+8.2f} sig")
        else:
            print(f"{row['label']:<22} {int(row['N']):>5d}  (too few SNe)")


tableP_V3 = bin_stats(zP, resid_P_V3, sigmaP, bin_edges)
tableD_V3 = bin_stats(zD, resid_D_V3, sigmaD, bin_edges)

# Combined Pantheon+ + DES
zPD = np.concatenate([zP, zD])
residPD_V3 = np.concatenate([resid_P_V3, resid_D_V3])
sigmaPD = np.concatenate([sigmaP, sigmaD])
tablePD_V3 = bin_stats(zPD, residPD_V3, sigmaPD, bin_edges)

# Compare LCDM to ensure we know what's catalog-intrinsic
tableP_L = bin_stats(zP, resid_P_L, sigmaP, bin_edges)
tableD_L = bin_stats(zD, resid_D_L, sigmaD, bin_edges)
residPD_L = np.concatenate([resid_P_L, resid_D_L])
tablePD_L = bin_stats(zPD, residPD_L, sigmaPD, bin_edges)

print_table(tableP_V3, "Pantheon+  under V_3 numerical KG")
print_table(tableD_V3, "DES        under V_3 numerical KG")
print_table(tablePD_V3, "Pantheon+ + DES combined  under V_3 numerical KG")

print_table(tableP_L, "Pantheon+  under LCDM (reference)")
print_table(tableD_L, "DES        under LCDM (reference)")
print_table(tablePD_L, "Pantheon+ + DES combined  under LCDM (reference)")


# ===================================================================
# Plot — bin RMS / expected as bar chart, with transition highlighted
# ===================================================================
fig, axes = plt.subplots(2, 1, figsize=(13, 9))

# Top: RMS / expected ratio per bin
ax = axes[0]
x = np.arange(len(bin_edges) - 1)
w = 0.27
mask_data = tablePD_V3['N'] >= 2
ratios_PD_V3 = tablePD_V3['ratio_rms_to_expected'].values
ratios_PD_L = tablePD_L['ratio_rms_to_expected'].values
ratios_P_V3 = tableP_V3['ratio_rms_to_expected'].values
ratios_D_V3 = tableD_V3['ratio_rms_to_expected'].values

ax.bar(x - w, ratios_P_V3, w, label='Pantheon+ V_3', color='C0')
ax.bar(x, ratios_D_V3, w, label='DES V_3', color='C2')
ax.bar(x + w, ratios_PD_V3, w, label='Combined V_3', color='C3', alpha=0.85)
ax.axhline(1.0, color='red', ls=':', alpha=0.7, label='null = 1')
# Mark transition bin (3-5 in the bin_edges sequence is 0.30-0.50, idx 3)
trans_idx = [i for i, lab in enumerate(bin_labels) if 'TRANSITION' in lab][0]
ax.axvspan(trans_idx - 0.5, trans_idx + 0.5, alpha=0.18, color='red',
           label='transition bin 0.28-0.32')
ax.set_xticks(x)
ax.set_xticklabels(bin_labels, rotation=30, ha='right')
ax.set_ylabel('rms_residual / expected_rms')
ax.set_title('Per-bin scatter ratio under V_3 numerical KG')
ax.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3, axis='y')

# Bottom: chi2 / SN per bin
ax = axes[1]
chi_P = tableP_V3['chi2_per_sn'].values
chi_D = tableD_V3['chi2_per_sn'].values
chi_PD = tablePD_V3['chi2_per_sn'].values
errs_PD = tablePD_V3['sigma_chi2'].values

ax.bar(x - w, chi_P, w, label='Pantheon+ V_3', color='C0')
ax.bar(x, chi_D, w, label='DES V_3', color='C2')
ax.bar(x + w, chi_PD, w, label='Combined V_3', color='C3', alpha=0.85)
# error bars on combined
ax.errorbar(x + w, chi_PD, yerr=errs_PD, fmt='none', color='black', alpha=0.6,
            capsize=3, label='1-sigma on combined')
ax.axhline(1.0, color='red', ls=':', alpha=0.7, label='null = 1')
ax.axvspan(trans_idx - 0.5, trans_idx + 0.5, alpha=0.18, color='red')
ax.set_xticks(x)
ax.set_xticklabels(bin_labels, rotation=30, ha='right')
ax.set_ylabel('chi^2 per SN  =  <(r/sigma)^2>')
ax.set_title('Standardized scatter per bin under V_3 numerical KG')
ax.legend(loc='upper right', fontsize=9)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
outdir = Path('reports/G8e')
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / 'continuity_scatter_broad_bins.png', dpi=140,
            bbox_inches='tight')

# Save bin-stats CSV
all_tables = pd.concat([
    tableP_V3.assign(catalog='Pantheon+', cosmology='V_3 numerical KG'),
    tableD_V3.assign(catalog='DES', cosmology='V_3 numerical KG'),
    tablePD_V3.assign(catalog='Pantheon+ + DES', cosmology='V_3 numerical KG'),
    tableP_L.assign(catalog='Pantheon+', cosmology='LCDM reference'),
    tableD_L.assign(catalog='DES', cosmology='LCDM reference'),
    tablePD_L.assign(catalog='Pantheon+ + DES', cosmology='LCDM reference'),
])
all_tables.to_csv(outdir / 'continuity_scatter_broad_bins.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 72)
