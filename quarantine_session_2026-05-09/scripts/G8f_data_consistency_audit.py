"""
G8f - Data consistency audit at z = 0.28 - 0.32

Don't compute residuals. Look at the catalogs' own per-SN quality and
spread metrics across z slices and check whether the data itself shows
unusual inconsistency in the z = 0.28 - 0.32 range relative to
neighboring z bins.

Quantities sampled per slice:
  Pantheon+:
    - m_b_corr_err_DIAG  (per-SN total error)
    - m_b_corr_err_RAW   (per-SN error before bias correction)
    - FITCHI2 / NDOF     (light-curve fit quality)
    - FITPROB            (fit probability)
    - x1, x1ERR          (SALT2 stretch + uncertainty)
    - c, cERR            (SALT2 color + uncertainty)
    - biasCor_m_b        (selection bias correction)
    - biasCorErr_m_b     (uncertainty on the correction)
    - biasCor_m_b_COVSCALE (cov scaling factor)
    - HOST_LOGMASS_ERR
  DES:
    - MUERR              (total per-SN error)
    - MUERR_RAW, MUERR_VPEC, MUERR_RENORM
    - MUPULL             (DES's standardized residual)
    - FITCHI2 / NDOF
    - FITPROB
    - x1, x1ERR, c, cERR
    - biasCor_mu, biasCorErr_mu, biasCor_muCOVSCALE, biasCor_muCOVADD
    - PROBCC_BEAMS

For each slice we report median, std, and (where appropriate) the
fraction of "anomalous" SNe (e.g. FITPROB < 0.01, |MUPULL| > 3).

Author: STAM Model-A research line, 2026-05-09.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from pathlib import Path

print("=" * 78)
print("G8f - Data consistency audit at z = 0.30 - 0.35")
print("=" * 78)

# Bin edges: standard 0.1 wide; narrow special bin at 0.28-0.32
bin_edges = [0.00, 0.10, 0.20, 0.28, 0.32, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
bin_labels = []
for i in range(len(bin_edges) - 1):
    lo, hi = bin_edges[i], bin_edges[i + 1]
    if lo == 0.28 and hi == 0.32:
        bin_labels.append(f"{lo:.2f}-{hi:.2f} *")
    else:
        bin_labels.append(f"{lo:.2f}-{hi:.2f}")


def slice_summary(df, z_col, columns, edges, robust=True):
    """For each (lo,hi) z slice, report N, median, std for each column.
    Robust=True uses median + IQR-based scale; False uses mean + std.
    """
    rows = []
    z = df[z_col].values
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = (z >= lo) & (z < hi)
        N = int(m.sum())
        row = {'z_lo': lo, 'z_hi': hi, 'label': bin_labels[i], 'N': N}
        if N >= 5:
            sub = df.loc[m, columns]
            for col in columns:
                vals = sub[col].values
                vals = vals[np.isfinite(vals)]
                if len(vals) >= 5:
                    if robust:
                        row[f'{col}_med'] = float(np.median(vals))
                        # robust std = MAD * 1.4826
                        mad = np.median(np.abs(vals - np.median(vals)))
                        row[f'{col}_madstd'] = float(mad * 1.4826)
                    else:
                        row[f'{col}_med'] = float(np.mean(vals))
                        row[f'{col}_madstd'] = float(np.std(vals))
                else:
                    row[f'{col}_med'] = np.nan
                    row[f'{col}_madstd'] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def fraction_above_threshold(df, z_col, col, edges, threshold,
                              comparison='abs_gt'):
    rows = []
    z = df[z_col].values
    vals = df[col].values
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = (z >= lo) & (z < hi)
        N = int(m.sum())
        if N >= 5:
            v = vals[m]
            v = v[np.isfinite(v)]
            if comparison == 'abs_gt':
                n_anom = int(np.sum(np.abs(v) > threshold))
            elif comparison == 'lt':
                n_anom = int(np.sum(v < threshold))
            else:
                n_anom = int(np.sum(v > threshold))
            frac = n_anom / max(len(v), 1)
        else:
            n_anom = 0
            frac = np.nan
        rows.append({'z_lo': lo, 'z_hi': hi, 'label': bin_labels[i],
                     'N': N, 'n_anom': n_anom, 'fraction_anomalous': frac})
    return pd.DataFrame(rows)


# ===================================================================
# Pantheon+ audit
# ===================================================================
print("\nLoading Pantheon+ ...")
pantheon = pd.read_csv('data/pantheon.csv')
mP = (pantheon['IS_CALIBRATOR'] == 0) & (pantheon['zCMB'] > 0.01)
pantheon_cosmo = pantheon[mP].copy()
print(f"  N = {len(pantheon_cosmo)}")

pcols = [
    'm_b_corr_err_DIAG', 'm_b_corr_err_RAW',
    'x1', 'x1ERR', 'c', 'cERR',
    'biasCor_m_b', 'biasCorErr_m_b', 'biasCor_m_b_COVSCALE',
    'HOST_LOGMASS_ERR',
]
# Compute FITCHI2/NDOF
pantheon_cosmo['FITCHI2_per_NDOF'] = (
    pantheon_cosmo['FITCHI2'] / pantheon_cosmo['NDOF'].replace(0, np.nan)
)
pcols_audit = pcols + ['FITCHI2_per_NDOF', 'FITPROB']
pantheon_summary = slice_summary(pantheon_cosmo, 'zCMB', pcols_audit, bin_edges)

print("\n=== Pantheon+: median & MAD-std per column per z slice ===")
# Print compactly
print(f"\n{'bin':<18} {'N':>5}", end='')
for c in pcols_audit:
    short = c[:14]
    print(f" {short:>14}", end='')
print()
for _, row in pantheon_summary.iterrows():
    if row['N'] >= 5:
        print(f"{row['label']:<18} {int(row['N']):>5}", end='')
        for c in pcols_audit:
            mv = row.get(f'{c}_med', np.nan)
            sv = row.get(f'{c}_madstd', np.nan)
            print(f" {mv:>+6.3f}±{sv:>5.3f}", end='')
        print()
    else:
        print(f"{row['label']:<18} {int(row['N']):>5}  (too few)")

# Anomaly fractions
print("\n=== Pantheon+: anomaly fractions per z slice ===")
print(f"{'bin':<18} {'N':>5} {'FITPROB<0.01':>14} "
      f"{'|x1|>3':>10} {'|c|>0.3':>10}")
af_fitprob = fraction_above_threshold(pantheon_cosmo, 'zCMB', 'FITPROB',
                                       bin_edges, 0.01, 'lt')
af_x1 = fraction_above_threshold(pantheon_cosmo, 'zCMB', 'x1',
                                  bin_edges, 3.0, 'abs_gt')
af_c = fraction_above_threshold(pantheon_cosmo, 'zCMB', 'c',
                                 bin_edges, 0.3, 'abs_gt')
for i in range(len(bin_edges) - 1):
    if af_fitprob.iloc[i]['N'] >= 5:
        print(f"{bin_labels[i]:<18} {int(af_fitprob.iloc[i]['N']):>5d}  "
              f"{af_fitprob.iloc[i]['fraction_anomalous']:>10.3%}  "
              f"{af_x1.iloc[i]['fraction_anomalous']:>10.3%}  "
              f"{af_c.iloc[i]['fraction_anomalous']:>10.3%}")

# ===================================================================
# DES audit
# ===================================================================
print("\nLoading DES ...")
des = pd.read_csv('data/des.csv')
mD = (des['zCMB'] > 0.01) & (des['MUERR'] > 0) & np.isfinite(des['MU'])
des_cosmo = des[mD].copy()
print(f"  N = {len(des_cosmo)}")

dcols = [
    'MUERR', 'MUERR_RAW', 'MUERR_VPEC', 'MUERR_RENORM',
    'MUPULL',
    'x1', 'x1ERR', 'c', 'cERR',
    'biasCor_mu', 'biasCorErr_mu', 'biasCor_muCOVSCALE',
    'PROBCC_BEAMS', 'LENSDMU',
]
des_cosmo['FITCHI2_per_NDOF'] = (
    des_cosmo['FITCHI2'] / des_cosmo['NDOF'].replace(0, np.nan)
)
dcols_audit = dcols + ['FITCHI2_per_NDOF', 'FITPROB']
des_summary = slice_summary(des_cosmo, 'zCMB', dcols_audit, bin_edges)

print("\n=== DES: median & MAD-std per column per z slice ===")
# Two passes for readability
group_a = ['MUERR', 'MUERR_RAW', 'MUERR_VPEC', 'MUERR_RENORM', 'MUPULL']
group_b = ['x1', 'c', 'biasCor_mu', 'PROBCC_BEAMS', 'LENSDMU',
           'FITCHI2_per_NDOF', 'FITPROB']

for grp_name, group in [('Errors + MUPULL', group_a),
                         ('Light-curve / host / quality', group_b)]:
    print(f"\n  ({grp_name})")
    print(f"  {'bin':<18} {'N':>5}", end='')
    for c in group:
        short = c[:13]
        print(f" {short:>13}", end='')
    print()
    for _, row in des_summary.iterrows():
        if row['N'] >= 5:
            print(f"  {row['label']:<18} {int(row['N']):>5}", end='')
            for c in group:
                mv = row.get(f'{c}_med', np.nan)
                sv = row.get(f'{c}_madstd', np.nan)
                if np.isfinite(mv):
                    print(f" {mv:>+5.2f}±{sv:>5.2f}", end='')
                else:
                    print(f"     n/a       ", end='')
            print()
        else:
            print(f"  {row['label']:<18} {int(row['N']):>5}  (too few)")

# DES anomaly fractions
print("\n=== DES: anomaly fractions per z slice ===")
print(f"{'bin':<18} {'N':>5} {'|MUPULL|>3':>12} "
      f"{'FITPROB<0.01':>14} {'PROBCC>0.05':>13}")
af_mupull = fraction_above_threshold(des_cosmo, 'zCMB', 'MUPULL',
                                      bin_edges, 3.0, 'abs_gt')
af_dfitprob = fraction_above_threshold(des_cosmo, 'zCMB', 'FITPROB',
                                        bin_edges, 0.01, 'lt')
af_probcc = fraction_above_threshold(des_cosmo, 'zCMB', 'PROBCC_BEAMS',
                                      bin_edges, 0.05, 'gt')
for i in range(len(bin_edges) - 1):
    if af_mupull.iloc[i]['N'] >= 5:
        print(f"{bin_labels[i]:<18} {int(af_mupull.iloc[i]['N']):>5d}  "
              f"{af_mupull.iloc[i]['fraction_anomalous']:>10.3%}    "
              f"{af_dfitprob.iloc[i]['fraction_anomalous']:>10.3%}    "
              f"{af_probcc.iloc[i]['fraction_anomalous']:>10.3%}")

# ===================================================================
# Save audit CSVs
# ===================================================================
outdir = Path('reports/G8f')
outdir.mkdir(parents=True, exist_ok=True)
pantheon_summary.to_csv(outdir / 'pantheon_audit_per_slice.csv', index=False)
des_summary.to_csv(outdir / 'des_audit_per_slice.csv', index=False)

# Anomaly tables combined
anom_P = pd.DataFrame({
    'bin': bin_labels,
    'N': af_fitprob['N'],
    'frac_FITPROB_lt_0.01': af_fitprob['fraction_anomalous'],
    'frac_abs_x1_gt_3': af_x1['fraction_anomalous'],
    'frac_abs_c_gt_0.3': af_c['fraction_anomalous'],
})
anom_P.to_csv(outdir / 'pantheon_anomaly_fractions.csv', index=False)

anom_D = pd.DataFrame({
    'bin': bin_labels,
    'N': af_mupull['N'],
    'frac_abs_MUPULL_gt_3': af_mupull['fraction_anomalous'],
    'frac_FITPROB_lt_0.01': af_dfitprob['fraction_anomalous'],
    'frac_PROBCC_gt_0.05': af_probcc['fraction_anomalous'],
})
anom_D.to_csv(outdir / 'des_anomaly_fractions.csv', index=False)

print(f"\nFiles saved to {outdir.resolve()}")
print("=" * 78)
