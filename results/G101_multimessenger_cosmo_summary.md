# G101 — compressed SN+BAO+CMB+chronometer stress test

This is a compressed toy likelihood, not a publication-grade cosmology result. It uses synthetic reference data generated from no-bias LCDM with lower `H0_ref`, then tests whether a STAM model with intrinsic `H0_true` plus photon-A distance bias can mimic it across probes.

## Locked convention

```text
D_C_obs(z) = D_C_intrinsic(z; H0_true) + f_los * b * shape(z)
D_L_obs(z) = (1+z) * D_C_obs(z)
Direct H(z) probes are not distance-biased.
```

- `A0 = 0.026525823849`
- `H0_true = 73.0400`
- `H0_ref = 67.4000`
- `b = 108.875163` Mpc = `355.103289` Mly

## Results

| Shape | Fit set | f_los best | chi2/N total | SN | BAO_DM | BAO_DH | CMB_D | H_z |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| linear | SN_only | 2.366 | 1.13 | 0.06 | 0.51 | 6.63 | 59.67 | 1.43 |
| linear | CMB_only | 1.577 | 0.85 | 0.37 | 1.95 | 6.63 | 0.00 | 1.43 |
| linear | distance_only_SN_BAO_CMB | 1.856 | 0.73 | 0.19 | 1.06 | 6.63 | 7.44 | 1.43 |
| linear | distance_plus_radial_BAO | 1.856 | 0.73 | 0.19 | 1.06 | 6.63 | 7.44 | 1.43 |
| linear | all_SN_BAO_CMB_Hz | 1.856 | 0.73 | 0.19 | 1.06 | 6.63 | 7.44 | 1.43 |
| stam_su | SN_only | 1.986 | 102124.90 | 0.13 | 1.23 | 6.63 | 10110303.71 | 1.43 |
| stam_su | CMB_only | 0.010 | 3.63 | 3.02 | 14.77 | 6.63 | 0.00 | 1.43 |
| stam_su | distance_only_SN_BAO_CMB | 0.010 | 3.63 | 3.02 | 14.76 | 6.63 | 0.01 | 1.43 |
| stam_su | distance_plus_radial_BAO | 0.010 | 3.63 | 3.02 | 14.76 | 6.63 | 0.01 | 1.43 |
| stam_su | all_SN_BAO_CMB_Hz | 0.010 | 3.63 | 3.02 | 14.76 | 6.63 | 0.01 | 1.43 |
| saturating | SN_only | 4.493 | 2.93 | 0.08 | 0.35 | 6.63 | 237.26 | 1.43 |
| saturating | CMB_only | 20.000 | 31.93 | 25.55 | 167.94 | 6.63 | 232.99 | 1.43 |
| saturating | distance_only_SN_BAO_CMB | 4.517 | 2.93 | 0.08 | 0.35 | 6.63 | 237.25 | 1.43 |
| saturating | distance_plus_radial_BAO | 4.517 | 2.93 | 0.08 | 0.35 | 6.63 | 237.25 | 1.43 |
| saturating | all_SN_BAO_CMB_Hz | 4.517 | 2.93 | 0.08 | 0.35 | 6.63 | 237.25 | 1.43 |
| log | SN_only | 3.446 | 2.79 | 0.00 | 0.00 | 6.63 | 231.86 | 1.43 |
| log | CMB_only | 20.000 | 57.55 | 46.90 | 339.90 | 6.63 | 201.24 | 1.43 |
| log | distance_only_SN_BAO_CMB | 3.488 | 2.79 | 0.00 | 0.00 | 6.63 | 231.78 | 1.43 |
| log | distance_plus_radial_BAO | 3.488 | 2.79 | 0.00 | 0.00 | 6.63 | 231.78 | 1.43 |
| log | all_SN_BAO_CMB_Hz | 3.488 | 2.79 | 0.00 | 0.00 | 6.63 | 231.78 | 1.43 |

## Interpretation

Distance-only probes can be shifted by `f_los`; direct expansion probes cannot. If a distance-bias shape improves SN/CMB transverse distances while BAO radial and chronometer H(z) remain badly off, the model is not yet a full Hubble-tension solution. It is a distance-bias mechanism that still needs a realistic multi-probe likelihood and likely a structure-dependent `f_los(z, sightline)` rather than one constant.

## Files

- CSV: `/mnt/data/G101_run/results/G101_multimessenger_cosmo_grid.csv`
- Plot: `/mnt/data/G101_run/plots/G101_chi2_by_probe.png`
- Plot: `/mnt/data/G101_run/plots/G101_best_residuals.png`
- Plot: `/mnt/data/G101_run/plots/G101_shape_comparison.png`
