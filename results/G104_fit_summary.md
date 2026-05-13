# G104 — SN + BAO + CMB + chronometer pipeline

- data mode: `synthetic`
- bias shape: `linear`
- f_los model: `one_plus_z`
- free r_d: `False`
- fit SN Mcal offset: `True`
- probes used: `SN,BAO,CMB,HZ`

## Best-fit parameters

| Model | H0 | Omega_m | f0 | f1 | r_d | Mcal | chi2 | N | chi2/N | AIC | BIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| lcdm | 67.4169 | 0.2998 | 0.0000 | 0.0000 | 147.100 | 0.0016 | 0.540 | 57 | 0.009 | 6.540 | 12.669 |
| stam | 67.3898 | 0.3000 | -0.0198 | 0.0199 | 147.100 | 0.0018 | 0.540 | 57 | 0.009 | 10.540 | 20.755 |

## Probe chi2 breakdown

| Model | SN | BAO | CMB | HZ | Total |
|---|---:|---:|---:|---:|---:|
| lcdm | 0.495 | 0.000 | 0.001 | 0.044 | 0.540 |
| stam | 0.495 | 0.001 | 0.000 | 0.044 | 0.540 |

## Interpretation

This pipeline tests whether the STAM photon-A distance layer can improve distance probes while direct expansion probes remain pressure tests. In synthetic mode, it should recover the lower-H0 reference model unless the STAM bias is strongly preferred by the generated data. In real-data mode, compare Δχ2, AIC, BIC, and probe-level residuals.

A positive STAM result requires more than better SN residuals: it must avoid wrecking BAO radial distances, chronometers, and CMB compressed distances under one locked convention.

## Files

- Best-fit params: `/mnt/data/G104_run/results/G104_best_fit_params.csv`
- Probe breakdown: `/mnt/data/G104_run/results/G104_probe_chi2_breakdown.csv`
- Model comparison: `/mnt/data/G104_run/results/G104_model_comparison.csv`
- Plots in `plots/`
