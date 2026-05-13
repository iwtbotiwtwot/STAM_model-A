# G107 — Identify L_mu and D_mu for decoherence channels

## Structural formula

```text
Gamma_res = sum_mu <L_mu^dagger L_mu> D_mu
```

STAM interpretation:

- `L_mu` is a physical interaction / write channel.
- `D_mu` is the distinguishability created by that channel.
- `Gamma_res` is the proper-time SU-write intensity.
- each write resolves `1 SU = A0 = 1/(12π)`.

`A0 = 0.026525823849`

## Chosen example

- separation Δx = `1e-06` m
- single photon wavelength λ = `5e-07` m
- single-channel rate γ = `1000` s^-1
- exposure time t = `0.01` s

## Results

| case | L_mu | D_mu | gamma | D | Gamma_res | P_resolved | N_write | A_ledger |
|---|---|---|---:|---:|---:|---:|---:|---:|
| single_photon_scattering | `sqrt(gamma_k) exp(i k x_hat)` | `1 - sinc(k Delta_x)` | 1.0000e+03 | 1.0000e+00 | 1.0000e+03 | 9.9995e-01 | 1.0000e+01 | 2.6526e-01 |
| thermal_photon_bath | `sqrt(gamma(lambda)) exp(i k(lambda) x_hat)` | `integral[1 - sinc(k Delta_x)]` | 1.0000e+03 | 6.1841e-02 | 6.1841e+01 | 4.6120e-01 | 6.1841e-01 | 1.6404e-02 |
| gas_collision | `sqrt(Gamma_q) exp(i q x_hat / hbar)` | `1 - cos(q Delta_x / hbar)` | 1.0000e+03 | 1.5611e+00 | 1.5611e+03 | 1.0000e+00 | 1.5611e+01 | 4.1409e-01 |
| detector_absorption | `sqrt(gamma_i) |i><i|` | `approximately 1` | 1.0000e+03 | 1.0000e+00 | 1.0000e+03 | 9.9995e-01 | 1.0000e+01 | 2.6526e-01 |

## Interpretation

This is the first concrete identification of `L_mu` and `D_mu` for ordinary environmental decoherence. Photon scattering uses `L_k = sqrt(gamma_k) exp(i k x_hat)` and `D_k = 1 - sinc(k Delta_x)`. Gas collision uses `L_q = sqrt(Gamma_q) exp(i q x_hat / hbar)` and a momentum-transfer distinguishability. Detector absorption is the high-D limit with macroscopically distinguishable records.

In STAM language, standard decoherence rates become SU-write rates: `dN_write = Gamma_res dt`, `dA_ledger = A0 Gamma_res dt`.

## Files

- CSV: `/mnt/data/G107_decoherence_run/results/G107_Lmu_Dmu_decoherence_grid.csv`
- Plot: `/mnt/data/G107_decoherence_run/plots/G107_distinguishability_vs_separation.png`
- Plot: `/mnt/data/G107_decoherence_run/plots/G107_resolution_probability.png`
- Plot: `/mnt/data/G107_decoherence_run/plots/G107_write_count_vs_time.png`
- Plot: `/mnt/data/G107_decoherence_run/plots/G107_gamma_vs_wavelength.png`
