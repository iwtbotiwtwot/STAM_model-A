# G100 — Hubble bias full distance-modulus target

## Purpose

G100 uses the G98/G99 locked convention but targets the full distance-modulus fit rather than a crude low-z slope. It asks: for true intrinsic `H0=73.04`, what `f_los` makes a biased photon-distance curve look like a no-bias `H0=67.4` curve over chosen redshift ranges?

## Locked convention

```text
D_C_obs(z) = D_C_intrinsic(z) + f_los * b * shape(z)
D_L_obs(z) = (1+z) * D_C_obs(z)
```

## Constants

- `A0 = 0.026525823849`
- `H0_true = 73.0400` km/s/Mpc
- `H0_target = 67.4000` km/s/Mpc
- `Omega_m = 0.3000`
- `b = 108.875163` Mpc = `355.103289` Mly

## Results

| Redshift range | Shape | H0 fit at f=0 | H0 fit at f=1 | f_los target | H0 target fit |
|---|---|---:|---:|---:|---:|
| low_SN_0p01_0p15 | linear | 73.040 | 71.118 | 3.097 | 67.400 |
| low_SN_0p01_0p15 | stam_su | 73.040 | 71.096 | 3.060 | 67.400 |
| low_SN_0p01_0p15 | saturating | 73.040 | 71.255 | 3.341 | 67.400 |
| pantheon_like_0p01_2p0 | linear | 73.040 | 70.593 | 2.415 | 67.400 |
| pantheon_like_0p01_2p0 | stam_su | 73.040 | 70.210 | 2.078 | 67.400 |
| pantheon_like_0p01_2p0 | saturating | 73.040 | 71.740 | 4.620 | 67.400 |
| low_mid_0p01_0p8 | linear | 73.040 | 70.953 | 2.845 | 67.400 |
| low_mid_0p01_0p8 | stam_su | 73.040 | 70.826 | 2.678 | 67.400 |
| low_mid_0p01_0p8 | saturating | 73.040 | 71.516 | 3.929 | 67.400 |
| high_z_0p1_2p0 | linear | 73.040 | 70.567 | 2.389 | 67.400 |
| high_z_0p1_2p0 | stam_su | 73.040 | 70.167 | 2.046 | 67.400 |
| high_z_0p1_2p0 | saturating | 73.040 | 71.764 | 4.709 | 67.400 |

## Interpretation

The full distance-modulus fit correctly recovers `H0_true` when `f_los=0`. Nonzero photon-A bias pushes the inferred no-bias H0 lower. The required `f_los` depends on redshift range and bias-shape model, so a real test must use the actual SN/BAO/CMB/chronometer likelihoods rather than a single low-z estimate.

## Status

G100 is still a scaffold, not a final cosmology result. It locks the distance convention and identifies the line-of-sight amplification scale required to turn local `H0≈73` into distance-inferred `H0≈67.4` under different bias shapes and redshift ranges. The next real test is to plug in actual SN and BAO/CMB data with covariance and fixed priors.

## Files

- Grid CSV: `/mnt/data/G100_run/results/G100_Hubble_bias_mu_fit_target_grid.csv`
- Plot: `/mnt/data/G100_run/plots/G100_required_flos_by_zrange.png`
- Plot: `/mnt/data/G100_run/plots/G100_H0_mu_fit_vs_flos.png`
- Plot: `/mnt/data/G100_run/plots/G100_best_fit_residuals.png`
