# G99 — Hubble tension as STAM distance bias, locked convention

## Convention

G99 uses the G98-approved convention:

```text
D_C_obs(z) = D_C_intrinsic(z) + f_los * b * shape(z)
D_L_obs(z) = (1+z) * D_C_obs(z)
```

This treats photon-A traversal as a path/comoving correction first, then applies the standard luminosity-distance factor.

## Constants

- `A0 = 0.026525823849`
- `H0_true = 73.0400` km/s/Mpc
- `H0_target = 67.4000` km/s/Mpc
- `Omega_m = 0.3000`
- `b = A0*c/H0_true = 108.875163` Mpc = `355.103289` Mly
- selected shape = `linear`

## Low-z estimate

All supported shapes are linear at very low z, so the leading low-z estimate is:

```text
H0_inferred ≈ H0_true / (1 + f_los*A0)
f_los_required = (H0_true/H0_target - 1)/A0
```

- `f_los_required = 3.154644`

## Selected-shape results

| f_los | H0 inferred low-z slope | H0 inferred full mu-fit z<=2.0 |
|---:|---:|---:|
| 0.000000 | 69.1200 | 73.0400 |
| 1.000000 | 67.3033 | 70.5927 |
| 3.154644 | 63.6962 | 65.8442 |

## Interpretation

The distance-bias layer can make a true local `H0_true` appear lower in distance-only fits. Direct H(z) probes are not changed by the photon-distance bias layer. This is the core STAM interpretation of the Hubble tension as a possible distance-inference tension.

This is still a scaffold, not a publication-grade cosmology likelihood. The next serious test must include SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, nuisance parameters, and locked priors.

## Files

- Grid CSV: `/mnt/data/G99_run/results/G99_Hubble_bias_locked_grid.csv`
- Plot: `/mnt/data/G99_run/plots/G99_inferred_H0_distance_vs_flos.png`
- Plot: `/mnt/data/G99_run/plots/G99_mu_residuals_locked_convention.png`
- Plot: `/mnt/data/G99_run/plots/G99_probe_separation.png`
- Plot: `/mnt/data/G99_run/plots/G99_bias_shapes.png`
