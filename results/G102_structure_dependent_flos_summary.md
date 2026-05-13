# G102 — structure-dependent f_los test

## Purpose

G102 replaces the single global `f_los` from G101 with simple redshift-dependent forms:

```text
constant:      f_los(z)=f0
one_plus_z:    f_los(z)=f0+f1 z/(1+z)
log:           f_los(z)=f0+f1 ln(1+z)
```

## Locked convention

```text
D_C_obs(z) = D_C_intrinsic(z; H0_true) + b*f_los(z)*shape(z)
D_L_obs(z) = (1+z)*D_C_obs(z)
Direct H(z) probes are not distance-biased.
```

## Constants

- `A0 = 0.026525823849`
- `H0_true = 73.0400`
- `H0_ref = 67.4000`
- `b = 108.875163` Mpc = `355.103289` Mly

## Best distance-only fit

- shape: `linear`
- f_los model: `one_plus_z`
- params: `3.2;-1.65`
- chi2/N total: `0.459`
- SN: `0.006`
- BAO_DM: `0.073`
- CMB_D: `0.064`

## Best all-probes fit

- shape: `linear`
- f_los model: `one_plus_z`
- params: `3.2;-1.65`
- chi2/N total: `0.459`
- SN: `0.006`
- BAO_DM: `0.073`
- BAO_DH: `6.625`
- CMB_D: `0.064`
- H_z: `1.429`

## Interpretation

If a structure-dependent `f_los(z)` improves SN/BAO_DM/CMB_D but BAO_DH and H_z remain strained, this supports the interpretation that the mechanism is a distance-bias layer, not an expansion-rate fix.

This remains a toy compressed likelihood. A real test must use actual SN covariance, BAO likelihoods, CMB acoustic angle, chronometers, structure growth, and locked priors.

## Files

- CSV: `/mnt/data/G102_fast_run/results/G102_structure_dependent_flos_grid.csv`
- Plot: `/mnt/data/G102_fast_run/plots/G102_chi2_comparison.png`
- Plot: `/mnt/data/G102_fast_run/plots/G102_best_flos_of_z.png`
- Plot: `/mnt/data/G102_fast_run/plots/G102_best_residuals.png`
