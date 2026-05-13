# G98 — Distance-bias units / definition audit

## Why G98 exists

G97 exposed a bookkeeping mismatch: the analytic low-z estimate gave `f_los_required ≈ 3.15`, while the numerical printout appeared to give `f_los ≈ 1`.

## Locked constants

- `A0 = 0.026525823849`
- `H0_true = 73.0400` km/s/Mpc
- `H0_target = 67.4000` km/s/Mpc
- `b = A0*c/H0_true = 108.875163` Mpc = `355.103289` Mly

## Correct low-z relation

If the photon-A bias adds to the same low-z fitted distance slope,

```text
D_obs ≈ (c/H0_true) z + f_los * (A0*c/H0_true) z
      ≈ (c/H0_true) z * (1 + f_los*A0)

H0_inferred ≈ H0_true / (1 + f_los*A0)
f_los_required = (H0_true/H0_target - 1) / A0
```

- `f_los_required = 3.154644`
- `H0(f=1) = 71.152618` km/s/Mpc
- `H0(f_required) = 67.400000` km/s/Mpc

## Main conclusion

The analytic value `f_los ≈ 3.15` is the correct low-z requirement if the bias term is `f_los*b*z` applied to the low-z fitted distance slope. If a script appears to return `f_los≈1` for the same target, it is mixing conventions or reporting the wrong row.

## Recommended convention for STAM going forward

Treat photon-A traversal as a path/comoving correction first:

```text
Intrinsic expansion distance:
    D_C(z) = ∫ c/H(z) dz

Photon-A traversal bias:
    ΔD_C(z) = f_los * b * shape(z)

Observed luminosity distance:
    D_L_obs(z) = (1+z) * [D_C(z) + ΔD_C(z)]
```

This makes the distance-bias layer explicit and keeps it separate from the expansion layer.

## Files

- CSV: `/mnt/data/G98_run/results/G98_distance_bias_units_grid.csv`
- Plot: `/mnt/data/G98_run/plots/G98_H0_inference_conventions.png`
- Plot: `/mnt/data/G98_run/plots/G98_lowz_residuals_conventions.png`
