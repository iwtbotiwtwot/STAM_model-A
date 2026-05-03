# Supernova b-sensitivity audit

This audit separates two questions:

1. Are Union3/Pantheon and DES discrepant in their `z -> Hubble ordinate` relation?
2. Can a slight adjustment in STAM Model-A's `b` bring Union3, Pantheon, and DES closer?

## Data columns used

The common Hubble-diagram ordinate is:

```text
Union3:   mb
Pantheon: MU_SH0ES
DES:      MU
```

Raw `mB` is not used for the main cross-catalog comparison because raw apparent magnitudes require nuisance corrections and absolute-magnitude normalization. The selected columns are distance-modulus-like corrected Hubble-diagram values.

## Locked baseline

```text
b = 354.95
```

STAM Model-A distance law:

```text
D_adj(z) = L z (1 + 0.5 z) + b z
```

A single global reference offset is calibrated from Pantheon and Union3 only:

```text
offset_ref(b) = catalog-balanced median[y_reference - mu_STAM(z; b)]
```

DES is then evaluated as a holdout.

## Fair b-sensitivity test

For each candidate `b`, recompute the Pantheon/Union3 reference offset and then evaluate all three catalogs.

This is the fair test because the model's nuisance/intercept normalization must remain tied to the reference family.

## Fixed-offset diagnostic

The script also reports a fixed-offset diagnostic where the `b = 354.95` offset is held fixed while `b` changes. This can make DES look better for large `b`, but it pushes Pantheon/Union3 away and should not be treated as a valid cross-catalog fit.

## Command

```powershell
python scripts\06_b_sensitivity_supernova.py `
  --union3 data\union3_bins.csv `
  --pantheon data\pantheon.csv `
  --des data\des.csv `
  --b 354.95 `
  --clean-output
```

Outputs are written to:

```text
results/b_sensitivity_supernova/
```

Key outputs:

```text
summary.json
catalog_residuals_locked_b.csv
binned_empirical_discrepancy.csv
des_pointwise_empirical_vs_stam_residuals.csv
b_sensitivity_reference_offset_fine_grid.csv
b_sensitivity_fixed_offset_diagnostic_fine_grid.csv
candidate_b_table_reference_offset.csv
des_top200_centered_outliers_locked_b.csv
des_top200_absolute_residuals_locked_b.csv
```

## Interpretation rule

If DES remains high after reference calibration over a broad range of `b`, then `b` is not the lever that solves DES.

If a large `b` only helps DES under the fixed-offset diagnostic while worsening Pantheon/Union3, that is evidence against treating DES as a simple retuning target.
