# Supernova discrepancy audit: Union3/Pantheon vs DES

This audit is designed to test whether the STAM Model-A residual against DES with locked `b = 354.95` is comparable to a catalog-level DES discrepancy relative to Union3/Pantheon.

It is **not** a DES fit.

## Default command

```powershell
python scripts\05_audit_supernova_catalog_discrepancy.py `
  --union3 data\union3_bins.csv `
  --pantheon data\pantheon.csv `
  --des data\des.csv `
  --b 354.95
```

Outputs are written to:

```text
results/supernova_discrepancy/
```

## Data columns used

The script uses a common distance-modulus-like Hubble-diagram ordinate:

```text
Union3:   mb
Pantheon: MU_SH0ES
DES:      MU
```

DES and Pantheon also contain raw/fitted `mB`, but Union3's `mb` values are already distance-modulus-like, so the primary cross-catalog comparison uses `MU`/`mb`-style values.

## Test logic

1. Fit a flexible empirical Pantheon reference curve using only Pantheon over the DES overlap range.
2. Check Union3 against that Pantheon curve.
3. Check DES against that Pantheon curve.
4. Compute STAM Model-A with locked `b = 354.95`.
5. Learn one global STAM magnitude offset from Pantheon and Union3 only.
6. Evaluate DES against locked STAM.
7. Compare:

```text
DES - Pantheon empirical curve
```

against:

```text
DES - STAM locked-b curve
```

If these are nearly the same, then the STAM-DES residual is largely the same as the independent DES-vs-Pantheon discrepancy.

## Outlier logic

The script reports DES outlier power in two ways:

1. absolute locked-b STAM residuals;
2. residuals after removing the DES median offset.

The second version isolates DES shape/outlier behavior after a global catalog offset is removed.
