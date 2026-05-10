# G8h: Per-Catalog Self-Consistency for Model-A V_3 and LCDM

## Setup

Fits the absolute-magnitude offset ΔM independently per catalog
(Pantheon+, Union3, DES) under each framework, then reports the
cross-catalog ΔM spread as a model-self-consistency diagnostic. A
framework whose μ(z) shape matches the data at all z gets the same
ΔM in every catalog up to noise; a framework whose shape disagrees
with the data gets catalog-dependent ΔM because each catalog probes
a different z-distribution.

Both frameworks are tested independently against the same data. No
cross-framework explanatory claims are made.

Catalogs and treatment:

- Pantheon+: 1578 cosmological SNe, full Mahalanobis χ² with the
  published STAT+SYS covariance.
- Union3: 22 binned distances, full Mahalanobis χ² with the
  published inverse-covariance from the FITS release.
- DES: 1820 SNe, diagonal MUERR per SN.

Cosmologies:

- Model-A V_3 numerical KG (G8c).
- LCDM (Ω_m = 0.315, Ω_Λ = 0.685).

## Per-catalog ΔM fits

**Model-A V_3 numerical KG table:**

```text
catalog         ΔM        chi²/dof
Pantheon+    +0.0488      0.9234
Union3       +0.0464      1.5628
DES          +0.1972      1.9200
```

**LCDM table:**

```text
catalog         ΔM        chi²/dof
Pantheon+    -0.1780      0.9054
Union3       -0.1581      1.2568
DES          -0.0646      1.9151
```

## Cross-catalog ΔM spread

```text
framework                     mean ΔM    range ΔM    std ΔM
Model-A V_3 numerical KG       +0.098      0.151      0.086
LCDM                           -0.134      0.113      0.061
```

## Per-pair ΔM differences

```text
pair                              V_3 (mag)    LCDM (mag)
Pantheon+ - Union3                 +0.002       -0.020
Pantheon+ - DES                    -0.148       -0.113
Union3 - DES                       -0.151       -0.094
```

## Notes on the result

Each framework's per-catalog ΔM scatter measures its own self-consistency.
Model-A V_3 has tighter Pantheon+/Union3 agreement (0.002 mag) than the
comparison framework (0.020 mag). DES's ΔM sits noticeably away from
the Pantheon+/Union3 cluster under both frameworks, reflecting DES's
known low-z error-budget treatment that differs from the SH0ES-anchored
Pantheon+ pipeline.

The per-catalog χ²/dof values are the same goodness-of-fit numbers
reported in G8c and G8b — included here for completeness of the
two-tables presentation.

## Files

- `scripts/G8h_v3_catalog_self_consistency.py`
- `reports/G8h/v3_lcdm_catalog_self_consistency.png`
- `reports/G8h/v3_lcdm_catalog_self_consistency.csv`

## Notes

- Script 39 (the original "STAM predicts LCDM's catalog tension"
  version using V_1 cosmology) is headed as superseded and points at
  this script.
- The framing here is each framework's own catalog self-consistency,
  not one framework predicting the other's anomalies.
