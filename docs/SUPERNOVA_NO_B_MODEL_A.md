# Supernova catalogs vs no-b Model-A distance prediction

## Purpose

This test runs the current supernova catalogs against the no-`b` Model-A distance law.

Primary Model-A prediction:

```text
D_adj,0(z) = L z (1 + 0.5 z)
```

Geometric spine comparison:

```text
D_geo(z) = L z (1 + 0.15 z)
```

Historical diagnostic only:

```text
D_adj,b(z) = L z (1 + 0.5 z) + 354.95 z
```

## Interpretation

This test does not treat `b` as core physics. It asks how current catalog-inferred distances sit relative to the no-`b` Model-A distance prediction.

A positive residual means:

```text
catalog-inferred distance is larger than the no-b Model-A prediction
```

This is not automatically a Model-A failure. The working STAM question is whether current supernova distance inference includes accumulation/catalog effects that should be separated from geometric distance.

## Run

```powershell
python scripts\11_supernova_no_b_model_a.py `
  --union3 data\union3_bins.csv `
  --pantheon data\pantheon.csv `
  --des data\des.csv
```

Outputs:

```text
results\no_b_supernova\
```
