# STAM Model-A high-z stress test

## Purpose

This test asks whether the current no-b Model-A formulas remain mathematically coherent when extrapolated toward high redshift.

It does not claim STAM explains:

```text
the origin of the Big Bang
recombination
CMB microphysics
the origin of redshift
```

## Formulas tested

```text
D_geo(z)      = Lz(1 + 0.15z)
D_adj,0(z)    = Lz(1 + 0.5z)
D_excess,0(z) = 0.35Lz²
```

Path accumulation:

```text
<A_path>(z)      = 0.35z/(1 + 0.15z)
A_path,local(z)  = 0.70z/(1 + 0.30z)
```

## Main high-z behavior

As:

```text
z → ∞
```

the path accumulation functions approach:

```text
7/3 ≈ 2.333333
```

rather than diverging.

## Run

```powershell
python scripts\17_high_z_stam_stress_test.py
```

Outputs:

```text
results\high_z_stress\summary.json
results\high_z_stress\high_z_model_a_grid.csv
results\high_z_stress\key_redshift_values.csv
results\high_z_stress\A_path_saturation_redshifts.csv
```
