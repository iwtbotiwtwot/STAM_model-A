# STAM Model-A CMB angular-scale toy test

## Purpose

This test asks whether current no-b Model-A distance layers can reproduce the observed CMB acoustic angular scale under simple mappings.

It does not claim to model:

```text
Big Bang origin
recombination
CMB microphysics
CMB power spectrum
```

## Observational reference

Planck reports:

```text
100θ* ≈ 1.0411
θ* ≈ 0.010411 rad
```

## Tested idea

For a sound horizon scale:

```text
r_s ≈ 144.4 Mpc
```

a simple angular relation would be:

```text
θ ≈ r_s / D
```

The script tests several simple candidates for `D` using STAM layers at `z≈1090`.

## Current result

No simple tested STAM layer reproduces the observed CMB acoustic angle at `z≈1090`.

This means Model-A needs a CMB-facing distance/visibility rule before CMB angular scales can be considered addressed.

## Run

```powershell
python scripts\18_cmb_angular_scale_toy.py
```

Outputs:

```text
results\cmb_angular_toy\summary.json
results\cmb_angular_toy\cmb_angular_mapping_results.csv
results\cmb_angular_toy\cmb_angular_best_mappings.csv
```
