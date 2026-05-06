# STAM Model-A black-hole decay and horizon recession

## Purpose

This test formalizes the STAM idea:

```text
A black hole loses mass if energy loss exceeds mass-energy accumulation.
As M decreases, the A=1 horizon recedes.
```

## Core formulas

```text
A(r,t) = 2GM(t)/(c²r)
```

```text
r_h(t) = 2GM(t)/c²
```

Mass balance:

```text
dM/dt = Ṁ_in - P_out/c²
```

Horizon recession/expansion:

```text
dr_h/dt = (2G/c²)dM/dt
```

Decay condition:

```text
P_out > Ṁ_in c²
```

## Interpretation

If net mass-energy decreases:

```text
dM/dt < 0
```

then:

```text
dr_h/dt < 0
```

and the horizon moves inward.

At a fixed old horizon radius, `A` drops below 1 after sufficient mass loss.

## Run

```powershell
python scripts\14_black_hole_decay_model_a.py
```

Outputs:

```text
results\black_hole_decay\summary.json
results\black_hole_decay\mass_fraction_horizon_recession.csv
results\black_hole_decay\mass_balance_rate_scenarios.csv
results\black_hole_decay\optional_hawking_reference_examples.csv
```

## Limitations

This is a spherical, non-rotating mass-balance threshold model. It does not derive a microscopic decay process, model Kerr rotation, charge, accretion disks, or black-hole thermodynamics from first principles.
