# STAM Model-A Shapiro delay test

## Purpose

This test moves away from cosmological distance and checks the local propagation side of Model-A.

Model-A uses:

```text
A(r) = Rs / r = 2GM / (c² r)
```

and:

```text
Delta t = (1/c) ∫ A(r) ds
```

For a straight weak-field path with impact parameter `b_imp`:

```text
Delta t = (2GM/c³) [asinh(x2/b_imp) - asinh(x1/b_imp)]
```

This is equivalent to the first-order logarithmic Shapiro structure:

```text
Delta t = (2GM/c³) ln[(r1+r2+R)/(r1+r2-R)]
```

## Interpretation

Passing this test means:

```text
STAM Model-A reproduces the first-order weak-field Shapiro delay structure.
```

It does not yet mean:

```text
full PPN equivalence
strong-field validation
complete solar-system precision validation
```

Those require separate tests.

## Run

```powershell
python scripts\09_shapiro_delay_model_a.py
```

Outputs:

```text
results\shapiro_delay\summary.json
results\shapiro_delay\shapiro_delay_cases.csv
```

## Headline solar-grazing Earth-Mars result

The script predicts a one-way delay of roughly:

```text
123.5 microseconds
```

and a two-way delay of roughly:

```text
247 microseconds
```

for a solar-grazing Earth-Mars-like path.
