# STAM Model-A gravitational-wave propagation test

## Purpose

This test checks the STAM statement:

```text
Gravitational waves propagate locally at c.
```

Apparent slowing near high accumulation is treated as accumulated traversal delay, not as a lower intrinsic wave speed.

## Propagation rule tested

```text
t_obs = ∫ds/c + k∫A(s)ds/c
```

Equal-coupling Model-A case:

```text
k_GW = k_EM
```

Then gravitational waves and electromagnetic waves receive the same accumulation traversal delay through the same `A` field.

## Horizon rule

Model-A also gives:

```text
A = (v_escape/c)^2
```

So:

```text
A < 1  ⇔ outward escape allowed in principle
A = 1  ⇔ horizon threshold
A > 1  ⇔ outward escape over-threshold
```

The wave still locally propagates at `c`; the problem in `A>1` is that outward escape would require more than `c`.

## Run

```powershell
python scripts\13_gw_propagation_model_a.py
```

Outputs:

```text
results\gw_propagation\summary.json
results\gw_propagation\gw_em_solar_accumulation_delay_cases.csv
results\gw_propagation\gw_em_near_horizon_apparent_delay_toy.csv
results\gw_propagation\gw_horizon_escape_threshold.csv
```

## Status

This is a first-pass consistency test. It is not a full wave-equation solution in dynamical curved spacetime.
