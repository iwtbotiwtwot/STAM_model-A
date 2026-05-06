# STAM Model-A horizon threshold test

## Purpose

This test checks the STAM horizon claim:

```text
A = 1
```

as both:

```text
r = Rs = 2GM/c²
```

and:

```text
v_escape = c
```

## Core identity

Model-A defines:

```text
A(r) = Rs/r = 2GM/(c²r)
```

Escape velocity is:

```text
v_escape² = 2GM/r
```

Therefore:

```text
A = v_escape²/c²
```

So:

```text
A < 1  ⇔  v_escape < c
A = 1  ⇔  v_escape = c
A > 1  ⇔  v_escape > c
```

## Interpretation

This is a clean non-distance Model-A result. It uses no `b`, no supernova catalogs, and no cosmological distance adjustment.

## Run

```powershell
python scripts\12_horizon_threshold_model_a.py
```

Outputs:

```text
results\horizon_threshold\summary.json
results\horizon_threshold\horizon_radius_from_mass_examples.csv
results\horizon_threshold\A_escape_threshold_scan_10Msun.csv
```

## Limitations

This is a spherical, non-rotating threshold test. It does not yet model rotating black holes, charged black holes, photon spheres, observed shadow size, accretion physics, or interior dynamics.
