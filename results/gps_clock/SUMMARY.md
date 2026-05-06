# STAM Model-A GPS satellite clock adjustment test

## Result

Using the STAM accumulation field:

```text
A(r)=2GM/(c²r)
```

and the weak-field clock mapping:

```text
dτ/dt ≈ 1 - A/2
```

the GPS geoid-reference calculation gives:

```text
STAM gravitational gain: 45.787467 μs/day
orbital kinematic loss: -7.213600 μs/day
net satellite gain:     38.573867 μs/day
```

Factory offset needed:

```text
Δf/f = -4.464567973354e-10
```

Adjusted frequency:

```text
f = 10.2299999954327 MHz
```

The published GPS factory offset is:

```text
Δf/f = -4.4647e-10
f = 10.2299999954326 MHz
```

## Meaning

This is a strong Model-A clock-rate alignment test. It does not use `b`, supernova distances, or the cosmological distance adjustment.

The simple spherical Earth surface calculation gives the same scale but is slightly low because the official GPS clock offset uses the geoid/effective potential reference, not a perfectly spherical Earth surface.
