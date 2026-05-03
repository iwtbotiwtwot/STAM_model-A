# STAM Model-A GPS satellite clock test

## Purpose

This test checks whether Model-A's local accumulation field can account for the gravitational part of GPS satellite clock adjustment.

It does not use supernova distances, `b`, or cosmological distance adjustment.

## STAM ingredients

```text
A(r) = 2GM/(c²r)
```

Weak-field clock mapping tested:

```text
dτ/dt ≈ 1 - A/2
```

Then the GPS gravitational clock-rate shift is:

```text
Δf/f_grav = (A_geoid - A_orbit)/2
```

For circular orbit kinematics:

```text
Δf/f_kin = -v²/(2c²)
```

## Headline result

At nominal GPS semimajor axis:

```text
r ≈ 26,560 km
```

using a GPS/geoid effective potential reference:

```text
STAM gravitational gain ≈ 45.787467 μs/day
kinematic loss          ≈ -7.213600 μs/day
net satellite gain      ≈ 38.573867 μs/day
```

Required factory clock offset:

```text
Δf/f ≈ -4.464567973354e-10
```

Adjusted nominal 10.23 MHz source:

```text
10.2299999954327 MHz
```

## Interpretation

This is a strong weak-field clock-rate alignment test:

```text
STAM accumulation difference between the geoid and GPS orbit
reproduces the gravitational clock-rate part of GPS.
```

Adding the standard kinematic term reproduces the published GPS net factory clock offset.

## Limitations

This is not a full GPS operations model. It does not include eccentricity, Sagnac correction, Earth multipole propagation, ionosphere, troposphere, relativistic broadcast clock correction terms, or satellite-specific ephemeris processing.
