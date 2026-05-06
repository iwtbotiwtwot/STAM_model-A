# STAM Model-A horizon threshold test

## Core result

Model-A defines:

```text
A(r) = Rs/r = 2GM/(c²r)
```

Escape speed is:

```text
v_escape² = 2GM/r
```

Therefore:

```text
A = v_escape²/c²
```

So:

```text
A = 1 ⇔ v_escape = c
A > 1 ⇔ v_escape > c
```

## Checks

```text
max |A(r_h)-1| = 0.000e+00
max |v_escape(r_h)/c - 1| = 0.000e+00
max |M_recovered/M - 1| = 0.000e+00
```

## Useful values

```text
1 solar mass horizon radius  ≈ 2.953339 km
10 solar mass horizon radius ≈ 29.533394 km
```

## Interpretation

This is a clean non-distance Model-A identity test.

The STAM threshold `A=1` is not arbitrary in the spherical case. It is exactly the point where the escape-speed requirement reaches the speed of light.
