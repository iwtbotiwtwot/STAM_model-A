# `b` as catalog bridge, not core Model-A physics

## Working position

STAM Model-A should preserve the historical `b` work, but moving forward the cleaner physical theory should not treat `b` as fundamental.

The core no-`b` Model-A distance structure is:

```text
D_geo(z)    = L z (1 + 0.15 z)

D_adj,0(z)  = L z (1 + 0.5 z)

D_excess,0(z) = 0.35 L z²
```

The historical/catalog comparison form is:

```text
D_catalog(z) ≈ D_adj,0(z) + b z
```

In this interpretation:

```text
b z
```

is a bridge term between current catalog-inferred supernova distances and STAM's no-`b` accumulation-distance law. It is not assumed to be a universal physical term.

## Why keep `b` documented?

`b` is still important because it measures how current supernova catalogs relate to the no-`b` STAM curve.

Historically useful values:

```text
Pantheon/Union-style bridge: b ≈ 354.95
Original retained bridge:    b ≈ 461.3626922
DES-style bridge:            b ≈ 1335.412792
```

These values should remain in the repository as historical/diagnostic results, not as silently retuned physical constants.

## Why move forward without `b`?

The current working hypothesis is:

```text
spacetime accumulation and distance are separate variables
```

Current luminosity-distance catalogs may include:

```text
geometric distance
+
accumulation/traversal effects
+
catalog calibration/inference effects
```

Therefore a mismatch between current catalog distances and the no-`b` STAM curve is not automatically a model failure. It is the main distance-interpretation question STAM is investigating.

## What weakens Model-A?

Model-A is not weakened merely because it predicts a flatter distance curve than current catalog distances.

Model-A is weakened if every mismatch requires an arbitrary after-the-fact bridge with no stable pattern, no independent explanation, and no observable separation between geometric and accumulation layers.

## What strengthens Model-A?

Model-A is strengthened if:

```text
1. no-b distance curves remain structured and meaningful;
2. b-like discrepancies are stable within similar catalog families;
3. catalog-to-catalog anomalies, such as DES vs Pantheon/Union3, are independently visible;
4. independent observables separate into STAM-predicted distance layers.
```

## Practical testing rule

Forward tests should report at least:

```text
no-b Model-A:
    D_adj,0(z) = L z (1 + 0.5 z)

historical bridge diagnostic:
    D_adj,b(z) = L z (1 + 0.5 z) + b z
```

The no-b result should be treated as the primary physical Model-A result. The `b` result should be treated as a catalog-bridge diagnostic.
