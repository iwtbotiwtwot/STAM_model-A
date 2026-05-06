# BAO mapping flip and the role of `b`

## Question

The first BAO geometric split test showed a flip:

```text
Full all-z BAO:
    D_geo / (1 + z) performed better.

Low/mid-z BAO, z <= 1.6:
    D_adj / (1 + z) performed better.
```

Does `b = 354.95` explain this flip?

## Short answer

No, not in the relevant physical sense.

The flip persists when:

```text
b = 0
```

That means the flip is mainly caused by the difference between the STAM geometric layer and the no-b accumulation-adjusted layer, plus the influence of the high-z BAO points. It is not caused by the Pantheon/Union bridge value `b = 354.95`.

## Quick numbers from the first compressed BAO check

Combined all-z BAO:

```text
D_geo/(1+z):
    chi2/dof ≈ 9.235

D_adj/(1+z), b = 0:
    chi2/dof ≈ 13.332

D_adj/(1+z), b = 354.95:
    chi2/dof ≈ 12.250
```

Combined z <= 1.6 BAO:

```text
D_geo/(1+z):
    chi2/dof ≈ 9.382

D_adj/(1+z), b = 0:
    chi2/dof ≈ 2.759

D_adj/(1+z), b = 354.95:
    chi2/dof ≈ 2.481
```

So the low-z preference for the adjusted layer is already present without `b`.

## Interpretation

`b` modestly changes the adjusted-layer curve. It does not produce the qualitative flip.

The flip is better described as:

```text
At low/mid redshift, BAO compressed D_M/r_d points prefer the accumulation-adjusted transverse mapping.

When the high-z BAO points are included, the geometric transverse mapping becomes better.
```

## Important caution

This was a first-pass compressed BAO test.

It used:

```text
D_M / r_d
```

with one nuisance scale:

```text
r_d_eff
```

It did not use full covariance matrices and did not test `D_H` or `D_V`.

So the correct status is:

```text
BAO result: meaningful but incomplete.
b explanation of flip: not supported.
```
