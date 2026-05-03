# STAM_model-A Theory Notes

## Working status

STAM_model-A is a proposed accumulation-field framework. The current repository is a research scaffold designed to separate:

```text
algebraic identities
phenomenological fits
physical derivations
independent predictions
falsification tests
```

Only the first category is implemented as passing smoke tests at initialization.

## Postulates

### 1. Accumulation field

Spacetime is assigned a dimensionless scalar accumulation field:

```text
A = A(x)
```

### 2. Source principle

Matter-energy sources accumulation. In the local, static, spherically symmetric weak-field limit:

```text
A(r) = Rs/r = 2GM/(c^2 r)
```

### 3. Local dynamics from gradients

The local gravitational acceleration field is identified operationally as:

```text
g_vec = (c^2/2) grad(A)
```

For `A(r)=Rs/r`:

```text
g_vec = -GM/r^2 r_hat
```

The minus sign means the acceleration points inward.

### 4. Propagation from path integrals

Accumulation contributes traversal delay through:

```text
Delta t = (1/c) int A ds
```

For a point-source weak-field path:

```text
Delta t = (2GM/c^3) int ds/r
```

This reproduces the logarithmic structure of weak-field Shapiro delay. The next requirement is coefficient-level comparison against standard Solar System tests.

### 5. Critical threshold

The model assigns a critical propagation threshold at:

```text
A = 1
```

Using `A(r)=Rs/r`, this occurs at:

```text
r = Rs
```

This is a horizon-like criterion. The current model still needs a complete invariant dynamical account of the `A >= 1` regime.

### 6. Cosmological path accumulation

STAM_model-A decomposes inferred distance into a geometric proxy and path-integrated excess:

```text
D_adj(z) = D_geo(z) + D_excess(z)
```

with:

```text
D_geo(z)    = L z (1 + 0.15z)
D_adj(z)    = L z (1 + 0.5z) + b z
D_excess(z) = b z + 0.35 L z^2
```

The implied average and local path accumulations are:

```text
<A_path>(z)      = D_excess / D_geo
A_path,local(z) = d(D_excess)/d(D_geo)
```

Explicitly:

```text
<A_path>(z)      = (b/L + 0.35z)/(1 + 0.15z)
A_path,local(z) = (b/L + 0.70z)/(1 + 0.30z)
```

The key open problem is deriving this redshift dependence from independent physical structure instead of defining it from the distance law.

## Major unresolved theoretical tasks

1. Derive `A_path,local(z)` independently.
2. Establish whether the local `A(r)` and cosmological `A_path(z)` are limits of one field equation.
3. Specify the `A >= 1` regime in coordinate-invariant language.
4. Close the normalization across acceleration, Shapiro delay, gravitational redshift, and lensing.
5. Test locked-parameter predictions against independent catalogs.
