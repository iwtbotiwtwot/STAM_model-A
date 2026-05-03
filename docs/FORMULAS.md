# Formula Reference

## Constants

```text
G = 6.67430e-11 m^3 kg^-1 s^-2
c = 299792458 m s^-1
C = 3.261563776
H = 0.000243635
L = C/H = 13387.090426252385
```

`C`, `H`, and `L` are preserved in the STAM distance scale used in the author's current notes. `L` is used as the scale in the cosmological formulas.

## Local accumulation

```text
Rs = 2GM/c^2
A(r) = Rs/r = 2GM/(c^2 r)
```

Radial derivative:

```text
dA/dr = -Rs/r^2 = -2GM/(c^2 r^2)
```

Gravity bridge:

```text
g_vec = (c^2/2) grad(A)
```

Spherical weak-field result:

```text
g_radial = (c^2/2) dA/dr
         = (c^2/2)(-2GM/(c^2 r^2))
         = -GM/r^2
```

Magnitude:

```text
|g| = GM/r^2
```

## Propagation delay

```text
Delta t = (1/c) int A(r) ds
```

Substitute the spherical weak-field form:

```text
Delta t = (2GM/c^3) int ds/r
```

For a straight path with impact parameter `impact_b` and coordinate `x`:

```text
r(x) = sqrt(x^2 + impact_b^2)
```

so:

```text
int dx/r = asinh(x/impact_b)
```

and:

```text
Delta t = (2GM/c^3) [asinh(x2/impact_b) - asinh(x1/impact_b)]
```

## Horizon threshold

```text
A = 1
Rs/r = 1
r = Rs
```

Black-hole mass from threshold radius:

```text
M = c^2 r_h / (2G)
```

## Mass inference formulas

```text
Horizon:              M = c^2 r_h / (2G)
Acceleration:         M = g r^2 / G
Orbital velocity:     M = v^2 r / G
Shapiro coefficient:  M = K c^3 / (2G)
Gravitational shift:  M ~= z_grav c^2 r / G
Lensing deflection:   M ~= alpha c^2 b / (4G)
```

These formulas are operational estimators under their stated approximations. The repo tests whether synthetic observables generated from one mass recover the same mass.

## Cosmology

```text
D_geo(z)    = L z (1 + 0.15z)
D_adj(z)    = L z (1 + 0.5z) + b z
D_excess(z) = D_adj - D_geo
             = b z + 0.35 L z^2
```

Average path accumulation:

```text
<A_path>(z) = D_excess/D_geo
            = (b/L + 0.35z)/(1 + 0.15z)
```

Local path accumulation:

```text
A_path,local(z) = d(D_excess)/d(D_geo)
                = (b + 0.70Lz)/(L(1 + 0.30z))
                = (b/L + 0.70z)/(1 + 0.30z)
```

Path-integral identity:

```text
D_adj(z) = int_0^z [1 + A_path,local(u)] dD_geo/du du
```

where:

```text
dD_geo/dz = L(1 + 0.30z)
```

and:

```text
[1 + A_path,local(z)] dD_geo/dz
= L(1 + 0.30z) + b + 0.70Lz
= L + b + Lz
```

Integrating from `0` to `z`:

```text
D_adj(z) = (L+b)z + 0.5Lz^2
         = Lz(1 + 0.5z) + bz
```
