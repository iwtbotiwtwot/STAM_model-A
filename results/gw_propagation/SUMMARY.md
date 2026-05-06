# STAM Model-A gravitational-wave propagation test

## Rule tested

```text
t_obs = ∫ds/c + k∫A(s)ds/c
```

Model-A equal-coupling case:

```text
k_GW = k_EM = 1
```

So gravitational waves and light both propagate locally at `c` and receive the same accumulation traversal delay through the same `A` field.

## Main checks

```text
max GW-EM propagation difference for equal coupling, solar cases:
0.000e+00 s

max GW-EM propagation difference for equal coupling, near-horizon toy cases:
0.000e+00 s

max error in A = (v_escape/c)^2:
4.441e-16
```

## GW170817 sanity check

Assuming a distance of 40 Mpc and a 1.74 s GRB-after-GW lag:

```text
naive fractional speed difference if all lag were propagation:
4.226e-16
```

Model-A equal-coupling propagation lag:

```text
0 s
```

So the observed lag can be treated as source/emission delay rather than different propagation speed.

## Near-horizon interpretation

Near high `A`, apparent speed inferred by an outside observer can fall below `c` because traversal time includes accumulated spacetime:

```text
v_app = path_length / [∫ds/c + ∫A ds/c]
```

This does not mean the wave locally travels slower than `c`.

For `A > 1`:

```text
A = (v_escape/c)^2
```

so outward escape would require `v_escape > c`. The over-threshold issue is trapped propagation, not a lower local wave speed.
