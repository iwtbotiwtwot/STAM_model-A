# STAM Model-A black-hole decay / horizon recession

## Core rule

```text
A(r,t)=2GM(t)/(c²r)
r_h(t)=2GM(t)/c²
```

Mass balance:

```text
dM/dt = Ṁ_in - P_out/c²
```

So:

```text
dr_h/dt = (2G/c²)dM/dt
```

Decay condition:

```text
P_out > Ṁ_in c²
```

## Headline 10 solar-mass example

Initial horizon radius:

```text
29.533394 km
```

After 5% net mass loss:

```text
28.056724 km
```

Horizon change:

```text
-1.476670 km
```

At the old horizon radius after 5% mass loss:

```text
A = 0.950000
```

So the old horizon point becomes sub-threshold and the `A=1` surface has moved inward.

## Checks

```text
horizon radius linear with mass max error: 2.220e-16
A at old horizon equals mass fraction max error: 2.220e-16
decay scenarios have dr_h/dt < 0: True
growth scenarios have dr_h/dt > 0: True
```
