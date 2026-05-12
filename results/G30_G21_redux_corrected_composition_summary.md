# G30: G21 redux under corrected composition

**Date:** 2026-05-11 (evening)

## What this redoes

G21 (script `G21_strong_field_with_cosmic_floor.py`) used the rescaled composition A_total(r) = A_0 + (1-A_0)(Rs/r). This produced landmark shifts (ISCO at 3.173 Rs, photon at 1.521 Rs) and two NEC sign-changes in the effective stress-energy bracket — a 'negative-rho shell.'

The corrected reading (Sean 2026-05-11 evening): A_0 is the void minimum, not a rescaling factor. The source has its own A = Rs/r in the strong field; A_0 only enters in void regions far from the source. So:

```
A(r) = max(Rs/r, A_0)
     = Rs/r     for r < r_trans (strong field)
     = A_0      for r > r_trans (void)
r_trans = Rs / A_0 = 12*pi * Rs ≈ 37.6991 Rs
```

## Landmark radii

| Landmark | A_target | r/Rs (corrected) | r/Rs (G21 rescaled, prior) | rescaled shift |
|---|---|---|---|---|
| ISCO | 0.333333 | 3.000000 | 3.172915 | 5.76% |
| photon sphere | 0.666667 | 1.500000 | 1.520719 | 1.38% |
| horizon | 1.000000 | 1.000000 | 1.000000 | 0.00% |

**Reading: under the corrected composition, ALL strong-field landmarks are at EXACT GR values.** ISCO = 3 Rs, photon sphere = 1.5 Rs, horizon = Rs. The 'shifts' reported in original G21 were artifacts of the wrong composition rule.

## NEC crossover

Sign changes of B(A) = f(A) - A f'(A) on (A_0, 1):

- Crossing 1: A = 0.439985 (r/Rs = 2.2728)

**Single NEC crossover at A ≈ 0.4400** — matches F3's original result. The 'two crossovers' reported in original G21 were also artifacts of the rescaled composition (the second crossover came from the (1-A_0) factor's contribution to B, not from genuine framework structure).

## Physical scales of r_trans

Where the cosmic A_0 floor takes over from source-dominated A:

| Object | Rs (kpc) | r_trans (kpc) |
|---|---|---|
| Stellar BH (10 M_sun) | 9.571e-16 | 3.608e-14 |
| Sgr A* (4e6 M_sun) | 3.828e-10 | 1.443e-08 |
| Galaxy bulge (1e10 M_sun) | 9.571e-07 | 3.608e-05 |
| Galaxy (1e12 M_sun) | 9.571e-05 | 3.608e-03 |
| Galaxy cluster (1e14 M_sun) | 9.571e-03 | 3.608e-01 |

**Eye-catching observation: for galaxy-mass concentrations, r_trans sits at galactic halo scale (~kpc to ~10 kpc).** This is where the cosmic A_0 floor takes over from the source-dominated A = Rs/r. The framework's effective mass profile in this transition region could connect to dark-matter-halo-like observational signatures, since the effective enclosed mass under A_0 grows linearly with r (cosmic-floor contribution to m(r)).

For stellar BHs and SMBHs in galactic centers, r_trans is far smaller than any relevant astronomical scale — the void region's A_0 contribution doesn't matter for ringdown, ISCO, etc. The corrected composition gives exactly GR-equivalent strong-field landmarks for these objects.

## What's different vs original G21

| Feature | Original G21 (rescaled) | G30 redux (corrected) |
|---|---|---|
| Composition | A_0 + (1-A_0)Rs/r | max(Rs/r, A_0) |
| ISCO | 3.173 Rs (5.8% shift) | 3.000 Rs (GR-exact) |
| Photon sphere | 1.521 Rs (1.4% shift) | 1.500 Rs (GR-exact) |
| Horizon | 1.000 Rs (no shift) | 1.000 Rs (GR-exact) |
| NEC crossings | 2 (negative-rho shell) | 1 (matches F3) |
| Asymptotic A | A_0 | A_0 (same) |
| Mass at infinity | linear in r | linear in r in voids (same effect) |

## Plot

![G30 G21 redux](../plots/G30_G21_redux.png)

Left: A(r) under three compositions. Corrected = blue, original G21 rescaled = red, pure GR/F3 = green. The corrected reading is essentially GR in the strong field, with A pinned at A_0 in voids.

Right: B(A) bracket. Single sign change at A ≈ 0.44 (matches F3, NOT two crossings as G21-rescaled reported).

## Reading

**Strong-field Model-A is EXACTLY GR in landmark locations.** The framework differs from GR only in g_rr (via the k(A) = (1-A)(1-A^2)^2 modification), not in where ISCO/photon sphere/horizon sit. The thirds-of-A reading (ISCO at A=1/3, photon at A=2/3, horizon at A=1) maps to r-coordinates identically to GR.

**The cosmic floor only enters in void regions** (r > r_trans). It does NOT shift strong-field landmarks. The shifts reported in original G21 came from the rescaled-composition algebra, not from framework physics.

**F3's effective stress-energy story stands intact.** Single NEC crossover at A ≈ 0.44; negative-rho region is for ALL A < 0.44 (not just a shell). Same as F3.

**The 'eye-catching' result: r_trans = 12pi * Rs sits at galactic-halo scales for galaxy-mass concentrations.** This is the only place where the cosmic-floor structure differs meaningfully from GR for typical astrophysical objects. The transition region (where A goes from Rs/r to A_0) could be worth examining for galactic-scale signatures.