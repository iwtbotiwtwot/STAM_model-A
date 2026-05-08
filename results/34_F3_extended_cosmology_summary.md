# F3-Extended Cosmology — Bold-STAM Effective Stress-Energy at Cosmic Scale

## Purpose

Compute the bold-STAM effective stress-energy (energy density AND pressure components) at cosmological scales. Determine whether the effective equation of state is dark-energy-like (w ≈ -1) and check whether this produces a derived value of `A_0 ≈ 0.0265` matching the historical bridge term b ≈ 354.95 Mly.

## Setup

- Cosmic A profile (Model-B): `A_cosmo(r) = (r/L)²` from critical-density Poisson.
- Bold-STAM metric: `g_tt = -(1-A)c²`, `g_rr = 1/[(1-A)(1-A²)²]`.
- Compute Einstein tensor components G^t_t, G^r_r, G^θ_θ for this metric.
- Interpret as effective stress-energy under Einstein gravity:
  ρ_eff ∝ -G^t_t, p_r_eff ∝ G^r_r, p_t_eff ∝ G^θ_θ.
- Volume-average over Hubble sphere to get cosmic-mean stress-energy.
- Effective equation of state: w_eff = <p_avg> / <ρc²>.

## Pointwise stress-energy at sample cosmic positions

|         u |          A |   G_t_t |   G_r_r |   G_theta_theta |   rho_eff |   p_avg_eff |     w_eff |
|----------:|-----------:|--------:|--------:|----------------:|----------:|------------:|----------:|
| 0.0999198 | 0.00998396 | 3.09844 | 3.01937 |        3.03854  |  -3.09844 |    3.03215  | -0.978606 |
| 0.299936  | 0.0899614  | 3.78048 | 3.13083 |        3.24435  |  -3.78048 |    3.20651  | -0.848176 |
| 0.499951  | 0.249951   | 4.52725 | 3.12114 |        3.10562  |  -4.52725 |    3.11079  | -0.687127 |
| 0.699967  | 0.489954   | 4.11412 | 2.59484 |        1.76247  |  -4.11412 |    2.03993  | -0.495836 |
| 0.899983  | 0.809969   | 1.86699 | 1.44346 |       -0.335924 |  -1.86699 |    0.257205 | -0.137765 |


## Cosmic volume averages

```text

<ρ_eff>      = -3.313973 (relative units)
<p_r_eff>    = +2.286891
<p_t_eff>    = +1.289983
<p_avg>      = +1.622286
w_cosmic_avg = -0.489529
u_max used   = 0.95 (avoiding singularity at u=1)
```

**Dark-energy-like reference**: w = -1 (cosmological constant).

**Verdict on equation of state**: NOT dark-energy-like (w_avg = -0.490, far from -1). Bold-STAM effective stress-energy at cosmic scales does not reduce to a simple dark-energy component.


## Bridge-term derivation status

Whether this cosmic stress-energy *exactly* reproduces `A_0 = 0.0265` requires solving modified Friedmann equations with the effective T^μ_ν as input, then extracting the distance-redshift relation, then identifying the constant-ambient-A equivalent. That last step is non-trivial because the modification doesn't always reduce to a simple constant-A picture.

What this script can say:
- The bold-STAM effective stress-energy at cosmic scale has a specific equation of state w (see above).
- The magnitude of the cosmic-average ρ_eff sets the SCALE of any dark-energy-like modification.
- Whether the resulting Friedmann modification produces A_0 = 0.0265 specifically is a separate, more elaborate calculation.

**Honest interpretation:**
- Cosmic-average ρ_eff magnitude: 3.313973 (relative units)
- Cosmic-average w: -0.489529
- abs_rho_cosmic_avg gives the relative magnitude of effective DE; exact A_0 derivation requires modified-Friedmann solution.


## Caveats

- Static-spherical-around-observer approximation. True FRW cosmology with bold-STAM would require deriving the Model-B field equations in expanding-universe context.
- Used A_cosmo(r) = (r/L)² as derived consequence of critical-density Poisson. If actual cosmological A profile differs, results differ.
- Volume averaging up to u_max < 1 to avoid singular endpoint. Results depend weakly on u_max choice (most volume is at moderate u).
- Pressure components computed from numerical second derivatives — some numerical noise possible at very small/large u.

## Generated plots

- `plots/34_einstein_stress_energy.png`
- `plots/34_cosmic_averages.png`
