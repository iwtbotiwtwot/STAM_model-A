# STAM Model-A Strong-Field SU Horizon Integral Test

## Purpose

This test evaluates the proposed strong-field traversal layer `S_h(A)=1/(1-A)` with `A(r)=Rs/r`. The goal is to compare ordinary radial distance against the strong-field SU measure as the inner radius approaches the horizon from outside.

## Core equations

```text
A(r) = Rs/r
S_h(A) = 1/(1-A)
dSU_A = dr/(1-Rs/r)
SU_A = ∫ r/(r-Rs) dr
SU_A = r + Rs ln|r-Rs| + constant
```

For an outward interval outside the horizon, using dimensionless units `x=r/Rs`:

```text
SU_A(x_inner → x_outer) = (x_outer - x_inner) + ln[(x_outer - 1)/(x_inner - 1)]
```

## Scale factor table

|           A |    S_h(A)=1/(1-A) |
|------------:|------------------:|
| 0           |          1        |
| 0.5         |          2        |
| 0.9         |         10        |
| 0.99        |        100        |
| 0.999       |       1000        |
| 0.999999    |     999999.999971 |
| 0.999999999 | 1000000028.28     |


## Interval comparison table

|   outer_r_over_Rs |   inner_r_over_Rs |        A_outer |        A_inner |   ordinary_delta_r_Rs |        SU_A_Rs |   SU_A_over_delta_r |
|------------------:|------------------:|---------------:|---------------:|----------------------:|---------------:|--------------------:|
|         10        |       2           | 0.1            | 0.5            |     8                 | 10.1972245773  |       1.27465307217 |
|          2        |       1.1         | 0.5            | 0.909090909091 |     0.9               |  3.20258509299 |       3.5584278811  |
|          1.1      |       1.01        | 0.909090909091 | 0.990099009901 |     0.09              |  2.39258509299 |      26.584278811   |
|          1.01     |       1.001       | 0.990099009901 | 0.999000999001 |     0.009             |  2.31158509299 |     256.84278811    |
|          1.001    |       1.000001    | 0.999000999001 | 0.999999000001 |     0.000999          |  6.90875427906 |    6915.66994901    |
|          1.000001 |       1.000000001 | 0.999999000001 | 0.999999999    |     9.98999999835e-07 |  6.90775619516 | 6914670.86717       |


## Interpretation

The ordinary radial interval can become very small near the horizon, while `SU_A` remains large and increases without bound as `r_inner/Rs → 1+`. This matches the intended strong-field behavior: the A=1 horizon is already derived from `A=(v_escape/c)^2`, and the candidate `S_h(A)` layer turns that horizon threshold into a divergent outside-comparison/traversal measure.

This does not by itself complete the strong-field STAM claim. It does establish that the proposed local SU layer produces the expected logarithmic horizon divergence while preserving finite ordinary radial distance outside the horizon.

## Generated plots

- `plots\30_sh_of_A_near_horizon.png`
- `plots\30_su_over_delta_r_vs_inner_radius.png`
- `plots\30_ordinary_distance_vs_strong_field_su.png`
- `plots\30_log_divergence_su_to_2Rs.png`
