# G22: Proper-time of radial infall in Model-A vs Schwarzschild

**Date:** 2026-05-11

**Setup.** Test particle released from rest at A = A_0 = 1/(12π) ≈ 0.026526. Free-fall radial infall. Compute the proper time τ the particle's own clock reads as it reaches various A values.

**Metric inputs (natural units c=1, Rs=1):**

- GR Schwarzschild: g_tt = -(1-A), g_rr = 1/(1-A); k/h = 1
- Model-A: g_tt = -(1-A), g_rr = 1/[(1-A)(1-A²)²]; k/h = (1-A²)²

**Geodesic (radial timelike free-fall):**

```
(dr/dτ)² = (k/h)·(A − A_release)
τ(A_end) = ∫_{A_release}^{A_end} dA / [A²·√((k/h)·(A − A_release))]
```

## Numerical results (proper time in units of Rs/c)

| landmark             |        A |      r/Rs |   tau_GR (Rs/c) |   tau_MA (Rs/c) |   ratio MA/GR |   diff MA-GR (Rs/c) |
|:---------------------|---------:|----------:|----------------:|----------------:|--------------:|--------------------:|
| release (A_0)        | 0.026526 | 37.699112 |        0.000000 |        0.000000 |    nan        |            0.000000 |
| A = 0.05             | 0.050000 | 20.000000 |      290.248627 |      290.555429 |      1.001057 |            0.306803 |
| A = 0.10             | 0.100000 | 10.000000 |      340.545217 |      341.089019 |      1.001597 |            0.543802 |
| A = 0.20             | 0.200000 |  5.000000 |      355.820670 |      356.662010 |      1.002365 |            0.841340 |
| A = 1/3 (ISCO)       | 0.333333 |  3.000000 |      360.043353 |      361.180615 |      1.003159 |            1.137262 |
| A = 0.50             | 0.500000 |  2.000000 |      361.677360 |      363.140095 |      1.004044 |            1.462735 |
| A = 2/3 (photon sph) | 0.666667 |  1.500000 |      362.354188 |      364.158898 |      1.004981 |            1.804710 |
| A = 0.80             | 0.800000 |  1.250000 |      362.652741 |      364.806231 |      1.005938 |            2.153490 |
| A = 0.90             | 0.900000 |  1.111111 |      362.806048 |      365.370119 |      1.007067 |            2.564072 |
| A = 0.99             | 0.990000 |  1.010101 |      362.911559 |      366.691141 |      1.010415 |            3.779581 |
| A = 0.999            | 0.999000 |  1.001001 |      362.920809 |      367.871909 |      1.013642 |            4.951100 |
| A = 0.9999           | 0.999900 |  1.000100 |      362.921722 |      369.040157 |      1.016859 |            6.118435 |
| A = 0.99999          | 0.999990 |  1.000010 |      362.921813 |      370.207167 |      1.020074 |            7.285354 |

## Divergence coefficient (Model-A near A=1)

Fit form: τ_MA(A) ≈ −slope · ln(1−A) + intercept (near A=1)

- Fitted slope:                       0.508834
- Analytical 1/(2·√(1−A_0)):          0.506766
- Ratio (fit/analytical):              1.004079
- Fitted intercept:                    364.351878

## Files

- `plots/G22_proper_time_infall.png`
- `plots/G22_ratio_and_diff.png`
- `results/G22_proper_time_table.csv`
