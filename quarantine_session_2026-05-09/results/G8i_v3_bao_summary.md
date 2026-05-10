# G8i: Model-A V_3 numerical KG vs DESI DR1 BAO

## Setup

Tests the V_3 numerical KG cosmology (G8c) against DESI DR1 BAO
distance measurements. The cosmology is fully derived from structural
commitments with no fitted cosmological parameters:

- A_0 = 1/(12π) (G7 commitment)
- α/β = [A_0/(1−A_0)]² (V_3 minimum at A_0)
- β_tilde = 0.4265 (Model-A internal closure: numerical-KG h²(today) = 1)
- Ω_m = 0.315, Ω_r = 9.2e-5

The framework predicts comoving distance D_M(z) = (c/H_0) ∫₀^z dz'/h(z').
A single nuisance r_d (sound horizon at the drag epoch) is fit such that
the predicted ratio D_M(z)/r_d matches the DESI observable.

## Result

```text
Best-fit r_d:                  128.40 Mpc
chi^2 (5 points)               = 8.03
chi^2 / dof (4 dof)            = 2.01
p-value                        ≈ 0.09
```

Per-point residuals:

```text
z         D_M/r_d obs    D_M/r_d pred    residual (sigma)
0.510         13.62          13.54          -0.32
0.706         16.85          17.58          +2.28
0.930         21.71          21.62          -0.33
1.317         27.79          27.45          -0.49
2.330         39.71          38.26          -1.54
```

## Reading the result

Four points at z = 0.510, 0.930, 1.317, and 2.330 sit between −0.32σ
and −1.54σ of Model-A's predicted curve. The point at z = 0.706 (LRG2)
deviates by +2.28σ. The combined chi² is 8.03; removing the LRG2 point
gives chi² = 2.82 across 3 dof (chi²/dof = 0.94), which quantifies
that the LRG2 residual contributes the majority of the total chi².

The +2σ-class deviation at LRG2 has been discussed in the DESI DR1
literature for multiple cosmological frameworks. Whether it reflects
a data-side feature or a real model tension remains open.

## Sound horizon prediction (open follow-up)

The BAO-aligned r_d = 128.40 Mpc is the sound horizon at the drag
epoch that the DESI data favors under V_3 cosmology. This is a specific
internal Model-A prediction that should be cross-checked against a
first-principles calculation: integrating c_s(z)/H(z) under Model-A's
early-universe content from pre-recombination through the drag epoch.

If the first-principles calculation produces ≈128 Mpc, the framework
is internally consistent across BAO and early-universe physics. If it
produces a different value, the gap identifies where Model-A's
pre-recombination physics needs further work.

This calculation has not been done.

## Multi-regime test record for Model-A V_3

| Regime | Test | χ²/dof |
|---|---|---:|
| Low-mid-z SN distance | G8c (Pantheon+ + Union3, full Mahalanobis) | 0.93 |
| Mid-z BAO standard ruler | G8i (DESI DR1 D_M/r_d, 5 points) | 2.01 |

## Files

- `scripts/G8i_v3_bao_test.py`
- `reports/G8i/v3_vs_lcdm_bao.png`
- `reports/G8i/v3_vs_lcdm_bao_summary.csv`
- `reports/G8i/v3_vs_lcdm_bao_per_point.csv`

## Open items

1. First-principles r_d calculation under V_3 cosmology, compared to
   the BAO-aligned value of 128.40 Mpc.
2. Re-run against DESI DR2 BAO when released.
