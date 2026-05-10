# G8g: Model-A V_3 vs LCDM — Distance Prediction Comparison

## Setup

Computes the distance modulus μ(z) under two independent frameworks —
Model-A V_3 numerical KG (G8c) and LCDM — and characterizes the shape
of the difference Δμ(z) = μ_V3(z) − μ_LCDM(z) as a function of z.

Cosmology parameters used:

- Model-A V_3 numerical KG (H_0 = 73.04 km/s/Mpc): A_0 = 1/(12π), α/β
  = [A_0/(1−A_0)]², β_tilde = 0.4265 (Model-A internal closure).
- LCDM (H_0 = 67.4 km/s/Mpc): Ω_m = 0.315, Ω_Λ = 0.685, Ω_r = 9.2e−5.

This is the V_3 numerical-KG analog of the original script 38, recast
under the two-tables framing: both frameworks are independent
predictions, neither is the reference for the other.

## Prediction difference

```text
z            Δμ = μ_V3 − μ_LCDM (mag)
0.05               -0.191
0.10               -0.206
0.30               -0.253
0.50               -0.286
1.00               -0.329
2.00               -0.346
```

## Functional form fits to Δμ(z)

```text
form                                              RMS residual    max |residual|
constant       Δμ = -0.3154                       0.0433          0.1406
linear in z    Δμ = -0.2525 - 0.0503 z            0.0235          0.0777
quadratic      Δμ = -0.2035 - 0.1681 z + 0.0471 z²  0.0083        0.0288
linear ln(1+z) Δμ = -0.2299 - 0.1135 ln(1+z)      0.0176          0.0551
```

The quadratic-in-z fit is closest (RMS 0.008 mag); the difference is
monotonic, saturating toward a roughly constant offset at high z. This
is consistent with the two cosmologies producing different h(z) shapes
that integrate to a smoothly diverging luminosity distance.

## Pantheon+ residuals under each framework

After per-catalog ΔM marginalization with diagonal weights:

```text
framework                        RMS residual (mag)
Model-A V_3 numerical KG               0.158
LCDM                                   0.157
```

Both frameworks reproduce the Pantheon+ distance modulus to comparable
RMS once the per-framework ΔM offset is absorbed.

## Files

- `scripts/G8g_v3_distance_prediction_comparison.py`
- `reports/G8g/v3_lcdm_distance_comparison.png`
- `reports/G8g/v3_lcdm_distance_difference.csv`

## Notes

- Script 38 (the original "derived adjustment" version using V_1
  cosmology and an LCDM-bridging framing) is now headed as superseded
  and points at this script.
- The Δμ shape is reported here as a model-comparison observable, not
  as Model-A's "adjustment to bridge to LCDM." The framing change
  follows the two-tables commitment.
