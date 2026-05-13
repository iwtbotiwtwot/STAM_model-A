# G86 redux — action-aware photon-sphere potential structure

This reruns the original G86 structural question after G90-G94.
It compares GR, STAM metric-proxy, and STAM action-aware axial diagnostic potentials near the photon sphere.

## Setup

```text
V_action = V_proxy + (sqrt(f))'' / sqrt(f)
```

where primes are STAM tortoise-coordinate derivatives and `f(A)` is reconstructed locally from the G70/G75 relation with `f(2/3)=1`.

## First nonzero derivative differences at r = 3M, ell = 2

| comparison | ordinary r derivative | tortoise derivative |
|---|---:|---:|
| proxy − GR | order 4: `-640/729` | order 4: `-640/59049` |
| action − GR | order 2: `80/81` | order 2: `80/729` |

## Interpretation

The metric-proxy and action-aware potentials both preserve the leading photon-sphere/eikonal match.
The finite-order derivative where the first mismatch appears explains why low-ell, non-eikonal QNM modes can shift even when the leading eikonal photon-orbit result is GR-exact.

The action-aware term tests whether the non-minimal coupling `f(Σ)R` introduces lower-order or larger photon-region structure than the metric proxy alone.
This script is structural only. It does not extract a QNM frequency; G91-G94 handle calibrated time-domain extraction.

## Files

- CSV: `C:\Users\drwho\OneDrive\Documents\STAM_model-A\STAM_model-A\results\G86_redux_action_aware_potential_derivatives.csv`
- Plot: `C:\Users\drwho\OneDrive\Documents\STAM_model-A\STAM_model-A\plots\G86_redux_action_aware_potential.png`
