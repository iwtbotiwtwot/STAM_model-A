# G103 — Greybody / scattering diagnostic from STAM throat

## Purpose

G103 computes real-frequency axial-wave scattering for GR, STAM metric-only proxy, and STAM action-aware potentials. It probes the near-horizon tortoise-distance enhancement from G90 without relying on complex-frequency QNM extraction.

## Method

Solve:

```text
ψ''(r*) + [ω² - V(r*)]ψ = 0
```

with a unit transmitted ingoing wave at the horizon and incident/reflected decomposition at infinity. Transmission is `T = 1/|A_in|²`.

## Representative results

| omega M | T_GR | T_proxy | T_action | proxy/GR | action/GR | T+R action |
|---:|---:|---:|---:|---:|---:|---:|
| 0.100 | 3.0443e-06 | 2.7680e-07 | 1.5197e-08 | 0.091 | 0.005 | 1.00000 |
| 0.200 | 8.6817e-04 | 5.8285e-04 | 1.7709e-04 | 0.671 | 0.204 | 1.00000 |
| 0.400 | 6.8981e-01 | 7.0527e-01 | 5.2700e-01 | 1.022 | 0.764 | 1.00000 |
| 0.800 | 1.0000e+00 | 9.9999e-01 | 9.9906e-01 | 1.000 | 0.999 | 1.00000 |
| 1.000 | 1.0000e+00 | 1.0000e+00 | 9.9980e-01 | 1.000 | 1.000 | 1.00000 |

## Interpretation

If the action-aware transmission differs systematically from GR while unitarity remains close to one, this is a cleaner observable-channel signature of the STAM throat than the early low-l QNM attempts. The result should still be treated as a diagnostic until the full constrained-Σ axial perturbation equation is derived.

## Files

- CSV: `/mnt/data/G103_run/results/G103_greybody_scattering.csv`
- Plot: `/mnt/data/G103_run/plots/G103_transmission_curves.png`
- Plot: `/mnt/data/G103_run/plots/G103_transmission_ratios.png`
- Plot: `/mnt/data/G103_run/plots/G103_potentials.png`
