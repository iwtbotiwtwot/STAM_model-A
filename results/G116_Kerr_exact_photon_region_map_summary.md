# G116 — Kerr exact photon-region STAM map

G116 converts the G85 exact Kerr photon-region surface into a reusable numerical map for STAM Kerr observables.

## Exact photon-region formulas

```text
lambda(r_p) = -(r_p^3 - 3r_p^2 + a^2 r_p + a^2) / [a(r_p - 1)]
eta(r_p) = r_p^2 [4a^2 Delta - (r_p^2 - 3r_p + 2a^2)^2] / [a^2 (r_p - 1)^2]
Theta(theta;r_p,a) = eta sin^2(theta) - cos^2(theta)(lambda^2 - a^2 sin^2(theta))
r_pr_exact(theta,a) = root of Theta=0 in [r_ph_prograde, r_polar]
```

## G84 interpolation error summary

| a | r_plus | r_eq | r_polar | max abs Δr | max rel Δr | max abs ΔΣ | max rel ΔΣ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 2.000000 | 3.000000 | 3.000000 | 0.000000 | 0.00% | 0.000000 | 0.00% |
| 0.30 | 1.953939 | 2.630026 | 2.959311 | 0.088762 | 3.19% | 0.064425 | 3.02% |
| 0.50 | 1.866025 | 2.347296 | 2.883218 | 0.155502 | 6.02% | 0.119034 | 5.30% |
| 0.70 | 1.714143 | 2.013334 | 2.757907 | 0.239608 | 10.37% | 0.191605 | 8.00% |
| 0.90 | 1.435890 | 1.557855 | 2.559997 | 0.381520 | 20.44% | 0.307188 | 11.77% |
| 0.99 | 1.141067 | 1.167642 | 2.430983 | 0.558077 | 40.31% | 0.418483 | 14.92% |

## Interpretation

G116 confirms the G85 result in map form: G84's sin^2(theta) interpolation is useful for intuition but introduces sizeable high-spin, intermediate-latitude errors. Future Kerr observable calculations should use r_pr_exact(theta,a) and Sigma_ph_exact(theta,a).

The structural STAM claims are unchanged: ghost-freedom, cubic horizon vanishing, eikonal recovery, and lambda_1/lambda_2 closure depend on the existence and topology of the photon-region surface. G116 supplies the precision geometry needed for EMRI/EHT/ringdown observables.

## Files

- Map CSV: `/mnt/data/G116_run/results/G116_Kerr_exact_photon_region_map.csv`
- Summary CSV: `/mnt/data/G116_run/results/G116_Kerr_exact_photon_region_summary_table.csv`
- Plot: `/mnt/data/G116_run/plots/G116_rpr_exact_vs_theta.png`
- Plot: `/mnt/data/G116_run/plots/G116_shell_width_vs_theta.png`
- Plot: `/mnt/data/G116_run/plots/G116_g84_error_vs_theta.png`
- Plot: `/mnt/data/G116_run/plots/G116_sigma_ph_vs_theta.png`
- Plot: `/mnt/data/G116_run/plots/G116_WK_factor_heatmap.png`
