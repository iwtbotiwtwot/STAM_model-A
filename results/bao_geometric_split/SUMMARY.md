# STAM Model-A BAO geometric split test results

Locked b: `354.95`

L: `13387.090426252385` Mly

One nuisance scale `r_d_eff` was fitted for each mapping and scenario.

Full BAO covariance matrices were not used in this first-pass compressed-data script.

## combined_all_z

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | geo_over_1pz | 9.235 | 107.5 | 0.8677 | 0.609 |
| 2 | adj_over_1pz | 12.25 | 138 | 0.9993 | 0.06783 |
| 3 | geo_raw | 171 | 211.7 | 3.734 | 1.442 |
| 4 | adj_raw | 386.2 | 301.7 | 5.611 | 3.218 |
| 5 | geo_times_1pz | 848.9 | 556.1 | 8.319 | 6.377 |

## combined_z_le_1p6

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | adj_over_1pz | 2.481 | 135 | 0.3998 | -0.04618 |
| 2 | geo_over_1pz | 9.382 | 108.2 | 0.7775 | 0.2157 |
| 3 | geo_raw | 58.06 | 189.5 | 1.934 | 0.5191 |
| 4 | adj_raw | 139 | 245.1 | 2.993 | 1.133 |
| 5 | geo_times_1pz | 334.4 | 365.6 | 4.642 | 2.495 |

## desi_dr1_bao_all_z

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | geo_over_1pz | 7.117 | 105.5 | 0.8257 | 0.3391 |
| 2 | adj_over_1pz | 12.47 | 139.1 | 1.093 | 0.006107 |
| 3 | geo_raw | 173.4 | 219.9 | 4.075 | 1.928 |
| 4 | adj_raw | 392.8 | 315.9 | 6.135 | 3.688 |
| 5 | geo_times_1pz | 878.5 | 581.6 | 9.174 | 7.275 |

## desi_dr1_bao_z_le_1p6

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | adj_over_1pz | 2.995 | 136.1 | 0.4703 | -0.1658 |
| 2 | geo_over_1pz | 5.815 | 106.4 | 0.6553 | 0.06471 |
| 3 | geo_raw | 45.03 | 196.4 | 1.824 | 0.05161 |
| 4 | adj_raw | 101.3 | 256.3 | 2.735 | 0.3778 |
| 5 | geo_times_1pz | 231.9 | 381.8 | 4.138 | 1.128 |

## sdss_dr17_bao_only_all_z

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | geo_over_1pz | 10.93 | 109.7 | 0.7979 | 0.8383 |
| 2 | adj_over_1pz | 14.54 | 136.7 | 0.9203 | 0.1688 |
| 3 | geo_raw | 198.4 | 201.2 | 3.4 | 2.035 |
| 4 | adj_raw | 456.9 | 282.2 | 5.16 | 3.506 |
| 5 | geo_times_1pz | 1015 | 518.7 | 7.692 | 5.844 |

## sdss_dr17_bao_only_z_le_1p6

| rank | mapping | chi2/dof | r_d_eff Mpc | weighted RMSE | median residual |
|---:|---|---:|---:|---:|---:|
| 1 | adj_over_1pz | 2.176 | 133.8 | 0.3098 | 0.01763 |
| 2 | geo_over_1pz | 13.59 | 110.1 | 0.7741 | 0.3568 |
| 3 | geo_raw | 76.58 | 181 | 1.838 | 0.6143 |
| 4 | adj_raw | 201.1 | 230.9 | 2.978 | 1.237 |
| 5 | geo_times_1pz | 527.6 | 343.7 | 4.824 | 2.808 |

## Plots

- `results/bao_geometric_split/bao_mapping_chi2_dof.png`
- `results/bao_geometric_split/combined_all_z_residuals_geo_vs_adj.png`
- `results/bao_geometric_split/combined_all_z_observed_vs_stam.png`