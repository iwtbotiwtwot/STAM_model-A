# G94 — QNM robustness sweep
**Purpose.** Test whether the G93 action-aware low-ell QNM shift is robust under fit-window, grid, throat-depth, outer-boundary, pulse, observer, and f-normalization choices.

Mode: quick sweep.

## Summary statistics for calibrated PASS rows
- PASS rows: 22 / 25
- Schwarzschild calibration error: mean 0.731%, range 0.244%–0.908%
- Proxy shift: mean 1.237%, range 1.132%–1.445%
- Action-aware shift: mean 5.384%, range 5.071%–5.632%
- Action-aware shift std: 0.118%

## Results table

| Group | Variant | Calib % | ω_GR | Proxy shift % | Action shift % | Status |
|---|---|---:|---:|---:|---:|---|
| baseline | baseline | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| fit_window | fit_90_230 | 8.873 | 0.339645-0.087009i | 0.972 | 5.879 | CALIB_FAIL |
| fit_window | fit_100_250 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| fit_window | fit_110_270 | 0.519 | 0.372243-0.090354i | 1.423 | 5.632 | PASS |
| fit_window | fit_120_290 | 0.244 | 0.373612-0.088027i | 1.445 | 5.372 | PASS |
| throat_depth | rsmin_-300 | 0.796 | 0.370633-0.089287i | 1.195 | 5.349 | PASS |
| throat_depth | rsmin_-500 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| throat_depth | rsmin_-700 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| outer_boundary | rsmax_300 | 0.797 | 0.370628-0.089277i | 1.195 | 5.351 | PASS |
| outer_boundary | rsmax_400 | 0.797 | 0.370628-0.089277i | 1.195 | 5.351 | PASS |
| outer_boundary | rsmax_500 | 0.797 | 0.370628-0.089277i | 1.195 | 5.351 | PASS |
| grid_resolution | Nstam_8001 | 0.796 | 0.370633-0.089287i | 1.195 | 5.348 | PASS |
| grid_resolution | Nstam_10001 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| grid_resolution | Nstam_12001 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| pulse | pulse_20_sig_3 | 0.373 | 0.372499-0.089788i | 1.417 | 5.575 | PASS |
| pulse | pulse_30_sig_3 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| pulse | pulse_40_sig_3 | 9.626 | 0.336698-0.089011i | 0.931 | 5.891 | CALIB_FAIL |
| pulse | pulse_30_sig_2 | 0.908 | 0.370223-0.088442i | 1.257 | 5.589 | PASS |
| pulse | pulse_30_sig_4 | 0.836 | 0.373022-0.092106i | 1.132 | 5.071 | PASS |
| observer | obs_40 | 0.459 | 0.372353-0.090131i | 1.423 | 5.605 | PASS |
| observer | obs_50 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| observer | obs_70 | 23.703 | 0.282948-0.081305i | 0.193 | 2.068 | CALIB_FAIL |
| f_normalization | fscale_0.1 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| f_normalization | fscale_1 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |
| f_normalization | fscale_10 | 0.796 | 0.370633-0.089287i | 1.196 | 5.350 | PASS |

## Interpretation guide
A robust action-aware QNM estimate should keep Schwarzschild calibration below 1% while leaving the action-aware shift in a narrow band under numerical/extraction variations. f-normalization should not change the action-aware result because (sqrt(f))''/sqrt(f) is invariant under constant rescaling of f.

## Files
- `scripts/G94_QNM_robustness_sweep.py`
- `results/G94_QNM_robustness_sweep.csv`
- `plots/G94_QNM_robustness_sweep.png`
