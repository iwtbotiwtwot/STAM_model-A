# G119b -- Kerr-base spin sweep redux (WKB stepping stone)

Stepping stone before full G120 Teukolsky.  Poschl-Teller / WKB
frequency-domain estimate on V_base, applied on the G119 Kerr base.

## V_base calibration audit (vs published Kerr l=m=2 n=0)
| a | WKB Re(om) on V_base | published Kerr Re(om) | WKB vs pub |
|---:|---:|---:|---:|
| 0.00 | 0.3783 | 0.3737 | 1.26% |
| 0.10 | 0.3787 | 0.3814 | 0.88% |
| 0.20 | 0.3802 | 0.3898 | 2.45% |
| 0.30 | 0.3827 | 0.3991 | 4.04% |
| 0.90 | 0.4336 | 0.4935 | 12.08% |
| 0.95 | 0.4443 | 0.5108 | 12.99% |
| 0.99 | 0.4552 | 0.5325 | 14.65% |

Reading: WKB-on-V_base differs from published Kerr by O(few %)
at low spin (standard WKB error) and grows at high spin (because
V_base lacks Kerr-specific m-coupling).  This audits the gap
between V_base and rigorous Kerr.

## WKB shift vs G118v2 time-domain shift on the same V_base
| a | G119b WKB shift | G118v2 TD shift | WKB - TD |
|---:|---:|---:|---:|
| 0.00 | 37.1003% | 4.5703% | +32.5300 |
| 0.10 | 37.3396% | 5.2665% | +32.0731 |
| 0.20 | 33.9051% | 5.8360% | +28.0691 |
| 0.30 | 30.7813% | 6.1631% | +24.6182 |
| 0.90 | 0.0000% | 14.0815% | -14.0815 |
| 0.95 | 0.0000% | 10.7923% | -10.7923 |
| 0.99 | 64.6839% | 12.1178% | +52.5661 |

If WKB and TD agree, the shift is robust to extraction method.
If they disagree, the difference reveals method sensitivity in
the time-domain fit window or higher-overtone contamination.

## Caveats
- V_base is Kerr-Delta-aware Schwarzschild-RW analog (G118v2 form).
  NOT the rigorous Kerr axial potential.
- WKB-PT is 0th-order WKB; absolute frequencies have few-% intrinsic
  error even at a = 0.
- This is NOT a final Kerr QNM prediction.  See G120 spec below.

## Next step:  G120 full Detweiler / Sasaki-Nakamura
G120 must:
1. Build the spin-weighted spheroidal angular eigenvalue solver A_lm(a*omega).
2. Build the Sasaki-Nakamura or Detweiler radial equation for s=-2.
3. Calibrate GR Kerr QNM to <0.5% at a in [0, 0.5, 0.9].
4. Insert G119's Delta_STAM only after calibration passes.
5. Track the l=m=2 n=0 mode branch across the spin sweep.
6. Report status flags per spin.

## Files
- [scripts/G119b_kerr_base_spin_sweep_redux.py](../scripts/G119b_kerr_base_spin_sweep_redux.py)
- [plots/G119b_kerr_base_spin_sweep_redux.png](../plots/G119b_kerr_base_spin_sweep_redux.png)
