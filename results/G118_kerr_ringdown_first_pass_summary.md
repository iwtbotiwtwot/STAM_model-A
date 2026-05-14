# G118 -- first-pass Kerr ringdown STAM correction (spin sweep)

## Approximations (honest)
- GR base potential = Schwarzschild RW.  Spinning Kerr Teukolsky
  is the rigorous follow-up (G119).
- Tortoise built with the STAM-modified equatorial Kerr substance
  density A_K = 2Mr/(r^2 + a^2).
- f(A_K) constructed by applying the G70/G75 spinless matched ODE
  along the Kerr A_K profile (extension, not Kerr-matched derivation).
- r_ph_prograde(a) from G116 sets the inside-shell upper boundary.
- (sqrt f_K)'' / sqrt f_K added analytically via chain rule.

## What the script reports
STAM shift % as a function of spin.  This is NOT the Kerr QNM frequency;
it is the framework's f-correction magnitude as a function of spin.

## Spin sweep results
| a | r_ph_prograde | Sigma_ph_eq | omega_STAM | shift |
|---:|---:|---:|---|---:|
| 0.00 | 3.0000 | 2.0000 | (0.3883795254523546-0.08257606229586545j) | 4.9773% |
| 0.10 | 2.8822 | 2.0792 | (0.38983596562554085-0.08290661714894798j) | 5.3084% |
| 0.20 | 2.7592 | 2.1632 | (0.3904367384132818-0.08302597692977369j) | 5.4486% |
| 0.30 | 2.6300 | 2.2520 | (0.38932861842216626-0.08237905543633255j) | 5.2286% |
| 0.50 | 2.3473 | 2.4452 | (0.37962634705527765-0.081183674224276j) | 3.1757% |
| 0.70 | 2.0133 | 2.6587 | (0.36940673762476095-0.09020409414212482j) | 0.4012% |
| 0.90 | 1.5579 | 2.8877 | (0.37063061183558677-0.08928653085859573j) | 0.0000% |
| 0.99 | 1.1676 | 2.9895 | (0.37063061183558677-0.08928653085859573j) | 0.0000% |

## Reading
At a = 0 the calculation recovers G108's spinless 4.977% shift.
As spin increases, the inside-shell region shrinks (r_ph_prograde
decreases) and the f(A_K) profile sharpens; the table shows how the
shift trends with spin.  Increasing or decreasing trend is the
first-pass spin-dependence prediction.

## What's still needed
- G119 (rigorous Kerr Teukolsky on STAM background) to get actual Kerr QNMs.
- Polar Kerr (G120) along the G114-rigorous lines.
- Higher-mode and overtone Kerr extensions (G121+).

## Files
- [scripts/G118_kerr_ringdown_first_pass.py](../scripts/G118_kerr_ringdown_first_pass.py)
- [plots/G118_kerr_ringdown_first_pass.png](../plots/G118_kerr_ringdown_first_pass.png)
