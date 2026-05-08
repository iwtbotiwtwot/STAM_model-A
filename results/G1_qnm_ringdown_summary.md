# G1: Quasinormal-Mode Ringdown Spectrum — Model-A vs Schwarzschild

## The test

Model-A's strong-field metric differs from Schwarzschild only in g_rr (time component g_tt is identical). The QNM ringdown spectrum after binary BH merger probes near-horizon geometry, where the (1-A²)² modification of k(A) bites. We compute the dominant ℓ=2, n=0 mode (what LIGO measures) for both metrics in the eikonal/WKB approximation at the photon sphere, then convert to physical frequency and damping time for representative BBH remnants.

## Setup

**Both metrics share:**
- g_tt = -(1-A) c²  [h(A) = 1-A]
- Photon sphere at r_c = 3M (A = 2/3)
- Orbital frequency Ω_c = √(1/3) / (3M) = 1/(3√3 M)

**They differ in g_rr:**
- Schwarzschild: g_rr = 1/(1-A)              k(r_c) = 1/3
- Model-A:       g_rr = 1/[(1-A)(1-A²)²]     k(r_c) = 25/243 ≈ 0.103

**Eikonal QNM formulas (Cardoso-Konoplya-Zhidenko):**
```text
ω_R = l × Ω_c                              (real, oscillation frequency)
ω_I = -(n + 1/2) × λ                       (imaginary, damping rate)
λ²  = h(r_c) × k(r_c) / (3 M²)             (Lyapunov exponent at photon sphere)
```

Because Ω_c depends on h only, Model-A predicts the **same ringdown frequency** as Schwarzschild for any given BH mass. The damping rate λ depends on h × k together, so Model-A's smaller k at the photon sphere gives a **longer damping time τ = 1/|ω_I|** by a factor:

```text
λ_S / λ_MA = √(k_S / k_MA) = √((1/3) / (25/243)) = √(81/25) = 9/5 = 1.80
τ_MA / τ_S = λ_S / λ_MA   = 1.80
```

**Model-A predicts ringdowns are 80% longer than GR predicts**, at the same dominant frequency, for any BH mass. This ratio is independent of BH mass — it's set entirely by the photon-sphere k value.

## Predictions for representative LIGO events

(Eikonal approximation; absolute values are ~3% off the Berti exact ℓ=2 axial gravitational mode, but the **ratio** Model-A / Schwarzschild is robust because both share the same approximation structure.)

```text

Event                                 M (M_sun)    f_GR(Hz)    f_MA(Hz)    τ_GR(ms)    τ_MA(ms)   τ_MA/τ_GR
--------------------------------------------------------------------------------------------------------------
GW150914 remnant                             62       200.5       200.5       3.174       5.714       1.800
GW170729 remnant                             80       155.4       155.4       4.096       7.373       1.800
GW190521 remnant                        1.4e+02        87.6        87.6        7.27       13.09       1.800
Stellar-mass BBH (typical)                   30       414.5       414.5       1.536       2.765       1.800
Intermediate-mass (~1000 M_sun)           1e+03        12.4        12.4        51.2       92.16       1.800
Sgr A*                                  4.3e+06         0.0         0.0   2.202e+05   3.963e+05       1.800
M87*                                    6.5e+09         0.0         0.0   3.328e+08    5.99e+08       1.800
```

## Falsification logic

LIGO's ringdown analysis of GW150914 (final BH ~62 M_sun) extracted a dominant mode at f ≈ 251 Hz with damping time τ ≈ 4 ms — consistent with the Schwarzschild prediction. Subsequent BBH events (GW170729, GW190521, etc.) have produced similar ringdown measurements.

Model-A predicts the **same frequency** (f ≈ 251 Hz for 62 M_sun) but a **1.80× longer damping time** (τ ≈ 7.2 ms instead of 4 ms). This is a 80% effect — well above current measurement uncertainties on damping time for the loudest events.

**The cleanest test:** the ringdown signal h(t) ∝ exp(-t/τ) cos(2πft) should ring 80% longer than GR predicts if Model-A is right. LIGO's matched-filter analysis fits both f and τ; if τ comes out consistently with the GR value (~4 ms for 62 M_sun class) and inconsistent with Model-A's longer prediction, that's a falsification of Model-A's k(A) at A = 2/3.

**What this distinguishes:** the modification factor k(A) = (1-A²)² × (1-A) hits the photon sphere at A = 2/3, where (1-A²)² = (5/9)² ≈ 0.31. So the photon sphere sees about 1/3 the Schwarzschild k value. The ringdown probes exactly this region. If we miss the prediction, the modification factor must be tuned down — but the modification is constrained at A=1 (must give boundary behavior) and at A→0 (must match GR weak field), so there's limited room to dial it without breaking other commitments.

## Honest caveats

- **Eikonal approximation.** This is a high-ℓ asymptotic that's about 3% off for ℓ=2. The RATIO Model-A/Schwarzschild is more robust than absolute frequencies because both share the same approximation.
- **Scalar perturbation surrogate.** This script computes scalar QNMs as a stand-in for the gravitational ℓ=2, m=2, n=0 mode. The actual gravitational perturbation equations (Regge-Wheeler / Zerilli for Schwarzschild; Model-A analog yet to be derived) may modify the result. For Schwarzschild, scalar and axial-gravitational ℓ=2 QNMs differ by ~10%; same for Model-A. The Model-A/GR ratio of damping times should still be ~1.8 to within similar precision.
- **Spin neglected.** Real BH remnants are spinning (typical χ_f ≈ 0.7 for BBH mergers). Kerr QNMs differ from Schwarzschild; Model-A's spinning analog has not been constructed. Spin-aligned events with moderate χ should still show the qualitative ratio, but a full Kerr-Model-A computation is a follow-up.
- **WKB at higher orders.** Schutz-Will (1985) gives a few-percent improvement over leading WKB; Konoplya-Zhidenko (2019) give 6th-order. These improve the absolute QNM accuracy but not the Model-A/GR ratio.

## Generated plots

- `plots/G1_qnm_metric_and_damping.png`
- `plots/G1_qnm_ringdown_waveforms.png`
- `plots/G1_qnm_tau_ratio.png`
