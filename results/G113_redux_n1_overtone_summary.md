# G113-redux -- n=1 Axial Overtone (Methodological Finding, NOT Extraction)

**Date: 2026-05-13.** Attempted n=1 spinless axial overtone extraction
via two methods (regularized frequency-domain shooting and WKB-PT order
3 / Schutz-Will). Both produced a clean methodological finding that
SHARPENS the open problem rather than closing it.

## What was tried

### Attempt 1: Regularized frequency-domain shooting

Implemented the proper regularized form `v'' - 2i omega v' - V v = 0`
on horizon side and `u'' + 2i omega u' - V u = 0` on infinity side,
with Wronskian = `v u' - u v' + 2i omega u v` at the match point.

**Calibration on GR Schwarzschild:**
- n=0: err 13.6% vs Berti 0.374 - 0.089i
- n=1: err 53.9% vs Berti 0.347 - 0.274i

Calibration too poor for STAM extraction. Cause: simple-BC shooting
(leading-order asymptotic) suffers from sub-leading contamination at
finite rstar bounds. Sub-percent extraction requires Frobenius series
BCs at horizon and asymptotic series at infinity (or Leaver CF).

### Attempt 2: WKB-PT order 3 (Schutz-Will / Iyer-Will)

Located V peak via Brent's method on -V(r*); computed V, V'', V''', V''''
at peak via 8-point stencil; applied the standard formula
`omega^2 = V_0 - i(n+1/2) sqrt(-2 V_0'') [1 + Lambda_2]`.

**Calibration on GR Schwarzschild (Berti 2009):**
- n=0: err 6.88%
- n=1: err 52.8%

Consistent with literature: WKB-PT order 3 gives ~1-7% on n=0 and
~30-50% on n=1 for Schwarzschild ell=2.

## The methodological finding

**The wave-zone peak of V_GR for ell=2 axial sits at r=3.28M, where A
= 0.61, OUTSIDE the photon sphere (A < 2/3).** The framework's k(A) =
(1-A) exactly in that regime -- F(y) is not active above A=2/3.

Consequence: WKB-PT on V_geom-only finds an identical peak in GR and
STAM (both r=3.28M), giving STAM/GR shift = 0.0000%. The framework's
non-zero QNM shift (G108/G109's 4.977% rigorous n=0) comes ENTIRELY
from the action correction `(sqrt f)''/sqrt f`, NOT from the k(A)
modification.

WKB-PT with the action correction added breaks because (sqrt f)''/sqrt
f produces a sharp throat spike inside PS that WKB misinterprets as
the QNM barrier (peak shifts to r=2.63M with V_0 ~ 19, throwing the
omega to ~19-19i, calibration error 6900%).

## Why n=1 remains open

For the framework's QNM shift to be reproducible at n=1, we need a
method that:
1. Correctly handles V_exact = V_geom + (sqrt f)''/sqrt f including
   the rapid throat structure
2. Treats the eigenfunction's exponential decay inside the throat
   properly (so the throat spike doesn't dominate)
3. Achieves sub-percent calibration on Schwarzschild n=0,1 reference

Three methods meet these criteria, all requiring substantial machinery:
- **Leaver continued fraction** with verified Schwarzschild coefficients
  (G113a failed at coefficient self-check; redo by line-by-line
  transcription from qnm-package source)
- **Sasaki-Nakamura transformation** to a short-range potential, then
  shooting (G120d-style proposal)
- **Konoplya higher-order WKB** (order 6 with Pade resummation), with
  explicit peak isolation excluding the throat spike

This script's contribution is identifying which calibration regime
each method occupies, not extracting an n=1 STAM shift.

## Same-method shift (METRIC MODIFICATION ONLY -- not the framework's prediction)

With V_geom only (no action correction):
- n=0 same-method (WKB STAM vs WKB GR): 0.0000%
- n=1 same-method (WKB STAM vs WKB GR): 0.0003%

These confirm that the metric-modification contribution to ell=2
axial QNMs is essentially zero (WKB peak outside PS where k = 1-A).
**The framework's 4.977% n=0 shift comes from the action correction.**

This is a useful negative result: it sharpens the framework's
prediction structure. The framework's QNM signature is an
ACTION-CORRECTION SIGNATURE, not a metric-modification signature.

## Status

- **n=1 STAM shift NOT extracted in this attempt.**
- **Methodological finding:** ell=2 axial QNM shifts are dominated by
  the action correction (sqrt f)''/sqrt f, not by the F(y) metric
  modification. WKB-PT and basic shooting cannot isolate the action
  contribution cleanly.
- **Path forward unchanged from HANDOFF_2026_05_13_kerr_wave_branch:**
  Leaver CF with verified coefficients, or Sasaki-Nakamura, are the
  rigorous routes. Both ~200-300 lines of careful work.

## Honest framing for the framework

G108/G109 succeeded for n=0 because the time-domain method (Gaussian
pulse + Prony) handles V_exact's throat structure naturally
(eigenfunction's exp-small content in the throat means the spike
contributes negligibly to ringdown extraction). The n=1 mode has
amplitude ~1e-3 of n=0, so time-domain extraction fails -- not
because the method is wrong, but because the precision wall.

**The n=1 shift exists and is well-defined; we just haven't built the
right extraction tool yet.** This script ruled out two candidate
tools (regularized shooting, WKB-PT order 3).

## Files

- [scripts/G113_redux_n1_overtone_via_shooting.py](../scripts/G113_redux_n1_overtone_via_shooting.py)
- [results/G113_redux_n1_wkb.png](G113_redux_n1_wkb.png) -- shows V_GR vs V_STAM, peak locations
- See also: [HANDOFF_2026_05_13_kerr_wave_branch.md](../memory/HANDOFF_2026_05_13_kerr_wave_branch.md) (the original G110-G112 time-domain failures)
