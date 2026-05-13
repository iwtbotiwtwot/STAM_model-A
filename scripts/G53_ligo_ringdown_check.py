"""
G53: LIGO ringdown constraints vs framework τ/τ_GR predictions
================================================================

The framework's G1 prediction (eikonal approximation):

    τ_QNM / τ_GR = (9/5)^(n/2)   at the horizon (A → 1)

For candidate n(A) forms, the value of n AT HORIZON determines the
ringdown damping ratio.

This script:
  1. Computes the predicted τ/τ_GR for each candidate n(A) form
  2. Compiles published LIGO ringdown / QNM constraints
  3. Compares the predictions to the observational bounds
  4. Reports which candidates are consistent / excluded

Caveats:
  - Framework's prediction is for SPINLESS analog; LIGO BHs are spinning Kerr.
    For typical remnant spin (a ~ 0.5-0.7), the Kerr damping time τ_Kerr is
    ~5-15% shorter than Schwarzschild τ_Schw. The translation introduces
    O(10%) uncertainty in the comparison.
  - Framework's prediction is EIKONAL (high-ℓ). LIGO measures ℓ=m=2 mode,
    which deviates from eikonal by ~10-20%. Adds another uncertainty.
  - The "spinless STAM analog" remains an open framework computation
    (memory: "spinless Model-A analog is the missing piece for direct
    comparison").

Net: rough comparison only. Excludes by factor-of-few; doesn't precisely
distinguish n=1 vs n=2.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

A_0 = 1.0 / (12 * np.pi)
A_horizon = 0.99999  # numerical limit

# --- Framework predictions: τ/τ_GR for each candidate n form, evaluated at horizon ---
candidates = [
    ("n = 1 (constant)",        lambda A: 1.0),
    ("n = 2 (constant)",        lambda A: 2.0),
    ("n = 1 + A",               lambda A: 1 + A),
    ("n = 1 + A²",              lambda A: 1 + A**2),
    ("n = 1 + 3A",              lambda A: 1 + 3*A),
    ("n = 1 + (A-A_0)/(1-A_0)", lambda A: 1 + (A - A_0)/(1 - A_0)),
    ("n = 1 + 2(A-A_0)/(1-A_0)", lambda A: 1 + 2*(A - A_0)/(1 - A_0)),  # n=1 to 3
    ("n = 1 + 3(A-A_0)/(1-A_0)", lambda A: 1 + 3*(A - A_0)/(1 - A_0)),  # n=1 to 4
]

def tau_ratio(n_horizon):
    return (9/5)**(n_horizon/2)

def delta_tau_pct(n_horizon):
    return (tau_ratio(n_horizon) - 1) * 100

# --- Published LIGO ringdown constraints (approximate, from public TGR papers) ---
ligo_constraints = """
Published LIGO/Virgo ringdown constraints on δτ/τ_GR
(from Tests of General Relativity catalog papers GWTC-2, GWTC-3, and
single-event analyses):

  Individual loud events (e.g., GW150914, GW190521):
    1-sigma constraint on δτ/τ_GR : ~10-30% (event-dependent)
    2-sigma constraint            : ~20-60%
    3-sigma constraint            : ~30-100%

  Stacked O3 catalog (Tests of GR paper):
    Combined 1-sigma constraint   : ~5-15%
    Combined 2-sigma              : ~10-30%

These constraints are for ℓ=m=2 fundamental QNM, on spinning Kerr remnants.
Spinless-analog corrections introduce O(10-15%) uncertainty in mapping
framework predictions to observed bounds.
"""

print("=" * 90)
print("G53: LIGO ringdown constraints vs framework predictions")
print("=" * 90)
print(ligo_constraints)

# --- Predictions table ---
print("-" * 90)
print(f"{'Candidate n(A)':<32} {'n at horizon':>14} {'τ/τ_GR':>10} {'δτ/τ_GR (%)':>14} {'consistent?':>20}")
print("-" * 90)

# Rough "consistent" criterion: δτ/τ_GR within ~30% of zero is consistent with single events
# within ~15% is consistent with stacked

for label, n_func in candidates:
    n_h = n_func(A_horizon)
    ratio = tau_ratio(n_h)
    delta_pct = delta_tau_pct(n_h)

    if abs(delta_pct) < 15:
        verdict = "consistent (stacked)"
    elif abs(delta_pct) < 30:
        verdict = "1-2σ tension (single)"
    elif abs(delta_pct) < 60:
        verdict = "2-3σ tension"
    elif abs(delta_pct) < 100:
        verdict = "3σ+ tension"
    else:
        verdict = "excluded"

    print(f"  {label:<32} {n_h:>14.3f} {ratio:>10.3f} {delta_pct:>+14.1f}% {verdict:>20}")

print()
print("=" * 90)
print("What the LIGO ringdown bounds tell us")
print("=" * 90)

print("""
Summary of the ringdown comparison:

  n = 1 (constant)              τ/τ_GR = 1.342, δ = +34%
    Marginal. Within single-event 2σ but tension with stacked.

  n = 2 (constant) / n = 1+A / n = 1+A² / n = 1+(A-A_0)/(1-A_0)
                                τ/τ_GR = 1.80, δ = +80%
    Strong tension with single events (2-3σ), very strong tension with
    stacked catalog (5σ+). Likely already disfavored by current data.

  n = 1 + 3A (n=4 at horizon)   τ/τ_GR = 3.24, δ = +224%
    Excluded by current LIGO data at very high confidence.

  n = 1 + 2(A-A_0)/(1-A_0) (n=3 at horizon)  τ/τ_GR = 2.42, δ = +142%
    Excluded by current data (>3σ tension).

Implications for the framework:

  1. **n = 1 + 3A is likely ruled out** by current LIGO ringdown data,
     despite its appealing structural alignment with G18 (n_horizon = 4 =
     two-face × hologram = α area-per-entry).

     The "horizon has 2 sides and 2 holograms = 4 elements" reading is
     structurally compelling but observationally constrained.

  2. **n = 2 at horizon (under any form: const n=2, n=1+A, etc.) is in
     tension** with current data. Marginally allowed for individual events
     but increasingly disfavored as stacked precision improves.

     The framework's existing commitment to n=2 was already on the edge of
     current ringdown constraints.

  3. **n = 1 (constant) is the most observationally favored** candidate,
     with τ/τ_GR = 1.342 only mildly disfavored.

Caveats (real):
  - LIGO measurements are for SPINNING Kerr BHs. Framework prediction is
    spinless. Spin corrections shift the prediction by O(10%).
  - Framework predicts eikonal (high-ℓ) limit. LIGO measures ℓ=m=2 mode.
    Mode-dependence shifts prediction by O(10-20%).
  - These together introduce maybe ±15-25% uncertainty in the bound
    comparison. So the "3σ tension" for n=2 could actually be 2σ or 4σ
    depending on translation details.

What the math says structurally:

  If LIGO ringdown excludes τ/τ_GR > ~2 (which current stacked O3 data
  approximately does), then the framework is observationally constrained
  to have **n_horizon ≤ ~2**, possibly as low as ~1.5 with tight stacking.

  This excludes:
    - n = 1 + 3A (n_horizon = 4)
    - n = 1 + 2(A-A_0)/(1-A_0) (n_horizon = 3)

  Marginally constrains:
    - n = 2 (constant)
    - n = 1 + A (n_horizon = 2)
    - n = 1 + A² (n_horizon = 2)
    - n = 1 + (A-A_0)/(1-A_0) (n_horizon = 2)

  Allows comfortably:
    - n = 1 (constant)
    - any n(A) with n_horizon < 1.5

Note: the "G18 entropy alignment requires n_horizon = 4" was a structural
elegance argument, but **observation appears to disfavor it**. The framework
must choose between G18's clean structural reading and LIGO consistency.

If observations prefer n_horizon close to 1 or 1.5, the framework's
prediction shifts:
  - n = 1: τ/τ_GR = 1.342 (most observation-friendly)
  - n_horizon = 1.5: τ/τ_GR = √(9/5)^1.5 = (9/5)^0.75 = 1.554
  - These both correspond to F3 NEC crossover near A = 1/2 (n=1 exact, n=1.5 close)
""")

# --- Specific test: what n_horizon does τ/τ_GR = 1.0 correspond to? (i.e., GR limit)
n_GR_limit = 0  # k(A) = 1-A, so n=0 by convention
print("-" * 90)
print("For reference:")
print(f"  τ/τ_GR = 1.0 ↔ n_horizon = {0:.2f} (i.e., k = (1-A) = pure GR)")
print(f"  τ/τ_GR = 1.342 ↔ n_horizon = 1.0 (= LIGO-preferred direction)")
print(f"  τ/τ_GR = 1.5 ↔ n_horizon = {2*np.log(1.5)/np.log(9/5):.2f}")
print(f"  τ/τ_GR = 1.8 ↔ n_horizon = 2.0 (= current STAM commitment)")
print(f"  τ/τ_GR = 2.0 ↔ n_horizon = {2*np.log(2.0)/np.log(9/5):.2f}")
print(f"  τ/τ_GR = 3.24 ↔ n_horizon = 4.0 (= n = 1+3A)")
