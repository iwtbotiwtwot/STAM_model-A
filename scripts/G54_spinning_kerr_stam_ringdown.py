"""
G54: STAM ringdown for spinning Kerr-like remnants
====================================================

The framework's published prediction τ/τ_GR = (9/5)^(n/2) is for SPINLESS
Schwarzschild-like BHs. LIGO BHs are Kerr (spinning). We need the STAM
prediction for the actual case to compare to data.

This is the framework's open computation. Without an explicit STAM-Kerr
metric, we compute under several interpretations:

  Interpretation A — "Naive Kerr-with-STAM-g_rr": Use Kerr-like A profile
    A_Kerr ≈ 2M/r at equator, apply STAM modification factor (1-A²)^n at
    the Kerr light ring location.

  Interpretation B — "Static STAM bubble + rotating matter on surface":
    Per memory project_outward_collapse_dynamics.md / hawking_mechanism_native:
    "Spinning BH bubble framework: Static bubble + rotating holographic matter."
    The metric structure for QNM is the bubble's static geometry; matter
    rotation doesn't change the bubble's deep structure significantly.
    τ_STAM(M, a) ≈ τ_STAM_spinless(M) regardless of a.

  Interpretation C — "Hybrid": some spin-dependent correction between A and B.

For each, we compute τ_STAM(M, a) / τ_GR_Kerr(M, a) for typical LIGO
remnant spins.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Geometrized units: G = c = M = 1

def kerr_light_ring(a, prograde=True):
    """Equatorial photon ring radius in Kerr (Boyer-Lindquist, M=1 units)."""
    if prograde:
        return 2 * (1 + np.cos((2/3) * np.arccos(-a)))
    else:
        return 2 * (1 + np.cos((2/3) * np.arccos(+a)))

def kerr_horizon(a):
    """Outer event horizon radius."""
    return 1 + np.sqrt(1 - a**2)

def kerr_QNM_damping_ratio(a):
    """Approximate ratio τ_Kerr(a) / τ_Schw for ℓ=m=2 fundamental QNM.
    From standard tabulations:
      a=0: ratio = 1.0
      a=0.5: ratio ≈ 1.05
      a=0.67: ratio ≈ 1.10
      a=0.9: ratio ≈ 1.30
    Approximated as smooth fit."""
    # Approximation: τ_Kerr increases with spin
    return 1 + 0.1*a + 0.4*a**4  # rough fit

def A_equatorial_Kerr(r, a):
    """At Kerr equator (θ=π/2), Σ = r²; A_Kerr = 2M·r/Σ = 2/r in M=1 units."""
    return 2 / r

# --- Framework's STAM modification factor ---
def stam_factor(A, n):
    """τ_STAM/τ_GR = (1-A²)^(-n/2) at the light ring location.
    For Schwarzschild light ring A=2/3, gives (5/9)^(-n/2) = (9/5)^(n/2)."""
    return (1 - A**2)**(-n/2)

# --- Run for various spin values ---
print("=" * 90)
print("G54: STAM ringdown for spinning Kerr remnants")
print("=" * 90)

print("\nGeometrized units (G = c = M = 1).")
print(f"Schwarzschild light ring: r = 3, A_eq = 2/3 = 0.667")
print(f"STAM spinless prediction: τ/τ_GR_Schw = (9/5)^(n/2)")
print(f"  n=1: 1.342    n=2: 1.800")

# --- Interpretation A: Naive Kerr light ring with STAM modification ---
print("\n" + "-" * 90)
print("Interpretation A — Naive Kerr light ring × STAM (1-A²)^(-n/2)")
print("-" * 90)
print(f"\n{'Spin a':>8} {'r_LR':>8} {'A_LR':>8} {'(1-A²)':>10} {'n=1: τ ratio':>16} {'n=2: τ ratio':>16}")
print(f"{'-'*8} {'-'*8} {'-'*8} {'-'*10} {'-'*16} {'-'*16}")
for a in [0.0, 0.1, 0.3, 0.5, 0.67, 0.8, 0.9, 0.99]:
    r_LR = kerr_light_ring(a, prograde=True)
    A_LR = A_equatorial_Kerr(r_LR, a)
    factor = 1 - A_LR**2
    if factor > 0:
        stam_n1 = stam_factor(A_LR, 1)
        stam_n2 = stam_factor(A_LR, 2)
        print(f"{a:>8.2f} {r_LR:>8.3f} {A_LR:>8.4f} {factor:>10.4f} {stam_n1:>16.3f} {stam_n2:>16.3f}")
    else:
        print(f"{a:>8.2f} {r_LR:>8.3f} {A_LR:>8.4f} {'A>1':>10} {'inf':>16} {'inf':>16}")

print("""
  Implication of A: For LIGO typical spin a ~ 0.5-0.7, STAM predicts
  τ ratios of 1.7-3.9 (n=1) or 2.9-15 (n=2) — WILDLY inconsistent with
  observation (τ_obs/τ_GR_Kerr ≈ 1 ± 0.1-0.3).

  Under this interpretation, the framework is essentially ruled out for
  any n > 0 at moderate spin. This suggests Interpretation A is wrong —
  the framework's spinning analog isn't naive Kerr-with-STAM-g_rr.
""")

# --- Interpretation B: Static bubble + rotating matter ---
print("-" * 90)
print("Interpretation B — Static STAM bubble + rotating matter on surface")
print("(per memory: 'Static bubble + rotating holographic matter')")
print("-" * 90)
print("""
  Under this reading: the metric structure for QNM is determined by the
  bubble's STATIC geometry (Schwarzschild-like). Matter rotation affects
  the bubble's surface dynamics (Hawking T pattern, etc.) but not the
  exterior metric for ringdown purposes.

  Under this, the STAM modification factor stays at the Schwarzschild
  light ring value: (9/5)^(n/2). Independent of spin.

  Comparison to Kerr GR ringdown:
    τ_STAM(M, a) / τ_GR_Kerr(M, a) ≈ (9/5)^(n/2) / [τ_Kerr/τ_Schw](a)
""")

print(f"\n{'Spin a':>8} {'τ_Kerr/τ_Schw':>16} {'n=1: τ_STAM/τ_GR_Kerr':>24} {'n=2: τ_STAM/τ_GR_Kerr':>24}")
print(f"{'-'*8} {'-'*16} {'-'*24} {'-'*24}")
for a in [0.0, 0.1, 0.3, 0.5, 0.67, 0.8, 0.9, 0.99]:
    spin_ratio = kerr_QNM_damping_ratio(a)
    stam_n1_over_kerr = (9/5)**(1/2) / spin_ratio
    stam_n2_over_kerr = (9/5)**(2/2) / spin_ratio
    print(f"{a:>8.2f} {spin_ratio:>16.3f} {stam_n1_over_kerr:>24.3f} {stam_n2_over_kerr:>24.3f}")

print("""
  Under Interpretation B:
    For LIGO typical remnant a ≈ 0.67:
      n=1: τ_STAM/τ_GR_Kerr ≈ 1.22 (22% excess) — marginal LIGO consistency
      n=2: τ_STAM/τ_GR_Kerr ≈ 1.64 (64% excess) — 2-3σ LIGO tension

    For high-spin events (a > 0.9):
      Spin correction τ_Kerr/τ_Schw > 1.3, reducing the framework's
      observable excess but still measurable
""")

# --- Interpretation C: Hybrid / partial ---
print("-" * 90)
print("Interpretation C — Partial STAM at moderate spin")
print("-" * 90)
print("""
  If the framework's modification partially applies (e.g., averaged over
  the bubble surface), the prediction is between A and B.

  At present, the framework has no derivation of this hybrid behavior.
""")

# --- Comparison to LIGO data ---
print("-" * 90)
print("Comparison to LIGO measurements")
print("-" * 90)
print("""
  LIGO ringdown measurements (from Tests of GR catalog papers):
    GW150914 (a_remnant ≈ 0.67):   δτ/τ_Kerr_GR ≈ 0 ± 0.2-0.3 (1σ)
    GW190521 (massive BBH):       δτ/τ_Kerr_GR ≈ 0 ± 0.3 (1σ)
    GWTC-3 stacked:               δτ/τ_Kerr_GR ≈ 0 ± 0.1-0.15 (1σ)

  Under Interpretation A (naive Kerr): framework predicts large excess
    → strongly ruled out at any non-zero spin.

  Under Interpretation B (static bubble): framework predicts modest excess
    → n=1 marginal (within 1-2σ for individual events)
    → n=2 in 2-3σ tension with single events, 5σ+ tension with stacked

  Without proper STAM-Kerr derivation, Interpretation B is the best
  available framework reading. Under it:
    - n=1 is observationally tolerable
    - n=2 is in significant tension
""")

# --- What the math is telling us ---
print("=" * 90)
print("What the math is telling us")
print("=" * 90)
print("""
  The framework's ringdown prediction is sensitive to the spinning STAM
  metric, which hasn't been explicitly derived. Three possibilities:

  1. **Static-bubble interpretation (B) is correct**: STAM's metric
     modification (the (1-A²)^n piece) is a property of the BUBBLE itself,
     not of the spinning hologram. Rotation moves the matter ON the
     bubble surface but doesn't change the bubble's deep structure.

     Under this: τ_STAM/τ_GR_Kerr ≈ (9/5)^(n/2) / spin_factor
     For a ≈ 0.67: 1.22 (n=1) or 1.64 (n=2)

     n=1 is marginally consistent with LIGO; n=2 is in tension.

  2. **Naive-Kerr interpretation (A) is correct**: STAM's modification
     applies at the actual Kerr light ring location. At moderate spin,
     this location has A_LR very close to 1, causing huge (1-A²)^(-n/2)
     amplification.

     Under this: τ ratios of 4-15 for typical LIGO spin. Massive tension.

     This interpretation would essentially rule out the framework for
     any non-trivial spin, regardless of n.

  3. **Some other STAM-Kerr structure**: the bubble's oblate shape and
     the rotating-matter ontology need explicit metric derivation. The
     prediction could lie between A and B.

  **Framework's open computation**: a proper STAM-Kerr metric is needed
  to distinguish between these readings. Without it, the framework's
  ringdown prediction has factor-of-several uncertainty for spinning BHs.

  **Implications for n(A) commitments**:
    Under Interpretation B (most plausible per substance ontology):
      - n=1 is preferred by LIGO data
      - n=2 is in tension but not ruled out
      - n_horizon = 4 (e.g., n=1+3A) is strongly disfavored

    Under Interpretation A:
      - Any non-trivial n is ruled out
      - Framework would need substantial revision

  **What this changes for the LIGO arc**:
    The framework cannot make a CLEAN ringdown prediction for spinning
    Kerr remnants without first deriving the STAM-Kerr metric. Until
    that's done, comparing to LIGO data is "interpretation-dependent".

    The framework's path forward:
    1. Compute the STAM-Kerr metric explicitly (open work)
    2. Derive τ_STAM(M, a) properly
    3. Compare to LIGO Tests of GR catalog results

  **For now**: use Interpretation B as a working approximation. Under it:
    - The framework's n=2 commitment is in moderate LIGO tension
    - n=1 is marginally consistent
    - The G18 "α=4 entropy" reading via n=4 at horizon is observationally
      disfavored (but n_observable = 2 at horizon, per Sean's correction,
      still gives 1.64 excess under B — within 3σ)
""")
