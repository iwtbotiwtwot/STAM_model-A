"""
G40: STAM Engine Time Prediction — Bubble Retraction Time vs Observed Engine Times
====================================================================================

Sean's picture predicts:
  Engine time (orbit-stop → explosion) = bubble retraction time
                                       = τ_critical (Peters-Mathews scaling)

where τ_critical is when orbital v² drops to c² · A_0 during the inspiral.

This script computes the structural prediction for various binary masses
and compares to GW170817 (the only BNS event with observed engine time).
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Constants
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30

# STAM
A_0 = 1.0 / (12 * np.pi)

print("="*78)
print("G40: STAM engine time prediction — bubble retraction model")
print("="*78)
print(f"\nSTAM prediction: engine time = τ_critical")
print(f"  where τ_critical = time for orbital v² to drop to c²·A_0 during inspiral")
print(f"  Equivalently: orbital separation at threshold = R_s / (2·A_0) = R_s × 6π")

# Peters-Mathews chirp time at orbital separation r:
#   t_chirp(r) = (5/256) c^5 r^4 / (G^3 M_total^2 μ)
# For equal masses (μ = M_total/4):
#   t_chirp(r) = (5/64) c^5 r^4 / (G^3 M_total^3)
#
# At threshold r = R_s × 6π = 12πGM/c²:
#   t_chirp = (5/64) × (12π)^4 × R_s / (2c)
#           = (5/128) × (12π)^4 × R_s / c

def tau_critical(M_total):
    """STAM structural prediction: time from orbital threshold to merger."""
    R_s = 2 * G * M_total / c**2
    prefactor = (5 / 128) * (12 * np.pi)**4
    return prefactor * R_s / c

# Test for various binary mass configurations
print(f"\n{'Binary type':<22} {'M_total (M_sun)':<18} {'R_s (km)':<12} {'τ_critical (s)':<20}")
print(f"{'-'*72}")

cases = [
    ("BNS-light (1.4+1.4)", 2.8),
    ("BNS-GW170817 (this)", 2.7),
    ("BNS-heavy (2+2)",     4.0),
    ("NS-BH (1.4+5)",       6.4),
    ("BBH-light (8+8)",     16.0),
    ("BBH-GW150914",        65.0),
    ("BBH-mid (30+30)",     60.0),
    ("BBH-heavy (50+50)",   100.0),
    ("BBH-GW190521",       150.0),
]

results = []
for name, M_solar in cases:
    M = M_solar * M_sun
    R_s = 2 * G * M / c**2
    tau_c = tau_critical(M)
    results.append((name, M_solar, R_s, tau_c))
    print(f"{name:<22} {M_solar:<18.2f} {R_s/1000:<12.2f} {tau_c:<20.4f}")

# Compare to observation
print(f"\n" + "="*78)
print("Comparison to GW170817")
print("="*78)
tau_observed = 1.74
tau_predicted = tau_critical(2.7 * M_sun)
ratio = tau_predicted / tau_observed
discrepancy = (tau_predicted - tau_observed) / tau_observed * 100

print(f"\n  Observed engine time:  {tau_observed} s")
print(f"  STAM predicted:        {tau_predicted:.4f} s")
print(f"  Ratio:                 {ratio:.3f}")
print(f"  Discrepancy:           {discrepancy:+.1f}%")

# Mass scaling
print(f"\n" + "="*78)
print("Mass scaling of the prediction (key falsifiability handle)")
print("="*78)
print(f"\n  τ_critical scales as: M_total^1 × (12π)^4 × 5/(128·c)")
print(f"                      ∝ M_total  (linear in total mass)")
print(f"\n  So engine time should scale LINEARLY with total binary mass.")
print(f"")

# Show the linear scaling
M_ref = 2.7
tau_ref = tau_critical(M_ref * M_sun)
print(f"  Predicted engine times (scaling from GW170817 anchor):")
print(f"  {'M_total':<18} {'τ_predicted (s)':<20} {'Scaling check':<20}")
print(f"  {'-'*60}")
for name, M_solar, _, tau_c in results:
    scaled = tau_ref * (M_solar / M_ref)
    print(f"  {M_solar:<18.2f} {tau_c:<20.4f} {scaled:.4f} (= linear × {M_solar/M_ref:.2f})")

# Important caveat for BBH
print(f"\n" + "="*78)
print("Caveat for BBH")
print("="*78)
print("""
BBH mergers do NOT produce gamma-ray bursts (no light emission typically).
There's no "engine time" to compare against for BBH.

The prediction is testable only for BNS events with EM counterparts.
Currently: GW170817 is the only such event. We need more BNS-EM events
to test the linear mass scaling.

If/when LIGO/Virgo/KAGRA detect more BNS+EM events, the engine times
should fall on a line: τ_observed = τ_critical × (M_total / 2.7 M_sun) × 1.74/2.05

Slope: 1.74 / 2.7 = 0.644 seconds per M_sun
""")

# Optional structural test: if engine time is NOT mass-scaling, picture falsified
print("="*78)
print("Falsification handle")
print("="*78)
print(f"""
STAM picture predicts engine time ∝ M_total.
Standard astrophysics (HMNS collapse + jet launch) may have different scaling.

If future BNS events show:
  - Engine time scales linearly with M_total       → STAM picture supported
  - Engine time mass-independent or different scaling → STAM picture refuted (revert to coincidence)
  - Engine time wildly different for similar masses → both mechanisms have noise

Current data: 1 point (GW170817). Cannot distinguish coincidence from prediction.
""")
