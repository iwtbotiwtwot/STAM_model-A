"""
G41: STAM-Corrected Inspiral Time (substance-velocity correction)
==================================================================

The substance ontology says: motion through elevated A is REAL slowdown,
not just observable. Applied to the inspiral, the binary's effective
orbital motion through elevated A is slower than v_Newton predicts.

This means GW emission (which depends on orbital frequency squared) is
reduced. But the orbit also takes longer per "Newtonian cycle" — the
binary moves through orbital phase more slowly.

We apply three candidate STAM corrections to the standard Peters-Mathews
inspiral integration and see what falls out for GW170817.

Each candidate is a candidate framework commitment for how substance
elevates A affects local coordinate motion.
"""

import numpy as np
from scipy.integrate import quad
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Constants
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30

# STAM
A_0 = 1.0 / (12 * np.pi)

# GW170817
M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2  # ~ 8 km

print("="*78)
print("G41: STAM-corrected inspiral merger time")
print("="*78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.3f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")

# --- Local A at the binary's position ---
# For equal-mass binary at separation r:
#   - Each NS sees A from the other: A_from_other = R_s_total/(2r)
#   - Each NS has kinematic A from its v²/c²: A_kin = R_s_total/(8r)
#   - Total: A_local = 5 R_s / (8r)
def A_local(r):
    return 5 * R_s / (8 * r)

# --- STAM correction candidates (modifies effective orbital motion) ---
def f_corr_propertime(A):
    """v_eff = v_Newton × √(1-A)  (proper-time factor)"""
    return np.sqrt(1 - A)

def f_corr_linear(A):
    """v_eff = v_Newton × (1-A)  (linear scaling)"""
    return 1 - A

def f_corr_metric(A):
    """v_eff = v_Newton × (1-A)² × (1+A)  (radial light-speed factor)"""
    return (1 - A)**2 * (1 + A)

# --- Peters-Mathews energy loss rate, with STAM correction ---
# GW emission rate ~ ω⁶ where ω is angular orbital frequency.
# If v_eff = v_Newton × f(A), then ω_eff = ω_Newton × f(A)
# So dE/dt_eff = dE/dt_Newton × f(A)⁶
#
# Peters-Mathews: dr/dt = -(64/5) G³ M³ / (c⁵ r³)  (equal mass)
# With STAM correction (less GW emission):
#   dr/dt_STAM = dr/dt_Newton × f(A_local(r))⁶
#
# Merger time = ∫_0^{r_initial} dr / |dr/dt|
#             = ∫_0^{r_initial} dr × (5 c⁵ r³)/(64 G³ M³ × f(A(r))⁶)

def merger_time(r_initial, correction_func, n_points=10000):
    """Compute STAM-corrected merger time from initial separation r_initial."""
    r_min = 2 * R_s  # geometric merger threshold

    # Integration of |dr/dt|^-1
    rs = np.linspace(r_min, r_initial, n_points)
    integrand_vals = np.zeros_like(rs)
    for i, r in enumerate(rs):
        A = A_local(r)
        if A >= 1:
            integrand_vals[i] = 0  # at saturation
            continue
        f = correction_func(A)
        # |dr/dt|^-1 = 5c^5 r^3 / (16 G^3 M^3 * f^6)  (for equal-mass binary)
        integrand_vals[i] = (5 * c**5 * r**3) / (16 * G**3 * M_total**3 * f**6)

    return np.trapezoid(integrand_vals, rs)

# --- Find threshold orbital separation (A_local = A_0) ---
# A_local(r) = 5 R_s / (8r) = A_0 → r = 5 R_s / (8 A_0)
r_threshold = 5 * R_s / (8 * A_0)
print(f"\nDoughnut threshold (A_local = A_0):")
print(f"  r_threshold = 5R_s/(8A_0) = {r_threshold/1000:.2f} km = {r_threshold/R_s:.2f} R_s")
print(f"  At threshold: A_local = {A_local(r_threshold):.5f}")

# --- Baseline (no STAM correction) ---
print("\n" + "="*78)
print("Baseline: Peters-Mathews leading-order (no STAM correction)")
print("="*78)
def f_none(A):
    return 1.0

t_baseline = merger_time(r_threshold, f_none)
print(f"  τ_merger from threshold (Newtonian PM) = {t_baseline:.4f} s")

# --- STAM-corrected times ---
print("\n" + "="*78)
print("STAM-corrected merger times")
print("="*78)

corrections = [
    ("v_eff = v_N × √(1-A)         (proper-time factor)", f_corr_propertime),
    ("v_eff = v_N × (1-A)           (linear)",            f_corr_linear),
    ("v_eff = v_N × (1-A)²(1+A)     (light radial factor)", f_corr_metric),
]

print(f"\n  {'Correction':<55} {'τ (s)':<12} {'vs observed (1.74s)':<20}")
print(f"  {'-'*100}")
print(f"  {'BASELINE (no correction)':<55} {t_baseline:<12.4f} {(t_baseline-1.74)/1.74*100:+.1f}%")
for name, fc in corrections:
    t = merger_time(r_threshold, fc)
    ratio = (t - 1.74) / 1.74 * 100
    print(f"  {name:<55} {t:<12.4f} {ratio:+.1f}%")

# --- Decomposition: how does the correction work? ---
print("\n" + "="*78)
print("Why this happens — the integrand profile")
print("="*78)
print(f"\n  At various orbital separations, what's A_local and correction factor?")
print(f"\n  {'r (km)':<12} {'r/R_s':<10} {'A_local':<10} {'(1-A)^6':<12} {'(1-A)^(6/2)':<14} {'(1-A)²(1+A)^6':<18}")
print(f"  {'-'*80}")
for r_factor in [r_threshold/R_s, 50, 20, 10, 5, 3, 2.1]:
    r = r_factor * R_s
    A = A_local(r)
    if A >= 1:
        continue
    f_pt = f_corr_propertime(A)**6
    f_lin = f_corr_linear(A)**6
    f_met = f_corr_metric(A)**6
    print(f"  {r/1000:<12.2f} {r_factor:<10.2f} {A:<10.4f} {f_lin:<12.5f} {f_pt:<14.5f} {f_met:<18.5f}")

# --- Summary ---
print("\n" + "="*78)
print("Summary")
print("="*78)
print(f"""
GW170817 observed engine time: 1.74 s
Peters-Mathews leading-order (no STAM): {t_baseline:.3f} s

STAM-corrected predictions (reducing GW emission by f(A)⁶):
  Proper-time factor √(1-A):  LONGER inspiral (correction weakens GW emission)
  Linear (1-A):                LONGER inspiral
  Light radial (1-A)²(1+A):    LONGER inspiral

⚠ All three corrections move in the WRONG direction (longer, not shorter).

This is structurally important:
  • The substance ontology's "motion through dense substance is slower"
    predicts REDUCED GW emission per unit orbital phase
  • This makes the inspiral take LONGER than Peters-Mathews leading-order
  • But observed is SHORTER than Peters-Mathews leading-order

So either:
  1. The substance correction has a different functional form than tested
  2. The substance correction is offset by another effect (e.g., enhanced GW
     emission via substance coupling, or substance drag adding energy loss)
  3. The 20% gap is standard GR post-Newtonian corrections, not STAM
  4. The match at all (2.10 s vs 1.74 s) is coincidence
""")
