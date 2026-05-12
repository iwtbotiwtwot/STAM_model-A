"""
G36: GW170817 Doughnut — Propagating Inspiral Model
=====================================================

Updated forward calculation using Sean's three commitments:
  1. Doughnut propagation is tied to orbital phase
     (A pulse emitted at time τ before merger is now at r = c·τ)
  2. A within the doughnut is to be discovered from the math
  3. Two contributions combine in the doughnut:
       - Mass: static R_s/r (current merger remnant)
       - Velocity: v_orb²/c² at the orbital phase that emitted from
         radius r (when system was at orbital separation matching τ = r/c)

STAM-only:
  - A_0 = 1/(12π) (V_4 derived)
  - A_static(r) = R_s/r (STAM weak-field A)
  - R_s = 2GM/c² (algebraic, STAM's mass-to-A bridge)
  - Keplerian + Peters-Mathews inspiral scaling r_orb ∝ τ^(1/4) used for
    inspiral evolution (STAM strong-field metric is GR-equivalent at the
    inspiral A values involved, so this scaling is STAM-consistent)

NO TARGET. Just running forward to see what falls out.
"""

import numpy as np
from scipy.integrate import quad
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Constants
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

# GW170817
M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2

print("=" * 75)
print("G36: GW170817 propagating doughnut (orbital-phase, mass + velocity)")
print("=" * 75)
print(f"  M_total = 2.7 M_sun")
print(f"  R_s = {R_s/1000:.3f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")

# --- Inspiral parameterization ---
# Calibration: at LIGO band entry (~100 s before merger), orbital separation
# corresponds to GW frequency ~25 Hz → orbital frequency ~12.5 Hz → period
# ~80 ms → Keplerian r_orb ≈ 400 km for M = 2.7 M_sun.
tau_band_entry = 100.0  # seconds before merger (LIGO band entry)
r_orb_band_entry = 400e3  # 400 km

def r_orb(tau):
    """Orbital separation at time τ before merger, Peters-Mathews scaling."""
    if tau <= 0:
        return 2 * R_s
    r = r_orb_band_entry * (tau / tau_band_entry)**(1/4)
    return max(r, 2 * R_s)

def v_orb_sq(tau):
    """Orbital velocity squared (Keplerian) at time τ before merger."""
    return G * M_total / r_orb(tau)

def A_velocity(tau):
    """A from orbital velocity at time τ before merger. Pure v²/c²."""
    return v_orb_sq(tau) / c**2

# Show inspiral evolution
print(f"\nInspiral evolution (Keplerian + Peters-Mathews):")
print(f"  {'τ (s before merger)':<25} {'r_orb (km)':<15} {'v_orb/c':<12} {'A_velocity':<15}")
print(f"  {'-'*70}")
for tau in [100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01, 0.001]:
    r = r_orb(tau)
    v = np.sqrt(v_orb_sq(tau))
    A_v = A_velocity(tau)
    above = "ABOVE A_0" if A_v > A_0 else "below A_0"
    print(f"  {tau:<25.4f} {r/1000:<15.2f} {v/c:<12.4f} {A_v:<15.5f} {above}")

# Find when A_velocity drops to A_0 (= structural threshold for doughnut)
# v_orb²/c² = A_0  →  GM/(r_orb c²) = A_0  →  r_orb = R_s/(2 A_0) = R_s × 6π
r_orb_critical = R_s / (2 * A_0)  # = R_s × 6π
# Solve for τ: r_orb_band_entry × (τ/τ_band)^(1/4) = r_orb_critical
tau_critical = tau_band_entry * (r_orb_critical / r_orb_band_entry)**4
print(f"\nThreshold (where v²/c² = A_0):")
print(f"  r_orb at threshold: R_s × 6π = {r_orb_critical/1000:.2f} km")
print(f"  τ at threshold:     {tau_critical:.3f} s before merger")
print(f"  Doughnut elevated-A extent at merger: r = c × τ = {c*tau_critical/1e6:.2f} × 10⁶ m")
print(f"                                      = {c*tau_critical/c:.3f} light-seconds")

# --- Doughnut A profile at merger time ---
# A pulse emitted at time τ during inspiral propagates outward at c, so at
# merger time, that pulse is at r = c × τ.
# Combined A at radius r at merger:
#   A_total(r) = R_s/r (static mass, current remnant)
#              + A_velocity(τ = r/c) (history contribution from inspiral)

def A_static(r):
    return R_s / r

def A_history(r):
    """A from inspiral velocity history: what was emitted at the system at
    time τ = r/c before merger, propagated outward at c."""
    tau = r / c
    if tau > tau_band_entry:
        return 0.0  # before LIGO band, assume no contribution
    return A_velocity(tau)

def A_total(r):
    return A_static(r) + A_history(r)

# Show the profile
print(f"\nA profile at merger time:")
print(f"  {'r':<22} {'r (km)':<15} {'A_static':<15} {'A_history':<15} {'A_total':<15}")
print(f"  {'-'*82}")
for r_m in [2*R_s, 10*R_s, 100*R_s, 1e6, 1e7, 1e8, 5e8, 1e9, 1e10, 3e10]:
    a_s = A_static(r_m)
    a_h = A_history(r_m)
    a_t = A_total(r_m)
    print(f"  {r_m:<22.3e} m {r_m/1000:<15.3e} {a_s:<15.5e} {a_h:<15.5e} {a_t:<15.5e}")

# --- Photon delay through the combined doughnut ---
# Integrate (A_total(r) - A_0) over the photon path from r_emit outward.
# Floor at A_0 (don't subtract baseline below itself; photons through baseline
# A_0 don't accumulate delay relative to baseline).

r_emit = 2 * R_s  # photon emitted just outside merger
r_outer_search = c * tau_band_entry  # = 3e10 m, where band-entry pulse is now

def integrand_above_baseline(r):
    """(A_total - A_0) when A_total > A_0, else 0."""
    A = A_total(r)
    return max(A - A_0, 0.0)

integral, err = quad(integrand_above_baseline, r_emit, r_outer_search, limit=500)
delay_total = integral / c
print(f"\n" + "=" * 75)
print("Photon delay calculation:")
print(f"  Integration: r_emit = {r_emit/1000:.2f} km → r_outer = {r_outer_search/1e9:.2f} × 10⁹ m")
print(f"  Integrand: max(A_static(r) + A_history(r) - A_0, 0)")
print(f"")
print(f"  ∫(A - A_0)|_above dr = {integral:.3e} m")
print(f"  Total photon delay  = {delay_total:.6f} s")
print(f"                       = {delay_total*1000:.3f} ms")

# Decompose into static + history contributions for understanding
def integrand_static_only(r):
    return max(A_static(r) - A_0, 0.0)

def integrand_history_only(r):
    return max(A_history(r) - A_0, 0.0)

int_static, _ = quad(integrand_static_only, r_emit, r_outer_search, limit=500)
int_history, _ = quad(integrand_history_only, r_emit, r_outer_search, limit=500)
print(f"\n  Decomposition (each above-baseline only):")
print(f"    Static-only:   ∫(R_s/r - A_0)+ dr        = {int_static:.3e} m → {int_static/c*1000:.3f} ms")
print(f"    History-only:  ∫(A_velocity - A_0)+ dr   = {int_history:.3e} m → {int_history/c*1000:.3f} ms")
print(f"    Note: combined ≠ sum because both add together where they overlap")

# --- Compare to observed and to G35 static-only result ---
print(f"\n" + "=" * 75)
print("Context (no targets, just for orientation):")
print(f"  G35 static doughnut alone: ~70 μs")
print(f"  G36 propagating doughnut (this run): {delay_total*1000:.3f} ms")
print(f"  GW170817 observed delay:   1.74 s")

# --- What the math says ---
print(f"\n" + "=" * 75)
print("What the math reveals:")
print("=" * 75)
print(f"")
print(f"  Doughnut elevated-A extent: out to r ≈ {c*tau_critical/1e6:.1f} × 10⁶ m")
print(f"    (set by orbital v² dropping to c² × A_0 at τ ≈ {tau_critical:.2f} s)")
print(f"")
print(f"  Photon delay through this region: {delay_total*1000:.3f} ms")
print(f"")
print(f"  Structural scaling:")
print(f"    The doughnut extent depends on when v_orb² drops to c²·A_0")
print(f"    This is r_orb_critical = R_s × 6π, then τ ~ r_orb_critical^4 (Peters-Mathews)")
print(f"    Delay scales as (R_s) × structural_factor (from integration)")
