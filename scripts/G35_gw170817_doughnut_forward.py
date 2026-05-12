"""
G35: GW170817 Doughnut — Forward Calculation
=============================================

Forward STAM calculation: given the merger system and STAM's A field
mechanism (∫A ds = photon path delay), compute what delay falls out for
light traversing the residual A profile around the merger after the
inspiral has built it up.

NO TARGET. We are not aiming for 1.74s. We just run the math forward
and see what comes out.

Per Sean's "scripts as listening instruments" methodology — the result is
a breadcrumb, not a verdict.

STAM-only commitments used:
  - A_0 = 1/(12π) (derived in V_4)
  - A(r) = R_s/r for static gravitational field (STAM weak-field def)
  - R_s = 2GM/c² (algebraically same as GR Schwarzschild radius)
  - ∫(A - A_0) ds along photon path gives delay above baseline

Inputs sandboxed FROM GR/LCDM where possible:
  - M = 2.7 M_sun: from LIGO GW170817 chirp-mass extraction. Uses GR templates
    in extraction, but STAM strong-field metric is GR-equivalent at the
    inspiral A values involved (A < 0.3), so this value is STAM-consistent.
  - We do NOT use LCDM distances. The doughnut calculation is local
    around the merger.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Physical constants (framework-agnostic)
c = 299792458.0   # m/s
G = 6.67430e-11   # m^3 kg^-1 s^-2
M_sun = 1.989e30  # kg

# STAM-only
A_0 = 1.0 / (12.0 * np.pi)

print("=" * 70)
print("STAM commitments:")
print(f"  A_0 = 1/(12π) = {A_0:.8f}")
print(f"  A(r) = R_s/r  (weak-field A definition)")
print(f"  R_s = 2GM/c²  (algebraic, STAM's mass-to-A bridge)")
print(f"  Δt_doughnut = (1/c) ∫(A - A_0) ds  (above baseline)")
print("=" * 70)

# ---------------------------------------------------------------------
# GW170817 observational input
# ---------------------------------------------------------------------
M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2

print(f"\nMerger system (STAM-derived from M = 2.7 M_sun):")
print(f"  R_s = {R_s:.2f} m = {R_s/1000:.3f} km")

# Where does static A drop to A_0?
r_doughnut_outer = R_s / A_0
print(f"\nNatural doughnut outer boundary (where A drops to A_0):")
print(f"  r* = R_s/A_0 = R_s × 12π = {r_doughnut_outer:.2f} m = {r_doughnut_outer/1000:.2f} km")
print(f"  This is the natural radial extent of the elevated-A region around")
print(f"  the merger, set entirely by R_s and A_0.")

# ---------------------------------------------------------------------
# Static doughnut: A(r) = R_s/r above baseline A_0
# Photon emitted at r_emit, traverses outward to r_doughnut_outer,
# then to infinity through baseline A_0.
#
# Delay above baseline within the doughnut region:
#   Δt = (1/c) ∫_{r_emit}^{r*} (R_s/r - A_0) dr
#      = (1/c) [R_s ln(r*/r_emit) - A_0 (r* - r_emit)]
# ---------------------------------------------------------------------

def doughnut_delay(r_emit):
    """Static-doughnut photon delay above baseline. Returns seconds."""
    if r_emit >= r_doughnut_outer:
        return 0.0
    integral = R_s * np.log(r_doughnut_outer / r_emit) - A_0 * (r_doughnut_outer - r_emit)
    return integral / c

print(f"\n" + "=" * 70)
print("Forward calculation: static doughnut delay vs. emission radius")
print("=" * 70)
print(f"  Photon emitted at r_emit, traverses outward through static A(r) = R_s/r")
print(f"  Delay above baseline is what falls out of the integral.\n")
print(f"  {'r_emit':<20} {'r_emit (km)':<15} {'∫(A-A_0)dr (m)':<20} {'Delay (s)':<15} {'Delay':<15}")
print(f"  {'-'*85}")

results = []
for r_emit_factor in [1.0, 1.5, 3.0, 10.0, 100.0, 1000.0]:
    r_emit = r_emit_factor * R_s
    dt = doughnut_delay(r_emit)
    integral_m = dt * c
    if dt > 1:
        dt_str = f"{dt:.3f} s"
    elif dt > 1e-3:
        dt_str = f"{dt*1000:.3f} ms"
    elif dt > 1e-6:
        dt_str = f"{dt*1e6:.3f} μs"
    else:
        dt_str = f"{dt*1e9:.3f} ns"
    print(f"  {r_emit_factor:.1f} R_s{'':<11} {r_emit/1000:<15.2f} {integral_m:<20.3e} {dt:<15.3e} {dt_str}")
    results.append((r_emit_factor, dt))

# ---------------------------------------------------------------------
# Include velocity-induced A buildup in the inspiral region
# At separation r during inspiral, orbital v² = GM/r = R_s c²/(2r)
# So A_velocity(r) = R_s/(2r) in the inspiral region.
#
# If the doughnut retains this velocity-induced buildup after merger
# (residual warping), then:
#   A_total(r) ≈ R_s/r + R_s/(2r) = 3R_s/(2r)  in inspiral region
#
# Doughnut outer radius for total A:
#   3R_s/(2r) = A_0  →  r* = 3R_s/(2 A_0) = R_s × 18π
# ---------------------------------------------------------------------

r_doughnut_outer_v = 1.5 * R_s / A_0  # 3R_s/(2A_0) = R_s × 18π
print(f"\n" + "=" * 70)
print("Forward calculation: with velocity-induced A in inspiral region")
print("=" * 70)
print(f"  A_total(r) = 3R_s/(2r)  (static + orbital velocity buildup)")
print(f"  Doughnut outer radius: {r_doughnut_outer_v/1000:.2f} km (= R_s × 18π)")

def doughnut_delay_with_velocity(r_emit):
    """Doughnut delay if A profile includes velocity buildup. Seconds."""
    if r_emit >= r_doughnut_outer_v:
        return 0.0
    integral = 1.5 * R_s * np.log(r_doughnut_outer_v / r_emit) - A_0 * (r_doughnut_outer_v - r_emit)
    return integral / c

print(f"\n  {'r_emit':<20} {'Delay (static)':<20} {'Delay (with v)':<20} {'Ratio':<10}")
print(f"  {'-'*70}")
for r_emit_factor in [1.0, 1.5, 3.0, 10.0, 100.0]:
    r_emit = r_emit_factor * R_s
    dt_static = doughnut_delay(r_emit)
    dt_v = doughnut_delay_with_velocity(r_emit)
    ratio = dt_v / dt_static if dt_static > 0 else float('inf')
    print(f"  {r_emit_factor:.1f} R_s{'':<11} {dt_static*1e6:<20.3f} μs {dt_v*1e6:<20.3f} μs {ratio:<10.3f}")

# ---------------------------------------------------------------------
# Reference scales for context (no target match)
# ---------------------------------------------------------------------
print(f"\n" + "=" * 70)
print("Reference scales for context (math output speaks for itself):")
print("=" * 70)
print(f"  STAM merger R_s: {R_s/1000:.3f} km")
print(f"  Static-doughnut outer radius (R_s × 12π): {r_doughnut_outer/1000:.2f} km")
print(f"  Velocity-doughnut outer radius (R_s × 18π): {r_doughnut_outer_v/1000:.2f} km")
print(f"")
print(f"  GW170817 observed GW-to-GRB delay: 1.74 s (for reference only,")
print(f"  not a target — listed last so the forward result speaks first).")

# ---------------------------------------------------------------------
# What the math says
# ---------------------------------------------------------------------
print(f"\n" + "=" * 70)
print("What the math reveals (forward calculation, no fit):")
print("=" * 70)
dt_inspiral_static = doughnut_delay(1.5 * R_s)
dt_inspiral_v = doughnut_delay_with_velocity(1.5 * R_s)
print(f"")
print(f"  For photon emitted near merger (r_emit = 1.5 R_s):")
print(f"    Static doughnut delay:        {dt_inspiral_static*1e6:.2f} μs")
print(f"    With orbital velocity buildup: {dt_inspiral_v*1e6:.2f} μs")
print(f"")
print(f"  Both are ~microseconds. The doughnut spatial extent is the size of")
print(f"  the merger system (hundreds of km, set by R_s/A_0 = R_s × 12π).")
print(f"")
print(f"  This is what comes out. We can read it from here.")
