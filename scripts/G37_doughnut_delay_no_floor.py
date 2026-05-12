"""
G37: GW170817 Doughnut — Photon Delay Through Interior High-A Region
======================================================================

Computes the delay for light to traverse the high-A interior of the
doughnut, FROM the merger event TO the doughnut edge.

Sean's question: how much delay does the elevated-A interior cause for
light to reach the doughnut edge?

Fix from G36: integrate A_system directly (no flooring at A_0).
The full A above zero is what produces extra path-stretching for light.
The A_0 baseline is shared by all paths and doesn't differentiate.

For the question "delay from event to doughnut edge," what we want is:
    Δt = (1/c) × ∫(A_static(r) + A_history(r)) dr from r_emit to r_d

This is the excess time light needs to traverse the elevated interior
compared to traversing baseline-only over the same distance.
"""

import numpy as np
from scipy.integrate import quad
import sys
sys.stdout.reconfigure(encoding='utf-8')

c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2

print("=" * 75)
print("G37: Doughnut interior delay (proper integration, no floor)")
print("=" * 75)
print(f"  M_total = 2.7 M_sun, R_s = {R_s/1000:.3f} km, A_0 = {A_0:.6f}")

# Inspiral parameterization (same as G36)
tau_band_entry = 100.0
r_orb_band_entry = 400e3

def r_orb(tau):
    if tau <= 0:
        return 2 * R_s
    r = r_orb_band_entry * (tau / tau_band_entry)**(1/4)
    return max(r, 2 * R_s)

def A_static(r):
    return R_s / r

def A_history(r):
    tau = r / c
    if tau <= 0 or tau > tau_band_entry:
        return 0.0
    return (G * M_total / r_orb(tau)) / c**2

def A_system(r):
    """System's contribution above cosmic baseline."""
    return A_static(r) + A_history(r)

# Define doughnut edge — where elevated A drops to baseline
# Where A_history = A_0: v² = c² × A_0 → r_orb = R_s × 6π
r_orb_critical = R_s / (2 * A_0)
tau_critical = tau_band_entry * (r_orb_critical / r_orb_band_entry)**4
r_doughnut_edge_strict = c * tau_critical  # ~600,000 km
r_doughnut_edge_full = c * tau_band_entry  # full inspiral extent ~30 million km

print(f"\nDoughnut edge definitions:")
print(f"  Strict (A_history drops to A_0):  r = {r_doughnut_edge_strict/1e6:.2f} × 10⁶ m")
print(f"                                       = {r_doughnut_edge_strict/c:.3f} light-seconds")
print(f"  Full propagation front:           r = {r_doughnut_edge_full/1e9:.2f} × 10⁹ m")
print(f"                                       = {r_doughnut_edge_full/c:.1f} light-seconds")

# Integrate A_system from r_emit to various edges
r_emit = 2 * R_s

def integrand(r):
    return A_system(r)

# Multiple integration ranges
print(f"\n" + "=" * 75)
print("Photon delay from merger to various radii (proper integration):")
print("=" * 75)
print(f"  Integrating ∫(A_static + A_history) dr from r_emit = {r_emit/1000:.2f} km")
print(f"  Delay = ∫A_system dr / c")
print()
print(f"  {'Integration to':<35} {'Integral (m)':<18} {'Delay':<20}")
print(f"  {'-'*73}")

ranges = [
    ("Strict doughnut edge (~600,000 km)", r_doughnut_edge_strict),
    ("1 million km", 1e9),
    ("10 million km", 1e10),
    ("Full propagation front (30 million km)", r_doughnut_edge_full),
]

for name, r_out in ranges:
    integral, _ = quad(integrand, r_emit, r_out, limit=500)
    delay = integral / c
    if delay > 0.1:
        delay_str = f"{delay:.4f} s"
    elif delay > 1e-3:
        delay_str = f"{delay*1000:.3f} ms"
    else:
        delay_str = f"{delay*1e6:.3f} μs"
    print(f"  {name:<35} {integral:<18.4e} {delay_str:<20}")

# Decompose: static-only vs history-only
print(f"\n" + "=" * 75)
print("Decomposition (proper, no floor):")
print("=" * 75)

for name, r_out in ranges:
    int_s, _ = quad(A_static, r_emit, r_out, limit=500)
    int_h, _ = quad(A_history, r_emit, r_out, limit=500)
    print(f"\n  Up to {name}:")
    print(f"    A_static contribution:  {int_s/c*1000:.3f} ms  (= R_s × ln(r/r_emit) / c)")
    print(f"    A_history contribution: {int_h/c*1000:.3f} ms  (= ∫A_velocity(τ=r/c) dr / c)")
    print(f"    Total:                  {(int_s+int_h)/c*1000:.3f} ms")

# What the math says
print(f"\n" + "=" * 75)
print("What the math reveals (forward, no target):")
print("=" * 75)
int_strict, _ = quad(A_system, r_emit, r_doughnut_edge_strict, limit=500)
int_full, _ = quad(A_system, r_emit, r_doughnut_edge_full, limit=500)
print(f"")
print(f"  Through strict doughnut edge (A above A_0): {int_strict/c*1000:.2f} ms")
print(f"  Through full propagation front:             {int_full/c:.3f} s")
print(f"")
print(f"  For reference, GW170817 observed delay: 1.74 s")
print(f"  (NOT a target — we are reading what the math says.)")
