"""
G47: Donut max radius, meeting point, A profile, light transit time
====================================================================

Sandboxed in STAM. Computes (per Sean's 2026-05-12 mechanism):

  1. Max orbital speed during inspiral (self-consistent: substance velocity-cap
     with A_static + A_kin both contributing to A_local at the orbit)

  2. Max donut radius at orbit-stop (substance disturbance extent)

  3. Meeting point of receding donut and expanding explosion (= GW birth point)

  4. A profile from r = R_s (merger surface) out to r_meet, combining
     A_static (from merged mass) and A_donut_loaded (carried out by the
     inspiral substance disturbance, frozen-in from orbital A_kin at launch time)

  5. Light transit time through this profile = signal delay

GW is generated at the meeting (already at r_meet outside the dense event
region). Light starts at r ≈ R_s (merger surface) and has to traverse the
elevated-A interior to reach the same r_meet. The light's slowed transit
through this region is the engine time.

Don't aim at 1.74 s. Report what falls out from the structural commitments.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# --- Constants ---
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

# --- GW170817 ---
M_total = 2.7 * M_sun
mu = M_total / 4.0
R_s = 2 * G * M_total / c**2

# --- Substance kinematics ---
def A_static_binary(r):
    """A from binary's mass (treating merged remnant as point at r=0 with mass M)."""
    return R_s / r if r > R_s else 1.0  # cap at 1 inside R_s

def A_orbit_total(r_orbit):
    """A at orbital position during inspiral: A_static + A_kin self-consistent."""
    x = R_s / (2 * r_orbit)
    return 2 * x / (1 + x)  # = x + x(1-x)/(1+x)

def v_orbital_local(r_orbit):
    """Local orbital speed with substance velocity-cap applied (relative v)."""
    A = A_orbit_total(r_orbit)
    v_Newton_sq = G * M_total / r_orbit
    return np.sqrt((1 - A) * v_Newton_sq)

def v_substance(A):
    return c * (1 - A)**2 * (1 + A)

v_sub_baseline = v_substance(A_0)

# --- Threshold ---
disc = (1 - A_0)**2 - 4 * A_0
x_threshold = (1 - A_0 - np.sqrt(disc)) / 2
r_threshold = R_s / (2 * x_threshold)

print("=" * 78)
print("G47: Donut radius, meeting, A profile, light transit")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  r_threshold/R_s = {r_threshold/R_s:.4f}")
print(f"  v_substance at A_0 baseline: {v_sub_baseline/c:.4f} c")

# --- Step 1: Max orbital speed ---
print("\n" + "-" * 78)
print("Step 1: Max orbital speed during inspiral (STAM-self-consistent)")
print("-" * 78)

r_grid = np.linspace(1.001*R_s, r_threshold, 10000)
v_local_grid = np.array([v_orbital_local(r) for r in r_grid])
v_max = v_local_grid.max()
r_at_vmax = r_grid[v_local_grid.argmax()]
print(f"\n  Max v_local: {v_max/c:.4f} c at r = {r_at_vmax/R_s:.3f} R_s")
print(f"  (substance velocity-cap factor at peak: {v_max/np.sqrt(G*M_total/r_at_vmax):.4f})")
print(f"  v at threshold: {v_orbital_local(r_threshold)/c:.4f} c")
print(f"  v at r = 1.1 R_s: {v_orbital_local(1.1*R_s)/c:.4f} c")

# --- Step 2: Max donut radius ---
print("\n" + "-" * 78)
print("Step 2: Max donut radius at orbit-stop")
print("-" * 78)

# Chirp time using STAM threshold + GR-equivalent (low A) emission rate
# τ_PM(r_thr) = (5/64) c⁵ r_thr⁴ / (G³ M³) for equal-mass binary
tau_chirp = (5/64) * c**5 * r_threshold**4 / (G**3 * M_total**3)
print(f"\n  Chirp time τ_chirp = {tau_chirp:.4f} s (STAM threshold, no cap on chirp rate)")

# Donut grows at v_substance (wave-propagation reading)
r_donut_max = v_sub_baseline * tau_chirp
print(f"  r_donut_max (wave propagation at v_substance):")
print(f"    = v_sub × τ_chirp = {v_sub_baseline/c:.4f}c × {tau_chirp:.4f}s")
print(f"    = {r_donut_max/c:.4f} light-seconds")
print(f"    = {r_donut_max:.3e} m")
print(f"    = {r_donut_max/R_s:.0f} R_s")

# --- Step 3: Meeting point ---
print("\n" + "-" * 78)
print("Step 3: Meeting point of receding donut and expanding explosion")
print("-" * 78)

# Both at v_substance(A_0) for first pass
v_recede = v_sub_baseline
v_explosion = v_sub_baseline
t_meet = r_donut_max / (v_recede + v_explosion)
r_meet = v_explosion * t_meet

print(f"\n  Donut recedes inward at v_substance(A_0) = {v_recede/c:.4f}c")
print(f"  Explosion expands at v_substance(A_0) = {v_explosion/c:.4f}c")
print(f"  Meeting time t_meet = {t_meet:.4f} s")
print(f"  Meeting radius r_meet = r_donut_max/2 = {r_meet/c:.4f} light-seconds")
print(f"                                       = {r_meet/R_s:.1f} R_s")

# --- Step 4: A profile from R_s to r_meet ---
print("\n" + "-" * 78)
print("Step 4: A profile from R_s (merger surface) to r_meet")
print("-" * 78)

# Substance at radius r in the donut at orbit-stop was launched at time
#   t_launch = τ_chirp - r/v_sub  (substance front propagated from binary to r in this time)
# At t_launch, orbital position was r_orbit(t_launch) (chirp law)
# Chirp law: r_orbit(t) = r_thr × (1 - t/τ_chirp)^(1/4)
# A_donut_loaded(r) = A_kin at orbit when launched = self-consistent value at r_orbit(t_launch)

def A_donut_loaded(r):
    """A carried into the donut by substance launched during inspiral."""
    if r <= 0 or r >= r_donut_max:
        return A_0
    t_launch = tau_chirp - r / v_sub_baseline
    if t_launch < 0:
        return A_0  # substance hasn't been launched yet (shouldn't happen)
    frac_remaining = max(1 - t_launch/tau_chirp, R_s/r_threshold)
    r_orbit_at_launch = r_threshold * frac_remaining**(1/4)
    x_orbit = R_s / (2 * r_orbit_at_launch)
    A_kin = x_orbit * (1 - x_orbit) / (1 + x_orbit)
    return max(A_kin, A_0)

def A_total_at_r(r):
    """Total A at radius r in the donut at orbit-stop.
    Uses additive composition for first-pass: A_static + A_donut_loaded.
    Caps at 1 inside R_s."""
    A_st = A_static_binary(r)
    A_dl = A_donut_loaded(r)
    return min(A_st + A_dl, 0.999)

# Sample the profile
print(f"\n  A profile (additive composition: A_static + A_donut_loaded):")
print(f"  {'r/R_s':>10} {'r (ls)':>10} {'A_static':>10} {'A_donut':>10} {'A_total':>10} {'v_light/c':>10}")
print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
sample_r_over_Rs = [1.0, 1.5, 2.0, 5.0, 10.0, 50.0, 500.0, 5000.0, 50000.0]
for u in sample_r_over_Rs:
    r = u * R_s
    if r > r_meet:
        break
    A_st = A_static_binary(r)
    A_dl = A_donut_loaded(r)
    A_t = A_total_at_r(r)
    v_l = v_substance(A_t)
    print(f"  {u:>10.1f} {r/c:>10.6f} {A_st:>10.5f} {A_dl:>10.5f} {A_t:>10.5f} {v_l/c:>10.5f}")

# --- Step 5: Light transit time from R_s_emit to r_meet ---
print("\n" + "-" * 78)
print("Step 5: Light transit time")
print("-" * 78)

# Light has to traverse from some emission point r_emit out to r_meet.
# r_emit can't be exactly R_s (A → 1 there, integral diverges).
# Try several emission points:

def light_transit(r_start, r_end, n_pts=50000):
    r_grid = np.geomspace(r_start, r_end, n_pts)
    A_grid = np.array([A_total_at_r(r) for r in r_grid])
    v_grid = v_substance(A_grid)
    # ∫ dr / v
    # Use trapezoidal on 1/v
    inv_v = 1.0 / v_grid
    return np.trapezoid(inv_v, r_grid)

print(f"\n  Light transit from various r_emit out to r_meet = {r_meet/c:.4f} ls:")
print(f"  {'r_emit/R_s':>12} {'A at r_emit':>14} {'transit time (s)':>18} {'gap vs 1.74':>14}")
print(f"  {'-'*12} {'-'*14} {'-'*18} {'-'*14}")
for r_emit_factor in [1.001, 1.01, 1.05, 1.1, 1.5, 2.0, 5.0, 10.0, 50.0, 500.0]:
    r_emit = r_emit_factor * R_s
    if r_emit > r_meet:
        continue
    A_emit = A_total_at_r(r_emit)
    t_transit = light_transit(r_emit, r_meet)
    gap = (t_transit - 1.74)/1.74 * 100
    print(f"  {r_emit_factor:>12.3f} {A_emit:>14.5f} {t_transit:>18.4f} {gap:>+13.2f}%")

# --- Listening ---
print("\n" + "=" * 78)
print("What the math is telling us")
print("=" * 78)
print(f"""
  Max donut radius: {r_donut_max/c:.3f} light-seconds ({r_donut_max/R_s:.0f} R_s)
  Meeting point:    {r_meet/c:.3f} light-seconds
  Max orbital v:    {v_max/c:.3f} c at r ≈ {r_at_vmax/R_s:.2f} R_s

  Light's transit time depends strongly on where light is emitted (r_emit):
  the A_static = R_s/r piece diverges at r → R_s, so transit time depends on
  how close to R_s the emission actually happens.

  Structural piece needed for a derived prediction:
    • Where does light first emerge from the merger event?
    • Equivalently: what is the merger remnant's "photospheric" r where
      A drops to a point that light can escape?

  Without this commitment, the signal delay is a function of r_emit, not
  a derived number. The framework needs a structural argument for r_emit
  before this mechanism becomes a falsifiable prediction.

  The A_donut_loaded piece (substance frozen-in from orbital A_kin during
  inspiral) contributes only weakly: peak ~0.17 near r=0, falling to A_0 at
  the outer donut. Most of the light transit time comes from the A_static
  region just outside R_s where A is close to 1.
""")
