"""
G49: STAM-native engine time — no double-dip on velocity
=========================================================

Correction from G48: the previous calculation applied the substance velocity-cap
with A_total = A_static + A_kin (self-consistent), which double-dips because
A_kin is itself v²/c² — so the velocity is used to compute the slowdown that
slows the velocity.

Fixed: apply cap with A_static (gravitational A) only.

  v_orbital² = (1 − A_static) × v_Newton²
  A_kin = v_orbital²/c²   (downstream consequence, not an input to the cap)
  Threshold: A_kin = A_0  →  x(1 − x) = A_0,  x = R_s/(2r)

Otherwise sequence unchanged from G48:
  - Donut max radius = v_sub × τ_chirp
  - Donut doesn't recede; explosion force reaches it at t_meet = r_donut/v_sub
  - Halfway condition: v_light = v_substance / 2 (matter slowdown factor 2)
  - Engine time = t_meet = τ_chirp

Don't aim. Report what falls out.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

M_total = 2.7 * M_sun
mu = M_total / 4.0
R_s = 2 * G * M_total / c**2

def A_static(r):
    return R_s / (2 * r)

def v_local(r):
    """Cap with A_static only — no double-dip."""
    x = A_static(r)
    v_Newton_sq = G * M_total / r
    return np.sqrt(max(1 - x, 0) * v_Newton_sq)

def A_kin(r):
    return v_local(r)**2 / c**2

def v_substance(A):
    return c * (1 - A)**2 * (1 + A)

v_sub_baseline = v_substance(A_0)

# Threshold: A_kin = A_0  →  x(1-x) = A_0
disc = 1 - 4 * A_0
x_thr = (1 - np.sqrt(disc)) / 2
r_threshold = R_s / (2 * x_thr)

tau_chirp = (5/64) * c**5 * r_threshold**4 / (G**3 * M_total**3)
r_donut_max = v_sub_baseline * tau_chirp

print("=" * 78)
print("G49: STAM-native engine time — no double-dip on velocity")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  v_substance baseline: {v_sub_baseline/c:.4f}c")

# Local speeds
print("\n" + "-" * 78)
print("Local orbital speeds (cap with A_static only)")
print("-" * 78)
print(f"\n  {'r/R_s':>8} {'A_static':>10} {'v_Newton/c':>14} {'v_local/c':>14} {'A_kin':>10}")
print(f"  {'-'*8} {'-'*10} {'-'*14} {'-'*14} {'-'*10}")
for u in [r_threshold/R_s, 12, 8, 5, 3, 2, 1.5, 1.1, 1.0]:
    r = u * R_s
    A_st = A_static(r)
    v_N = np.sqrt(G * M_total / r)
    v_L = v_local(r)
    A_k = A_kin(r)
    print(f"  {u:>8.3f} {A_st:>10.4f} {v_N/c:>14.4f} {v_L/c:>14.4f} {A_k:>10.4f}")

# Max v_local
r_search = np.linspace(1.001 * R_s, r_threshold, 10000)
v_search = np.array([v_local(r) for r in r_search])
v_max = v_search.max()
r_at_vmax = r_search[v_search.argmax()]
print(f"\n  Max v_local: {v_max/c:.4f}c at r = {r_at_vmax/R_s:.3f} R_s")
print(f"  v at threshold: {v_local(r_threshold)/c:.4f}c")

# Threshold
print("\n" + "-" * 78)
print("Threshold (A_kin = A_0 with no double-dip)")
print("-" * 78)
print(f"\n  Equation: x(1-x) = A_0  →  x² − x + A_0 = 0")
print(f"  x_threshold = {x_thr:.5f}")
print(f"  r_threshold/R_s = {r_threshold/R_s:.4f}")

# Chirp + donut
print("\n" + "-" * 78)
print("Chirp + donut")
print("-" * 78)
print(f"\n  τ_chirp (PM at threshold) = {tau_chirp:.4f} s")
print(f"  r_donut_max = v_sub × τ_chirp = {r_donut_max/c:.4f} light-seconds")

# Engine time
t_meet = r_donut_max / v_sub_baseline  # = τ_chirp
engine_time = t_meet
print("\n" + "-" * 78)
print("Engine time (donut doesn't recede, halfway condition: v_light = v_sub/2)")
print("-" * 78)
print(f"\n  t_meet = r_donut_max / v_substance = {t_meet:.4f} s")
print(f"  Engine time = t_meet = τ_chirp = {engine_time:.4f} s")
print(f"\n  Observed engine time (GW170817): 1.74 s")
print(f"  Gap: {(engine_time - 1.74)/1.74 * 100:+.2f}%")

# Comparison
print("\n" + "-" * 78)
print("Comparison: cap-with-A_total (G48) vs cap-with-A_static-only (G49)")
print("-" * 78)
print(f"""
  Framing                             r_thr/R_s   τ_chirp   gap vs 1.74
  G48 (cap with A_static + A_kin):    17.821      1.677 s    −3.62%
  G49 (cap with A_static only):       {r_threshold/R_s:.3f}      {tau_chirp:.4f} s   {(engine_time-1.74)/1.74*100:+.2f}%
  Naive Newton (no cap at all):       18.850      2.099 s   +20.7%
""")

# Mass scaling
print("-" * 78)
print("Mass-scaling prediction (G49, no double-dip)")
print("-" * 78)
print(f"\n  Slope: {engine_time/2.7:.4f} s/M_sun")
print(f"\n  Predicted engine times:")
print(f"  {'M_total (M_sun)':>16} {'τ_engine (s)':>14}")
for M in [2.5, 2.6, 2.7, 2.74, 2.8, 3.0, 4.0]:
    M_kg = M * M_sun
    R_s_M = 2 * G * M_kg / c**2
    r_thr_M = R_s_M / (2 * x_thr)
    tau_M = (5/64) * c**5 * r_thr_M**4 / (G**3 * M_kg**3)
    print(f"  {M:>16.2f} {tau_M:>14.4f}")

# What M_total would give exactly 1.74 s
# τ ∝ M, so M_for_1.74 = 2.7 × 1.74/engine
M_for_obs = 2.7 * 1.74 / engine_time
print(f"\n  M_total that lands at exactly 1.74 s: {M_for_obs:.3f} M_sun")

# Listening
print("\n" + "=" * 78)
print("What the math is telling us")
print("=" * 78)
print(f"""
  Removing the double-dip moves the threshold outward (18.33 R_s vs 17.82),
  which moves τ_chirp from 1.677 s to {tau_chirp:.3f} s — and shifts the prediction
  from -3.6% under observation to {(engine_time-1.74)/1.74*100:+.1f}% over.

  Observation sits between the two framings:
    G48 (double-dip):     1.677 s  (-3.6% under)
    G49 (no double-dip):  {engine_time:.3f} s  ({(engine_time-1.74)/1.74*100:+.1f}% over)
    Observed:             1.74 s

  The structural identity τ_engine = τ_chirp holds in both framings.
  What changes is τ_chirp itself, via the threshold definition.

  G49 needs M_total ≈ {M_for_obs:.2f} M_sun to land exactly at observed.
  GW170817's chirp-mass-derived M_total is 2.74±~0.04 M_sun. The +8%
  gap doesn't fit inside this mass uncertainty (which would shift by ~1.5%).

  The math is pointing at one of:
    (a) The no-double-dip framing overstates τ_chirp slightly
    (b) Some piece between Newton-cap-with-A_static and self-consistent
        feedback is the right structural reading
    (c) Observation has another contribution (conventional astrophysics)
        beyond the STAM-geometric piece
""")
