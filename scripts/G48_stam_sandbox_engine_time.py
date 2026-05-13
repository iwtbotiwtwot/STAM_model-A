"""
G48: STAM-native engine time — Sean's full sequence
=====================================================

Sandboxed in STAM. Sean's mechanism, step by step:

  1. Inspiral: bodies orbit with self-consistent A_static + A_kin.
     Compute v_local at each separation r along inspiral.

  2. Donut grows during inspiral.
     r_donut_max = v_substance(A_0) × τ_chirp at orbit-stop.

  3. At merger (t=0):
     • Light starts moving outward (slowed by A along its path)
     • Donut starts retracting inward at v_substance(A_0)
     • Explosion force expands outward at v_substance(A_0)
     • Net effect on donut: retracts slightly but is pushed back out by the
       explosion force. Take the donut's position as roughly stationary at
       r_donut_max during the relevant interval.

  4. Force of the explosion reaches donut at t_meet:
       t_meet = r_donut_max / v_substance(A_0)
     GW is born here, at r = r_donut_max.

  5. Light has been propagating since t=0 through the elevated-A event region.
     Light is slower than the explosion force because it traverses through
     A-elevated substance ("light must travel along A to leave the event").

  6. Engine time at observer = time light reaches r_donut_max minus
     time GW is born there.

The "halfway condition" Sean stated: at t_meet, light still needs the same
amount of time to reach the meeting that it has already traveled. This forces
v_light = v_substance / 2 (matter slowdown factor exactly 2). Under this
condition, engine time = t_meet = r_donut_max / v_substance.

Don't aim. Report what falls out.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# Constants
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

# GW170817
M_total = 2.7 * M_sun
mu = M_total / 4.0
R_s = 2 * G * M_total / c**2

# Self-consistent A at orbital position
def x_of(r):
    return R_s / (2 * r)

def A_static_orbit(r):
    return x_of(r)

def A_kin_orbit(r):
    x = x_of(r)
    return x * (1 - x) / (1 + x)

def A_total_orbit(r):
    return 2 * x_of(r) / (1 + x_of(r))

def v_local_orbit(r):
    """Body's local orbital speed with substance velocity-cap applied."""
    A = A_total_orbit(r)
    v_Newton_sq = G * M_total / r
    return np.sqrt((1 - A) * v_Newton_sq)

def v_substance(A):
    return c * (1 - A)**2 * (1 + A)

v_sub_baseline = v_substance(A_0)

# Threshold
disc = (1 - A_0)**2 - 4 * A_0
x_thr = (1 - A_0 - np.sqrt(disc)) / 2
r_threshold = R_s / (2 * x_thr)

# Chirp time (PM with STAM threshold; STAM ≡ GR at low A by commitment)
tau_chirp = (5/64) * c**5 * r_threshold**4 / (G**3 * M_total**3)

# Donut max radius
r_donut_max = v_sub_baseline * tau_chirp

print("=" * 78)
print("G48: STAM-native engine time — full sequence")
print("=" * 78)

# Step 1: Local orbital speeds
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  v_substance at baseline: {v_sub_baseline/c:.4f} c")

print("\n" + "-" * 78)
print("Step 1: Local orbital speeds (substance velocity-cap, A_static + A_kin)")
print("-" * 78)
print(f"\n  {'r/R_s':>8} {'A_static':>10} {'A_kin':>10} {'A_total':>10} {'v_Newton/c':>14} {'v_local/c':>14}")
print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*14} {'-'*14}")
for u in [17.82, 12, 8, 5, 3, 2, 1.5, 1.1, 1.0]:
    r = u * R_s
    A_st = A_static_orbit(r)
    A_kin = A_kin_orbit(r)
    A_t = A_total_orbit(r)
    v_N = np.sqrt(G*M_total/r)
    v_L = v_local_orbit(r)
    print(f"  {u:>8.3f} {A_st:>10.4f} {A_kin:>10.4f} {A_t:>10.4f} {v_N/c:>14.4f} {v_L/c:>14.4f}")

# Find max v_local
r_search = np.linspace(1.001*R_s, r_threshold, 10000)
v_search = np.array([v_local_orbit(r) for r in r_search])
v_max = v_search.max()
r_at_vmax = r_search[v_search.argmax()]
print(f"\n  Max v_local: {v_max/c:.4f} c at r = {r_at_vmax/R_s:.3f} R_s")
print(f"  v at threshold: {v_local_orbit(r_threshold)/c:.4f} c")

# Step 2: Donut radius
print("\n" + "-" * 78)
print("Step 2: Donut grows during inspiral, max at orbit-stop")
print("-" * 78)
print(f"\n  r_threshold (substance threshold) = {r_threshold/R_s:.4f} R_s")
print(f"  τ_chirp (Peters-Mathews from threshold) = {tau_chirp:.4f} s")
print(f"  r_donut_max = v_sub × τ_chirp = {r_donut_max/c:.4f} light-seconds")
print(f"             = {r_donut_max:.3e} m")
print(f"             = {r_donut_max/R_s:.0f} R_s")

# Step 3: Meeting (stationary-donut framing per Sean: retracts slightly, pushed back)
print("\n" + "-" * 78)
print("Step 3: Force of explosion reaches donut at r_donut_max")
print("-" * 78)
print(f"""
  Per Sean's framing: donut retracts slightly at orbit-stop but is pushed
  back out by the explosion force. Net: donut stays roughly at r_donut_max.
  Force of explosion propagates from R_s outward at v_substance(A_0).
""")
t_meet = r_donut_max / v_sub_baseline
r_meet = r_donut_max
print(f"  t_meet = r_donut_max / v_substance = {t_meet:.4f} s")
print(f"  GW born at r_meet = r_donut_max = {r_meet/c:.4f} ls")

# Step 4: Light traversal — halfway condition
print("\n" + "-" * 78)
print("Step 4: Light traverses elevated-A region, slowed by matter+A")
print("-" * 78)
print(f"""
  Sean's halfway condition at the meeting:
    At t_meet, light has traveled X seconds AND still needs X seconds.
    This forces v_light_avg × t_meet = r_donut_max / 2.
    Combined with t_meet = r_donut_max / v_substance:
      v_light_avg = v_substance / 2

  Under this condition, light takes 2 × t_meet to reach r_donut_max.
""")
v_light_required = v_sub_baseline / 2
print(f"  v_light required (halfway condition): v_substance/2 = {v_light_required/c:.4f} c")

t_light_to_r_donut = r_donut_max / v_light_required
print(f"  t_light to reach r_donut_max at v_light = v_sub/2: {t_light_to_r_donut:.4f} s")

# Step 5: Engine time
print("\n" + "-" * 78)
print("Step 5: Engine time at observer")
print("-" * 78)

engine_time = t_light_to_r_donut - t_meet
print(f"\n  Engine time = t_light_at_r_donut_max − t_meet")
print(f"              = {t_light_to_r_donut:.4f} − {t_meet:.4f}")
print(f"              = {engine_time:.4f} s")
print(f"\n  Observed GW170817 engine time: 1.74 s")
print(f"  STAM prediction: {engine_time:.4f} s")
print(f"  Gap: {(engine_time - 1.74)/1.74 * 100:+.2f}%")

# --- Structural identity ---
print("\n" + "-" * 78)
print("Structural identity")
print("-" * 78)
print(f"""
  Under the halfway condition (v_light = v_substance/2):
    engine_time = r_donut_max / v_substance = τ_chirp

  For GW170817:
    τ_chirp = {tau_chirp:.4f} s
    engine_time (this calculation) = {engine_time:.4f} s
    observed = 1.74 s
    gap = {(engine_time - 1.74)/1.74 * 100:+.2f}%

  The framework's prediction:
    τ_engine = r_donut_max / v_substance = τ_chirp

  Scales linearly with M_total:
    Slope = 1 / (M_sun × c² × 5 / (64 G³ M_total³) × ... ) — let me just compute numerically

  Slope d(engine)/d(M_total) = engine / M_total = {engine_time/2.7:.4f} s per M_sun
""")

# Mass scaling
print("-" * 78)
print("Mass-scaling prediction (falsification handle)")
print("-" * 78)
print(f"\n  τ_engine ∝ M_total via τ_chirp ∝ M_total")
print(f"  Slope: {engine_time/2.7:.4f} s/M_sun")
print(f"\n  Predicted engine times for various BNS systems:")
print(f"  {'M_total (M_sun)':>16} {'τ_engine (s)':>14}")
for M in [1.4*2, 1.4+1.3, 2.7, 1.5+1.5, 1.7+1.7, 2.0+2.0]:
    M_kg = M * M_sun
    R_s_M = 2 * G * M_kg / c**2
    # Recompute threshold and chirp for this mass
    r_thr_M = R_s_M / (2 * x_thr)
    tau_M = (5/64) * c**5 * r_thr_M**4 / (G**3 * M_kg**3)
    r_donut_M = v_sub_baseline * tau_M
    engine_M = r_donut_M / v_sub_baseline  # = tau_M under halfway
    print(f"  {M:>16.2f} {engine_M:>14.4f}")

print("\n" + "=" * 78)
print("What the math is telling us")
print("=" * 78)
print(f"""
  STAM-native sequence:
    • Bodies orbit at v_local determined by self-consistent A_static + A_kin
      (peak v_local = {v_max/c:.3f}c at r ≈ {r_at_vmax/R_s:.2f} R_s)
    • Donut builds to r_donut_max = {r_donut_max/c:.3f} ls = v_sub × τ_chirp
    • Force of explosion reaches donut at t_meet = {t_meet:.3f} s
    • Light, slowed by ~factor 2 (the halfway condition), reaches donut at
      2 × t_meet = {t_light_to_r_donut:.3f} s
    • Engine time = differential = t_meet = τ_chirp = {engine_time:.3f} s

  Match to observed 1.74 s: {(engine_time - 1.74)/1.74 * 100:+.2f}%

  Structural identity: τ_engine = r_donut_max / v_substance = τ_chirp.
  Scales linearly with M_total. Falsifiable with future BNS+EM events.

  The factor-2 matter slowdown (v_light = v_substance/2) is the open piece.
  If structurally derivable from kilonova physics in STAM, the prediction
  is closed. If it's conventional astrophysics input, STAM provides the
  geometry (r_donut_max) and conventional physics provides the magnitude.
""")
