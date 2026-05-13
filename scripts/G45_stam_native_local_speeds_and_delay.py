"""
G45: STAM-native — local orbital speeds, doughnut extent, signal delay
=======================================================================

Sandboxed entirely in STAM. No LIGO-frame translation, no GR import.
Single consistent frame (the metric's coord time t, which is STAM's bookkeeping
clock; we're not converting to any observer's proper time).

Three questions, in order:
  1. How fast are the two objects moving locally at each separation r?
  2. Where does this push the doughnut/bubble outer edge by orbit-stop?
  3. What is the signal delay (GW front leads light by how much)?

STAM commitments used (all framework-internal, no GR):
  - Substance velocity-cap: v_local² = (1 - A_local) × v_Newton²
  - Self-consistent local A at orbital position: A_local = 2x/(1+x), x = R_s/(2r)
  - Threshold (doughnut formation): A_kin = A_0 ⇒ r_threshold/R_s = 17.82
  - Substance wave speed at baseline: v_sub = c(1-A_0)²(1+A_0)
  - STAM ≡ GR emission rate at low A (framework's principle #1 commitment), so
    P_GW ∝ μ²r⁴ω⁶ with the STAM-corrected ω.

Don't aim at 1.74 s. Report what falls out.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# --- Constants ---
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

# --- GW170817 binary ---
M_total = 2.7 * M_sun
mu = M_total / 4.0
R_s = 2 * G * M_total / c**2

# --- STAM substance kinematics ---
def A_total(r):
    """Self-consistent A at the orbital position (static + kinematic, both substance-corrected)."""
    x = R_s / (2 * r)
    return 2 * x / (1 + x)  # = x + x(1-x)/(1+x)

def v_local(r):
    """Local orbital speed with substance velocity-cap applied. Each body's
    speed in the equal-mass case is v_rel/2; here we report v_rel."""
    A = A_total(r)
    v_Newton_sq = G * M_total / r       # Keplerian relative velocity squared
    return np.sqrt((1 - A) * v_Newton_sq)

def omega_coord(r):
    """STAM-corrected orbital angular frequency in metric coord time."""
    return v_local(r) / r

v_sub_baseline = c * (1 - A_0)**2 * (1 + A_0)

# --- Threshold ---
disc = (1 - A_0)**2 - 4 * A_0
x_threshold = (1 - A_0 - np.sqrt(disc)) / 2
r_threshold = R_s / (2 * x_threshold)

print("=" * 78)
print("G45: STAM-native — local speeds, doughnut extent, signal delay")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  r_threshold/R_s = {r_threshold/R_s:.4f}")
print(f"  v_substance at baseline A_0: {v_sub_baseline/c:.4f} c")

# --- Question 1: How fast are they moving locally? ---
print("\n" + "-" * 78)
print("Q1: Local orbital speeds along the inspiral (substance velocity-cap applied)")
print("-" * 78)

sample_r_over_Rs = [17.82, 12.0, 8.0, 5.0, 3.0, 2.0, 1.5, 1.1, 1.0]
print(f"\n  {'r/R_s':>8} {'A_total':>10} {'v_Newton/c':>14} {'v_local/c':>12} {'cap factor':>14}")
print(f"  {'-'*8} {'-'*10} {'-'*14} {'-'*12} {'-'*14}")
for u in sample_r_over_Rs:
    r = u * R_s
    A = A_total(r)
    v_N = np.sqrt(G * M_total / r)
    v_L = v_local(r)
    print(f"  {u:>8.3f} {A:>10.4f} {v_N/c:>14.4f} {v_L/c:>12.4f} {v_L/v_N:>14.4f}")

print(f"\n  Threshold v_local: {v_local(r_threshold)/c:.4f} c "
      f"({v_local(r_threshold)/np.sqrt(G*M_total/r_threshold):.3f} of Newton)")
print(f"  Just outside R_s (r=1.1 R_s) v_local: {v_local(1.1*R_s)/c:.4f} c "
      f"(velocity cap is reducing motion by {(1 - v_local(1.1*R_s)/np.sqrt(G*M_total/(1.1*R_s)))*100:.1f}%)")

# --- Question 2: Where does this push the doughnut to? ---
print("\n" + "-" * 78)
print("Q2: Doughnut/bubble outer edge at orbit-stop")
print("-" * 78)

# Energy loss rate with STAM-corrected ω:
#   P_coord ∝ μ²r⁴ω_coord⁶ = μ²r⁴((1-A)·GM/r³)³ = (1-A)³ × P_PM
# So dr/dt_STAM = (1-A)³ × dr/dt_PM_coord
#   τ_chirp_STAM = ∫_{R_s}^{r_thr} r³ (1-A(r))^(-3) dr / K
#                = τ_PM × ⟨(1-A)^(-3)⟩_{r³-weighted}

n_pts = 20000
r_grid = np.linspace(R_s, r_threshold, n_pts)
A_grid = np.array([A_total(r) for r in r_grid])
integrand_STAM = r_grid**3 / (1 - A_grid)**3
integrand_PM   = r_grid**3
I_STAM = np.trapezoid(integrand_STAM, r_grid)
I_PM   = np.trapezoid(integrand_PM, r_grid)
ratio  = I_STAM / I_PM

# τ_PM with STAM threshold
tau_PM_at_threshold = (5/8) * (r_threshold/R_s)**4 * (R_s/c)
tau_chirp_STAM = tau_PM_at_threshold * ratio

print(f"\n  Chirp duration in coord time:")
print(f"    τ_PM with STAM threshold (no velocity-cap on chirp): {tau_PM_at_threshold:.4f} s")
print(f"    Velocity-cap correction factor ⟨(1-A)^(-3)⟩:         {ratio:.4f}")
print(f"    τ_chirp with substance cap on orbital motion:        {tau_chirp_STAM:.4f} s")

# Bubble extent
d_bubble = v_sub_baseline * tau_chirp_STAM
print(f"\n  Bubble outer edge at orbit-stop:")
print(f"    d_bubble = v_substance × τ_chirp = {d_bubble/c:.4f} light-seconds")
print(f"             = {d_bubble:.3e} m")
print(f"             = {d_bubble/1000:.0f} km")

# For comparison, no-cap version
tau_no_cap = tau_PM_at_threshold
d_bubble_no_cap = v_sub_baseline * tau_no_cap
print(f"\n  (without velocity-cap on chirp rate: d_bubble = {d_bubble_no_cap/c:.4f} light-seconds)")

# --- Question 3: Signal delay ---
print("\n" + "-" * 78)
print("Q3: Signal delay (GW head start over light)")
print("-" * 78)

# Mechanism (Sean's head-start framing):
#   GW front = bubble outer edge at orbit-stop, propagating outward at v_substance
#   Light = electromagnetic emission from r=0 at orbit-stop
#   Both propagate at v_substance through cosmic baseline (F5c symmetric coupling)
#   Δt_observer = d_bubble / v_substance

delay = d_bubble / v_sub_baseline
print(f"\n  Δt_observer = d_bubble / v_substance = {delay:.4f} s")
print(f"  (Equivalent to τ_chirp, since bubble grew at v_sub during chirp.)")

# Naive Newton chirp for comparison (different threshold)
r_naive = 6 * np.pi * R_s
tau_naive = (5/8) * (r_naive/R_s)**4 * (R_s/c)
print(f"\n  Compare:")
print(f"    Naive Newton chirp (r_thr = 6π R_s, no STAM):  {tau_naive:.4f} s")
print(f"    STAM threshold + STAM-cap on chirp rate:       {tau_chirp_STAM:.4f} s")
print(f"    STAM threshold only (no cap on chirp rate):    {tau_PM_at_threshold:.4f} s")
print(f"    Observed engine time (GW170817):               1.74 s")

# --- Listening to the math ---
print("\n" + "=" * 78)
print("What the math is telling us")
print("=" * 78)

print(f"""
  Substance velocity-cap applied to orbital motion at every r:
    • Local speeds reduced from Newton by factor √(1-A) at each r —
      threshold v: {v_local(r_threshold)/c:.3f} c → close to merger: {v_local(1.1*R_s)/c:.3f} c
    • Chirp rate reduced by (1-A)³ at each r → τ_chirp grows by factor {ratio:.3f}
    • τ_STAM_chirp ≈ τ_naive_Newton: substance correction at threshold
      (smaller r_thr) and substance correction at chirp rate (slower at each r)
      cancel almost exactly. STAM's two corrections sum to no net change
      from naive Newton at the engine-time level.
    • Doughnut extends to {d_bubble/c:.2f} light-seconds; signal delay = {delay:.2f} s.

  Compared to observed 1.74 s: STAM under this reading overshoots by
  {(delay - 1.74)/1.74*100:+.1f}%.

  Three honest readings of this:
    (a) Substance velocity-cap (1-A) is too strong on orbital motion;
        a softer form (e.g. √(1-A)) would land closer to observation.
    (b) Substance velocity-cap doesn't apply to circular orbital motion
        the same way it does to translational motion. The framework needs
        a structural commitment on what "encountering substance" means
        for circular orbits.
    (c) The (1-A)³ → cancellation against threshold reduction is structural
        and STAM truly does predict ≈ naive Newton at this scale.
        Observation is then telling us about a piece we haven't surfaced.

  The math doesn't choose between these — but it shows the framework has
  a real conceptual question at the substance-velocity-cap-for-orbits level,
  with measurable consequences for engine-time predictions.
""")
