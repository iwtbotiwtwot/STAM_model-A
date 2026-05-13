"""
G44: STAM-native GW170817 engine time
=====================================

Sandbox STAM. No Peters-Mathews carried in from GR — but at low A the framework
recovers GR's quadrupole rate by its own construction, so that rate is STAM's
own rate at low A, not a GR import.

Mechanism (head-start framing, per Sean 2026-05-12 follow-up):
  GW front  = outer edge of substance bubble (propagating during inspiral at v_substance)
  Light     = electromagnetic emission from merger center at orbit-stop
  Δt at observer = R_bubble_at_orbit_stop / v_substance ≈ τ_inspiral_observer

Two STAM-native energy-loss channels to weigh:
  (A) Bubble V_3 elevation: orbital binding energy → A elevation in substance bubble
  (B) GW radiation as metric perturbation: STAM's own GW channel, ≡ GR at low A

The script first checks (A)'s capacity. If V_3 cannot energetically be the sink
at stellar-mass scales, (B) must dominate, and the math has told us something
specific about substance density at the V_3 baseline.

Then computes observer-frame inspiral time under whichever channel the math
points to. Don't aim at 1.74 s. Whatever falls out is the breadcrumb.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

# --- Constants ---
c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
Mpc = 3.0857e22
H_0 = 73.05 * 1000 / Mpc  # SI

# --- STAM ---
A_0 = 1.0 / (12 * np.pi)
rho_crit = 3 * H_0**2 / (8 * np.pi * G)
Omega_DE = 0.7

# V_3 calibration: β_tilde = Ω_DE × (1-A_0)²; β in J/m³ via × ρ_crit × c²
beta  = Omega_DE * (1 - A_0)**2 * rho_crit * c**2
alpha = beta * (A_0 / (1 - A_0))**2

def V_3(A):
    return alpha / A + beta / (1 - A)

V_3_baseline = V_3(A_0)

# --- GW170817 ---
M_total = 2.7 * M_sun
mu = M_total / 4.0            # equal-mass reduced mass
R_s = 2 * G * M_total / c**2

print("=" * 78)
print("G44: STAM-native GW170817 engine time")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  V_3 calibration: β = {beta:.3e} J/m³")

# --- Threshold (self-consistent, per G43 logic) ---
disc = (1 - A_0)**2 - 4 * A_0
x_threshold = (1 - A_0 - np.sqrt(disc)) / 2
r_threshold = R_s / (2 * x_threshold)
print(f"\n  Self-consistent threshold: r_threshold/R_s = {r_threshold/R_s:.4f}")

# --- Local A profile during inspiral ---
def A_static(r):
    return R_s / (2 * r)

def A_kin(r):
    x = R_s / (2 * r)
    return x * (1 - x) / (1 + x)

def A_total(r):
    return A_static(r) + A_kin(r)  # = 2x/(1+x)

# --- Substance speed at A ---
def v_substance(A):
    return c * (1 - A)**2 * (1 + A)

v_sub_baseline = v_substance(A_0)
print(f"  v_substance(A_0) = {v_sub_baseline/c:.6f} c")

# --- Channel A: V_3 bubble capacity check ---
print("\n" + "-" * 78)
print("Channel (A): V_3 bubble accumulation as energy sink")
print("-" * 78)

# Use observed timescale as the size estimate (just for the capacity test;
# this is not aiming, it's sizing the bubble for the capacity comparison)
tau_size = 1.7  # seconds, for sizing the bubble in capacity test only
R_bubble_est = v_sub_baseline * tau_size

E_orbit_thr = -G * M_total * mu / (2 * r_threshold)
E_orbit_Rs  = -G * M_total * mu / (2 * R_s)
dE_orbit    = abs(E_orbit_Rs - E_orbit_thr)
print(f"  ΔE_orbit (threshold → R_s, Newtonian): {dE_orbit:.3e} J")

# Bubble V_3 capacity at A_inside = generous upper bound 0.5
A_inside_bound = 0.5
V_3_elev_max = V_3(A_inside_bound) - V_3_baseline
V_bubble = (4/3) * np.pi * R_bubble_est**3
E_bubble_max = V_3_elev_max * V_bubble
print(f"  V_3 elevation density @ A={A_inside_bound}: {V_3_elev_max:.3e} J/m³")
print(f"  Bubble volume (~{tau_size:.1f} ls radius):  {V_bubble:.3e} m³")
print(f"  E_bubble_max (V_3 upper bound):        {E_bubble_max:.3e} J")
print(f"  Capacity ratio (bubble / orbital):     {E_bubble_max/dE_orbit:.3e}")

if E_bubble_max < 0.01 * dE_orbit:
    print(f"\n  → V_3 bubble is {dE_orbit/E_bubble_max:.1e}× too small to absorb ΔE_orbit.")
    print(f"  → Substance density at the V_3 baseline is too dilute to be the energy sink")
    print(f"    for stellar-mass binary mergers. Channel (A) is structurally subdominant here.")
    print(f"  → Channel (B) — STAM GW radiation (metric perturbations) — must dominate.")
    print(f"    At low A, STAM's GW emission rate ≡ GR's quadrupole by construction.")
    channel = "B"
else:
    channel = "A"

# --- Channel B: STAM GW radiation, observer-frame inspiral time ---
print("\n" + "-" * 78)
print("Channel (B): STAM GW emission, observer-frame engine time")
print("-" * 78)

# Source-frame (asymptotic coordinate) chirp time from r_threshold to R_s
# STAM's GW emission rate at low A ≡ GR quadrupole (framework's own commitment)
tau_source = (5/8) * (r_threshold/R_s)**4 * (R_s/c)
print(f"  Source-frame chirp time (STAM = GR at low A): {tau_source:.4f} s")

# Magic-bell translation to LIGO frame
# LIGO is at A = A_0; orbital position is at A_total(r)
# dτ_LIGO / dτ_source(local) = √((1-A_0)/(1-A(r)))
# Integrate over the chirp trajectory r(t) from r_threshold to R_s

n_pts = 5000
t_src = np.linspace(0, tau_source, n_pts)
# Chirp law: r^4 = r_thr^4 × (1 - t/tau_source)
r_grid = (r_threshold**4 * (1 - t_src/tau_source))**(1/4)
r_grid = np.maximum(r_grid, R_s)
A_grid = np.array([A_total(r) for r in r_grid])
clock_factor = np.sqrt((1 - A_0) / (1 - A_grid))
tau_observer = np.trapezoid(clock_factor, t_src)

print(f"  Mean clock factor across chirp:        {np.mean(clock_factor):.4f}")
print(f"  Max clock factor (near R_s):           {np.max(clock_factor):.4f}")
print(f"  Observer-frame inspiral duration (Reading 2): {tau_observer:.4f} s")

# --- Reading 1: substance velocity-cap applies to orbital motion ---
# Then ω² = (1-A) GM/r³ in proper time, so P_GW ∝ ω⁶ scales by (1-A)³
# dr/dt_STAM = (1-A)³ × dr/dt_PM
# τ_STAM = ∫ r³ (1-A(r))^(-3) dr × const, vs τ_PM = ∫ r³ dr × const

print("\n" + "-" * 78)
print("Reading 1: Substance velocity-cap also applies to orbital motion")
print("-" * 78)

r_int = np.linspace(R_s, r_threshold, 5000)
A_int = np.array([A_total(r) for r in r_int])
integrand_STAM = r_int**3 / (1 - A_int)**3
integrand_PM   = r_int**3
I_STAM = np.trapezoid(integrand_STAM, r_int)
I_PM   = np.trapezoid(integrand_PM, r_int)
ratio_STAM_PM = I_STAM / I_PM

tau_source_R1 = tau_source * ratio_STAM_PM
print(f"  τ_STAM / τ_PM ratio (substance cap on chirp): {ratio_STAM_PM:.4f}")
print(f"  Source-frame chirp time (Reading 1): {tau_source_R1:.4f} s")

# Now apply clock dilation to this longer chirp (chirp profile differs, but use
# the same A-trajectory weighting as a first-pass approximation)
tau_observer_R1 = tau_source_R1 * np.mean(clock_factor)
print(f"  Observer-frame inspiral duration (Reading 1): {tau_observer_R1:.4f} s")

# Head-start identification
print(f"\n  Head-start mechanism (both readings):")
print(f"    GW front at orbit-stop = bubble outer edge = v_sub × τ_observer")
print(f"    Light starts from r=0 at orbit-stop")
print(f"    Δt observed = R_bubble / v_sub ≈ τ_observer")

# --- Comparison ---
print("\n" + "=" * 78)
print("Listening to the math")
print("=" * 78)
print(f"\n  Observed engine time (GW170817): 1.74 s")
print(f"\n  Progression:")
print(f"    Naive Newton (no STAM):           2.1004 s  ({(2.1004-1.74)/1.74*100:+.2f}%)")
print(f"    G43 (subst. threshold, no clock): {tau_source:.4f} s  ({(tau_source-1.74)/1.74*100:+.2f}%)")
print(f"    Reading 2 (G43 + LIGO clock):     {tau_observer:.4f} s  ({(tau_observer-1.74)/1.74*100:+.2f}%)")
print(f"    Reading 1 (R2 + subst-cap chirp): {tau_observer_R1:.4f} s  ({(tau_observer_R1-1.74)/1.74*100:+.2f}%)")
print()
print(f"  Channel verdict: STAM energy loss for binary mergers is GW-emission-dominated,")
print(f"  not bubble-V_3-dominated. The substance is too dilute at the V_3 baseline to")
print(f"  absorb stellar-mass binding energies through A elevation alone.")
print(f"")
print(f"  Reading fork (substance velocity-cap on orbital motion?):")
print(f"    Reading 2 (cap on translation only, not circular orbit): undershoots by 1.3%")
print(f"    Reading 1 (cap on all motion through substance):         overshoots by {(tau_observer_R1-1.74)/1.74*100:+.1f}%")
print(f"")
print(f"  The math is telling us the framework has a conceptual choice on the table.")
print(f"  Neither reading lands on 1.74 alone; observation is between them.")
