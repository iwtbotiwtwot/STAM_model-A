"""
G50: Target 1.74 s — what does the math need for STAM-A only?
==============================================================

Sean's targeting request. Compute what r_emit (the radius where light first
emerges from the merger event) gives engine time exactly 1.74 s under
STAM-A only (no matter), assuming:

  - Force from explosion propagates at v_substance(A_0) baseline from R_s to
    r_donut_max, takes t_meet = r_donut_max / v_substance(A_0).
  - Light propagates at v_substance(A_static(r)) through the A profile
    A_static(r) = R_s/r near the merger remnant.
  - GW is born at r_donut_max at t_meet, propagates outward at baseline.

Engine time = light's transit from r_emit to r_donut_max minus t_meet.

Targeting: solve for r_emit such that engine time = 1.74 s.

This is a fitting exercise, NOT a derivation. The framework's goal is to
identify a STRUCTURAL meaning for whatever r_emit value falls out.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

def brentq(f, a, b, xtol=1e-12, maxiter=200):
    """Simple bisection method as a stand-in for scipy.optimize.brentq."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError(f"Function has same sign at endpoints: f({a})={fa}, f({b})={fb}")
    for _ in range(maxiter):
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < xtol or (b - a) / 2 < xtol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return (a + b) / 2

c = 299792458.0
G = 6.67430e-11
M_sun = 1.989e30
A_0 = 1.0 / (12 * np.pi)

M_total = 2.7 * M_sun
R_s = 2 * G * M_total / c**2

target_engine = 1.74

def v_substance(A):
    return c * (1 - A)**2 * (1 + A)

v_sub_baseline = v_substance(A_0)

# Two threshold framings
disc_G48 = (1 - A_0)**2 - 4 * A_0
x_thr_G48 = (1 - A_0 - np.sqrt(disc_G48)) / 2
r_thr_G48 = R_s / (2 * x_thr_G48)
tau_chirp_G48 = (5/64) * c**5 * r_thr_G48**4 / (G**3 * M_total**3)
r_donut_G48 = v_sub_baseline * tau_chirp_G48

disc_G49 = 1 - 4 * A_0
x_thr_G49 = (1 - np.sqrt(disc_G49)) / 2
r_thr_G49 = R_s / (2 * x_thr_G49)
tau_chirp_G49 = (5/64) * c**5 * r_thr_G49**4 / (G**3 * M_total**3)
r_donut_G49 = v_sub_baseline * tau_chirp_G49

def light_transit(r_emit, r_donut_max, n_pts=200000):
    """Integrate light's transit from r_emit to r_donut_max through A_static profile."""
    r_grid = np.geomspace(r_emit, r_donut_max, n_pts)
    A_grid = R_s / r_grid
    A_grid = np.clip(A_grid, 0, 0.999999)
    v_grid = v_substance(A_grid)
    return np.trapezoid(1.0 / v_grid, r_grid)

def engine_time(r_emit, r_donut_max):
    t_meet = r_donut_max / v_sub_baseline
    t_light = light_transit(r_emit, r_donut_max)
    return t_light - t_meet

print("=" * 78)
print("G50: Target engine time = 1.74 s, find r_emit (STAM-A only)")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km = {R_s:.2f} m")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  v_substance(A_0) = {v_sub_baseline/c:.6f}c")
print(f"\n  Target engine time: {target_engine} s")

for label, r_donut_max, tau_chirp in [
    ("G48 (cap with A_total)", r_donut_G48, tau_chirp_G48),
    ("G49 (cap with A_static)", r_donut_G49, tau_chirp_G49),
]:
    print("\n" + "-" * 78)
    print(f"Framing: {label}")
    print("-" * 78)
    print(f"  r_donut_max = {r_donut_max/c:.4f} ls = {r_donut_max:.3e} m = {r_donut_max/R_s:.1f} R_s")
    print(f"  τ_chirp = {tau_chirp:.4f} s")

    t_meet = r_donut_max / v_sub_baseline
    print(f"  t_meet (force at baseline) = {t_meet:.4f} s")

    # Sweep r_emit to find the one that gives target engine time
    # r_emit must be > R_s; sweep ε = r_emit/R_s - 1 from very small to ~0.5
    def f(eps):
        r_emit = R_s * (1 + eps)
        return engine_time(r_emit, r_donut_max) - target_engine

    # Find bracket
    try:
        eps_low = 1e-8
        eps_high = 0.1
        f_low = f(eps_low)
        f_high = f(eps_high)
        if f_low * f_high > 0:
            print(f"  No bracket found: f({eps_low}) = {f_low}, f({eps_high}) = {f_high}")
            continue
        eps_sol = brentq(f, eps_low, eps_high, xtol=1e-12)
        r_emit_sol = R_s * (1 + eps_sol)
        A_emit = R_s / r_emit_sol
        v_emit = v_substance(A_emit)

        print(f"\n  Target solution:")
        print(f"    ε = r_emit/R_s - 1 = {eps_sol:.6e}")
        print(f"    r_emit = R_s × (1 + ε) = {r_emit_sol:.6e} m")
        print(f"    r_emit - R_s = {r_emit_sol - R_s:.4f} m ({(r_emit_sol - R_s)*100:.2f} cm)")
        print(f"    A at r_emit = {A_emit:.7f}")
        print(f"    1 − A at r_emit = {1 - A_emit:.4e}")
        print(f"    v_light at r_emit = {v_emit/c:.4e} c = {v_emit:.3e} m/s")

        # Check structural ratios
        u_max = r_donut_max / R_s
        ratio_to_inverse_2umax = eps_sol * 2 * u_max
        ratio_to_A0 = eps_sol / A_0
        ratio_to_RsCEngine = eps_sol / (R_s/(c * target_engine))

        print(f"\n  Structural checks for ε = {eps_sol:.4e}:")
        print(f"    ε × 2u_max = {ratio_to_inverse_2umax:.4f}   (1.0 if ε = R_s/(2 r_donut_max))")
        print(f"    ε / A_0 = {ratio_to_A0:.4e}    (1.0 if ε = A_0)")
        print(f"    ε × c × τ_engine / R_s = {(eps_sol * c * target_engine / R_s):.4f}")
        print(f"    ε × c × τ_chirp / R_s = {(eps_sol * c * tau_chirp / R_s):.4f}")
        print(f"    R_s/(2 r_donut_max) = {R_s/(2*r_donut_max):.4e}")
        print(f"    R_s/(2 × c × τ_engine) = {R_s/(2*c*target_engine):.4e}")

        # Verify
        engine_check = engine_time(r_emit_sol, r_donut_max)
        print(f"\n  Verification: engine_time(r_emit_sol, r_donut_max) = {engine_check:.6f} s")

    except Exception as e:
        print(f"  Error: {e}")

# Sweep across r_emit to show the function shape
print("\n" + "=" * 78)
print("Sweep r_emit to show engine time landscape (G49 r_donut_max)")
print("=" * 78)
print(f"\n  {'ε':>14} {'r_emit/R_s':>14} {'A at r_emit':>14} {'engine (s)':>14}")
for eps in [1e-8, 1e-7, 1e-6, 7.45e-6, 1e-5, 5e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5]:
    r_emit = R_s * (1 + eps)
    A_emit = R_s / r_emit
    try:
        engine = engine_time(r_emit, r_donut_G49)
        print(f"  {eps:>14.4e} {(1+eps):>14.8f} {A_emit:>14.7f} {engine:>14.4f}")
    except Exception as e:
        print(f"  {eps:>14.4e}  ERROR: {e}")

print("\n" + "=" * 78)
print("Reading the math")
print("=" * 78)
print(f"""
  The targeted r_emit is very close to R_s — within ~6 cm above the merger
  surface for GW170817 (where R_s = {R_s:.0f} m = {R_s/1000:.2f} km).

  Physically: light effectively emerges from a thin shell just outside the
  saturation surface where A is still essentially 1 but light can begin to
  propagate (v_substance > 0).

  The structural question: what determines ε? Candidates checked above —
  see which ratio (if any) lands cleanly on 1.0. If ε = R_s/(2 r_donut_max)
  is the right structural relationship, then engine time = r_donut_max × (1−A_0)/c
  approximately, scaling linearly with r_donut_max → linearly with M_total.

  This is what the targeting reveals — STAM-A alone CAN produce 1.74 s
  for GW170817 if r_emit lives in this specific thin near-saturation shell.
  The question is whether this shell has structural meaning.
""")
