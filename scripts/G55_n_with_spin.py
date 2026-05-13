"""
G55: Combining spin into the n(A) functional form
====================================================

Idea: n depends on both A and spin a. Possible forms:

  Form 1: Direct spin damping multiplied onto Sean's form
    n(A, a) = (1+3A)(1-A/2)(1-a²)

  Form 2: Linear spin damping
    n(A, a) = (1+3A)(1-A/2)(1-a)

  Form 3: Spin in observability factor
    n(A, a) = (1+3A)(1 - A/2 - a²/4)  # spin adds to "hidden" fraction

  Form 4: Quartic spin damping (slower onset)
    n(A, a) = (1+3A)(1-A/2)(1-a^4)

Test: for each form, compute the framework's prediction for
τ_STAM/τ_GR_Kerr at LIGO-relevant spins, see which gives consistency.

Physical intuition: at higher spin, maybe more of the framework's
"hidden" structure becomes observable (frame-dragging brings inner
content outward), OR the structural count reduces (matter rotation
disperses elements). Either way, n should decrease with spin.
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

A_0 = 1.0 / (12 * np.pi)
A_lightring = 2/3  # Schwarzschild light ring under Interpretation B

# --- Spin-dependent Kerr GR ringdown ratio (from earlier) ---
def kerr_QNM_damping_ratio(a):
    return 1 + 0.1*a + 0.4*a**4

# --- Candidate forms for n(A, a) ---
def n_form_1(A, a):
    """(1+3A)(1-A/2)(1-a²)"""
    return (1 + 3*A) * (1 - A/2) * (1 - a**2)

def n_form_2(A, a):
    """(1+3A)(1-A/2)(1-a)"""
    return (1 + 3*A) * (1 - A/2) * (1 - a)

def n_form_3(A, a):
    """(1+3A)(1 - A/2 - a²/4)"""
    return (1 + 3*A) * (1 - A/2 - a**2/4)

def n_form_4(A, a):
    """(1+3A)(1-A/2)(1-a^4)"""
    return (1 + 3*A) * (1 - A/2) * (1 - a**4)

def n_form_baseline(A, a):
    """No spin dependence — Sean's original form"""
    return (1 + 3*A) * (1 - A/2)

forms = [
    ("baseline (Sean): (1+3A)(1-A/2)",          n_form_baseline),
    ("Form 1: × (1-a²)",                         n_form_1),
    ("Form 2: × (1-a)",                          n_form_2),
    ("Form 3: (1+3A)(1-A/2-a²/4)",               n_form_3),
    ("Form 4: × (1-a^4)",                        n_form_4),
]

print("=" * 100)
print("G55: Spin combined into n(A, a)")
print("=" * 100)

print("\nTesting forms at the Schwarzschild light ring (A = 2/3, where eikonal QNMs probe).")
print("Computing n(2/3, a) and the resulting τ_STAM/τ_GR_Kerr for typical LIGO remnant spins.")
print()

# Header
header = f"{'Spin a':>8} {'τ_Kerr/τ_Schw':>14}"
for label, n_func in forms:
    header += f" {label[:24]:>26}"
print(header)
header2 = f"{'':>8} {'':>14}"
for label, n_func in forms:
    header2 += f" {'n_LR | τ_S/τ_K':>26}"
print(header2)
print("-" * 100)

spin_values = [0.0, 0.1, 0.3, 0.5, 0.67, 0.8, 0.9, 0.99]

for a in spin_values:
    spin_ratio = kerr_QNM_damping_ratio(a)
    row = f"{a:>8.2f} {spin_ratio:>14.3f}"
    for label, n_func in forms:
        n_LR = max(n_func(A_lightring, a), 0)
        if n_LR < 0:
            row += f" {'—':>26}"
            continue
        tau_spinless = (9/5)**(n_LR/2)
        tau_over_kerr = tau_spinless / spin_ratio
        row += f" {f'{n_LR:.2f}|{tau_over_kerr:.2f}':>26}"
    print(row)

# --- Compare to LIGO constraints ---
print("\n" + "=" * 100)
print("LIGO consistency check")
print("=" * 100)
print("""
LIGO ringdown constraint: τ_obs/τ_GR_Kerr(observed a) = 1 ± ~10-30% (depending on event/catalog)
For framework to be consistent: τ_STAM/τ_GR_Kerr should be within ~1.0-1.3 at typical spin.

Looking at the typical LIGO remnant spin a ≈ 0.67:""")

a_typical = 0.67
spin_ratio = kerr_QNM_damping_ratio(a_typical)
print(f"\n  At a = {a_typical} (typical LIGO remnant):  τ_Kerr/τ_Schw = {spin_ratio:.3f}")
print()
print(f"{'Form':<36} {'n_LR':>8} {'τ_STAM/τ_Kerr':>16} {'LIGO consistency':>20}")
print("-" * 90)
for label, n_func in forms:
    n_LR = max(n_func(A_lightring, a_typical), 0)
    tau_spinless = (9/5)**(n_LR/2)
    tau_over_kerr = tau_spinless / spin_ratio
    delta_pct = (tau_over_kerr - 1) * 100

    if abs(delta_pct) < 15:
        verdict = "consistent"
    elif abs(delta_pct) < 30:
        verdict = "marginal"
    elif abs(delta_pct) < 60:
        verdict = "tension"
    elif abs(delta_pct) < 100:
        verdict = "strong tension"
    else:
        verdict = "excluded"

    print(f"{label:<36} {n_LR:>8.2f} {tau_over_kerr:>16.3f} {verdict:>20}")

# --- Cross check: predictions at low and high spin ---
print("\n" + "-" * 100)
print("Cross-check at extremes:")
print("-" * 100)
print()
print("  Low-spin remnant (a = 0.1):")
print(f"  τ_Kerr/τ_Schw = {kerr_QNM_damping_ratio(0.1):.3f}")
for label, n_func in forms:
    n_LR = max(n_func(A_lightring, 0.1), 0)
    tau = (9/5)**(n_LR/2) / kerr_QNM_damping_ratio(0.1)
    print(f"    {label:<40} n_LR = {n_LR:.2f}, τ_STAM/τ_Kerr = {tau:.3f}")
print()
print("  High-spin remnant (a = 0.9):")
print(f"  τ_Kerr/τ_Schw = {kerr_QNM_damping_ratio(0.9):.3f}")
for label, n_func in forms:
    n_LR = max(n_func(A_lightring, 0.9), 0)
    tau = (9/5)**(n_LR/2) / kerr_QNM_damping_ratio(0.9)
    print(f"    {label:<40} n_LR = {n_LR:.2f}, τ_STAM/τ_Kerr = {tau:.3f}")

# --- Solve for the spin function that makes prediction match LIGO ---
print("\n" + "=" * 100)
print("Reverse-engineering: what spin damping function gives LIGO-consistent predictions?")
print("=" * 100)
print()
print("For τ_STAM/τ_GR_Kerr ≈ 1 at every spin (no detectable STAM ringdown signature):")
print("  Required n_LR(a) = 2·log(τ_Kerr/τ_Schw)/log(9/5)")
print()
print(f"{'Spin a':>8} {'τ_Kerr/τ_Schw':>16} {'Required n_LR':>16} {'ratio to baseline':>20}")
print("-" * 60)
for a in spin_values:
    spin_ratio = kerr_QNM_damping_ratio(a)
    n_required = 2 * np.log(spin_ratio) / np.log(9/5)
    baseline_n = n_form_baseline(A_lightring, a)  # = 2 always
    ratio = n_required / baseline_n if baseline_n > 0 else 0
    print(f"{a:>8.2f} {spin_ratio:>16.3f} {n_required:>16.3f} {ratio:>20.4f}")

print()
print("""
  Reading: to make STAM match GR exactly, n_LR at a=0 must be 0 (no STAM
  at all), but the framework has structural reasons to want n_LR(a=0) ≥ 1.

  For STAM to give τ_STAM ≈ τ_GR_Kerr at all spins, the function (9/5)^(n_LR/2)
  must equal τ_Kerr/τ_Schw at each spin. This requires:

    n_LR(0) = 0    (no STAM at Schwarzschild — contradicts framework)
    n_LR(0.67) ≈ 0.48
    n_LR(0.9) ≈ 1.02

  The "n_LR rises with spin" behavior of "required n_LR" is opposite to
  what physical motivation (spin damping reducing n) would suggest. So
  no simple n(A,a) form can make STAM match GR exactly.

  Instead, the framework needs to either:
  (a) Accept that LIGO will measure some τ excess and predict the shape
  (b) Find a spin-dependent form where the framework's prediction is small
      but non-zero across the LIGO observable range
""")

# --- What the math is telling us ---
print("=" * 100)
print("What the math is telling us")
print("=" * 100)
print(f"""
  Combining spin into n cleanly:

  Form 1: n = (1+3A)(1-A/2)(1-a²)  [spin damping quadratic in a]
    At a=0.67: n_LR = {n_form_1(2/3, 0.67):.2f}, τ_STAM/τ_GR_Kerr = {(9/5)**(n_form_1(2/3,0.67)/2)/kerr_QNM_damping_ratio(0.67):.3f}
    At a=0:    n_LR = 2.0, τ_STAM/τ_GR_Kerr = 1.8 (Schwarzschild — tension persists)
    At a=0.9:  n_LR = {n_form_1(2/3, 0.9):.2f}, τ_STAM/τ_GR_Kerr ≈ {(9/5)**(n_form_1(2/3,0.9)/2)/kerr_QNM_damping_ratio(0.9):.3f}

  Form 4: n × (1-a^4)  [slower spin onset]
    At a=0.67: n_LR = {n_form_4(2/3, 0.67):.2f}, τ_STAM/τ_GR_Kerr = {(9/5)**(n_form_4(2/3,0.67)/2)/kerr_QNM_damping_ratio(0.67):.3f}
    Comparison-friendly only at high spin.

  Form 3: (1+3A)(1-A/2-a²/4)  [spin in observability factor]
    At a=0.67: n_LR = {n_form_3(2/3, 0.67):.2f}, τ_STAM/τ_GR_Kerr = {(9/5)**(n_form_3(2/3,0.67)/2)/kerr_QNM_damping_ratio(0.67):.3f}

  None of these forms gives LIGO consistency across the full spin range
  while keeping the framework's Schwarzschild commitment (n=2 at A=2/3, a=0).

  The structural tension:
    - Framework wants n_LR(A=2/3, a=0) = 2 (Schwarzschild light ring at n=2)
    - LIGO wants the spin-corrected prediction to be ~1 (no STAM excess)
    - These are incompatible at a=0 (where there's no spin to "damp" away n=2)

  The framework's options:
    1. Commit to the Schwarzschild prediction (τ/τ_GR = 1.8) and accept LIGO
       tension. Wait for tighter data to resolve.
    2. Move to n_LR(A=2/3, a=0) < 1 (e.g., constant n=1 or n(A) with
       n at light ring close to 1). Loses some structural elegance.
    3. Argue that the framework's "spinless τ_STAM" computation is naive
       and the full STAM-Kerr derivation gives smaller departures (open work).

  Sean's spin-combined form helps at HIGH spin but doesn't resolve the
  tension at Schwarzschild/low-spin.
""")
