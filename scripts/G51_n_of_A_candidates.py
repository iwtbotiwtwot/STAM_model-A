"""
G51: Scratch test for n(A) candidate forms
============================================

Sean's concern: n=A directly fails because n must satisfy n ≥ 1 across the
framework's k(A) family. So n(A) must be a function that stays ≥ 1.

Goals:
  - Explore candidate functional forms for n(A) that satisfy n ≥ 1 everywhere
  - Compute the framework's predictions at key A landmarks under each
  - Look for forms with structurally clean values (matching STAM commitments)

Candidates explored:
  1. n = 1 (constant, simplest)
  2. n = 2 (constant, current commitment)
  3. n = 1 + A (linear bridge from 1 at weak field to 2 at horizon)
  4. n = 1 + A² (quadratic bridge)
  5. n = 2 - A (decreasing from 2 to 1)
  6. n = 1/(1-A) (diverges at horizon)
  7. n = 1 + ln(1/(1-A)) (logarithmic divergence)
  8. n = 1 + 3A (linear, steeper)
  9. n = 1 + 12π × A × A_0 = 1 + A (same as #3 in structural form)
 10. n = 2 - A² (decreasing quadratically)
 11. n = 1 + A^3 (very flat at low A)
 12. n = 1 + A/(1-A) (diverges at horizon)

Tests at structural landmarks:
  - A = 0       (weak field; recovers GR)
  - A = A_0     (cosmic baseline)
  - A = 1/3     (ISCO landmark)
  - A = 1/2     (NEC crossover for n=1)
  - A = 2/3     (photon sphere landmark)
  - A = 1       (horizon)
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

A_0 = 1.0 / (12 * np.pi)

# --- Candidate n(A) forms ---
def n_const1(A):       return 1.0
def n_const2(A):       return 2.0
def n_linear(A):       return 1 + A
def n_quad(A):         return 1 + A**2
def n_descending(A):   return 2 - A
def n_descend_quad(A): return 2 - A**2
def n_inv(A):          return 1.0 / max(1 - A, 1e-12)
def n_log(A):          return 1 + np.log(1.0 / max(1 - A, 1e-12))
def n_3A(A):           return 1 + 3 * A
def n_cubic(A):        return 1 + A**3
def n_horizon_div(A):  return 1 + A / max(1 - A, 1e-12)
def n_12piA(A):        return 1 + 12 * np.pi * A * A_0   # = 1 + A
def n_n2_at_A0(A):
    # n = 1 + A/A_0 gives n=1 at A=0, n=2 at A=A_0, n=12π+1 at A=1
    return 1 + A / A_0

candidates = [
    ("n = 1 (const)",                n_const1),
    ("n = 2 (const)",                n_const2),
    ("n = 1 + A",                    n_linear),
    ("n = 1 + A²",                   n_quad),
    ("n = 2 - A",                    n_descending),
    ("n = 2 - A²",                   n_descend_quad),
    ("n = 1 + 3A",                   n_3A),
    ("n = 1 + A³",                   n_cubic),
    ("n = 1 + A/(1-A)",              n_horizon_div),
    ("n = 1 + ln(1/(1-A))",          n_log),
    ("n = 1 + A/A_0",                n_n2_at_A0),
    ("n = 1/(1-A)",                  n_inv),
]

# --- Predictions ---
def G1_ringdown_ratio(n):
    """G1 prediction: τ/τ_GR = (9/5)^(n/2). Evaluated at the relevant A (horizon)."""
    return (9/5)**(n/2)

def g_rr_2nd_order_ratio(n):
    """G17: STAM:GR second-order coefficient ratio = n+1 (at small A)."""
    return n + 1

def k_value(A, n):
    """k(A) = (1-A)(1-A²)^n"""
    return (1 - A) * (1 - A**2)**n

def NEC_crossover_check(n_func):
    """F3 NEC crossover: B(A) = A²(... depending on n). For constant n=1, A*=1/2.
    Approximate: find where d/dA[1/k] changes structure."""
    # Sample 1/k and look for inflection-like feature
    A_grid = np.linspace(0.01, 0.99, 200)
    inv_k = np.array([1.0 / k_value(A, n_func(A)) for A in A_grid])
    # Compute d²/dA²
    d2 = np.gradient(np.gradient(inv_k, A_grid), A_grid)
    # NEC crossover is where d²/dA² hits some specific value (varies by n)
    # For now just report the inflection point (zero of d²/dA²)
    sign_changes = np.where(np.diff(np.sign(d2)))[0]
    if len(sign_changes) > 0:
        return A_grid[sign_changes[0]]
    return None

print("=" * 90)
print("G51: Scratch test — n(A) candidate forms")
print("=" * 90)
print(f"\nA_0 = 1/(12π) = {A_0:.6f}")

# Header table
print("\n" + "-" * 90)
print("n(A) values at structural landmarks")
print("-" * 90)
print(f"{'Candidate':<24} {'n(0)':>8} {'n(A_0)':>10} {'n(1/3)':>10} {'n(1/2)':>10} {'n(2/3)':>10} {'n(1)':>10}")
for label, n_func in candidates:
    try:
        n_vals = [n_func(A) for A in [0, A_0, 1/3, 1/2, 2/3, 0.9999]]
        print(f"{label:<24} {n_vals[0]:>8.3f} {n_vals[1]:>10.4f} {n_vals[2]:>10.3f} {n_vals[3]:>10.3f} {n_vals[4]:>10.3f} {n_vals[5]:>10.3f}")
    except Exception as e:
        print(f"{label:<24} ERROR: {e}")

# G1 ringdown predictions
print("\n" + "-" * 90)
print("G1 ringdown prediction: τ/τ_GR = (9/5)^(n_horizon/2) — evaluated at A=1")
print("(LIGO precision currently bracket-level; constraints tightening)")
print("-" * 90)
print(f"{'Candidate':<24} {'n at horizon':>14} {'τ/τ_GR':>12}")
for label, n_func in candidates:
    try:
        n_h = n_func(0.9999)
        ratio = G1_ringdown_ratio(n_h)
        print(f"{label:<24} {n_h:>14.3f} {ratio:>12.3f}")
    except:
        pass

# Weak-field departure (g_rr 2nd order)
print("\n" + "-" * 90)
print("G17 g_rr 2nd-order ratio at A→0: ratio = 1 + n(0)")
print("-" * 90)
print(f"{'Candidate':<24} {'n(0)':>10} {'STAM:GR ratio':>16}")
for label, n_func in candidates:
    try:
        n_w = n_func(0)
        print(f"{label:<24} {n_w:>10.3f} {1 + n_w:>16.3f}")
    except:
        pass

# Bridge term structure: at A_0
print("\n" + "-" * 90)
print("n at cosmic baseline A_0 — relates to V_3, bridge term, and SU-ruler reading")
print("-" * 90)
print(f"{'Candidate':<24} {'n(A_0)':>10} {'(1-A_0²)^n':>14} {'k(A_0)':>14}")
for label, n_func in candidates:
    try:
        n_A0 = n_func(A_0)
        suppression = (1 - A_0**2)**n_A0
        k_val = (1 - A_0) * suppression
        print(f"{label:<24} {n_A0:>10.4f} {suppression:>14.6f} {k_val:>14.6f}")
    except:
        pass

# Proper-time divergence behavior at horizon
print("\n" + "-" * 90)
print("Proper-time integral ∫dA/k near A=1 — divergence structure")
print("-" * 90)
print(f"{'Candidate':<24} {'n at A=1':>12} {'1/k(0.99)':>14} {'1/k(0.999)':>14} {'1/k(0.9999)':>14}")
for label, n_func in candidates:
    try:
        n_h = n_func(0.9999)
        vals = [1.0 / k_value(A, n_func(A)) for A in [0.99, 0.999, 0.9999]]
        print(f"{label:<24} {n_h:>12.3f} {vals[0]:>14.3e} {vals[1]:>14.3e} {vals[2]:>14.3e}")
    except:
        pass

# Identify "clean" forms — those with integer values at landmarks
print("\n" + "=" * 90)
print("Structural cleanliness check")
print("=" * 90)
print("""
Candidates with structurally clean values at landmarks:

  n = 1 + A:
    n(0) = 1, n(A_0) = 1 + 1/(12π), n(1/3) = 4/3, n(1/2) = 3/2,
    n(2/3) = 5/3, n(1) = 2
    → Bridges from n=1 at weak field to n=2 at horizon.
    → Gives τ/τ_GR = √(9/5)^2 = 9/5 at horizon (matches n=2 commitment)
    → Gives τ/τ_GR = √(9/5) at weak field (matches n=1)
    → "Both n=1 and n=2 simultaneously, depending on regime"

  n = 1 + A²:
    Same horizon-limit τ/τ_GR = 9/5
    Slower transition: stays close to 1 at low A

  n = 2 - A:
    n(0) = 2, n(1) = 1. Reverse of n = 1 + A.

  n = 1 + A/A_0:
    n(0) = 1, n(A_0) = 2, n(1) = 1 + 12π ≈ 38.7
    Extreme growth toward horizon — possibly unphysical
    BUT: clean structural meaning at A_0 (n=2 = current commitment)

Most interesting candidate: **n = 1 + A**
  - Continuously interpolates between n=1 (weak field) and n=2 (horizon)
  - Both committed values pop out as endpoint behaviors
  - Resolves the n=1 vs n=2 ambiguity as "both — at different A regimes"
  - G1 ringdown at A=1 → n=2 → τ/τ_GR = 1.8 (current commitment, falsifiable target)
  - Second-order coefficient at weak field → 1+1 = 2 → STAM:GR = 2 (matches n=1)
  - F3 NEC behavior: needs computation
  - G18 entropy: at A=1, n=2, so (2×2)=4 decomposition works
""")
