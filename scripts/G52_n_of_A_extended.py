"""
G52: Extended n(A) exploration with A in [A_0, 1]
==================================================

Sean's correction: A can't be zero anywhere. The minimum A in any physical
spacetime is A_0 (cosmic baseline). So the framework's k(A) family lives
on A ∈ [A_0, 1].

Re-running the n(A) candidate exploration with this corrected domain.
Looking for forms that:
  - Stay n ≥ 1 across [A_0, 1]
  - Have clean structural values at landmarks (thirds-of-A, A_0, A=1)
  - Connect the framework's two committed values (n=1 weak, n=2 strong)
  - Have a structural meaning (count, depth, smooth interpolation)

Candidates explored beyond G51:
  - Forms tied to A_0 (n = 1 at A=A_0, n = 2 at A=1)
  - Forms tied to thirds-of-A
  - Forms motivated by "depth" or "count" interpretation
"""

import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8')

A_0 = 1.0 / (12 * np.pi)
A_min = A_0  # Sean's correction
A_max = 0.99999  # avoid singular A=1

# --- Candidate n(A) forms (all designed to give n ≥ 1 on [A_0, 1]) ---

def n_const1(A):       return 1.0
def n_const2(A):       return 2.0
def n_lin(A):          return 1 + A
def n_quad(A):         return 1 + A**2
def n_sqrt(A):         return 1 + np.sqrt(A)
def n_norm_lin(A):     return 1 + (A - A_0)/(1 - A_0)  # n=1 at A_0, n=2 at A=1
def n_norm_quad(A):    return 1 + ((A - A_0)/(1 - A_0))**2
def n_norm_sqrt(A):    return 1 + np.sqrt((A - A_0)/(1 - A_0))
def n_norm_smooth(A):  # smooth s-curve
    x = (A - A_0)/(1 - A_0)
    return 1 + 3*x**2 - 2*x**3  # smoothstep from 0 to 1
def n_thirds(A):       return 1 + 3 * (A - A_0)/(1 - A_0)  # n=1 at A_0, n=4 at A=1
def n_anchor_A0_n2(A): return 1 + A/A_0  # n=2 at A=A_0
def n_two_minus(A):    return 2 - (1-A)/(1-A_0)  # mirror

candidates = [
    ("n = 1 (const)",           n_const1),
    ("n = 2 (const)",           n_const2),
    ("n = 1 + A",               n_lin),
    ("n = 1 + A²",              n_quad),
    ("n = 1 + √A",              n_sqrt),
    ("n = 1 + (A-A_0)/(1-A_0)", n_norm_lin),
    ("n = 1 + ((A-A_0)/(1-A_0))²", n_norm_quad),
    ("n = 1 + √((A-A_0)/(1-A_0))", n_norm_sqrt),
    ("n = smoothstep(A_0,1)",   n_norm_smooth),
    ("n = 2 - (1-A)/(1-A_0)",   n_two_minus),
]

def k_value(A, n):
    return (1 - A) * (1 - A**2)**n

# --- Predictions at landmarks ---
landmarks = [
    ("A_0 (baseline)",   A_0),
    ("1/3 (ISCO)",       1/3),
    ("1/2 (NEC region)", 1/2),
    ("2/3 (photon sph)", 2/3),
    ("0.99 (near horiz)", 0.99),
    ("1 (horizon)",      A_max),
]

print("=" * 100)
print("G52: n(A) candidate exploration on A ∈ [A_0, 1]")
print("=" * 100)
print(f"\nA_0 = 1/(12π) = {A_0:.6f}")

# --- n values at landmarks ---
print("\n" + "-" * 100)
print("n(A) at structural landmarks:")
print("-" * 100)
header = f"{'Candidate':<32}"
for name, A in landmarks:
    header += f" {name:>14}"
print(header)
for label, n_func in candidates:
    row = f"{label:<32}"
    for name, A in landmarks:
        try:
            row += f" {n_func(A):>14.4f}"
        except:
            row += f" {'ERROR':>14}"
    print(row)

# --- G1 ringdown prediction at horizon ---
print("\n" + "-" * 100)
print("G1 ringdown: τ/τ_GR = (9/5)^(n_horizon/2)")
print("LIGO target: spinless analog of observed BBH ringdown; n=1 gives 1.342, n=2 gives 1.800")
print("-" * 100)
for label, n_func in candidates:
    n_h = n_func(A_max)
    print(f"  {label:<32} n_horizon = {n_h:>7.3f}   τ/τ_GR = {(9/5)**(n_h/2):>7.3f}")

# --- Weak-field departure at A_0 ---
print("\n" + "-" * 100)
print("Weak-field test (G17 g_rr 2nd-order ratio): evaluated at A = A_0 (baseline)")
print("Effective ratio = 1 + n(A_0)")
print("-" * 100)
for label, n_func in candidates:
    n_baseline = n_func(A_0)
    print(f"  {label:<32} n(A_0) = {n_baseline:>7.4f}   STAM:GR (2nd order) = {1 + n_baseline:>7.4f}")

# --- k(A) shape comparison ---
print("\n" + "-" * 100)
print("k(A) values across the range (lower k = larger g_rr = stronger STAM departure)")
print("-" * 100)
A_samples = [A_0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
header = f"{'Candidate':<32}"
for A in A_samples:
    header += f" {f'k(A={A:.2f})':>10}"
print(header)
for label, n_func in candidates:
    row = f"{label:<32}"
    for A in A_samples:
        try:
            n = n_func(A)
            k = k_value(A, n)
            row += f" {k:>10.4e}"
        except:
            row += f" {'ERROR':>10}"
    print(row)

# --- Look for structural alignment: thirds-of-A check ---
print("\n" + "=" * 100)
print("Structural alignment check — does n take 'clean' values at thirds-of-A landmarks?")
print("=" * 100)
print(f"""
For each candidate, list (n at A=1/3, n at A=1/2, n at A=2/3, n at A=1):

  n = 1 + A:           (4/3, 3/2, 5/3, 2)
    Values follow thirds: 4/3, 5/3, 6/3=2. NEC at 1/2 gives n=3/2.
    Connects n=1 (weak) and n=2 (horizon) endpoints.
    Thirds-of-A landmarks → thirds-of-n values. Suggestive structural alignment.

  n = 1 + (A−A_0)/(1−A_0):  (1.31, 1.49, 1.66, 2)
    Anchored at n=1 at the cosmic baseline (not at fictitious A=0).
    Endpoints: n(A_0)=1, n(1)=2. Same horizon prediction.
    More "framework-respecting" since it doesn't reach below A_0.
    But thirds-of-A landmarks no longer give clean thirds-of-n.

  n = 1 + ((A−A_0)/(1−A_0))²:  (1.10, 1.22, 1.43, 2)
    Quadratic version. Stays close to 1 at low A (stealthy).
    Cleanest weak-field GR-like behavior.

  n = 1 + 3(A−A_0)/(1−A_0):  (1.94, 2.46, 2.98, 4)
    n=2 reached at A near 1/3, n=3 at A near 2/3, n=4 at horizon.
    Cleanest thirds-of-n alignment at thirds-of-A.
    BUT: gives τ/τ_GR = (9/5)² = 3.24 at horizon — outside current LIGO bracket for n=2.

  n = 1 + A_0/(1-A): bridges from 1 at A=A_0 to ∞ near horizon. Diverges.
""")

# --- Look at specific physically motivated test: at A = A_0, what is the V_3 connection? ---
print("\n" + "-" * 100)
print("V_3 potential connection: V_3(A) = α/A + β/(1-A) with α/β = (A_0/(1-A_0))²")
print("V_3' at A_0 = 0 (potential minimum). V_3'' at A_0 = curvature.")
print("Does any n(A) form connect to V_3 structure?")
print("-" * 100)
# V_3'' at A_0
# V_3 = α/A + β/(1-A)
# V_3' = -α/A² + β/(1-A)²
# V_3'' = 2α/A³ + 2β/(1-A)³
# With α/β = (A_0/(1-A_0))² → α = β(A_0/(1-A_0))²
# V_3''(A_0) = 2 β (A_0/(1-A_0))² / A_0³ + 2β/(1-A_0)³
#            = 2β / [A_0 (1-A_0)²] + 2β/(1-A_0)³
#            = 2β/(1-A_0)² × [1/A_0 + 1/(1-A_0)]
#            = 2β/(1-A_0)² × 1/[A_0(1-A_0)]
#            = 2β/[A_0 (1-A_0)³]

V_3pp_factor = 2 / (A_0 * (1-A_0)**3)
print(f"  V_3''(A_0) = 2β/(A_0(1-A_0)³)")
print(f"  Coefficient (without β): {V_3pp_factor:.4f}")
print(f"  Ratio (1-A_0)²/A_0 = {(1-A_0)**2/A_0:.4f} (appears in α/β)")

# --- Listening ---
print("\n" + "=" * 100)
print("What the math is telling us")
print("=" * 100)
print(f"""
  Several n(A) candidates satisfy the framework's constraints. The most
  structurally clean by different criteria:

  1. **n = 1 + A** — thirds-of-A in n. Clean horizon=2, clean weak-field
     coefficient = 2. Bridges n=1 and n=2 commitments smoothly. Best
     structural alignment but technically reaches "n=1" at A=0 which
     doesn't exist physically.

  2. **n = 1 + (A−A_0)/(1−A_0)** — anchored at A_0 (the physical baseline).
     Same horizon behavior (n=2 at A=1). Linear interpolation from baseline
     to saturation. Most "framework-respecting" since it acknowledges A_0
     as the floor.

  3. **n = 1 + ((A−A_0)/(1−A_0))²** — quadratic anchored form. Stays close
     to n=1 at low A for the longest, then ramps up to 2 at horizon.
     Most "stealthy" — weak-field tests pass trivially.

  4. **n = 1 + 3(A−A_0)/(1−A_0)** — thirds alignment with n values 2, 3,
     4 at thirds. But gives τ/τ_GR = 3.24 at horizon, likely outside
     current LIGO bracket.

  The framework's choice depends on which structural features are
  load-bearing:
    - If weak-field stealthiness is most important → quadratic form (#3)
    - If thirds-of-A alignment is most important → linear (#1) or
      thirds form (#4)
    - If A_0-anchoring is most important → normalized linear (#2)

  All candidates 1, 2, 3 give horizon n=2 → G1 ringdown prediction 1.800.
  Distinguishing them requires precision at intermediate A values:
    - Late-inspiral GW chirp shape sensitive to n at A ~ 0.1-0.5
    - QNM higher-order overtones sensitive to n profile
    - PSR-binary Shapiro at high precision sensitive to n(small A)
""")
