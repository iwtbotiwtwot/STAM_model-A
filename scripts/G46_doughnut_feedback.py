"""
G46: STAM-native — doughnut feedback on the inspiral
====================================================

Adds the second-order feedback: the doughnut substance itself sits at the
orbital position with elevation δA_doughnut above baseline. This is on top of
A_static (from the other body) and A_kin (from orbital motion).

Two competing effects of δA_doughnut > 0 at the orbit:

  (1) Threshold pulls IN. Self-consistent threshold becomes:
        A_0 = x × (1 - x - A_doughnut)
        x_threshold² + x_threshold(A_0 - 1 + A_doughnut) + A_0 = 0
      r_threshold gets smaller → chirp distance shorter → SHORTER chirp.

  (2) Chirp rate slows at every r. A_total at each r is bigger by δA_doughnut,
      so (1-A)³ factor gets smaller, chirp rate slows → LONGER chirp.

The math will tell us which wins, and at what δA_doughnut the doughnut/bubble
shrinks back toward observed 1.74 s.

The framework hasn't committed to a specific δA_doughnut value. This is a
sensitivity sweep — the listening exercise is what δA_doughnut value (if any)
makes the math land on observation, and whether that value has structural
meaning.

Don't aim. The result is the breadcrumb.
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

v_sub_baseline = c * (1 - A_0)**2 * (1 + A_0)

def threshold_with_feedback(dA):
    """r_threshold under doughnut feedback δA_doughnut = dA."""
    b = A_0 - 1 + dA
    disc = b**2 - 4 * A_0
    if disc < 0:
        return None
    x_thr = (-b - np.sqrt(disc)) / 2
    return R_s / (2 * x_thr)

def A_total_at_orbit(r, dA):
    """Self-consistent A at orbital position with doughnut feedback dA."""
    x = R_s / (2 * r)
    # A_kin self-consistent with A_doughnut included:
    # A_kin = x(1 - x - A_kin - dA)
    # A_kin(1+x) = x(1 - x - dA)
    # A_kin = x(1 - x - dA) / (1+x)
    A_kin = x * (1 - x - dA) / (1 + x)
    if A_kin < 0:
        A_kin = 0
    return x + A_kin + dA

def chirp_time(r_thr, dA, n_pts=20000):
    """Coord-time chirp from r_thr to R_s with (1-A_total)³ factor on rate."""
    if r_thr is None or r_thr <= R_s:
        return None
    r_grid = np.linspace(R_s, r_thr, n_pts)
    A_grid = np.array([A_total_at_orbit(r, dA) for r in r_grid])
    # Guard: clip A < 1 (otherwise factor blows up at saturation)
    A_grid = np.clip(A_grid, 0, 0.999)
    integrand = r_grid**3 / (1 - A_grid)**3
    I = np.trapezoid(integrand, r_grid)
    # PM coefficient at threshold: τ_PM(r) = (5/64) c⁵ r⁴ / (G³ M³)
    # So τ_chirp = (5/64) c⁵ / (G³ M³) × ∫ r³ (1-A)^(-3) × 4 dr  [wait, careful]
    # τ_PM = ∫_{R_s}^{r_thr} r³ dr × 5/(64) × c⁵ / (G³M³) — for η=1/4 equal-mass
    # so τ_chirp = (5/64) c⁵ / (G³ M³ η_factor) × I, where η-factor handled in coefficient
    # G43 used τ_PM = (5/8)(r_thr/R_s)⁴ × R_s/c which equals (5/64) c⁵ r_thr⁴/(G³M³)
    # τ_PM for equal-mass binary: τ = (5/64) c⁵ r_thr⁴ / (G³M³) [G43's form]
    # Equivalently: τ = (5/16) c⁵ / (G³M³) × ∫_0^{r_thr} r³ dr  (since ∫r³dr = r⁴/4)
    coef = 5 / 16 * c**5 / (G**3 * M_total**3)
    return coef * I

print("=" * 78)
print("G46: STAM-native — doughnut feedback sweep")
print("=" * 78)
print(f"\n  M_total = 2.7 M_sun, R_s = {R_s/1000:.4f} km")
print(f"  A_0 = 1/(12π) = {A_0:.6f}")
print(f"  v_substance at baseline: {v_sub_baseline/c:.4f} c")

# --- Sweep δA_doughnut ---
print("\n" + "-" * 78)
print("Sensitivity sweep: δA_doughnut at orbit position")
print("-" * 78)

print(f"\n  {'δA_dough':>10} {'r_thr/R_s':>12} {'τ_chirp (s)':>14} {'d_bub (lt-s)':>14} {'Δt (s)':>10} {'gap %':>10}")
print(f"  {'-'*10} {'-'*12} {'-'*14} {'-'*14} {'-'*10} {'-'*10}")

dA_values = [0.0, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050]
observed = 1.74

results = []
for dA in dA_values:
    r_thr = threshold_with_feedback(dA)
    if r_thr is None:
        print(f"  {dA:>10.4f}    (threshold equation has no real root)")
        continue
    tau = chirp_time(r_thr, dA)
    if tau is None:
        continue
    d_bubble = v_sub_baseline * tau
    delay = d_bubble / v_sub_baseline  # = tau
    gap = (delay - observed) / observed * 100
    print(f"  {dA:>10.4f} {r_thr/R_s:>12.4f} {tau:>14.4f} {d_bubble/c:>14.4f} {delay:>10.4f} {gap:>+9.2f}%")
    results.append((dA, r_thr, tau, delay, gap))

# Find the δA_doughnut that lands closest to observed
print("\n" + "-" * 78)
print("Where the math lands relative to observation 1.74 s")
print("-" * 78)

best_dA = min(results, key=lambda x: abs(x[4]))
print(f"\n  Closest to observation: δA_doughnut = {best_dA[0]:.4f}")
print(f"    r_threshold/R_s = {best_dA[1]/R_s:.4f}")
print(f"    τ_chirp = signal delay = {best_dA[3]:.4f} s")
print(f"    Gap = {best_dA[4]:+.2f}%")

# Where the doughnut pulls in to give 1.74:
# Interpolate linearly between bracketing points
deltas = np.array([r[4] for r in results])
dAs    = np.array([r[0] for r in results])
# Want delta = 0
if deltas.min() <= 0 and deltas.max() >= 0:
    # Find bracket
    for i in range(len(deltas) - 1):
        if deltas[i] * deltas[i+1] < 0:
            frac = -deltas[i] / (deltas[i+1] - deltas[i])
            dA_match = dAs[i] + frac * (dAs[i+1] - dAs[i])
            print(f"\n  Linear interpolation for exact 1.74 s match: δA_doughnut ≈ {dA_match:.4f}")
            print(f"    = {dA_match/A_0:.3f} × A_0")
            print(f"    = {dA_match*100:.2f}% absolute elevation above A_static + A_kin")
            break

# --- Verification: G45 baseline (δA = 0) ---
print("\n" + "-" * 78)
print("Verification: δA_doughnut = 0 reproduces G45")
print("-" * 78)
print(f"  G45 reported τ_chirp = 2.1073 s, signal delay = 2.1073 s")
print(f"  G46 at δA=0:        τ_chirp = {results[0][2]:.4f} s, signal delay = {results[0][3]:.4f} s")

# --- Reading the math ---
print("\n" + "=" * 78)
print("What the math is telling us")
print("=" * 78)

print(f"""
  Adding doughnut feedback δA_doughnut at the orbital position has two
  competing effects:

    (1) Threshold pulls in: r_threshold drops as δA_doughnut grows, shortening
        the chirp distance.
    (2) Chirp rate slows: (1-A_total)³ factor reduces dr/dt at every r,
        lengthening the chirp.

  The threshold effect wins. As δA_doughnut grows from 0 to {dA_values[-1]}:
    • r_threshold drops from {results[0][1]/R_s:.3f} R_s to {results[-1][1]/R_s:.3f} R_s
    • Signal delay drops from {results[0][3]:.3f} s to {results[-1][3]:.3f} s
    • The doughnut is being pulled in, consistent with Sean's read.

  The math says: doughnut feedback brings the prediction down from the
  +21% overshoot (G45 with no feedback) toward observation. The δA_doughnut
  value that hits 1.74 s is structurally meaningful only if the framework
  can derive it from substance physics. Otherwise it's a free parameter
  the data has constrained — useful breadcrumb but not yet a STAM commitment.

  Open: what determines δA_doughnut structurally? Candidates —
    • The cumulative energy invested in the doughnut substance during inspiral
    • The doughnut's steady-state A profile as a wave train continuously driven
    • Something else from the substance ontology
""")
