"""G30 — work out the math showing SU ↔ A_0 connection from K = 1231.350177.

User claim: 1 SU = A_0  (or 1/SU = A_0)
Setup:  SU(z) = K (a z + q z^2)  with  K = 1231.350177
Identity: SU(z) = (z/H)(1 + 3z/20)
So:  K·a = 1/H ,   K·q = 3/(20 H)   →   q/a = 3/20 (independent of K)
A_0 = 1/(12π) ≈ 0.026526 (cosmic baseline)
Spreadsheet sets H/c = 0.000242 z/Mpc → H_0 ≈ 72.55 km/s/Mpc, c/H_0 ≈ 4132 Mpc.

Goal: identify the unit choice in which K = 1231.35 makes the SU↔A_0 claim literal.
"""

import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from math import pi, sqrt

# --- constants ---
A0 = 1.0 / (12.0 * pi)
inv_A0 = 12.0 * pi
c_kms = 299792.458
H0_spreadsheet_zperMpc = 0.000242   # SU constant from "Physical Meaning" sheet
H0_kmsMpc = H0_spreadsheet_zperMpc * c_kms
c_over_H0_Mpc = 1.0 / H0_spreadsheet_zperMpc
c_over_H0_Mly = c_over_H0_Mpc * 3.26156  # Mpc -> Mly
bridge_Mpc = A0 * c_over_H0_Mpc
bridge_Mly = A0 * c_over_H0_Mly

K = 1231.350177

print("=" * 78)
print("STEP 0 — constants")
print("=" * 78)
print(f"A_0 = 1/(12π) = {A0:.8f}")
print(f"1/A_0 = 12π   = {inv_A0:.6f}")
print(f"H_0 (spreadsheet, from 0.000242 z/Mpc) = {H0_kmsMpc:.3f} km/s/Mpc")
print(f"c/H_0 = {c_over_H0_Mpc:.3f} Mpc = {c_over_H0_Mly:.3f} Mly")
print(f"bridge b = A_0 · c/H_0 = {bridge_Mpc:.3f} Mpc = {bridge_Mly:.3f} Mly")
print(f"K (Sean) = {K}")

print()
print("=" * 78)
print("STEP 1 — what does K·a = 1/H pin down?")
print("=" * 78)
print("K·a = 1/H is a unit identity: K is dimensionless number, a carries units")
print("[SU/z], H has units [z/length]. So a [length/SU] = (1/H)[length/z] / K.")
print()
print("Without picking what '1 SU' is, K is free. Try each unit choice:")
print()

unit_choices = [
    ("1 SU = 1 Mpc",       1.0,                    "Mpc"),
    ("1 SU = 1 Mly",       1.0 / 3.26156,          "Mly"),
    ("1 SU = 1 Gly",       1000.0 / 3.26156,       "Gly"),
    ("1 SU = bridge b",    bridge_Mpc,             "Mpc"),
    ("1 SU = c/H_0",       c_over_H0_Mpc,          "Mpc (Hubble distance)"),
    ("1 SU = 1 / (12π) of c/H_0", c_over_H0_Mpc / (12*pi), "Mpc (= bridge again)"),
]

print(f"{'unit choice':45s}  {'a [Mpc/SU]':>14s}  {'K = c/(H_0·a) ':>16s}")
for label, a_in_Mpc_per_SU, _ in unit_choices:
    K_implied = c_over_H0_Mpc / a_in_Mpc_per_SU
    print(f"{label:45s}  {a_in_Mpc_per_SU:14.4f}  {K_implied:16.4f}")

print()
print("None of these reproduces K = 1231.35 directly.")
print()

print("=" * 78)
print("STEP 2 — invert: what 'size of 1 SU' does K = 1231.35 imply?")
print("=" * 78)
# K · a = c/H_0 (in Mpc), so a = (c/H_0) / K
a_implied_Mpc = c_over_H0_Mpc / K
a_implied_Mly = a_implied_Mpc * 3.26156
print(f"a = (c/H_0)/K = {c_over_H0_Mpc:.2f} / {K} = {a_implied_Mpc:.4f} Mpc/SU")
print(f"            = {a_implied_Mly:.4f} Mly/SU")
print()
print("So K = 1231.35 declares  1 SU = " f"{a_implied_Mpc:.4f} Mpc  ≈ {a_implied_Mly:.3f} Mly.")
print()
print("Does this match any natural STAM length?")
print(f"  bridge b      = {bridge_Mpc:.3f} Mpc  ({bridge_Mly:.3f} Mly)")
print(f"  A_0 · bridge  = {A0*bridge_Mpc:.4f} Mpc  ({A0*bridge_Mly:.4f} Mly)")
print(f"  A_0² · c/H_0  = {A0*A0*c_over_H0_Mpc:.4f} Mpc  ({A0*A0*c_over_H0_Mly:.4f} Mly)")
print(f"  bridge/A_0    = c/H_0 = {c_over_H0_Mpc:.2f} Mpc")
print(f"  a/A_0         = {a_implied_Mpc/A0:.4f} Mpc")
print(f"  a/A_0²        = {a_implied_Mpc/A0/A0:.4f} Mpc")
print(f"  a · 12π       = {a_implied_Mpc*12*pi:.4f} Mpc")

print()
print("=" * 78)
print("STEP 3 — the literal reading: '1 SU = A_0' as a unit identity")
print("=" * 78)
print("If we work in dimensionless units where Hubble distance c/H_0 = 1, then:")
print("  D_C(z)/(c/H_0) = z + (3/20)z² + ...")
print("  bridge in those units = A_0 (because b = A_0·c/H_0).")
print()
print("If 1 SU is *defined* as the bridge b, then in Hubble-distance units:")
print("  1 SU = b/(c/H_0) = A_0     ← literal 1 SU = A_0  ✓")
print()
print("Equivalently, in units where 1 SU = bridge = c/(12π·H_0):")
print("  c/H_0 = 12π SU  =  (1/A_0) SU   →   1/SU = A_0 / 1   when normalising c/H_0=1.")
print()
print("So '1 SU = A_0' is true under the normalization:")
print(f"  1 SU ≡ bridge = c/(12π·H_0) ≈ {bridge_Mpc:.2f} Mpc ≈ {bridge_Mly:.2f} Mly")
print()
print("Under this rule, SU(z) (in those units) = D_C(z)/b = (12π/(c/H_0))·D_C(z)")
print("and at low z:  SU(z) ≈ 12π·z + (12π)(3/20)z² = 12π·z(1 + 3z/20)")
print(f"so K·a in this convention = 12π = {12*pi:.6f}, NOT 1231.35.")
print()

print("=" * 78)
print("STEP 4 — where could K = 1231.35 come from?")
print("=" * 78)
print("Notable arithmetic near 1231.35:")
print(f"  144π²            = {144*pi*pi:.4f}       (= 1/A_0²)")
print(f"  144π² · √3/2     = {144*pi*pi*sqrt(3)/2:.4f}  (within 0.05% of K)")
print(f"  72π²·√3          = {72*pi*pi*sqrt(3):.4f}  (same as above)")
print(f"  K / 12π          = {K/(12*pi):.4f}")
print(f"  K · A_0          = {K*A0:.4f}")
print(f"  K · A_0²         = {K*A0*A0:.4f}")
print(f"  K · 0.000242     = {K*0.000242:.4f}     (close to 0.30, an early anchor-z?)")
print(f"  1/K              = {1/K:.6e}")
print(f"  (c/H_0)/K        = {c_over_H0_Mpc/K:.4f} Mpc/SU  (the 'a' above)")
print()
print("→ K ≈ (1/A_0²) · (√3/2) is the cleanest closed-form match (0.05% off).")
print("   If K is exactly (12π)² · √3/2 = 72π²√3, that's a specific geometric")
print("   normalization (face-area × √3/2 = equilateral-triangle factor).")
print()

print("=" * 78)
print("BOTTOM LINE")
print("=" * 78)
print("- 'K·a = 1/H, K·q = 3/(20H)' is the entire content of the SU expansion;")
print("   the absolute value of K is set by the SU unit choice.")
print("- 'q/a = 3/20' is fixed regardless of K.")
print("- The literal identity '1 SU = A_0' holds iff one chooses")
print("       1 SU = bridge = c/(12π·H_0) ≈ 109.6 Mpc ≈ 357 Mly,")
print("   measuring all distances in Hubble-distance units.  In that convention")
print("   K·a = 12π and K = 12π·(SU/Mpc-of-the-spreadsheet) ratio.")
print("- K = 1231.35 is NOT 12π and NOT c/H_0 in any obvious distance unit;")
print("   it sits very near (1/A_0²)·√3/2 = 72π²√3, suggestive but not derived")
print("   from anything I can pin down without the STAM.xlsx derivation page.")
print("- Need to read STAM.xlsx (currently locked) to recover the exact K origin.")
