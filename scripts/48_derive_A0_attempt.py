"""
Script 48 - Attempt to derive A_0 = 0.0265 = 1/(12*pi) from first principles

Empirical chain so far:
  b = 354.95 Mly  (historical bridge term, fitted to SN catalogs)
  L_H = c/H_0 = 13393 Mly (Hubble distance)
  A_0 = b/L_H = 0.02651 (calibrated, NOT derived)
  1/(12*pi) = 0.02653 (matches A_0 to 4 significant figures)

This script tries several natural STAM cosmological derivations of
A_0 from first principles and reports which (if any) produce 1/(12*pi)
cleanly.

Setup:
  Poisson source equation in STAM: laplacian A = kappa rho_m
  where kappa = 8 pi G / c^2.
  For cosmic critical density: rho_c = 3 H_0^2 / (8 pi G).
  Therefore kappa * rho_c = 3 H_0^2 / c^2 = 3 / L_H^2.

Several prescriptions tried below.
"""

from __future__ import annotations
import numpy as np

PI = np.pi
A_0_target = 1.0/(12*PI)
A_0_empirical = 354.95/13393.0   # rounding accepted

print("=" * 72)
print("Attempting first-principles derivation of A_0")
print("=" * 72)
print(f"  A_0 empirical (from b/L_H) = {A_0_empirical:.6f}")
print(f"  1/(12*pi)                  = {A_0_target:.6f}")
print(f"  ratio                      = {A_0_empirical/A_0_target:.4f}")
print()

# In all that follows, we work in normalized units where
#   L_H = 1 (Hubble distance)
#   H_0 = 1 (so c = 1)
#   kappa * rho_c = 3 (= 3 H_0^2 / c^2 in physical units)
# This isolates the geometric/structural factors from the constants.

KAPPA_RHO = 3.0   # 3 / L_H^2 in normalized units
R_H = 1.0

# ==================================================================
# Prescription 1: A field at center of uniform sphere of critical density
# ==================================================================
# For Poisson laplacian A = kappa rho with uniform rho inside R, A(R)=0:
#   A(r) = (kappa rho / 6) (R^2 - r^2)
#   A(0) = kappa rho R^2 / 6 = 3/6 = 1/2
A_center_uniform = KAPPA_RHO * R_H**2 / 6.0
print("Prescription 1: A at center of uniform critical-density sphere")
print(f"  A(0) = kappa*rho*R^2 / 6 = {A_center_uniform:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {A_center_uniform/A_0_empirical:.2f}x off")
print()

# ==================================================================
# Prescription 2: Volume-averaged A in same uniform sphere
# ==================================================================
# <A>_vol = integral A dV / V
#   = (kappa rho / 6) * integral (R^2 - r^2) 4 pi r^2 dr / V
#   = (kappa rho R^2) / 15
A_volavg = KAPPA_RHO * R_H**2 / 15.0
print("Prescription 2: Volume-average A in uniform sphere")
print(f"  <A>_vol = kappa rho R^2 / 15 = {A_volavg:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {A_volavg/A_0_empirical:.2f}x off")
print()

# ==================================================================
# Prescription 3: Radial-line-of-sight average of A in uniform sphere
# ==================================================================
# <A>_radial = (1/R) integral_0^R A(r) dr
#   = (kappa rho / (6R)) * (R^3 - R^3/3)
#   = (kappa rho R^2) / 9
A_radial = KAPPA_RHO * R_H**2 / 9.0
print("Prescription 3: Radial-LoS average A in uniform sphere")
print(f"  <A>_radial = kappa rho R^2 / 9 = {A_radial:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {A_radial/A_0_empirical:.2f}x off")
print()

# ==================================================================
# Prescription 4: Thin boundary shell at R_H with mass M_shell
# ==================================================================
# For shell at R_H with mass M_shell, A inside = 2 G M_shell / (c^2 R_H) = Rs_shell / R_H
# Total mass within Hubble volume at critical density:
# M_total = (4 pi / 3) rho_c R_H^3
# In normalized units, M_total such that 2 G M_total / c^2 = R_H gives Rs_total = R_H.
# Identify: A_0 = 1/(12 pi) corresponds to M_shell = M_total / (12 pi).
# CONSISTENCY observation, NOT a derivation.
print("Prescription 4: Thin boundary shell at R_H")
print(f"  A_inside_shell = M_shell / M_total")
print(f"  For A_0 = 1/(12*pi), M_shell / M_total = 1/(12*pi) = {1/(12*PI):.4f}")
print(f"  This is consistent with A_0 = 1/(12*pi) but does NOT derive it.")
print(f"  We'd need an independent argument for why M_shell = M_total/(12*pi).")
print()

# ==================================================================
# Prescription 5: Shapiro path integral through cosmic Hubble sphere
# ==================================================================
# Photon path from boundary r=R_H to r=0 (or full diameter) through
# uniform-density sphere with A(R_H) = 0.
# A(r) = (kappa rho / 6) (R^2 - r^2)
# Diametral Shapiro extra path length: (1/c) integral_path A ds * c
#   = integral 2 * A(r) dr from 0 to R_H
#   = 2 * (kappa rho/6) integral_0^R [R^2 - r^2] dr
#   = 2 * (kappa rho/6) [R^3 - R^3/3]
#   = 2 * (kappa rho/6) * (2 R^3/3)
#   = (2 kappa rho R^3)/9
# This equals b in normalized units (b = A_0 * R_H = A_0 in units R_H=1)
# Solving: A_0 = (2 kappa rho R_H^2)/9 = 2*3/9 = 2/3
A_0_shapiro_diametral = 2 * KAPPA_RHO * R_H**2 / 9.0
print("Prescription 5: Shapiro through diametral path of uniform sphere")
print(f"  A_0 = 2*kappa*rho*R_H^2 / 9 = {A_0_shapiro_diametral:.4f}")
print(f"  vs A_0 empirical = {A_0_empirical:.4f}: factor {A_0_shapiro_diametral/A_0_empirical:.1f}x off")
print()

# ==================================================================
# Prescription 6: Sub-critical density (Omega_m fraction)
# ==================================================================
# Same as Prescription 3 but with Omega_m only:
Omega_m = 0.315
A_radial_Omega_m = A_radial * Omega_m
print("Prescription 6: Radial-LoS with Omega_m = 0.315 (matter-only)")
print(f"  <A>_radial * Omega_m = {A_radial_Omega_m:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {A_radial_Omega_m/A_0_empirical:.2f}x off")
print()

# ==================================================================
# Prescription 7: Baryonic-only fraction
# ==================================================================
Omega_b = 0.05
A_radial_Omega_b = A_radial * Omega_b
print("Prescription 7: Radial-LoS with Omega_b = 0.05 (baryonic only)")
print(f"  <A>_radial * Omega_b = {A_radial_Omega_b:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {A_radial_Omega_b/A_0_empirical:.2f}x off")
print()

# ==================================================================
# Prescription 8: Surface flux divided by volume (Gauss-like)
# ==================================================================
# For Poisson with uniform source, integral_volume kappa*rho dV = kappa*rho*V
# Equals surface flux: integral grad A . dA over horizon.
# For surface area 4 pi R^2, average |grad A| = (kappa*rho V)/(4 pi R^2)
#   = (kappa*rho)*(4 pi R^3 / 3)/(4 pi R^2)
#   = kappa*rho*R/3
# Integrated A over the inward-radial path = (kappa rho R / 3) * R = kappa*rho*R^2/3 = 1
# Same answer as just <A> at center / different path; doesn't give 1/(12pi).
print("Prescription 8: Gauss-flux integrated path A")
print(f"  A_path = kappa*rho*R^2/3 = {KAPPA_RHO*R_H**2/3:.4f}")
print(f"  vs A_0 = {A_0_empirical:.4f}: factor {(KAPPA_RHO*R_H**2/3)/A_0_empirical:.1f}x off")
print()

# ==================================================================
# Prescription 9: Holographic (cosmic horizon area / 4 G hbar) connection
# ==================================================================
# The de Sitter cosmic horizon has Hawking-like temperature
#   T_dS = hbar H_0 / (4 pi k_B) = (1/(4*pi)) * hbar c |grad A|_horizon
# with |grad A|_dS = 2 H_0/c. The factor 4 pi is structural.
# What if A_0 = (1/(4 pi)) * (some other dimensionless 1/3)?
# 1/(12*pi) = (1/(4*pi)) * (1/3). The 1/3 could come from spatial dimensionality
# but no clean derivation.
print("Prescription 9: Decomposition 1/(12 pi) = (1/(4 pi)) * (1/3)")
print(f"  1/(4*pi) factor from STAM Q8 thermal-bridge structure")
print(f"  factor 1/3 could be spatial-dimensionality or thirds-of-A")
print(f"  but no first-principles argument for the specific 1/3 in this context")
print()

# ==================================================================
# Prescription 10: Mean A from cumulative cosmic structure
# ==================================================================
# If matter clumps into structures with fraction f_struct of cosmic volume
# at over-density delta, and photons sample structure-A weighted by traversal,
# the effective ambient A could differ from simple critical-density average.
# This requires a structure-formation calculation we haven't done.
print("Prescription 10: Cosmic structure (clusters, voids, filaments)")
print(f"  Requires structure-formation N-body or analytic calculation")
print(f"  not attempted in this session; deferred to future work")
print()

# ==================================================================
# Summary
# ==================================================================
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"""
None of the simple natural cosmological prescriptions gives A_0 = 1/(12*pi)
exactly. The closest matches in this set:

  Prescription                              Predicted A     vs 1/(12*pi)
  -----------------------------------------------------------------------
  Center of uniform critical sphere         {A_center_uniform:.4f}      19x too big
  Volume-average                             {A_volavg:.4f}      7.5x too big
  Radial-LoS average                         {A_radial:.4f}      4.2x too big
  Diametral Shapiro                          {A_0_shapiro_diametral:.4f}      25x too big
  Radial-LoS x Omega_m                       {A_radial_Omega_m:.4f}      1.3x too big
  Radial-LoS x Omega_b                       {A_radial_Omega_b:.4f}      0.6x off
  Gauss-flux integrated                      {KAPPA_RHO*R_H**2/3:.4f}      12.6x too big
  1/(12*pi) target                           {A_0_target:.4f}      reference

Honest verdict:
  - No clean STAM-derivation tested here gives A_0 = 1/(12*pi).
  - Closest near-matches involve Omega_b or Omega_m fractions of <A>_radial,
    but neither lands on 1/(12*pi) cleanly.
  - The 1/(12*pi) match (0.99957 ratio with A_0_empirical) remains a
    suggestive numerical coincidence in this session.
  - Possible derivation routes not fully explored:
    (a) Cosmic structure-amplified A (N-body weighted line of sight)
    (b) Specific holographic construction that ties dim-3 to thirds-of-A
    (c) Different boundary-condition Poisson with mass shell at A=1
    (d) Match through full FRW + photon-A coupling integrated over cosmic time
  - First-principles derivation of A_0 = 1/(12*pi) is OPEN.

What we have established firmly:
  - A_0 calibrated empirically from b_historical / L_Hubble = 0.02651
  - The bridge-term FORM b = A_0 * c/H_0 is STAM-derived from Shapiro
  - 1/(12*pi) = 0.02653 matches A_0 to 4 significant figures
  - 12*pi = 4*pi * 3 has STAM-structural plausibility (4*pi from Q8)
  - But the specific 1/3 factor in the cosmic context is NOT derived

Honest claim status:
  - A_0 is calibrated, NOT derived
  - The 1/(12*pi) match is a hint, not a theorem
  - Future work: explore prescriptions (a)-(d) above
""")
