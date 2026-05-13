#!/usr/bin/env python3
"""
G73_box_A_and_routes.py

Computes the covariant d'Alembertian box A on the committed metric --
the building block for every higher-derivative route (DHOST, multi-field
kinetic mixing, or the substance-ontology Gamma_res formulation).

For static spherically symmetric ds^2 = -h dt^2 + dr^2/k + r^2 dOmega^2
with A = A(r):

    box A = (1 / sqrt(-g)) d_r ( sqrt(-g) g^{rr} d_r A )
          = (1 / (r^2 sqrt(h/k))) d_r ( r^2 sqrt(h/k) * k * A'(r) )

For A = 2/r (M = 1, spinless), and the committed k(A) = (1-A) F(y),
this is a specific closed-form rational function of r.

Also computes:
  - (nabla A)^2 = g^{mu nu} d_mu A d_nu A
  - (box A)^2  (the leading DHOST building block)
  - The mass scale set by box A near PS and near horizon
"""

import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

r = sp.symbols('r', positive=True)
A = 2 / r
y = 3 * A - 2
F = 1 - 5 * y**4 + 4 * y**5
h = 1 - A
k = h * F

Ap = sp.diff(A, r)
sqrt_minus_g_radial = r**2 * sp.sqrt(h / k)  # angular part 4 pi sin(theta) absorbed

inner = sqrt_minus_g_radial * k * Ap
box_A = sp.diff(inner, r) / sqrt_minus_g_radial
box_A = sp.simplify(box_A)

grad_A_sq = k * Ap**2
grad_A_sq = sp.simplify(grad_A_sq)

print("=" * 80)
print("Covariant d'Alembertian box A on committed metric (M = 1)")
print("=" * 80)
print()
print("(grad A)^2 = g^{mu nu} d_mu A d_nu A:")
print(f"  {sp.factor(grad_A_sq)}")
print()
print("box A:")
print(f"  {sp.factor(box_A)}")
print()
print("(box A)^2:")
box_A_sq = sp.expand(box_A**2)
print(f"  {sp.factor(box_A_sq)}")
print()

# Sample at key points
print("=" * 80)
print("Sample values across the final shell:")
print("=" * 80)
print(f"{'A':>6}{'r':>8}{'(grad A)^2':>16}{'box A':>16}{'(box A)^2':>16}")
print("-" * 62)
for A_val in [0.67, 0.70, 0.80, 0.90, 0.953, 0.99]:
    r_val = 2 / A_val
    gA2 = float(grad_A_sq.subs(r, r_val))
    bA = float(box_A.subs(r, r_val))
    bA2 = float(box_A_sq.subs(r, r_val))
    print(f"{A_val:>6.3f}{r_val:>8.3f}{gA2:>16.6e}{bA:>16.6e}{bA2:>16.6e}")
print()

# Boundary behavior
print("=" * 80)
print("Boundary behavior:")
print("=" * 80)
print()
print("Near PS (A -> 2/3, r -> 3):")
box_A_at_PS = sp.limit(box_A, r, 3, '-')
print(f"  box A | A=2/3 = {sp.simplify(box_A_at_PS)}")
print()
print("Near horizon (A -> 1, r -> 2):")
# Series expansion around r = 2
box_A_series = sp.series(box_A, r, 2, 4).removeO()
print(f"  box A near r=2:  {sp.simplify(box_A_series)}")
print()

# Outside PS (A < 2/3), the metric is Schwarzschild (F = 1, k = h = 1-A).
# Compute box A in that limit to confirm continuity.
print("=" * 80)
print("Outside PS sanity check (k = h = 1 - 2/r, Schwarzschild vacuum):")
print("=" * 80)
A_sw = 2 / r
h_sw = 1 - A_sw
k_sw = 1 - A_sw  # F = 1 outside
grad_A_sq_sw = k_sw * sp.diff(A_sw, r)**2
inner_sw = r**2 * sp.sqrt(h_sw / k_sw) * k_sw * sp.diff(A_sw, r)
box_A_sw = sp.simplify(sp.diff(inner_sw, r) / (r**2 * sp.sqrt(h_sw / k_sw)))
print(f"  box A (Schwarzschild) = {sp.factor(box_A_sw)}")
print()
print("For Schwarzschild, A = 2M/r is harmonic with respect to the spatial part:")
print("  Indeed box A involves only k h' / 2 + ... pieces that simplify cleanly.")
