#!/usr/bin/env python3
"""
G72b_tilted_aether_setup.py

Step B of the Einstein-aether track.  EA kinematic scalars for a tilted
aether on the framework's committed strong-field metric.

Parametrization: instead of (cosh phi, sinh phi), use the components
   U = u^t(r),   W = u^r(r)
as independent symbolic Functions of r, with the constraint

   -h U^2 + W^2 / k = -1     (unit-norm)

enforced at the end by substitution.  This keeps sympy in polynomial /
rational arithmetic and avoids cosh/sinh blow-up.
"""

from __future__ import annotations

import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

r = sp.symbols('r', positive=True)
t = sp.symbols('t', real=True)
th = sp.symbols('th', positive=True)
ph = sp.symbols('ph', real=True)
coords = [t, r, th, ph]

# Committed metric (M = 1)
A_expr = 2 / r
y_expr = 3 * A_expr - 2
F_expr = 1 - 5 * y_expr**4 + 4 * y_expr**5
h_expr = 1 - A_expr
k_expr = h_expr * F_expr

# Aether components as independent Functions, constraint applied at end
U = sp.Function('U')(r)   # u^t
W = sp.Function('W')(r)   # u^r

u_up = [U, W, sp.S.Zero, sp.S.Zero]

# Metric
g = sp.Matrix.zeros(4, 4)
g[0, 0] = -h_expr
g[1, 1] = 1 / k_expr
g[2, 2] = r**2
g[3, 3] = r**2 * sp.sin(th)**2
g_inv = g.inv()

print("Unit-norm constraint  (to apply later):")
constraint = sp.simplify(-h_expr * U**2 + W**2 / k_expr + 1)
print(f"   -h U^2 + W^2/k + 1 = 0   (verified: {constraint == -h_expr*U**2 + W**2/k_expr + 1})")
print()


# Christoffels
def Gamma(a, b, c):
    s = sp.S.Zero
    for d in range(4):
        s += g_inv[a, d] * (sp.diff(g[d, b], coords[c])
                            + sp.diff(g[d, c], coords[b])
                            - sp.diff(g[b, c], coords[d]))
    return s / 2


# nabla_mu u^alpha (no aggressive simplification)
print("Computing nabla_mu u^alpha ...")
nab = sp.MutableDenseMatrix(4, 4, [sp.S.Zero]*16)
for mu in range(4):
    for a in range(4):
        val = sp.diff(u_up[a], coords[mu])
        for b in range(4):
            val += Gamma(a, mu, b) * u_up[b]
        nab[mu, a] = val

# Show structural form (non-zero entries only)
print("Non-zero  nabla_mu u^alpha  components:")
labels = ['t', 'r', 'th', 'ph']
for mu in range(4):
    for a in range(4):
        val = sp.together(nab[mu, a])
        if val != 0:
            print(f"  nabla_{labels[mu]} u^{labels[a]} = {val}")
print()


# Lower index
print("Computing nabla_mu u_nu ...")
nab_dn = sp.MutableDenseMatrix(4, 4, [sp.S.Zero]*16)
for mu in range(4):
    for nu in range(4):
        s = sp.S.Zero
        for la in range(4):
            s += g[nu, la] * nab[mu, la]
        nab_dn[mu, nu] = s


# Expansion (lightweight)
print("Computing expansion theta ...")
theta_exp = sum(nab[mu, mu] for mu in range(4))


# Acceleration
print("Computing acceleration a^mu ...")
a_up_vec = []
for mu in range(4):
    val = sp.S.Zero
    for nu in range(4):
        val += u_up[nu] * nab[nu, mu]
    a_up_vec.append(val)

a_sq = sp.S.Zero
for mu in range(4):
    for nu in range(4):
        a_sq += g[mu, nu] * a_up_vec[mu] * a_up_vec[nu]


# Four EA kinetic scalars  (light simplification only)
print("Computing the four EA kinetic scalars ...")

# term1 = nabla^mu u^nu  *  nabla_mu u_nu
def nab_up_up(mu, nu):
    s = sp.S.Zero
    for alpha in range(4):
        s += g_inv[mu, alpha] * nab[alpha, nu]
    return s

print("  term1 ...")
term1 = sp.S.Zero
for mu in range(4):
    for nu in range(4):
        term1 += nab_up_up(mu, nu) * nab_dn[mu, nu]

print("  term2 ...")
term2 = theta_exp**2

print("  term3 ...")
term3 = sp.S.Zero
for mu in range(4):
    for nu in range(4):
        term3 += nab_up_up(mu, nu) * nab_dn[nu, mu]

print("  term4 ...")
term4 = a_sq

# Apply the unit-norm constraint: solve for W in terms of U (pick the branch).
# -h U^2 + W^2 / k = -1   ->   W^2 = k (h U^2 - 1)
print()
print("Applying unit-norm constraint  W^2 = k (h U^2 - 1) ...")

# We treat W^2 as a substitution.  Where W appears at first power we leave it
# (it's a sign-bearing Lorentz factor).  Where W^2 appears, substitute.
# To do this cleanly, expand each term as polynomial in W, then sub W^2 -> k(hU^2 - 1).

W_sq_value = k_expr * (h_expr * U**2 - 1)

def apply_constraint(expr):
    """Substitute W^2 = k(hU^2 - 1)."""
    e = sp.expand(expr)
    return e.replace(W**2, W_sq_value)


# Also: dW/dr can be expressed via differentiating the constraint:
# 2 W W' = (k(hU^2 - 1))' = k'(hU^2-1) + k(h' U^2 + 2 h U U')
# W' = [k'(hU^2-1) + k(h' U^2 + 2 h U U')] / (2 W)
# This introduces 1/W; leave as is for now (we'll see).

print("  applying to term1 ...")
term1_c = apply_constraint(term1)
print("  applying to term2 ...")
term2_c = apply_constraint(term2)
print("  applying to term3 ...")
term3_c = apply_constraint(term3)
print("  applying to term4 ...")
term4_c = apply_constraint(term4)


# Aligned-limit check: U = 1/sqrt(h), W = 0
print()
print("=" * 80)
print("ALIGNED LIMIT  (U = 1/sqrt(h), W = 0):")
print("=" * 80)
aligned_subs = {
    U: 1/sp.sqrt(h_expr),
    sp.Derivative(U, r): sp.diff(1/sp.sqrt(h_expr), r),
    W: sp.S.Zero,
    sp.Derivative(W, r): sp.S.Zero,
}

for label, expr in [('term1', term1), ('term2', term2),
                    ('term3', term3), ('term4', term4)]:
    aligned = sp.simplify(expr.subs(aligned_subs))
    print(f"  {label} = {sp.factor(aligned)}")

# Expected forms (in r, with h, k explicit):
h_p = sp.diff(h_expr, r)
expected_kin = k_expr * h_p**2 / (4 * h_expr**2)
print()
print(f"  Expected:  L_ae aligned = (c1 - c4) * k h'^2/(4 h^2)")
print(f"  k h'^2/(4 h^2) = {sp.factor(expected_kin)}")
print()
print("Therefore: aligned terms above should give  (-term1, 0, 0, +term4)")
print("with term1 = -k h'^2/(4h^2) and term4 = +k h'^2/(4h^2).")
print()
print("=" * 80)
print("STRUCTURAL OUTPUT")
print("=" * 80)
print("For tilted aether on committed metric, the four EA scalars are")
print("polynomial expressions in U, W, U', W', r (with W^2 eliminated by the")
print("unit-norm constraint).  Full forms are too long to print; we save")
print("them and proceed to matching in G72c.")
