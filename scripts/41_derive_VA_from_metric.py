"""
Script 41 - Derive V(A) from the bold STAM metric (geometric derivation)

Sean's commitment (2026-05-07):
  "The current V(A) looks clean but that doesn't mean it is correct."
  V(A) = beta/(1-A) is a guess. The honest test is whether the bold STAM
  metric ansatz, treated as a complete scalar-tensor theory, REQUIRES a
  specific V(A) for self-consistency.

Setup:
  Bold STAM static spherically symmetric metric (in c = 1 units):
      ds^2 = -(1 - A) dt^2 + dr^2 / [(1 - A)(1 - A^2)^2] + r^2 dOmega^2
  where A = A(r) is a scalar field.

  Treat this as Einstein gravity coupled to a canonical scalar:
      S = integral d^4x sqrt(-g) [R/(16 pi G) - (1/2) g^{mu nu} d_mu A d_nu A - V(A)]

  Then Einstein equations G_{mu nu} = 8 pi G T_{mu nu} with the scalar
  field stress-energy tensor must be self-consistent.

  T^t_t = -[ (1/2) g^rr (A')^2 + V(A) ]   (energy density-like)
  T^r_r =  [ (1/2) g^rr (A')^2 - V(A) ]   (radial pressure)
  T^theta_theta = T^phi_phi = -[ (1/2) g^rr (A')^2 + V(A) ]

  Compute G^mu_nu for the metric, set G^mu_nu = 8 pi G T^mu_nu, and see
  whether a consistent V(A) emerges.

  Three independent equations:
    (G^t_t equation) gives one relation between A, A', and V
    (G^r_r equation) gives another
    (G^theta_theta equation) gives a third

  Self-consistency: all three must be compatible. Subtracting two gives
  V(A) directly; subtracting the third gives a consistency check.

  If a consistent V(A) emerges, the bold STAM metric is a valid solution
  of scalar-tensor gravity with that V. If no consistent V exists, the
  metric ansatz needs modification (or non-canonical kinetic term, etc).
"""

from __future__ import annotations
import sympy as sp
from pathlib import Path

print("=" * 72)
print("Geometric derivation of V(A) from bold STAM metric")
print("=" * 72)

# Symbols
r = sp.Symbol('r', positive=True)
A_func = sp.Function('A')
A = A_func(r)
Ap = sp.diff(A, r)
App = sp.diff(A, r, r)

# Metric functions: f = -g_tt and h = g_rr
f = 1 - A
h = 1 / ((1 - A) * (1 - A**2)**2)

print(f"f(r) = -g_tt = {sp.simplify(f)}")
print(f"h(r) = g_rr = {sp.simplify(h)}")
print()
# Verify reduction to GR Schwarzschild when (1-A^2)^2 -> 1, i.e., when A -> 0
# (but Schwarzschild has A = Rs/r so A is not small near horizon)
# The key point: the EXTRA factor in g_rr vs Schwarzschild is (1-A^2)^(-2)

# ==================================================================
# Compute Christoffels, Ricci tensor, Einstein tensor symbolically
# Using formulas for a static spherically symmetric metric:
#   ds^2 = -f dt^2 + h dr^2 + r^2 dOmega^2
# the nonzero Einstein tensor components (mixed indices) are:
# G^t_t   = -1/r^2 + 1/(r^2 h) - h'/(r h^2)
# G^r_r   = -1/r^2 + 1/(r^2 h) + f'/(r f h)
# G^theta_theta = G^phi_phi
#               = (1/(2 h)) [ f''/f - (1/2)(f'/f)^2 - (1/2)(f' h')/(f h)
#                             + (1/r)(f'/f - h'/h) ]
# ==================================================================
fp = sp.diff(f, r)
fpp = sp.diff(f, r, r)
hp = sp.diff(h, r)

G_t_t_raw = -1/r**2 + 1/(r**2 * h) - hp/(r * h**2)
G_r_r_raw = -1/r**2 + 1/(r**2 * h) + fp/(r * f * h)
G_theta_raw = (1/(2*h))*(fpp/f - sp.Rational(1,2)*(fp/f)**2
                        - sp.Rational(1,2)*(fp*hp)/(f*h)
                        + (fp/f - hp/h)/r)

print("Computing Einstein tensor components (this may take a moment)...")
G_t_t = sp.simplify(G_t_t_raw)
G_r_r = sp.simplify(G_r_r_raw)
G_theta = sp.simplify(G_theta_raw)
print("Done.")
print()

# ==================================================================
# Scalar field stress-energy in mixed-index form for canonical scalar:
# T^t_t   = -[ (1/2) g^rr (A')^2 + V(A) ] = -[ (A')^2/(2h) + V ]
# T^r_r   =  [ (1/2) g^rr (A')^2 - V(A) ] = [ (A')^2/(2h) - V ]
# T^theta_theta = T^phi_phi = -[ (1/2) g^rr (A')^2 + V(A) ]
# ==================================================================
V = sp.Function('V')(A)

T_t_t = -(Ap**2 / (2*h) + V)
T_r_r =  (Ap**2 / (2*h) - V)
T_theta = -(Ap**2 / (2*h) + V)

# Einstein equations: G^mu_nu = kappa T^mu_nu, with kappa = 8 pi G (set to 1 in geom. units)
kappa = sp.Symbol('kappa', positive=True)

eq_t = sp.Eq(G_t_t, kappa * T_t_t)
eq_r = sp.Eq(G_r_r, kappa * T_r_r)
eq_theta = sp.Eq(G_theta, kappa * T_theta)

# ==================================================================
# Strategy: subtract eq_t from eq_r to eliminate V (since T^t_t and T^r_r
# differ in V sign). This gives a relation between A, A', A''.
# Then use eq_t alone to solve for V.
# ==================================================================
diff_r_t = sp.simplify(G_r_r - G_t_t - kappa * (T_r_r - T_t_t))
print("Constraint from G^r_r - G^t_t equation (must be zero for consistency):")
print(f"  G^r_r - G^t_t - kappa(T^r_r - T^t_t) = {diff_r_t}")
print()

# T^r_r - T^t_t = (A')^2/h
# So G^r_r - G^t_t = kappa * (A')^2/h
relation = sp.simplify(G_r_r - G_t_t)
print(f"G^r_r - G^t_t (geometric) = {relation}")
print()

# This must equal kappa * (A')^2 / h. So:
#   (A')^2 = (h/kappa) * (G^r_r - G^t_t)
# This DETERMINES A'(r) from the metric structure alone.
A_prime_sq_required = sp.simplify(h * relation / kappa)
print(f"(A')^2 required by EFE = {sp.simplify(A_prime_sq_required)}")
print()

# Now use G^t_t equation to solve for V:
#   G^t_t = -kappa [(A')^2/(2h) + V]
#   V = -G^t_t/kappa - (A')^2/(2h)
V_from_eq_t = sp.simplify(-G_t_t/kappa - Ap**2/(2*h))
print("V(r) from G^t_t equation (with A' replaced via consistency):")

# Substitute A_prime_sq in V using consistency: (A')^2 = h(G^r_r - G^t_t)/kappa
V_subbed = sp.simplify(V_from_eq_t.subs(Ap**2, A_prime_sq_required))
print(f"  V from t-eq, after substituting A'^2 from r-t consistency: ")
print(f"     V(r) = {V_subbed}")
print()

# Now check theta-theta equation as a consistency condition
# The theta equation should be automatically satisfied if t and r are.
# Compute residual:
V_in_theta = -G_theta/kappa - Ap**2/(2*h)
V_in_theta_subbed = sp.simplify(V_in_theta.subs(Ap**2, A_prime_sq_required))
residual = sp.simplify(V_in_theta_subbed - V_subbed)
print(f"Theta-equation consistency residual = {residual}")
print("(should be zero for the metric to be a valid scalar-tensor solution)")
print()

# ==================================================================
# Express V as a function of A (parametrically): we have V(r) and we
# need A(r) to invert. The (A')^2 relation gives A(r) implicitly.
# ==================================================================
print("="*72)
print("Summary of the geometric derivation")
print("="*72)

# Print the key relation: (A')^2 in terms of A and r
print("\n(A')^2 from EFE (G^r_r - G^t_t):")
expr = sp.together(A_prime_sq_required)
print(f"  (dA/dr)^2 = {expr}")

# And V from the t-equation
print("\nV expressed via metric (still has dA/dr in it):")
expr_V = sp.together(V_from_eq_t)
print(f"  V = {expr_V}")

# Also express V purely in terms of A and r after substituting (A')^2:
print("\nV after using EFE consistency to eliminate A':")
expr_V_clean = sp.together(V_subbed)
print(f"  V(r) = {expr_V_clean}")

# Save full report
outdir = Path('reports/script_41')
outdir.mkdir(parents=True, exist_ok=True)
with open(outdir/'derivation_report.txt', 'w', encoding='utf-8') as fh:
    fh.write("Geometric derivation of V(A) from bold STAM metric\n")
    fh.write("=" * 72 + "\n\n")
    fh.write(f"Metric:\n  f = -g_tt = {sp.pretty(f)}\n  h = g_rr = {sp.pretty(h)}\n\n")
    fh.write(f"G^t_t = {sp.pretty(G_t_t)}\n\n")
    fh.write(f"G^r_r = {sp.pretty(G_r_r)}\n\n")
    fh.write(f"G^theta_theta = {sp.pretty(G_theta)}\n\n")
    fh.write(f"(A')^2 required for consistency = {sp.pretty(A_prime_sq_required)}\n\n")
    fh.write(f"V from G^t_t (raw) = {sp.pretty(V_from_eq_t)}\n\n")
    fh.write(f"V from G^t_t (with A'^2 substituted) = {sp.pretty(V_subbed)}\n\n")
    fh.write(f"Theta-equation residual = {sp.pretty(residual)}\n")
    fh.write("(must equal zero for the metric to be a valid scalar-tensor solution)\n")

print(f"\nFull report saved to {outdir.resolve()}/derivation_report.txt")
print("="*72)
