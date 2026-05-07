# STAM Model-B — Formal Specification

**Author**: Sean Brady

**Status**: refinement of Model-A; introduces ONE category distinction (GWs are space, light is matter on space) that resolves F5 and dissolves several other ambiguities.



Model-B *refines* Model-A; it does not replace it. All Model-A weak-field, thermodynamic, and strong-field derivations carry over unchanged. Model-B specifies the gravitational-wave sector and the cosmological A field structure more cleanly than Model-A did, and in doing so resolves the F5 fatality.


---

## Core commitments

```text

STAM MODEL-B — CORE COMMITMENTS
================================

1. TWO-FIELD STRUCTURE.
   Model-B has two fundamental fields:
       (a) The metric g_μν: a (0,2) tensor. The geometric structure of space.
       (b) The accumulation field A: a scalar. A feature of space, sourced
           by matter.
   Both are independent dynamical objects, coupled through the construction
   rule (item 3) and the source equation (item 4).

2. CATEGORY DISTINCTION (the new commitment).
   GWs ARE space — propagating perturbations of g_μν itself. Tensor character
   gives h_+, h_× polarizations natively. Speed = c intrinsically.
   Light IS matter on space — photons propagating along null geodesics of
   the metric. Slowed by A through the standard Shapiro/SU mechanism.
   These are categorically different; STAM does not predict them to behave
   the same way.

3. METRIC CONSTRUCTION RULE (carried from Model-A Framework C).
   For static spherically symmetric configurations, the metric components
   are determined by the A field via:
       g_tt = -(1-A) c²
       g_rr = 1 / [(1-A)(1-A²)²]
       g_θθ = r²
       g_φφ = r² sin²θ
   For non-static / non-spherical configurations, the metric has its own
   degrees of freedom (giving GWs as tensor perturbations). The construction
   rule applies to the "static averaged" structure; perturbations on top
   of it carry the tensor GW modes.

4. SOURCE EQUATION FOR A (carried from Model-A Framework C).
   The scalar A field is sourced by matter density:
       ∇²A = (8π G / c²) ρ_matter   (static / weak-field)
       □A  = (8π G / c²) ρ_source    (relativistic generalization)
   In cosmic voids where ρ_matter ≈ 0, A satisfies Laplace's equation, so
   A ≈ 0 (the only bounded solution). NO smooth cosmological background;
   the de Sitter A_cosmo(r) ansatz used in Model-A Q8 was incorrect for
   the actual cosmology.

5. BUBBLE PICTURE FOR HORIZONS (carried from Model-A).
   A=1 surfaces are the 2D phase boundaries between spacetime (A<1) and
   not-spacetime (A>1 doesn't exist as a manifold region). Black holes
   are bubbles. Matter accumulates asymptotically on the boundary.
   No interior, no singularity, no information paradox.

6. RESOLUTION RULE (carried from Model-A).
   Beyond A=1, the manifold ends. Inward-going paths from A=1 do not
   exist. Vacuum fluctuations at the boundary can resolve only outward,
   producing thermal emission (Hawking radiation) at temperature
   k_B T = (1/(4π)) ℏ c |∇A| via phase-boundary equilibrium.

7. NO-CROSSING (carried from Model-A bold extension).
   For an infaller approaching A=1, proper time to reach the boundary
   diverges logarithmically. The traveler "falls forever" in their own
   frame; the observer sees asymptotic freezing. Both frames register
   no-crossing event in the universal ledger.

8. THIRDS-OF-A (carried from Model-A).
   ISCO at A=1/3, photon sphere at A=2/3, horizon at A=1.
   Preserved exactly because g_tt is identical to GR.

```

## What carries over from Model-A

```text

ALL MODEL-A RESULTS CARRIED INTO MODEL-B (UNCHANGED):

Local weak-field (5 derivations):
    1. Horizon threshold A = 1 ↔ r = Rs.
    2. Newtonian gravity from gravity bridge g = (c²/2) ∇A.
    3. GPS clock correction (+38.57 μs/day satellite gain).
    4. Shapiro propagation delay (123.6 μs Earth-Mars solar grazing).
    5. Mass-estimator consistency.

Black hole thermodynamics (6 derivations, scripts Q8-Q12):
    6. Hawking T (Schwarzschild + Unruh + de Sitter from one rule).
    7. Planckian spectrum from phase-boundary equilibrium.
    8. Bekenstein entropy S = k_B A / (4 ell_P²).
    9. First law dE = T dS (machine precision).
    10. Smarr formula M c² = 2 T S (machine precision).
    11. Hawking evaporation lifetime (matches textbook).
    12. Generalized second law (+1/3 surplus).

Strong-field (4 derivations, scripts SF3-SF6):
    13. Bold-STAM g_rr metric ansatz (preserves thirds-of-A).
    14. No-crossing infall with logarithmic proper-time divergence.
    15. GR exterior recovery (passes solar system / pulsar tests by 10⁴ - 10¹⁴).
    16. Multi-source linear superposition; binary merger d_crit = 4 Rs exact.

Fatality survivals (Model-A status, all carried):
    F1 (Penrose-Hawking): manifold class evades premise.
    F2 (Equivalence principle): collapses into F1.
    F3 (Birkhoff): bold STAM is not Einstein vacuum gravity.

NONE of the above depends on GW polarization structure or cosmological-A
profile, so all are preserved unchanged in Model-B.

```

## What Model-B refines

```text

WHAT MODEL-B REFINES vs MODEL-A:

(R1) F5 RESOLVED.
    Model-A Framework C predicted GWs as scalar A perturbations → falsified
    by LIGO (no transverse-tensor modes).
    Model-B treats GWs as propagating perturbations of the metric tensor
    g_μν itself, with intrinsic tensor character. h_+, h_× polarizations
    are produced natively. F5 is no longer a falsification.

(R2) GW170817 CONSISTENCY MADE EXPLICIT.
    GWs travel at c intrinsically (they ARE space).
    Light travels at c through cosmic voids (where A ≈ 0 by Framework C
    source equation), slowed only near matter where A > 0.
    Both arrive at Earth from cosmic-distance sources within ~astrophysical
    emission timing differences. Consistent with the |Δv|/c < 10⁻¹⁵
    constraint from GW170817.

(R3) DE SITTER A_cosmo(r) = (r/R_dS)² ANSATZ ABANDONED.
    Model-A Q8 used this ansatz to derive the de Sitter horizon temperature
    alongside Schwarzschild and Unruh. Model-B recognizes this ansatz was
    not derived from Framework C — it was a postulate. The actual STAM
    cosmology (Framework C source equation with no cosmological background)
    has A ≈ 0 in cosmic voids.

    Consequences:
        - Q8's Schwarzschild and Unruh derivations are unaffected (they're
          local, depending on |grad A| at the boundary, not on cosmological A).
        - Q8's de Sitter case becomes more speculative; whether STAM
          predicts the de Sitter horizon temperature in this cleaner
          cosmology is an open computation.
        - Script 31 (cosmological distance) was based on the wrong ansatz.
          Bold-STAM cosmological distance prediction needs redoing under
          Framework C/D's actual matter-sourced A.

(R4) CRITICAL FOLLOW-UP RECONTEXTUALIZED.
    The F3-extended cosmology calculation flagged earlier (potential dark-
    energy-like behavior from bold-STAM effective stress-energy) needs to
    be redone with Model-B's actual A profile (matter-sourced), not the
    de Sitter ansatz. Whether STAM derives the bridge term b is still an
    open computation, but now it's a properly specified one.

```

## What Model-B leaves open

```text

WHAT MODEL-B LEAVES OPEN:

(O1) Metric field equations.
    Model-B has the metric construction rule for static spherical
    configurations and treats the metric as having tensor degrees of
    freedom for GW perturbations. A unified field equation for the
    metric (analog of Einstein equations or scalar-tensor field
    equations) that recovers BOTH the construction rule AND the GW
    propagation, ideally from a single Lagrangian, is open work.

(O2) Cosmological distance / bridge term b.
    The matter-sourced A field (Framework C) over cosmological scales
    needs to be computed properly. With the correct A profile, the
    distance-redshift relation can be derived. Whether it predicts
    catalog observations directly or with a derivable bridge term is
    open.

(O3) Galaxy rotation curves.
    Marked exploratory in Model-A. Same status in Model-B; the matter-
    sourced A field around a galaxy might give different rotation
    curves than Newton or LCDM-with-DM. Open computation.

(O4) Quasinormal mode frequencies for binary BH ringdown.
    Bold-STAM g_rr modification near horizon could shift QNM frequencies
    by a few percent. Specific predictions for LIGO ringdown comparisons
    are open.

(O5) Dynamical-collapse rigorous proof.
    F1 was clean for the static manifold. Whether gravitational collapse
    smoothly produces a Model-B-class manifold without ever transiently
    forming a strict trapped surface is open dynamical work.

```

## Internal consistency verifications

- **Field structure: g_uv (tensor) and A (scalar) as independent fields**: PASS
  - Both fields specified. Coupling via construction rule + source equation.
- **Model-B predicts |v_GW - v_EM|/c << 1e-15 for cosmic-void propagation**: PASS
  - Both messengers at c through cosmic voids. Multi-messenger consistent.
  - Δv/c = 0.000e+00, constraint = 1.000e-15
- **Model-B GWs: tensor polarizations h_+, h_x native to metric structure**: PASS
  - GWs are perturbations of (0,2) tensor g_uv, hence tensor character. F5 falsification dissolved.
- **All Q8-Q12 thermodynamics and SF3-SF6 strong-field results carry over**: PASS
  - Local results depend on g_tt (unchanged) and on |grad A| (unchanged at boundaries). Model-B preserves all 16 derived results from Model-A.
- **Model-B predicts A ~= 0 in cosmic voids (Laplace + bounded BC)**: PASS
  - Replaces the Model-A de Sitter ansatz A_cosmo = (r/R_dS)^2. Voids: rho_matter = 0 -> grad^2 A = 0 -> A = 0.


## Generated plots

- `plots/STAM_Model_B_vs_Model_A.png`
- `plots/STAM_Model_B_field_structure.png`


---

## Bottom line

STAM Model-B is the cleanest single statement of what bold STAM actually IS. It has two independent fields (g_μν tensor and A scalar), one source equation for A from matter density, one construction rule for the static metric, and a clear category distinction between GWs (which ARE space, with intrinsic tensor character) and light (which is matter on space, slowed by A).

Every result derived in Model-A carries over. F5 (the only fatality Model-A failed) is resolved by the category distinction — GWs naturally have tensor character because the metric naturally has tensor character. Multi-messenger consistency (GW170817) is automatic in cosmic voids where the matter-sourced A field is essentially zero.

Model-B is *not* a complete theory — it leaves open the unified field equations for the metric, the cosmological-distance derivation, and the galaxy-rotation problem. But it is a substantially stronger framework than Model-A, free of the F5 falsification and free of the de Sitter ansatz that didn't follow from Framework C.
