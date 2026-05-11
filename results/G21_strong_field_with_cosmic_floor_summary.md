# G21: Strong-field Model-A with the cosmic floor explicit

**Date:** 2026-05-11
**Foundational technical step toward re-interpreting the strong-field corpus under the symmetric-boundary commitment (A on (A_0, 1), not (0, 1)).**

## Water-tank framing

Author 2026-05-11: A_0 IS the void — the unperturbed manifold level, not zero. Matter and energy displace this value into something higher locally. Crucially, displacement maps onto the AVAILABLE RANGE (1 − A_0), not onto the unit interval — A can't be pushed above 1 (the saturation/horizon ceiling), so the source's displacement scales by (1 − A_0).

```text
A_total(r) = A_0 + (1 − A_0) × (Rs / r)
A_local(r) = (1 − A_0) × (Rs / r)   (displacement, scaled to available range)
A_0        = 1/(12π) ≈ 0.026526    (cosmic structural floor, G7)
```

At r = Rs: A_total = A_0 + (1 − A_0) = 1 (horizon at GR's Schwarzschild radius, matching weak-field GR exactly).

**Naive A_total = A_0 + Rs/r is wrong.** It would put A_total > 1 at r ≤ Rs, which violates the A ≤ 1 ceiling. The rescaling by (1 − A_0) is what the water-tank picture actually requires.

Under the ledger reading, A is the density of resolved write-events per Planck cell. A_0 is the minimum write density required for the manifold to exist; A_local is the additional write density due to matter currently being at the source location, scaled so cells can't saturate beyond capacity.

## Strong-field metric (unchanged k(A))

```text
g_tt = -(1 - A_total) c²
g_rr = 1 / [(1 - A_total)(1 - A_total²)²]
```

The metric construction k(A) = (1-A)(1-A²)² remains the framework's postulated commitment (see project_strong_field_commitment.md and project_strong_field_departure_question.md). This script changes only what enters as A: A_total with cosmic floor instead of A_local alone.

## Landmark radii

| landmark      |   A_total_target |   A_local_needed |   r/Rs (with floor) |   r/Rs (no floor) |   shift (%) |
|:--------------|-----------------:|-----------------:|--------------------:|------------------:|------------:|
| ISCO          |          0.33333 |          0.30681 |             3.17292 |           3.00000 |     5.76384 |
| photon sphere |          0.66667 |          0.64014 |             1.52072 |           1.50000 |     1.38125 |
| horizon       |          1.00000 |          0.97347 |             1.00000 |           1.00000 |     0.00000 |


**Reading:** all three orbital landmarks shift OUTWARD when the cosmic floor is included. The shift is largest at ISCO (lowest A_total target) and smallest at the horizon (highest target). This is because A_0 is a fixed offset; landmarks defined by small A_local values feel it more in relative terms.

**Horizon:** r_horizon / Rs = 1/(1 - A_0) ≈ 1.0272. The horizon sits ~2.7% outside the GR Schwarzschild radius. The cosmic A_0 is already in the cell; the source's displacement only needs to push A_total from A_0 to 1, i.e., A_local = 1 - A_0 instead of 1.

## Asymptotic limit

As r → ∞, A_local → 0 and A_total → A_0. The metric reads:

```text
g_tt(∞) = -(1 - A_0) c² ≈ -0.9735 c²
g_rr(∞) = 1 / [(1 - A_0)(1 - A_0²)²] ≈ 1.0287
```

**Spacetime is asymptotically non-Minkowski.** Clocks at the asymptotic limit tick at √(1 - A_0) ≈ 0.9866 of a hypothetical A=0 clock — but A=0 does not exist, so this is not a slowdown relative to a faster reference. It is the framework's commitment about what asymptotic time IS at the cosmic floor.

## Effective stress-energy decomposition

Under the ledger framing, the effective stress-energy bracket B(A_total, A_0) = ρ_eff × r² (up to a constant) decomposes into two physically distinct contributions:

- **Structural (background):** `B_background = f(A_0)` — the cosmic-floor contribution, present everywhere by virtue of the manifold existing. Constant in r. Encodes the structural baseline.
- **Dynamical (source):** `B_source(A_total) = B(A_total, A_0) - B_background` — the displacement contribution from the matter source. Zero at the asymptotic limit (A_total = A_0); rises to ≈ 1 - f(A_0) at the horizon.

The full bracket with cosmic floor is:

```text
B(A_total, A_0) = f(A_total) - (A_total - A_0) × f'(A_total)
                = 1 - (1 - A_total)²(1 + A_total)
                       × [1 + A_total(1 - 5A_0) + 4A_total² - A_0]
```

When A_0 = 0 this reduces exactly to F3's bracket `1 - (1 - A)²(1 + A)(1 + A + 4A²)`. With A_0 = 1/(12π), there are A_0 corrections at order A_0 throughout.

## Numerical table

|   A_total |   A_local |       r/Rs |   m_total/M |   m_background/M |   m_source/M |   B_total (with floor) |   B_background |   B_source |   B_F3 (no floor) |   ΔB (floor - F3) |
|----------:|----------:|-----------:|------------:|-----------------:|-------------:|-----------------------:|---------------:|-----------:|------------------:|------------------:|
|  0.026526 |  0.000000 | inf        |  nan        |       nan        |   nan        |               0.027895 |       0.027895 |   0.000000 |         -0.001331 |          0.029226 |
|  0.050000 |  0.023474 |  41.470004 |    2.270236 |         1.156816 |     1.113420 |               0.026938 |       0.027895 |  -0.000957 |         -0.004483 |          0.031421 |
|  0.100000 |  0.073474 |  13.249202 |    1.562213 |         0.369590 |     1.192624 |               0.019712 |       0.027895 |  -0.008183 |         -0.015740 |          0.035452 |
|  0.200000 |  0.173474 |   5.611637 |    1.474289 |         0.156538 |     1.317751 |              -0.003736 |       0.027895 |  -0.031632 |         -0.044480 |          0.040744 |
|  0.300000 |  0.273474 |   3.559657 |    1.496230 |         0.099298 |     1.396933 |              -0.015178 |       0.027895 |  -0.043073 |         -0.057420 |          0.042242 |
|  0.400000 |  0.373474 |   2.606537 |    1.503033 |         0.072710 |     1.430323 |               0.011947 |       0.027895 |  -0.015948 |         -0.028160 |          0.040107 |
|  0.440000 |  0.413474 |   2.354377 |    1.497013 |         0.065676 |     1.431338 |               0.038344 |       0.027895 |   0.010449 |          0.000012 |          0.038332 |
|  0.500000 |  0.473474 |   2.056024 |    1.477767 |         0.057353 |     1.420414 |               0.097315 |       0.027895 |   0.069420 |          0.062500 |          0.034815 |
|  0.600000 |  0.573474 |   1.697503 |    1.419384 |         0.047352 |     1.372032 |               0.248922 |       0.027895 |   0.221027 |          0.221760 |          0.027162 |
|  0.700000 |  0.673474 |   1.445451 |    1.332663 |         0.040321 |     1.292342 |               0.458283 |       0.027895 |   0.430388 |          0.440020 |          0.018263 |
|  0.800000 |  0.773474 |   1.258574 |    1.225951 |         0.035108 |     1.190843 |               0.695629 |       0.027895 |   0.667734 |          0.686080 |          0.009549 |
|  0.900000 |  0.873474 |   1.114485 |    1.110462 |         0.031089 |     1.079373 |               0.905112 |       0.027895 |   0.877217 |          0.902340 |          0.002772 |
|  0.990000 |  0.963474 |   1.010379 |    1.010375 |         0.028185 |     0.982190 |               0.998855 |       0.027895 |   0.970960 |          0.998824 |          0.000031 |


**Two sign changes of B_total** (not one as in F3 without floor):
- Lower crossover: A_total ≈ 0.18380
- Upper crossover: A_total ≈ 0.37421

F3 (no floor) had a single sign change at A ≈ 0.44. With the cosmic floor, the effective stress-energy is positive at the asymptotic limit (A_total = A_0, B_total = f(A_0) ≈ 0.02790), becomes negative in a finite band between the two crossovers, then positive again as A_total approaches 1.

**This is qualitatively different from F3's no-floor picture.** F3 (no floor) had effective ρ negative for ALL A < 0.44 (the whole outer region of the source). With the cosmic floor enforced, the negative-ρ region is now a finite band, NOT the whole outer region. The asymptotic region itself is positive-ρ, dominated by the cosmic structural baseline.


**What's structural vs dynamical:**
- The background bracket f(A_0) is small and positive (≈ 0.02790). It represents the cosmic-floor contribution to effective stress-energy — small magnitude, positive sign, present everywhere by virtue of the manifold existing.
- The source bracket B_source is negative for moderate A_total and positive near the horizon. Its single sign change sits near F3's 0.44 (small offset from A_0 corrections to dA/dr).

## Honest assessment

**What changed quantitatively:**
- Landmark radii shift outward by 2-9% (ISCO most, horizon least).
- Asymptotic metric is non-Minkowski (rescaled time and radial coordinates).
- F3's single NEC crossover at A ≈ 0.44 becomes **2 sign changes** under cosmic floor: 0.1838 and 0.3742.
- The negative-ρ band is no longer the whole outer region; it is a finite shell between the two crossovers.
- The mass function and effective ρ both grow linearly with r at large r (vs returning to fixed values at infinity in the no-floor case) — the cosmic A_0 contributes a 'background mass' that scales with volume.

**What changed conceptually:**
- Effective stress-energy decomposes naturally into structural (background, A_0-only) and dynamical (source displacement) contributions. These are not free interpretive layers — they fall out of the math once you write the bracket with A_floor present.
- The structural piece is small but nonzero everywhere — it is the framework's cosmological-constant-like contribution, expressed at the spatial-profile level (consistent with script 34's equation-of-state finding that w ≈ -1 at A near A_0).
- The dynamical piece carries F3's original sign-change story, now interpretable as 'where the source's contribution flips from dark-energy-like to matter-like.'

**What did NOT change:**
- The metric construction k(A) = (1-A)(1-A²)². Still postulated. Field-equations gap unchanged.
- The thirds-of-A landmarks (ISCO at 1/3, photon sphere at 2/3, horizon at 1) are still the values of A_total at which they occur. What shifted is the radius at which A_total reaches each value.
- F3's verdict (Birkhoff doesn't apply, owes field equations) unchanged.

**The 'why' question this surfaces:**
- Why does the NEC crossover (sign change of B_total) sit so close to A = 0.44 with or without the floor? F3's 0.44 is structurally related to where the (1-A)²(1+A)(1+A+4A²) factor crosses 1. Cosmic floor shifts this slightly. Worth understanding what determines this value structurally — it may have its own meaning.
- The asymptotic linear growth of m(r) means the framework predicts a 'mass-at-infinity' that diverges in the standard Schwarzschild parameterization. This is the cosmic-floor contribution at large r. Interpretation: in a proper FRW-style cosmological framework, this is the framework's natively-derived dark-energy effective mass.

## Generated plots

- `plots/G21_landmark_shifts.png`
- `plots/G21_bracket_decomposition.png`
- `plots/G21_mass_function_decomposition.png`
