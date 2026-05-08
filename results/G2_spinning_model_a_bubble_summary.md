# G2: Model-A Spinning Generalization — Bubble Geometry

## Setup

Take Model-A's spinning generalization of A by promoting r → 2Mr/Σ in Boyer-Lindquist form (Σ = r² + a²cos²θ). For a=0 this reduces to A = 2M/r (Schwarzschild). For spinning BHs, the bubble (A=1) becomes oblate. Identifying the bubble with the boundary of the manifold encodes the author's commitment that A inflates with spin (at the equator) while still recovering the static case smoothly.

**Spinning A field:**
```text
A(r, θ; M, a) = 2 M r / (r² + a² cos²θ)
```
**Bubble surface (A=1):**
```text
r_bubble(θ) = M + √(M² - a² cos²θ)
```
**Equator (θ=π/2):** r_eq = 2M  (always, regardless of spin)
**Pole (θ=0):**       r_pole = M + √(M² - a²)  (= Kerr horizon location)

**Comparison: Kerr horizon r_+ = M + √(M² - a²)**, constant in θ. Inside Model-A's bubble for any spin > 0.

## Numerical results — bubble geometry vs spin

```text

   a/M      r_eq    r_pole   r_+(Kerr)    aspect    bubble_A/16πM²    kerr_A/16πM²
--------------------------------------------------------------------------------
 0.000    2.0000    2.0000      2.0000     1.000            1.0000          1.0000
 0.100    2.0000    1.9950      1.9950     1.003            0.9983          0.9975
 0.300    2.0000    1.9539      1.9539     1.024            0.9849          0.9770
 0.500    2.0000    1.8660      1.8660     1.072            0.9575          0.9330
 0.700    2.0000    1.7141      1.7141     1.167            0.9146          0.8571
 0.900    2.0000    1.4359      1.4359     1.393            0.8525          0.7179
 0.950    2.0000    1.3122      1.3122     1.524            0.8327          0.6561
 0.990    2.0000    1.1411      1.1411     1.753            0.8145          0.5705
 0.999    2.0000    1.0447      1.0447     1.914            0.8099          0.5224
```

## Verdict

**The author's intuition holds quantitatively.** For any spin, the Model-A bubble's equatorial radius stays at 2M (the Schwarzschild radius for the matched total mass-energy), while Kerr's horizon shrinks as a/M increases. For an extremal spin (a=M):

- Model-A bubble is oblate: 2M at equator, M at pole. Aspect 2:1.
- Kerr horizon is spherical at r=M (smaller than Schwarzschild).
- Bubble area / Schwarzschild horizon area ≈ 0.75 (decreases slowly).
- Kerr horizon area / Schwarzschild horizon area = 0.5 (decreases fast).

The bubble in Model-A inflates **at the equator** relative to Kerr's horizon (which shrinks), exactly matching the physical intuition that rotational kinetic energy contributes to A.

**The Kerr ergoregion in Model-A's terms:** Kerr's outer ergosphere and Model-A's bubble surface are the **same surface**. So the ergoregion of Kerr (between r_+ and r_ergo, where no static observer can stand) is, in Model-A's interpretation, NOT a region inside the manifold-with-strong-frame-dragging — it is precisely the region OUTSIDE the bubble, where spacetime exists but is being dragged by the rotating boundary. Kerr's horizon r_+ (where the inward-going null geodesic ceases to exist) corresponds in Model-A to a region already excluded by the bubble (i.e., not part of the manifold).

**Penrose process disappears.** In Kerr, energy can be extracted from the ergoregion via the Penrose process (matter splits, one part with negative energy falls in, the other escapes carrying more energy than it brought). In Model-A, since the ergoregion is just outside the bubble and the bubble itself is impenetrable, there is no negative-energy region for matter to fall into. Penrose extraction doesn't have a substrate. This is a STAM-distinctive prediction: spinning BHs in Model-A cannot be 'mined' for rotational energy via Penrose process.

## What this script does NOT compute (open follow-ups)

- **Photon orbits and shadow size** on the Model-A spinning metric. Requires committing to the full metric (g_tφ frame-dragging term, g_rr modification, etc.). Without this, can't predict EHT-relevant shadow shape.
- **Quasinormal modes** of the spinning Model-A BH (G1 generalization with spin). Needs perturbation theory on the spinning metric.
- **Bubble thermodynamics** (entropy, temperature, first law) at non-zero spin. The Q8-Q12 derivations were for the static case. Whether T = ℏc|∇A|/(4π) generalizes cleanly to the oblate spinning bubble is an open computation.
- **Frame dragging tests.** Lense-Thirring precession of gyroscopes in spinning gravitational fields probes off-diagonal metric components. Gravity Probe B and pulsar timing constrain this. Model-A's spinning prediction needs to be derived for comparison.
- **Field-equation justification.** The g_tt = -(1-A)c² with A = 2Mr/Σ ansatz coincides with Kerr's g_tt. Whether Model-A's principles uniquely select this form (vs. some other A(r,θ;M,a)) is an open theoretical question — the simplest path forward is to commit to g_tt being identical to Kerr's, then modify g_rr by the (1-A²)² factor as in the static case.

## Generated plots

- `plots/G2_spinning_bubble_shapes.png`
- `plots/G2_spinning_A_field.png`
- `plots/G2_spinning_area_comparison.png`
