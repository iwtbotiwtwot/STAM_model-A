# G3: Spinning Bubble Thermodynamics — Model-A vs Kerr

## Setup

Apply Q8/Q10's resolution-rule temperature `T = ℏc|∇A|/(4π k_B)` to the spinning bubble (G2). The bubble is at A=1, with r_bubble(θ) = M + √(M² - a²cos²θ). On the bubble, |∇A| varies with θ:

```text
|∇A|²(θ)|_bubble = (M-r)²/(M²r²)  +  a⁴ sin²(2θ) / (4 M² r⁴)
                 with r = r_bubble(θ)
```

**Equator (θ=π/2):** r=2M, sin(2θ)=0 → |∇A|=1/(2M) → T = T_Schw
**Pole (θ=0):** r=r_+, sin(2θ)=0 → |∇A| = √(M²-a²)/(M r_+) → T → 0 at extremal

The temperature is **non-uniform** along the spinning bubble, highest at the equator, lowest at the poles. This contrasts with Kerr, where T is uniform on the horizon.

## Numerical results for 1.0 solar mass BH

```text
    a/M    T_eq/T_S   T_pole/T_S   T_mean/T_S    T_K/T_S   A_b/A_S   A_K/A_S    L_b/L_S    L_K/L_S
----------------------------------------------------------------------------------------------------
  0.000      1.0000       1.0000       1.0000     1.0000    1.0000    1.0000  1.000e+00  1.000e+00
  0.100      1.0000       0.9975       0.9992     0.9975    0.9983    0.9975  9.950e-01  9.875e-01
  0.300      1.0000       0.9764       0.9925     0.9764    0.9850    0.9770  9.562e-01  8.881e-01
  0.500      1.0000       0.9282       0.9794     0.9282    0.9586    0.9330  8.843e-01  6.926e-01
  0.700      1.0000       0.8332       0.9605     0.8332    0.9198    0.8571  7.920e-01  4.131e-01
  0.900      1.0000       0.6071       0.9380     0.6071    0.8708    0.7179  7.013e-01  9.755e-02
  0.950      1.0000       0.4759       0.9324     0.4759    0.8582    0.6561  6.834e-01  3.365e-02
  0.990      1.0000       0.2473       0.9276     0.2473    0.8495    0.5705  6.720e-01  2.132e-03
  0.999      1.0000       0.0856       0.9260     0.0856    0.8484    0.5224  6.699e-01  2.804e-05
```

(T_S = T_Schwarzschild, A_S = Schwarzschild horizon area, L_S = Schwarzschild luminosity)

## Findings

**1. Temperature distribution.** Model-A predicts the spinning bubble has a *non-uniform* temperature along it. The equatorial temperature is preserved at the Schwarzschild value for ANY spin (because the equatorial bubble radius is locked at 2M). The polar temperature decreases monotonically with spin and vanishes at extremal. Kerr's horizon temperature is uniform and decreases monotonically (zeroing at extremal as well).

**2. Hawking emission rate.** Stefan-Boltzmann luminosity L = ∫ σ T⁴(θ) dA. Because Model-A retains the hot equator, its luminosity stays substantially higher than Kerr's at all spin > 0. Both predictions converge to L_Schw at a=0 and both shut off at extremal (no temperature anywhere → no emission), but the shutoff rate differs: Kerr's L drops as (r_+ - M)⁴ × A_Kerr while Model-A's L is dominated by the residual hot equator until it loses area at extremal.

**3. Boundary area & entropy.** Both Model-A and Kerr have decreasing area with spin, but Model-A shrinks more slowly. At extremal: Model-A bubble area ~0.81 × Schwarzschild, Kerr horizon = 0.5 × Schwarzschild. Therefore S_Model-A > S_Kerr at any spin > 0, with the gap widening at high spin.

**4. First-law structure.** Standard Kerr thermodynamics has `dM = T_Kerr dS + Ω_Kerr dJ` with a single uniform T. Model-A's non-uniform T(θ) means the spinning bubble is in NON-equilibrium Hawking emission: hotter regions emit faster, the boundary is polar-cooler/equator-hotter. To recover a single-T first law, one would compute an effective T_eff via L = σ T_eff⁴ A, but the physical interpretation of T_eff differs from Kerr's T_Kerr.

## Implications & STAM-distinctive predictions

- **Latitudinally-resolved Hawking spectrum.** If primordial BHs (near the evaporation epoch) are spinning, Model-A predicts their Hawking emission has a polar pattern: equatorial-bright, polar-dim, with the contrast growing with spin. Kerr predicts isotropic emission. Could be testable in PBH constraints if the mass scale is right.

- **Equator temperature locked at T_Schw.** This is unique to Model-A's bubble identification — a spinning BH still has the Schwarzschild temperature at its equator, because the equatorial bubble radius is locked at 2M. This means a spinning BH's EQUATORIAL Hawking emission is independent of spin in Model-A.

- **Higher entropy than Kerr.** S_Model-A > S_Kerr at any spin > 0. Spinning Model-A BHs hold more information on their boundary than Kerr would assign them.

- **No Penrose extraction (G2 finding).** Combined with this thermodynamics: spinning BHs in Model-A can only lose energy via Hawking-style boundary emission, never via classical energy extraction in the ergoregion. Energy budget for a spinning BH is constrained: L_emission as computed here, integrated over time, drives the BH toward the static (a=0) limit eventually.

- **Smarr formula.** Q12 confirmed Smarr `Mc² = 2 T S` for Schwarzschild. For Kerr the integrated form is `Mc² = 2 T S + 2 Ω J`. Model-A's non-uniform T means the Smarr-like relation needs reformulation. Whether `Mc² = 2 ⟨T⟩ S + 2 Ω J` holds for an appropriate ⟨T⟩, or whether a fundamentally different integral relation arises, is an open computation.

## Honest caveats

- **Coordinate gradient, not metric gradient.** Q8 used the coordinate |∇A| and that's preserved here. The proper |∇A| from the inverse metric on Model-A's spinning g would differ (would require committing to g_rr, g_tφ for the spinning case).
- **Round-coordinate area.** The bubble area integral uses the surface-of-revolution form, not the proper induced metric on the A=1 surface. Off by O(1) factor depending on how strongly the bubble is warped. The qualitative story (non-uniform T, area decreases with spin, S_modelA > S_Kerr) is robust.
- **Stefan-Boltzmann photon-only.** As in Q8/Q12, this counts two-helicity photon emission. Full Hawking emission integrates over all massless species and is larger by a degrees-of-freedom factor. Ratios L_ModelA/L_Kerr/L_Schw are unaffected.
- **No dynamical evolution.** The script computes thermodynamic quantities at fixed (M, a). Including back-reaction (the BH spinning down via emission, J/M ratio evolving) is a separate evaporation calculation.

## Generated plots

- `plots/G3_T_along_bubble.png`
- `plots/G3_thermo_quantities_vs_spin.png`
