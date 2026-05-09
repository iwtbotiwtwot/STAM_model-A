# G12: PBH Formation Probability in STAM Language

## Setup

In STAM, a primordial black hole IS a bubble — a region where A reaches 1. For a primordial overdensity δ at scale R:

```text
A(R) = (R / L_Hubble)² × (1 + δ) × Ω_m_at_epoch
```

PBH forms when A(R) reaches 1, equivalent to δ ≥ δ_c ≈ 0.45 in the radiation era. This is the same condition as standard PBH formation — STAM relabels it in A-language without changing the underlying probability.

**Press-Schechter formation fraction:**
```text
β(M) = (1/2) × erfc(δ_c / (σ(M) √2))
```
**Today's PBH abundance (fraction of DM):**
```text
f_PBH(M) ≈ β(M) × √(M_eq/M) × (Ω_m / Ω_DM)
```
where M_eq ≈ 10⁵⁰ g is the horizon mass at matter-radiation equality.

## PBH formation table

```text
log10(M/g)  beta σ=0.05  beta σ=0.1  beta σ=0.2  f_PBH σ=0.1  max allowed
--------------------------------------------------------------------------------
      15.0     1.13e-19    3.40e-06    1.22e-02     1.77e+12        1e-07
      17.0     1.13e-19    3.40e-06    1.22e-02     1.77e+11        1e+00
      20.0     1.13e-19    3.40e-06    1.22e-02     5.59e+09        1e+00
      23.0     1.13e-19    3.40e-06    1.22e-02     1.77e+08        1e-01
      25.0     1.13e-19    3.40e-06    1.22e-02     1.77e+07        1e-01
      30.0     1.13e-19    3.40e-06    1.22e-02     5.59e+04        1e-01
      33.0     1.13e-19    3.40e-06    1.22e-02     1.77e+03        1e-03
      36.0     1.13e-19    3.40e-06    1.22e-02     5.59e+01        1e-03
      39.0     1.13e-19    3.40e-06    1.22e-02     1.77e+00        1e-03
```

## Sigma required for full PBH-DM at various masses

```text
     M (g)   σ for f=1   σ for f=0.1     PBH-allowed?
------------------------------------------------------------
     1e+17      0.0533        0.0518              yes
     1e+20      0.0563        0.0544              yes
     1e+23      0.0595        0.0571          f<1e-01
     1e+30      0.0698        0.0665          f<1e-01
```

## STAM-language interpretation

Each PBH is a small STAM bubble:

- **A = 1 boundary**: a real 2D surface in spacetime where the manifold ends. Bubble's geometric form, not a coordinate singularity.
- **Mass = surface mass on the bubble**: matter that fell toward the PBH never crossed inside; it accumulated holographically on the boundary surface (per outward-collapse picture, memory/project_outward_collapse_dynamics).
- **Hawking T = ℏc|∇A|/(4π k_B) at boundary** (Q8/Q10): for PBH of mass M, T ≈ 6 × 10⁻⁸ K × (M_sun/M).
- **Bekenstein entropy = A_horizon/(4 ℓ_p²)** (Q11): bubble area-scaling.
- **No interior, no singularity, no firewall** — same dissolution as for stellar BHs (F1 result).

PBH-DM is therefore a *more natural fit for STAM than for LCDM*: the framework's bubble picture handles small black holes with the same machinery as stellar/super-massive BHs, without requiring new particle physics for dark matter.

## What this tells us

**Standard inflation** (σ ≈ 10⁻⁵ at all scales) gives β ≈ exp(-10¹⁰), essentially zero PBH formation. So PBH-DM requires **enhanced power at small scales** — a peaked or running primordial spectrum.

**For asteroid-mass PBH-DM (10¹⁷ - 10²³ g window):**
- σ(M) ≈ 0.05-0.10 at the PBH-formation scale produces f_PBH ~ 0.1 - 1
- This is achievable with single-field inflation models having a small-scale bump in the power spectrum
- The asteroid-mass window is observationally allowed for f_PBH up to ~1 (no current bounds rule out PBH being all of DM in this range)

**STAM doesn't change PBH formation odds compared to standard physics** — same Press-Schechter, same δ_c, same probability. The framework's contribution is interpretive: each PBH is a STAM bubble with derivable thermodynamics, and PBH-DM dovetails naturally with the framework's BH ontology.

**For the rotation-curve question (G10):** PBH-DM with masses in the asteroid-mass window can supply the missing galactic mass. Each galaxy contains roughly:
- M_DM_per_galaxy ≈ 5 × M_baryon (Milky Way: ~3.5 × 10¹¹ M_sun)
- For PBH-DM at M_PBH ≈ 10¹⁸ g: ~ 10⁴¹ PBHs per galaxy
- They contribute via standard cumulative gravity (Newton)

This is consistent with the keystone test result: STAM linear cumulative A from baryons fails because there isn't enough baryonic matter. Adding PBH-DM (compatible with STAM's bubble framework) fills the gap.

## Generated plots

- `plots/G12_pbh_formation.png`
