# G14: BH Horizon Growth vs Big Bang Expansion (Wild Theory)

## Setup

User's wild theory: a stellar collapse to singularity 'resets' matter and inverts space, with the horizon expansion as something that might be analogous to the Big Bang. Connects to mainstream ideas (Smolin's CNS, Penrose's CCC, Big Bounce, white-hole cosmologies) without committing to any of them.

STAM commits to: matter accumulates on the A=1 bubble surface, no interior. The bubble grows from r=0 to r=Rs as matter joins.

## Numerical comparison

**BH horizon growth during stellar collapse:**

```text
Case                                                           Rs (m)    tau_ff (s)
------------------------------------------------------------------------------------------
Stellar BH (10 M_sun, R = 5×10^9 m progenitor)               2.95e+04      1.08e+04
Massive BH (40 M_sun, R = 1×10^10 m progenitor)              1.18e+05      1.52e+04
PBH formation (10^15 g, horizon-scale at t = 10^-23 s)       1.49e-15      2.23e-23
```

**Big Bang expansion eras:**

```text
Era              Functional form    Expansion factor
----             ---------------    ----------------
Inflation        a ~ exp(H_inf t)   ~10^26 (60 e-folds)
Radiation        a ~ t^(1/2)        ~10^4 (BBN to recomb)
Matter           a ~ t^(2/3)        ~10^3 (recomb to today)
Dark energy      a ~ exp(H_0 t)     1+ (today onward)
```

## Mathematical similarities and differences

**Similarities:**
- Both involve a 'horizon-like' boundary that grows or recedes.
- Both are characterized by a single mass-energy scale (BH mass M; universe content ρ_total).
- Both have a characteristic timescale set by the mass-energy (τ_BH ~ 2GM/c³; τ_universe ~ 1/H).
- For BH formation in time-reverse (white hole) and Big Bang looked at from outside, the mathematics partially aligns.

**Differences:**
- BH horizon SATURATES at Rs; Big Bang scale factor has no upper saturation (continues growing).
- BH formation timescale: minutes to hours for stellar collapse; BB cosmology unfolds over ~13.8 billion years.
- BH horizon expansion is driven by mass *accreting onto* a fixed location; BB expansion is the *spatial fabric stretching* between fixed-content matter.
- The 'inside-out' picture only works if you redefine the time direction (white-hole analog).
- BH horizon growth is sub-luminal (matter must fall in); BB expansion can superluminally separate distant points (no info transfer, but space itself stretches).

## STAM's take

STAM Model-A doesn't commit to BH-as-baby-universe. The framework's specific claim is that A=1 is the manifold edge — there is no 'inside' to be a separate universe. Matter accumulates on the boundary, full stop.

However, the framework is *compatible* with the wild-theory extension if one wanted to add it:
- The bubble's outer surface (where matter accumulates from outside) could be one hologram
- A complementary 'inner' surface (logically distinct, but ontologically denied in current STAM) would be the baby-universe interior
- The two-hologram parking-lot idea (your 2-side-of-horizon intuition) sits on this edge.

Mathematically, STAM doesn't predict any specific BH-to-BB mapping. The two processes have different drivers (matter accretion vs. cosmic expansion) and different functional forms (bubble saturation vs. monotonic scale-factor growth).

## Mainstream connections

**Smolin's Cosmological Natural Selection (1992):** every BH spawns a baby universe with slightly different physical constants. This is the closest mainstream proposal to the wild-theory intuition. Predicts that the constants we observe should be near-locally-optimized for BH production. Untested but interesting.

**Penrose's Conformal Cyclic Cosmology (CCC):** at the universe's max extent (when only photons remain after all BHs evaporate), the conformal structure resets and becomes the next Big Bang. Specific predictions for CMB anisotropy patterns; modest empirical support, contested.

**Big Bounce / loop quantum cosmology:** quantum gravity at the BB singularity replaces it with a non-singular bounce — the universe before BB was a contracting cosmos (which could be the inside of a BH from a parent universe).

**White hole cosmologies:** the Big Bang as the time-reverse of BH collapse. Various proposals (Markov, Frolov, etc.).

## Honest framing

The wild-theory direction has interesting mainstream company and STAM is structurally compatible with it (the two-hologram idea sits naturally there). But STAM doesn't currently make specific BB-from-BH predictions, and the mathematical analogies break down quickly under detailed comparison.

Worth keeping in the parking-lot file — alongside the multiverse and dark-energy-as-outer-hologram intuitions — but not in the active research line. The current framework's distinguishing predictions are at strong-field (G1), quantum (F6), and cosmological (G6/G11/G13) scales, not at BB-from-BH scales.

## Generated plots

- `plots/G14_bh_bigbang_expansion.png`
