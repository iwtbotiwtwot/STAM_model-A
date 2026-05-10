# G15: A_0 Sensitivity Test

## Motivation

The framework currently commits to A_0 = 1/(12π) on structural grounds: the (4π × 3) decomposition where 4π is the Q8/Q10 thermal prefactor and 3 is hypothesized to be spatial dimensionality. Because the framework author has flagged a self-recognized bias toward seeing 'thirds' structure, this script tests directly whether the precise value 1/(12π) is doing real work in the framework, or whether any A_0 in roughly [0.02, 0.04] would produce equivalent observational predictions.

If A_0 = 1/(12π) sits in a sharp peak where downstream observables are tightly fit, the structural commitment is load-bearing. If it sits on a broad shelf where any value in the range works about equally well, the commitment is aesthetic and we should not get attached to the specific decomposition.

## Setup

- H_0 = 73.04 km/s/Mpc (SH0ES local distance ladder)
- Ω_m_target = 0.315 (PBH-DM-compatible matter content)
- Historical empirical bridge term = 354.95 Mly (STAM's own prior catalog calibration)
- Observed θ_⋆ = 0.010410 rad (Planck)
- No-LoS θ_⋆ at Ω_m = 0.315: 0.010813 rad (+3.87% from observed)
- f_LoS required for full closure (data-driven, A_0-independent): 1.0387× (i.e. A_LoS_effective = 0.0387)

Candidate A_0 values: 1/(8π), 1/(10π), 1/(12π), 1/(14π), 1/(16π).

## Results

```text
A_0 candidate   A_0 value    b (Mly)    b offset   closure %   amp / A_0   physical?
------------------------------------------------------------------------------------------
      1/(8π)    0.039789     532.65     +50.06%      97.41%       0.974  NO (amp<1)
     1/(10π)    0.031831     426.12     +20.05%      82.71%       1.217         yes
     1/(12π)    0.026526     355.10      +0.04%      69.28%       1.461         yes
     1/(14π)    0.022736     304.37     -14.25%      59.60%       1.704         yes
     1/(16π)    0.019894     266.33     -24.97%      52.30%       1.948         yes
```

## Reading the table

**Bridge term sensitivity.** `b = A_0 × c/H_0` is strictly linear in A_0. Only A_0 = 1/(12π) lands within ±5% of STAM's own historical empirical bridge term. Adjacent candidates (1/(10π), 1/(14π)) miss by ~20%; outer candidates (1/(8π), 1/(16π)) miss by 25-50%.

**CMB closure sensitivity.** The data-required `f_LoS` does not depend on A_0 — at Ω_m = 0.315, the data wants `f_LoS ≈ 1.0387` regardless of where A_0 sits. The A_0 sensitivity enters as 'how much structure amplification above the A_0 baseline is needed?'. If the required amplification factor is < 1, the framework would need to anti-amplify A below A_0 along the line of sight — physically impossible (amplification by cosmic structure can only raise A above the void-dominated background). So `amp < 1` rules a candidate out on the CMB side.

**Joint constraint.** Candidates passing both observables (bridge match within ±5% AND amp > 1) are the only ones the framework can sit on.

## Verdict

- Bridge passers (|offset| < 5%): **1/(12π)**
- CMB passers (amp > 1, physical): **1/(10π), 1/(12π), 1/(14π), 1/(16π)**
- Joint passers (both constraints): **1/(12π)**


**Reading.** Of the candidates tested, only **1/(12π)** survives both the bridge-term constraint and the CMB-closure physical-amplification constraint. The structural commitment to A_0 = 1/(12π) is doing real work — it's not arbitrary, it's not interchangeable with nearby (mathematically similar) values, and it sits at the unique point where both internal STAM observables are consistent. The bridge term is the tighter constraint (linear in A_0); the CMB constraint is the sign check (amp > 1 required). The two together pin A_0 to within a narrow band, and 1/(12π) is the structurally-elegant value inside that band.

**What this does NOT establish.** This sensitivity test does not derive the '3' factor from first principles — it shows that 1/(12π) is observationally distinguished, not that it's theoretically necessary. Different mathematical decompositions could give the same numerical value (e.g. 0.0265 doesn't have to factor as (4π × 3)). What it does establish is that the *value* 0.0265 is load-bearing — it's not interchangeable with 0.020 or 0.040 in a way that the framework would tolerate.


## Generated plots

- `plots/G15_a0_sensitivity.png`
