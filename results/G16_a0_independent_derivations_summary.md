# G16: Independent A_0 Derivation Candidates

## Question

Do any independent physical-principle derivations converge on A_0 = 1/(12π) ≈ 0.0265 *without* referencing the value 0.0265 in their setup?

If two or more independent routes land on the same value, the structural commitment is *overdetermined* — multiple physical arguments converge on it. If all routes give different values, the (4π × 3) decomposition stands alone as the structural anchor.

## Methodology

Each candidate must:
1. Start from a physical principle (Newton-shell, uniform sphere, holographic, etc.) that does NOT reference A_0 = 0.0265 in its setup.
2. Compute the value the principle predicts.
3. Compare to the target.

Two candidates are flagged as **tautological** for completeness:
- The (4π × 3) decomposition itself — that's the argument under test.
- V_3's minimum at A_0 — α/β was tuned by hand to put it there.

These don't count as independent confirmations.

## Results

```text
                                         Candidate   Predicted A_0     Ratio      Offset                     Verdict
--------------------------------------------------------------------------------------------------------------------
              (4π × 3) decomposition  [structural]      0.02652582    1.0000      +0.00%                tautological
                       V_3 minimum  [tautological]      0.02652582    1.0000      +0.00%                tautological
             Newton-shell with Ω_b mass at horizon      0.04930000    1.8586     +85.86%       no — off by tens of %
             Newton-shell with Ω_m mass at horizon      0.31500000   11.8752   +1087.52%          no — off by factor
            Newton-shell with Ω_DE mass at horizon      0.68500000   25.8239   +2482.39%          no — off by factor
                 Center of uniform critical sphere      0.50000000   18.8496   +1784.96%          no — off by factor
         Volume-average of uniform critical sphere      0.20000000    7.5398    +653.98%          no — off by factor
        Radial-LoS through uniform critical sphere      0.33333333   12.5664   +1156.64%          no — off by factor
              Radial-LoS × Ω_b  [near-match probe]      0.01643333    0.6195     -38.05%       no — off by tens of %
              Radial-LoS × Ω_m  [near-match probe]      0.10500000    3.9584    +295.84%          no — off by factor
```

## Verdict

- Independent candidates tested: **8**
- Converged within ±5%:  **0**
- Near match (±5-20%):   **0**
- Off by ≥20%:           **8**


**No independent derivation lands on A_0 = 1/(12π) within ±5%.** This is the result, reported honestly.

What this means for the framework's status on A_0:

- The numerical match to STAM's historical bridge term (0.04%, G15) remains real, sharp, and load-bearing.
- The (4π × 3) structural decomposition remains the single structural argument for *why* the value is 1/(12π) specifically.
- V_3's minimum at A_0 is by α/β construction, so it doesn't constitute independent evidence.

What this *doesn't* mean:

- A_0 = 1/(12π) is not invalidated by this result. The bridge match and the decomposition are real evidence on their own.
- Future independent derivations (a proper Lagrangian-for-A derivation, a holographic argument, an F3-extended Friedmann first-principles calculation) might still converge. This script tests only the candidates available with current framework machinery.

Closest independent candidates from the set tested:

- **Radial-LoS × Ω_b  [near-match probe]**: A_0 = 0.016433 (-38.0% from target). Script 48 Prescription 7 — closest existing near-match

- **Newton-shell with Ω_b mass at horizon**: A_0 = 0.049300 (+85.9% from target). M_shell = Ω_b × M_H (baryonic fraction at Hubble boundary)

- **Radial-LoS × Ω_m  [near-match probe]**: A_0 = 0.105000 (+295.8% from target). Script 48 Prescription 6


## How to read this honestly

A negative result on multi-route convergence is not the same as the framework being wrong about A_0. The G15 sensitivity test established that A_0 = 1/(12π) is observationally distinguished — only this value satisfies both the bridge-term and CMB-closure constraints. That result stands independently of whether multiple theoretical derivations also converge on the value.

What this G16 result tells us: the *theoretical* support for A_0 = 1/(12π) currently rests on the (4π × 3) decomposition argument, not on multiple independent first-principles derivations. The open theoretical work — the Lagrangian for A, the first-principles derivation of the spatial-dimensionality factor 3 — would convert the single argument into multiple, if successful. Until then, A_0 is observationally distinguished and theoretically anchored on one argument.


## Generated plots

- `plots/G16_a0_independent_derivations.png`
