# STAM Model-A F5: Gravitational Wave Polarization Structure

## The attack

GR predicts gravitational waves with two transverse-tensor polarizations: `h_+` and `h_×`. LIGO/Virgo/KAGRA observations have confirmed these tensor modes and constrained any non-tensor polarizations (scalar / vector / longitudinal) to less than ~10% of the total GW power for typical compact-binary mergers (representative value; specific bounds vary by analysis).

Pure scalar theories of gravity (Nordstrom 1913 and similar) predict only a longitudinal scalar 'breathing mode' and were ruled out historically. The opponent asks: bold STAM in Framework C is scalar gravity (with bold-STAM metric construction rule). Doesn't this place bold STAM in the same falsified class as Nordstrom?

## The honest analysis

Linearizing the bold-STAM construction rule around flat space (small `δA` perturbation):
```text
g_tt = -(1-A) c²       =>   δg_tt = +c² δA       (time-time component)
g_rr = 1/[(1-A)(1-A²)²] =>   δg_rr ≈ +δA          (radial component)
g_θθ = r²              =>   δg_θθ = 0            (angular — UNCHANGED)
g_φφ = r² sin²θ        =>   δg_φφ = 0
```

For a GW propagating in the z direction, in Cartesian coordinates at the detector:
```text
h_tt = +c² δA          (scalar / time-time mode)
h_zz = +δA             (longitudinal — in direction of propagation)
h_xx = h_yy = 0        (transverse components vanish — NO h_+ MODE!)
h_xy = 0               (NO h_× MODE!)
```

**The TT-gauge tensor modes (h_+, h_×) are exactly zero in bold STAM Framework C.** Bold STAM, in its simplest generalization to time-dependent perturbations, predicts only longitudinal scalar GW modes. This is the same problem Nordstrom's gravity had.

## Comparison with LIGO

- Bold STAM tensor fraction: 0.00  (LIGO requires ≥ 0.50). **FAIL** (False).
- Bold STAM scalar fraction: 1.00  (LIGO requires ≤ 0.10). **FAIL** (False).

**Verdict on F5: bold STAM in Framework C, with the simplest scalar generalization, is FALSIFIED by LIGO/Virgo/KAGRA observations of GW polarization.**

## Path forward — what bold STAM has to do to survive

F5 does not refute bold STAM as a *concept*. It refutes Framework C in its simplest form. The construction rule must be amended so that the metric has additional structure beyond what the scalar A provides directly. The natural amendments are:

**Option I: Scalar-tensor amendment.** Promote the metric to an independent dynamical field, coupled to A through the construction rule but with its own tensor degrees of freedom. Bold STAM becomes a scalar-tensor theory of gravity, structurally similar to Brans-Dicke or its generalizations. The TT modes (h_+, h_×) come from the metric's own dynamics; the scalar A provides an additional mode that LIGO has bounded but not ruled out at the few-percent level.

**Option II: A as a higher-rank field.** Promote A from a scalar to a tensor-valued field, where the construction rule maps a tensor A_μν to the metric. This is more elaborate and less STAM-native; the original framework treats A as a scalar.

**Option III: Amended construction rule with derivative couplings.** Allow the metric to depend on derivatives of A (∂A, ∂²A) in addition to A itself. Could potentially produce TT modes from how A propagates, even with A as a scalar.

Option I is the most direct fix. It does require giving up some of bold STAM's current cleanness — the metric is no longer a derived consequence of A alone. But it preserves all the local results (Q8-Q12 thermodynamics, SF3-SF5 strong-field structure, F1 evasion, F3 evasion) because those don't depend on the tensor mode structure of GWs. The amendment shows up only when GWs are computed.

## Honest verdict

**F5 is the first fatality bold STAM does not survive in its current form.** Framework C, as we specified it in SF6, predicts the wrong GW polarization structure. LIGO observations directly contradict this prediction.

However, F5 is *not* a refutation of bold STAM as a research program. It is a specific failure of the simplest scalar generalization of the bold-STAM construction rule. The fix (scalar-tensor amendment) is well-understood theoretically — it just hasn't been written down for bold STAM specifically. Doing so is real research work, comparable in scope to what was needed to write down Q8-Q12 or SF3-SF5.

The framework's status, after F5:
- *Survived cleanly:* F1 (Penrose-Hawking), F2 (equivalence principle).
- *Survived after Framework C commitment:* F3 (Birkhoff).
- *Failed in current form, with clear path to amend:* F5 (GW polarization).

**Bold STAM is unrefuted as a conceptual framework**, but it is not yet a complete theory of gravity. The completion requires committing to scalar-tensor structure (Option I above) and rebuilding the perturbation theory accordingly.

## Generated plots

- `plots/F5_polarization_patterns.png`
- `plots/F5_polarization_fractions.png`
