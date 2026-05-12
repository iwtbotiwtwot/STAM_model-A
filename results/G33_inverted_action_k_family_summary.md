# G33: Inverted action approach -- k(A) family from A-fundamental view

**Date:** 2026-05-11 (evening)

**Swing in the dark:** treat A as fundamental, derive what k(A) can be.

## The setup

Standard scalar-tensor approach (G25-G27): assume g_munu fundamental, A is scalar field, derive Z(A) and V(A) from action variation. Result: ghost regions and Rs-dependent V.

Inverted approach (this script): A is fundamental, g_munu constructed from A via g_tt = -(1-A)c^2 and g_rr = 1/k(A). Question: what k(A) is consistent with the framework's structural commitments?

Constraints from framework:
- Weak-field GR at first order in A: k(0)=1, k'(0)=-1
- Horizon at A=1 with proper-time divergence: k(A) ~ (1-A)^p near A=1, p>=2
- Pair structure (1-A^2) appearing throughout framework results

Natural family satisfying all three: **k_n(A) = (1-A)(1-A^2)^n** for n>=1.

## Family analysis

| n | k_n(A) | weak-field GR | proper time div | tau_n/tau_GR (QNM) | NEC crossover |
|---|---|---|---|---|---|
| 0 | (1-A)(1-A^2)^0 | 1, -1 | no | 1.0000 | none |
| 1 | (1-A)(1-A^2)^1 | 1, -1 | YES | 1.3416 | 0.5000 |
| 2 | (1-A)(1-A^2)^2 | 1, -1 | YES | 1.8000 | 0.4400 |
| 3 | (1-A)(1-A^2)^3 | 1, -1 | YES | 2.4150 | 0.3992 |
| 4 | (1-A)(1-A^2)^4 | 1, -1 | YES | 3.2400 | 0.3689 |
| 5 | (1-A)(1-A^2)^5 | 1, -1 | YES | 4.3469 | 0.3452 |

**All n >= 1 satisfy the framework's stated constraints.** Model-A's n=2 is ONE MEMBER of a 1-parameter family. The 'three principles force k(A)' claim in PDF Section 4 is overstated -- three principles narrow the form to this family but don't uniquely specify n=2.

## QNM eikonal damping ratio

For static spherical metric with the framework's g_tt and g_rr = 1/k_n:

```
tau_n / tau_GR = (9/5)^(n/2)
```

This is a concrete framework prediction depending on n. Specific values:

- n=0 (Schwarzschild): tau/tau_GR = 1 (no departure)
- n=1: tau/tau_GR ≈ 1.342
- n=2 (Model-A): tau/tau_GR = 9/5 = 1.800 -- PDF's G1 prediction
- n=3: tau/tau_GR ≈ 2.415
- n=4: tau/tau_GR = 81/25 = 3.240
- n=5: tau/tau_GR ≈ 4.350

LIGO current precision is tau_obs/tau_GR ~ 1.0 +/- 0.2 (mass-independent).
This nominally excludes n>=1 in the eikonal approximation. However:
- LIGO BHs are spinning (Kerr-like); Model-A spinless analog is not yet derived
- Exact (non-eikonal) Regge-Wheeler computation has not been done for the family
- So this comparison is indicative, not falsifying

## What picks out n=2 (current Model-A commitment)?

Heuristic structural arguments for n=2 (NOT derivations):
- Two-face refinement: 'two faces x pair structure' -> (1-A^2)^2
- G18 ledger entropy decomposition: (2x2)=4 area-per-entry   (but G18 assumed n=2 to start, so this is potentially circular)
- Bekenstein-Hawking S = A/(4*L_P^2) is reproduced under n=2 via thermodynamic argument
- (4π × 3) decomposition of A_0 has integer 3 but doesn't directly fix n

None of these UNIQUELY force n=2 vs n=1 or n=3. The framework's commitment to n=2 is currently structural-aesthetic, not derived.

## What we learn from the miss

**1. The Lagrangian-gap remains.** Inverting the action principle doesn't close it -- it just reveals that the framework's constraints admit a 1-parameter family of k_n forms. The Lagrangian (if it exists) would need to fix n through some additional structural principle we haven't yet identified.

**2. Model-A's specific n=2 needs a fixing point.** Two candidates:
- (a) Observational: a clean measurement of tau/tau_GR for a spinless BH ringdown
- (b) Theoretical: an additional structural principle that uniquely picks n=2

**3. Pair structure (1-A^2) IS load-bearing.** It appears in all viable k_n forms; n=0 (no pair factor) gives pure GR, no STAM departure. The number of pair factors n is the open question.

**4. The (9/5)^(n/2) scaling is a concrete framework prediction.** Different n give measurably different QNM signatures. This is the framework's most direct exposure point to falsification once the spinless-analog issue is sorted out.

## Where this points us next

- **Push G1 toward exact (non-eikonal) QNM** for the family {k_n}. The exact spectrum may discriminate more sharply than the eikonal estimate.
- **Push G18 toward area-per-entry derivation WITHOUT assuming n=2.** If a first-principles ledger count gives 4 area/entry, n=2 is forced; if it gives a different number, n is what that count picks out.
- **Update PDF v0.4** to reflect that k(A) is a 1-parameter family, not a uniquely-forced form. Model-A's n=2 is a working commitment within the family.

## Plot

![G33 k family](../plots/G33_k_family.png)

Left: k_n(A) shapes for n in {0, 1, 2, 3, 4, 5}. Model-A's n=2 is bold blue. All n>=1 satisfy the framework's stated constraints; n=2 is one choice among them.

Right: QNM eikonal damping ratio vs n. The (9/5)^(n/2) scaling parameterizes a continuum of framework predictions. Model-A's n=2 commitment gives the 1.80 value cited in the PDF.

## Reading

**The swing missed in an informative way.** We didn't derive k(A) = (1-A)(1-A^2)^2 from first principles -- we discovered that the framework's constraints admit a 1-parameter family, and Model-A's commitment to n=2 is structural-aesthetic rather than forced. The miss tells us where the actual fixing point would have to come from: either an exact QNM observation, or a non-circular derivation of the entropy decomposition.

This is a sharpening of the strong-field-departure question (see project_strong_field_departure_question.md). The question is no longer 'what makes k(A) concrete?' but 'what picks n=2 out of the family?'.