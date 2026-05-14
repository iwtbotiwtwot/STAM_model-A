---
name: Microscopic SU-write fragments already in the repo
description: The framework has several committed microscopic-scale fragments (1 SU = A_0, 1 SU = 4 Planck cells on horizon, 1 SU = 12π Planck cells in bulk, Γ_res Lindblad form, saturation-pair Hawking) but no formal stochastic process or coarse-graining theorem. Any future microscopic SU-write model must reproduce these fragments.
metadata:
  type: project
---

# Microscopic SU-write fragments

After a thorough sweep (Explore agent + grep, 2026-05-13 evening) for any formal microscopic SU-write process or coarse-graining theorem, the finding is: **no formal stochastic process (Poisson/Markov/lattice/spin-foam) and no coarse-graining theorem (RG/hydrodynamic/path-integral) exists in the repo.** But there are several microscopic structural fragments that any future formal model must reproduce.

## Foundational identity: A_0 = 1 SU

The three "1 SU = ?" readings below are NOT independent. They are scale-conjugate views of ONE primitive: **A_0 = 1 SU** is the natural-unit identity (Sean confirmed 2026-05-13 evening).

| Context | Statement | Reading |
|---|---|---|
| Natural units | A_0 = 1 SU | SU is the ruler; A_0 is one tick |
| Bulk volumetric | 1 SU = 12π Planck cells | Min cell = 4πD = 12π Planck cells at D=3 |
| Horizon holographic | 1 SU = 4 Planck cells | α_H × gravity-bridge = 2 × 2 (G59) |
| Density normalization | A_0 = 1/(12π) | Same primitive, Planck-cell-density units |

The 12π → 4 dimensional reduction is geometric (bulk D-volume → horizon (D-1)-area). Ratio 12π/4 = 3π. Bekenstein-Hawking entropy reads as S_BH = total SU count on horizon when each SU occupies its 4 Planck areas. G59 numerical verification: ratios = 1.000000 across all masses.

Do NOT treat the three readings as separate facts. They are scale-faces of the single primitive A_0 = 1 SU.

## The fragments (in order of structural strength)

### 1 SU = 4 Planck cells on the horizon (G59)
Strongest microscopic identification. See [results/G59_entropy_under_ledger_channel_summary.md](../results/G59_entropy_under_ledger_channel_summary.md).

  alpha = alpha_H * gravity-bridge = 2 * 2 = 4 Planck cells per manifold-support entry

Decomposition:
- alpha_H = 2: two-face horizon-pair (1 cell on inner face + 1 cell on outer face)
- gravity-bridge = 2: 2 Planck-cell depth per channel-cell (from A = 2GM/c^2 r factor of 2)

Reproduces Bekenstein-Hawking area-per-bit factor numerically (ratio = 1.000000 across all masses tested). Both factors are independently load-bearing elsewhere (alpha_H in Beta(D+1, 2), gravity-bridge in A definition + resolution rule + bridge term), so this isn't curve-fitting.

**Each SU write has a well-defined Planck-scale footprint on the holographic boundary: 4 cells.**

### 1 SU = 12π Planck cells in the bulk vacuum
See [[project_v04_priority_SU_equals_A0]]. A_0 = 1/(12π) interpreted as "1 SU per minimum cell, where 1 minimum cell = (4π × D) = 12π Planck cells at D=3". Volumetric, bulk-density reading.

### 1 SU = A_0 per quantum interaction
The atomic unit commitment, [[HANDOFF_2026_05_13_unified_framework]]. Each quantum physical interaction resolves exactly 1 SU = A_0 of substance. Distinct from the geometric "what does an SU look like" questions above.

### Γ_res Lindblad-analog form
README Open Problem #5a; [[HANDOFF_2026_05_13_unified_framework]] section on Γ_res structural form:

  Gamma_res = sum_mu <L_mu^dagger L_mu> * D_mu        (dynamical)
            = -d ln(C)/d tau                          (diagnostic)
  0 <= Gamma_res <= (A/A_0) / tau_P                   (bounds)

Form committed; specific L_mu and D_mu for particular interactions open. The upper bound (A/A_0)/tau_P is the carrying-capacity ceiling on the write rate.

### N_SU support vs dN_write rate distinction
G68, [[HANDOFF_2026_05_13_unified_framework]]. The structural carrying capacity N_SU = A/A_0 is DIFFERENT from the actual write rate per Planck tick dN_write = Gamma_res * d tau. Earth surface: N_SU ~ 3e-8 (thin support). Near horizon: N_SU -> 38 (saturation). Hawking rate at horizon: ~10^-122 per cell per Planck tick (far below "1 write per tick per cell").

Do NOT claim "one write per Planck tick per cell." Support is structural; rate is physical.

### Hawking saturation-pair mechanism
[[project_hawking_mechanism_native]]. At A=1, full cell forces exchange-shaped pair writes: each event is outward-write + content-reduction. Pair structure is the saturation requirement itself, not a Bogoliubov artifact.

## What's missing (Target 1 + Target 2 status)

**No formal stochastic process exists.** Absent: Poisson point process, Markov chain / master equation, spin foam, spin network, lattice / cellular automaton, explicit "writes form a Poisson process" formalism.

**No coarse-graining theorem exists.** Absent: RG flow, hydrodynamic limit, path integral over write configurations, mean-field derivation, scaling-limit theorem with error bounds.

The README's Open Problem #1 explicitly notes "microscopic derivation of the constrained shell-count action from underlying SU-write dynamics and Gamma_res coarse-graining remains a deeper foundational project."

## Concrete candidate microscopic model

The fragments above are consistent with (but don't uniquely select) the following structural proposal:

> Spacetime Poisson point process of SU-write events.
> Local rate density Gamma_res(x) with bound (A/A_0)/tau_P.
> Each event deposits 1 A_0 of substance in a 12pi-Planck-cell bulk neighborhood
> OR registers 4 Planck cells of horizon area at the boundary.
> Events at A=1 are forced into exchange pairs (saturation rule).
> Coarse-graining: Sigma(x) = path-integrated count over Planck volume
> -> smooth field at scales >> ell_P.

This proposal is *consistent* with every fragment above but is not yet *derived*. Next-step script would be **G124_su_write_poisson_model.py**: instantiate this Poisson process, derive Gamma_res hazard from p(y)/F(y), and check the coarse-grained limit reproduces the constrained shell-count action.

## How to use this memory

- When future-Claude searches for "microscopic SU-write model," this memory consolidates what exists vs what's open.
- The 1 SU = 4 Planck cells result (G59) is the strongest microscopic structural identification — any future formalism must reproduce it.
- Two different "SU footprint" readings coexist: 4 Planck cells (horizon area) vs 12pi Planck cells (bulk volume). Both are correct in their respective contexts.
- Sean has surfaced the "4 Planck cells" reading multiple times (recall on 2026-05-13 evening) — it is load-bearing in his intuition.

## See also
- [[project_three_strain_points_fix_strategy.md]] — the broader fix strategy this memory supports
- [[project_v04_priority_SU_equals_A0]] — the bulk volumetric reading
- [[project_hawking_mechanism_native]] — the saturation-pair mechanism
- [[HANDOFF_2026_05_13_unified_framework]] — original Γ_res commitment
- [[project_five_primitive_architecture]] — where the SU-write process would close the deepest primitives
