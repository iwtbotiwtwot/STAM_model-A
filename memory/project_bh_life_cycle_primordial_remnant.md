---
name: BH life cycle — outer-face dynamics and primordial-mass stable remnants
description: Full BH life cycle picture under two-face refinement: formation presses primordial mass to inner face (permanent); accretion builds outer face; Hawking radiation drains outer face; horizon size tracks outer-face content; final state is a stable primordial-mass remnant. BHs do NOT fully evaporate. Worked out 2026-05-11 with Sean.
type: project
---

Worked out with Sean 2026-05-11. End-to-end picture of black hole formation, growth, decay, and final state under Model-A's two-face refinement.

## The life cycle, end to end

**Stage 1 — Formation.**
Matter concentrates such that A reaches 1 *beyond* its physical surface. The horizon (A=1 surface) forms outside the matter. The matter that triggered formation is pressed onto the **inner face** of the horizon as the A=1 surface establishes itself. This becomes the **primordial mass M_0** — the mass associated with the formation event. The inner face holds this content **permanently**.

**Stage 2 — Active growth.**
Subsequent matter accretes onto the **outer face** of the horizon (from outside). Outer-face content M_acc(t) accumulates. The BH's *total* horizon size grows because the boundary's total content is M_0 + M_acc(t). From an outside observer, the BH gets larger as it accretes — this is just the outer face growing.

**Stage 3 — Steady state (optional).**
If accretion rate equals decay rate, the BH size is constant. Outer-face content flows in (from accretion) and out (as Hawking radiation) at the same rate. Horizon size is stable.

**Stage 4 — Decay.**
When the accretion rate drops below the decay rate, the outer face drains faster than it's replenished. The **horizon recedes** — the BH gets smaller from the outside. This continues as the outer-face content decreases.

**Stage 5 — Stable primordial remnant.**
Eventually the outer face is empty. The boundary has only inner-face content, M_0. Hawking radiation stops (no outer-face content available for pair-exchange). The horizon size stabilizes at the size determined by M_0 alone. The BH is now a **stable primordial-mass remnant** — it doesn't change size, doesn't emit Hawking radiation, and doesn't fully evaporate.

**Prediction: NO black hole ever fully evaporates.** Every BH eventually becomes a stable remnant equal to its formation mass.

## Horizon size as real-time accounting

A key insight: the horizon's size from outside reflects the *current* outer-face content. It is not a fixed geometric quantity tied to the BH's "underlying" mass. The horizon size is dynamic and tracks:

- Current accretion rate (additions)
- Current decay rate (subtractions)
- Cumulative outer-face content (the running total)

A BH actively accreting in a dense environment will appear to have a large horizon. The same BH in a quiet environment will gradually shrink its horizon as outer-face content drains. The inner face (M_0) sets a hard lower limit on horizon size, but the *current* horizon size is determined by total outer-face content plus M_0.

## Connection to other framework commitments

- **Two-face refinement (working commitment as of 2026-05-10):** the inner/outer distinction is what allows this life-cycle picture to work. Without it, full evaporation would be inevitable.
- **Information conservation:** inner face permanence means information about formation is never emitted. The information paradox dissolves cleanly: formation information stays put forever; only accretion information (outer-face content) gets radiated as Hawking radiation.
- **PBH-DM compatibility (G12, G13):** primordial black holes formed in the early universe have primordial mass M_0 set by formation. They have minimal outer-face accretion in cosmic voids. So they are stable to their full formation mass for all time — exactly what's needed for them to be dark matter.
- **Hawking mechanism (project_hawking_mechanism_native.md):** the pair-exchange mechanism at the saturated boundary drains the outer face. Inner face is not accessible to the pair-exchange (different category of content; held permanently by the formation event).

## Observational implications

- **Stellar and supermassive BHs:** Hawking timescales are so long that we wouldn't observe full evaporation in current astronomy regardless of whether BHs leave remnants. The remnant prediction is consistent with all current observations.
- **Asteroid-mass PBHs (10^17–10^23 g):** if these are dark matter (PBH-DM hypothesis), they need to be stable on cosmological timescales. Model-A says they are — they don't evaporate below their primordial mass. Stronger PBH-DM compatibility than standard physics.
- **No evidence of BH "death bursts":** standard Hawking predicts a final burst of energy as a BH's mass drops to zero. Model-A predicts no such burst — emission stops at the primordial remnant stage. Absence of such bursts in observation would be weak support for Model-A; non-detection is consistent with timescales being too long anyway.
- **BH size doesn't directly read mass.** The accretion-decay accounting picture means a BH's horizon size reflects its current outer-face state, not a fixed geometric mass. This is consistent with observed BH variability (active vs quiescent phases) but doesn't yet give a sharp testable prediction.

## How to apply

- When discussing BH "evaporation": clarify that Model-A's evaporation is outer-face drainage, not full disappearance. Use "evaporates to primordial remnant" or similar — never just "evaporates."
- When discussing BH "mass": distinguish primordial mass (M_0, set at formation, permanent) from accreted mass (M_acc(t), dynamic). Total mass = M_0 + M_acc(t).
- When discussing horizon size: it tracks current outer-face content; not a fixed property of the BH.
- For PBH-DM scenarios: PBHs are stable to primordial mass; this is feature, not bug.
- Information conservation: inner face is permanent, no info loss; framework-native resolution of the paradox.
- Don't claim "BHs are stable forever" without qualification — they may have ongoing outer-face dynamics; what's stable is the floor.

## Open quantitative pieces

- Accretion rate vs decay rate quantitatively: when does outer face balance, grow, or shrink? Needs specific input on accretion physics.
- Final remnant size for known BH types: stellar-mass formations from supernovae have specific M_0 values; supermassive BHs have less clear primordial mass histories (mergers of many smaller BHs?).
- Whether mergers preserve the inner faces of both progenitors (inner faces combine? remain distinct?) — open question for binary BH merger physics under two-face refinement. SF5 needs to be redone under this picture.
