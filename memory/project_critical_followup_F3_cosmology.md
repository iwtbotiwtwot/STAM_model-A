---
name: CRITICAL — F3-extended cosmology calculation
description: Author flagged this as critical to come back to. F3 exposed bold-STAM effective stress-energy (negative ρ in outer regions, positive near horizon). Cosmologically integrated, this might have dark-energy-like effects relevant to the bridge term b. Has NOT been computed yet.
type: project
---

**Critical follow-up flagged by author 2026-05-07.** Must come back to.

**The setup F3 left open:**

F3 (Birkhoff fatality) exposed that bold-STAM metric `g_rr = 1/[(1-A)(1-A²)²]` has non-zero effective stress-energy when interpreted in Einstein gravity:
- Outer exterior (A < 0.44): effective ρ < 0 — looks like dark energy
- Inner exterior (A > 0.44): effective ρ > 0 — looks like ordinary matter
- Mass function m(r)/M peaks at ~1.44 around A ≈ 0.5

**Why this is critical for distance / bridge term b:**

Script 31 (cosmological distance) computed only **Shapiro propagation** through cosmological A field — light traversal time. It found bold STAM is flatter than LCDM, doesn't help match catalogs.

But script 31 did NOT include **gravitational dynamics from the bold-STAM effective stress-energy**. Those are different calculations. The F3 finding suggests bold STAM might have a built-in dark-energy-like cosmological component that script 31 missed.

**What the calculation would do:**

1. Compute G^r_r and G^θ_θ for bold-STAM metric (radial and tangential effective pressure). G^t_t (energy density) already done in F3.
2. Set up cosmological-scale bold-STAM A field profile (probably linear superposition of A from many distributed sources, à la SF5).
3. Solve modified Friedmann equations with bold-STAM effective stress-energy.
4. Extract distance-redshift prediction.
5. Compare to D_adj_no_b, LCDM, and to catalog data.

**Caveat: must specify field equations first.**

Two interpretations of bold STAM give different cosmological predictions:
- "Einstein gravity + effective matter" interpretation: one answer.
- "Modified gravity (scalar-tensor or f(R))" interpretation: different answer.

Until bold STAM commits to a field-equation class, the calculation has to either be done both ways, or be done conditionally on one interpretation.

**Possible outcomes (all informative):**
- Match catalogs better than current STAM: major positive finding; bold STAM derives dark-energy-like effect, partially explains catalog discrepancy / b term.
- Match catalogs about the same: script 31's negative result confirmed.
- Match catalogs worse: bold-STAM cosmology broken on this front; need to rethink cosmological A profile.

**Estimated effort:** 2-3 hours of focused script writing.

**How to apply:**
- This is the highest-priority unfinished cosmological piece. When the author returns to STAM work, this should be near the top.
- Best done after committing to bold-STAM field equations (close F3 properly first), but a preliminary "Einstein-gravity + effective matter" version could be done independently for fast feedback.
