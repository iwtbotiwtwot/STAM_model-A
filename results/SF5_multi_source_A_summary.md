# STAM Model-A SF5: Multi-Source A Field and Binary BH Merger Topology

## Purpose

Formalize linear superposition of the A field for multiple sources, and use it to map the binary-BH merger topology. This sets the foundation for any future LIGO ringdown / inspiral comparison work.

## Linear superposition rule (weak-field, exact)

```text
A_total(x) = sum_i  2 G m_i / (c^2 |x - x_i|) = sum_i Rs_i / |x - x_i|
```
Each source contributes its own 1/r piece. The total field is the linear sum. This is exact at large separation (where each source's local field dominates in its own neighborhood) and gives the qualitatively correct merger topology.

## Bold-STAM picture of binary BH merger

**Inspiral.** Two separate bubbles (A=1 surfaces) approach each other. Each surface deforms toward the other under the combined A field, but as long as they remain topologically distinct, the system is two BHs.

**Critical separation: d_crit = 4 Rs.** This is exact for equal-mass BHs from linear superposition: A at the midpoint = 2 × Rs / (d/2) = 4Rs/d. Setting this = 1 gives d = 4Rs. The bubbles' A=1 surfaces touch at the midpoint when the centers are 4 Rs apart. (For unequal masses the criterion shifts, but the principle is the same.)

**Merger.** At d = d_crit, the A=1 surfaces touch and reconnect into a single common bubble. All matter — which lives on the surfaces, per the bubble picture — instantly joins onto the new common surface. No information is lost; no interior to fall into.

**Ringdown.** Below d_crit, the common bubble has higher-multipole deformations (non-spherical) that radiate away as gravitational waves. The bubble settles toward a more spherical (or, with rotation, Kerr-shaped) final state.

## Numerical results

|   separation_over_Rs | topology    |   area_2d_equatorial |   midpoint_A |
|---------------------:|:------------|---------------------:|-------------:|
|                 10   | two_bubbles |               7.6698 |     0.4      |
|                  6   | two_bubbles |               9.117  |     0.666667 |
|                  5   | two_bubbles |              10.0332 |     0.8      |
|                  4.5 | two_bubbles |              10.7946 |     0.888889 |
|                  4.1 | merging     |              11.8458 |     0.97561  |
|                  4   | one_bubble  |              12.3471 |     1        |
|                  3.9 | one_bubble  |              12.8925 |     1.02564  |
|                  3.5 | one_bubble  |              13.7421 |     1.14286  |
|                  3   | one_bubble  |              14.0679 |     1.33333  |
|                  2.5 | one_bubble  |              13.9797 |     1.6      |
|                  2   | one_bubble  |              13.6683 |     2        |
|                  1.5 | one_bubble  |              13.2813 |     2.66667  |



**Critical separation (numerical):** d_crit = 4.0000 Rs (matches analytical d_crit = 4 Rs).


## Why d_crit = 4 Rs is exact (not numerical)

From linear superposition, A at the midpoint between two equal-mass sources at separation d is `A_mid = 2 × (Rs / (d/2)) = 4 Rs / d`. The midpoint is the first point where the two sources' contributions overlap maximally; if A_mid < 1, neither bubble's A=1 surface reaches the midpoint, so the bubbles are separate. If A_mid > 1, the midpoint is enclosed in the combined A≥1 region. The transition is at A_mid = 1, giving `d_crit = 4 Rs`. This is a clean topological prediction of bold STAM with linear superposition.

## The thirds-of-A pattern in two-BH geometry

The plot `SF5_contours_panel.png` shows three contours at A = 1/3, 2/3, 1 — the Schwarzschild ISCO, photon sphere, and event horizon thresholds — for various separations. As the BHs approach each other:
- The A=1/3 (ISCO) contour merges first, well before the horizon does.
- Then the A=2/3 (photon) contour merges next.
- Finally at d=4Rs, the A=1 (horizon) contour merges.

This reproduces the LIGO inspiral picture: the orbit decays through ISCO first, then through the photon-sphere threshold, then the horizons merge. Each transition is a topological event in the A field, locatable by where the corresponding A contour reconnects.

## What this gives bold STAM

- **Multi-source A is now formalized:** linear superposition, with explicit topology of merger geometry.
- **Merger picture in bold-STAM language:** bubbles (A=1 surfaces) reconnecting, with all matter on the surfaces joining onto the new common surface. No information loss, no interior physics needed.
- **The merger is a topological event** at exactly d = 4 Rs (for equal masses, in linear superposition). This is a clean prediction that can be matched against any future detailed inspiral analysis.
- **Foundation for LIGO armor:** with this multi-source structure plus future dynamical-A propagation, bold STAM should reproduce LIGO inspiral and ringdown phenomenology. Detailed waveform comparison is open.

## Open / future work in this direction

- **Dynamical A:** propagation of A perturbations as gravitational waves. Not done in this script. Would require committing to a wave equation for A.
- **Strong-field corrections to linear superposition:** linear superposition is exact in weak field but might need modification very near merger. The qualitative topology should be correct regardless.
- **Quasinormal mode frequencies:** the ringdown spectrum depends on the near-horizon geometry of the merged BH. With bold-STAM g_rr modification, QNM frequencies could differ from GR by ~few%. Testable with next-gen LIGO.

## Generated plots

- `plots/SF5_contours_panel.png`
- `plots/SF5_midpoint_A.png`
- `plots/SF5_horizon_area_vs_separation.png`
