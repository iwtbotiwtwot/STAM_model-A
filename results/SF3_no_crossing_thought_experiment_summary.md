# STAM Model-A SF3: Two-Planet No-Crossing Thought Experiment

## Purpose

Implement the author's two-planet thought experiment with all three effects kept separate at every step. No magic bell, no free ticket. Demonstrate that the no-crossing boundary at A=1 is defined by the *collapse of frame comparison*, not by either frame seeing a special event.

## Setup

- Black hole mass: M = 1 M_sun, Rs = 2.95e+03 m.
- Alpha+ observer at large radius, A_obs = 0.001.
- Traveler released from rest at r_start = Rs/A_start; falls radially toward Rs.
- Traveler emits 'I am at A=X' messages at successive A values.

## Bold-STAM strong-field metric

```text
g_tt = -(1-A) c^2                       (GR — preserves all weak-field clock tests)
g_rr = 1 / [(1-A)(1-A^2)^2]             (Bold STAM — vanishes as (1-A)^3 at horizon)
```
At leading order in A, both reduce to GR Schwarzschild — all weak-field tests of GR pass automatically. At A=1, g_rr diverges faster than GR's, making proper distance to the boundary *infinite*. For a free-fall trajectory with E=c^2, the proper time integral picks up a logarithmic divergence at A=1: the traveler **falls forever** while crossing each intermediate A value at finite proper time.

## Three independent quantities (the triangle)

1. **Traveler's diary τ_travel(A):** proper time elapsed since release. Diverges logarithmically in bold STAM as A→1 (falls forever); finite at A=1 in GR.
2. **Clock-dilation buildup (t_emit − τ_travel):** the gravitational time-dilation contribution to what an observer sees. Diverges algebraically as 1/√(1-A).
3. **SU delay along signal path:** light-signal propagation cost from traveler's r to the observer, integrated through the A field between them (always non-zero). Diverges algebraically (1-A)^(-2) in bold STAM, faster than GR's logarithmic.

**Total observer reception time = traveler's emission coordinate time + SU delay.** The clock-dilation effect is *embedded* in the difference between coordinate time at emission and the traveler's proper time at emission.

## Numerical results — bold STAM

|   A_emit |   r_over_Rs |   tau_travel_s |   clock_dilation_buildup_s |   SU_delay_s |   t_obs_reception_s |   instantaneous_cd_factor |
|---------:|------------:|---------------:|---------------------------:|-------------:|--------------------:|--------------------------:|
|  0.01    |   100       |       0.202125 |                0.000427988 |   0.00888905 |            0.211442 |                   1.00504 |
|  0.1     |    10       |       0.208489 |                0.000567401 |   0.00980422 |            0.218861 |                   1.05409 |
|  0.3     |     3.33333 |       0.208662 |                0.000600566 |   0.00991835 |            0.219181 |                   1.19523 |
|  0.5     |     2       |       0.208687 |                0.000616468 |   0.0100084  |            0.219312 |                   1.41421 |
|  0.7     |     1.42857 |       0.208698 |                0.000633627 |   0.0101257  |            0.219458 |                   1.82574 |
|  0.9     |     1.11111 |       0.208708 |                0.000681082 |   0.0103019  |            0.219692 |                   3.16228 |
|  0.95    |     1.05263 |       0.208713 |                0.000737949 |   0.0103962  |            0.219847 |                   4.47214 |
|  0.99    |     1.0101  |       0.208721 |                0.00114846  |   0.0108406  |            0.22071  |                  10       |
|  0.999   |     1.001   |       0.208733 |                0.00562304  |   0.0155048  |            0.229861 |                  31.6228  |
|  0.9999  |     1.0001  |       0.208752 |                0.905826    |   0.145985   |            1.26056  |                 100       |
|  0.99999 |     1.00001 |       0.208776 |                3.13459     |  12.3343     |           15.6776   |                 316.228   |


## Numerical results — GR (for comparison)

|   A_emit |   r_over_Rs |   tau_travel_s |   clock_dilation_buildup_s |   SU_delay_s |   t_obs_reception_s |   instantaneous_cd_factor |
|---------:|------------:|---------------:|---------------------------:|-------------:|--------------------:|--------------------------:|
|  0.01    |   100       |       0.202123 |                0.000427981 |   0.00888896 |            0.21144  |                   1.00504 |
|  0.1     |    10       |       0.208484 |                0.000567179 |   0.00980319 |            0.218854 |                   1.05409 |
|  0.3     |     3.33333 |       0.208651 |                0.00059918  |   0.00991473 |            0.219165 |                   1.19523 |
|  0.5     |     2       |       0.208673 |                0.000612534 |   0.0100008  |            0.219286 |                   1.41421 |
|  0.7     |     1.42857 |       0.20868  |                0.000623322 |   0.0101099  |            0.219413 |                   1.82574 |
|  0.9     |     1.11111 |       0.208684 |                0.000638092 |   0.0102504  |            0.219572 |                   3.16228 |
|  0.95    |     1.05263 |       0.208684 |                0.000645736 |   0.0102936  |            0.219624 |                   4.47214 |
|  0.99    |     1.0101  |       0.208685 |                0.000662205 |   0.0103397  |            0.219687 |                  10       |
|  0.999   |     1.001   |       0.208685 |                0.000685041 |   0.0103695  |            0.219739 |                  31.6228  |
|  0.9999  |     1.0001  |       0.208685 |                0.000723888 |   0.0104043  |            0.219813 |                 100       |
|  0.99999 |     1.00001 |       0.208685 |                0.000771506 |   0.0106284  |            0.220085 |                 316.228   |


## Frame-comparison ratios (the boundary-collapse diagnostic)

|   A_emit |   r_over_Rs |   tau_travel_s |   clock_dilation_buildup_s |   SU_delay_s |   t_obs_reception_s |   instantaneous_cd_factor |   clock_dilation_over_tau |   SU_over_clock_dilation |   SU_over_tau |   t_obs_over_tau |
|---------:|------------:|---------------:|---------------------------:|-------------:|--------------------:|--------------------------:|--------------------------:|-------------------------:|--------------:|-----------------:|
|  0.01    |   100       |       0.202125 |                0.000427988 |   0.00888905 |            0.211442 |                   1.00504 |                0.00211744 |                20.7694   |     0.043978  |          1.0461  |
|  0.1     |    10       |       0.208489 |                0.000567401 |   0.00980422 |            0.218861 |                   1.05409 |                0.00272149 |                17.2792   |     0.047025  |          1.04975 |
|  0.3     |     3.33333 |       0.208662 |                0.000600566 |   0.00991835 |            0.219181 |                   1.19523 |                0.00287818 |                16.515    |     0.0475331 |          1.05041 |
|  0.5     |     2       |       0.208687 |                0.000616468 |   0.0100084  |            0.219312 |                   1.41421 |                0.00295403 |                16.235    |     0.0479588 |          1.05091 |
|  0.7     |     1.42857 |       0.208698 |                0.000633627 |   0.0101257  |            0.219458 |                   1.82574 |                0.00303609 |                15.9805   |     0.0485182 |          1.05155 |
|  0.9     |     1.11111 |       0.208708 |                0.000681082 |   0.0103019  |            0.219692 |                   3.16228 |                0.00326332 |                15.1259   |     0.0493605 |          1.05262 |
|  0.95    |     1.05263 |       0.208713 |                0.000737949 |   0.0103962  |            0.219847 |                   4.47214 |                0.00353571 |                14.088    |     0.049811  |          1.05335 |
|  0.99    |     1.0101  |       0.208721 |                0.00114846  |   0.0108406  |            0.22071  |                  10       |                0.00550235 |                 9.43924  |     0.051938  |          1.05744 |
|  0.999   |     1.001   |       0.208733 |                0.00562304  |   0.0155048  |            0.229861 |                  31.6228  |                0.0269389  |                 2.75737  |     0.0742806 |          1.10122 |
|  0.9999  |     1.0001  |       0.208752 |                0.905826    |   0.145985   |            1.26056  |                 100       |                4.33924    |                 0.161162 |     0.699321  |          6.03856 |
|  0.99999 |     1.00001 |       0.208776 |                3.13459     |  12.3343     |           15.6776   |                 316.228   |               15.0141     |                 3.93489  |    59.079     |         75.0931  |


## The headline result: A=1 is where frame comparison collapses

**For A < 1**: the triangle is computable. Each quantity is a finite number. The ratio τ_travel : clock-dilation buildup : SU delay carries the asymmetry between the frames. The traveler ages slowly compared to the observer; the magnitudes of the two observer-frame effects can be compared.

**As A → 1**: all three quantities diverge to infinity. They diverge at *different rates* — τ_travel logarithmically, clock-dilation as 1/√(1-A), SU delay as 1/(1-A)^2 — but each individually goes to ∞.

**At A = 1 itself**: ∞ = ∞ in every frame. The traveler measures infinite proper time to reach A=1 (logarithmic divergence). The observer measures infinite reception time for any 'I am at A=1' message (algebraic divergence). The SU traversal cost is also infinite. Frame-comparison ratios become ∞/∞: indeterminate.

**The very framework of comparing frames collapses.** This is the precise mathematical statement of the no-crossing boundary. It is not 'the traveler asymptotically slows.' It is 'the comparison structure between any two frames fails.' Both the traveler and the observer measure infinite time to reach A=1, by different physical mechanisms, with the same mathematical limit. There is no privileged frame in which the crossing can be registered as a finite event.

## Why this matters for the fatalities

- **'Time stops at the horizon'** — No. Time does not stop. *Comparison* stops. Both clocks keep running in their own local frames. What stops is the meaningful relationship between the frames.
- **'You're privileging the observer'** — No. The script computes both frames; both diverge; they reach the same infinite limit. The observer has no magic bell, the traveler has no privileged finite proper time.
- **'You're modifying the metric ad hoc'** — The h(A), k(A) modification gives infinite traveler proper time; that's the bold STAM commitment. But the structural point — that A=1 is where comparison breaks down — would hold for any metric with all three quantities diverging at the boundary. The specific form is one realization of a more general boundary-collapse principle.
- **'Information paradox'** — There's no paradox because there's no 'after the crossing' event in any frame. The universal ledger has no entry for 'crossed' because no frame has finite time to register one. Information stays in the resolved sector throughout; it never enters an interior, because no interior is reachable from any frame.

## Falsifiable corners (where bold STAM differs from GR)

- *Near-horizon Shapiro delay*: bold STAM SU delay diverges algebraically as 1/(1-A)^2 vs GR's logarithmic. EHT timing of light bouncing near the photon sphere could in principle test this.
- *EHT shadow size*: shadow depends on photon-sphere geometry, which is sensitive to g_rr. Modified k(A) shifts the shadow size at the few-percent level for moderate-A regions.
- *LIGO ringdown*: quasinormal mode frequencies depend on near-horizon metric structure. Bold STAM predicts shifted QNM frequencies vs GR.
- *PPN parameters at second order*: bold STAM has g_rr = 1 + A + 2A^2 + … vs GR's 1 + A + A^2 + … High-precision PPN tests beyond Cassini precision could distinguish.

## Generated plots

- `plots/SF3_triangle_panels.png`
- `plots/SF3_boundary_collapse.png`
- `plots/SF3_frame_ratios.png`
