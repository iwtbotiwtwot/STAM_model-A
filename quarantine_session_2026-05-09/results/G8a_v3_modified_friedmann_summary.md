# G8a: Model-A V_3 Modified Friedmann

## Setup

V_3 = α/A + β/(1−A), selected by G8 under the symmetric-boundary commitment.
Constants:

```text
A_0 = 1/(12π) = 0.02652582                 (G7 commitment)
α/β = [A_0/(1−A_0)]² = 7.4249e-4            (V_3 minimum at A_0)
Ω_m = 0.315, Ω_r = 9.2e-5                   (Model-A internal content)
```

β_tilde is calibrated by Model-A's own internal closure h²(today) = 1
(not anchored to any external Ω_DE target). Tracking ansatz:
A_eq(a) is the instantaneous equilibrium V'(A_eq) = κρ_m(a), and
ω_tracking = −3h × Y'(A_eq)/Y''(A_eq).

## Results — tracking-ansatz cosmology

Model-A V_3 tracking table at H_0 = 73.04 km/s/Mpc:

```text
β_tilde          = 0.647255   (Model-A internal closure)
α_tilde          = 4.806e-4   (= β_tilde × α/β)
A_eq(today)      = 0.035524    (slightly above A_0 because matter is non-zero)
ω(today)         = -0.0413     (small — field near its minimum)
Ω_DE_potential   = 0.6846      (V_3(A_eq) at a=1)
Ω_DE_kinetic     = 0.0003      (ω²/6)
Ω_DE_total       = 0.6849
h²(today)        = 1.000000    (exact, by construction)
```

w_eff(z) (effective dark-energy equation of state):

```text
w_eff(z=0) = -0.991   (essentially cosmological-constant today)
w_eff(z=1) = -0.684   (quintessence-like at intermediate z)
w_eff(z=3) = -0.544   (approaches slow-roll value)
```

H(z) tabulation:

```text
z = 0.0:  H = 73.04 km/s/Mpc
z = 0.1:  H = 76.86
z = 0.5:  H = 113.17
z = 1.0:  H = 149.75
z = 2.0:  H = 244.20
z = 3.0:  H = 358.47
z = 5.0:  H = 632.81
```

## Tracking-vs-numerical caveat (V_3 specific)

The full Klein-Gordon + Friedmann numerical integration disagrees
substantially with the tracking ansatz at late times:

```text
                       tracking      numerical KG
A(today)               0.0355        0.2134
ω(today)               -0.041        -0.754
h(today)               1.0000        1.1113
```

**Why:** V_3 has a true minimum at A_0 (V'(A_0) = 0). A field with inertia
will overshoot the minimum and oscillate around it under Hubble damping.
The slow-roll tracking ansatz ω = −3h × Y'/Y'' assumes the field stays at
the instantaneous matter-driven equilibrium, which V_3 does not honor in
full dynamics. V_1 had the same issue at smaller magnitude, because V_1
has no minimum (V'(A) > 0 everywhere) so the tracking attractor is the
natural late-time solution.

**For the SN distance test (G8b):** uses the tracking h²(a) for direct
methodological comparison with the V_1 SN test (script 37), which also
used tracking. The numerical-KG cosmology of V_3 is a separate study.

## Files

- `scripts/G8a_v3_modified_friedmann.py`
- `reports/G8a/v3_modified_friedmann_summary.png`
- `reports/G8a/v3_tracking_attractor_check.png`
- `reports/G8a/v3_today_budget.csv`
