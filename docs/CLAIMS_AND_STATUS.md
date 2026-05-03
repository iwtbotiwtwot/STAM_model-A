# Claims and Status

This document separates what is currently an algebraic identity, what is a phenomenological fit, what is a hypothesis, and what remains open.

## A. Algebraic / definitional identities

These are useful checks but not empirical validation.

| Item | Status |
|---|---|
| `A(r)=Rs/r` with `Rs=2GM/c^2` | Definition in the local spherical weak-field limit |
| `g=(c^2/2)grad(A)` recovers `GM/r^2` for `A=Rs/r` | Algebraic identity |
| `A=1` gives `r=Rs` | Algebraic identity |
| `Delta t=(1/c)int A ds` gives `(2GM/c^3)int ds/r` | Algebraic substitution |
| `D_adj=D_geo+D_excess` | Definitional decomposition |
| `A_path,local=dD_excess/dD_geo` gives the stated formula | Algebraic derivative |

## B. Synthetic consistency tests

These are useful regression tests. They check that code and formulas agree.

| Item | Status |
|---|---|
| Multiple mass estimators recover a synthetic input mass | Implemented as smoke/unit tests |
| Numerical path integral reproduces analytic cosmology formula | Implemented as smoke/unit tests |
| Straight-line weak-field Shapiro integral has logarithmic/asinh structure | Implemented as smoke/unit tests |

## C. Phenomenological claims requiring data tests

| Item | Current status |
|---|---|
| STAM distance law can fit SNe data with catalog-dependent `b` | To be rerun in this repo |
| Pantheon/Union-style catalogs prefer `b ~= 355` | Historical candidate value, not yet rerun here |
| DES-style catalogs prefer `b ~= 1335` | Historical candidate value, not yet rerun here |
| Catalog-family difference is physically meaningful rather than calibration/preprocessing | Open |
| `A_path,local(z)` represents real path accumulation | Hypothesis |

## D. Independent theoretical derivations still needed

| Item | Status |
|---|---|
| Derive `A_path,local(z)` from physical structure | Open |
| Connect local `A(r)` and cosmological `A_path(z)` through one field equation | Open |
| Derive lensing normalization without hand-tuned factors | Open |
| Specify invariant horizon dynamics for `A>=1` | Open |
| Reproduce GPS/gravitational clock corrections through accumulation language | Open |

## E. Hard falsification targets

The model becomes weaker or fails if:

1. locked `b` values cannot predict independent supernova catalogs within predeclared tolerance;
2. the inferred `b` drifts systematically with redshift in a way not predicted by the model;
3. BAO/angular-diameter/CMB distance tests require incompatible geometry;
4. local acceleration, time delay, redshift, and lensing cannot be normalized by one accumulation field;
5. Solar System tests require deviations from GR/Newtonian weak-field results outside observational bounds;
6. the horizon criterion `A=1` cannot be made coordinate-invariant or dynamically meaningful.
