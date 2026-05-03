# STAM_model-A

**STAM_model-A** is a clean research scaffold for the variable-accumulation version of the Spacetime Accumulation Model (STAM). It treats `A` as a dimensionless accumulation field, not as a universal constant.

Author: **Sean Brady**

Status: **proposed theoretical framework / testbed**. This repository is for organizing formulas, scripts, outputs, falsification tests, and versioned claims. It is not a claim of established physics.

## Core idea

STAM_model-A proposes that spacetime admits a dimensionless accumulation field:

```text
A = A(x)
```

Matter-energy sources `A`. In local spherical weak-field form:

```text
A(r) = Rs / r = 2GM / (c^2 r)
```

The model then interprets several observables as different projections of the same accumulation structure:

```text
local gravity            -> gradients of A
propagation delay        -> path integrals through A
redshift-distance excess -> cosmological path accumulation
horizon criterion        -> A = 1
mass inference           -> source strength recovered through multiple observables
```

## Core equations

Local accumulation field:

```text
A(r) = Rs / r = 2GM / (c^2 r)
```

Gravity bridge:

```text
g_vec = (c^2 / 2) grad(A)
```

For the spherical weak-field case, `grad(A) = -Rs/r^2 r_hat`, so:

```text
g_vec = -GM/r^2 r_hat
```

Propagation delay:

```text
Delta t = (1/c) int A(r) ds
```

For `A(r)=2GM/(c^2 r)`:

```text
Delta t = (2GM/c^3) int ds/r
```

Critical threshold:

```text
A = 1 <=> r = Rs
```

Cosmological distance decomposition:

```text
D_adj(z)    = D_geo(z) + D_excess(z)
D_geo(z)    = L z (1 + 0.15 z)
D_adj(z)    = L z (1 + 0.5 z) + b z
D_excess(z) = b z + 0.35 L z^2
```

Average path accumulation:

```text
<A_path>(z) = D_excess / D_geo
            = (b/L + 0.35 z) / (1 + 0.15 z)
```

Local differential path accumulation:

```text
A_path,local(z) = d(D_excess) / d(D_geo)
                = (b/L + 0.70 z) / (1 + 0.30 z)
```

Path-integral representation:

```text
D_adj(z) = int [1 + A_path,local(z)] dD_geo
```

## Current parameter values recorded for testing

These values are carried forward as historical fitted/calibration candidates. They should be locked before prediction tests and not silently retuned.

```text
Pantheon/Union-style b: 354.95
Original retained b:    461.3626922
DES-style b:            1335.412792
```

Model scale:

```text
C = 3.261563776
H = 0.000243635
L = C/H = 13387.090426252385
```

## What has been implemented here

This scaffold includes scripts and tests for:

1. Spherical local accumulation `A(r)=Rs/r`.
2. Newtonian weak-field identity through `g=(c^2/2)grad(A)`.
3. Horizon threshold identity `A=1 <=> r=Rs`.
4. Weak-field Shapiro path-integral logarithmic structure.
5. Operational mass estimator consistency checks.
6. STAM cosmological distance decomposition and path-accumulation identities.
7. Template scripts for future catalog fitting and locked-parameter falsification tests.

## What this repository does not claim yet

This repository does **not** yet claim:

- a full derivation of `A_path,local(z)` from independent physical structure;
- validated agreement with independent supernova, BAO, CMB, or lensing datasets;
- a complete dynamical theory for the `A >= 1` horizon regime;
- a replacement for general relativity;
- catalog-level results unless generated and stored under `results/` with scripts and parameters recorded.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .

python scripts/00_validate_model_a.py
pytest
```

Generate first tables and synthetic checks:

```bash
python scripts/01_generate_stam_distance_tables.py
python scripts/02_mass_inference_synthetic.py
```

Outputs are written to `results/`.

## Repo structure

```text
STAM_model-A/
├── README.md
├── CHANGELOG.md
├── PRIORITY_RECORD.md
├── pyproject.toml
├── requirements.txt
├── docs/
│   ├── FORMULAS.md
│   ├── THEORY.md
│   ├── CLAIMS_AND_STATUS.md
│   ├── FALSIFICATION_TESTS.md
│   ├── PARAMETER_LOCKING.md
│   └── archive/
├── src/stam_model_a/
│   ├── constants.py
│   ├── local.py
│   ├── propagation.py
│   ├── cosmology.py
│   └── mass_estimators.py
├── scripts/
├── tests/
├── data/
├── results/
└── notebooks/
```

## Development rule

Supportive demonstrations and falsification tests must be kept separate.

A script that proves an algebraic identity is useful, but it is not empirical validation. A script that fits `b` to one dataset is useful, but it is not prediction. A real model test must lock parameters first and then evaluate an independent holdout.
