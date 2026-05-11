# G23: QNM damping ratio scan over metric forms

**Date:** 2026-05-11

Exploratory scan. Goal: see whether the (m=2, n=2) connection to the framework's structural '4' (per Sean) shows anything clean in the QNM landscape.

## What's scanned

k(A) = (1−A)^m · (1+A)^n for integer m ∈ [1,5], n ∈ [0,5]. h(A) = 1−A fixed. Photon sphere at A = 2/3.

## Damping ratio formula (photon sphere, A=2/3)

```
τ_alt / τ_GR = 3^((m−1)/2) · (3/5)^(n/2)
```

## Weak-field constraint

k(A) recovers GR weak-field only when **n = m − 1**. This is the 'diagonal' in the (m, n) grid.

### Weak-field-respecting family

Along n = m−1, write k = h · (1−A²)^p with p = n, m = p+1.

|   m |   n | k(A)            |   tau_ratio |
|----:|----:|:----------------|------------:|
|   1 |   0 | (1-A)^1·(1+A)^0 |      1.0000 |
|   2 |   1 | (1-A)^2·(1+A)^1 |      1.3416 |
|   3 |   2 | (1-A)^3·(1+A)^2 |      1.8000 |
|   4 |   3 | (1-A)^4·(1+A)^3 |      2.4150 |
|   5 |   4 | (1-A)^5·(1+A)^4 |      3.2400 |

Scaling: **τ_ratio = (9/5)^(p/2)**. Each (1−A²) factor multiplies by √(9/5) ≈ 1.342.

## The (m=2, n=2) connection to '4'

At (m=2, n=2): k = (1−A)²(1+A)² = (1−A²)². τ_ratio = 1.0392.

**FAILS GR weak-field recovery** (n=2, but m−1=1). Light speed at small A goes as √(1−A) instead of (1−A) — factor-of-2 error, detectable by Cassini-class tests.

**However** — Model-A's actual k = h · (1−A²)². The MODIFICATION beyond h is (1−A²)², which has m_mod=2, n_mod=2 in its own exponents. **Total modification factors: 4.**

This is the structural '4' Sean noted. The (m=2, n=2) reading applies to the modification factor, not the full k(A).

**Reading**: Model-A's modification factor (1−A²)² has exactly the same 2×2 structure as Bekenstein-Hawking's α = 2-face × 2-gravity-bridge = 4 area-per-entry. Both are '2 independent pair axes meeting at the boundary.'

## Damping ratio along the weak-field-respecting family

| p (pair²-factors in mod) | (m, n) | k(A) | τ_ratio | Comment |
|---|---|---|---|---|
| 0 | (1, 0) | (1−A) | 1.0000 | GR (no modification) |
| 1 | (2, 1) | (1−A)²(1+A) | 1.3416 | minimal pair |
| 2 | (3, 2) | (1−A)³(1+A)² = (1−A)(1−A²)² | 1.8000 | **Model-A; mod has 4 factors** |
| 3 | (4, 3) | (1−A)⁴(1+A)³ | 2.4150 | |
| 4 | (5, 4) | (1−A)⁵(1+A)⁴ | 3.2404 | |
| 5 | (6, 5) | (1−A)⁶(1+A)⁵ | 4.3478 | |

**Model-A's p=2 is the first case along the diagonal with a fully 'doubled-pair' modification (4 factors).** p=1 has 2 modification factors (one (1−A) × one (1+A)) — a single pair, not the structural 4. p=2 is the first 'two-independent-pair-axes' case.

## Full scan grid

|   m |   n | k(A)            | modification beyond h   |   mod factor count |   k at A=2/3 |   tau_ratio | WF OK   |
|----:|----:|:----------------|:------------------------|-------------------:|-------------:|------------:|:--------|
|   1 |   0 | (1-A)^1·(1+A)^0 | 1 (no mod)              |                  0 |       0.3333 |      1.0000 | True    |
|   1 |   1 | (1-A)^1·(1+A)^1 | (1-A)^0·(1+A)^1         |                  1 |       0.5556 |      0.7746 | False   |
|   1 |   2 | (1-A)^1·(1+A)^2 | (1-A)^0·(1+A)^2         |                  2 |       0.9259 |      0.6000 | False   |
|   1 |   3 | (1-A)^1·(1+A)^3 | (1-A)^0·(1+A)^3         |                  3 |       1.5432 |      0.4648 | False   |
|   1 |   4 | (1-A)^1·(1+A)^4 | (1-A)^0·(1+A)^4         |                  4 |       2.5720 |      0.3600 | False   |
|   1 |   5 | (1-A)^1·(1+A)^5 | (1-A)^0·(1+A)^5         |                  5 |       4.2867 |      0.2789 | False   |
|   2 |   0 | (1-A)^2·(1+A)^0 | (1-A)^1·(1+A)^0         |                  1 |       0.1111 |      1.7321 | False   |
|   2 |   1 | (1-A)^2·(1+A)^1 | (1-A)^1·(1+A)^1         |                  2 |       0.1852 |      1.3416 | True    |
|   2 |   2 | (1-A)^2·(1+A)^2 | (1-A)^1·(1+A)^2         |                  3 |       0.3086 |      1.0392 | False   |
|   2 |   3 | (1-A)^2·(1+A)^3 | (1-A)^1·(1+A)^3         |                  4 |       0.5144 |      0.8050 | False   |
|   2 |   4 | (1-A)^2·(1+A)^4 | (1-A)^1·(1+A)^4         |                  5 |       0.8573 |      0.6235 | False   |
|   2 |   5 | (1-A)^2·(1+A)^5 | (1-A)^1·(1+A)^5         |                  6 |       1.4289 |      0.4830 | False   |
|   3 |   0 | (1-A)^3·(1+A)^0 | (1-A)^2·(1+A)^0         |                  2 |       0.0370 |      3.0000 | False   |
|   3 |   1 | (1-A)^3·(1+A)^1 | (1-A)^2·(1+A)^1         |                  3 |       0.0617 |      2.3238 | False   |
|   3 |   2 | (1-A)^3·(1+A)^2 | (1-A)^2·(1+A)^2         |                  4 |       0.1029 |      1.8000 | True    |
|   3 |   3 | (1-A)^3·(1+A)^3 | (1-A)^2·(1+A)^3         |                  5 |       0.1715 |      1.3943 | False   |
|   3 |   4 | (1-A)^3·(1+A)^4 | (1-A)^2·(1+A)^4         |                  6 |       0.2858 |      1.0800 | False   |
|   3 |   5 | (1-A)^3·(1+A)^5 | (1-A)^2·(1+A)^5         |                  7 |       0.4763 |      0.8366 | False   |
|   4 |   0 | (1-A)^4·(1+A)^0 | (1-A)^3·(1+A)^0         |                  3 |       0.0123 |      5.1962 | False   |
|   4 |   1 | (1-A)^4·(1+A)^1 | (1-A)^3·(1+A)^1         |                  4 |       0.0206 |      4.0249 | False   |
|   4 |   2 | (1-A)^4·(1+A)^2 | (1-A)^3·(1+A)^2         |                  5 |       0.0343 |      3.1177 | False   |
|   4 |   3 | (1-A)^4·(1+A)^3 | (1-A)^3·(1+A)^3         |                  6 |       0.0572 |      2.4150 | True    |
|   4 |   4 | (1-A)^4·(1+A)^4 | (1-A)^3·(1+A)^4         |                  7 |       0.0953 |      1.8706 | False   |
|   4 |   5 | (1-A)^4·(1+A)^5 | (1-A)^3·(1+A)^5         |                  8 |       0.1588 |      1.4490 | False   |
|   5 |   0 | (1-A)^5·(1+A)^0 | (1-A)^4·(1+A)^0         |                  4 |       0.0041 |      9.0000 | False   |
|   5 |   1 | (1-A)^5·(1+A)^1 | (1-A)^4·(1+A)^1         |                  5 |       0.0069 |      6.9714 | False   |
|   5 |   2 | (1-A)^5·(1+A)^2 | (1-A)^4·(1+A)^2         |                  6 |       0.0114 |      5.4000 | False   |
|   5 |   3 | (1-A)^5·(1+A)^3 | (1-A)^4·(1+A)^3         |                  7 |       0.0191 |      4.1828 | False   |
|   5 |   4 | (1-A)^5·(1+A)^4 | (1-A)^4·(1+A)^4         |                  8 |       0.0318 |      3.2400 | True    |
|   5 |   5 | (1-A)^5·(1+A)^5 | (1-A)^4·(1+A)^5         |                  9 |       0.0529 |      2.5097 | False   |

## What this empirically surfaces

1. **The weak-field constraint n = m−1 is a sharp line in the grid.** Off-diagonal cases fail observational tests (e.g., (m=2, n=2) gives a factor-2 error in radial light speed at small A).

2. **Along the WF-respecting diagonal, the damping ratio grows monotonically with p.** Each additional (1−A²) factor in the modification multiplies the ratio by ~1.342.

3. **p=2 (Model-A) is the smallest WF-respecting case with 4 modification factors.** The next case (p=3) has 6 factors; p=1 has 2. p=2 is uniquely the '4-factor' case.

4. **If pair structure forces '4 modification factors at the boundary'** (matching Bekenstein's 2-face × 2-gravity-bridge structure), then p=2 is structurally forced, and Model-A's specific k(A) = (1−A)(1−A²)² is the natural minimal form. The 1.80 damping ratio at the photon sphere becomes a structural prediction, not a free choice.

5. **LIGO tension stays.** The 1.80 ratio is significantly above observed (GR-consistent) ringdown damping times. If the '4-factor' structural reading is right, the tension is a real concern for the framework. Either (a) the eikonal overshoots and exact Regge-Wheeler gives a smaller ratio, or (b) the structural 4-reading is wrong, or (c) Model-A's k(A) needs revision.

