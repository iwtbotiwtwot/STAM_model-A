# STAM Model-A no-b supernova distance run

Primary physical Model-A no-b prediction:

```text
D_adj,0(z) = L z (1 + 0.5 z)
```

Comparison geometric spine:

```text
D_geo(z) = L z (1 + 0.15 z)
```

Common range:

```text
0.05 <= z <= 1.14418
```

## No-b Model-A residuals: observed - predicted

| catalog | median mag | mean mag | RMSE mag | median obs/model distance ratio |
|---|---:|---:|---:|---:|
| Union3 | 0.0512 | 0.0429 | 0.0539 | 1.0238 |
| Pantheon | 0.0808 | 0.0800 | 0.1664 | 1.0379 |
| DES | 0.2068 | 0.2383 | 0.3626 | 1.0999 |

## Geometric spine residuals: observed - D_geo/SU

| catalog | median mag | mean mag | RMSE mag | median obs/model distance ratio |
|---|---:|---:|---:|---:|
| Union3 | 0.4119 | 0.3657 | 0.3987 | 1.2089 |
| Pantheon | 0.2815 | 0.2951 | 0.3497 | 1.1384 |
| DES | 0.5344 | 0.5615 | 0.6320 | 1.2791 |

## b-like bridge implied by no-b Model-A

| catalog | median b-like bridge [Mly] | mean b-like bridge [Mly] | std [Mly] |
|---|---:|---:|---:|
| Union3 | 387.94 | 326.77 | 269.26 |
| Pantheon | 563.54 | 620.91 | 1103.35 |
| DES | 1642.50 | 2095.93 | 2938.70 |

## Early read

The no-b Model-A curve is coherent and structured. Catalog observations remain above the no-b prediction, with Pantheon/Union3 closer to one another and DES still elevated relative to that family.
