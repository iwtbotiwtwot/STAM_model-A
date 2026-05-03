# data/

Place input datasets here.

Recommended layout:

```text
data/
├── pantheon/
├── union3/
├── des/
├── bao/
└── clocks/
```

Do not commit large or license-restricted datasets unless their license allows it. If datasets are not committed, add exact download instructions and checksums.

For supernova template scripts, expected columns are one of:

```text
z, D_obs
```

or:

```text
z, mu
```

where `mu` is distance modulus. If `mu` is used, scripts convert to luminosity distance in Mpc and then to Mly using:

```text
D_Mpc = 10 ** ((mu - 25)/5)
D_Mly = D_Mpc * 3.261563776
```
