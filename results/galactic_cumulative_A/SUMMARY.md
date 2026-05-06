# STAM Model-A galactic cumulative A scratch test

## Question

Can individually tiny nonzero `A` contributions accumulate into a galaxy-scale field?

## Linear cumulative Model-A

```text
A_total(x)=Σ 2Gm_i/(c²|x-x_i|)
```

At 10 kpc:

```text
A from one solar mass      ≈ 9.571121e-18
A from 1e11 solar masses   ≈ 9.571121e-07
compact circular speed     ≈ 207.390 km/s
```

So many individually tiny `A` contributions can sum into a galaxy-scale quantity.

## Toy visible galaxy result

Toy components:

```text
disk  = 6.00e+10 Msun
gas   = 1.00e+10 Msun
bulge = 1.00e+10 Msun
```

Visible linear accumulation gives:

```text
outer 15-35 kpc median speed ≈ 121.076 km/s
outer slope                  ≈ -2.822 km/s/kpc
```

## Collective-envelope exploratory toy

If a target edge speed of `220.0 km/s` is imposed at `20.0 kpc`, the required envelope parameter is:

```text
epsilon ≈ 5.775256e-07
```

This is a small dimensionless `A`-gradient strength.

The collective-envelope result gives:

```text
outer 15-35 kpc median speed ≈ 198.689 km/s
outer slope                  ≈ -1.708 km/s/kpc
```

## Interpretation

The linear test supports the idea that many tiny `A` contributions can accumulate into a galaxy-scale field.

The "greater than the linear sum" idea requires an additional collective/coherence/envelope term:

```text
A_total = A_linear + A_collective
```

That envelope part is exploratory only.
