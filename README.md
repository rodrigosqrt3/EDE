# EDE

[![CRAN status](https://www.r-pkg.org/badges/version/EDE)](https://CRAN.R-project.org/package=EDE) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/EDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rodrigosqrt3/EDE/actions/workflows/R-CMD-check.yaml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/EDE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/EDE)

Extinction date estimation from sighting records.

Given a time-ordered table of sighting counts, EDE estimates when a
species most likely went extinct, or tests whether it can still be
considered extant. Eleven procedures from the sighting-record literature are
implemented under one interface, including endpoint estimators and
constant-rate, declining-rate, sighting-rate, Weibull and combinatorial
persistence tests, a sighting-trend index, and reliability-adjusted inference.

## Installation

You can install the released version of EDE from [CRAN](https://CRAN.R-project.org) with:

```r
install.packages("EDE")
```
And the development version from [GitHub](https://github.com/rodrigosqrt3/EDE) with:

```r
# install.packages("remotes")
remotes::install_github("rodrigosqrt3/EDE")
```
## Input

Every estimator takes a `sighting_data` object: a table with one row per
time unit and the number of sightings recorded in it.

The earliest supplied time defines the observation origin. If monitoring
began before the first sighting, include an initial row with count zero.
Methods differ in their treatment of multiple sightings in one time unit:
OLE and the Weibull test treat them as independent events, Solow's 1993
constant-rate test uses occupied times, and Burgman's discrete method uses
the full frequencies.

```r
library(EDE)

years <- c(1900, 1902, 1903, 1905, 1907, 1908, 1910, 1912, 1915, 1918,
           1920, 1923, 1925, 1928, 1930, 1933, 1936)
sightings <- c(4, 3, 5, 2, 3, 4, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1)

sd <- sighting_data(data.frame(year = years, sightings = sightings))
```

The last confirmed sighting is 1936. No sightings were recorded afterward.

## Estimators

```r
ole(sd)
#> <OLE (Roberts & Solow 2003)>
#>   estimate: 1941.97
#>   95% CI: [1937.48, 1956.75]

robson1964(sd)
#> <Robson & Whitlock (1964)>
#>   estimate: 1939
#>   95% one-sided CI: [1936, 1993]

strauss1989(sd)
#> <Strauss & Sadler (1989)>
#>   estimate: 1938.25
#>   95% one-sided CI: [1936, 1943.41]

solow1993(sd, test_year = 2000)
#> <Solow (1993)>
#>   estimate: 1944

solow1993b(sd, test_year = 2000)
#> <Solow (1993b)>
#>   estimate: NA
#> Warning message:
#> p-value never falls to alpha before `test_year`.

solow2005(sd, test_year = 2000)
#> <Solow (2005) Weibull test>
#>   estimate: 1954

mcinerny2006(sd, test_year = 2000)
#> <McInerny, Roberts, Davy & Cribb (2006)>
#>   estimate: 1942

burgman1995(sd, test_year = 1945)
#> <Burgman, Grimson & Ferson (1995)>
#>   estimate: 1944

solow_roberts2003(sd, test_year = 2000)
jaric2010(sd, test_year = 2000)
```

`solow1993()`, `solow1993b()`, `solow2005()`, `mcinerny2006()`, and
`burgman1995()`, `solow_roberts2003()`, and `jaric2010()` test a grid of candidate
years up to `test_year` and return the first year at which the p-value drops
to or below `alpha`. Pass `data_out = TRUE` to get the full curve instead of
the first crossing. For compatibility with version 0.1.0, that p-value is
stored in a column currently named `chance`.

## Sighting record and persistence curves

```r
curve_1993 <- solow1993(sd, test_year = 2000, data_out = TRUE)
curve_2005 <- solow2005(sd, test_year = 2000, data_out = TRUE)

plot(curve_1993$time, curve_1993$chance, type = "l",
     xlab = "candidate extinction year", ylab = "p-value")
lines(curve_2005$time, curve_2005$chance, col = "firebrick")
```

The two curves diverge because they encode different models. `solow1993()`
assumes a stationary Poisson process and uses the number of distinct
sighting times after the observation origin. `solow2005()` instead uses the
spacing of the most recent sightings under a Weibull extreme-value model,
the same model used by `ole()`.

For uncertain records, supply one reliability per occupied time (or one per
row):

```r
uncertain <- sighting_data(data.frame(
  year = c(1900, 1910, 1920, 1930),
  sightings = c(1, 1, 1, 1)
))

jaric_roberts2014(uncertain, reliability = c(1, 0.9, 0.6, 0.3))
```

## Methods

| Function | Reference | Output |
|---|---|---|
| `ole()` | Roberts & Solow (2003) | point estimate + CI |
| `robson1964()` | Robson & Whitlock (1964) | point estimate + one-sided CI |
| `strauss1989()`, `strauss1989_curve()` | Strauss & Sadler (1989) | point estimate + one-sided CI / curve |
| `solow1993()` | Solow (1993) | first rejection time or p-value curve |
| `solow1993b()` | Solow (1993b) | first rejection time or p-value curve |
| `solow2005()` | Solow (2005) | first rejection time or p-value curve |
| `mcinerny2006()` | McInerny, Roberts, Davy & Cribb (2006) | first rejection time or p-value curve |
| `burgman1995()` | Burgman, Grimson & Ferson (1995) | first rejection time or p-value curve |
| `solow_roberts2003()` | Solow & Roberts (2003) | first rejection time or p-value curve |
| `jaric2010()` | Jarić & Ebenhard (2010) | first rejection time or p-value curve |
| `jaric_roberts2014()` | Jarić & Roberts (2014) | reliability-adjusted estimate + upper bound |

See `vignette("EDE")` for the statistical background and a worked
comparison across all eleven procedures.

## License

GPL (>= 3)
