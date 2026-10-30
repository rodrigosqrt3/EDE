# EDE

[![CRAN status](https://www.r-pkg.org/badges/version/EDE)](https://CRAN.R-project.org/package=EDE) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/EDE/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/EDE/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/EDE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/EDE)

Extinction date estimation from sighting records.

Given a time-ordered record of sightings, EDE estimates when a species
most likely went extinct, or how likely it is to still be extant. It
collects methods from the sighting-record and stratigraphic-range
literature under one interface: endpoint estimators, frequentist persistence
tests, Bayesian models, methods for sightings of uncertain reliability,
models that combine records with surveys (including the model behind the
IUCN Red List extinction-probability tool), and confidence intervals on the
end-points of fossil ranges. Every method is tested against the published
examples of its source paper.

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
Methods differ in their treatment of multiple sightings in one time unit,
following their source papers: OLE, the Weibull test, both Solow (1993)
tests and Burgman's tests count every sighting event, while the
Robson-Whitlock, Strauss-Sadler, McInerny, Solow-Roberts and Jarić-Ebenhard
methods use only the distinct sighting times.

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
#> <Solow (1993a)>
#>   95% upper confidence bound: 1939.32

solow1993b(sd, test_year = 2000)
#> <Solow (1993b)>
#>   95% upper confidence bound: NA
#> Warning message:
#> p-value never falls to alpha before `test_year`; returning NA.

solow2005(sd, test_year = 2000)
#> <Solow (2005) Weibull test>
#>   95% upper confidence bound: 1953.99

mcinerny2006(sd, test_year = 2000)
#> <McInerny, Roberts, Davy & Cribb (2006)>
#>   95% upper confidence bound: 1942

burgman1995(sd, test_year = 2000)
#> <Burgman, Grimson & Ferson (1995), equation 2>
#>   95% upper confidence bound: 1940

solow_roberts2003(sd, test_year = 2000)
#> <Solow & Roberts (2003)>
#>   95% upper confidence bound: 1993

jaric2010(sd, test_year = 2000)
#> <Jaric & Ebenhard (2010)>
#>   95% upper confidence bound: 1980.02
```

`solow1993()`, `solow1993b()`, `solow2005()`, `mcinerny2006()`,
`burgman1995()`, `solow_roberts2003()` and `jaric2010()` are persistence
tests. Their source papers define a test, not a point estimate, so they
return `estimate = NA` and the one-sided `1 - alpha` upper confidence bound
obtained by inverting the test: the earliest time at which persistence is
rejected. For the continuous-time methods this bound is computed exactly
(for example, equation 17 of Solow 2005); for the discrete-period methods of
McInerny et al. and Burgman et al. it is the first rejected time unit. A
bound beyond `test_year`, or one that does not exist, is returned as `NA`
with a warning. Pass `data_out = TRUE` to get the p-value curve at integer
steps up to `test_year`. For compatibility with version 0.1.0, that p-value
is stored in a column currently named `chance`.

## Comparing methods

`ede_compare()` applies every method that needs only the sighting record and
returns one table. Methods that cannot be applied stay in the table with the
reason in `note`.

```r
ede_compare(sighting_data(monk_seal), test_year = 2050)
```

## Example data

EDE bundles classic records from the literature, each reproducing its
published example: `dodo` (Roberts & Solow 2003), `monk_seal` (Solow 1993a),
`ferret` (Solow 1993b; Jarić & Ebenhard 2010), `orchid` (Solow & Roberts
2003) and `extinct_birds`, four possibly extinct birds with sightings
classified by evidence quality (Roberts et al. 2010; Jarić & Roberts 2014).

```r
ole(sighting_data(dodo))
#> <OLE (Roberts & Solow 2003)>
#>   estimate: 1690.42
#>   95% CI: [1669.13, 1798.88]
```

## Sighting record and persistence curves

```r
curve_1993 <- solow1993(sd, test_year = 2000, data_out = TRUE)
curve_2005 <- solow2005(sd, test_year = 2000, data_out = TRUE)

plot(curve_1993$time, curve_1993$chance, type = "l",
     xlab = "candidate extinction year", ylab = "p-value")
lines(curve_2005$time, curve_2005$chance, col = "firebrick")
```

The two curves diverge because they encode different models. `solow1993()`
assumes a stationary Poisson process and uses the number of sighting events
after the observation origin. `solow2005()` instead uses the spacing of the
most recent sightings under a Weibull extreme-value model, the same model
used by `ole()`.

For uncertain records, supply one reliability per sighting (one per occupied
row, or one per row):

```r
uncertain <- sighting_data(data.frame(
  year = c(1900, 1910, 1920, 1930),
  sightings = c(1, 1, 1, 1)
))

jaric_roberts2014(uncertain, reliability = c(1, 0.9, 0.6, 0.3))
```

Pass `test_year` to `jaric_roberts2014()` to also get its reliability-adjusted
p-value, and `data_out = TRUE` for the p-value curve.

For a Bayesian view of the constant-rate model, `solow1993_bayes()` gives the
Bayes factor in favour of persistence and the posterior probability of
persistence for a chosen prior:

```r
solow1993_bayes(sd, test_year = 2000, prior = 0.5)
```

## Methods

**Endpoint estimators**

| Function | Reference | Output |
|---|---|---|
| `ole()` | Roberts & Solow (2003) | point estimate + two-sided CI |
| `robson1964()` | Robson & Whitlock (1964) | point estimate + one-sided CI |
| `strauss1989()`, `strauss1989_curve()` | Strauss & Sadler (1989) | point estimate + one-sided CI / curve |

**Persistence tests**

| Function | Reference | Output |
|---|---|---|
| `solow1993()` | Solow (1993a) | upper bound or p-value curve |
| `solow1993b()` | Solow (1993b); Solow (2005) | unbiased estimate + upper bound, or p-value curve |
| `solow2005()` | Solow (2005) | upper bound or p-value curve |
| `mcinerny2006()` | McInerny, Roberts, Davy & Cribb (2006) | upper bound or p-value curve |
| `burgman1995()` | Burgman, Grimson & Ferson (1995) | upper bound or p-value curve (extinction, eq. 2; change in status, eqs. 3-4) |
| `solow_roberts2003()` | Solow & Roberts (2003) | upper bound or p-value curve |
| `jaric2010()` | Jarić & Ebenhard (2010) | upper bound or index curve |

**Bayesian models**

| Function | Reference | Output |
|---|---|---|
| `solow1993_bayes()` | Solow (1993a) | Bayes factor + posterior probability of persistence |
| `rout2009()`, `rout2009_absences()` | Rout, Salomon & McCarthy (2009) | Bayes factor under a declining sighting rate; surveys needed to declare eradication |
| `caley_barry2014()` | Caley & Barry (2014) | posterior probability of persistence from yearly detections (exact constant model, or MCMC) |
| `strauss1989_bayes()` | Strauss & Sadler (1989) | posterior mean + credible interval of the end-point |

**Uncertain sightings**

| Function | Reference | Output |
|---|---|---|
| `jaric_roberts2014()` | Jarić & Roberts (2014) | reliability-adjusted estimate + upper bound, p-value |
| `thompson2013()` | Thompson et al. (2013) | probability of persistence with several observation types |
| `brook2019()` | Brook, Buettel & Jarić (2019) | resampling by reliability around any estimator |

**Records and surveys**

| Function | Reference | Output |
|---|---|---|
| `thompson2017()` | Thompson et al. (2017) | year-by-year probability extant (IUCN tool), with uncertainty bands |
| `thompson2019()` | Thompson et al. (2019) | Bayesian year-by-year updating |

**Stratigraphic ranges**

| Function | Reference | Output |
|---|---|---|
| `marshall1994()`, `marshall1994_trend()` | Marshall (1994) | distribution-free bounds on the CI; test for a trend in gap size |
| `marshall1997()` | Marshall (1997) | one-sided CI under non-uniform recovery potential |
| `weiss_marshall1999()` | Weiss & Marshall (1999) | Bayesian end-point for discretely sampled sections |

For fossil data, the end-point estimated is always the one at the largest
`time`: code ages as `time = -age`, and first appearances as negative
positions.

See `vignette("EDE")` for the statistical background and guidance on
choosing a method.

## License

GPL (>= 3)
