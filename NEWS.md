# EDE 0.3.0

A consistency review of every method against its source paper.

## New features

Every new method reproduces the published examples of its source paper;
where a paper's own numbers are inconsistent, the help page says so.

### Bayesian models

- `solow1993_bayes()`: Bayes factor and posterior probability of persistence
  of Solow (1993a), equation 3.
- `rout2009()` and `rout2009_absences()`: Bayes factor under a declining
  sighting rate and the number of absent surveys needed to declare
  eradication (Rout, Salomon & McCarthy 2009).
- `caley_barry2014()`: probability of persistence from yearly detections
  (Caley & Barry 2014). The constant model is computed exactly; the
  non-constant model ports the authors' MCMC sampler.
- `strauss1989_bayes()`: Bayesian posterior for the end of a range (Strauss &
  Sadler 1989, equations 22-27).

### Uncertain sightings and surveys

- `thompson2013()`: probability of persistence with several observation types
  of uncertain detection probability (Thompson et al. 2013), with the
  annealed and quenched averages computed exactly.
- `lee2014()`: probability of persistence from certain and uncertain sightings
  with false positives (Lee et al. 2014), following the authors' OpenBUGS
  model but computed exactly.
- `kodikara2020()`: posterior of the extinction year from yearly certain and
  uncertain sightings, Models 1-3 of Kodikara et al. (2020), computed exactly
  instead of by MCMC.
- `kodikara2021()`: extinction time when sightings follow non-homogeneous
  Poisson processes with Weibull rates and a change point (Kodikara et al.
  2021), by numerical integration. The default reproduces the paper with the
  authors' likelihood; equation 9 as printed is available as an option.
- `brook2019()`: resampling with sighting reliabilities as inclusion
  probabilities, around any estimator (Brook, Buettel & Jarić 2019).
- `thompson2017()`: the iterative model of records and surveys behind the
  IUCN Red List extinction-probability tool (Thompson et al. 2017).
- `thompson2019()`: Bayesian year-by-year updating (Thompson et al. 2019).
- `jaric_roberts2014()` gains `test_year` and `data_out` for the
  reliability-adjusted p-value of equations 7-8.

### Stratigraphic ranges

- `marshall1994()` and `marshall1994_trend()`: distribution-free confidence
  intervals on a range end-point and the test for a trend in gap size.
- `marshall1997()`: confidence intervals under a non-uniform fossil recovery
  potential.
- `weiss_marshall1999()` and `weiss_marshall1999_prior()`: Bayesian end-point
  for sections sampled at discrete horizons.

### Extensions of existing methods

- `solow1993b()` returns the minimum variance unbiased estimate of the
  extinction time (Solow 2005, equation 8) in `estimate`.
- `burgman1995(equation = 3)`: the runs test on occupied and empty periods.

### Tools and data

- `ede_compare()` applies every method that needs only the sighting record
  and returns one table.
- Example records from the literature: `dodo`, `monk_seal`, `ferret`,
  `orchid` and `extinct_birds`.
- Bayesian intervals print as credible intervals.

## Breaking changes

- The persistence tests (`solow1993()`, `solow1993b()`, `solow2005()`,
  `mcinerny2006()`, `burgman1995()`, `solow_roberts2003()` and `jaric2010()`)
  no longer report the year of first rejection as a point estimate, which
  none of the source papers define. They return `estimate = NA` and the
  one-sided upper confidence bound from inverting the test in `upper`. For
  continuous-time methods the bound is now exact rather than rounded up to a
  whole year: equation 5 of Solow (2005) for `solow1993()`, equation 17 of
  Solow (2005) for `solow2005()`, equations 7-8 of Jarić & Ebenhard (2010)
  for `jaric2010()`, the Robson-Whitlock bound for `solow_roberts2003()`, and
  a numerical root for `solow1993b()`.
- `burgman1995()` now applies equation 2 of Burgman et al. (1995), the test
  for extinction, by default. The former calculation, equation 4, tests for a
  change in conservation status using the longest empty run anywhere in the
  record; it remains available as `equation = 4`.
- `solow1993()` and `solow1993b()` count every sighting event, as in Solow
  (1993a), instead of distinct sighting times (`solow1993()`) or dropping all
  sightings in the origin time unit (both). Only the sighting that defines the
  origin is omitted.
- `solow1993()` is now labelled Solow (1993a).

## Bug fixes

- `ole()` returned an upper confidence limit earlier than the last sighting
  when `-log(alpha / 2) >= k`. The limit does not exist in that case and is
  now returned as `Inf` with a warning.
- `mcinerny2006()` ignored a zero-count first row. It now uses it as the known
  start date and counts every sighting time, as the source paper allows.
- `burgman1995(equation = 4)` returned `NaN` for large numbers of sightings.
  The inclusion-exclusion sum is now computed on the log scale.
- Candidate times no longer extend past `test_year` when it is less than one
  time unit after the last sighting.

## Documentation

- Corrected the equations claimed by `jaric_roberts2014()` and stated that
  its binary-count requirement is a package restriction.
- Documented how every method treats repeated sightings and the observation
  origin, the reduced coverage of the Robson-Whitlock bound under declining
  sighting rates, the one-sided relationship between `ole()` and
  `solow2005()`, and that `jaric2010()` is an index rather than a p-value.
- Credited Cooke (1980) as the origin of optimal linear estimation.
- Added regression tests from the published examples in each paper.

# EDE 0.2.0

## Scientific corrections

- Corrected the Robson--Whitlock point estimator and separated it from its
  one-sided upper confidence bound.
- Added the unbiased Strauss--Sadler endpoint estimate instead of reporting
  its confidence bound as the point estimate.
- Corrected the OLE covariance matrix and added a configurable `k` for the
  most recent sighting events.
- Corrected conditioning on the observation origin in the Solow (1993),
  declining-rate Solow (1993), and Burgman et al. (1995) procedures.
- Replaced the unsupported former `solow2005()` calculation with the Weibull
  extreme-value p-value from equation 16 of Solow (2005).
- Clarified that frequentist outputs are p-values rather than posterior
  probabilities of persistence.
- Added regression tests based on published examples for the dodo, Caribbean
  monk seal, and black-footed ferret.

## New methods

- Added `solow_roberts2003()` for the nonparametric two-endpoint persistence
  test.
- Added `jaric2010()` for average-interval and sighting-trend inference.
- Added `jaric_roberts2014()` for records with individual sighting
  reliabilities.
