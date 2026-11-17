# EDE 0.4.0

An extreme-value toolkit, methods for radiometrically dated records and
several further methods from the literature, and S3 methods for the result
classes.

## New features

### Documentation

- New vignette `extreme-value-methods`, which takes over the extreme-value
  material of `vignette("EDE")` and adds `vogel2009()` and `mccarthy1998()`.
- New vignette `choosing-a-method`: the sighting model assumed by each method
  and a simulation study (15 scenario-size combinations, 1000 replicates) of
  coverage, bias and the existence of upper limits. The script is
  `data-raw/simulation.R`.

### Extreme-value estimators

- `hall_wang1999()`: the minimum-distance end-point estimator of Hall & Wang
  (1999) for the `k` most recent sightings, with its exact confidence
  interval. It does not depend on the extremal exponent, so it is defined
  more often than maximum likelihood and, unlike `ole()`, its interval does
  not rely on a plug-in shape. It reproduces the 100 m example of the paper
  (estimate 9.91, extremal exponent 0.82, 90% interval 9.80 to 9.92).
- `fraga_alves_neves2014()`: the closed-form end-point estimator of Fraga
  Alves & Neves (2014) for the Gumbel domain of attraction, which is never
  earlier than the last sighting.
- `ede_select_k()`: automatic choice of the number of recent sightings `k`
  by ordered goodness-of-fit tests with ForwardStop false-discovery-rate
  control, adapting the threshold selection of Bader, Yan & Zhang (2018). It
  is a package construct, not a published method for sighting records.
- `dekkers1989()` and `ede_evi()`: the moment estimator of the extreme-value
  index (Dekkers, Einmahl & de Haan, 1989) and the end-point estimate that
  follows from it, for one `k` or over a range of `k` (a stability diagram with
  `plot()`). Unlike the Weibull-domain estimators it also says whether a finite
  end-point is plausible. It is not invariant to the origin of the time axis;
  `origin` is an argument and the help page explains the dependence.
- `smith1985()`: maximum likelihood end-point and profile-likelihood interval
  for the same model (Smith, 1985). The likelihood is regular only when the
  extremal exponent exceeds 2, which excludes a constant sighting rate; the
  function warns when the fit is outside that range, and the help page reports
  simulated coverage against `hall_wang1999()`, which is exact for any shape.
- The simulated null distribution of Greenwood's statistic, used by
  `hall_wang1999()` and `ede_endpoint_test()`, is checked against the exact
  points of Burrows (1979) quoted by Stephens (1981).
- `vogel2009()` and `ede_field_test()`: the probability plot correlation test
  of a constant sighting rate (uniform sighting times) of Vogel et al. (2009),
  and the four field-significance procedures they use to combine the tests of
  many records. The simulated null reproduces the published significance
  levels of the field tests, and the rejections of their Table 2.
- `ede_endpoint_test()`: tests whether the record is compatible with a finite
  end-point (the Weibull domain assumed by `ole()` and `hall_wang1999()`)
  against the Gumbel null, using the statistics of Neves, Picek & Fraga Alves
  (2006) and Hasofer & Wang (1992). The Neves p-value is exact (Fisher's
  distribution of the maximum over the sum of exponentials), because the
  Gumbel approximation of the paper is very conservative for small `k`. A
  third test, `test = "rosario"`, is the statistic G* of Fraga Alves, Neves &
  Rosario (2016), built on `fraga_alves_neves2014()`, which has a Gumbel limit
  under the null; its p-value under an exponential tail is simulated exactly
  as for the Hasofer & Wang test, and the asymptotic Gumbel p-value of the
  paper is given next to it. Its size is checked by simulation and its power
  against a finite end-point and a heavy tail is verified; the paper's own
  application (supercentenarian women) is not reproduced because its data are
  not published with it.
- `ede_compare()` now also reports `hall_wang1999()` and
  `fraga_alves_neves2014()`.

### Persistence tests

- `mccarthy1998()`: the effort-adjusted Solow test ("partial Solow equation")
  of McCarthy (1998), also used by Roberts & Jaric (2020), with the sampling
  effort given per time unit; it returns the upper bound or the p-value curve.

### Dated records

- `solow2006()`, `solow2006_region()` and `solow2006_check()`: inference
  about the extinction time, or the settlement time, from dated remains with
  the dating error taken into account (Solow, 1997; Solow, Roberts & Robbirt,
  2006). The true ages are uniform between a lower end-point and the
  endpoint of interest and each recorded age carries a normal error with known
  standard deviation; the endpoint is estimated by maximum likelihood with a
  profile-likelihood interval (Solow, 1997), `solow2006_region()` gives the
  joint confidence region of the extinction times of two taxa, and
  `solow2006_check()` the simulation check of the 2006 paper. The settlement
  example of Solow (1997) is reproduced (651 yr BP, interval 617 to 687, and
  the profile log-likelihood of -113.26), and so are the maximum likelihood
  estimates and the conclusions about the joint region of the mammoth and
  horse example of Solow et al. (2006), to the 0.01 ka resolution of the
  published figures. The ages can be given as a data frame with one row per
  remain, which allows tied dates with different errors.
- `bradshaw2012()`: the GRIWM method of Bradshaw et al. (2012) for dated
  records with dating error. The sighting-rate model of McInerny et al. (2006)
  is applied to the `k` most recent records for every `k`, the results are
  averaged with weights inversely proportional to the span of the series, and
  the dates are resampled from their normal errors; the estimate is the median
  and the interval the quantiles of the resampled values. The estimates for
  each `k` and the weighted estimate of the paper's Figure 2 are reproduced for
  the mammoth and horse ages, and the resampled medians agree with the
  extinction dates of the paper's Figure 4 (12,200 and 11,150 radiocarbon years
  BP, to 30 years). Without errors it returns the inverse-weighted McInerny
  estimate. The Neanderthal example is not reproduced because its dates are not
  published in the paper.

### Bayesian models

- `alroy2014()`: the iterative Bayesian method of Alroy (2014), which turns the
  frequency of sightings and the length of the observed range into a posterior
  probability of extinction for each time unit after the last sighting (the
  "creeping-shadow-of-a-doubt" values), with the best estimate of the
  extinction time (the shadow value reaches 0.5) and credible limits. The
  standard and the run-based versions are included. The worked example of the
  paper is reproduced (the posterior is 0.0183 after one unit, 0.072 after
  three, 0.506 after ten and passes 0.99 after 25), and so is its example of
  the run-based sampling probabilities. The real-data examples of the paper use
  museum records that are not published with it, and the sampling-control and
  maximum-likelihood variants are not implemented.

### Stratigraphic ranges

- `solow1996()`: the likelihood ratio test of Solow (1996) of whether several
  taxa share an upper (or lower) end-point, with an exactly known significance
  level, and the confidence interval for the common end-point. It reproduces
  the worked examples of the paper for two, three and ten ammonite taxa of
  Seymour Island, and for one taxon it equals `strauss1989()`.

### Data sets

- New data sets `alaska_megafauna` (Table 1 of Solow, Roberts & Robbirt, 2006),
  `nz_shells` (Table 1 of Solow, 1997) and `seymour_ammonites` (Table 1 of
  Solow, 1996).

### S3 methods

- `plot()` methods for `sighting_data`, `ede_estimate`, `ede_compare` (a
  forest plot) and `ede_select_k` objects.
- `as.data.frame()` method for `ede_estimate` objects.
- `ede_compare()` returns an object of class `ede_compare`, which still
  inherits from `data.frame`.

## Deprecations

- The p-value curves returned with `data_out = TRUE` by `solow1993()`,
  `solow1993b()`, `solow2005()`, `mcinerny2006()`, `burgman1995()`,
  `solow_roberts2003()`, `jaric2010()` and `jaric_roberts2014()` now have a
  column `p_value`. The former name, `chance`, is kept as an alias and will be
  removed in a future release.

## Documentation

- The `r = ...` examples in three help pages no longer trigger roxygen's
  inline-code parser.

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
