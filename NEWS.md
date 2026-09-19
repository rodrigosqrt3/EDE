# EDE 0.2.0 (development)

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
