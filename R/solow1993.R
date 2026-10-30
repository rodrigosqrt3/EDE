#' Solow (1993a) constant-rate persistence test
#'
#' Parametric test of the null hypothesis that a species was still extant at a
#' candidate test year, under a stationary Poisson sighting process. For
#' candidate time `T`, the p-value is `(t[n] / T)^n` (equation 2 of Solow,
#' 1993a), with times measured from the observation origin.
#'
#' The earliest supplied time defines the observation origin. If that time
#' has a sighting, one sighting defines the origin and is omitted, as in the
#' Caribbean monk seal example of the source paper; any further sightings in
#' the same time unit are counted. Every sighting event after that is counted,
#' so a time unit with several sightings contributes several events to `n`.
#'
#' Solow (1993a) defines only the test. Inverting it gives the one-sided
#' upper confidence bound `t[n] / alpha^(1/n)` for the extinction time
#' (Solow, 2005, equation 5), which is returned as `upper`.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). The returned upper bound is the
#'   earliest candidate time at which persistence is rejected at this level.
#' @param test_year Latest time to test. Must be later than the last
#'   sighting. A bound beyond `test_year` is reported as `NA` with a warning.
#' @param data_out If `TRUE`, return the p-value curve at integer steps after
#'   the last sighting instead of the upper bound.
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `chance`. `chance` contains p-values; the column name
#'   is retained for compatibility with EDE 0.1.0.
#'
#' @references
#' Solow, A. R. (1993a). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @export
solow1993 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)

  origin <- min(sd$time)
  last_sight <- max(sd$time[sd$count > 0])
  if (last_sight <= origin) {
    stop("Solow (1993a) needs at least 1 sighting after the observation origin.", call. = FALSE)
  }
  check_test_year(test_year, last_sight)

  if (data_out) {
    full <- expand_record(sd, test_year)
    candidates <- full$time[full$time > last_sight]
    chance <- vapply(candidates, function(t) {
      solow1993_chance(full[full$time <= t, , drop = FALSE])
    }, numeric(1))
    return(data.frame(time = candidates, chance = chance))
  }

  n <- solow_event_count(sd)
  bound <- origin + (last_sight - origin) / alpha^(1 / n)
  test_bound_result(bound, test_year, "Solow (1993a)", alpha)
}

#' Solow (1993a) Bayes factor for persistence
#'
#' Bayesian counterpart of [solow1993()] under the same stationary Poisson
#' model. With a uniform prior on the extinction time over the observation
#' period and the noninformative prior `dm / m` on the sighting rate, the Bayes
#' factor in favour of persistence to `T` is
#' `B = (n - 1) / ((T / t[n])^(n - 1) - 1)` (equation 3 of Solow, 1993a),
#' which tends to `1 / log(T / t[n])` when `n = 1`. Given the prior
#' probability `prior` of persistence, the posterior probability is
#' `1 / (1 + (1 - prior) / (prior * B))`.
#'
#' Times are measured from the observation origin and sighting events are
#' counted as in [solow1993()].
#'
#' @inheritParams solow1993
#' @param test_year Candidate time `T`. Must be later than the last sighting.
#' @param prior Prior probability that the species persisted to `test_year`,
#'   in (0, 1).
#' @param data_out If `TRUE`, return the Bayes factor and posterior
#'   probability at integer steps after the last sighting up to `test_year`.
#'
#' @return A list of class `ede_bayes` with components `bayes_factor`,
#'   `posterior`, `prior`, `test_year` and `method`; or (if
#'   `data_out = TRUE`) a data frame with columns `time`, `bayes_factor` and
#'   `posterior`.
#'
#' @references
#' Solow, A. R. (1993a). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' @export
solow1993_bayes <- function(sd, test_year, prior = 0.5, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (length(prior) != 1L || !is.numeric(prior) || is.na(prior) ||
      prior <= 0 || prior >= 1) {
    stop("`prior` must be in (0, 1).", call. = FALSE)
  }
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)

  origin <- min(sd$time)
  last_sight <- max(sd$time[sd$count > 0])
  if (last_sight <= origin) {
    stop("Solow (1993a) needs at least 1 sighting after the observation origin.", call. = FALSE)
  }
  check_test_year(test_year, last_sight)

  n <- solow_event_count(sd)
  tn <- last_sight - origin
  bayes_factor <- function(t) {
    ratio <- (t - origin) / tn
    if (n == 1) 1 / log(ratio) else (n - 1) / (ratio^(n - 1) - 1)
  }
  posterior <- function(b) 1 / (1 + (1 - prior) / (prior * b))

  if (data_out) {
    times <- candidate_times(last_sight, test_year)
    b <- vapply(times, bayes_factor, numeric(1))
    return(data.frame(time = times, bayes_factor = b, posterior = posterior(b)))
  }

  b <- bayes_factor(test_year)
  structure(
    list(bayes_factor = b, posterior = posterior(b), prior = prior,
         test_year = test_year, method = "Solow (1993a) Bayes factor"),
    class = "ede_bayes"
  )
}

#' @export
print.ede_bayes <- function(x, ...) {
  cat(sprintf("<%s>\n", x$method))
  cat(sprintf("  Bayes factor for persistence to %s: %s\n",
              format(x$test_year), format(signif(x$bayes_factor, 4))))
  cat(sprintf("  posterior probability of persistence (prior %g): %s\n",
              x$prior, format(signif(x$posterior, 4))))
  invisible(x)
}

#' Number of sighting events after the observation origin
#'
#' One sighting at the origin defines it and is omitted (Solow 1993a); any
#' further sightings in the origin time unit are counted.
#' @keywords internal
#' @noRd
solow_event_count <- function(d) {
  origin <- min(d$time)
  sum(d$count) - (sum(d$count[d$time == origin]) > 0)
}

#' @keywords internal
#' @noRd
solow1993_chance <- function(d) {
  origin <- min(d$time)
  positive <- d$time[d$count > 0]
  if (length(positive) == 0L || max(positive) <= origin) return(1.0)
  tn <- max(positive) - origin
  tmax <- max(d$time) - origin
  n <- solow_event_count(d)
  (tn / tmax)^n
}
