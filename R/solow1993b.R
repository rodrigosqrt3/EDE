#' Solow (1993b) declining-population persistence test
#'
#' Nonparametric test of the null hypothesis that a declining species was still
#' extant at a candidate test year, based on the ratio of Fisher's gap
#' distribution for a non-stationary Poisson process with an exponential
#' declining rate.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1).
#' @param test_year Latest year to test. Must be later than the last sighting.
#' @param data_out If `TRUE`, return the full chance-of-persistence curve
#'   instead of the single first-rejection year.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Solow, A. R. (1993b). Inferring extinction in a declining population.
#'   Journal of Mathematical Biology, 32(1), 79-82.
#'
#' @export
solow1993b <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  full <- expand_record(sd, test_year)

  n_total <- sum(sd$count[sd$count > 0])
  if (n_total < 2) {
    stop("Solow (1993b) needs at least 2 sightings with count > 0.", call. = FALSE)
  }

  last_sight <- max(sd$time[sd$count > 0])
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  candidates <- full$time[full$time > last_sight]

  chance <- vapply(candidates, function(t) {
    solow1993b_chance(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "chance of persistence never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Solow (1993b)", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Solow (1993b)", alpha = alpha)
}

#' @keywords internal
#' @noRd
#' @keywords internal
#' @noRd
solow1993b_fisher_series <- function(u, s, n) {
  if (is.na(u) || is.nan(u) || u <= 0 || s <= 0) return(0)
  y <- u / s
  n_terms <- floor(1 / y)
  if (n_terms <= 0 || !is.finite(n_terms)) return(1)

  i <- seq_len(min(as.integer(n_terms), as.integer(n)))
  terms <- (-1)^(i - 1) * choose(n, i) * (1 - i * y)^(n - 1)
  1 - sum(terms)
}

#' @keywords internal
#' @noRd
solow1993b_chance <- function(d) {
  year <- d$time - min(d$time)
  s <- sum(year * d$count)
  tn <- max(year[d$count > 0])
  tmax <- max(year)
  n <- sum(d$count)

  if (s == 0 || n < 2) return(1.0)

  fs1 <- solow1993b_fisher_series(tn, s, n)
  fs2 <- solow1993b_fisher_series(tmax, s, n)

  p <- fs1 / fs2
  max(0, min(1, p))
}
