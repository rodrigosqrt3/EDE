#' Solow (1993) constant-rate persistence test
#'
#' Parametric test of the null hypothesis that a species was still extant at a
#' candidate test year, under a stationary Poisson sighting process.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). Persistence is rejected for
#'   the first candidate year at which the p-value falls to
#'   or below `alpha`.
#' @param test_year Latest year to test. Must be supplied as a number.
#' @param data_out If `TRUE`, return the full p-value curve
#'   instead of the single first-rejection year.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Solow, A. R. (1993). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' @export
solow1993 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  full <- expand_record(sd, test_year)

  times <- sort(unique(sd$time[sd$count > 0]))
  origin <- min(sd$time)
  n_used <- sum(times > origin)
  if (n_used < 1) {
    stop("Solow (1993) needs at least 2 distinct sighting times when the first sighting defines the observation origin.", call. = FALSE)
  }

  t1 <- times[1]
  last_sight <- times[length(times)]
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  candidates <- full$time[full$time > last_sight]
  chance <- vapply(candidates, function(t) {
    solow1993_chance(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "p-value never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Solow (1993)", alpha = alpha))
  }

  new_ede_estimate(below[1], method = "Solow (1993)", alpha = alpha)
}

#' @keywords internal
#' @noRd
solow1993_chance <- function(d) {
  origin <- min(d$time)
  event_times <- sort(unique(d$time[d$count > 0 & d$time > origin]))
  if (length(event_times) == 0L) return(1.0)
  tn <- max(event_times) - origin
  tmax <- max(d$time) - origin
  n <- length(event_times)
  (tn / tmax)^n
}
