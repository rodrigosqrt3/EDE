#' Jaric & Ebenhard (2010) sighting-trend index
#'
#' Index of persistence from the average interval `I` between distinct
#' sighting times, optionally adjusted by the trend coefficient `c`. For
#' candidate time `T`, the index is `I / (I + T - t[n])` (equations 3-4 of
#' Jaric & Ebenhard, 2010) or, with `trend = TRUE`,
#' `(I + c) / (I + c + T - t[n])` (equations 5-6). The paper calls this the
#' probability that the species is still present; it is an index rather than
#' a p-value from a fitted probability model.
#'
#' Several sightings in one time unit count as one sighting, as in the paper.
#' `I` is the mean interval between the first and last sighting times, so the
#' observation origin is not used. `c` is the mean change in consecutive
#' interval lengths, which reduces to the difference between the last and
#' first intervals divided by `n - 2`; it is positive when sightings become
#' less frequent.
#'
#' The returned upper bound is `t[n] + (I + c) (1 - alpha) / alpha`, the upper
#' bound of the confidence interval in equations 7-8 of the paper (`c = 0`
#' when `trend = FALSE`).
#'
#' @inheritParams solow1993
#' @param trend Logical. If `TRUE` (the default), adjust the average interval
#'   by the trend coefficient (equations 5-6 and 8).
#'
#' @return An [ede_estimate] object with `estimate = NA`, the one-sided upper
#'   bound in `upper`, and additional components `average_interval` and
#'   `trend_coefficient`; or (if `data_out = TRUE`) a data frame with columns
#'   `time` and `chance` (index values) carrying the same two quantities as
#'   attributes.
#'
#' @references
#' Jaric, I., & Ebenhard, T. (2010). A method for inferring extinction based
#'   on sighting records that change in frequency over time. Wildlife Biology,
#'   16(3), 267-275.
#'
#' @export
jaric2010 <- function(sd, alpha = 0.05, test_year, data_out = FALSE,
                      trend = TRUE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  check_test_year(test_year)
  if (length(trend) != 1L || !is.logical(trend) || is.na(trend)) {
    stop("`trend` must be TRUE or FALSE.", call. = FALSE)
  }

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  minimum_n <- if (trend) 3L else 2L
  if (n < minimum_n) {
    stop(sprintf("Jaric & Ebenhard (2010) needs at least %d distinct sighting times when `trend = %s`.",
                 minimum_n, trend), call. = FALSE)
  }

  last_sight <- times[n]
  check_test_year(test_year, last_sight)

  intervals <- diff(times)
  average_interval <- mean(intervals)
  trend_coefficient <- if (trend) mean(diff(intervals)) else 0
  adjusted_interval <- average_interval + trend_coefficient
  if (!is.finite(adjusted_interval) || adjusted_interval <= 0) {
    stop("the average interval plus its trend must be positive.", call. = FALSE)
  }

  candidates <- candidate_times(last_sight, test_year)
  p_value <- adjusted_interval /
    (adjusted_interval + candidates - last_sight)

  if (data_out) {
    out <- data.frame(time = candidates, chance = p_value)
    attr(out, "average_interval") <- average_interval
    attr(out, "trend_coefficient") <- trend_coefficient
    return(out)
  }

  bound <- last_sight + adjusted_interval * (1 - alpha) / alpha
  result <- test_bound_result(bound, test_year, "Jaric & Ebenhard (2010)", alpha)
  result$average_interval <- average_interval
  result$trend_coefficient <- trend_coefficient
  result
}
