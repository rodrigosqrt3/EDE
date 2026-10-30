#' Solow & Roberts (2003) nonparametric persistence test
#'
#' Tests the null hypothesis that a species persisted to a candidate time
#' using only the two most recent distinct sighting times. For candidate time
#' `T`, the p-value is `(t[n] - t[n-1]) / (T - t[n-1])` (equation 4 of the
#' source paper). It is exact for a constant sighting rate and asymptotically
#' correct otherwise.
#'
#' The paper derives the test from the Robson & Whitlock (1964) interval;
#' inverting it gives the one-sided upper confidence bound
#' `t[n] + (t[n] - t[n-1]) (1 - alpha) / alpha`, returned as `upper`. It is
#' the same bound as [robson1964()].
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `chance` (p-values).
#'
#' @references
#' Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for extinction
#'   based on a sighting record. Ecology, 84(5), 1329-1332.
#'
#' @export
solow_roberts2003 <- function(sd, alpha = 0.05, test_year,
                              data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  check_test_year(test_year)

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 2L) {
    stop("Solow & Roberts (2003) needs at least 2 distinct sighting times.",
         call. = FALSE)
  }
  last_sight <- times[n]
  previous_sight <- times[n - 1L]
  check_test_year(test_year, last_sight)

  gap <- last_sight - previous_sight
  candidates <- candidate_times(last_sight, test_year)
  p_value <- gap / (candidates - previous_sight)

  if (data_out) return(data.frame(time = candidates, chance = p_value))

  bound <- last_sight + gap * (1 - alpha) / alpha
  test_bound_result(bound, test_year, "Solow & Roberts (2003)", alpha)
}
