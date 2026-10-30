#' Solow (2005) Weibull extreme-value persistence test
#'
#' Tests persistence using the Weibull extreme-value model for the `k` most
#' recent sighting events, the model behind [ole()]. For a candidate time
#' `T`, the approximate p-value is equation 16 of Solow (2005), with the
#' shape parameter estimated by equation 15.
#'
#' Inverting the test gives the one-sided upper confidence bound of equation
#' 17, `(t[1] - c t[k]) / (1 - c)` with `c = (-log(alpha) / k)^m`, returned as
#' `upper`. The bound does not exist when `-log(alpha) >= k`. Note that
#' [ole()] reports the two-sided interval of Roberts & Solow (2003), so its
#' upper limit at level `alpha` equals this bound at level `alpha / 2`.
#'
#' @inheritParams solow1993
#' @param k Number of most recent sighting events to use. Counts greater than
#'   one are treated as independent events at the recorded time; Solow (2005)
#'   recommends combining sightings that are obviously related. By default,
#'   all sighting events are used. Must be at least 3; Solow (2005) suggests
#'   at least 5.
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `chance` (p-values).
#'
#' @references
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @export
solow2005 <- function(sd, alpha = 0.05, test_year, data_out = FALSE, k = NULL) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  check_test_year(test_year)

  times <- sort(rep(sd$time[sd$count > 0], sd$count[sd$count > 0]), decreasing = TRUE)
  n_times <- length(times)
  if (n_times < 3L) {
    stop("Solow (2005) needs at least 3 sighting events with count > 0.", call. = FALSE)
  }
  if (is.null(k)) k <- n_times
  if (length(k) != 1L || !is.numeric(k) || is.na(k) ||
      k != as.integer(k) || k < 3L || k > n_times) {
    stop("`k` must be an integer between 3 and the number of sighting events.", call. = FALSE)
  }
  k <- as.integer(k)
  times <- times[seq_len(k)]
  if (times[1] == times[2]) {
    stop("Solow (2005) requires the most recent sighting time to be strictly greater than the second most recent.", call. = FALSE)
  }

  last_sight <- times[1]
  check_test_year(test_year, last_sight)

  candidates <- candidate_times(last_sight, test_year)
  shape <- ole_shape(times)
  oldest <- times[k]
  chance <- exp(-k * ((candidates - last_sight) /
                        (candidates - oldest))^(1 / shape))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  c_alpha <- (-log(alpha) / k)^shape
  bound <- if (c_alpha < 1) (last_sight - c_alpha * oldest) / (1 - c_alpha) else NA_real_
  test_bound_result(bound, test_year, "Solow (2005) Weibull test", alpha)
}
