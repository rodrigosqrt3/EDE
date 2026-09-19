#' Solow (2005) Weibull extreme-value persistence test
#'
#' Tests persistence using the Weibull extreme-value model for the `k` most
#' recent sighting events. This is the hypothesis-test counterpart of
#' [ole()]. For a candidate time `T`, the p-value is equation 16 of Solow
#' (2005). Earlier EDE versions incorrectly described a Fisher-gap calculation
#' as a sighting-effort-weighted method from this paper.
#'
#' @inheritParams solow1993
#' @param k Number of most recent sighting events to use. Counts greater than
#'   one are treated as independent events at the recorded time. By default,
#'   all sighting events are used. Must be at least 3.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @export
solow2005 <- function(sd, alpha = 0.05, test_year, data_out = FALSE, k = NULL) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year) || length(test_year) != 1L || is.na(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

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
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  candidates <- seq(last_sight + 1, test_year)
  shape <- ole_shape(times)
  oldest <- times[k]
  chance <- exp(-k * ((candidates - last_sight) /
                        (candidates - oldest))^(1 / shape))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "p-value never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Solow (2005) Weibull test", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Solow (2005) Weibull test", alpha = alpha)
}
