#' Jaric & Ebenhard (2010) sighting-trend index
#'
#' Tests persistence from the average interval between distinct sighting times,
#' optionally adjusted by the mean trend in consecutive interval lengths. The
#' trend-adjusted calculation implements equations 5-6 of Jaric & Ebenhard
#' (2010); `trend = FALSE` implements equations 3-4.
#'
#' @inheritParams solow1993
#' @param trend Logical. If `TRUE` (the default), adjust the average interval
#'   by the mean change in consecutive interval lengths.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`. The latter contains index p-values; the
#'   column name is retained for compatibility with EDE 0.1.0.
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
  if (missing(test_year) || !is.numeric(test_year) ||
      length(test_year) != 1L || is.na(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }
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
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  intervals <- diff(times)
  average_interval <- mean(intervals)
  trend_coefficient <- if (trend) mean(diff(intervals)) else 0
  adjusted_interval <- average_interval + trend_coefficient
  if (!is.finite(adjusted_interval) || adjusted_interval <= 0) {
    stop("the average interval plus its trend must be positive.", call. = FALSE)
  }

  candidates <- seq(last_sight + 1, test_year)
  p_value <- adjusted_interval /
    (adjusted_interval + candidates - last_sight)

  if (data_out) {
    out <- data.frame(time = candidates, chance = p_value)
    attr(out, "average_interval") <- average_interval
    attr(out, "trend_coefficient") <- trend_coefficient
    return(out)
  }

  below <- candidates[p_value <= alpha]
  if (length(below) == 0L) {
    warning("p-value never falls to alpha before `test_year`; returning NA.",
            call. = FALSE)
    return(new_ede_estimate(NA_real_, method = "Jaric & Ebenhard (2010)",
                            alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Jaric & Ebenhard (2010)", alpha = alpha)
}
