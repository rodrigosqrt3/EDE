#' Solow & Roberts (2003) nonparametric persistence test
#'
#' Tests the null hypothesis that a species persisted to a candidate time
#' using only the two most recent distinct sighting times. For candidate time
#' `T`, the p-value is `(t[n] - t[n-1]) / (T - t[n-1])`.
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`. The latter contains p-values; the column
#'   name is retained for compatibility with EDE 0.1.0.
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
  if (missing(test_year) || !is.numeric(test_year) ||
      length(test_year) != 1L || is.na(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 2L) {
    stop("Solow & Roberts (2003) needs at least 2 distinct sighting times.",
         call. = FALSE)
  }
  last_sight <- times[n]
  previous_sight <- times[n - 1L]
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  candidates <- seq(last_sight + 1, test_year)
  p_value <- (last_sight - previous_sight) /
    (candidates - previous_sight)

  if (data_out) return(data.frame(time = candidates, chance = p_value))

  below <- candidates[p_value <= alpha]
  if (length(below) == 0L) {
    warning("p-value never falls to alpha before `test_year`; returning NA.",
            call. = FALSE)
    return(new_ede_estimate(NA_real_, method = "Solow & Roberts (2003)",
                            alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Solow & Roberts (2003)", alpha = alpha)
}
