#' Robson & Whitlock (1964) truncation point estimator
#'
#' Estimates the extinction date by jackknife bias correction of the most
#' recent sighting. If the two most recent distinct sighting times are
#' `t[n-1]` and `t[n]`, the estimate is `t[n] + (t[n] - t[n-1])`.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1), for the approximate one-sided
#'   confidence interval.
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Robson, D. S., & Whitlock, J. H. (1964). Estimation of a truncation
#'   point. Biometrika, 51(1/2), 33-39.
#'
#' @export
robson1964 <- function(sd, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 2) {
    stop("Robson & Whitlock (1964) needs at least 2 distinct sighting times with count > 0.", call. = FALSE)
  }

  gap <- times[n] - times[n - 1]
  estimate <- times[n] + gap
  upper <- times[n] + gap * (1 - alpha) / alpha

  new_ede_estimate(estimate, lower = times[n], upper = upper,
                   method = "Robson & Whitlock (1964)", alpha = alpha,
                   interval_type = "one-sided")
}
