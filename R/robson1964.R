#' Robson & Whitlock (1964) truncation point estimator
#'
#' Estimates the extinction date by a Quenouille (jackknife-type) bias
#' reduction of the most recent sighting. If the two most recent distinct
#' sighting times are `t[n-1]` and `t[n]`, the estimate is
#' `t[n] + (t[n] - t[n-1])` and the one-sided upper confidence bound is
#' `t[n] + (t[n] - t[n-1]) (1 - alpha) / alpha`.
#'
#' The estimator is mean-unbiased to order `n^-2` in general and exactly
#' unbiased for a uniform distribution. The bound is exact for a uniform
#' distribution and asymptotically valid when the density is continuous and
#' positive at the endpoint; its coverage is below `1 - alpha` when the
#' density decreases towards the endpoint, as with a declining sighting rate.
#'
#' Only distinct sighting times are used; counts are ignored, so the paper's
#' order statistics are applied to distinct times.
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
