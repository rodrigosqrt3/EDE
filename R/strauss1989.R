#' Strauss & Sadler (1989) confidence interval for the end of a range
#'
#' Classical confidence interval for the true endpoint of a temporal range,
#' derived from the distribution of the sample range under a uniform
#' occurrence model. Originally developed for stratigraphic ranges, applied
#' here to sighting records.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1).
#'
#' @return An [ede_estimate] object. Only `upper` is defined: the method
#'   gives a one-sided bound on how much later than the last sighting
#'   extinction could plausibly have occurred, not a point estimate with a
#'   two-sided interval.
#'
#' @references
#' Strauss, D., & Sadler, P. M. (1989). Classical confidence intervals and
#'   Bayesian probability estimates for ends of local taxon ranges.
#'   Mathematical Geology, 21(4), 411-427.
#'
#' @export
strauss1989 <- function(sd, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sd$time[sd$count > 0]
  h <- length(times)
  if (h < 2) {
    stop("Strauss & Sadler (1989) needs at least 2 sighting times with count > 0.", call. = FALSE)
  }

  range_r <- max(times) - min(times)
  b <- alpha^(-1 / (h - 1)) - 1
  estimate <- max(times) + b * range_r

  new_ede_estimate(estimate, method = "Strauss & Sadler (1989)", alpha = alpha)
}

#' Full confidence curve for the Strauss & Sadler (1989) estimator
#'
#' Same estimator as [strauss1989()], evaluated over a grid of alpha values
#' from 0.01 to 1, for plotting the confidence curve instead of a single
#' bound.
#'
#' @param sd A [sighting_data] object.
#' @return A data frame with columns `time` and `chance` (1 - alpha).
#' @export
strauss1989_curve <- function(sd) {
  stopifnot(inherits(sd, "sighting_data"))
  times <- sd$time[sd$count > 0]
  h <- length(times)
  if (h < 2) {
    stop("Strauss & Sadler (1989) needs at least 2 sighting times with count > 0.", call. = FALSE)
  }

  alphas <- seq(0.01, 1, by = 0.01)
  range_r <- max(times) - min(times)
  b <- alphas^(-1 / (h - 1)) - 1
  time_est <- max(times) + b * range_r

  data.frame(time = rev(time_est), chance = rev(1 - alphas))
}
