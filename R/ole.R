#' Optimal Linear Estimation of extinction date
#'
#' Estimates the extinction date by optimal linear estimation (Cooke, 1980),
#' as applied by Roberts & Solow (2003). The joint distribution of the `k`
#' most recent sighting times (the extreme order statistics) is assumed to
#' have the Weibull extreme-value form; the estimate is the weighted sum of
#' those times given by equations 12-15 of Solow (2005).
#'
#' The confidence interval is the two-sided `1 - alpha` interval of Roberts &
#' Solow (2003). Its upper limit does not exist when `-log(alpha / 2) >= k`;
#' it is then returned as `Inf` with a warning. [solow2005()] inverts the
#' corresponding test one-sidedly (Solow, 2005, equation 17), so its bound at
#' level `alpha` equals the upper limit here at level `2 * alpha`.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level for the confidence interval, in (0, 1).
#' @param k Number of most recent sighting events to use. Counts greater than
#'   one are treated as independent events at the recorded time; Solow (2005)
#'   recommends combining sightings that are obviously related. By default,
#'   all sighting events are used. Must be at least 3; Solow (2005) suggests
#'   at least 5, and notes that the extreme-value approximation may fail if
#'   `k` is too large.
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Cooke, P. (1980). Optimal linear estimation of bounds of random variables.
#'   Biometrika, 67(1), 257-258.
#'
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: When did the
#'   dodo become extinct? Nature, 426(6964), 245.
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @export
ole <- function(sd, alpha = 0.05, k = NULL) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sort(rep(sd$time[sd$count > 0], sd$count[sd$count > 0]), decreasing = TRUE)
  n_times <- length(times)
  if (n_times < 3) {
    stop("OLE needs at least 3 sighting times with count > 0.", call. = FALSE)
  }
  if (is.null(k)) k <- n_times
  if (length(k) != 1L || !is.numeric(k) || is.na(k) ||
      k != as.integer(k) || k < 3L || k > n_times) {
    stop("`k` must be an integer between 3 and the number of sighting events.", call. = FALSE)
  }
  k <- as.integer(k)
  times <- times[seq_len(k)]

  t1 <- times[1]
  tk <- times[k]

  if (t1 == times[2]) {
    stop("OLE requires the most recent sighting time (T1) to be strictly greater than T2.", call. = FALSE)
  }

  v <- ole_shape(times)

  cov_ij <- function(i, j, v) {
    i_max <- pmax(i, j)
    j_min <- pmin(i, j)
    exp(lgamma(2 * v + i_max) + lgamma(v + j_min) -
          lgamma(v + i_max) - lgamma(j_min))
  }

  lambda <- outer(seq_len(k), seq_len(k), cov_ij, v = v)
  e <- rep(1, k)
  x <- solve(lambda, e)
  weights <- as.vector(x / sum(x))

  estimate <- sum(weights * times)

  span <- t1 - tk
  sl <- (-log(1 - alpha / 2) / k)^(-v)
  su <- (-log(alpha / 2) / k)^(-v)
  lower <- t1 + span / (sl - 1)
  upper <- if (su > 1) t1 + span / (su - 1) else Inf
  if (is.infinite(upper)) {
    warning(
      "the upper confidence limit does not exist for k = ", k,
      " at this alpha (it requires -log(alpha / 2) < k); returning Inf.",
      call. = FALSE
    )
  }

  new_ede_estimate(estimate, lower, upper, method = "OLE (Roberts & Solow 2003)", alpha = alpha)
}

#' Estimate the Weibull extreme-value shape parameter
#' @keywords internal
#' @noRd
ole_shape <- function(times) {
  times <- sort(times, decreasing = TRUE)
  k <- length(times)
  if (k < 3L || times[1] <= times[k]) return(NA_real_)
  sum(log((times[1] - times[k]) /
            (times[1] - times[2:(k - 1)]))) / (k - 1)
}
