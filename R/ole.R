#' Optimal Linear Estimation of extinction date
#'
#' Estimates the extinction date from the sighting times with a positive
#' count, using the best linear unbiased estimator (BLUE) of Roberts & Solow
#' (2003) under a Weibull-type record-value model for the spacing of the
#' largest order statistics.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level for the confidence interval, in (0, 1).
#' @param k Number of most recent sighting events to use. Counts greater than
#'   one are treated as independent events at the recorded time. By default,
#'   all sighting events are used. Must be at least 3.
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: When did the
#'   dodo become extinct? Nature, 426(6964), 245.
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
  upper <- t1 + span / (su - 1)

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
