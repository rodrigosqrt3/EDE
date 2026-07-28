#' Optimal Linear Estimation of extinction date
#'
#' Estimates the extinction date from the sighting times with a positive
#' count, using the best linear unbiased estimator (BLUE) of Roberts & Solow
#' (2003) under a Weibull-type record-value model for the spacing of the
#' largest order statistics.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level for the confidence interval, in (0, 1).
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: When did the
#'   dodo become extinct? Nature, 426(6964), 245.
#'
#' @export
ole <- function(sd, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sort(rep(sd$time[sd$count > 0], sd$count[sd$count > 0]), decreasing = TRUE)
  k <- length(times)
  if (k < 3) {
    stop("OLE needs at least 3 sighting times with count > 0.", call. = FALSE)
  }

  t1 <- times[1]
  tk <- times[k]
  ti <- times[2:(k - 1)]

  if (t1 == times[2]) {
    stop("OLE requires the most recent sighting time (T1) to be strictly greater than T2.", call. = FALSE)
  }

  v <- sum(log((t1 - tk) / (t1 - ti))) / (k - 1)

  cov_ij <- function(i, j, v) {
    i_min <- pmin(i, j)
    j_max <- pmax(i, j)
    (gamma(2 * v + i_min) * gamma(v + j_max)) / (gamma(v + i_min) * gamma(j_max))
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
