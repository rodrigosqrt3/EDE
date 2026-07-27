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

  times <- sort(sd$time[sd$count > 0], decreasing = TRUE)
  k <- length(times)
  if (k < 3) {
    stop("OLE needs at least 3 sighting times with count > 0.", call. = FALSE)
  }

  # shape parameter of the assumed record-value model, from the spacing
  # of the k largest sighting times (Roberts & Solow 2003, eq. 1).
  # NOTE: normalised by (k - 1), not by the number of summed terms (k - 2).
  v <- sum(log((times[1] - times[k]) / (times[1] - times[2:(k - 1)]))) / (k - 1)

  # BLUE weights under that model
  cov_ij <- function(i, j, v) (gamma(2 * v + i) * gamma(v + j)) / (gamma(v + i) * gamma(j))
  lambda <- outer(seq_len(k), seq_len(k), cov_ij, v = v)
  lambda[upper.tri(lambda)] <- t(lambda)[upper.tri(lambda)]
  lambda_inv <- solve(lambda)
  e <- rep(1, k)
  weights <- as.vector(1 / (t(e) %*% lambda_inv %*% e)) * (lambda_inv %*% e)
  estimate <- sum(weights * times)

  span <- times[1] - times[k]
  sl <- (-log(1 - alpha / 2) / k)^(-v)
  su <- (-log(alpha / 2) / k)^(-v)
  lower <- times[1] + span / (sl - 1)
  upper <- times[1] + span / (su - 1)

  new_ede_estimate(estimate, lower, upper, method = "OLE (Roberts & Solow 2003)", alpha = alpha)
}
