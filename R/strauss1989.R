#' Strauss & Sadler (1989) confidence interval for the end of a range
#'
#' Unbiased point estimator and classical one-sided confidence interval for
#' the true endpoint of a temporal range, assuming finds are uniformly
#' distributed between the true endpoints. With `n` finds spanning the range
#' `y` to `z`, the estimate is `z + (z - y) / (n - 1)` (equation 8 of the
#' source paper) and the interval is `(z, z + c (z - y))` with
#' `c = alpha^(-1 / (n - 1)) - 1` (equations 13, 16 and 20).
#'
#' Each distinct sighting time is one find, as each fossil horizon is in the
#' paper; counts are ignored. The paper's model is continuous, so with dense
#' yearly records the unbiasedness and the coverage hold only approximately.
#' Both endpoints are estimated from the data, so the observation origin (a
#' leading zero-count row) is not used.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1).
#'
#' @return An [ede_estimate] object containing the unbiased point estimate and
#'   a one-sided confidence interval.
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

  times <- sort(unique(sd$time[sd$count > 0]))
  h <- length(times)
  if (h < 2) {
    stop("Strauss & Sadler (1989) needs at least 2 distinct sighting times with count > 0.", call. = FALSE)
  }

  range_r <- max(times) - min(times)
  b <- alpha^(-1 / (h - 1)) - 1
  estimate <- max(times) + range_r / (h - 1)
  upper <- max(times) + b * range_r

  new_ede_estimate(estimate, lower = max(times), upper = upper,
                   method = "Strauss & Sadler (1989)", alpha = alpha,
                   interval_type = "one-sided")
}

#' Confidence curve for the Strauss & Sadler (1989) upper bound
#'
#' Evaluates the upper confidence bound of [strauss1989()],
#' `z + c (z - y)` with `c = alpha^(-1 / (n - 1)) - 1`, for alpha from 0.01
#' to 1 in steps of 0.01, for plotting instead of a single bound.
#'
#' @param sd A [sighting_data] object.
#' @return A data frame with columns `time` and `chance`, where `chance` is
#'   the confidence level `1 - alpha` of the bound `time`, i.e. the
#'   probability that the true endpoint is at or before `time`.
#' @export
strauss1989_curve <- function(sd) {
  stopifnot(inherits(sd, "sighting_data"))
  times <- sort(unique(sd$time[sd$count > 0]))
  h <- length(times)
  if (h < 2) {
    stop("Strauss & Sadler (1989) needs at least 2 distinct sighting times with count > 0.", call. = FALSE)
  }

  alphas <- seq(0.01, 1, by = 0.01)
  range_r <- max(times) - min(times)
  b <- alphas^(-1 / (h - 1)) - 1
  time_est <- max(times) + b * range_r

  data.frame(time = rev(time_est), chance = rev(1 - alphas))
}

#' Strauss & Sadler (1989) Bayesian estimate of the end of a range
#'
#' Bayesian counterpart of [strauss1989()] (Section 4.3 of the source paper).
#' Finds are independent and uniform between the true endpoints
#' `theta1 < theta2` (equation 23), and the prior on the endpoints is uniform
#' over `a < theta1 < theta2 < b` (equation 22, there scaled to `(0, 1)`),
#' where `(a, b) = prior_range` is the interval that the endpoints are known
#' to lie in. With `n` finds spanning `y` to `z`, the posterior density of
#' the extinction time `theta2` is
#' `h(theta2) = (n - 2) [(theta2 - y)^(1-n) - (theta2 - a)^(1-n)] / u_n` on
#' `z < theta2 < b` (equation 26), with
#' `u_n = (z-y)^(2-n) - (b-y)^(2-n) - (z-a)^(2-n) + (b-a)^(2-n)` (equation
#' 25). The estimate is the posterior mean (equation 27), the Bayes point
#' estimate under quadratic loss. Since the density decreases from its
#' maximum at `z`, the shortest `1 - alpha` posterior probability interval is
#' `(z, q)`, with `q` the `1 - alpha` posterior quantile (as shaded in
#' Fig. 3 of the paper); it is returned in `lower` and `upper`.
#'
#' The paper writes the equations for `(a, b) = (0, 1)`; they are applied
#' here after the corresponding linear change of scale, which leaves the
#' posterior unchanged. The paper suggests `(0, Inf)` for taxa that may be
#' extant; `b = Inf` is allowed (the posterior is proper for `n >= 3`, but
#' its mean is infinite for `n = 3`). The paper's example takes `(0, 1500)`
#' m for the ammonite section of Seymour Island.
#'
#' As in [strauss1989()], each distinct sighting time is one find and counts
#' are ignored; the observation origin (a leading zero-count row) is not
#' used, but it may be a natural choice for `a`. The paper requires
#' `n > 2`.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Posterior probability outside the interval, in (0, 1).
#' @param prior_range Numeric vector `c(a, b)`: the prior is uniform for the
#'   endpoints over `a < theta1 < theta2 < b`. `a` must be finite and earlier
#'   than the first find; `b` must be later than the last find and may be
#'   `Inf`.
#' @param data_out If `TRUE`, return the posterior on a grid of 200 times
#'   from `z` to the 0.999 posterior quantile instead of the estimate.
#'
#' @return An [ede_estimate] object with the posterior mean in `estimate`
#'   and the one-sided posterior probability interval `(z, q)` in `lower` and
#'   `upper`, plus components `density` and `cdf` (functions of time giving
#'   the posterior density and distribution function of the extinction time)
#'   and `prior_range`. If `data_out = TRUE`, a data frame with columns
#'   `time`, `density` and `chance`, where `chance` is the posterior
#'   probability that the extinction time is at or before `time`.
#'
#' @references
#' Strauss, D., & Sadler, P. M. (1989). Classical confidence intervals and
#'   Bayesian probability estimates for ends of local taxon ranges.
#'   Mathematical Geology, 21(4), 411-427.
#'
#' @export
strauss1989_bayes <- function(sd, alpha = 0.05, prior_range, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 3) {
    stop("Strauss & Sadler (1989) Bayesian estimation needs at least 3 distinct sighting times with count > 0.",
         call. = FALSE)
  }
  y <- min(times)
  z <- max(times)
  if (missing(prior_range) || !is.numeric(prior_range) ||
      length(prior_range) != 2L || anyNA(prior_range)) {
    stop("`prior_range` must be supplied as c(lower, upper).", call. = FALSE)
  }
  a <- prior_range[1]
  b <- prior_range[2]
  if (!is.finite(a) || a >= y) {
    stop("`prior_range[1]` must be finite and earlier than the first sighting.", call. = FALSE)
  }
  if (b <= z) {
    stop("`prior_range[2]` must be later than the last sighting.", call. = FALSE)
  }

  post <- strauss1989_posterior(n, y, z, a, b)
  q <- post$quantile(1 - alpha)

  if (data_out) {
    grid <- seq(z, post$quantile(0.999), length.out = 200)
    return(data.frame(time = grid, density = post$density(grid),
                      chance = post$cdf(grid)))
  }

  res <- new_ede_estimate(post$mean, lower = z, upper = q,
                          method = "Strauss & Sadler (1989), Bayesian",
                          alpha = alpha, interval_type = "one-sided")
  res$interval_kind <- "credible"
  res$density <- post$density
  res$cdf <- post$cdf
  res$prior_range <- c(a, b)
  res
}

#' Posterior of the upper endpoint, Strauss & Sadler (1989) eqs. 25-27.
#' Powers are taken relative to r = z - y to avoid overflow for large n.
#' @keywords internal
#' @noRd
strauss1989_posterior <- function(n, y, z, a, b) {
  r <- z - y
  # u_k scaled by r^k: (z-y)^-k - (b-y)^-k - (z-a)^-k + (b-a)^-k
  u_scaled <- function(k) 1 - (r / (b - y))^k - (r / (z - a))^k + (r / (b - a))^k
  u <- u_scaled(n - 2)

  density <- function(t) {
    inside <- t >= z & t <= b
    out <- numeric(length(t))
    ti <- t[inside]
    out[inside] <- (n - 2) / r * ((r / (ti - y))^(n - 1) - (r / (ti - a))^(n - 1)) / u
    out
  }
  cdf <- function(t) {
    tt <- pmin(pmax(t, z), b)
    p <- (1 - (r / (tt - y))^(n - 2) - (r / (z - a))^(n - 2) + (r / (tt - a))^(n - 2)) / u
    pmin(pmax(p, 0), 1)
  }
  quantile <- function(p) {
    upper0 <- if (is.finite(b)) b else z + r
    stats::uniroot(function(t) cdf(t) - p, lower = z, upper = upper0,
                   extendInt = if (is.finite(b)) "no" else "upX",
                   tol = 1e-10 * max(1, abs(z)))$root
  }

  mean <- if (n > 3) {
    # equation 27, in the original scale
    (n - 2) / u * (r * u_scaled(n - 3) / (n - 3) +
                     (y * (1 - (r / (b - y))^(n - 2)) -
                        a * ((r / (z - a))^(n - 2) - (r / (b - a))^(n - 2))) / (n - 2))
  } else if (is.finite(b)) {
    z + stats::integrate(function(t) 1 - cdf(t), z, b)$value
  } else {
    Inf
  }

  list(density = density, cdf = cdf, quantile = quantile, mean = mean)
}
