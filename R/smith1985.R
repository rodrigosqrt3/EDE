#' Maximum likelihood end-point and profile-likelihood interval (Smith, 1985)
#'
#' Fits the generalised Pareto model to the excesses of the `k` latest
#' sightings over the `(k + 1)`-th latest, in the Weibull domain, and returns
#' the maximum likelihood estimate of the end-point (the extinction time)
#' with a profile-likelihood interval. This is the likelihood treatment of
#' non-regular end-point problems of Smith (1985, Section 7), applied to
#' sighting times.
#'
#' With `u = t[k + 1]` the threshold (sightings ordered latest first) and
#' `theta` the end-point, the excesses `y = t - u` satisfy
#' `y / (theta - u) ~ Beta(1, a)`, where the shape `a = -1 / gamma` is the
#' extremal exponent of [ole()] and [hall_wang1999()] (`gamma` is the
#' extreme-value index, negative in the Weibull domain). Maximising over `a`
#' for fixed `theta` gives the closed-form profile log-likelihood
#' `k * log(a) - k * log(theta - u) - k - sum(log((theta - t) / (theta - u)))`
#' with `a = -k / sum(log((theta - t) / (theta - u)))`, which is maximised over
#' `theta`. The interval is `{theta : 2 * (l(theta_hat) - l(theta)) <=
#' qchisq(1 - alpha, 1)}` on the profile `l`.
#'
#' **Regularity.** The likelihood is unbounded as `theta` approaches the last
#' sighting, so the maximum likelihood estimate is a local interior maximum,
#' and it exists, and is asymptotically normal, only when the shape is larger
#' than 2 (`gamma > -1/2`; Smith, 1985, Theorems 1 and 3). For `a < 2` the
#' problem is non-regular, and for `a < 1` (`gamma < -1`) no consistent
#' maximum likelihood estimator exists. A sighting rate that is constant until
#' extinction gives `a = 1`, so in the setting that motivates most extinction
#' models the likelihood treatment is outside its regular range. When the
#' fitted shape is not larger than 2 the function warns that the chi-squared
#' calibration of the interval is not justified; [hall_wang1999()], whose
#' interval is exact for any shape, is then the better choice. The estimate is
#' `NA`, with a warning, when the profile likelihood has no interior maximum.
#'
#' In simulations (300 records each, nominal 95% interval) the estimate
#' existed in 61% to 97% of records with shape 3 to 4 and `k` of 15 to 40,
#' in 53% with shape 2 and in 6% with shape 1 (`k` 10 to 15). Where it existed
#' with shape above 2 the interval covered the true end-point in 91% to 96% of
#' records, against 93% to 97% for [hall_wang1999()] in all cases.
#'
#' An upper limit that does not exist is returned as `Inf`, as in [ole()].
#' By default each distinct sighting time counts once; set `events = TRUE` to
#' count every sighting event.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1); the interval has nominal
#'   coverage `1 - alpha`.
#' @param k Number of latest sightings, an integer of at least 4 with `k + 1`
#'   at most the number of sighting times, or `NULL` for the 10 latest (fewer
#'   if the record is shorter).
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#'
#' @return An [ede_estimate] object with `estimate` (possibly `NA`), `lower`,
#'   `upper`, and components `shape` (the extremal exponent `a`), `gamma`
#'   (`-1 / a`), `k`, `loglik` (maximised profile log-likelihood) and
#'   `regular` (`TRUE` when the shape exceeds 2).
#'
#' @references
#' Smith, R. L. (1985). Maximum likelihood estimation in a class of
#'   nonregular cases. Biometrika, 72(1), 67-90. doi:10.1093/biomet/72.1.67
#'
#' @seealso [hall_wang1999()] and [ole()] for the same model with exact
#'   intervals; [ede_endpoint_test()] to test the Gumbel null.
#'
#' @examples
#' # the dodo record: shape near 2, so the interval is not justified (warning)
#' smith1985(sighting_data(dodo), k = 8)
#'
#' @export
smith1985 <- function(sd, alpha = 0.05, k = NULL, events = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  t <- evt_times(sd, events)
  n <- length(t)
  k <- evt_k(k, n, kmin = 4L, default = min(n - 1L, 10L),
             label = "Smith (1985)", kmax = n - 1L)
  u <- t[k + 1L]
  x <- t[seq_len(k)]
  if (!(x[1L] > u)) {
    stop("the `k` latest sightings must span more than one time.",
         call. = FALSE)
  }
  g <- x[1L] - u   # span; the end-point is theta = t[1] + s, s > 0
  fit <- smith_fit(x, u, g)
  crit <- stats::qchisq(1 - alpha, 1) / 2

  if (is.na(fit$s)) {
    warning("the profile likelihood has no interior maximum for this record ",
            "(the data favour the exponential model or a spike at the last ",
            "sighting); returning NA.", call. = FALSE)
    ci <- c(NA_real_, NA_real_)
  } else {
    ci <- smith_interval(x, u, g, fit, crit)
    if (is.infinite(ci[2L])) {
      warning("the upper confidence limit does not exist at this alpha and ",
              "record; returning Inf.", call. = FALSE)
    }
    if (fit$shape <= 2) {
      warning("the fitted shape (", format(fit$shape, digits = 3),
              ") is not larger than 2: the problem is non-regular (Smith, ",
              "1985), and the chi-squared calibration of the interval is not ",
              "justified. Consider hall_wang1999().", call. = FALSE)
    }
  }
  out <- new_ede_estimate(
    estimate = if (is.na(fit$s)) NA_real_ else x[1L] + fit$s,
    lower = if (anyNA(ci)) NA_real_ else x[1L] + ci[1L],
    upper = if (anyNA(ci)) NA_real_ else x[1L] + ci[2L],
    method = "Smith (1985)", alpha = alpha, interval_type = "two-sided"
  )
  out$shape <- fit$shape
  out$gamma <- if (is.na(fit$shape)) NA_real_ else -1 / fit$shape
  out$k <- k
  out$loglik <- fit$ll
  out$regular <- isTRUE(fit$shape > 2)
  out
}

#' Profile log-likelihood of the end-point `theta = x[1] + s`
#' @keywords internal
#' @noRd
smith_profile <- function(s, x, u, g) {
  k <- length(x)
  r <- (x[1L] - x + s) / (x[1L] - u + s)   # (theta - t_i) / (theta - u)
  slr <- sum(log(r))
  a <- -k / slr
  k * log(a) - k * log(x[1L] - u + s) - k - slr
}

#' Limit of the profile log-likelihood as the end-point goes to infinity
#' (exponential model for the excesses)
#' @keywords internal
#' @noRd
smith_profile_inf <- function(x, u) {
  y <- x - u
  -length(x) * log(mean(y)) - length(x)
}

#' Interior local maximum of the profile likelihood
#' @keywords internal
#' @noRd
smith_fit <- function(x, u, g) {
  k <- length(x)
  grid <- g * exp(seq(log(1e-4), log(1e5), length.out = 600L))
  ll <- vapply(grid, smith_profile, numeric(1), x = x, u = u, g = g)
  ll_inf <- smith_profile_inf(x, u)
  i <- which(ll[-c(1L, length(ll))] > ll[-c(length(ll) - 1L, length(ll))] &
             ll[-c(1L, length(ll))] > ll[-c(1L, 2L)]) + 1L
  if (!length(i)) return(list(s = NA_real_, shape = NA_real_, ll = NA_real_))
  i <- i[which.max(ll[i])]
  o <- stats::optimize(function(ls) smith_profile(exp(ls), x, u, g),
                       log(grid[c(i - 1L, i + 1L)]), maximum = TRUE,
                       tol = 1e-12)
  s <- exp(o$maximum)
  if (!(o$objective > ll_inf)) {
    return(list(s = NA_real_, shape = NA_real_, ll = NA_real_))
  }
  r <- (x[1L] - x + s) / (x[1L] - u + s)
  list(s = s, shape = -k / sum(log(r)), ll = o$objective)
}

#' Profile-likelihood interval around the interior maximum
#' @keywords internal
#' @noRd
smith_interval <- function(x, u, g, fit, crit) {
  target <- fit$ll - crit
  f <- function(ls) smith_profile(exp(ls), x, u, g) - target
  # upper limit: first crossing beyond the maximum; Inf if none up to 1e5 spans
  up <- fit$s * exp(seq(0, log(1e5 * g / fit$s), length.out = 400L))
  fu <- vapply(log(up), f, numeric(1))
  below <- which(fu < 0)
  upper <- if (length(below)) {
    j <- below[1L]
    exp(stats::uniroot(f, log(up[c(j - 1L, j)]), tol = 1e-12)$root)
  } else if (smith_profile_inf(x, u) >= target) Inf else {
    up[length(up)]
  }
  # lower limit: the nearest crossing below the maximum, ignoring the spike at
  # the last sighting; the end-point cannot be earlier than the last sighting
  lo <- fit$s * exp(seq(0, log(1e-6), length.out = 400L))
  fl <- vapply(log(lo), f, numeric(1))
  belowl <- which(fl < 0)
  lower <- if (length(belowl)) {
    j <- belowl[1L]
    exp(stats::uniroot(f, log(lo[c(j, j - 1L)]), tol = 1e-12)$root)
  } else 0
  c(lower, upper)
}
