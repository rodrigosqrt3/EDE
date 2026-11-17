#' Extinction time from dated remains, allowing for dating error (Solow et al. 2006)
#'
#' Maximum likelihood estimate of the extinction time of a taxon from the
#' radiometric ages of its most recent remains, with the dating error taken
#' into account (Solow, Roberts & Robbirt, 2006). The true ages of the remains
#' are assumed independent and uniformly distributed between a lower end-point
#' `gamma` and the extinction time `theta`; each recorded age is the true age
#' plus a normal dating error with a known standard deviation. The recorded
#' age then has density `(pnorm((x - gamma) / s) - pnorm((x - theta) / s)) /
#' (theta - gamma)` (equation 2 of the paper, written for a time axis that
#' increases towards the present). Because the dating error is normal and
#' positive the usual likelihood theory applies, which is not the case for the
#' same model without dating error ([strauss1989()], [robson1964()]).
#'
#' The uniform distribution is a strong assumption and holds, if at all, only
#' for the most recent part of the record; `k` selects the number of most
#' recent ages used, and [solow2006_check()] implements the simulation check of
#' the paper. For two taxa, [solow2006_region()] gives the joint confidence
#' region of their extinction times (equations 5 to 7 of the paper).
#'
#' The model for one taxon, the profile log-likelihood and the interval of
#' the extinction time (the set of values at which twice the drop of the
#' profile log-likelihood is at most the upper `alpha` quantile of the
#' chi-square distribution on one degree of freedom) are those of Solow (1997),
#' who applies them to the settlement time of New Zealand; Solow et al. (2006)
#' use the same likelihood for two taxa and a joint region, see
#' [solow2006_region()].
#'
#' The endpoint estimated is the one at the largest time. For an extinction,
#' ages in years before present enter as minus the age, so that the most
#' recent end of the range is the largest time; for a settlement or other first
#' appearance, as in Solow (1997), they enter as the age, so that the earliest
#' end is the largest time. The ages are on the scale of the source, without
#' calibration, which is outside this function.
#'
#' The dated remains are given in one of two ways. With `error`, the rows of a
#' [sighting_data] object with a positive count are the remains (a count above
#' one repeats the age with the same error). Without `error`, `sd` is a data
#' frame with one row per remain and columns `time` and `error`, which allows
#' ties in the times with different errors, as in radiocarbon data.
#'
#' @param sd A [sighting_data] object whose positive-count rows are the dated
#'   remains, or, without `error`, a data frame with columns `time` and `error`.
#' @param error Standard deviation of the dating error, positive, in the unit
#'   of the time axis: one value, or one per row of `sd` in the order of its
#'   rows (`sd` is sorted by time). Omit it when `sd` is a data frame with an
#'   `error` column.
#' @param alpha Significance level for the profile-likelihood interval, in
#'   (0, 1).
#' @param k Number of most recent ages to use; at least 3. By default all.
#'
#' @return An [ede_estimate] object (class `ede_solow2006`) with the maximum
#'   likelihood estimate of the extinction time in `estimate`, the interval in
#'   `lower` and `upper` (`Inf` if the upper limit does not exist), and the
#'   components `gamma` (estimated lower end-point), `loglik`, `k`, and `data`
#'   (the ages used, most recent first, with their errors).
#'
#' @references
#' Solow, A. R. (1997). Estimating settlement time. Radiocarbon, 39(3),
#'   351-354. doi:10.1017/S0033822200053327
#'
#' Solow, A. R., Roberts, D. L., & Robbirt, K. M. (2006). On the Pleistocene
#'   extinctions of Alaskan mammoths and horses. Proceedings of the National
#'   Academy of Sciences, 103(19), 7351-7353. doi:10.1073/pnas.0509480103
#'
#' @seealso [solow2006_region()], [solow2006_check()], [alaska_megafauna],
#'   [nz_shells].
#'
#' @examples
#' mammoth <- alaska_megafauna[alaska_megafauna$taxon == "mammoth", ]
#' # extinction: ka before present enters as minus the age
#' dates <- data.frame(time = -mammoth$age_ka_bp, error = mammoth$sd_ka)
#' fit <- solow2006(dates)
#' fit
#' fit$gamma
#'
#' # settlement of New Zealand (Solow 1997): the earliest end is the one wanted,
#' # so the age enters as it is
#' solow2006(data.frame(time = nz_shells$age_bp, error = nz_shells$sd_yr))
#'
#' @export
solow2006 <- function(sd, error = NULL, alpha = 0.05, k = NULL) {
  check_alpha(alpha)
  d <- sr_data(sd, error, k)
  fit <- sr_fit(d$x, d$s)
  lim <- sr_limits(fit, d$x, d$s, stats::qchisq(1 - alpha, 1))
  out <- new_ede_estimate(
    estimate = fit$theta, lower = lim[1L], upper = lim[2L],
    method = "Solow (1997); Solow, Roberts & Robbirt (2006)", alpha = alpha
  )
  if (is.infinite(lim[2L])) {
    warning("the upper confidence limit does not exist; returning Inf.",
            call. = FALSE)
  }
  out$gamma <- fit$gamma
  out$loglik <- fit$ll
  out$k <- length(d$x)
  out$data <- data.frame(time = d$x, error = d$s)
  class(out) <- c("ede_solow2006", class(out))
  out
}

#' Joint confidence region for the extinction times of two taxa
#'
#' Approximate joint confidence region for the extinction times of two taxa
#' from their dated remains, by the profile likelihood of the two-taxon model
#' of Solow, Roberts & Robbirt (2006, equations 5 to 7). Each taxon follows the
#' model of [solow2006()] and the two are independent, so the log-likelihood is
#' the sum and the profile log-likelihood of the pair of extinction times is
#' the sum of the single-taxon profiles. The region is the set of pairs for
#' which twice the drop of the profile log-likelihood from its maximum is at
#' most the upper `alpha` quantile of the chi-square distribution on two
#' degrees of freedom.
#'
#' The region answers questions about the relation of the two extinction times,
#' for instance whether it contains pairs in which the second taxon went
#' extinct after the first, or after a given event (the arrival of humans in
#' the paper), through the columns of `grid`.
#'
#' @param sd1,sd2 The dated remains of each taxon, as the `sd` argument of
#'   [solow2006()]: [sighting_data] objects, or data frames with columns `time`
#'   and `error`.
#' @param error1,error2 Dating error standard deviations for `sd1` and `sd2`,
#'   as the `error` argument of [solow2006()]; omit them for data frames.
#' @param alpha Significance level, in (0, 1).
#' @param k1,k2 Number of most recent ages of each taxon; by default all.
#' @param n Number of grid points on each axis.
#'
#' @return An object of class `solow2006_region`: a list with `fit1` and `fit2`
#'   (the single-taxon fits, with unadjusted 1-df intervals), `grid` (a data
#'   frame with `theta1`, `theta2`, `statistic` and `inside`), `critical` and
#'   `alpha`. Use `plot()` to draw the region.
#'
#' @references
#' Solow, A. R., Roberts, D. L., & Robbirt, K. M. (2006). On the Pleistocene
#'   extinctions of Alaskan mammoths and horses. Proceedings of the National
#'   Academy of Sciences, 103(19), 7351-7353. doi:10.1073/pnas.0509480103
#'
#' @seealso [solow2006()], [alaska_megafauna].
#'
#' @examples
#' mam <- alaska_megafauna[alaska_megafauna$taxon == "mammoth", ]
#' hor <- alaska_megafauna[alaska_megafauna$taxon == "horse", ]
#' hor <- hor[-which.max(hor$age_ka_bp), ] # the paper omits the oldest horse age
#' mk <- function(d) data.frame(time = -d$age_ka_bp, error = d$sd_ka)
#' reg <- solow2006_region(mk(mam), mk(hor), n = 40)
#' reg
#' # pairs in the region with the horse extinct later than 12.0 ka BP
#' range(reg$grid$theta2[reg$grid$inside])
#'
#' @export
solow2006_region <- function(sd1, sd2, error1 = NULL, error2 = NULL,
                             alpha = 0.05, k1 = NULL, k2 = NULL, n = 100L) {
  check_alpha(alpha)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 10) {
    stop("`n` must be a whole number of at least 10.", call. = FALSE)
  }
  n <- as.integer(n)
  d1 <- sr_data(sd1, error1, k1)
  d2 <- sr_data(sd2, error2, k2)
  f1 <- sr_fit(d1$x, d1$s)
  f2 <- sr_fit(d2$x, d2$s)
  crit2 <- stats::qchisq(1 - alpha, 2)
  axis <- function(f, d) {
    lim <- sr_limits(f, d$x, d$s, crit2)
    rng <- diff(range(d$x))
    lo <- lim[1L]
    hi <- if (is.finite(lim[2L])) lim[2L] else f$theta + 2 * rng
    seq(lo, hi, length.out = n)
  }
  g1 <- axis(f1, d1)
  g2 <- axis(f2, d2)
  s1 <- vapply(g1, sr_stat, 0, fit = f1, x = d1$x, s = d1$s)
  s2 <- vapply(g2, sr_stat, 0, fit = f2, x = d2$x, s = d2$s)
  grid <- expand.grid(theta1 = g1, theta2 = g2)
  grid$statistic <- as.vector(outer(s1, s2, "+"))
  grid$inside <- grid$statistic <= crit2
  mk <- function(f, d) {
    lim <- sr_limits(f, d$x, d$s, stats::qchisq(1 - alpha, 1))
    out <- new_ede_estimate(f$theta, lim[1L], lim[2L],
                            method = "Solow, Roberts & Robbirt (2006)",
                            alpha = alpha)
    out$gamma <- f$gamma
    out
  }
  structure(
    list(fit1 = mk(f1, d1), fit2 = mk(f2, d2), grid = grid, critical = crit2,
         alpha = alpha, n = n),
    class = "solow2006_region"
  )
}

#' @export
print.solow2006_region <- function(x, ...) {
  cat("<Joint confidence region, Solow, Roberts & Robbirt (2006)>\n")
  cat(sprintf("  maximum likelihood: theta1 = %s, theta2 = %s\n",
              format(round(x$fit1$estimate, 3)),
              format(round(x$fit2$estimate, 3))))
  ins <- x$grid[x$grid$inside, , drop = FALSE]
  cat(sprintf("  %g%% region (grid of %d x %d): theta1 in [%s, %s], theta2 in [%s, %s]\n",
              100 * (1 - x$alpha), x$n, x$n,
              format(round(min(ins$theta1), 3)),
              format(round(max(ins$theta1), 3)),
              format(round(min(ins$theta2), 3)),
              format(round(max(ins$theta2), 3))))
  invisible(x)
}

#' @export
plot.solow2006_region <- function(x, xlab = "theta1", ylab = "theta2", ...) {
  g1 <- sort(unique(x$grid$theta1))
  g2 <- sort(unique(x$grid$theta2))
  z <- matrix(x$grid$statistic, nrow = length(g1))
  graphics::contour(g1, g2, z, levels = x$critical, drawlabels = FALSE,
                    xlab = xlab, ylab = ylab, ...)
  graphics::abline(0, 1, lty = 2)
  graphics::points(x$fit1$estimate, x$fit2$estimate, pch = 3)
  invisible(x)
}

#' Simulation check of the fitted model of Solow et al. (2006)
#'
#' The check of the paper (its Figure 2). The fitted model is used to simulate
#' `nsim` records: `k` true ages are drawn uniformly between the estimated
#' lower end-point and extinction time, ordered, and the `j`-th is given the
#' normal dating error with the standard deviation of the `j`-th recorded age
#' (the ages ordered as recorded); the simulated ages are then ordered again.
#' If the model fits, the ordered recorded ages are close to the expected
#' ordered simulated ages and within their lower and upper `0.025` quantiles.
#'
#' @param fit A result of [solow2006()].
#' @param nsim Number of simulated records.
#' @param level Coverage of the band; the paper uses the `0.025` and `0.975`
#'   quantiles, that is `0.95`.
#'
#' @return A data frame of class `solow2006_check` with one row per age and
#'   columns `observed` (ordered recorded ages, oldest first), `expected` (mean
#'   of the ordered simulated ages), `lower` and `upper`. `plot()` draws the
#'   observed against the expected ages with the band and the 45 degree line.
#'
#' @references
#' Solow, A. R., Roberts, D. L., & Robbirt, K. M. (2006). On the Pleistocene
#'   extinctions of Alaskan mammoths and horses. Proceedings of the National
#'   Academy of Sciences, 103(19), 7351-7353. doi:10.1073/pnas.0509480103
#'
#' @seealso [solow2006()].
#'
#' @examples
#' mam <- alaska_megafauna[alaska_megafauna$taxon == "mammoth", ]
#' fit <- solow2006(data.frame(time = -mam$age_ka_bp, error = mam$sd_ka))
#' chk <- solow2006_check(fit, nsim = 200)
#' head(chk, 3)
#'
#' @export
solow2006_check <- function(fit, nsim = 1000, level = 0.95) {
  if (!inherits(fit, "ede_solow2006")) {
    stop("`fit` must be a result of solow2006().", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 10) {
    stop("`nsim` must be a whole number of at least 10.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  d <- fit$data
  o <- order(d$time) # oldest first
  obs <- d$time[o]
  s <- d$error[o]
  m <- length(obs)
  sims <- vapply(seq_len(nsim), function(i) {
    u <- sort(stats::runif(m, fit$gamma, fit$estimate))
    sort(u + stats::rnorm(m, 0, s))
  }, numeric(m))
  a <- (1 - level) / 2
  out <- data.frame(
    observed = obs,
    expected = rowMeans(sims),
    lower = apply(sims, 1L, stats::quantile, a, names = FALSE),
    upper = apply(sims, 1L, stats::quantile, 1 - a, names = FALSE)
  )
  class(out) <- c("solow2006_check", "data.frame")
  out
}

#' @export
plot.solow2006_check <- function(x, xlab = "expected ordered age",
                                 ylab = "observed ordered age", ...) {
  rng <- range(x$observed, x$expected, x$lower, x$upper)
  graphics::plot(x$expected, x$observed, type = "l", xlim = rng, ylim = rng,
                 xlab = xlab, ylab = ylab, ...)
  graphics::lines(x$expected, x$lower, lty = 3)
  graphics::lines(x$expected, x$upper, lty = 3)
  graphics::abline(0, 1, lty = 2)
  invisible(x)
}

# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

#' Dated remains for the Solow et al. (2006) model, most recent first
#' @keywords internal
#' @noRd
sr_data <- function(sd, error = NULL, k = NULL, kmin = 3L) {
  if (is.null(error)) {
    if (!is.data.frame(sd) || !all(c("time", "error") %in% names(sd))) {
      stop("without `error`, `sd` must be a data frame with columns `time` ",
           "and `error`, one row per dated remain.", call. = FALSE)
    }
    x <- sd$time
    s <- sd$error
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
      stop("`time` must be finite numbers.", call. = FALSE)
    }
    if (!is.numeric(s) || anyNA(s) || any(!is.finite(s)) || any(s <= 0)) {
      stop("`error` must be positive and finite.", call. = FALSE)
    }
  } else {
    stopifnot(inherits(sd, "sighting_data"))
    pos <- sd$count > 0
    if (!any(pos)) stop("no sightings with count > 0.", call. = FALSE)
    if (!is.numeric(error) || anyNA(error) || any(!is.finite(error)) ||
        any(error <= 0)) {
      stop("`error` must be positive and finite.", call. = FALSE)
    }
    if (length(error) == 1L) error <- rep(error, nrow(sd))
    if (length(error) != nrow(sd)) {
      stop("`error` must have one value, or one per row of `sd` (", nrow(sd),
           ").", call. = FALSE)
    }
    x <- rep(sd$time[pos], sd$count[pos])
    s <- rep(error[pos], sd$count[pos])
  }
  o <- order(x, decreasing = TRUE)
  x <- x[o]
  s <- s[o]
  n <- length(x)
  if (is.null(k)) k <- n
  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k != as.integer(k) ||
      k < kmin || k > n) {
    stop("`k` must be an integer between ", kmin, " and the number of ages (",
         n, ").", call. = FALSE)
  }
  list(x = x[seq_len(k)], s = s[seq_len(k)])
}

#' log(pnorm(a) - pnorm(b)) for a > b, stable in both tails
#' @keywords internal
#' @noRd
sr_logdiff <- function(a, b) {
  out <- numeric(length(a))
  up <- b > 0
  if (any(!up)) {
    la <- stats::pnorm(a[!up], log.p = TRUE)
    lb <- stats::pnorm(b[!up], log.p = TRUE)
    out[!up] <- la + log1p(-exp(lb - la))
  }
  if (any(up)) {
    ua <- stats::pnorm(-a[up], log.p = TRUE)
    ub <- stats::pnorm(-b[up], log.p = TRUE)
    out[up] <- ub + log1p(-exp(ua - ub))
  }
  out
}

#' Log-likelihood of the extinction time `theta` and lower end-point `gamma`
#' @keywords internal
#' @noRd
sr_ll <- function(theta, gamma, x, s) {
  sum(sr_logdiff((x - gamma) / s, (x - theta) / s)) -
    length(x) * log(theta - gamma)
}

#' Maximise over `gamma` at a fixed `theta`
#' @keywords internal
#' @noRd
sr_gamma_hat <- function(theta, x, s) {
  rng <- diff(range(x))
  lo <- min(x) - 5 * rng - 10 * max(s)
  hi <- min(min(x) + 5 * max(s), theta - 1e-6 * rng)
  o <- stats::optimize(function(g) sr_ll(theta, g, x, s), c(lo, hi),
                       maximum = TRUE, tol = 1e-10)
  list(gamma = o$maximum, ll = o$objective)
}

#' Joint maximum likelihood estimate
#' @keywords internal
#' @noRd
sr_fit <- function(x, s) {
  rng <- diff(range(x))
  if (!(rng > 0)) {
    stop("the ages must not all be equal.", call. = FALSE)
  }
  grid <- seq(max(x) - 5 * max(s), max(x) + 2 * rng + 5 * max(s),
              length.out = 80L)
  v <- vapply(grid, function(t) sr_gamma_hat(t, x, s)$ll, 0)
  i <- which.max(v)
  o <- stats::optimize(function(t) sr_gamma_hat(t, x, s)$ll,
                       c(grid[max(i - 1L, 1L)], grid[min(i + 1L, length(grid))]),
                       maximum = TRUE, tol = 1e-10)
  list(theta = o$maximum, gamma = sr_gamma_hat(o$maximum, x, s)$gamma,
       ll = o$objective)
}

#' Twice the drop of the profile log-likelihood at `theta`
#' @keywords internal
#' @noRd
sr_stat <- function(theta, fit, x, s) {
  max(0, 2 * (fit$ll - sr_gamma_hat(theta, x, s)$ll))
}

#' Profile-likelihood limits of the extinction time at the critical value `crit`
#' @keywords internal
#' @noRd
sr_limits <- function(fit, x, s, crit) {
  f <- function(t) sr_stat(t, fit, x, s) - crit
  step <- max(max(s), diff(range(x)) / 20)
  find <- function(dir) {
    b <- fit$theta + dir * step
    for (i in seq_len(60L)) {
      if (f(b) > 0) {
        return(stats::uniroot(f, sort(c(fit$theta, b)), tol = 1e-9)$root)
      }
      b <- fit$theta + dir * step * 2^i
    }
    if (dir > 0) Inf else -Inf
  }
  c(find(-1), find(1))
}
