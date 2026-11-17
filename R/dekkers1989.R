#' Moment estimator of the extreme-value index and of the end-point
#'
#' Estimates the extreme-value index `gamma` of the sighting-time
#' distribution, and from it the extinction time, with the moment estimator of
#' Dekkers, Einmahl & de Haan (1989) applied to the `k` latest sightings. The
#' index separates the three domains of attraction: `gamma < 0` is the
#' Weibull domain, a finite end-point, the setting of [ole()] and
#' [hall_wang1999()] (the Weibull extremal exponent of those functions is
#' `-1 / gamma`); `gamma = 0` is the Gumbel domain, and `gamma > 0` an
#' unbounded upper tail. Unlike those estimators it does not assume the
#' Weibull domain, so it also says whether the assumption is plausible.
#'
#' With `x[i]` the sighting times measured from `origin`, latest first, and
#' the threshold `x[k + 1]`, let
#' `M1 = mean(log(x[i] / x[k + 1]))` and `M2 = mean(log(x[i] / x[k + 1])^2)`
#' for `i = 1, ..., k`. The index estimate is
#' `gamma = M1 + 1 - 0.5 / (1 - M1^2 / M2)`. For `gamma < 0` the end-point
#' estimate is `x[k + 1] - a / gamma`, with scale
#' `a = x[k + 1] * M1 * (1 - gamma_minus)` and
#' `gamma_minus = 1 - 0.5 / (1 - M1^2 / M2)` (de Haan & Ferreira, 2006,
#' Section 4.5). The end-point estimate is `NA`, with a warning, when the
#' estimated index is not negative.
#'
#' The estimator uses logarithms of the observations, so it requires the times
#' to be positive and, unlike the other estimators in this package, it is not
#' invariant to the origin of the time axis (de Haan & Ferreira, 2006, p.
#' 102). The default `origin` is the earliest time in the record, the start of
#' observation; the estimate depends on that choice, so check it with
#' different origins. The estimator needs `k + 1` sighting times and the
#' threshold must be later than `origin`.
#'
#' The asymptotic standard error of the index estimate is
#' `sqrt(v / k)` with `v = 1 + gamma^2` for `gamma >= 0` and
#' `v = (1 - gamma)^2 (1 - 2 gamma) (1 - gamma + 6 gamma^2) /
#' ((1 - 3 gamma) (1 - 4 gamma))` for `gamma < 0`. It is a large-sample
#' result, and with the few sightings typical of extinction records the
#' interval is only a guide. No interval is given for the end-point.
#'
#' By default each distinct sighting time counts once, as in
#' [hall_wang1999()]; set `events = TRUE` to count every sighting event.
#'
#' @param sd A [sighting_data] object.
#' @param k Number of latest sightings, an integer of at least 3 with `k + 1`
#'   at most the number of sighting times later than `origin`, or `NULL` for
#'   the 10 latest (fewer if the record is shorter). With the default
#'   `origin`, the earliest sighting cannot serve as threshold, so `k` is at
#'   most the number of sighting times minus two.
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#' @param origin Origin of the time axis, earlier than the threshold
#'   sighting `k + 1`; `NULL` (default) uses the earliest time in `sd`.
#' @param alpha Significance level, in (0, 1), of the normal interval for
#'   `gamma`.
#'
#' @return An [ede_estimate] object with `estimate` (possibly `NA`, see
#'   Details), no interval, and components `gamma`, `gamma_se`,
#'   `gamma_ci` (a `1 - alpha` normal interval), `k` and `origin`.
#'
#' @references
#' Dekkers, A. L. M., Einmahl, J. H. J., & de Haan, L. (1989). A moment
#'   estimator for the index of an extreme-value distribution. Annals of
#'   Statistics, 17(4), 1833-1855. doi:10.1214/aos/1176347397
#'
#' de Haan, L., & Ferreira, A. (2006). Extreme Value Theory: An Introduction.
#'   Springer, Sections 3.5 and 4.5.
#'
#' @seealso [ede_evi()] for the estimates over a range of `k`;
#'   [ede_endpoint_test()] to test the Gumbel null; [hall_wang1999()] and
#'   [ole()] for the Weibull domain.
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' dekkers1989(sd, k = 6)
#'
#' @export
dekkers1989 <- function(sd, k = NULL, events = FALSE, origin = NULL,
                        alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  t <- evt_times(sd, events)
  n <- length(t)
  origin <- moment_origin(sd, origin)
  n <- sum(t > origin)   # the threshold must be later than the origin
  k <- evt_k(k, n, kmin = 3L, default = min(n - 1L, 10L),
             label = "Dekkers et al. (1989)", kmax = n - 1L)
  fit <- moment_fit(t - origin, k)
  if (is.na(fit$gamma)) {
    stop("the `k` latest sightings must span more than one time.",
         call. = FALSE)
  }
  se <- sqrt(moment_var(fit$gamma) / k)
  z <- stats::qnorm(1 - alpha / 2)
  est <- if (is.na(fit$endpoint)) NA_real_ else origin + fit$endpoint
  if (is.na(est)) {
    warning("the estimated extreme-value index is not negative (",
            format(fit$gamma, digits = 3), "), so the data suggest no finite ",
            "end-point; returning NA.", call. = FALSE)
  }
  out <- new_ede_estimate(est, method = "Dekkers et al. (1989)",
                          alpha = alpha, interval_type = "two-sided")
  out$gamma <- fit$gamma
  out$gamma_se <- se
  out$gamma_ci <- c(fit$gamma - z * se, fit$gamma + z * se)
  out$k <- k
  out$origin <- origin
  out
}

#' Moment estimates over a range of k (stability diagram)
#'
#' Applies the moment estimator of [dekkers1989()] for every `k` in a range,
#' the usual diagnostic for choosing `k` in extreme-value analysis: look for a
#' range of `k` over which the estimate of the extreme-value index is stable.
#'
#' @inheritParams dekkers1989
#' @param kmin,kmax Smallest and largest `k`; `kmax = NULL` uses the largest
#'   admissible value (the threshold needs one more sighting later than
#'   `origin`).
#'
#' @return An object of class `ede_evi`, a data frame with columns `k`,
#'   `gamma`, `se`, `lower`, `upper` (normal interval for `gamma`) and
#'   `endpoint` (`NA` where the index is not negative), with attributes
#'   `alpha` and `origin`. `plot()` draws the estimates against `k`.
#'
#' @references
#' Dekkers, A. L. M., Einmahl, J. H. J., & de Haan, L. (1989). A moment
#'   estimator for the index of an extreme-value distribution. Annals of
#'   Statistics, 17(4), 1833-1855. doi:10.1214/aos/1176347397
#'
#' @seealso [dekkers1989()]
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' ede_evi(sd)
#' plot(ede_evi(sd))
#'
#' @export
ede_evi <- function(sd, kmin = 3L, kmax = NULL, events = FALSE, origin = NULL,
                    alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  t <- evt_times(sd, events)
  n <- length(t)
  origin <- moment_origin(sd, origin)
  n <- sum(t > origin)   # the threshold must be later than the origin
  kmax <- if (is.null(kmax)) n - 1L else kmax
  kmin <- evt_k(kmin, n, kmin = 3L, default = 3L, label = "Dekkers et al. (1989)",
                kmax = n - 1L)
  kmax <- evt_k(kmax, n, kmin = kmin, default = n - 1L,
                label = "Dekkers et al. (1989)", kmax = n - 1L)
  ks <- kmin:kmax
  z <- stats::qnorm(1 - alpha / 2)
  rows <- lapply(ks, function(k) {
    f <- moment_fit(t - origin, k)
    se <- if (is.na(f$gamma)) NA_real_ else sqrt(moment_var(f$gamma) / k)
    data.frame(k = k, gamma = f$gamma, se = se, lower = f$gamma - z * se,
               upper = f$gamma + z * se,
               endpoint = if (is.na(f$endpoint)) NA_real_
                          else origin + f$endpoint)
  })
  out <- do.call(rbind, rows)
  attr(out, "alpha") <- alpha
  attr(out, "origin") <- origin
  class(out) <- c("ede_evi", "data.frame")
  out
}

#' @export
print.ede_evi <- function(x, digits = 4L, ...) {
  cat(sprintf("<Moment estimates of the extreme-value index (origin %g)>\n",
              attr(x, "origin")))
  print.data.frame(as.data.frame(lapply(x, function(v)
    if (is.numeric(v)) signif(v, digits) else v)), row.names = FALSE)
  invisible(x)
}

#' @export
plot.ede_evi <- function(x, ...) {
  args <- utils::modifyList(
    list(xlab = "k", ylab = expression(hat(gamma)), type = "n",
         ylim = range(c(x$lower, x$upper, 0, -0.5), na.rm = TRUE),
         main = "Moment estimator of the extreme-value index"),
    list(...)
  )
  do.call(graphics::plot, c(list(x$k, x$gamma), args))
  graphics::polygon(c(x$k, rev(x$k)), c(x$lower, rev(x$upper)),
                    col = "grey88", border = NA)
  graphics::abline(h = 0, lty = 2)
  graphics::abline(h = -0.5, lty = 3, col = "grey50")
  graphics::lines(x$k, x$gamma, type = "b", pch = 19)
  invisible(x)
}

#' Origin of the time axis for the moment estimator
#' @keywords internal
#' @noRd
moment_origin <- function(sd, origin) {
  if (is.null(origin)) return(min(sd$time))
  if (!is.numeric(origin) || length(origin) != 1L || !is.finite(origin)) {
    stop("`origin` must be a single finite number.", call. = FALSE)
  }
  origin
}

#' Moment estimate of the index and end-point from positive data
#'
#' `x` are positive observations, latest (largest) first.
#' @keywords internal
#' @noRd
moment_fit <- function(x, k) {
  thr <- x[k + 1L]
  if (!(thr > 0)) {
    stop("the threshold sighting (the `k + 1`-th latest) must be later than ",
         "`origin`; use a smaller `k` or an earlier `origin`.", call. = FALSE)
  }
  l <- log(x[seq_len(k)] / thr)
  m1 <- mean(l)
  m2 <- mean(l^2)
  if (!(m2 > 0) || !is.finite(m1 / m2)) {
    return(list(gamma = NA_real_, endpoint = NA_real_))
  }
  gm <- 1 - 0.5 / (1 - m1^2 / m2)
  gamma <- m1 + gm
  endpoint <- if (gamma < 0) {
    thr - thr * m1 * (1 - gm) / gamma
  } else NA_real_
  list(gamma = gamma, endpoint = endpoint)
}

#' Asymptotic variance of the moment estimator (Dekkers et al., 1989)
#' @keywords internal
#' @noRd
moment_var <- function(g) {
  if (g >= 0) return(1 + g^2)
  (1 - g)^2 * (1 - 2 * g) * (1 - g + 6 * g^2) / ((1 - 3 * g) * (1 - 4 * g))
}
