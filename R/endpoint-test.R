#' Test whether a sighting record is compatible with a finite end-point
#'
#' Extinction dates are the finite upper end-point of the sighting-time
#' distribution, but a finite end-point is not guaranteed: if the sighting
#' rate declines smoothly to zero (an exponential tail, the Gumbel domain of
#' attraction) or decays like a power law (the Frechet domain), the end-point
#' estimators for the Weibull domain, [ole()] and [hall_wang1999()], rest on a
#' model the data do not support. This function tests, from the `k` latest
#' sightings, the null hypothesis that the sighting times lie in the Gumbel
#' domain of attraction against the alternative that they lie in the Weibull
#' domain (a finite end-point, the setting of those estimators) or in the
#' Frechet domain (an unbounded upper tail).
#'
#' Three tests are available. All are invariant to the origin and the unit of
#' time, and all depend only on the latest sightings (the `k` latest, or the
#' `2k + 1` latest for `"rosario"`), taken in the same way as for
#' [hall_wang1999()].
#'
#' * `"neves"`: the statistic of Neves, Picek & Fraga Alves (2006), the ratio
#'   of the maximum to the mean of the excesses over the `(k + 1)`-th latest
#'   sighting, minus `log(k)`. With `t[1] > t[2] > ...` the sighting times,
#'   latest first, it is
#'   `T* = (t[1] - t[k + 1]) / mean(t[1:k] - t[k + 1]) - log(k)`. It needs
#'   `k + 1` sighting times. The paper rejects when `T*` falls outside the
#'   Gumbel quantiles. Here the p-value is exact for an exponential tail:
#'   `T* + log(k)` is then `k` times the maximum over the sum of `k`
#'   independent exponentials, whose distribution is Fisher's (1929) statistic
#'   for periodogram analysis. The column `p_asymptotic` gives the Gumbel
#'   approximation of the paper.
#' * `"hasofer_wang"`: the statistic of Hasofer & Wang (1992),
#'   `W = k (mean(x) - x[k])^2 / ((k - 1) sum((x - mean(x))^2))`, with `x` the
#'   `k` latest sightings. Under the null the distances `t[i] - t[k]` of the
#'   other `k - 1` sightings are distributed as `k - 1` independent
#'   exponentials, and `W = 1 / ((k - 1) (k G - 1))` with `G` Greenwood's
#'   statistic of those distances. The p-value is therefore taken from the
#'   same simulated null as [hall_wang1999()].
#' * `"rosario"`: the statistic of Fraga Alves, Neves & Rosario (2016), built on
#'   the end-point estimator of [fraga_alves_neves2014()]. With `xhat` that
#'   estimate for `k` and `t[1] > t[2] > ...` the sighting times, latest first,
#'   `G = (xhat - t[k + 1]) / (t[k + 1] - t[2 k + 1])` and the statistic is
#'   `G* = log(2) G - (log(k) + log(2) / 2)`, which has a Gumbel limit under
#'   the null (their Theorem 6) and goes to minus infinity for a finite
#'   end-point and plus infinity for an unbounded tail. It needs `2k + 1`
#'   sighting times. The paper rejects when `G*` is below the Gumbel quantile of
#'   level `alpha` (finite end-point); `p_asymptotic` is that Gumbel p-value.
#'   `p_value` is the p-value under an exponential tail, simulated from the
#'   same representation as for the other tests (given the `(2k + 1)`-th latest
#'   time, the `2k` later ones are the order statistics of `2k` independent
#'   exponentials), which is more reliable than the asymptotic one when `k` is
#'   small. Like `"neves"` it does not estimate the extreme value index.
#'
#' Small values of `T*` and `G*` and large values of `W` point to a finite
#' end-point.
#' The p-value is for `alternative`: `"finite"` (default) is the Weibull
#' domain, `"infinite"` the Frechet domain and `"two.sided"` either (twice the
#' smaller one-sided p-value, as in the paper's two-sided rule, capped at 1).
#'
#' A non-significant result does not show that the end-point is infinite: with
#' only a handful of sightings the tests have little power against a finite
#' end-point (see the simulations in Neves et al., 2006, and in Hasofer &
#' Wang, 1992). Use it to flag records for which the finite end-point model is
#' contradicted, not to confirm it.
#'
#' By default each distinct sighting time counts once, as in
#' [hall_wang1999()]; set `events = TRUE` to count every sighting event.
#'
#' @param sd A [sighting_data] object.
#' @param k Number of latest sightings, an integer of at least 4, or `NULL`
#'   for the 10 latest (fewer if the record is shorter; the Neves test needs
#'   one more sighting than `k`, and the Rosario test needs `2k + 1`
#'   sightings). [ede_select_k()] chooses `k` for the estimators of the
#'   Weibull domain.
#' @param test Which tests to run: any of `"neves"`, `"hasofer_wang"` and
#'   `"rosario"`; the first two by default.
#' @param alternative `"finite"` (default), `"infinite"` or `"two.sided"`; see
#'   Details.
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#' @param nsim Number of draws used for the null of Greenwood's statistic in
#'   the Hasofer & Wang test.
#'
#' @return An object of class `ede_endpoint_test`, a list with a data frame
#'   `table` (columns `test`, `k`, `statistic`, `p_value`, `p_asymptotic`),
#'   `alternative`, `k`, and the number of sighting times `n`.
#'
#' @references
#' Neves, C., Picek, J., & Fraga Alves, M. I. (2006). The contribution of the
#'   maximum to the sum of excesses for testing max-domains of attraction.
#'   Journal of Statistical Planning and Inference, 136(4), 1281-1301.
#'   doi:10.1016/j.jspi.2004.09.008
#'
#' Hasofer, A. M., & Wang, Z. (1992). A test for extreme value domain of
#'   attraction. Journal of the American Statistical Association, 87(417),
#'   171-177. doi:10.2307/2290466
#'
#' Fraga Alves, I., Neves, C., & Rosario, P. (2016). A general estimator for
#'   the right endpoint with an application to supercentenarian women's
#'   records. arXiv:1412.3972 (version 3).
#'
#' Fisher, R. A. (1929). Tests of significance in harmonic analysis.
#'   Proceedings of the Royal Society A, 125(796), 54-59.
#'
#' @seealso [hall_wang1999()] and [ole()] for the Weibull domain,
#'   [fraga_alves_neves2014()] for the Gumbel domain.
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' ede_endpoint_test(sd)
#' ede_endpoint_test(sd, k = 6, test = "neves", alternative = "two.sided")
#'
#' @export
ede_endpoint_test <- function(sd, k = NULL, test = c("neves", "hasofer_wang"),
                              alternative = c("finite", "infinite",
                                              "two.sided"),
                              events = FALSE, nsim = 1e5) {
  stopifnot(inherits(sd, "sighting_data"))
  test <- match.arg(test, c("neves", "hasofer_wang", "rosario"),
                    several.ok = TRUE)
  alternative <- match.arg(alternative)
  check_nsim(nsim)
  t <- evt_times(sd, events)
  n <- length(t)
  kmax <- n
  if ("neves" %in% test) kmax <- min(kmax, n - 1L)
  if ("rosario" %in% test) kmax <- min(kmax, (n - 1L) %/% 2L)
  k <- evt_k(k, n, kmin = 4L,
             default = min(kmax, 10L), label = "the end-point tests",
             kmax = kmax)
  rows <- list()
  if ("neves" %in% test) {
    r <- endpoint_neves(t, k, alternative)
    rows$neves <- data.frame(test = "Neves, Picek & Fraga Alves (2006)", k = k,
                             statistic = r$stat, p_value = r$p,
                             p_asymptotic = r$p_asy,
                             stringsAsFactors = FALSE)
  }
  if ("hasofer_wang" %in% test) {
    r <- endpoint_hasofer_wang(t, k, alternative, nsim)
    rows$hw <- data.frame(test = "Hasofer & Wang (1992)", k = k,
                          statistic = r$stat, p_value = r$p,
                          p_asymptotic = NA_real_, stringsAsFactors = FALSE)
  }
  if ("rosario" %in% test) {
    r <- endpoint_rosario(t, k, alternative, nsim)
    rows$rosario <- data.frame(test = "Fraga Alves, Neves & Rosario (2016)",
                               k = k, statistic = r$stat, p_value = r$p,
                               p_asymptotic = r$p_asy,
                               stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  structure(list(table = tab, alternative = alternative, k = k, n = n),
            class = "ede_endpoint_test")
}

#' @export
print.ede_endpoint_test <- function(x, digits = 4L, ...) {
  alt <- switch(x$alternative,
    finite = "finite end-point (Weibull domain)",
    infinite = "unbounded upper tail (Frechet domain)",
    two.sided = "not the Gumbel domain")
  cat("<Test of the Gumbel domain of attraction>\n")
  cat(sprintf("  alternative: %s\n  k = %d latest sightings, %d sighting times\n",
              alt, x$k, x$n))
  tab <- x$table
  tab$statistic <- signif(tab$statistic, digits)
  tab$p_value <- signif(tab$p_value, digits)
  tab$p_asymptotic <- signif(tab$p_asymptotic, digits)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' Neves, Picek & Fraga Alves (2006) statistic and p-values
#'
#' `t` are the sighting times, latest first.
#' @keywords internal
#' @noRd
endpoint_neves <- function(t, k, alternative) {
  thr <- t[k + 1L]
  exc <- t[seq_len(k)] - thr
  if (!(mean(exc) > 0)) {
    stop("the `k` latest sightings must span more than one time.",
         call. = FALSE)
  }
  tn <- exc[1L] / mean(exc)
  lo <- 1 - fisher_g_upper(tn / k, k)   # P(T <= t): small T means finite
  hi <- fisher_g_upper(tn / k, k)
  g <- tn - log(k)
  lo_a <- exp(-exp(-g))
  hi_a <- 1 - lo_a
  pick <- function(lo, hi) switch(alternative, finite = lo, infinite = hi,
                                  two.sided = min(1, 2 * min(lo, hi)))
  list(stat = g, p = pick(lo, hi), p_asy = pick(lo_a, hi_a))
}

#' P(max / sum > x) for `k` independent exponentials (Fisher, 1929)
#' @keywords internal
#' @noRd
fisher_g_upper <- function(x, k) {
  if (x >= 1) return(0)
  if (x <= 1 / k) return(1)
  j <- seq_len(floor(1 / x))
  p <- sum((-1)^(j - 1) * choose(k, j) * (1 - j * x)^(k - 1))
  min(1, max(0, p))
}

#' Statistic G* of Fraga Alves, Neves & Rosario (2016, eq. 17)
#'
#' `t` are the sighting times, latest first, at least `2 k + 1` of them.
#' @keywords internal
#' @noRd
rosario_stat <- function(t, k) {
  i <- 0:(k - 1L)
  a <- (log(k + i + 1L) - log(k + i)) / log(2)
  xhat <- t[1L] + sum(a * (t[k + 1L] - t[k + i + 1L]))
  g <- (xhat - t[k + 1L]) / (t[k + 1L] - t[2L * k + 1L])
  log(2) * g - (log(k) + log(2) / 2)
}

rosario_cache <- new.env(parent = emptyenv())

#' Null distribution of G* for an exponential tail
#'
#' Given the `(2 k + 1)`-th latest time, the `2 k` later ones are the order
#' statistics of `2 k` independent exponentials, so with that time at zero the
#' statistic is a function of exponentials alone. Simulated with a fixed seed;
#' the global random number stream is restored.
#' @keywords internal
#' @noRd
rosario_null <- function(k, nsim) {
  key <- paste(k, nsim)
  if (!is.null(rosario_cache[[key]])) return(rosario_cache[[key]])
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20261113L)
  m <- 2L * k
  e <- matrix(stats::rexp(nsim * m), nsim, m)
  out <- sort(vapply(seq_len(nsim), function(i) {
    rosario_stat(c(sort(e[i, ], decreasing = TRUE), 0), k)
  }, numeric(1)))
  rosario_cache[[key]] <- out
  out
}

#' Fraga Alves, Neves & Rosario (2016) statistic and p-values
#' @keywords internal
#' @noRd
endpoint_rosario <- function(t, k, alternative, nsim) {
  if (!(t[k + 1L] > t[2L * k + 1L])) {
    stop("the sightings from the (k + 1)-th to the (2k + 1)-th latest must ",
         "span more than one time.", call. = FALSE)
  }
  g <- rosario_stat(t, k)
  null <- rosario_null(k, nsim)
  lo <- mean(null <= g)   # small G* means a finite end-point
  hi <- mean(null >= g)
  lo_a <- exp(-exp(-g))   # Gumbel d.f.
  hi_a <- 1 - lo_a
  pick <- function(lo, hi) switch(alternative, finite = lo, infinite = hi,
                                  two.sided = min(1, 2 * min(lo, hi)))
  list(stat = g, p = pick(lo, hi), p_asy = pick(lo_a, hi_a))
}

#' Hasofer & Wang (1992) statistic and p-values
#' @keywords internal
#' @noRd
endpoint_hasofer_wang <- function(t, k, alternative, nsim) {
  x <- t[seq_len(k)]
  xbar <- mean(x)
  ss <- sum((x - xbar)^2)
  if (!(ss > 0)) {
    stop("the `k` latest sightings must span more than one time.",
         call. = FALSE)
  }
  w <- k * (xbar - x[k])^2 / ((k - 1) * ss)
  d <- x[-k] - x[k]
  gobs <- sum(d^2) / sum(d)^2
  null <- greenwood_null(k - 1L, nsim)
  lo <- mean(null <= gobs)   # small G = large W = finite end-point
  hi <- mean(null >= gobs)
  p <- switch(alternative, finite = lo, infinite = hi,
              two.sided = min(1, 2 * min(lo, hi)))
  list(stat = w, p = p)
}
