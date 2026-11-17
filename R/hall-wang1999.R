#' Hall & Wang (1999) minimum-distance end-point estimator
#'
#' Estimates the extinction time as the end-point of the sighting-time
#' distribution, from the `k` most recent sightings, by the minimum-distance
#' method of Hall & Wang (1999). The `k` latest sightings are taken to follow
#' the limiting Weibull extreme-value model, the same model that underlies
#' [ole()]. If `t[1] > t[2] > ... > t[k]` are the sighting times, latest
#' first, and `theta` is a candidate end-point, then
#' `Y[i] = i * (log(theta - t[i + 1]) - log(theta - t[i]))` for
#' `i = 1, ..., k - 1`. At the true end-point the `Y[i]` are independent
#' exponential variables with a common mean. The estimate is the `theta` that
#' minimises Greenwood's statistic `G = sum(Y^2) / sum(Y)^2` of the `Y[i]`.
#'
#' Unlike maximum likelihood, the criterion does not involve the extremal
#' exponent (the shape of the Weibull model), which makes the estimator
#' well defined much more often than maximum likelihood when `k` is small.
#' For the same reason the distribution of `G` at the true end-point does not
#' depend on any nuisance parameter, so the confidence interval
#' `{theta : G(theta) < g}`, with `g` the `1 - alpha` quantile of that
#' distribution, has exact coverage when the model holds. The package
#' computes the quantile by simulation (`nsim` draws, with a fixed internal
#' seed, so the result is reproducible and the global random number stream is
#' left untouched). The distribution of Greenwood's statistic has no closed
#' form; with the default `nsim` the simulated quantiles agree to three
#' decimals with the exact points of Burrows (1979) for 10 uniform spacings
#' (as quoted by Stephens, 1981), and are more accurate than the Pearson-curve
#' points of Stephens (1981) in the lower tail for small samples. The Monte
#' Carlo error of a tail probability near 0.05 is about 0.0007.
#'
#' The estimator does not exist when `G` has no minimum below the end-point
#' at infinity, which happens for some records; `estimate` is then `NA` with
#' a warning. The interval is still reported. An upper limit that does not
#' exist is returned as `Inf` with a warning, as in [ole()].
#'
#' The model assumes distinct (continuous) sighting times, so by default each
#' distinct sighting time counts once, as in [robson1964()]; set
#' `events = TRUE` to count every sighting event instead, as [ole()] does.
#' With few sightings Hall & Wang recommend a small `k`, because the
#' distribution of the sighting times may change over a long record; by
#' default the 10 latest sightings are used, following the usual practice for
#' [ole()] (Solow 2005). [ede_select_k()] chooses `k` from the data.
#'
#' Hall & Wang (1999, section 3) report their example with the 90% interval;
#' it is reproduced with `alpha = 0.1`.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). The interval has coverage
#'   `1 - alpha`.
#' @param k Number of most recent sightings to use, an integer of at least 4,
#'   or `NULL` for the 10 latest (all, if there are fewer).
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#' @param nsim Number of draws used to compute the critical values of
#'   Greenwood's statistic.
#' @param data_out If `TRUE`, return instead the curve of the statistic and of
#'   the p-value of candidate end-points between the last sighting and
#'   `test_year` (the confidence interval is the set where `p_value > alpha`).
#' @param test_year Upper end of the curve, later than the last sighting;
#'   required when `data_out = TRUE`.
#'
#' @return An [ede_estimate] object with components `estimate`, `lower`,
#'   `upper`, plus `shape` (estimated extremal exponent), `k`, `statistic`
#'   (minimum of Greenwood's statistic) and `critical` (its critical value).
#'   With `data_out = TRUE`, a data frame with columns `time`, `statistic` and
#'   `p_value`.
#'
#' @references
#' Hall, P., & Wang, J. Z. (1999). Estimating the end-point of a probability
#'   distribution using minimum-distance methods. Bernoulli, 5(1), 177-189.
#'   doi:10.2307/3318618
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @seealso [ole()], which fits the same model by optimal linear estimation;
#'   [fraga_alves_neves2014()] for the Gumbel domain of attraction;
#'   [ede_select_k()] to choose `k`.
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' hall_wang1999(sd)
#' hall_wang1999(sd, k = 6, alpha = 0.1)
#'
#' @export
hall_wang1999 <- function(sd, alpha = 0.05, k = NULL, events = FALSE,
                          nsim = 1e5, data_out = FALSE, test_year = NULL) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  check_nsim(nsim)
  t <- evt_times(sd, events)
  k <- evt_k(k, length(t), kmin = 4L, default = min(length(t), 10L),
             label = "Hall & Wang (1999)")
  t <- t[seq_len(k)]
  if (t[1L] <= t[k]) {
    stop("Hall & Wang (1999) needs the `k` latest sightings to span more ",
         "than one time.", call. = FALSE)
  }
  g <- t[1L] - t
  null <- greenwood_null(k - 1L, nsim)

  if (isTRUE(data_out)) {
    check_test_year(test_year, t[1L])
    cand <- t[1L] + (test_year - t[1L]) * seq_len(500L) / 500
    stat <- hw_stat(cand - t[1L], g)
    p <- 1 - findInterval(stat, null, left.open = TRUE) / length(null)
    return(data.frame(time = cand, statistic = stat, p_value = p))
  }

  fit <- hw_fit(g)
  crit <- stats::quantile(null, 1 - alpha, names = FALSE)
  ci <- hw_interval(g, fit, crit)

  if (is.na(fit$s)) {
    warning("the minimum-distance estimator does not exist for this record ",
            "(Greenwood's statistic has no minimum before infinity); ",
            "returning NA.", call. = FALSE)
  }
  if (is.infinite(ci[2L])) {
    warning("the upper confidence limit does not exist at this alpha and ",
            "record; returning Inf.", call. = FALSE)
  }
  if (anyNA(ci)) {
    warning("no end-point is compatible with the record at this alpha; ",
            "the confidence interval is empty.", call. = FALSE)
  }
  out <- new_ede_estimate(
    estimate = t[1L] + fit$s, lower = t[1L] + ci[1L], upper = t[1L] + ci[2L],
    method = "Hall & Wang (1999)", alpha = alpha, interval_type = "two-sided"
  )
  out$shape <- fit$shape
  out$k <- k
  out$statistic <- fit$stat
  out$critical <- crit
  out
}

# ---------------------------------------------------------------------------
# Internals shared by the extreme-value estimators
# ---------------------------------------------------------------------------

#' Sighting times, latest first
#' @keywords internal
#' @noRd
evt_times <- function(sd, events = FALSE) {
  pos <- sd$count > 0
  if (!any(pos)) stop("no sightings with count > 0.", call. = FALSE)
  t <- if (isTRUE(events)) rep(sd$time[pos], sd$count[pos]) else sd$time[pos]
  sort(t, decreasing = TRUE)
}

#' Validate the number `k` of extreme sightings
#' @keywords internal
#' @noRd
evt_k <- function(k, n, kmin, default, label, kmax = n) {
  if (is.null(k)) k <- default
  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k != as.integer(k) ||
      k < kmin) {
    stop("`k` must be an integer of at least ", kmin, ".", call. = FALSE)
  }
  if (n < kmin) {
    stop(label, " needs at least ", kmin, " sighting times; the record has ",
         n, ".", call. = FALSE)
  }
  if (k > kmax) {
    stop("`k` cannot exceed ", kmax, " for this record and method.",
         call. = FALSE)
  }
  as.integer(k)
}

check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }
  invisible(alpha)
}

check_nsim <- function(nsim) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1000) {
    stop("`nsim` must be a number of at least 1000.", call. = FALSE)
  }
  invisible(nsim)
}

#' Greenwood's statistic as a function of the candidate end-point
#'
#' `s` is the distance of the candidate end-point beyond the latest sighting
#' (a vector); `g` are the gaps `t[1] - t[i]` of the `k` latest sightings.
#' @keywords internal
#' @noRd
hw_stat <- function(s, g) {
  k <- length(g)
  logs <- log(outer(s, g, "+"))
  y <- logs[, -1L, drop = FALSE] - logs[, -k, drop = FALSE]
  y <- sweep(y, 2L, seq_len(k - 1L), "*")
  rowSums(y^2) / rowSums(y)^2
}

#' Limit of the statistic as the end-point goes to infinity
#' @keywords internal
#' @noRd
hw_stat_inf <- function(g) {
  y <- seq_len(length(g) - 1L) * diff(g)
  sum(y^2) / sum(y)^2
}

#' Grid of end-point distances, relative to the span of the `k` sightings
#' @keywords internal
#' @noRd
hw_grid <- function(g) {
  g[length(g)] * exp(seq(log(1e-4), log(1e5), length.out = 400L))
}

#' Minimum of the statistic over the grid and the limit at infinity
#' @keywords internal
#' @noRd
hw_min_stat <- function(g) min(hw_stat(hw_grid(g), g), hw_stat_inf(g))

#' Fit the minimum-distance estimator
#' @keywords internal
#' @noRd
hw_fit <- function(g) {
  grid <- hw_grid(g)
  stat <- hw_stat(grid, g)
  inf <- hw_stat_inf(g)
  i <- which.min(stat)
  if (i > 1L && i < length(grid) && stat[i] < inf) {
    o <- stats::optimize(function(ls) hw_stat(exp(ls), g),
                         log(c(grid[i - 1L], grid[i + 1L])), tol = 1e-12)
    s <- exp(o$minimum)
    gmin <- o$objective
    logs <- log(s + g)
    shape <- (length(g) - 1L) /
      sum(seq_len(length(g) - 1L) * diff(logs))
  } else {
    s <- NA_real_
    gmin <- min(stat[i], inf)
    shape <- NA_real_
  }
  list(s = s, stat = gmin, shape = shape)
}

#' Confidence interval {s : G(s) < crit} as distances beyond the last sighting
#' @keywords internal
#' @noRd
hw_interval <- function(g, fit, crit) {
  grid <- hw_grid(g)
  stat <- hw_stat(grid, g)
  ok <- which(stat < crit)
  inf_ok <- hw_stat_inf(g) < crit
  if (!length(ok) && !inf_ok) return(c(NA_real_, NA_real_))
  # only the limit at infinity is below the critical value: the interval
  # starts beyond the end of the grid, which is returned as its lower end
  if (!length(ok)) return(c(grid[length(grid)], Inf))
  f <- function(ls) hw_stat(exp(ls), g) - crit
  first <- ok[1L]
  last <- ok[length(ok)]
  lower <- if (first == 1L) 0 else
    exp(stats::uniroot(f, log(grid[c(first - 1L, first)]), tol = 1e-12)$root)
  upper <- if (last == length(grid) || (inf_ok && last >= length(grid))) Inf else
    exp(stats::uniroot(f, log(grid[c(last, last + 1L)]), tol = 1e-12)$root)
  c(lower, upper)
}

greenwood_cache <- new.env(parent = emptyenv())

#' Null distribution of Greenwood's statistic for `m` exponential spacings
#'
#' `sum(E^2) / sum(E)^2` for `m` independent unit exponentials, simulated
#' with a fixed seed. The global random number stream is restored.
#' @keywords internal
#' @noRd
greenwood_null <- function(m, nsim) {
  key <- paste(m, nsim)
  if (!is.null(greenwood_cache[[key]])) return(greenwood_cache[[key]])
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20261112L)
  e <- matrix(stats::rexp(nsim * m), nsim, m)
  out <- sort(rowSums(e^2) / rowSums(e)^2)
  greenwood_cache[[key]] <- out
  out
}
