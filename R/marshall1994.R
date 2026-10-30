#' Marshall (1994) distribution-free confidence intervals on a range end-point
#'
#' Confidence interval on the end-point of a stratigraphic (or sighting) range
#' that does not assume randomly placed finds, only that gap sizes come from a
#' continuous distribution with no trend along the section. The `C`
#' confidence interval equals the `C`th percentile of the gap-size
#' distribution (Marshall, 1994, p. 460), which is bracketed by order
#' statistics of the `N` observed gaps. With `F(x)` the binomial
#' `Binomial(N, C)` distribution function and `gamma = (1 - conf_prob) / 2`:
#' the lower limit on the interval length is the `(x + 1)`th smallest gap,
#' `x` the largest integer `< N` with `F(x) <= gamma` (equation 1); the upper
#' limit is the `(x + 1)`th smallest gap, `x` the smallest integer with
#' `F(x) >= 1 - gamma` (equation 2). With probability `conf_prob` (the paper's
#' confidence probability `1 - 2 gamma`) the length of the `C` interval lies
#' between the two. If no `x < N` satisfies equation 1 there is no lower
#' limit; if only `x = N` satisfies equation 2 there is no upper limit and the
#' interval is unbounded (e.g. a 95% interval needs at least 72 gaps at
#' `conf_prob = 0.95`).
#'
#' The method gives a range of possible interval lengths rather than one
#' interval. The returned `upper` is the conservative end, last find plus the
#' upper-limit gap (`Inf` if there is none); `upper_min` is the last find
#' plus the lower-limit gap, so the `C` bound lies in `[upper_min, upper]`
#' with probability `conf_prob`. This mapping onto [ede_estimate] is a
#' package construct.
#'
#' Gaps are the distances between successive find horizons. By default
#' (`horizons = "distinct"`) each distinct time with a positive count is one
#' horizon, as in [strauss1989()]. With `horizons = "count"` a row with count
#' `k` stands for `k` horizons at the same (binned) time, i.e. `k - 1` gaps of
#' size zero, as in the paper's Metrarabdotos ages binned to 0.05 Myr.
#' Zero-count rows are not used. The end point is the one at the largest
#' `time`: for a first appearance, or for ages in Ma, code `time = -age` (or
#' `-height`). The method assumes continuous sampling; it is not meant for
#' discretely sampled cores (use [weiss_marshall1999()] for those). Check for
#' a trend in gap size with [marshall1994_trend()] first.
#'
#' @param sd A [sighting_data] object.
#' @param alpha One minus the confidence level `C` of the interval, in (0, 1).
#' @param conf_prob Confidence probability `1 - 2 gamma` with which the
#'   length of the `C` interval is bracketed, in (0, 1). The paper uses 0.95
#'   or 0.99.
#' @param horizons `"distinct"` (each distinct positive-count time is one
#'   horizon) or `"count"` (a count of `k` is `k` horizons at that time).
#'
#' @return An [ede_estimate] object with `estimate = NA`, `lower = NA`,
#'   `upper` (last find plus upper-limit gap, `Inf` if not available) and
#'   `interval_type = "one-sided"`, plus components `upper_min` (last find
#'   plus lower-limit gap, `NA` if not available), `ci_length` (named
#'   `lower`/`upper` limits on the interval length), `gap_ranks` (ranks `l`
#'   and `u` of the gaps used, `NA` if not available), `n_gaps`,
#'   `conf_prob` and `last_sighting`.
#'
#' @references
#' Marshall, C. R. (1994). Confidence intervals on stratigraphic ranges:
#'   partial relaxation of the assumption of randomly distributed fossil
#'   horizons. Paleobiology, 20(4), 459-469.
#'
#' @export
marshall1994 <- function(sd, alpha = 0.05, conf_prob = 0.95,
                         horizons = c("distinct", "count")) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (!is.numeric(conf_prob) || length(conf_prob) != 1L || is.na(conf_prob) ||
      conf_prob <= 0 || conf_prob >= 1) {
    stop("`conf_prob` must be in (0, 1).", call. = FALSE)
  }
  horizons <- match.arg(horizons)
  gaps <- marshall1994_gaps(sd, horizons)
  n <- length(gaps)
  conf_level <- 1 - alpha
  gamma <- (1 - conf_prob) / 2

  ranks <- marshall1994_ranks(n, conf_level, gamma)
  sorted <- sort(gaps)
  len_lo <- if (is.na(ranks[["lower"]])) NA_real_ else sorted[ranks[["lower"]]]
  len_hi <- if (is.na(ranks[["upper"]])) Inf else sorted[ranks[["upper"]]]
  lad <- max(sd$time[sd$count > 0])

  result <- new_ede_estimate(NA_real_, upper = lad + len_hi,
                             method = "Marshall (1994)", alpha = alpha,
                             interval_type = "one-sided")
  result$upper_min <- lad + len_lo
  result$ci_length <- c(lower = len_lo, upper = len_hi)
  result$gap_ranks <- ranks
  result$n_gaps <- n
  result$conf_prob <- conf_prob
  result$last_sighting <- lad
  result
}

#' Marshall (1994) transposition test for a trend in gap size
#'
#' Tests the assumption of [marshall1994()] (and of the classical intervals)
#' that gap size is not correlated with stratigraphic position. Following
#' Marshall (1994, p. 462-463), the gaps are listed in stratigraphic order and
#' the number of transpositions `T` (pairs in which a larger gap precedes a
#' smaller one) is counted; each group of `m` equal gaps adds `m (m - 1) / 4`.
#' `T` near 0 indicates gaps increasing in size, `T` near
#' `T_max = N (N - 1) / 2` gaps decreasing in size.
#'
#' The paper compares `min(T, T_max - T)` with tabulated critical values
#' (its Table 3, from Kreyszig 1979). Here a p-value is computed instead from
#' the exact null distribution of the number of transpositions in a random
#' ordering of `N` distinct gaps: `p = min(1, 2 P(I <= min(T, T_max - T)))`,
#' with a non-integer `T` (from ties) rounded down. This two-sided p-value is
#' a package construct. Table 3's critical values are two-sided (each tail
#' `p / 2`): the exact distribution reproduces 67 of its 72 entries, and the
#' other five differ by one at borderline probabilities.
#'
#' @inheritParams marshall1994
#'
#' @return A list of class `marshall1994_trend` with components
#'   `transpositions` (`T`, including the tie correction), `t_max`,
#'   `n_gaps`, `direction` (`"increasing"` if `T < T_max / 2`, `"decreasing"`
#'   if `T > T_max / 2`, otherwise `"none"`), `p_value` and `gaps` (in
#'   stratigraphic order).
#'
#' @references
#' Marshall, C. R. (1994). Confidence intervals on stratigraphic ranges:
#'   partial relaxation of the assumption of randomly distributed fossil
#'   horizons. Paleobiology, 20(4), 459-469.
#'
#' @export
marshall1994_trend <- function(sd, horizons = c("distinct", "count")) {
  stopifnot(inherits(sd, "sighting_data"))
  horizons <- match.arg(horizons)
  gaps <- marshall1994_gaps(sd, horizons)
  n <- length(gaps)
  if (n < 2L) stop("the transposition test needs at least 2 gaps.", call. = FALSE)

  # pairs i < j with gaps[i] > gaps[j]
  larger_first <- outer(gaps, gaps, ">")
  tp <- sum(larger_first[upper.tri(larger_first)])
  ties <- table(gaps)
  tp <- tp + sum(ties * (ties - 1) / 4)

  t_max <- n * (n - 1) / 2
  dist <- marshall1994_inversions(n)
  k <- floor(min(tp, t_max - tp))
  p_value <- min(1, 2 * sum(dist[seq_len(k + 1L)]))
  direction <- if (tp < t_max / 2) "increasing" else if (tp > t_max / 2) "decreasing" else "none"

  structure(
    list(transpositions = tp, t_max = t_max, n_gaps = n,
         direction = direction, p_value = p_value, gaps = gaps),
    class = "marshall1994_trend"
  )
}

#' @export
print.marshall1994_trend <- function(x, ...) {
  cat("<Marshall (1994) transposition test>\n")
  cat(sprintf("  %d gaps, T = %s (T_max = %s), suggested trend: %s\n",
              x$n_gaps, format(x$transpositions), format(x$t_max), x$direction))
  cat(sprintf("  two-sided p-value: %s\n", format(signif(x$p_value, 4))))
  invisible(x)
}

#' Gaps between successive find horizons, in stratigraphic order
#' @keywords internal
#' @noRd
marshall1994_gaps <- function(sd, horizons) {
  occ <- sd[sd$count > 0, , drop = FALSE]
  times <- if (horizons == "count") rep(occ$time, occ$count) else unique(occ$time)
  times <- sort(times)
  if (length(times) < 2L) {
    stop("Marshall (1994) needs at least 2 find horizons (1 gap).", call. = FALSE)
  }
  gaps <- diff(times)
  # gaps equal up to floating-point error (e.g. 1.40 - 1.35 vs 0.05 - 0)
  # are made exactly equal, so ties are recognized
  tol <- sqrt(.Machine$double.eps) * max(1, abs(times))
  o <- order(gaps)
  group <- cumsum(c(TRUE, diff(gaps[o]) > tol))
  gaps[o] <- stats::ave(gaps[o], group, FUN = min)
  gaps
}

#' Ranks of the gaps bounding the C-th percentile (equations 1 and 2)
#' @keywords internal
#' @noRd
marshall1994_ranks <- function(n, conf_level, gamma) {
  cdf <- stats::pbinom(0:n, n, conf_level)
  # tolerance guards against rounding at exact boundaries
  eps <- 1e-12
  lower_x <- which(cdf[seq_len(n)] <= gamma + eps) - 1L
  lower <- if (length(lower_x)) max(lower_x) + 1L else NA_integer_
  upper_x <- min(which(cdf >= 1 - gamma - eps)) - 1L
  upper <- if (upper_x < n) upper_x + 1L else NA_integer_
  c(lower = as.integer(lower), upper = as.integer(upper))
}

#' Null distribution of the number of inversions in a random permutation
#' @keywords internal
#' @noRd
marshall1994_inversions <- function(n) {
  d <- 1
  if (n >= 2L) {
    for (k in 2:n) {
      nd <- numeric(length(d) + k - 1L)
      for (s in 0:(k - 1L)) {
        idx <- seq_along(d) + s
        nd[idx] <- nd[idx] + d
      }
      d <- nd / k
    }
  }
  d
}
