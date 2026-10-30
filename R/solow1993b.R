#' Solow (1993b) declining-population persistence test
#'
#' Conditional test of the null hypothesis that a declining species was still
#' extant at a candidate test year, under a non-stationary Poisson process with
#' an exponentially declining rate. The test conditions on `n` and on the sum
#' `S` of the sighting times, so the rate parameters are not estimated; the
#' p-value is `F(t[n]) / F(T)`, where `F` is Fisher's distribution of the
#' largest gap (equations 3-4 of the source paper).
#'
#' Sighting events are counted as in [solow1993()]: one sighting at the
#' observation origin defines it and is omitted, and every other sighting
#' event is counted, including repeated sightings in one time unit. Include a
#' zero-count first row when observation began before the first sighting, as
#' in the black-footed ferret example of the source paper.
#'
#' Solow (1993b) defines only the test. Under the assumption that the species
#' is extinct, the minimum variance unbiased estimate of the extinction time
#' under the same truncated exponential model (Beg, 1982) is equation 8 of
#' Solow (2005), also listed in Table 1 of McCrea et al. (2024):
#' `T_E = t[n] + sum_{i=0}^{[s/t[n]]} (-1)^i C(n, i) (s - i t[n])^(n-1) /
#' (n (n-1) sum_{i=0}^{[s/t[n]]-1} (-1)^i C(n-1, i) (s - (i+1) t[n])^(n-2))`,
#' which equals `t[n] + F(t[n]) / F'(t[n])`. It is returned as `estimate`
#' (in the original time units, i.e. shifted by the origin). Both sums are
#' evaluated as Fisher's largest-gap probabilities through the stable
#' B-spline (Irwin-Hall) recursion rather than the alternating series, which
#' loses all precision when the sighting times are concentrated near `t[n]`.
#' The estimate is `NA` when it is undefined (all sightings other than the
#' last at the origin with `n > 2`, or all at `t[n]`).
#'
#' Inverting the test numerically gives a one-sided upper confidence bound
#' for the extinction time, returned as `upper` (a package construct; Solow,
#' 2005, Section 3, describes the numerical inversion). Under this model the
#' bound may not exist: the p-value cannot fall below `F(t[n])`, however late
#' the candidate time (Solow, 2005, Section 3).
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object with the minimum variance unbiased
#'   estimate in `estimate` and the one-sided upper bound in `upper` (`NA`,
#'   with a warning, if it does not exist before `test_year`), or (if
#'   `data_out = TRUE`) a data frame with columns `time` and `chance`
#'   (p-values).
#'
#' @references
#' Solow, A. R. (1993b). Inferring extinction in a declining population.
#'   Journal of Mathematical Biology, 32(1), 79-82.
#'
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' Beg, M. A. (1982). Optimal tests and estimators for truncated exponential
#'   families. Metrika, 29, 103-113.
#'
#' Solow, A. R., & Helser, T. (2000). Detecting extinction in sighting data.
#'   In S. Ferson & M. Burgman (Eds.), Quantitative Methods for Conservation
#'   Biology (pp. 1-6). Springer, New York.
#'
#' McCrea, R. S., Cheale, T., Campillo-Funollet, E., & Roberts, D. L.
#'   (2024). Inferring species extinction from sighting data. Cambridge
#'   Prisms: Extinction, 2, e19.
#'
#' @export
solow1993b <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)

  last_sight <- max(sd$time[sd$count > 0])
  check_test_year(test_year, last_sight)

  origin <- min(sd$time)
  n <- solow_event_count(sd)
  s <- sum((sd$time - origin) * sd$count)
  if (n < 2 || last_sight <= origin) {
    stop("Solow (1993b) needs at least 2 sighting events after the observation origin.", call. = FALSE)
  }

  if (data_out) {
    full <- expand_record(sd, test_year)
    candidates <- full$time[full$time > last_sight]
    chance <- vapply(candidates, function(t) {
      solow1993b_chance(full[full$time <= t, , drop = FALSE])
    }, numeric(1))
    return(data.frame(time = candidates, chance = chance))
  }

  # F(T) reaches 1 at T = S, so the p-value is bounded below by F(t[n])
  tn <- last_sight - origin
  floor_p <- solow1993b_fisher_series(tn, s, n)
  bound <- if (floor_p > alpha) {
    NA_real_
  } else {
    origin + stats::uniroot(
      function(u) floor_p / solow1993b_fisher_series(u, s, n) - alpha,
      lower = tn, upper = s, tol = 1e-8
    )$root
  }
  res <- test_bound_result(bound, test_year, "Solow (1993b)", alpha)
  res$estimate <- origin + solow1993b_mvue(tn, s, n)
  res
}

#' Minimum variance unbiased estimate of the extinction time (Solow 2005,
#' eq. 8; Beg 1982), with times measured from the origin
#' @keywords internal
#' @noRd
solow1993b_mvue <- function(tn, s, n) {
  if (n == 2) return(tn + (2 * tn - s) / 2)
  if (s <= tn) return(NA_real_)
  # numerator sum = s^(n-1) F_n(tn / s); denominator sum =
  # (s - tn)^(n-2) F_{n-1}(tn / (s - tn)), with F Fisher's largest-gap cdf
  log_ratio <- (n - 2) * log(s / (s - tn)) +
    solow1993b_log_fisher(tn / s, n) -
    solow1993b_log_fisher(tn / (s - tn), n - 1)
  est <- tn + s / (n * (n - 1)) * exp(log_ratio)
  if (is.finite(est)) est else NA_real_
}

#' Log of Fisher's largest-gap probability F_m(y), the probability that the
#' largest of the m spacings of m - 1 uniform points on (0, 1) is at most y.
#' Uses F_m(y) = (m-1)! y^(m-1) M_m(1/y), with M_m the Irwin-Hall density
#' (cardinal B-spline), evaluated by the positive Cox-de Boor recursion
#' N_k(x) = x N_{k-1}(x) + (k - x) N_{k-1}(x - 1), N_k = (k-1)! M_k,
#' rescaled at each step to avoid overflow.
#' @keywords internal
#' @noRd
solow1993b_log_fisher <- function(y, m) {
  if (y >= 1) return(0)
  if (y * m <= 1) return(-Inf)
  z <- 1 / y - 0:(m - 1)
  v <- as.numeric(z >= 0 & z < 1)
  log_scale <- 0
  if (m >= 2) {
    for (k in 2:m) {
      idx <- seq_len(m - k + 1)
      v <- z[idx] * v[idx] + (k - z[idx]) * v[idx + 1]
      # 1 < 1/y < m puts some z in [0, 1), and the recursion keeps the
      # B-spline positive on its support, so top > 0
      top <- max(v)
      v <- v / top
      log_scale <- log_scale + log(top)
    }
  }
  (m - 1) * log(y) + log(v[1]) + log_scale
}

#' @keywords internal
#' @noRd
solow1993b_fisher_series <- function(u, s, n) {
  if (is.na(u) || is.nan(u) || u <= 0 || s <= 0) return(0)
  y <- u / s
  n_terms <- floor(1 / y)
  if (n_terms <= 0 || !is.finite(n_terms)) return(1)

  i <- seq_len(min(as.integer(n_terms), as.integer(n)))
  terms <- (-1)^(i - 1) * choose(n, i) * (1 - i * y)^(n - 1)
  1 - sum(terms)
}

#' @keywords internal
#' @noRd
solow1993b_chance <- function(d) {
  origin <- min(d$time)
  year <- d$time - origin
  s <- sum(year * d$count)
  positive <- year > 0 & d$count > 0
  if (!any(positive)) return(1.0)
  tn <- max(year[positive])
  tmax <- max(year)
  n <- solow_event_count(d)

  if (s == 0 || n < 2) return(1.0)

  fs1 <- solow1993b_fisher_series(tn, s, n)
  fs2 <- solow1993b_fisher_series(tmax, s, n)

  if (!is.finite(fs1) || !is.finite(fs2) || fs2 <= 0) return(1.0)

  p <- fs1 / fs2
  max(0, min(1, p))
}
