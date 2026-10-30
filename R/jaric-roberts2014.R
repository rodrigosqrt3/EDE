#' Reliability-adjusted extinction inference
#'
#' Modifies the constant-rate Solow (1993a) model by assigning a reliability
#' `p[i]` in `[0, 1]` to each sighting. The number of sightings is replaced by
#' the effective number `r = sum(p[i])` (equation 4 of Jaric & Roberts, 2014)
#' and the last sighting time by its expectation `t[r]` over which record is
#' the last valid one (equations 5-6). The point estimate is
#' `(r + 1) / r * t[r]` (equation 9) and the upper confidence bound is
#' `t[r] / alpha^(1 / r)` (equation 10). For a candidate time `T`, the
#' modified Solow p-value is `(t[r] / T)^r` (equations 7-8), with times
#' measured from the observation origin.
#'
#' The earliest supplied time defines the observation origin. A sighting at
#' that time fixes the origin and does not enter `r` or `t[r]`, so its
#' reliability is not used; prepend a zero-count row to include it.
#' Reliabilities given for zero-count rows are ignored.
#'
#' @param sd A [sighting_data] object. Counts must be zero or one. This is a
#'   restriction of the package, not of the paper: reliability belongs to an
#'   individual record, and the paper does not say how to combine several
#'   records of different reliability in one time unit.
#' @param reliability Numeric reliabilities in `[0, 1]`, supplied either for
#'   every row of `sd` or for the positive-count rows only.
#' @param alpha Significance level, in (0, 1), for the upper confidence bound.
#' @param test_year Optional candidate time at which to compute the p-value
#'   of equations 7-8. Must be later than the last sighting. Required when
#'   `data_out = TRUE`.
#' @param data_out If `TRUE`, return the p-value curve at integer steps
#'   after the last sighting up to `test_year` instead of the estimate.
#'
#' @return An [ede_estimate] object with additional components
#'   `effective_sightings` (`r`), `effective_endpoint` (`t[r]`, on the
#'   original time scale) and, if `test_year` is given, `p_value`; or (if
#'   `data_out = TRUE`) a data frame with columns `time` and `chance`
#'   (p-values).
#'
#' @references
#' Jaric, I., & Roberts, D. L. (2014). Accounting for observation reliability
#'   when inferring extinction based on sighting records. Biodiversity and
#'   Conservation, 23(11), 2801-2815.
#'
#' @export
jaric_roberts2014 <- function(sd, reliability, alpha = 0.05,
                              test_year = NULL, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (data_out && is.null(test_year)) {
    stop("`test_year` is required when `data_out = TRUE`.", call. = FALSE)
  }
  if (!is.null(test_year)) check_test_year(test_year)
  if (missing(reliability) || !is.numeric(reliability) || anyNA(reliability) ||
      any(!is.finite(reliability)) || any(reliability < 0 | reliability > 1)) {
    stop("`reliability` must contain finite values in [0, 1].", call. = FALSE)
  }
  if (any(sd$count > 1)) {
    stop("jaric_roberts2014() requires binary sighting counts (one record per row).",
         call. = FALSE)
  }

  occupied <- which(sd$count > 0)
  if (length(reliability) == nrow(sd)) {
    event_reliability <- reliability[occupied]
  } else if (length(reliability) == length(occupied)) {
    event_reliability <- reliability
  } else {
    stop("`reliability` must have one value per row or per occupied sighting time.",
         call. = FALSE)
  }

  event_times <- sd$time[occupied]
  origin <- min(sd$time)
  used <- event_times > origin
  event_times <- event_times[used] - origin
  event_reliability <- event_reliability[used]
  if (length(event_times) == 0L) {
    stop("at least one sighting is required after the observation origin.",
         call. = FALSE)
  }

  r <- sum(event_reliability)
  if (r <= 0) {
    stop("the sum of reliabilities after the observation origin must be positive.",
         call. = FALSE)
  }

  n <- length(event_times)
  last_valid_probability <- vapply(seq_len(n), function(i) {
    later_false <- if (i == n) 1 else prod(1 - event_reliability[(i + 1L):n])
    event_reliability[i] * later_false
  }, numeric(1))
  tr <- sum(event_times * last_valid_probability)

  if (!is.null(test_year)) {
    last_sight <- origin + max(event_times)
    check_test_year(test_year, last_sight)
    candidates <- candidate_times(last_sight, test_year)
    if (data_out) {
      return(data.frame(time = candidates,
                        chance = (tr / (candidates - origin))^r))
    }
  }

  estimate <- origin + ((r + 1) / r) * tr
  upper <- origin + tr / alpha^(1 / r)
  result <- new_ede_estimate(
    estimate,
    lower = NA_real_,
    upper = upper,
    method = "Jaric & Roberts (2014)",
    alpha = alpha,
    interval_type = "one-sided"
  )
  result$effective_sightings <- r
  result$effective_endpoint <- origin + tr
  if (!is.null(test_year)) result$p_value <- (tr / (test_year - origin))^r
  result
}
