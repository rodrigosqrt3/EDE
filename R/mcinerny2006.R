#' McInerny, Roberts, Davy & Cribb (2006) sighting-rate persistence test
#'
#' A modification of the Solow (1993a) persistence test that conditions on
#' the sighting rate observed up to the last sighting, `s = n / t[n]`, instead
#' of on the length of the whole observation window. For candidate time `T`,
#' the p-value is `(1 - s)^(T - t[n])` (equation 2 of the source paper). This
#' makes the test comparable across records with very different total
#' observation periods.
#'
#' The record is treated as binary: `n` counts time units with at least one
#' sighting, as in the source paper. The earliest supplied time defines the
#' observation origin. If that time has a sighting, it is used as the start
#' date and not counted, so `n` is reduced by one, as in the paper. If the
#' record starts with a zero-count row, that row is the known start date and
#' every sighting time is counted.
#'
#' The model is in whole time units, so the returned upper bound is the first
#' time unit after the last sighting at which persistence is rejected, which
#' is how the paper itself reports the test (Fig. 1).
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1).
#' @param test_year Latest time to test. Must be later than the last
#'   sighting. A bound beyond `test_year` is reported as `NA` with a warning.
#' @param data_out If `TRUE`, return the p-value curve instead of the upper
#'   bound.
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `chance` (p-values).
#'
#' @references
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#'   Significance of sighting rate in inferring extinction and threat.
#'   Conservation Biology, 20(2), 562-567.
#'
#' @export
mcinerny2006 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  check_test_year(test_year)

  times <- sort(unique(sd$time[sd$count > 0]))
  origin <- min(sd$time)
  # a sighting at the origin is the start date, not a tested sighting
  # (McInerny et al. 2006, p. 564)
  n_used <- sum(times > origin)
  if (n_used < 1) {
    stop(
      "McInerny et al. (2006) needs at least 1 sighting time after the ",
      "observation origin.", call. = FALSE
    )
  }

  last_sight <- max(times)
  check_test_year(test_year, last_sight)

  tn <- last_sight - origin
  rate <- n_used / tn
  if (rate >= 1) {
    warning(
      "estimated sighting rate n / tn >= 1 (a sighting in essentially ",
      "every time unit); the sighting-rate equation is degenerate for ",
      "this record and will reject persistence immediately.", call. = FALSE
    )
  }

  candidates <- candidate_times(last_sight, test_year)
  chance <- if (rate >= 1) rep(0, length(candidates)) else
    (1 - rate)^(candidates - last_sight)

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  bound <- if (length(below) > 0L) below[1] else NA_real_
  test_bound_result(bound, test_year, "McInerny, Roberts, Davy & Cribb (2006)", alpha)
}
