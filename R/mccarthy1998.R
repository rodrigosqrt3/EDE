#' McCarthy (1998) persistence test adjusted for sampling effort
#'
#' Persistence test of Solow (1993a) in which the time axis is replaced by the
#' cumulative sampling effort, the "partial Solow equation" of McCarthy (1998,
#' equation 2). If the probability that a record falls in a given time unit is
#' proportional to the effort `e[i]` spent in it, then, under the hypothesis
#' that the species persisted to the candidate time `T`, the probability that
#' all `N` records fall no later than the last record `t[n]` is
#' `p = (sum(e[0:t[n]]) / sum(e[0:T]))^N`, where the sums run over the time
#' units from the observation origin. With constant effort it is the equation
#' of [solow1993()] up to the discretisation (`(t + 1) / (T + 1)` in place of
#' `t / T`, because the origin unit counts). Roberts & Jaric (2020) use the
#' same test for species known from a single specimen, with the number of
#' specimens collected as the effort index.
#'
#' The count `N` and the origin are those of [solow1993()]: the earliest time
#' of `sd` defines the origin, one sighting there is omitted, and every other
#' sighting event is counted. The time unit is one whole number of the time
#' scale (typically a year), as in McCarthy's data; `sd$time` must be integers
#' and effort must be given for every time unit from the origin to `test_year`.
#' Inverting the test gives the one-sided upper bound reported in `upper`: the
#' earliest time unit after the last sighting at which `p <= alpha`.
#'
#' Effort is only a proxy, and McCarthy and Roberts & Jaric report that the
#' method is well calibrated when the extinction is sudden and effort is stable
#' or rising, but biased when the decline is gradual or effort falls. The
#' choice of the effort index (for instance the number of specimens of all
#' species, or of related species, collected per year) is part of the analysis.
#'
#' @param sd A [sighting_data] object whose times are whole numbers.
#' @param effort The sampling effort in each time unit, non-negative and not
#'   all zero: either a numeric vector with one value for each time unit from
#'   the origin (`min(sd$time)`) to at least `test_year`, in order, or a
#'   vectorised function of time returning the effort in each time unit.
#' @param alpha Significance level, in (0, 1). The upper bound is the earliest
#'   time unit at which persistence is rejected at this level.
#' @param test_year Latest time to test, a whole number later than the last
#'   sighting. A bound beyond `test_year` is reported as `NA` with a warning;
#'   extend `effort` and `test_year` to project further.
#' @param data_out If `TRUE`, return the p-value at each time unit after the
#'   last sighting up to `test_year` instead of the upper bound.
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `p_value` (plus the deprecated alias `chance`).
#'
#' @references
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#'   museum data. Biological Conservation, 83(1), 9-17.
#'   doi:10.1016/S0006-3207(97)00048-7
#'
#' Roberts, D. L., & Jaric, I. (2020). Inferring the extinction of species
#'   known only from a single specimen. Oryx, 54(2), 161-166.
#'
#' Solow, A. R. (1993a). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' @seealso [solow1993()], the test with constant effort; [marshall1997()] for
#'   stratigraphic ranges with a recovery function.
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925),
#'                                sightings = 1))
#' # effort doubles after 1915
#' effort <- ifelse(1900:1950 < 1915, 1, 2)
#' mccarthy1998(sd, effort, test_year = 1950)
#' solow1993(sd, test_year = 1950)
#' mccarthy1998(sd, effort, test_year = 1950, data_out = TRUE)[1:3, ]
#'
#' @export
mccarthy1998 <- function(sd, effort, alpha = 0.05, test_year,
                         data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)
  if (any(sd$time != round(sd$time))) {
    stop("`sd$time` must be whole numbers (time units).", call. = FALSE)
  }
  origin <- min(sd$time)
  last_sight <- max(sd$time[sd$count > 0])
  if (last_sight <= origin) {
    stop("McCarthy (1998) needs at least 1 sighting after the observation ",
         "origin.", call. = FALSE)
  }
  check_test_year(test_year, last_sight)
  if (test_year != round(test_year)) {
    stop("`test_year` must be a whole number.", call. = FALSE)
  }
  units <- origin:test_year
  e <- if (is.function(effort)) effort(units) else effort
  if (!is.numeric(e) || anyNA(e) || any(!is.finite(e)) || any(e < 0)) {
    stop("`effort` must be non-negative and finite.", call. = FALSE)
  }
  if (length(e) < length(units)) {
    stop("`effort` must cover every time unit from ", origin, " to ",
         test_year, " (", length(units), " values); it has ", length(e), ".",
         call. = FALSE)
  }
  e <- e[seq_along(units)]
  cum <- cumsum(e)
  e_last <- cum[units == last_sight]
  if (!(e_last > 0)) {
    stop("no sampling effort up to the last sighting.", call. = FALSE)
  }
  n <- solow_event_count(sd)
  cand <- units[units > last_sight]
  p <- (e_last / cum[units > last_sight])^n

  if (isTRUE(data_out)) return(p_curve(cand, p))

  hit <- which(p <= alpha)
  bound <- if (length(hit)) cand[hit[1L]] else NA_real_
  test_bound_result(bound, test_year, "McCarthy (1998)", alpha)
}
