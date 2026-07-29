
#' Solow (1993) nonparametric persistence test
#'
#' Nonparametric test of the null hypothesis that a species was still extant
#' at a candidate test year, based on the ratio of the time since the last
#' sighting to the total observation window, under a homogeneous sighting
#' process.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). Persistence is rejected for
#'   the first candidate year at which the chance of persistence falls to
#'   or below `alpha`.
#' @param test_year Latest year to test. Must be supplied as a number.
#' @param data_out If `TRUE`, return the full chance-of-persistence curve
#'   instead of the single first-rejection year.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Solow, A. R. (1993). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' @export
solow1993 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  full <- expand_record(sd, test_year)

  times <- sort(unique(sd$time[sd$count > 0]))
  n_total <- sum(sd$count[sd$count > 0])
  if (n_total < 2) {
    stop("Solow (1993) needs at least 2 sightings with count > 0.", call. = FALSE)
  }

  t1 <- times[1]
  last_sight <- times[length(times)]
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  tn <- last_sight - t1

  candidates <- seq(last_sight + 1, test_year)
  chance <- (tn / (candidates - t1))^n_total

  if (data_out) return(data.frame(time = candidates, chance = chance))

  t_needed <- ceiling(tn / (alpha^(1 / n_total)))
  est_year <- t1 + t_needed

  if (est_year > test_year) {
    warning(
      "chance of persistence never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Solow (1993)", alpha = alpha))
  }

  new_ede_estimate(est_year, method = "Solow (1993)", alpha = alpha)
}

#' @keywords internal
#' @noRd
solow1993_chance <- function(d) {
  tmin <- min(d$time)
  tn <- max(d$time[d$count > 0]) - tmin
  tmax <- max(d$time) - tmin
  n <- sum(d$count)
  if (tmax <= 0) return(1.0)
  (tn / tmax)^n
}
