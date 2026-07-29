#' Solow (2005) sighting-effort-weighted persistence test
#'
#' Parametric extension of [solow1993()] that weights the observation
#' window by cumulative sighting effort (sighting count x time since first
#' sighting) instead of raw elapsed time, giving a more realistic null model
#' when sighting effort was not constant over the record.
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Solow, A. R. (2005). Inferring extinction from a sighting record.
#'   Mathematical Biosciences, 195(1), 47-55.
#'
#' @export
solow2005 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  full <- expand_record(sd, test_year)
  last_sight <- max(sd$time[sd$count > 0])
  candidates <- full$time[full$time > last_sight]
  if (length(candidates) == 0L) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  chance <- vapply(candidates, function(t) {
    solow2005_chance(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "chance of persistence never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Solow (2005)", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Solow (2005)", alpha = alpha)
}

#' @keywords internal
#' @noRd
solow2005_survival_series <- function(y, n) {
  if (y <= 0) return(0)
  n_terms <- floor(1 / y)
  if (n_terms <= 0) return(1)

  i <- seq_len(min(n_terms, n))
  terms <- (-1)^(i - 1) * choose(n, i) * (1 - i * y)^(n - 1)
  1 - sum(terms)
}

#' @keywords internal
#' @noRd
solow2005_chance <- function(d) {
  year <- d$time - min(d$time)
  s <- sum(year * d$count)
  tn <- max(year[d$count > 0])
  tmax <- max(year)

  n_total <- sum(d$count[d$count > 0])
  n_used <- n_total - 1

  if (s == 0 || n_used < 1) return(1.0)

  fs1 <- solow2005_survival_series(tn / s, n_used)
  fs2 <- solow2005_survival_series(tmax / s, n_used)

  if (fs2 <= 0) return(1.0)
  p <- fs1 / fs2
  max(0, min(1, p))
}
