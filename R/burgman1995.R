#' Burgman, Grimson & Ferson (1995) combinatorial persistence test
#'
#' Computes the probability that, if sighting events were distributed
#' uniformly at random over the candidate observation window, the largest
#' observed gap between sightings would not exceed the gap actually seen in
#' the data. Uses an inclusion-exclusion (Stirling-number) argument on the
#' occupancy of time bins by sighting events.
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring threat
#'   from scientific collections. Conservation Biology, 9(4), 923-928.
#'
#' @export
burgman1995 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
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
    burgman_chance(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  # sExtinct additionally requires the chance curve to be strictly
  # decreasing at the cut point (`diff(chance) < 0`), since this
  # combinatorial probability is not monotone in general.
  decreasing <- c(FALSE, diff(chance) < 0)
  below <- candidates[decreasing & chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "chance of persistence never falls to alpha (on a decreasing run) ",
      "before `test_year`; returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Burgman, Grimson & Ferson (1995)", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Burgman, Grimson & Ferson (1995)", alpha = alpha)
}

#' @keywords internal
#' @noRd
stirling2 <- function(n, j) {
  i <- 0:j
  sum((-1)^i * choose(j, i) * (j - i)^n) / factorial(j)
}

#' @keywords internal
#' @noRd
burgman_chance <- function(d) {
  n <- sum(d$count)   # total sighting events
  ct <- nrow(d)       # number of time bins from first sighting to cutoff

  nz_times <- d$time[d$count > 0]
  last_time <- d$time[nrow(d)]
  r <- max(diff(sort(c(nz_times, last_time))))

  total <- 0
  for (j in seq_len(n)) {
    s2 <- stirling2(n, j)
    k_max <- min(j + 1, floor(ct / r))
    if (k_max >= 1) {
      for (k in seq_len(k_max)) {
        fallfac <- prod(ct - r * k - (0:(j - 1)))
        total <- total + (-1)^(k + 1) * choose(j + 1, k) * fallfac * s2
      }
    }
  }
  ct^(-n) * total
}
